package org.metagene.ftdbexp.store;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import org.metagene.genestrip.DefaultExecutionContext;
import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResourceListStream;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.match.FastqKMerMatcher;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.match.CountsPerTaxid;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

/**
 * Compares the two k-mer stores of Genestrip on one and the same database: the radix store a
 * database is filled into, and a {@link KMerSortedArray} copied from it, which is the binary search
 * layout of the first study. Both are held in one JVM and both classify the same fastq file, so they
 * differ in the data structure a lookup walks and in nothing else.
 * <p>
 * Reported per store, and written as a row of {@code results/storebench.csv}:
 * <ul>
 * <li>the entries and the capacity the store reserved for them,</li>
 * <li>what the store and its Bloom filter take by construction, eight bytes per k-mer for the radix
 *     store and ten for the sorted array, which needs a second array for the taxon index,</li>
 * <li>the heap the store took when it was built, measured after a collection,</li>
 * <li>the wall time and the read throughput of classifying the file.</li>
 * </ul>
 * The memory figures are the point of the comparison. The throughput depends on how much of a run is
 * lookup at all: a database that classifies two percent of a short-read file spends its time reading,
 * not searching.
 */
public class StoreBenchMain {
    private static final File BASE_DIR = new File("./data");
    private static final File RESULTS_DIR = new File("./results");

    /** The store the database file brings. */
    private static final String STORE_DB = "db";
    /** The sorted array copied from it. */
    private static final String STORE_SORTED_ARRAY = "sortedarray";
    /** System property for the number of consumer threads. */
    private static final String THREADS_PROP = "ftdbexp.storebench.threads";
    /** System property for the discarded warm-up pass; {@code false} turns it off. */
    private static final String WARMUP_PROP = "ftdbexp.storebench.warmup";

    /**
     * Runs the comparison.
     *
     * @param args the database project name, the database file or an empty string for the project's
     *             refined one, and the fastq files to classify as one comma separated argument
     * @throws Exception if the database cannot be read or a fastq file cannot be classified
     */
    public static void main(String[] args) throws Exception {
        if (args.length < 3) {
            System.err.println("Usage: StoreBenchMain <db> <database file> <fastq file>[,<fastq file>...]");
            System.err.println("  <db>             name of the database project under data/projects");
            System.err.println("  <database file>  empty for data/projects/<db>/db/<db>_ftdb.zip");
            System.err.println("  <fastq file>     the files to classify, absolute or under data/fastq,");
            System.err.println("                   comma separated; each is classified with either store");
            System.exit(1);
        }
        String db = args[0];
        File dbFile = args[1].isEmpty() ? new File(BASE_DIR, "projects/" + db + "/db/" + db + "_ftdb.zip")
                : new File(args[1]);
        List<File> fastqs = new ArrayList<File>();
        for (String name : args[2].split(",")) {
            if (name.trim().isEmpty()) {
                continue;
            }
            File fastq = resolveFastq(db, name.trim());
            if (!fastq.exists()) {
                throw new IllegalArgumentException("No fastq file at " + fastq);
            }
            fastqs.add(fastq);
        }
        if (fastqs.isEmpty()) {
            throw new IllegalArgumentException("No fastq file named.");
        }
        if (!dbFile.exists()) {
            throw new IllegalArgumentException("No database at " + dbFile);
        }

        GSProject project = new GSProject(new GSCommon(BASE_DIR), db, true);
        project.initConfigParam(GSConfigKey.THREADS, -1);
        // One consumer thread per processor less one by default, as the goals do it. Fewer is what
        // shows the store: with many consumers the single producer thread - reading, inflating and
        // parsing - is the limit of the pipeline, and a faster lookup then barely moves the total.
        int configured = Integer.getInteger(THREADS_PROP, -1);
        int threads = configured < 0 ? Runtime.getRuntime().availableProcessors() - 1 : configured;
        if (threads < 1) {
            throw new IllegalArgumentException("At least one consumer thread is needed, got " + configured + ".");
        }
        boolean useFilter = project.booleanConfigValue(GSConfigKey.USE_BLOOM_FILTER_FOR_MATCH);
        // Each file is classified twice per store and only the second pass is measured. The first one
        // warms the JIT and the page cache, and it does so for either store alike: without it the
        // store measured first reads a cold file and compiles the matching loop, and the one measured
        // second gets both for free.
        boolean warmup = !"false".equalsIgnoreCase(System.getProperty(WARMUP_PROP, "true"));

        System.out.println("Genestrip k-mer store comparison");
        System.out.println("  project:  " + db);
        System.out.println("  database: " + dbFile);
        for (File fastq : fastqs) {
            System.out.println("  input:    " + fastq + " (" + String.format(Locale.US, "%.0f", toMB(fastq.length()))
                    + " MB compressed)");
        }
        System.out.println("  threads:  " + threads + " consumer and one producer, of "
                + Runtime.getRuntime().availableProcessors() + " processors");
        System.out.println("  warm-up:  " + (warmup ? "one discarded pass per file and store" : "none"));
        System.out.println();

        long heapBeforeLoad = usedHeapAfterGc();
        System.out.println("Loading " + dbFile + " ...");
        Database database = Database.load(dbFile, useFilter);
        // The converted store shares its arrays with the loaded one, so what the load added to the
        // heap is the store, its filter and the taxonomy tree together.
        KMerStore<SmallTaxIdNode> loaded = database.convertKMerStore();
        loaded.setUseFilter(useFilter);
        long loadHeap = usedHeapAfterGc() - heapBeforeLoad;

        List<Row> rows = new ArrayList<Row>();
        for (String kind : new String[] { STORE_DB, STORE_SORTED_ARRAY }) {
            KMerStore<SmallTaxIdNode> store;
            long heap;
            if (STORE_SORTED_ARRAY.equals(kind)) {
                long before = usedHeapAfterGc();
                long entries = loaded.getEntries();
                System.out.println("Copying " + entries + " entries into a sorted array ...");
                store = KMerSortedArray.copyOf(loaded,
                        project.doubleConfigValue(GSConfigKey.OPT_BLOOM_FILTER_FPP));
                // The copy shares nothing with the source, so this delta is the copy alone.
                heap = usedHeapAfterGc() - before;
                if (store.getEntries() != entries) {
                    throw new IllegalStateException("The copy holds " + store.getEntries()
                            + " entries against the source's " + entries + ".");
                }
                System.out.println("Copied: " + store.getEntries() + " entries.");
            } else {
                store = loaded;
                heap = loadHeap;
            }
            Row row = new Row(db, kind, store, heap);
            row.print();
            // Every file with every store, and the row holds the totals: the comparison is of the two
            // stores on one input set, so what counts is the time that set took, not a single file's.
            for (File fastq : fastqs) {
                classify(row, fastq, project, database, store, threads, warmup);
            }
            if (fastqs.size() > 1) {
                System.out.println("    " + fastqs.size() + " files together: "
                        + String.format(Locale.US, "%.2f", row.seconds) + " s, "
                        + String.format(Locale.US, "%.0f", row.readsPerSecond()) + " reads / s");
            }
            rows.add(row);
        }

        write(rows, db);
    }

    /**
     * Classifies the file with the given store and records the timing of the measured pass in the row.
     *
     * @param warmup whether to classify the file once before the measured pass and discard that
     */
    private static void classify(Row row, File fastq, GSProject project, Database database,
            KMerStore<SmallTaxIdNode> store, int threads, boolean warmup) throws Exception {
        if (warmup) {
            long discarded = runOnce(fastq, project, database, store, threads)[0];
            System.out.println("    " + fastq.getName() + " warm-up: "
                    + String.format(Locale.US, "%.2f", discarded / 1000d) + " s");
        }
        long[] measured = runOnce(fastq, project, database, store, threads);
        double seconds = measured[0] / 1000d;
        row.threads = threads;
        row.seconds += seconds;
        row.reads += measured[1];
        row.kmers += measured[2];
        System.out.println("    " + fastq.getName() + ": " + String.format(Locale.US, "%.2f", seconds)
                + " s, " + String.format(Locale.US, "%.0f", measured[1] / seconds) + " reads / s");
    }

    /**
     * Classifies the file once.
     *
     * @return the milliseconds it took, the reads and the k-mers it saw
     */
    private static long[] runOnce(File fastq, GSProject project, Database database,
            KMerStore<SmallTaxIdNode> store, int threads) throws Exception {
        // The marks of a previous pass would make the next one count fewer distinct k-mers, and they
        // live in the store, which both passes share.
        if (store.isMarkVisited()) {
            store.clearVisitedMarks();
        }
        SmallTaxTree taxTree = database.getTaxTree();
        ExecutionContext bundle = new DefaultExecutionContext(null, threads,
                project.longConfigValue(GSConfigKey.LOG_PROGRESS_UPDATE_CYCLE));
        FastqKMerMatcher matcher = newMatcher(store, taxTree, project, bundle, database);
        try {
            long start = System.currentTimeMillis();
            MatchingResult result = matcher.runMatcher(streamOf(fastq), null, null);
            long millis = System.currentTimeMillis() - start;
            CountsPerTaxid stats = result.getGlobalStats();
            return new long[] { millis, stats.getReads(), stats.getKMers() };
        } finally {
            matcher.dump();
            bundle.dump();
        }
    }

    /**
     * Creates a matcher configured as the project's match goal would configure it.
     */
    private static FastqKMerMatcher newMatcher(KMerStore<SmallTaxIdNode> store, SmallTaxTree taxTree,
            final GSProject project, ExecutionContext bundle, Database database) {
        return new FastqKMerMatcher(store, project.intConfigValue(GSConfigKey.INITIAL_READ_SIZE_BYTES),
                project.intConfigValue(GSConfigKey.THREAD_QUEUE_SIZE), bundle,
                project.booleanConfigValue(GSConfigKey.WITH_PROBS), taxTree,
                project.intConfigValue(GSConfigKey.MAX_CLASSIFICATION_PATHS),
                project.doubleConfigValue(GSConfigKey.MAX_READ_TAX_ERROR_COUNT),
                project.doubleConfigValue(GSConfigKey.MAX_READ_CLASS_ERROR_COUNT),
                project.booleanConfigValue(GSConfigKey.WRITE_ALL),
                project.intConfigValue(GSConfigKey.MIN_KMERS_FOR_CLASS),
                database.getConfigInfo().getProperty(GSProject.DB_MD5)) {
            @Override
            protected boolean isOwnUniqueKMerBits() {
                return project.booleanConfigValue(GSConfigKey.PARALLEL_DB_MATCHING);
            }

            @Override
            protected boolean isProgressBar() {
                return false;
            }
        };
    }

    private static StreamingResourceStream streamOf(File fastq) {
        return new StreamingResourceListStream(new StreamingFileResource(fastq));
    }

    private static File resolveFastq(String db, String name) {
        File file = new File(name);
        return file.isAbsolute() ? file : new File(BASE_DIR, "projects/" + db + "/fastq/" + name);
    }

    /**
     * Returns the used heap after asking the collector to run, which is as close to the live set as a
     * JVM lets one get without a heap dump. Repeated because one collection need not finish the job.
     *
     * @return the used heap in bytes
     */
    private static long usedHeapAfterGc() {
        for (int i = 0; i < 3; i++) {
            System.gc();
            try {
                Thread.sleep(200);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                break;
            }
        }
        Runtime runtime = Runtime.getRuntime();
        return runtime.totalMemory() - runtime.freeMemory();
    }

    private static double toMB(long bytes) {
        return bytes / 1000d / 1000d;
    }

    /**
     * Writes one row per store, which is what the paper's appendix reads. The file is named after the
     * database, so that a run against another one does not overwrite it and the appendix can hold the
     * two side by side.
     */
    private static void write(List<Row> rows, String db) throws Exception {
        if (!RESULTS_DIR.exists() && !RESULTS_DIR.mkdirs()) {
            throw new IllegalStateException("Could not create " + RESULTS_DIR);
        }
        File out = new File(RESULTS_DIR, "storebench_" + db + ".csv");
        PrintWriter writer = new PrintWriter(out, "UTF-8");
        try {
            // One row per measured figure, with the sorted array first: it is what the radix store is
            // held against, and the paper's table gives the radix column as a percentage of it. The
            // measurement order is the other way round and stays that way, since the first
            // classification of a run pays the JIT warm-up and that should not fall on one store only.
            Row sortedArray = rowOf(rows, STORE_SORTED_ARRAY);
            Row radix = rowOf(rows, STORE_DB);
            writer.println("project;parameter;sorted array;radix store;");
            // The first three follow from the entry count and the bytes per entry, the rest is measured.
            writer.println(line(db, "Arrays (MB)", toMB(sortedArray.arrays), toMB(radix.arrays), 0));
            writer.println(line(db, "Bloom filter (MB)", toMB(sortedArray.filter), toMB(radix.filter), 0));
            writer.println(line(db, "Store (MB)", toMB(sortedArray.arrays + sortedArray.filter),
                    toMB(radix.arrays + radix.filter), 0));
            writer.println(line(db, "Heap (MB)", toMB(sortedArray.heap), toMB(radix.heap), 0));
            writer.println(line(db, "Wall time (s)", sortedArray.seconds, radix.seconds, 2));
            writer.println(line(db, "Speed (reads / s)", sortedArray.readsPerSecond(),
                    radix.readsPerSecond(), 0));
        } finally {
            writer.close();
        }
        System.out.println();
        System.out.println("Wrote " + out);
    }

    /** Returns the row of the given store, which must be there. */
    private static Row rowOf(List<Row> rows, String store) {
        for (Row row : rows) {
            if (store.equals(row.store)) {
                return row;
            }
        }
        throw new IllegalStateException("No row for store '" + store + "'.");
    }

    /** One CSV line: the figure, the sorted array's value and the radix store's. */
    private static String line(String project, String parameter, double sortedArray, double radix, int decimals) {
        String format = "%." + decimals + "f";
        return project + ';' + parameter + ';' + String.format(Locale.US, format, sortedArray) + ';'
                + String.format(Locale.US, format, radix) + ';';
    }

    /** One store's figures. */
    private static final class Row {
        private final String project;
        private final String store;
        private final String clazz;
        private final long entries;
        private final long capacity;
        private final long perEntry;
        private final long arrays;
        private final long filter;
        private final long heap;
        private int threads;
        private double seconds;
        private long reads;
        private long kmers;

        private Row(String project, String store, KMerStore<SmallTaxIdNode> s, long heap) {
            this.project = project;
            this.store = store;
            this.clazz = s.getClass().getSimpleName();
            this.entries = s.getEntries();
            this.capacity = s.getSize();
            // A radix bucket holds one 64-bit word per k-mer. The sorted array needs a second array
            // for the taxon index, a short per k-mer, hence ten bytes instead of eight.
            this.perEntry = s instanceof KMerSortedArray ? 10 : 8;
            this.arrays = capacity * perEntry;
            this.filter = s.getFilter() == null ? 0 : s.getFilter().getBitSize() / 8;
            this.heap = heap;
        }

        private double readsPerSecond() {
            return seconds > 0 ? reads / seconds : 0;
        }

        private void print() {
            System.out.println("Store '" + store + "': " + clazz + ", " + entries + " entries, capacity "
                    + capacity);
            System.out.println("    " + perEntry + " bytes per entry: " + String.format(Locale.US, "%.0f", toMB(arrays))
                    + " MB of store, " + String.format(Locale.US, "%.0f", toMB(filter))
                    + " MB of Bloom filter, " + String.format(Locale.US, "%.0f", toMB(arrays + filter))
                    + " MB together");
            System.out.println("    heap it took: " + String.format(Locale.US, "%.0f", toMB(heap)) + " MB");
        }

    }
}
