package org.metagene.ftdbexp.stquality;

import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.FinerTreeMaker;
import org.metagene.genestrip.goals.MatchResultGoal;
import org.metagene.genestrip.match.FastqKMerMatcher;
import org.metagene.genestrip.match.CountsPerTaxid;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayDeque;
import java.util.Deque;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

/**
 * Scores a refined database on isolates it has never seen, at the sequence-type level.
 * <p>
 * This is the held-out counterpart of {@link STQualityCountsGoal}. That goal re-reads the genomes
 * the database was filled from, so every k-mer it looks for is there by construction and its recall
 * is one whatever the refinement did; it measures placement, not diagnosis. Here the reads come from
 * isolates outside the database, and the question is the clinical one: does the classification name
 * the lineage the isolate belongs to?
 * <p>
 * The rule is the one Genestrip already classifies a read with, applied to a whole isolate: see
 * {@link PathVoteSTModel}. It places the isolate at a node of the taxonomy, and the answer is the
 * types that node still leaves in question -- one where the node holds a single lineage, several
 * where it holds several. An answer of several is not a failure to answer but a weaker answer, and
 * it is scored as one: {@code 1/n} where the isolate's own type is among the {@code n}, nothing
 * where it is not, which is the paper's {@code p_st} applied to an isolate rather than to a k-mer.
 * Placing an isolate at the species therefore earns almost nothing instead of counting as right.
 * <p>
 * The same rule is scored in three currencies -- one vote per read, per matched k-mer, and per
 * distinct matched k-mer -- and against the majority call as a control. {@link NaiveBayesSTModel} is
 * the weighted alternative and can be selected instead; it is not the default, since it answers a
 * different question and lost the comparison it was built for.
 * <p>
 * Invoked directly:
 * <pre>
 *   java ... IsolateSTAccuracyMain cdiff data/fastq/cdiff_isolates.txt results/cdiff_isolate_st.csv [ftdb|db]
 * </pre>
 * Run it against both database variants: what the refinement is worth is the difference, and the
 * unrefined run is expected to make almost no call at all, since without refined nodes there is
 * hardly a node below the species that leaves a single type in question.
 */
public class IsolateSTAccuracyMain {
    private static final File BASE_DIR = new File("./data");
    /** The rule Genestrip itself classifies a read with, applied to the isolate; the default. */
    private static final String PATH = "path";
    /** The weighted alternative. Kept because it is a different question, not because it is better. */
    private static final String BAYES = "bayes";
    /** The currencies a rule is scored in, in the order their columns are written. */
    private static final String[] CURRENCIES = { "read", "kmer", "unique" };

    /**
     * Returns one isolate's counts in the currency of the given column.
     *
     * @param call  the isolate's tally
     * @param index the index into {@link #CURRENCIES}
     * @return the counts by dense node position
     */
    private static Map<Integer, Long> currency(IsolateSTCall call, int index) {
        switch (index) {
            case 0: return call.getReadsPerNode();
            case 1: return call.getKMersPerNode();
            default: return call.getUniqueKMersPerNode();
        }
    }

    /**
     * Creates the requested rule.
     *
     * @param name        {@link #PATH} or {@link #BAYES}
     * @param composition the counted database
     * @return the classifier
     * @throws IllegalArgumentException if the name is neither
     */
    private static STClassifier createClassifier(String name, STComposition composition) {
        if (PATH.equals(name)) {
            return new PathVoteSTModel(composition);
        }
        if (BAYES.equals(name)) {
            return new NaiveBayesSTModel(composition);
        }
        throw new IllegalArgumentException("Not a classifier: '" + name + "'. Expected '" + PATH
                + "' or '" + BAYES + "'.");
    }

    /**
     * Runs the experiment.
     *
     * @param args {@code <db> <fastq map> <out csv> [ftdb|db]}
     * @throws Exception if the database, the typing or the reads cannot be read
     */
    public static void main(String[] args) throws Exception {
        if (args.length < 3) {
            System.err.println("Usage: IsolateSTAccuracyMain <db> <fastq map> <out csv> [ftdb|db]");
            System.err.println("  <db>          name of the database project under data/projects");
            System.err.println("  <fastq map>   the map of isolate key to fastq file, e.g.");
            System.err.println("                data/fastq/cdiff_isolates.txt");
            System.err.println("  <out csv>     where the per-isolate result is written");
            System.err.println("  ftdb|db       which database variant to classify against; default 'ftdb'");
            System.err.println("  path|bayes    which rule to classify with; default 'path', the rule");
            System.err.println("                Genestrip itself uses on a read, over the isolate's");
            System.err.println("                pooled counts. 'bayes' is the weighted alternative.");
            System.err.println();
            System.err.println("The isolates' sequence types are read from <fastq map>.st.csv, i.e. the");
            System.err.println("map's path with '.st.csv' appended - see IsolateSTTruth for the format.");
            System.exit(1);
        }
        String db = args[0];
        File fqMap = new File(args[1]);
        File outCsv = new File(args[2]);
        boolean refined = args.length <= 3 || args[3].isEmpty() || "ftdb".equals(args[3]);
        String classifier = args.length <= 4 || args[4].isEmpty() ? PATH : args[4].trim();

        if (!fqMap.exists()) {
            System.err.println("No such fastq map: " + fqMap.getAbsolutePath());
            System.exit(1);
        }
        File truthFile = new File(fqMap.getPath() + ".st.csv");
        if (!truthFile.exists()) {
            System.err.println("No such typing: " + truthFile.getAbsolutePath());
            System.err.println("It maps each key of " + fqMap.getName() + " to the isolate's sequence type.");
            System.err.println("For the C. difficile isolates the types are in Table 2 of the source study,");
            System.err.println("keyed by the same identifiers the archive carries as sample_alias.");
            System.exit(1);
        }
        IsolateSTTruth truth = new IsolateSTTruth(truthFile);
        System.out.println("Read " + truth.size() + " typed isolate(s) from " + truthFile);
        if (truth.size() == 0) {
            System.err.println("Nothing to score against.");
            System.exit(1);
        }

        Classified classified = classify(db, fqMap, refined, truth, classifier);
        write(classified.calls, classified.model, outCsv, db, refined);
    }

    /**
     * Classifies every isolate of the map and returns one tally per fastq key.
     *
     * @param db      the database project
     * @param fqMap   the fastq map
     * @param refined whether to classify against the refined database
     * @param truth   the isolates' sequence types
     * @return the tallies keyed by fastq key, and the model to score them with
     * @throws IOException if the database or the reads cannot be read
     */
    private static Classified classify(String db, File fqMap, boolean refined,
                                       IsolateSTTruth truth, String classifier) throws IOException {
        FTProject project = new FTProject(new GSCommon(BASE_DIR), db, null, null, fqMap.getPath(),
                null, null, null, null, null, null, false);
        project.initConfigParam(GSConfigKey.THREADS, -1);
        FinerTreeMaker<FTProject> maker = new FinerTreeMaker<>(project);

        GoalKey matchKey = refined ? FTGoalKey.FTMATCHRES : org.metagene.genestrip.GSGoalKey.MATCHRES;
        GoalKey loadKey = refined ? FTGoalKey.LOAD_FTDB : org.metagene.genestrip.GSGoalKey.LOAD_DB;

        @SuppressWarnings("unchecked")
        ObjectGoal<Database, FTProject> dbGoal = (ObjectGoal<Database, FTProject>) maker.getGoal(loadKey);
        SmallTaxTree tree = dbGoal.get().getTaxTree();
        File mlst = new File("results/" + db + "_mlst.csv");
        if (!mlst.exists()) {
            throw new IOException("No such typing of the database's genomes: " + mlst.getAbsolutePath()
                    + ". Produce it with bin/" + db + "_eval.sh; without it no node predicts anything.");
        }
        // Counted once, whichever rule is scored against it, so that two result files differ by the
        // rule and never by two readings of the same database.
        STComposition composition = new STComposition(tree, new STGroundTruth(mlst));
        STClassifier model = createClassifier(classifier, composition);
        System.out.println("Classifying with '" + model.getName() + "'.");
        System.out.println("Model over " + composition.getTypes().size() + " sequence type(s) and "
                + composition.getCollectionTotal() + " genome(s).");

        Map<String, IsolateSTCall> result = new LinkedHashMap<>();
        // One tally per matcher thread rather than one behind a lock: the callback runs for every
        // read, and guarding it globally serialises the whole run - AccuracyEvaluator carries the
        // same construction and the same reason. The merge in afterKey() is exact because every
        // field is a count, and it sees the worker threads' writes for the reason set out there.
        final Queue<IsolateSTCall> threadTallies = new ConcurrentLinkedQueue<>();
        final ThreadLocal<IsolateSTCall> local = new ThreadLocal<IsolateSTCall>() {
            @Override
            protected IsolateSTCall initialValue() {
                // The key is filled in at afterKey(), where it is known; a per-thread tally is only
                // ever a set of counters to be merged and is never read by key.
                IsolateSTCall fresh = new IsolateSTCall(null, null);
                threadTallies.add(fresh);
                return fresh;
            }
        };

        MatchResultGoal<?> matchGoal = (MatchResultGoal<?>) maker.getGoal(matchKey);
        matchGoal.setAfterMatchCallback(new MatchResultGoal.AfterMatchCallback() {
            @Override
            public void afterMatch(FastqKMerMatcher.MatcherReadEntry entry, boolean found) {
                IsolateSTCall tally = local.get();
                SmallTaxTree.SmallTaxIdNode node = entry.classNode;
                if (node == null) {
                    tally.recordUnclassified();
                } else {
                    int pos = node.getPosition();
                    tally.recordClassified(pos, composition.getMajorityST(pos), composition.getCandidates(pos));
                }
            }

            @Override
            public void afterKey(String key, MatchingResult res) {
                IsolateSTCall merged = new IsolateSTCall(key, truth.getST(key));
                // The k-mer level, taken from the matcher's own per-taxon tally rather than from the
                // reads: both counts hold every k-mer the database keeps at that node and the isolate
                // carries, whether or not the read carrying it was classified there. `kmers' counts
                // them with multiplicity and `unique kmers' once each; see IsolateSTCall for what
                // each is robust against.
                for (Map.Entry<String, CountsPerTaxid> e : res.getTaxid2Stats().entrySet()) {
                    SmallTaxTree.SmallTaxIdNode node = tree.getNodeByTaxId(e.getKey());
                    if (node == null) {
                        continue;
                    }
                    int pos = node.getPosition();
                    merged.recordKMers(pos, e.getValue().getKMers(), e.getValue().getUniqueKMers(),
                            composition.getMajorityST(pos));
                }
                for (IsolateSTCall threadTally : threadTallies) {
                    merged.add(threadTally);
                    // Reset, not discard: the thread still holds this instance in its ThreadLocal
                    // and would keep writing into an object no longer in the queue, so the next
                    // file's votes from that thread would go missing.
                    threadTally.reset();
                }
                // merge() rather than put(): a paired isolate arrives as two files under one key,
                // and the second must add to the first instead of replacing it.
                IsolateSTCall previous = result.put(key, merged);
                if (previous != null) {
                    merged.add(previous);
                }
            }
        });
        try {
            matchGoal.make();
        } finally {
            // Shuts the shared thread pool down, which is what lets the virtual machine exit. Its
            // threads are created with plain `new Thread', so they are not daemons and keep the
            // process alive on their own once main() has returned: the run then prints its results,
            // reports success and hangs, and under Maven it is the forked JVM that never ends. Only
            // dumpAll() reaches ExecutionContext.dump() and hence executorService.shutdown().
            maker.dumpAll();
        }
        return new Classified(result, model);
    }

    /** The tallies and the model they are to be scored with. */
    private static final class Classified {
        private final Map<String, IsolateSTCall> calls;
        private final STClassifier model;

        Classified(Map<String, IsolateSTCall> calls, STClassifier model) {
            this.calls = calls;
            this.model = model;
        }
    }

    /**
     * Writes one row per isolate and a summary to standard output.
     *
     * @param calls   the tallies
     * @param model   the classifier the tallies are scored with
     * @param outCsv  the file to write
     * @param db      the database project
     * @param refined whether the refined database was classified against
     * @throws IOException if the file cannot be written
     */
    private static void write(Map<String, IsolateSTCall> calls, STClassifier model, File outCsv,
                              String db, boolean refined) throws IOException {
        File parent = outCsv.getParentFile();
        if (parent != null && !parent.isDirectory() && !parent.mkdirs()) {
            throw new IOException("Cannot create " + parent.getAbsolutePath());
        }
        int typed = 0;
        int untyped = 0;
        int majorityCorrect = 0;
        long kMersCorrect = 0;
        long kMersMatched = 0;
        long readsCorrect = 0;
        long readsClassified = 0;
        // One accumulator per currency, in the order the columns are written.
        double[] precSum = new double[CURRENCIES.length];
        double[] gainSum = new double[CURRENCIES.length];
        double[] sizeSum = new double[CURRENCIES.length];
        int[] hits = new int[CURRENCIES.length];
        int[] exact = new int[CURRENCIES.length];
        try (PrintStream out = new PrintStream(outCsv)) {
            // The rule's own name in the header, so that two result files scored under different
            // rules cannot be mistaken for one another.
            String n = model.getName();
            StringBuilder header = new StringBuilder("db;variant;isolate;true st;");
            for (String currency : CURRENCIES) {
                header.append(n).append(' ').append(currency).append(" st;")
                        .append(n).append(' ').append(currency).append(" cand;")
                        .append(n).append(' ').append(currency).append(" hit;")
                        .append(n).append(' ').append(currency).append(" prec;")
                        .append(n).append(' ').append(currency).append(" gain;");
            }
            header.append("majority st;majority verdict;majority share;")
                    .append("read accuracy;kmer accuracy;unique kmer accuracy;mean candidates;reads;")
                    .append("classified;correct reads;matched kmers;matched unique kmers;");
            out.println(header);
            for (IsolateSTCall call : calls.values()) {
                String trueST = call.getTrueST();
                if (trueST == null) {
                    untyped++;
                } else {
                    typed++;
                    kMersCorrect += Math.round(call.getKMerAccuracy() * call.getMatchedKMers());
                    kMersMatched += call.getMatchedKMers();
                    readsCorrect += call.getCorrectReads();
                    readsClassified += call.getClassified();
                }
                IsolateSTCall.Verdict majority = call.getMajorityVerdict();
                if (majority == IsolateSTCall.Verdict.CORRECT) {
                    majorityCorrect++;
                }
                out.print(db);
                out.print(';');
                out.print(refined ? "refined" : "unrefined");
                out.print(';');
                out.print(call.getKey());
                out.print(';');
                out.print(trueST == null ? "" : trueST);
                out.print(';');
                for (int i = 0; i < CURRENCIES.length; i++) {
                    // The same rule in every currency: one vote per read, per matched k-mer, and per
                    // distinct matched k-mer.
                    Map<Integer, Long> counts = currency(call, i);
                    String st = model.classify(counts);
                    int cand = model.getCandidates(counts);
                    boolean hit = trueST != null && model.isHit(counts, trueST);
                    double prec = trueST == null ? 0 : model.getPrecision(counts, trueST);
                    double gain = trueST == null ? 0 : model.getInformationGain(counts, trueST);
                    if (trueST != null) {
                        precSum[i] += prec;
                        gainSum[i] += gain;
                        sizeSum[i] += cand;
                        if (hit) {
                            hits[i]++;
                            if (cand == 1) {
                                exact[i]++;
                            }
                        }
                    }
                    out.print(st == null ? "" : st);
                    out.print(';');
                    out.print(cand);
                    out.print(';');
                    out.print(hit);
                    out.print(';');
                    out.printf("%.6f", prec);
                    out.print(';');
                    out.printf("%.6f", gain);
                    out.print(';');
                }
                out.print(call.getCalledST() == null ? "" : call.getCalledST());
                out.print(';');
                out.print(majority);
                out.print(';');
                out.printf("%.8f", call.getCalledShare());
                out.print(';');
                out.printf("%.8f", call.getReadAccuracy());
                out.print(';');
                out.printf("%.8f", call.getKMerAccuracy());
                out.print(';');
                out.printf("%.8f", call.getUniqueKMerAccuracy());
                out.print(';');
                out.printf("%.4f", call.getMeanCandidates());
                out.print(';');
                out.print(call.getReads());
                out.print(';');
                out.print(call.getClassified());
                out.print(';');
                out.print(call.getCorrectReads());
                out.print(';');
                out.print(call.getMatchedKMers());
                out.print(';');
                out.print(call.getMatchedUniqueKMers());
                out.print(';');
                out.println();
            }
        }
        System.out.println("Wrote " + outCsv.getAbsolutePath());
        System.out.println();
        System.out.println("  isolates typed          : " + typed
                + (untyped > 0 ? " (" + untyped + " untyped, not scored)" : ""));
        System.out.printf("  isolate accuracy (major): %.4f%n",
                typed == 0 ? Double.NaN : ((double) majorityCorrect) / typed);
        System.out.println();
        System.out.printf("  %-8s %8s %8s %8s %8s %8s%n",
                "currency", "gain", "prec", "hit", "exact", "answer");
        for (int i = 0; i < CURRENCIES.length; i++) {
            System.out.printf("  %-8s %8.4f %8.4f %8.4f %8.4f %8.2f%n", CURRENCIES[i],
                    typed == 0 ? Double.NaN : gainSum[i] / typed,
                    typed == 0 ? Double.NaN : precSum[i] / typed,
                    typed == 0 ? Double.NaN : ((double) hits[i]) / typed,
                    typed == 0 ? Double.NaN : ((double) exact[i]) / typed,
                    typed == 0 ? Double.NaN : sizeSum[i] / typed);
        }
        System.out.println();
        System.out.println("  gain    how much likelier the isolate's own type became, over how much");
        System.out.println("          likelier it would have to become to be named outright. Zero for an");
        System.out.println("          answer that narrowed nothing, one for an answer that pinned the");
        System.out.println("          lineage, and it does not read a node of 35 genomes of one lineage");
        System.out.println("          and one of another as a coin toss the way 1/n does.");
        System.out.println("  prec    the paper's p_st on an isolate: 1/n where the answer leaves n types");
        System.out.println("          in question and the true one is among them, else nothing.");
        System.out.println("  hit     the share of isolates whose own type is in the answer at all, and");
        System.out.println("  exact   the share where it is the only one left. An average alone cannot");
        System.out.println("          tell a vague answer from a wrong one; these two can.");
        System.out.println("  answer  the mean number of types left in question.");
        System.out.println();
        System.out.printf("  read accuracy (pooled)  : %.4f  over %,d classified read(s)%n",
                readsClassified == 0 ? Double.NaN : ((double) readsCorrect) / readsClassified,
                readsClassified);
        System.out.printf("  k-mer accuracy (pooled) : %.4f  over %,d matched k-mer(s)%n",
                kMersMatched == 0 ? Double.NaN : ((double) kMersCorrect) / kMersMatched, kMersMatched);
        System.out.println();
        System.out.println("The k-mer level sits ON TOP of the read level and is not a substitute for it:"
                + " a read is classified to the lowest node its k-mers agree on, so a single specific"
                + " k-mer is outvoted by its own read and never reaches the read-level tally.");
        System.out.println();
        System.out.println("The majority call counts every read alike and is therefore carried by whichever"
                + " type the collection has most of; it is the control the rule has to beat. Run both"
                + " database variants: what the refinement is worth is the difference between them.");
    }
}
