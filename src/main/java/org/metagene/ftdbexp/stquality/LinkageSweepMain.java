package org.metagene.ftdbexp.stquality;

import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.finertree.FTConfigKey;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.cluster.SimpleAggloClustering;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Rebuilds the refined database under each of several linkage methods and files the resulting
 * sequence-type quality under its own name, so that the linkage question is settled by comparing the
 * files rather than by keeping three runs apart by hand.
 * <p>
 * The linkages are given as one comma-separated argument and default to the three the paper reports:
 * <pre>
 *   java ... org.metagene.ftdbexp.stquality.LinkageSweepMain cdiff results/cdiff_mlst.csv
 *   java ... org.metagene.ftdbexp.stquality.LinkageSweepMain cdiff results/cdiff_mlst.csv UPGMA
 * </pre>
 * The result lands in {@code results/} as {@code <db>_ftstquality_<linkage>.csv}, the linkage
 * appended so the files sort together and the unsuffixed name stays free for whatever the project's
 * own configuration produces. The unrefined measure does not depend on the linkage; take it once
 * with {@code STQualityMain}.
 * <p>
 * All linkages run in ONE process, which is what makes the sweep worth having: the k-mer index and
 * the pairwise intersection counts are built from the genomes once and every further linkage only
 * re-cuts the dendrogram. What that requires is that {@link #reset} drops the previous dendrogram
 * without poisoning the goals that read -- see the comment there, which records how it was got wrong.
 * A linkage that fails does not stop the others; the exit status says whether any did.
 * <p>
 * WHAT IS RESET, and nothing else: the goals whose output encodes the dendrogram. Never
 * {@code storekmerindex}, {@code kmerindexbloom} or the unrefined database. Doing that by hand means
 * deleting exactly the right files and no others, and the make targets are no help --
 * {@code clean} reaches only the internal goal a request is aggregated into, while {@code cleanall}
 * descends into every dependency that permits it and takes the unrefined database with it.
 */
public class LinkageSweepMain {
    private static final File BASE_DIR = new File("./data");
    private static final File RESULTS_DIR = new File("./results");

    /**
     * The goals whose output encodes the dendrogram, and hence the ones that have to go before the
     * next linkage is made. Ordered from the top of the graph down, so that a file goal is gone
     * before anything that could consider it up to date is asked.
     * <p>
     * {@code FT_ST_QUALITY} and its counting goal are here because a report describing the previous
     * dendrogram is worse than none. {@code FTDB}, {@code FTDBINFO} and {@code FT_SVG_TAX_TREE} are
     * the files that hold the refined tree. {@code LOAD_FTDB}, {@code UPDATE_STORE_GOAL} and
     * {@code DENDROGRAM} are in-memory values, dropped so that a run started against an existing
     * refined database rebuilds it rather than measuring the previous linkage under this one's name.
     */
    private static final GoalKey[] DENDROGRAM_DEPENDENT = {
            STQualityMaker.FT_ST_QUALITY,
            STQualityMaker.FT_ST_QUALITY_COUNTS,
            FTGoalKey.FTDBINFO,
            FTGoalKey.FT_SVG_TAX_TREE,
            FTGoalKey.LOAD_FTDB,
            FTGoalKey.FTDB,
            FTGoalKey.UPDATE_STORE_GOAL,
            FTGoalKey.DENDROGRAM
    };

    /** The linkages the sweep takes when the call names none. */
    private static final SimpleAggloClustering.Method[] DEFAULT_METHODS = {
            SimpleAggloClustering.Method.SINGLE_LINKAGE,
            SimpleAggloClustering.Method.UPGMA,
            SimpleAggloClustering.Method.COMPLETE_LINKAGE
    };

    /**
     * Rebuilds and measures under each requested linkage.
     *
     * @param args {@code <db> <sequence type csv> [<linkage>[,<linkage>...]]}
     * @throws Exception if the database or the typing cannot be read
     */
    public static void main(String[] args) throws Exception {
        if (args.length < 2) {
            System.err.println("Usage: LinkageSweepMain <db> <sequence type csv> [<linkage>[,<linkage>...]]");
            System.err.println("  <db>                  name of the database project under data/projects");
            System.err.println("  <sequence type csv>   the CSV written by bin/mlst_assemblies.sh");
            System.err.println("  <linkage>             SINGLE_LINKAGE, COMPLETE_LINKAGE, UPGMA or WPGMA;");
            System.err.println("                        several comma separated, all in one process.");
            System.err.println("                        Left out: SINGLE_LINKAGE, UPGMA, COMPLETE_LINKAGE.");
            System.exit(1);
        }
        String db = args[0];
        File stCsv = new File(args[1]);
        if (!stCsv.exists()) {
            System.err.println("No such sequence type file: " + stCsv.getAbsolutePath());
            System.err.println("Produce it with bin/cdiff_eval.sh.");
            System.exit(1);
        }
        List<SimpleAggloClustering.Method> methods = parseMethods(args);

        STGroundTruth groundTruth = new STGroundTruth(stCsv);
        System.out.println("Read " + groundTruth.getTypedLeaves() + " typed genome(s) in "
                + groundTruth.getSTCount() + " distinct sequence type(s) from " + stCsv
                + "; " + groundTruth.getUntypedLeaves() + " listed without a type.");
        if (groundTruth.getTypedLeaves() == 0) {
            System.err.println("Nothing to measure against.");
            System.exit(1);
        }
        if (!RESULTS_DIR.isDirectory() && !RESULTS_DIR.mkdirs()) {
            System.err.println("Cannot create " + RESULTS_DIR.getAbsolutePath());
            System.exit(1);
        }

        FTProject project = new FTProject(new GSCommon(BASE_DIR), db, null, null, null,
                null, null, null, null, null, null, false);
        project.initConfigParam(GSConfigKey.THREADS, -1);

        STQualityMaker maker = new STQualityMaker(project, groundTruth);
        List<String> written = new ArrayList<>();
        List<SimpleAggloClustering.Method> failed = new ArrayList<>();
        try {
            for (SimpleAggloClustering.Method method : methods) {
                System.out.println();
                System.out.println("=== refined database, clusterMethod=" + method + " ===");
                try {
                    if (!project.initConfigParam(FTConfigKey.CLUSTER_METHOD, method)) {
                        throw new IllegalStateException("The configuration refused the linkage " + method
                                + ". Without it taking effect the run would cluster the way the project's"
                                + " own configuration says and be filed under a name it does not have.");
                    }
                    reset(maker);
                    reportIndexReuse(maker);

                    STQualityCSVGoal goal = (STQualityCSVGoal) maker.getCSVGoal(true);
                    goal.make();
                    written.addAll(fileInto(goal, db, method));
                } catch (Exception e) {
                    // One linkage failing is no reason to lose the others: they are independent
                    // re-cuts of the same dendrogram input, and the next pass resets everything the
                    // failed one could have left behind.
                    failed.add(method);
                    System.out.println();
                    System.err.println("Linkage " + method + " failed, going on with the rest:");
                    e.printStackTrace();
                }
            }
        } finally {
            maker.dumpAll();
        }

        System.out.println();
        System.out.println("Wrote " + written.size() + " file(s) to " + RESULTS_DIR.getAbsolutePath() + ":");
        for (String name : written) {
            System.out.println("  " + name);
        }
        System.out.println();
        System.out.println("The unrefined baseline is linkage-independent and is NOT taken here; produce"
                + " it once with STQualityMain. p_st can only be judged as a change from what the"
                + " taxonomy already gave, so the comparison needs it.");
        if (!failed.isEmpty()) {
            System.err.println();
            System.err.println("FAILED: " + failed + " - the file(s) above do not cover " + failed.size()
                    + " of the " + methods.size() + " requested linkage(s).");
            System.exit(1);
        }
    }

    /**
     * Reads the requested linkages off the command line, accepting them either as separate arguments
     * or comma separated within one, since the Maven property that carries them is a single string.
     *
     * @param args the command line, the linkages starting at index 2
     * @return the linkages to sweep, {@link #DEFAULT_METHODS} if none were named
     */
    private static List<SimpleAggloClustering.Method> parseMethods(String[] args) {
        List<SimpleAggloClustering.Method> methods = new ArrayList<>();
        for (int i = 2; i < args.length; i++) {
            for (String name : args[i].split(",")) {
                name = name.trim();
                if (name.isEmpty()) {
                    continue;
                }
                try {
                    SimpleAggloClustering.Method method =
                            SimpleAggloClustering.Method.valueOf(name.toUpperCase(Locale.ROOT));
                    if (!methods.contains(method)) {
                        methods.add(method);
                    }
                } catch (IllegalArgumentException e) {
                    System.err.println("Not a linkage method: " + name);
                    System.err.println("Expected one of: SINGLE_LINKAGE, COMPLETE_LINKAGE, UPGMA, WPGMA");
                    System.exit(1);
                }
            }
        }
        if (methods.isEmpty()) {
            for (SimpleAggloClustering.Method method : DEFAULT_METHODS) {
                methods.add(method);
            }
        }
        System.out.println("Sweeping " + methods.size() + " linkage(s) in this process: " + methods);
        return methods;
    }

    /**
     * Drops everything that encodes the previous dendrogram, and nothing else. A file goal has its
     * own output removed and an object goal its value, both without touching dependencies - which is
     * the whole point, since the k-mer index is a dependency and rebuilding it would cost the pass
     * over the genomes this sweep exists to avoid.
     *
     * @param maker the maker holding the goals
     */
    private static void reset(STQualityMaker maker) {
        for (GoalKey key : DENDROGRAM_DEPENDENT) {
            Goal<FTProject> goal = maker.getGoal(key);
            if (goal == null) {
                // Not a reason to stop, but not something to pass over in silence either: a goal that
                // has been renamed leaves its predecessor's output in place, and the next pass would
                // then measure the previous linkage under this one's name.
                System.out.println("  WARNING: no goal for '" + key.getName() + "' - if it exists under"
                        + " another key now, its output survives this reset and the result is not this"
                        + " linkage's.");
                continue;
            }
            // cleanThis(), not dump(): for an object goal both discard the value, but
            // FastaReaderGoal.dump() additionally ends the consumer threads of a pass that is not
            // running, which this reset has no reason to do.
            //
            // Neither of them, though, is what made the second linkage read nothing. A pass ends by
            // setting the flag its consumers watch -- every reader goal does that itself once it has
            // its result -- and nothing put the flag back, so a goal that had read once could not read
            // again whatever it was reset with. That is fixed where it belongs, in
            // FastaReaderGoal.readFastas(), which now readies itself; see the comment there.
            goal.cleanThis();
            System.out.println("  reset " + key.getName());
        }
    }

    /**
     * Says whether the k-mer index survived the previous pass, which is the one thing that decides
     * whether this sweep is cheaper than three separate runs.
     *
     * @param maker the maker holding the goals
     */
    private static void reportIndexReuse(STQualityMaker maker) {
        report(maker, FTGoalKey.STORE_KMER_INDEX, "the k-mer index on disk");
        report(maker, FTGoalKey.INTERSECT_COUNT, "the pairwise intersection counts");
    }

    private static void report(STQualityMaker maker, GoalKey key, String what) {
        Goal<FTProject> goal = maker.getGoal(key);
        if (goal != null) {
            System.out.println("  " + (goal.isMade() ? "reusing " : "rebuilding ") + what
                    + " (" + key.getName() + ")");
        }
    }

    /**
     * Copies a goal's output into {@code results/}, appending the linkage to the file's stem.
     *
     * @param goal    the goal whose files are to be filed
     * @param db      the database name, only used to check the file belongs to it
     * @param method  the linkage to append, or {@code null} for a result that has none
     * @return the names written
     * @throws IOException if a file cannot be copied
     */
    private static List<String> fileInto(STQualityCSVGoal goal, String db,
                                         SimpleAggloClustering.Method method) throws IOException {
        List<String> written = new ArrayList<>();
        for (File file : goal.getFiles()) {
            if (!file.exists()) {
                throw new IOException("The goal reported " + file + " but did not write it.");
            }
            String name = file.getName();
            if (method != null) {
                int dot = name.lastIndexOf('.');
                String stem = dot < 0 ? name : name.substring(0, dot);
                String ext = dot < 0 ? "" : name.substring(dot);
                name = stem + "_" + method.name().toLowerCase(Locale.ROOT) + ext;
            }
            File target = new File(RESULTS_DIR, name);
            Files.copy(file.toPath(), target.toPath(), StandardCopyOption.REPLACE_EXISTING);
            System.out.println("  " + file.getAbsolutePath() + "  ->  " + target.getPath());
            written.add(name);
        }
        return written;
    }
}
