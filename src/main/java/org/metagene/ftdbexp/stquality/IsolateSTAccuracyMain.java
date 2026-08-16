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
 * Every node of the database is a predictor: the types in question at it are the ones its genomes
 * carry, and it predicts the most frequent of them. A read therefore always yields a prediction and
 * nothing has to be thresholded -- what would otherwise be a convention about how many reads make a
 * call becomes a plain count of how many reads were predicted right. See {@link IsolateSTCall} for
 * what follows from that, and why the per-read accuracy and the per-isolate call answer different
 * questions.
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
            System.err.println();
            System.err.println("The isolates' sequence types are read from <fastq map>.st.csv, i.e. the");
            System.err.println("map's path with '.st.csv' appended - see IsolateSTTruth for the format.");
            System.exit(1);
        }
        String db = args[0];
        File fqMap = new File(args[1]);
        File outCsv = new File(args[2]);
        boolean refined = args.length <= 3 || args[3].isEmpty() || "ftdb".equals(args[3]);

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

        Classified classified = classify(db, fqMap, refined, truth);
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
                                       IsolateSTTruth truth) throws IOException {
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
        NaiveBayesSTModel model = new NaiveBayesSTModel(tree, new STGroundTruth(mlst));
        System.out.println("Model over " + model.getTypes().size() + " sequence type(s) and "
                + model.getPriorTotal() + " genome(s).");

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
                    tally.recordClassified(pos, model.getMajorityST(pos), model.getCandidates(pos));
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
                            model.getMajorityST(pos));
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
        matchGoal.make();
        return new Classified(result, model);
    }

    /** The tallies and the model they are to be scored with. */
    private static final class Classified {
        private final Map<String, IsolateSTCall> calls;
        private final NaiveBayesSTModel model;

        Classified(Map<String, IsolateSTCall> calls, NaiveBayesSTModel model) {
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
    private static void write(Map<String, IsolateSTCall> calls, NaiveBayesSTModel model, File outCsv,
                              String db, boolean refined) throws IOException {
        File parent = outCsv.getParentFile();
        if (parent != null && !parent.isDirectory() && !parent.mkdirs()) {
            throw new IOException("Cannot create " + parent.getAbsolutePath());
        }
        int correct = 0;
        int wrong = 0;
        int noCall = 0;
        int untyped = 0;
        int majorityCorrect = 0;
        int kmerCorrect = 0;
        int uniqueCorrect = 0;
        long kMersCorrect = 0;
        long kMersMatched = 0;
        long readsCorrect = 0;
        long readsClassified = 0;
        try (PrintStream out = new PrintStream(outCsv)) {
            out.println("db;variant;isolate;true st;bayes st;bayes verdict;kmer st;kmer verdict;"
                    + "unique kmer st;unique kmer verdict;majority st;majority verdict;majority share;"
                    + "read accuracy;kmer accuracy;unique kmer accuracy;mean candidates;reads;"
                    + "classified;correct reads;matched kmers;matched unique kmers;");
            for (IsolateSTCall call : calls.values()) {
                String bayes = model.classify(call.getReadsPerNode());
                String kmer = model.classify(call.getKMersPerNode());
                String uniqueKmer = model.classify(call.getUniqueKMersPerNode());
                IsolateSTCall.Verdict verdict = call.getVerdict(bayes);
                IsolateSTCall.Verdict kmerVerdict = call.getVerdict(kmer);
                IsolateSTCall.Verdict uniqueVerdict = call.getVerdict(uniqueKmer);
                IsolateSTCall.Verdict majority = call.getMajorityVerdict();
                switch (verdict) {
                    case CORRECT: correct++; break;
                    case WRONG: wrong++; break;
                    case NO_CALL: noCall++; break;
                    default: untyped++; break;
                }
                if (majority == IsolateSTCall.Verdict.CORRECT) {
                    majorityCorrect++;
                }
                if (kmerVerdict == IsolateSTCall.Verdict.CORRECT) {
                    kmerCorrect++;
                }
                if (uniqueVerdict == IsolateSTCall.Verdict.CORRECT) {
                    uniqueCorrect++;
                }
                if (kmerVerdict != IsolateSTCall.Verdict.NOT_SCORED) {
                    kMersCorrect += Math.round(call.getKMerAccuracy() * call.getMatchedKMers());
                    kMersMatched += call.getMatchedKMers();
                }
                if (verdict != IsolateSTCall.Verdict.NOT_SCORED) {
                    readsCorrect += call.getCorrectReads();
                    readsClassified += call.getClassified();
                }
                out.print(db);
                out.print(';');
                out.print(refined ? "refined" : "unrefined");
                out.print(';');
                out.print(call.getKey());
                out.print(';');
                out.print(call.getTrueST() == null ? "" : call.getTrueST());
                out.print(';');
                out.print(bayes == null ? "" : bayes);
                out.print(';');
                out.print(verdict);
                out.print(';');
                out.print(kmer == null ? "" : kmer);
                out.print(';');
                out.print(kmerVerdict);
                out.print(';');
                out.print(uniqueKmer == null ? "" : uniqueKmer);
                out.print(';');
                out.print(uniqueVerdict);
                out.print(';');
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
                out.println(';');
            }
        }
        int typed = correct + wrong + noCall;
        System.out.println("Wrote " + outCsv.getAbsolutePath());
        System.out.println();
        System.out.println("  isolates typed          : " + typed
                + (untyped > 0 ? " (" + untyped + " untyped, not scored)" : ""));
        System.out.println("  correct / wrong / none  : " + correct + " / " + wrong + " / " + noCall);
        System.out.printf("  isolate accuracy (Bayes): %.4f%n",
                typed == 0 ? Double.NaN : ((double) correct) / typed);
        System.out.printf("  isolate accuracy (k-mer): %.4f%n",
                typed == 0 ? Double.NaN : ((double) kmerCorrect) / typed);
        System.out.printf("  isolate accuracy (uniq) : %.4f%n",
                typed == 0 ? Double.NaN : ((double) uniqueCorrect) / typed);
        System.out.printf("  isolate accuracy (major): %.4f%n",
                typed == 0 ? Double.NaN : ((double) majorityCorrect) / typed);
        System.out.printf("  read accuracy (pooled)  : %.4f  over %,d classified read(s)%n",
                readsClassified == 0 ? Double.NaN : ((double) readsCorrect) / readsClassified,
                readsClassified);
        System.out.printf("  k-mer accuracy (pooled) : %.4f  over %,d matched k-mer(s)%n",
                kMersMatched == 0 ? Double.NaN : ((double) kMersCorrect) / kMersMatched, kMersMatched);
        System.out.println();
        System.out.println("The k-mer level sits ON TOP of the read level and is not a substitute for it:"
                + " a read is classified to the lowest node its k-mers agree on, so a single specific"
                + " k-mer is outvoted by its own read and never reaches the read-level classifier."
                + " Where the k-mer accuracy exceeds the read accuracy, that is the evidence the read"
                + " classification is discarding.");
        System.out.println();
        System.out.println("The read accuracy measures the database, since it is the nodes that predict."
                + " The two isolate accuracies differ only in how a read is weighted: the majority call"
                + " counts every read alike and is therefore carried by whichever type the collection"
                + " has most of, while the Bayes call weighs a read by what its node adds over that"
                + " prior and ignores one that adds nothing. Run both database variants; what the"
                + " refinement is worth is the difference between them.");
    }
}
