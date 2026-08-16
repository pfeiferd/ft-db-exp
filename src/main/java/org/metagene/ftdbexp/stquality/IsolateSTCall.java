package org.metagene.ftdbexp.stquality;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * The sequence-type result for one isolate, accumulated read by read.
 * <p>
 * Every node of the database is treated as a predictor: the types still in question at it are the
 * ones its genomes carry, and it predicts the most frequent of them. A read therefore always yields
 * a prediction --- the majority type of wherever it was classified to --- and whether that
 * prediction is right is a plain comparison against the isolate's type. Nothing has to be thresholded
 * and no read abstains, which is what makes the resulting rate a per-read accuracy rather than a
 * convention.
 * <p>
 * Two numbers follow and they answer different questions. The <em>read accuracy</em> is the share of
 * classified reads whose node predicts the isolate's type; it measures the database, since it is the
 * nodes that predict. The <em>isolate call</em> is the type most reads predicted, and whether it is
 * right is the clinical question --- an isolate is reported once, not once per read.
 * <p>
 * A refinement acts on both through the same mechanism. Before it, most reads rest on a node
 * spanning every lineage, whose majority type is simply the most frequent type of the collection: a
 * prior, right for isolates of that lineage and wrong for the rest. Refined nodes narrow the set the
 * majority is taken over, and the accuracy rises exactly as far as that narrowing carries. The
 * candidate count is kept beside the counts so the two can be read together.
 */
public class IsolateSTCall {
    /** What became of one isolate. */
    public enum Verdict {
        /** The type most reads predicted is the isolate's. */
        CORRECT,
        /** The type most reads predicted is another one. */
        WRONG,
        /** No read was classified, so nothing was predicted. */
        NO_CALL,
        /** The isolate carries no ground truth, so nothing can be said about it. */
        UNTYPED
    }

    private final String key;
    private final String trueST;
    private final Map<String, Long> predictions = new LinkedHashMap<>();
    /** Reads per node, by dense node position; what the naive Bayes classifier is scored from. */
    private final Map<Integer, Long> readsPerNode = new LinkedHashMap<>();
    private long reads;
    private long classified;
    private long correctReads;
    /** Sum of the candidate counts of the nodes the reads landed at, for the mean below. */
    private long candidateSum;

    /**
     * Creates the tally.
     *
     * @param key    the fastq map key, which identifies the isolate
     * @param trueST the isolate's sequence type, or {@code null} if it is not typed
     */
    public IsolateSTCall(String key, String trueST) {
        this.key = key;
        this.trueST = trueST;
    }

    /** Records a read the matcher did not classify at all. */
    public void recordUnclassified() {
        reads++;
    }

    /**
     * Records a classified read.
     *
     * @param nodePos     the dense position of the node it was classified to
     * @param predictedST the majority sequence type of that node, or {@code null} if no genome
     *                    below it carries a type at all
     * @param candidates  how many types are still in question at that node
     */
    public void recordClassified(int nodePos, String predictedST, int candidates) {
        reads++;
        classified++;
        candidateSum += candidates;
        readsPerNode.merge(nodePos, 1L, Long::sum);
        if (predictedST != null) {
            predictions.merge(predictedST, 1L, Long::sum);
            if (predictedST.equals(trueST)) {
                correctReads++;
            }
        }
    }

    /**
     * Adds another tally of the same isolate, which is how the per-thread tallies of one fastq file
     * are merged. Every field is a count, so the sum is exactly what one thread would have produced.
     *
     * @param other the tally to add
     */
    public void add(IsolateSTCall other) {
        reads += other.reads;
        classified += other.classified;
        correctReads += other.correctReads;
        candidateSum += other.candidateSum;
        for (Map.Entry<String, Long> e : other.predictions.entrySet()) {
            predictions.merge(e.getKey(), e.getValue(), Long::sum);
        }
        for (Map.Entry<Integer, Long> e : other.readsPerNode.entrySet()) {
            readsPerNode.merge(e.getKey(), e.getValue(), Long::sum);
        }
    }

    /**
     * Empties the tally for reuse. The per-thread tallies are reset rather than discarded once a
     * file is merged: a matcher thread holds its tally in a {@code ThreadLocal} and would go on
     * writing into a discarded one, so the next file's reads from that thread would be lost.
     */
    public void reset() {
        predictions.clear();
        readsPerNode.clear();
        reads = 0;
        classified = 0;
        correctReads = 0;
        candidateSum = 0;
    }

    /**
     * Returns the type most reads predicted, which is the isolate's call.
     * <p>
     * A tie goes to the type seen first, which is deterministic for a given run but arbitrary
     * between runs; a tie at the top means the reads did not decide and the result should not be
     * read as if they had.
     *
     * @return the called type, or {@code null} if no read predicted anything
     */
    public String getCalledST() {
        String best = null;
        long bestVotes = -1;
        for (Map.Entry<String, Long> e : predictions.entrySet()) {
            if (e.getValue() > bestVotes) {
                best = e.getKey();
                bestVotes = e.getValue();
            }
        }
        return best;
    }

    /**
     * Returns how many reads predicted the called type.
     *
     * @return the called type's reads, or 0 if none
     */
    public long getCalledReads() {
        String called = getCalledST();
        return called == null ? 0 : predictions.get(called);
    }

    /**
     * Returns the called type's share of the reads that predicted anything, i.e. how decided the
     * call was. Reported rather than thresholded.
     *
     * @return the share in {@code [0, 1]}, or {@link Double#NaN} if nothing was predicted
     */
    public double getCalledShare() {
        long total = 0;
        for (long v : predictions.values()) {
            total += v;
        }
        return total == 0 ? Double.NaN : ((double) getCalledReads()) / total;
    }

    /**
     * Returns the share of classified reads whose node predicts the isolate's type.
     *
     * @return the read accuracy, or {@link Double#NaN} if no read was classified
     */
    public double getReadAccuracy() {
        return classified == 0 ? Double.NaN : ((double) correctReads) / classified;
    }

    /**
     * Returns the mean number of types still in question at the nodes the reads landed at. It is
     * what a refinement reduces, and it explains a read accuracy without having to guess at one.
     *
     * @return the mean candidate count, or {@link Double#NaN} if no read was classified
     */
    public double getMeanCandidates() {
        return classified == 0 ? Double.NaN : ((double) candidateSum) / classified;
    }

    /**
     * Decides the isolate on the given call.
     *
     * @param calledST the type the classifier settled on, or {@code null} if it could not
     * @return the verdict
     */
    public Verdict getVerdict(String calledST) {
        if (trueST == null) {
            return Verdict.UNTYPED;
        }
        if (calledST == null) {
            return Verdict.NO_CALL;
        }
        return calledST.equals(trueST) ? Verdict.CORRECT : Verdict.WRONG;
    }

    /**
     * Decides the isolate on the plain majority of the reads' predictions, which is the classifier
     * without the weighting -- reported beside the Bayes call so that what the weighting is worth
     * can be read off the difference.
     *
     * @return the verdict of the majority call
     */
    public Verdict getMajorityVerdict() {
        return getVerdict(getCalledST());
    }

    /**
     * Returns the fastq map key identifying the isolate.
     *
     * @return the key
     */
    public String getKey() {
        return key;
    }

    /**
     * Returns the isolate's sequence type as the ground truth has it.
     *
     * @return the type, or {@code null} if it is not typed
     */
    public String getTrueST() {
        return trueST;
    }

    /**
     * Returns how many reads were seen.
     *
     * @return the number of reads
     */
    public long getReads() {
        return reads;
    }

    /**
     * Returns how many reads the matcher classified anywhere in the database.
     *
     * @return the number of classified reads
     */
    public long getClassified() {
        return classified;
    }

    /**
     * Returns how many classified reads landed on a node predicting the isolate's type.
     *
     * @return the number of correctly predicted reads
     */
    public long getCorrectReads() {
        return correctReads;
    }

    /**
     * Returns how many reads landed at each node, by dense node position. It is what
     * {@link NaiveBayesSTModel#classify} is given, and keeping it per node rather than scoring each
     * read as it arrives is what makes the classifier cost one pass over the types at the end
     * instead of one per read.
     *
     * @return the reads per node position
     */
    public Map<Integer, Long> getReadsPerNode() {
        return readsPerNode;
    }

    /**
     * Returns the reads per predicted type.
     *
     * @return the prediction map, in the order the types were first predicted
     */
    public Map<String, Long> getPredictions() {
        return predictions;
    }
}
