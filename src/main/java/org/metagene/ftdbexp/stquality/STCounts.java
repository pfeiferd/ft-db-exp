package org.metagene.ftdbexp.stquality;

import java.io.Serializable;
import java.util.BitSet;

/**
 * The per-node tally behind the sequence-type level precision of the Genestrip-FT paper.
 * <p>
 * It is the counterpart of {@code DBQualityCountsGoal.Counts} with the unit replaced, which is all
 * the paper's definition asks for: let {@code S_n} be the sequence types represented by the genomes
 * below a node {@code n} ({@code n} included), and for a k-mer {@code a} stored at {@code v(a)} let
 * {@code c_st(a)} be the number of types in {@code S_v(a)} at least one of whose genomes contains
 * {@code a}; then
 * <pre>
 *     p_st(a) = c_st(a) / |S_v(a)| .
 * </pre>
 * Node precision, subtree precision and the restricted variant {@code sp*} follow from it by the same
 * averages as at the data-taxon level.
 * <p>
 * One thing does not carry over from that level, and it is the reason this is a separate class rather
 * than the same one with a different number in it: {@code |D_n|} is a count of leaves and aggregates
 * up the tree by addition, while {@code |S_n|} is the size of a <em>union</em> - two children may
 * carry the same type - and does not. Hence {@link #stSet}, which is unioned into the parent instead
 * of summed, and only turned into a number once the tree has been walked.
 */
public class STCounts implements Serializable {
    private static final long serialVersionUID = 1L;

    /** Whether this tally belongs to a leaf, i.e. to the deepest artificial node on its branch. */
    private final boolean forLeaf;

    /** Number of k-mers the database stores under exactly this tax id, not counting descendants. */
    private final long kmerSumForNode;

    /**
     * The sequence types represented below this node, this node included - the set {@code S_n}. Bits
     * are the dense indices of {@link STGroundTruth}. A leaf sets its own genome's type, an inner node
     * receives the union of its descendants'.
     */
    private final BitSet stSet = new BitSet();

    /**
     * Sum over the k-mers {@code a} stored at exactly this node of {@code c_st(a)}, i.e. of the number
     * of sequence types at least one of whose genomes carries {@code a}. Counted once per
     * {@code (k-mer, sequence type)} pair, which is what makes it {@code c_st} rather than a count of
     * genomes.
     */
    private long tpForNodePrecision;

    /**
     * True positives pooled over the sequence types below this node: {@code (k-mer, type)} pairs whose
     * k-mer the database stores on the path from the genome it was read in. Filled once, after the
     * reading pass, from the per-type tallies - see {@code STQualityCountsGoal}.
     */
    private long tp;
    /** The same pairs whether or not they are stored on that path, i.e. {@code tp + fn}. */
    private long tpPlusFn;
    /** Sum of the per-type recalls over the types below this node, for the unweighted average. */
    private double recallSum;
    /** Number of types that contributed to {@link #recallSum}. */
    private int recallTypes;

    /** Sum of {@code p_st(a)} over the k-mers of the subtree rooted here for which it is defined. */
    private double subtreePrecisionSum;
    /** Number of k-mers of that subtree, i.e. the denominator of {@link #getSubtreePrecision()}. */
    private long subtreeKmerSum;
    /** Sum of {@code p_st(a)} over the k-mers of that subtree stored above the data taxa. */
    private double subtreePrecisionSumAboveData;
    /** Number of k-mers of that subtree stored above the data taxa. */
    private long subtreeKmerSumAboveData;

    /**
     * Creates an empty tally.
     *
     * @param forLeaf        whether the node is the deepest artificial node on its branch
     * @param kmerSumForNode the number of k-mers the database stores under exactly this tax id
     */
    public STCounts(boolean forLeaf, long kmerSumForNode) {
        this.forLeaf = forLeaf;
        this.kmerSumForNode = kmerSumForNode;
    }

    /**
     * Returns whether this tally belongs to a leaf node.
     *
     * @return whether the node is the deepest artificial node on its branch
     */
    public boolean isForLeaf() {
        return forLeaf;
    }

    /**
     * Records that the genomes below this node represent the given sequence type.
     *
     * @param stIndex the dense index of the sequence type
     */
    void addST(int stIndex) {
        stSet.set(stIndex);
    }

    /**
     * Adds the sequence types of a descendant to this node's set, which is how {@code S_n} is built:
     * by union and not by addition, since two descendants may share a type.
     *
     * @param descendant the tally of a node below this one
     */
    void unionSTs(STCounts descendant) {
        stSet.or(descendant.stSet);
    }

    /**
     * Returns {@code |S_n|}, the number of distinct sequence types represented below this node.
     *
     * @return the number of sequence types below this node, this node included
     */
    public int getSTCount() {
        return stSet.cardinality();
    }

    /**
     * Returns the index of the first sequence type below this node at or after {@code from}, or
     * {@code -1} if there is none. Lets a caller walk {@code S_n} without a copy of the set.
     *
     * @param from the index to start looking at
     * @return the next sequence type index, or {@code -1}
     */
    public int nextST(int from) {
        return stSet.nextSetBit(from);
    }

    /**
     * Adds one sequence type's read counts to this node, which is how the recall of a node is pooled
     * from the lineages below it.
     *
     * @param typeTp       the type's true positives
     * @param typeTpPlusFn the type's true positives plus false negatives
     */
    void addTypeRecall(long typeTp, long typeTpPlusFn) {
        tp += typeTp;
        tpPlusFn += typeTpPlusFn;
        if (typeTpPlusFn > 0) {
            recallSum += ((double) typeTp) / typeTpPlusFn;
            recallTypes++;
        }
    }

    /**
     * Returns the number of {@code (k-mer, sequence type)} pairs found below this node whose k-mer the
     * database stores on the path from the genome it was read in.
     *
     * @return the pooled true positives
     */
    public long getTp() {
        return tp;
    }

    /**
     * Returns the number of {@code (k-mer, sequence type)} pairs found below this node, on the path or
     * not.
     *
     * @return the pooled true positives plus false negatives
     */
    public long getTpPlusFn() {
        return tpPlusFn;
    }

    /**
     * Returns the recall pooled over the lineages below this node: of the database's k-mers that their
     * genomes carry, the share stored where those genomes can claim them.
     * <p>
     * Pooled and therefore weighted by how many k-mers each lineage contributes, which is the same
     * convention the data-taxon level uses. Path correctness makes this one for a database whose fill
     * saw every genome; it falls below one where it did not - a genome outside the reference set, or
     * one the per-taxon limits dropped - and that is what it is worth reporting for.
     *
     * @return the pooled recall, or {@link Double#NaN} if no pair was found below this node
     */
    public double getRecall() {
        if (tpPlusFn == 0) {
            return Double.NaN;
        }
        return ((double) tp) / tpPlusFn;
    }

    /**
     * Returns the mean of the per-lineage recalls below this node, each lineage counting once however
     * many genomes or k-mers it brings. It is the measure to read when the collection is as lopsided
     * as a genome collection usually is - here ten of 196 types hold well over half the genomes.
     *
     * @return the unweighted average recall, or {@link Double#NaN} if no lineage below has a recall
     */
    public double getAvgRecall() {
        if (recallTypes == 0) {
            return Double.NaN;
        }
        return recallSum / recallTypes;
    }

    /**
     * Records one more {@code (k-mer, sequence type)} pair for a k-mer stored at exactly this node.
     */
    void incTpForNodePrecision() {
        tpForNodePrecision++;
    }

    /**
     * Adds one reader's tally for this node.
     *
     * @param count the number of {@code (k-mer, sequence type)} pairs that reader found for this node
     */
    void addTpForNodePrecision(long count) {
        tpForNodePrecision += count;
    }

    /**
     * Returns the accumulated {@code c_st} over this node's own k-mers.
     *
     * @return the number of distinct {@code (k-mer, sequence type)} pairs found for this node
     */
    public long getTpForNodePrecision() {
        return tpForNodePrecision;
    }

    /**
     * Returns the number of k-mers stored at exactly this node.
     *
     * @return the k-mers the database stores under this tax id, not counting descendants
     */
    public long getKmerSumForNode() {
        return kmerSumForNode;
    }

    /**
     * Returns this node's precision at the sequence-type level: the mean of {@code p_st(a)} over the
     * k-mers stored at this node alone, which is {@code tp / (|S_n| * kmers at node)}.
     * <p>
     * Undefined, and hence {@code NaN}, for a node that holds no k-mer or has no typed genome
     * underneath. Deliberately not extended to those by a convention such as {@code 1}: a node
     * storing no k-mer makes no statement about k-mer placement, and a perfect score would raise the
     * subtree precision above it.
     *
     * @return the node precision, or {@link Double#NaN} where it is undefined
     */
    public double getNodePrecision() {
        int sts = getSTCount();
        if (sts == 0 || kmerSumForNode == 0) {
            return Double.NaN;
        }
        return ((double) tpForNodePrecision) / (((double) sts) * kmerSumForNode);
    }

    /**
     * Returns the subtree precision at the sequence-type level: the mean of {@code p_st(a)} over every
     * k-mer residing in the subtree rooted here.
     *
     * @return the subtree precision, or {@link Double#NaN} if the subtree holds no such k-mer
     */
    public double getSubtreePrecision() {
        if (subtreeKmerSum == 0) {
            return Double.NaN;
        }
        return subtreePrecisionSum / subtreeKmerSum;
    }

    /**
     * Returns the number of k-mers {@link #getSubtreePrecision()} averages over.
     *
     * @return the k-mers of this subtree for which {@code p_st} is defined
     */
    public long getSubtreeKmerSum() {
        return subtreeKmerSum;
    }

    /**
     * Returns the subtree precision restricted to the k-mers stored above the data taxa, which is
     * where a refinement can act at all: a k-mer sitting at a single genome has {@code |S_n| = 1} and
     * hence {@code p_st = 1} by construction, and such k-mers are the bulk of a database, so they
     * dominate the unrestricted average while being incapable of improvement.
     *
     * @return the restricted subtree precision, or {@link Double#NaN} if there are no such k-mers
     */
    public double getRestrictedSubtreePrecision() {
        if (subtreeKmerSumAboveData == 0) {
            return Double.NaN;
        }
        return subtreePrecisionSumAboveData / subtreeKmerSumAboveData;
    }

    /**
     * Returns the number of k-mers {@link #getRestrictedSubtreePrecision()} averages over.
     *
     * @return the k-mers of this subtree stored above its data taxa
     */
    public long getSubtreeKmersAboveData() {
        return subtreeKmerSumAboveData;
    }

    /**
     * Adds the contribution of one node of this node's subtree - possibly this node itself - to the
     * pooled sums behind the two subtree averages.
     * <p>
     * Every k-mer of the database sits at exactly one node, so summing per node over the subtree
     * averages over precisely the subtree's k-mers, each counted once. A node holding no k-mer
     * contributes nothing, which is what makes the average insensitive to how many nodes the taxonomy
     * or the refinement happens to provide.
     *
     * @param counts the tally of a node of this subtree, which must have at least one typed genome
     *               underneath
     */
    void aggregateSubtree(STCounts counts) {
        if (counts.kmerSumForNode > 0) {
            // Sum of p_st(a) = c_st(a) / |S_n| over the k-mers a stored at that node. The caller
            // guarantees |S_n| > 0, so p_st is defined for all of them.
            double precisionSum = ((double) counts.tpForNodePrecision) / counts.getSTCount();
            subtreePrecisionSum += precisionSum;
            subtreeKmerSum += counts.kmerSumForNode;
            if (!counts.forLeaf) {
                subtreePrecisionSumAboveData += precisionSum;
                subtreeKmerSumAboveData += counts.kmerSumForNode;
            }
        }
    }
}
