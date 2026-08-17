package org.metagene.ftdbexp.stquality;

import org.metagene.genestrip.tax.SmallTaxTree;

import java.util.ArrayList;
import java.util.List;
import java.util.Collections;
import java.util.Map;
import java.util.Set;

/**
 * The classifier Genestrip uses on a single read, applied to a whole isolate at once.
 * <p>
 * Genestrip classifies a read the way Kraken does: the read's {@code k}-mers are counted at the nodes
 * they are stored at, each node that received one is a candidate, and the candidate whose root path
 * carries the largest total wins -- {@code FastqKMerMatcher#sumCounts} sums a node together with all
 * of its ancestors. Candidates that tie are resolved to their lowest common ancestor, and where a
 * minimum is configured a candidate is first raised to the lowest ancestor at which the running total
 * reaches it. This class runs exactly that rule, only over counts pooled across every read of one
 * isolate rather than within one read, and it takes those counts in either of two currencies:
 * <ul>
 * <li>the matched {@code k}-mers of the isolate per node, and</li>
 * <li>the reads of the isolate per node, i.e. one vote per read where {@link NaiveBayesSTModel}
 *     would have taken one factor per read.</li>
 * </ul>
 * The two answer different questions and are reported side by side, as the {@code k}-mer and read
 * levels of the Bayes classifier are.
 * <p>
 * Two properties of the rule make it worth trying on this data, and both are consequences of summing
 * along the root path rather than at the node itself.
 * <p>
 * The first is that the uninformative mass cancels. Nearly every {@code k}-mer of a sample sits at
 * the node holding what the species shares -- for cdiff that node holds ≈ 95 percent of everything
 * above the genomes -- and that node is an ancestor of every candidate below it. Its count therefore
 * enters every candidate's total identically and cannot affect which is larger. Where the naive Bayes
 * classifier had to be arranged so that such a node contributes a factor of exactly one, here it
 * falls out of the comparison on its own.
 * <p>
 * The second is that a tie is answered with a lineage. Two reference genomes that attract the same
 * total resolve to their lowest common ancestor: in a refined database that is the node grouping
 * them, which carries a sequence type; in an unrefined one there is nothing between a genome and the
 * species, so the same tie resolves to the species and says nothing. The refinement's contribution to
 * this classifier is therefore structural rather than a matter of degree, which is what the case
 * study needs to establish and what comparing the two databases under it should show.
 * <p>
 * What the rule does not do is weigh a class against how frequent it is in the collection. It reports
 * where the isolate's sequence sits in the tree, and the lineage follows from the node. A node
 * holding a single reference genome therefore answers with that genome's type on the strength of one
 * genome, which is nearest-neighbour reasoning and is meant here: it is the baseline the weighted
 * classifier has to beat, and in an unrefined database it is all there is.
 */
public class PathVoteSTModel extends STClassifier {
    /** Take every candidate at the node it was counted at, i.e. no minimum. */
    public static final long NO_MINIMUM = 1;

    /**
     * The total a candidate's root path must reach before it is taken. A rule with a minimum and one
     * without are two rules, not one rule configured twice, so it is fixed when the classifier is
     * made and the results are labelled accordingly.
     */
    private final long minimum;

    /**
     * Creates the classifier over an already counted composition, taking every candidate where it
     * stands.
     *
     * @param composition the type composition of the database's nodes
     */
    public PathVoteSTModel(STComposition composition) {
        this(composition, NO_MINIMUM);
    }

    /**
     * Creates the classifier over an already counted composition.
     *
     * @param composition the type composition of the database's nodes
     * @param minimum     the total a candidate's root path must reach before it is taken, raising it
     *                    towards the root until it does; {@link #NO_MINIMUM} to take every candidate
     *                    where it stands
     */
    public PathVoteSTModel(STComposition composition, long minimum) {
        super(composition);
        this.minimum = minimum;
    }

    @Override
    public String getName() {
        return minimum > NO_MINIMUM ? "path/min" + minimum : "path";
    }

    /**
     * Names the lineage of one isolate, i.e. what the node it is classified to predicts.
     *
     * @param countsPerNode how much the isolate contributed at each node, by dense node position
     * @return the sequence type, or {@code null} if the isolate is not classified
     */
    @Override
    public String classify(Map<Integer, Long> countsPerNode) {
        SmallTaxTree.SmallTaxIdNode node = classifyNode(countsPerNode, minimum);
        return node == null ? null : composition.getMajorityST(node.getPosition());
    }

    /**
     * Names every lineage the winning node leaves in question.
     * <p>
     * The node an isolate lands on need not hold a single type, and where it holds several the
     * honest answer is all of them: the rule has placed the isolate there and no further. The
     * isolate is then classified to each of them at once and scored accordingly -- see
     * {@link #getPrecision}.
     *
     * @param countsPerNode how much the isolate contributed at each node, by dense node position
     * @return the types below the winning node, empty if the isolate is not classified
     */
    @Override
    public Set<String> classifyAll(Map<Integer, Long> countsPerNode) {
        SmallTaxTree.SmallTaxIdNode node = classifyNode(countsPerNode, minimum);
        return node == null ? Collections.<String>emptySet()
                : composition.getTypesAt(node.getPosition());
    }

    /**
     * Classifies an isolate to a node, by the rule described for this class.
     *
     * @param countsPerNode how much the isolate contributed at each node, by dense node position;
     *                      reads or {@code k}-mers, whichever currency is being scored
     * @param minimum       the total a candidate's root path must reach before it is taken, raising
     *                      it towards the root until it does; {@link #NO_MINIMUM} to take every
     *                      candidate where it stands
     * @return the node the isolate is classified to, or {@code null} if nothing was counted or no
     *         candidate reaches the minimum anywhere
     */
    public SmallTaxTree.SmallTaxIdNode classifyNode(Map<Integer, Long> countsPerNode, long minimum) {
        if (countsPerNode.isEmpty()) {
            return null;
        }
        long best = Long.MIN_VALUE;
        List<SmallTaxTree.SmallTaxIdNode> winners = new ArrayList<>();
        for (Map.Entry<Integer, Long> e : countsPerNode.entrySet()) {
            if (e.getValue() <= 0) {
                continue;
            }
            SmallTaxTree.SmallTaxIdNode node = composition.getNode(e.getKey());
            if (node == null) {
                // A count for a node this composition does not know is not something to pass over
                // quietly: it means the counts and the tree come from different databases, and the
                // result would be a classification of one against the shape of the other.
                throw new IllegalStateException("There is a count for node position " + e.getKey()
                        + ", which this database's taxonomy does not have. The counts and the tree"
                        + " come from different databases.");
            }
            long sum = sumOnPath(node, countsPerNode);
            if (sum > best) {
                best = sum;
                winners.clear();
                winners.add(node);
            } else if (sum == best) {
                winners.add(node);
            }
        }
        if (winners.isEmpty()) {
            return null;
        }
        if (minimum > NO_MINIMUM) {
            for (int i = 0; i < winners.size(); i++) {
                SmallTaxTree.SmallTaxIdNode raised = lowestNodeReaching(winners.get(i), countsPerNode, minimum);
                if (raised == null) {
                    // Not even the root gathers enough: the isolate is not classified rather than
                    // classified to whatever happened to be nearest.
                    return null;
                }
                winners.set(i, raised);
            }
        }
        SmallTaxTree.SmallTaxIdNode node = winners.get(0);
        for (int i = 1; i < winners.size() && node != null; i++) {
            node = composition.getTree().getLowestCommonAncestor(node, winners.get(i));
        }
        return node;
    }

    /**
     * Sums what the isolate contributed at the given node and at every ancestor of it.
     *
     * @param node          the candidate
     * @param countsPerNode the isolate's counts by dense node position
     * @return the total on the root path
     */
    private long sumOnPath(SmallTaxTree.SmallTaxIdNode node, Map<Integer, Long> countsPerNode) {
        long sum = 0;
        for (SmallTaxTree.SmallTaxIdNode n = node; n != null; n = n.getParent()) {
            Long c = countsPerNode.get(n.getPosition());
            if (c != null) {
                sum += c;
            }
        }
        return sum;
    }

    /**
     * Walks from the given node towards the root and returns the first node at which the running
     * total reaches the minimum.
     *
     * @param node          the candidate to raise
     * @param countsPerNode the isolate's counts by dense node position
     * @param minimum       the total to reach
     * @return the lowest such node, or {@code null} if the root is reached without getting there
     */
    private SmallTaxTree.SmallTaxIdNode lowestNodeReaching(SmallTaxTree.SmallTaxIdNode node,
                                                           Map<Integer, Long> countsPerNode, long minimum) {
        long sum = 0;
        for (SmallTaxTree.SmallTaxIdNode n = node; n != null; n = n.getParent()) {
            Long c = countsPerNode.get(n.getPosition());
            if (c != null) {
                sum += c;
                if (sum >= minimum) {
                    return n;
                }
            }
        }
        return null;
    }
}
