package org.metagene.ftdbexp.taxquality;

import org.metagene.genestrip.tax.SmallTaxTree;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Names a sample's data taxon the way Kraken names a read's: by the heaviest root-to-node path.
 * <p>
 * Each node the sample contributed at is scored with the sum of the contributions along its path to
 * the root, and the heaviest wins; ties are resolved to the lowest common ancestor of the tied
 * nodes, so a tie widens the answer instead of picking one of its halves arbitrarily. The rule is
 * the C.~difficile study's {@code PathVoteSTModel} with one simplification: there the winning node
 * still had to be turned into a sequence type, by taking the majority type below it, because a
 * sequence type is not a place in the tree. A data taxon is, so the winning node <em>is</em> the
 * answer wherever it is one -- and where it is not, the rule says so rather than guessing.
 * <p>
 * A minimum may be required of the winner. Where the heaviest path does not gather it, the answer is
 * raised to the lowest ancestor that does; where not even the root gathers it, the sample is left
 * unclassified rather than classified to whatever happened to be nearest. That is the knob by which
 * this rule trades specificity for confidence, and it is the reason it is in the comparison at all:
 * the naive Bayes rule of {@link NaiveBayesTaxonModel} always commits.
 */
public class PathVoteTaxonModel extends TaxonClassifier {
    /** Value of {@code minimum} meaning that no minimum is required. */
    public static final long NO_MINIMUM = 0;

    private final long minimum;

    /**
     * @param composition the data-taxon composition of the database's nodes
     * @param minimum     how much the winning path must gather, or {@link #NO_MINIMUM}
     */
    public PathVoteTaxonModel(TaxonComposition composition, long minimum) {
        super(composition);
        this.minimum = minimum;
    }

    @Override
    public String classify(Map<Integer, Long> countsPerNode) {
        SmallTaxTree.SmallTaxIdNode node = classifyNode(countsPerNode, minimum);
        if (node == null) {
            return null;
        }
        // The winning node is the answer when it is a data taxon or lies below one. Above one it is
        // not: a node spanning several species does not name a species, and saying which of them is
        // the majority would answer a question the reads did not settle. The majority is available
        // through the composition for a caller that wants it, but this rule does not take it.
        String taxon = composition.getMajorityClass(node.getPosition());
        return composition.getCandidates(node.getPosition()) == 1 ? taxon : null;
    }

    /**
     * The node the heaviest root-to-node path ends at, before it is turned into a taxon.
     * <p>
     * Exposed because it is the more informative answer of the two: a caller comparing the unrefined
     * database with the refined one wants to see <em>where</em> the vote landed, and after a
     * refinement that may be a node the unrefined tree did not have.
     *
     * @param countsPerNode how much the sample contributed at each node, by dense node position
     * @param minimum       how much the winning path must gather, or {@link #NO_MINIMUM}
     * @return the winning node, or {@code null} where the rule makes no call
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
     * @param node          the node whose path is summed
     * @param countsPerNode the sample's contributions
     * @return the sum of the contributions from {@code node} up to the root
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
     * @param node          the node to start from
     * @param countsPerNode the sample's contributions
     * @param minimum       the amount the path must gather
     * @return the lowest ancestor of {@code node}, itself included, whose path gathers
     *         {@code minimum}, or {@code null} if not even the root does
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

    @Override
    public String getName() {
        return minimum > NO_MINIMUM ? "path vote (min " + minimum + ")" : "path vote";
    }
}
