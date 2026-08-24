package org.metagene.ftdbexp.taxquality;

import org.metagene.genestrip.tax.SmallTaxTree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * Names a sample's species by a plain majority over what the matcher assigned: the node most of the
 * sample's reads -- or k-mers, depending on the currency -- were classified to.
 * <p>
 * Each node is worth its own count and nothing more. Ties are resolved to the lowest common ancestor
 * of the tied nodes, so a tie widens the answer instead of picking one of its halves arbitrarily.
 * <p>
 * Not the sum along the node's path to the root, which this rule carried over from the C.~difficile
 * study at first and which is wrong here. That sum is Kraken's, and Kraken's is right -- for a
 * <em>read</em>. There the weights on a path are the k-mers of one sequence, so a k-mer resting at
 * the genus really is evidence about the same read as a k-mer resting at the species. A sample is not
 * a read: its counts are already aggregated over millions of reads, and those that stopped at a
 * coarse node are not those that reached a fine one. Kraken scores paths within one sequence and
 * aggregates afterwards, never the other way round.
 * <p>
 * What the path sum cost here is worth recording. The {@code strepto} database carries an
 * {@code unclassified Streptococcus} node with 468 species below it, where the LCA update parks every
 * read it cannot place further. Its whole weight was credited to each of those 468 in full, so a
 * draft genome with a single read of its own outscored {@code S. pneumoniae} with fifty, and 178 of
 * 243 calls landed somewhere in that one branch. Under a majority the node wins for itself, and the
 * answer is the honest "some unnamed Streptococcus" rather than an arbitrary one of its children.
 * <p>
 * A minimum may be required of the winner. Where it does not gather that on its own, the answer is
 * raised to the lowest ancestor whose path does -- upwards, where the counts of the coarser nodes are
 * the reads that were classified there. Where not even the root gathers it, the sample is left
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
        // The winning node is the answer when it is a species or lies below one. Above one it is
        // not: a node spanning several species does not name a species, and saying which of them is
        // the majority would answer a question the reads did not settle. The majority is available
        // through the composition for a caller that wants it, but this rule does not take it.
        String taxon = composition.getMajorityClass(node.getPosition());
        return composition.getCandidates(node.getPosition()) == 1 ? taxon : null;
    }

    /**
     * The node the majority names, before it is turned into a taxon.
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
        List<SmallTaxTree.SmallTaxIdNode> best = classifyNodes(countsPerNode, minimum, 1);
        return best.isEmpty() ? null : best.get(0);
    }

    /**
     * The {@code howMany} nodes with the most reads behind them, heaviest first.
     * <p>
     * One winner is the wrong shape of answer for a sputum sample. Such a sample carries a community
     * of streptococci -- one with more than a thousand reads spreads them over some 320 nodes of the
     * genus -- and the reference standard says so too, naming several organisms for nine of the
     * thirteen samples it decides. A rule that crowns one of them answers a question the sample does
     * not settle, and the runner-up may be the organism actually asked about. Reporting the first few
     * leaves that judgement to the reader instead of making it silently.
     * <p>
     * Nodes that score alike are one answer, not several: each group of equal weight is resolved to
     * the lowest common ancestor of its members, exactly as a tie is resolved for a single winner, so
     * a rank widens rather than picking arbitrarily from a tie. Groups are then taken from the top
     * until {@code howMany} <em>distinct</em> nodes have been collected -- raising to a minimum can
     * carry two groups onto the same ancestor, and repeating it would suggest support that is not
     * there.
     *
     * @param countsPerNode the sample's contributions by node position
     * @param minimum       how much a winner's path must gather, or {@link #NO_MINIMUM}
     * @param howMany       how many nodes to return at most
     * @return the nodes, heaviest first, possibly fewer than asked for and possibly empty
     */
    public List<SmallTaxTree.SmallTaxIdNode> classifyNodes(Map<Integer, Long> countsPerNode, long minimum,
                                                           int howMany) {
        List<SmallTaxTree.SmallTaxIdNode> result = new ArrayList<>();
        if (countsPerNode.isEmpty() || howMany <= 0) {
            return result;
        }
        // Score every node the sample contributed at, then group the nodes by score.
        SortedMap<Long, List<SmallTaxTree.SmallTaxIdNode>> byScore = new TreeMap<>(Collections.reverseOrder());
        for (Map.Entry<Integer, Long> e : countsPerNode.entrySet()) {
            if (e.getValue() <= 0) {
                continue;
            }
            SmallTaxTree.SmallTaxIdNode node = composition.getNode(e.getKey());
            if (node == null) {
                throw new IllegalStateException("There is a count for node position " + e.getKey()
                        + ", which this database's taxonomy does not have. The counts and the tree"
                        + " come from different databases.");
            }
            // The node's own count and nothing else: a plain majority over what the matcher assigned.
            byScore.computeIfAbsent(e.getValue(), k -> new ArrayList<>()).add(node);
        }
        for (List<SmallTaxTree.SmallTaxIdNode> tied : byScore.values()) {
            if (result.size() >= howMany) {
                break;
            }
            SmallTaxTree.SmallTaxIdNode node = resolve(tied, countsPerNode, minimum);
            if (node != null && !result.contains(node)) {
                result.add(node);
            }
        }
        return result;
    }

    /**
     * Turns one group of equally weighted nodes into the single node that stands for it: each is
     * first raised to where it gathers the minimum, and what remains is folded into its lowest
     * common ancestor.
     *
     * @param tied          the nodes of one score group
     * @param countsPerNode the sample's contributions
     * @param minimum       how much the path must gather, or {@link #NO_MINIMUM}
     * @return the node standing for the group, or {@code null} if the minimum is out of reach
     */
    private SmallTaxTree.SmallTaxIdNode resolve(List<SmallTaxTree.SmallTaxIdNode> tied,
                                                Map<Integer, Long> countsPerNode, long minimum) {
        List<SmallTaxTree.SmallTaxIdNode> winners = new ArrayList<>(tied);
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
