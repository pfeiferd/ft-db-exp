package org.metagene.ftdbexp.eval;

import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.util.IdentityHashMap;
import java.util.Map;

/**
 * Counts, for a node of a database's taxonomy tree, how many species taxa a read classified to that
 * node could still belong to.
 * <p>
 * Plain precision at the genus rank cannot tell an unrefined database from a refined one: a read
 * that ends up at its genus counts as correct in both cases, no matter how many species that genus
 * comprises. What the refinement improves is exactly how far a classification narrows the species
 * down, and that is what this class quantifies. A classification to a node covering three species is
 * worth much more than one to a genus covering forty, and the reciprocal of the candidate count
 * expresses that.
 * <p>
 * The candidates of a node are the species taxa <em>comparable</em> with it, i.e. those below it and
 * -- if there are none -- the species above it, which applies when a read is classified to a strain
 * or another node below the species rank. Species are counted over the entire subtree rather than
 * among the direct children, because a taxonomy may well place intermediate nodes such as a subgenus
 * between a genus and its species. The counts are memoized per node, so the tree is walked once.
 */
public class SpeciesCandidates {
    private final Map<SmallTaxTree.SmallTaxIdNode, Integer> speciesBelowCache =
            new IdentityHashMap<SmallTaxTree.SmallTaxIdNode, Integer>();

    /**
     * Returns how many species taxa a classification to the given node leaves in question.
     *
     * @param node the node a read was classified to
     * @return the number of candidate species, or {@code 0} if the node's subtree holds no species
     * taxon and no species lies above it
     */
    public int candidatesFor(SmallTaxTree.SmallTaxIdNode node) {
        int below = speciesBelow(node);
        return below > 0 ? below : speciesAbove(node);
    }

    /**
     * Returns the weight a read classified to the given node contributes, i.e. the reciprocal of its
     * candidate count. A classification that pins down a single species is worth one, a
     * classification to a genus of forty species is worth a fortieth.
     *
     * @param node the node a read was classified to
     * @return the weight in {@code (0, 1]}, or {@code 0} if the node leaves no species in question
     */
    public double weightFor(SmallTaxTree.SmallTaxIdNode node) {
        int candidates = candidatesFor(node);
        return candidates == 0 ? 0 : 1.0 / candidates;
    }

    /**
     * Returns whether the two nodes lie on a common path to the root, i.e. whether one is an
     * ancestor of the other or they are identical. Only then does a classification say anything
     * about the read's true species.
     *
     * @param a the first node
     * @param b the second node
     * @return whether the nodes are comparable in the tree order
     */
    public static boolean areComparable(SmallTaxTree.SmallTaxIdNode a, SmallTaxTree.SmallTaxIdNode b) {
        return isAncestorOrSelf(a, b) || isAncestorOrSelf(b, a);
    }

    /**
     * Returns whether the first node is an ancestor of the second one or identical to it.
     *
     * @param ancestor  the presumed ancestor
     * @param candidate the node to check
     * @return whether {@code ancestor} lies on the path from {@code candidate} to the root
     */
    private static boolean isAncestorOrSelf(SmallTaxTree.SmallTaxIdNode ancestor,
                                            SmallTaxTree.SmallTaxIdNode candidate) {
        for (SmallTaxTree.SmallTaxIdNode node = candidate; node != null; node = node.getParent()) {
            if (node == ancestor) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns the number of species taxa in the subtree rooted at the given node, including the node
     * itself. Intermediate ranks between the node and the species are passed through, so a subgenus
     * or any other node in between does not hide the species underneath it.
     *
     * @param node the subtree root
     * @return the number of species taxa in the subtree
     */
    private int speciesBelow(SmallTaxTree.SmallTaxIdNode node) {
        Integer cached = speciesBelowCache.get(node);
        if (cached != null) {
            return cached;
        }
        int count = Rank.SPECIES.equals(node.getRank()) ? 1 : 0;
        SmallTaxTree.SmallTaxIdNode[] subNodes = node.getSubNodes();
        if (subNodes != null) {
            for (SmallTaxTree.SmallTaxIdNode subNode : subNodes) {
                if (subNode != null) {
                    count += speciesBelow(subNode);
                }
            }
        }
        speciesBelowCache.put(node, count);
        return count;
    }

    /**
     * Returns the number of species taxa on the path from the given node to the root, which is one
     * for a node below the species rank and zero otherwise.
     *
     * @param node the node to start from
     * @return the number of species taxa above the node
     */
    private int speciesAbove(SmallTaxTree.SmallTaxIdNode node) {
        for (SmallTaxTree.SmallTaxIdNode current = node.getParent(); current != null; current = current.getParent()) {
            if (Rank.SPECIES.equals(current.getRank())) {
                return 1;
            }
        }
        return 0;
    }
}
