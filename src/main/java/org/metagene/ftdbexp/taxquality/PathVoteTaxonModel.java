package org.metagene.ftdbexp.taxquality;

import org.metagene.genestrip.tax.SmallTaxTree;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

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
     * The {@code howMany} species with the most reads behind them, heaviest first.
     * <p>
     * Every count is first carried up to the species it belongs to -- a read assigned to a strain, or
     * to one of the artificial nodes the fill inserts below it, is a read for that strain's species --
     * and the species are then ranked by what they gathered. Rolling up is what keeps a well
     * represented species from being beaten by an obscure one: {@code S. pneumoniae} enters the
     * database with 153 genomes and would otherwise field 153 separate candidates, each with a
     * fraction of the reads, against a draft genome that fields one.
     * <p>
     * Counts that belong to no species are not carried anywhere. A read the classifier left at the
     * genus, at a refined node, or in a bucket such as {@code unclassified Streptococcus} says that
     * the evidence did not reach a species, and crediting it to one would invent a specificity the
     * read does not have. {@link #unplaced} reports how much of the sample that is; it is the figure
     * a refinement should lower, since pushing k-mers down is exactly what lets a read reach a
     * species it could not reach before.
     * <p>
     * Ties are broken by tax id so that a repeated run answers the same. There is no folding into a
     * common ancestor here, unlike a vote over arbitrary nodes: the answer has to be a species, and
     * the ancestor of two species is not one.
     *
     * @param countsPerNode the sample's contributions by node position
     * @param minimum       how much a species must gather to be reported at all, or {@link #NO_MINIMUM}
     * @param howMany       how many species to return at most
     * @return the species nodes, heaviest first, possibly fewer than asked for and possibly empty
     */
    public List<SmallTaxTree.SmallTaxIdNode> classifyNodes(Map<Integer, Long> countsPerNode, long minimum,
                                                           int howMany) {
        List<SmallTaxTree.SmallTaxIdNode> result = new ArrayList<>();
        if (countsPerNode.isEmpty() || howMany <= 0) {
            return result;
        }
        Map<String, Long> perSpecies = countsPerSpecies(countsPerNode);
        List<Map.Entry<String, Long>> ranked = new ArrayList<>(perSpecies.entrySet());
        ranked.sort((a, b) -> {
            int byCount = Long.compare(b.getValue(), a.getValue());
            return byCount != 0 ? byCount : a.getKey().compareTo(b.getKey());
        });
        for (Map.Entry<String, Long> e : ranked) {
            if (result.size() >= howMany) {
                break;
            }
            if (minimum > NO_MINIMUM && e.getValue() < minimum) {
                // Ranked by size, so nothing further down can reach it either.
                break;
            }
            SmallTaxTree.SmallTaxIdNode node = composition.getTree().getNodeByTaxId(e.getKey());
            if (node != null) {
                result.add(node);
            }
        }
        return result;
    }

    /**
     * The sample's counts summed per species, leaving out what belongs to no species.
     *
     * @param countsPerNode the sample's contributions by node position
     * @return the count each species gathered, keyed by tax id
     */
    public Map<String, Long> countsPerSpecies(Map<Integer, Long> countsPerNode) {
        Map<String, Long> perSpecies = new LinkedHashMap<>();
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
            String species = TaxonComposition.speciesOf(node);
            if (!TaxonComposition.UNPLACED.equals(species)) {
                perSpecies.merge(species, e.getValue(), Long::sum);
            }
        }
        return perSpecies;
    }

    /**
     * How much of the sample reached no species at all: the counts at the genus, at refined nodes and
     * in the unranked buckets between them.
     *
     * @param countsPerNode the sample's contributions by node position
     * @return the sum of the counts that belong to no species
     */
    public long unplaced(Map<Integer, Long> countsPerNode) {
        long total = 0;
        long placed = 0;
        for (Long v : countsPerNode.values()) {
            if (v > 0) {
                total += v;
            }
        }
        for (Long v : countsPerSpecies(countsPerNode).values()) {
            placed += v;
        }
        return total - placed;
    }

    @Override
    public String getName() {
        return minimum > NO_MINIMUM ? "species vote (min " + minimum + ")" : "species vote";
    }
}
