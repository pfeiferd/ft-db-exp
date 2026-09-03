package org.metagene.ftdbexp.taxquality;

import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.util.*;

/**
 * Which species sit below each node of a database, counted once and shared by the classifiers that
 * need it.
 * <p>
 * This is the taxonomic counterpart of the {@code STComposition} the C.~difficile study used. There
 * the class of a genome was its seven-locus sequence type, which had to be established outside the
 * database and carried in from a typing file. Here the class is the genome's <em>species</em>, which
 * the tree already knows, so nothing has to be read in beside it. The rest is unchanged, and
 * deliberately so: every classifier weighs the same composition, counted once, so that what separates
 * two answers is the rule and never two readings of the database.
 * <p>
 * The species and not the data taxon the genome is filed at, although the fill files it deeper. Two
 * things ask for the species: the candidate precision of the paper averages {@code 1/|Sigma(k)|} over
 * the species of a node's subtree, and the reference standard a call is held against names a species.
 * A composition of strains answers neither -- see {@link #speciesOf}.
 * <p>
 * A leaf whose lineage carries no rank at or below species is counted as {@link #UNPLACED} rather
 * than dropped, for the same reason the ST composition carried its untyped genomes: left out, they
 * would silently inflate the confidence of whatever classes remain, since a node holding ten placed
 * and ten unplaced genomes would speak for the placed ones as though nothing else sat below it.
 */
public class TaxonComposition {
    /** The class of a genome whose lineage offers no taxon at or below the species rank. */
    public static final String UNPLACED = "*";

    private final Map<Integer, Map<String, Integer>> countsByPos = new HashMap<>();
    private final Map<Integer, String> majorityByPos = new HashMap<>();
    private final Map<Integer, Integer> candidatesByPos = new HashMap<>();
    private final Map<Integer, SmallTaxTree.SmallTaxIdNode> nodeByPos = new HashMap<>();
    private final Map<String, Integer> collectionCounts = new LinkedHashMap<>();
    private final List<String> classes;
    private final SmallTaxTree tree;
    private int collectionTotal;

    /**
     * Counts the species below every node of the given tree.
     *
     * @param tree the taxonomy of the database being classified against
     */
    public TaxonComposition(SmallTaxTree tree) {
        this.tree = tree;
        Deque<SmallTaxTree.SmallTaxIdNode> order = new ArrayDeque<>();
        Deque<SmallTaxTree.SmallTaxIdNode> stack = new ArrayDeque<>();
        for (SmallTaxTree.SmallTaxIdNode node : tree) {
            if (node.getParent() == null) {
                stack.push(node);
            }
        }
        while (!stack.isEmpty()) {
            SmallTaxTree.SmallTaxIdNode node = stack.pop();
            order.push(node);
            nodeByPos.put(node.getPosition(), node);
            SmallTaxTree.SmallTaxIdNode[] subNodes = node.getSubNodes();
            if (subNodes != null) {
                for (SmallTaxTree.SmallTaxIdNode sub : subNodes) {
                    stack.push(sub);
                }
            }
        }
        while (!order.isEmpty()) {
            SmallTaxTree.SmallTaxIdNode node = order.pop();
            Map<String, Integer> own = new LinkedHashMap<>();
            SmallTaxTree.SmallTaxIdNode[] subNodes = node.getSubNodes();
            if (subNodes == null || subNodes.length == 0) {
                String leafClass = classOfLeaf(node);
                if (leafClass != null) {
                    own.put(leafClass, 1);
                }
            } else {
                for (SmallTaxTree.SmallTaxIdNode sub : subNodes) {
                    Map<String, Integer> subCounts = countsByPos.get(sub.getPosition());
                    if (subCounts != null) {
                        for (Map.Entry<String, Integer> e : subCounts.entrySet()) {
                            own.merge(e.getKey(), e.getValue(), Integer::sum);
                        }
                    }
                }
            }
            if (own.isEmpty()) {
                continue;
            }
            countsByPos.put(node.getPosition(), own);
            candidatesByPos.put(node.getPosition(), own.size());
            String majority = null;
            int best = -1;
            for (Map.Entry<String, Integer> e : own.entrySet()) {
                if (e.getValue() > best) {
                    best = e.getValue();
                    majority = e.getKey();
                }
            }
            majorityByPos.put(node.getPosition(), majority);
        }
        for (SmallTaxTree.SmallTaxIdNode node : tree) {
            if (node.getParent() == null) {
                Map<String, Integer> rootCounts = countsByPos.get(node.getPosition());
                if (rootCounts != null) {
                    for (Map.Entry<String, Integer> e : rootCounts.entrySet()) {
                        collectionCounts.merge(e.getKey(), e.getValue(), Integer::sum);
                        collectionTotal += e.getValue();
                    }
                }
            }
        }
        List<String> l = new ArrayList<>(collectionCounts.keySet());
        Collections.sort(l);
        classes = Collections.unmodifiableList(l);
    }

    /**
     * Resolves the species a leaf belongs to: the <em>highest</em> ancestor, itself included, whose
     * rank is species or below it.
     * <p>
     * The highest and not the first one met, which is what this did before and what made every one of
     * the 153 strains of {@code Streptococcus pneumoniae} a class of its own. That is the wrong unit
     * here twice over. The candidate precision of the paper is {@code 1/|Sigma(k)|} over the
     * <em>species</em> of a node's subtree -- "an assignment to a genus comprising 40 species scores
     * 1/40" -- so counting strains inflates the denominator; and the reference standard a call is
     * held against names a species, so a candidate set of strain tax ids could never contain it. A
     * correct call on {@code Streptococcus pyogenes} scored zero for that reason alone.
     * <p>
     * Ranks between {@link Rank#SPECIES} and {@link Rank#FORMA_SPECIALIS} are the ones that count as
     * "species or below"; {@link Rank#STRAIN} sits among them, which is why the walk cannot stop at
     * the first hit. Everything else is stepped over rather than treated as an end: the artificial
     * DATA, FILE, ID and REFINED nodes, {@link Rank#ISOLATE}, the {@code no rank} and {@code clade}
     * buckets such as {@code unclassified Streptococcus}, and a rank the enum does not know at all.
     * None of them says anything about the species. The walk ends at the first real rank above the
     * species -- a genus, a family -- and answers with whatever it last remembered.
     *
     * @param leaf the leaf to resolve
     * @return the tax id of its species, or {@link #UNPLACED} if its lineage names none
     */
    /**
     * The class a leaf contributes to the candidate sets above it, or {@code null} where it
     * contributes none.
     * <p>
     * A childless node of rank {@link Rank#REFINED} is an OTHER bucket: the refinement creates it to
     * stand for "none of the siblings", it holds no genome, and in every database of this paper no
     * k-mer either. It is not an organism a k-mer could have come from, so it must not enter the
     * candidate set. Counting it added exactly one to every refined node above it and to no unrefined
     * node at all, since the unrefined tree has none: the genus of {@code strepto} reported 158
     * species against the unrefined 157, and the node holding the mitis complex reported nine against
     * eight. The measure understated the refinement it was there to assess.
     * <p>
     * {@link #UNPLACED} stays a candidate where a real leaf produces it. A genome under
     * {@code unclassified Streptococcus} names no species, and a k-mer there genuinely could have
     * come from an organism the taxonomy cannot place; that is a fact about the reference and not an
     * artefact of the refinement.
     *
     * @param leaf a node without sub-nodes
     * @return its species, {@link #UNPLACED}, or {@code null} for an OTHER bucket
     */
    static String classOfLeaf(SmallTaxTree.SmallTaxIdNode leaf) {
        return leaf.getRank() == Rank.REFINED ? null : speciesOf(leaf);
    }

    static String speciesOf(SmallTaxTree.SmallTaxIdNode leaf) {
        String species = UNPLACED;
        for (SmallTaxTree.SmallTaxIdNode n = leaf; n != null; n = n.getParent()) {
            int r = n.getRankOrdinal();
            if (r >= Rank.SPECIES.ordinal() && r <= Rank.FORMA_SPECIALIS.ordinal()) {
                species = n.getTaxId();
            } else if (r >= 0 && r < Rank.SPECIES.ordinal()) {
                // A rank the enum knows and that sits above the species. Note the r >= 0: an unknown
                // rank is -1, which would otherwise end the walk here and lose the species above.
                break;
            }
        }
        return species;
    }

    /** @param nodePos dense node position
     *  @return the species most genomes below that node belong to, or {@code null} */
    public String getMajorityClass(int nodePos) {
        return majorityByPos.get(nodePos);
    }

    /** @param nodePos dense node position
     *  @return how many distinct species sit below that node */
    public int getCandidates(int nodePos) {
        Integer c = candidatesByPos.get(nodePos);
        return c == null ? 0 : c;
    }

    /** @param nodePos dense node position
     *  @return the species below that node */
    public Set<String> getClassesAt(int nodePos) {
        Map<String, Integer> m = countsByPos.get(nodePos);
        return m == null ? Collections.<String>emptySet() : Collections.unmodifiableSet(m.keySet());
    }

    /** @param nodePos dense node position
     *  @param taxon a species
     *  @return how many genomes of that taxon sit below the node */
    public int getCount(int nodePos, String taxon) {
        Map<String, Integer> m = countsByPos.get(nodePos);
        if (m == null) {
            return 0;
        }
        Integer c = m.get(taxon);
        return c == null ? 0 : c;
    }

    /** @param nodePos dense node position
     *  @return how many genomes sit below the node in total */
    public int getCount(int nodePos) {
        Map<String, Integer> m = countsByPos.get(nodePos);
        if (m == null) {
            return 0;
        }
        int sum = 0;
        for (int v : m.values()) {
            sum += v;
        }
        return sum;
    }

    /** @param taxon a species
     *  @return how many genomes of it the whole database holds */
    public int getCollectionCount(String taxon) {
        Integer c = collectionCounts.get(taxon);
        return c == null ? 0 : c;
    }

    /** @return how many genomes the whole database holds */
    public int getCollectionTotal() {
        return collectionTotal;
    }

    /** @return every species the database holds a genome of, sorted */
    public List<String> getClasses() {
        return classes;
    }

    /** @param nodePos dense node position
     *  @return the node at that position, or {@code null} if the tree has none */
    public SmallTaxTree.SmallTaxIdNode getNode(int nodePos) {
        return nodeByPos.get(nodePos);
    }

    /** @return the taxonomy this composition was counted over */
    public SmallTaxTree getTree() {
        return tree;
    }
}
