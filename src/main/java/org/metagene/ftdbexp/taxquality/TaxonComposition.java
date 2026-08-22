package org.metagene.ftdbexp.taxquality;

import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.util.*;

/**
 * Which data taxa sit below each node of a database, counted once and shared by the classifiers that
 * need it.
 * <p>
 * This is the taxonomic counterpart of the {@code STComposition} the C.~difficile study used. There
 * the class of a genome was its seven-locus sequence type, which had to be established outside the
 * database and carried in from a typing file. Here the class is the genome's <em>data taxon</em> in
 * the sense of the paper -- the taxon its genome is filed at -- so the tree already knows it and
 * nothing has to be read in beside it. The rest is unchanged, and deliberately so: every classifier
 * weighs the same composition, counted once, so that what separates two answers is the rule and
 * never two readings of the database.
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
     * Counts the data taxa below every node of the given tree.
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
                own.put(dataTaxonOf(node), 1);
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
     * Resolves the data taxon a leaf belongs to: the first ancestor, itself included, whose rank is
     * species or below it. The artificial nodes the fill and the refinement insert -- DATA, FILE, ID
     * and REFINED -- have no taxonomic rank of their own and are skipped, which is what makes the
     * answer a taxon of the reference taxonomy rather than a node of this particular database.
     *
     * @param leaf the leaf to resolve
     * @return the tax id of its data taxon, or {@link #UNPLACED} if its lineage has none
     */
    private static String dataTaxonOf(SmallTaxTree.SmallTaxIdNode leaf) {
        for (SmallTaxTree.SmallTaxIdNode n = leaf; n != null; n = n.getParent()) {
            int r = n.getRankOrdinal();
            if (r >= Rank.SPECIES.ordinal() && r <= Rank.FORMA_SPECIALIS.ordinal()) {
                return n.getTaxId();
            }
        }
        return UNPLACED;
    }

    /** @param nodePos dense node position
     *  @return the data taxon most genomes below that node belong to, or {@code null} */
    public String getMajorityClass(int nodePos) {
        return majorityByPos.get(nodePos);
    }

    /** @param nodePos dense node position
     *  @return how many distinct data taxa sit below that node */
    public int getCandidates(int nodePos) {
        Integer c = candidatesByPos.get(nodePos);
        return c == null ? 0 : c;
    }

    /** @param nodePos dense node position
     *  @return the data taxa below that node */
    public Set<String> getClassesAt(int nodePos) {
        Map<String, Integer> m = countsByPos.get(nodePos);
        return m == null ? Collections.<String>emptySet() : Collections.unmodifiableSet(m.keySet());
    }

    /** @param nodePos dense node position
     *  @param taxon a data taxon
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

    /** @param taxon a data taxon
     *  @return how many genomes of it the whole database holds */
    public int getCollectionCount(String taxon) {
        Integer c = collectionCounts.get(taxon);
        return c == null ? 0 : c;
    }

    /** @return how many genomes the whole database holds */
    public int getCollectionTotal() {
        return collectionTotal;
    }

    /** @return every data taxon the database holds a genome of, sorted */
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
