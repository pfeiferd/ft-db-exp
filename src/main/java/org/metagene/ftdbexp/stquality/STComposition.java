package org.metagene.ftdbexp.stquality;

import org.metagene.genestrip.tax.SmallTaxTree;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Which sequence types sit below each node of a database, counted once and shared by the classifiers
 * that need it.
 * <p>
 * Every classifier here answers the same question from a different angle -- which lineage do an
 * isolate's reads point at -- and every one of them needs the same two things: the type composition
 * of each node, and a way from a node's dense position back to the node. Counting that twice would
 * not merely be wasteful; it would make it possible for two classifiers to disagree about the
 * database rather than about the reads, which is the one difference their comparison must not
 * contain.
 * <p>
 * The class of a genome the typing scheme does not cover is {@link #STAR}. See there for why an
 * untyped genome is carried as a class of its own rather than dropped, and {@link #getTypes()} for
 * why a class the collection holds no genome of is not carried at all.
 */
public class STComposition {
    /**
     * The class of a genome the typing scheme does not cover: the paper's star in
     * {@code S' = S union {star}}.
     * <p>
     * Written {@code *} and not {@code -}, which is what {@code mlst} writes in its input and what
     * {@link STGroundTruth} recognises there. The two are deliberately different symbols: one is an
     * absence in a data file, the other a class a classifier may call an isolate for, and a reader of
     * the output should be able to tell a lineage that was settled on from a field nobody filled in.
     * <p>
     * A genome without a sequence type is not a genome without a lineage. The seven-locus scheme
     * fails on the cryptic clades of C. difficile outright, and a novel allele combination yields no
     * type until PubMLST issues one; a broken assembly yields none either, for a reason that has
     * nothing to do with the organism. Left out, all of them would silently inflate the confidence of
     * whatever types remain -- a node holding ten typed and ten untyped genomes would speak for the
     * typed ones as though nothing else sat below it.
     */
    public static final String STAR = "*";

    /** Per node position, the number of genomes of each type below it; absent where none. */
    private final Map<Integer, Map<String, Integer>> countsByPos = new HashMap<>();
    /** Per node position, the type most of its genomes carry. */
    private final Map<Integer, String> majorityByPos = new HashMap<>();
    /** Per node position, how many types are in question there. */
    private final Map<Integer, Integer> candidatesByPos = new HashMap<>();
    /** Per node position, the node, so that a classifier can walk from a count back up the tree. */
    private final Map<Integer, SmallTaxTree.SmallTaxIdNode> nodeByPos = new HashMap<>();
    /** Genomes of each type in the whole collection, i.e. {@code g(s)}. */
    private final Map<String, Integer> collectionCounts = new LinkedHashMap<>();
    /** Genomes in the whole collection, i.e. {@code g}. */
    private int collectionTotal;
    /** The classes, i.e. the types the collection holds a genome of. */
    private final List<String> types;
    private final SmallTaxTree tree;

    /**
     * Counts the composition of every node of the given database.
     *
     * @param tree        the database's taxonomy
     * @param genomeTypes the type of every genome the database was filled from
     * @throws IllegalStateException if no genome of the database carries a type the file knows of,
     *                               which means the two are keyed differently rather than that the
     *                               genomes are untyped
     */
    public STComposition(SmallTaxTree tree, STGroundTruth genomeTypes) {
        this.tree = tree;
        // One post-order pass: a node's counts are the sum of its children's, a leaf contributes its
        // own genome. Iterative because a refined subtree is routinely hundreds of levels deep, which
        // is the very shape these models are built to evaluate.
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
                // Every leaf counts, typed or not: a leaf left out would be evidence quietly removed
                // rather than evidence absent.
                int index = genomeTypes.getSTIndex(node.getName());
                own.put(index < 0 ? STAR : genomeTypes.getSTName(index), 1);
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
                    majority = e.getKey();
                    best = e.getValue();
                }
            }
            majorityByPos.put(node.getPosition(), majority);
            // The root of the requested subtree carries every genome exactly once, so g(s) is read
            // off the largest node rather than counted separately.
            int total = 0;
            for (int v : own.values()) {
                total += v;
            }
            if (total > collectionTotal) {
                collectionTotal = total;
                collectionCounts.clear();
                collectionCounts.putAll(own);
            }
        }

        // The classes are those the database actually holds a genome of, read off the collection and
        // not off the typing file. The file types every genome the wider collection knows of -- for
        // cdiff 196 sequence types against the 58 that were built into this database -- and a class
        // with no genome here is not one a model can weigh: it has no distribution to estimate, so
        // whatever convention were chosen for it would decide rather than describe. STAR is a class
        // here exactly when some genome of the database carries no type, and in cdiff none does.
        List<String> classes = new ArrayList<>();
        for (String name : genomeTypes.getSTNames()) {
            if (collectionCounts.containsKey(name)) {
                classes.add(name);
            }
        }
        if (collectionCounts.containsKey(STAR)) {
            classes.add(STAR);
        }
        if (classes.isEmpty()) {
            throw new IllegalStateException("The database holds no genome of any class the typing"
                    + " knows of, so there is nothing to predict. Most likely the typing keys its rows"
                    + " differently from the way the database names its file nodes, in which case every"
                    + " genome looks untyped and even the untyped class stays empty.");
        }
        this.types = classes;
    }

    /**
     * Returns the type most genomes below the given node carry, i.e. what that node predicts on its
     * own.
     *
     * @param nodePos the node's dense position
     * @return the majority type, or {@code null} if no genome sits below the node
     */
    public String getMajorityST(int nodePos) {
        return majorityByPos.get(nodePos);
    }

    /**
     * Returns how many types are still in question at the given node.
     *
     * @param nodePos the node's dense position
     * @return the number of candidate types, or 0 if none
     */
    public int getCandidates(int nodePos) {
        Integer c = candidatesByPos.get(nodePos);
        return c == null ? 0 : c;
    }

    /**
     * Returns the types below the given node, i.e. the ones still in question there.
     *
     * @param nodePos the node's dense position
     * @return the types, empty if no genome sits below the node
     */
    public Set<String> getTypesAt(int nodePos) {
        Map<String, Integer> counts = countsByPos.get(nodePos);
        return counts == null ? Collections.<String>emptySet()
                : Collections.unmodifiableSet(counts.keySet());
    }

    /**
     * Returns the number of genomes of the given type below the given node.
     *
     * @param nodePos the node's dense position
     * @param type    the sequence type
     * @return the number of genomes, 0 if none
     */
    public int getCount(int nodePos, String type) {
        Map<String, Integer> counts = countsByPos.get(nodePos);
        return counts == null ? 0 : counts.getOrDefault(type, 0);
    }

    /**
     * Returns the number of genomes below the given node, whatever their type.
     *
     * @param nodePos the node's dense position
     * @return the number of genomes, 0 if none
     */
    public int getCount(int nodePos) {
        Map<String, Integer> counts = countsByPos.get(nodePos);
        if (counts == null) {
            return 0;
        }
        int sum = 0;
        for (int v : counts.values()) {
            sum += v;
        }
        return sum;
    }

    /**
     * Returns the number of genomes of the given type in the whole collection, i.e. {@code g(s)}.
     *
     * @param type the sequence type
     * @return the number of genomes, 0 if none
     */
    public int getCollectionCount(String type) {
        return collectionCounts.getOrDefault(type, 0);
    }

    /**
     * Returns the number of genomes the collection holds, i.e. {@code g}.
     *
     * @return the number of genomes
     */
    public int getCollectionTotal() {
        return collectionTotal;
    }

    /**
     * Returns the classes a model built on this composition can predict.
     *
     * @return the sequence types, in the order the typing file names them, {@link #STAR} last
     */
    public List<String> getTypes() {
        return Collections.unmodifiableList(types);
    }

    /**
     * Returns the node at the given dense position.
     *
     * @param nodePos the position
     * @return the node, or {@code null} if the tree has none there
     */
    public SmallTaxTree.SmallTaxIdNode getNode(int nodePos) {
        return nodeByPos.get(nodePos);
    }

    /**
     * Returns the tree this composition was counted over.
     *
     * @return the taxonomy
     */
    public SmallTaxTree getTree() {
        return tree;
    }
}
