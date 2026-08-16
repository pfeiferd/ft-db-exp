package org.metagene.ftdbexp.stquality;

import org.metagene.genestrip.tax.SmallTaxTree;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Every node of a database as a predictor of sequence type, and the naive Bayes classifier that
 * follows from taking the reads of one isolate as independent observations of such nodes.
 * <p>
 * The node a read is classified to is what is observed; the isolate's class is what is to be
 * inferred. The classes are the sequence types of the collection together with {@link #STAR},
 * for the genomes the scheme does not cover. Writing {@code P(s)} for the share of the database's
 * genomes in class {@code s} and {@code P(s | v)} for that share among the genomes below node
 * {@code v}, taking the reads as independent observations given the class gives
 * <pre>
 *   P(s | v_1..v_n)  ~  P(s) * PRODUCT_i  P(s | v_i) / P(s)
 * </pre>
 * up to a factor common to every class, which cannot affect which one wins. Evaluated as a sum of
 * logarithms, since a product over millions of reads underflows long before it is compared.
 * <p>
 * The form is worth reading twice: a read classified to a node whose class distribution is the
 * collection's own contributes a factor of one. It is not evidence, and it is not counted as any.
 * Only a node that shifts the distribution away from the prior moves the estimate, and it moves it
 * the further the rarer the class it favours -- which is what a majority vote over reads gets wrong
 * when most reads rest on a node that has narrowed nothing.
 * <p>
 * Both probabilities are Laplace-smoothed over the classes. Without it a single read classified to a
 * node whose subtree happens to hold no genome of the true class would send that class's estimate to
 * zero and no number of later reads could recover it, which is not a statement the data supports: a
 * read may come from sequence the isolate shares with another lineage.
 * <p>
 * The independence assumption is false here and is worth saying so. Reads of one isolate come from
 * one genome, so they are heavily correlated, and the resulting scores are far too confident to be
 * read as probabilities. What survives the assumption is the ordering, which is all the classifier
 * uses: the type with the highest score is the one the reads point at, and that is the decision
 * being scored.
 */
public class NaiveBayesSTModel {
    /**
     * The class of a genome the typing scheme does not cover: the paper's star in
     * {@code S' = S union {star}}.
     * <p>
     * Written {@code *} and not {@code -}, which is what {@code mlst} writes in its input and what
     * {@link STGroundTruth} recognises there. The two are deliberately different symbols: one is an
     * absence in a data file, the other a class the model predicts and may call an isolate for, and
     * a reader of the output should be able to tell a lineage the classifier settled on from a
     * field nobody filled in.
     * <p>
     * A genome without a sequence type is not a genome without a lineage. The seven-locus scheme
     * fails on the cryptic clades of C. difficile outright, and a novel allele combination
     * yields no type until PubMLST issues one; a broken assembly yields none either, for a reason
     * that has nothing to do with the organism. Left out of the model, all of them would silently
     * inflate the confidence of whatever types remain -- ten typed and ten untyped genomes below a
     * node would read as {@code P(s | v) = 1} where half the evidence says nothing of the sort.
     * Carried as a class, the estimate is diluted correctly, and the classifier may call an isolate
     * untyped, which for this organism is the informative answer that it looks like nothing the
     * scheme covers.
     * <p>
     * The two reasons for being untyped need not be told apart beforehand. Where untyped genomes
     * scatter, as a broken assembly does, {@code P(*|v)} stays at the prior and the class contributes
     * a factor of one wherever it appears; where they cluster, as a cryptic clade does, it
     * contributes exactly as much as the clustering warrants.
     */
    public static final String STAR = "*";

    /** The classes, i.e. the types of the collection followed by {@link #STAR}. */
    private final List<String> types;
    /** Per node position, the number of genomes of each type below it; {@code null} where none. */
    private final Map<Integer, Map<String, Integer>> countsByPos = new HashMap<>();
    /** Per node position, the type most of its genomes carry. */
    private final Map<Integer, String> majorityByPos = new HashMap<>();
    /** Per node position, how many types are in question there. */
    private final Map<Integer, Integer> candidatesByPos = new HashMap<>();
    /** Genomes of each type in the whole collection, i.e. what the prior is taken over. */
    private final Map<String, Integer> priorCounts = new LinkedHashMap<>();
    private int priorTotal;

    /**
     * Builds the model from a database's taxonomy and the typing of the genomes it was filled from.
     *
     * @param tree        the database's taxonomy
     * @param genomeTypes the type of every genome the database was filled from
     */
    public NaiveBayesSTModel(SmallTaxTree tree, STGroundTruth genomeTypes) {
        List<String> classes = new ArrayList<>(genomeTypes.getSTNames());
        classes.add(STAR);
        this.types = classes;

        // One post-order pass: a node's counts are the sum of its children's, a leaf contributes its
        // own genome. Iterative because a refined subtree is routinely hundreds of levels deep, which
        // is the very shape this model is built to evaluate.
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
                // from the denominator rather than evidence absent.
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
            // The root of the requested subtree carries every genome exactly once, so the prior is
            // read off the largest node rather than counted separately.
            int total = 0;
            for (int v : own.values()) {
                total += v;
            }
            if (total > priorTotal) {
                priorTotal = total;
                priorCounts.clear();
                priorCounts.putAll(own);
            }
        }
    }

    /**
     * Returns the type most genomes below the given node carry, i.e. what that node predicts on its
     * own.
     *
     * @param nodePos the node's dense position
     * @return the majority type, or {@code null} if no typed genome sits below the node
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
     * Returns the smoothed {@code log ( P(s | v) / P(s) )} that one read classified to the given node
     * contributes to the type's score.
     *
     * @param nodePos the node's dense position
     * @param type    the sequence type
     * @return the log ratio, or 0 where the node carries no typed genome and hence no information
     */
    public double getLogRatio(int nodePos, String type) {
        Map<String, Integer> counts = countsByPos.get(nodePos);
        if (counts == null) {
            return 0;
        }
        int total = 0;
        for (int v : counts.values()) {
            total += v;
        }
        double atNode = (counts.getOrDefault(type, 0) + 1.0) / (total + types.size());
        double prior = (priorCounts.getOrDefault(type, 0) + 1.0) / (priorTotal + types.size());
        return Math.log(atNode / prior);
    }

    /**
     * Classifies an isolate from the nodes its reads were classified to.
     *
     * @param readsPerNode how many reads landed at each node, by dense node position
     * @return the type with the highest score, or {@code null} if no read carried information
     */
    public String classify(Map<Integer, Long> readsPerNode) {
        if (readsPerNode.isEmpty()) {
            return null;
        }
        String best = null;
        double bestScore = Double.NEGATIVE_INFINITY;
        for (String type : types) {
            // Uniform over the classes, hence a constant and hence omitted: it cannot change which
            // class wins. The prior belonging here is the one over the ISOLATE's class, and the
            // reference collection is the wrong population to read it from -- it records what has
            // been deposited, not what walks into a ward. In this collection ST 1 holds 71 genomes
            // and none of the 37 isolates is ST 1, while ST 11 is the most frequent among them.
            //
            // Not to be confused with the prior in getLogRatio's denominator, which is a different
            // quantity and stays: that one is the prior already contained in the node's own estimate
            // P(s | v), and dividing it out is what turns a posterior back into a likelihood ratio.
            // Removing it as well would multiply the collection's composition in n times over
            // instead of cancelling it.
            double score = 0;
            for (Map.Entry<Integer, Long> e : readsPerNode.entrySet()) {
                score += e.getValue() * getLogRatio(e.getKey(), type);
            }
            if (score > bestScore) {
                bestScore = score;
                best = type;
            }
        }
        return best;
    }

    /**
     * Returns the number of genomes the prior is taken over.
     *
     * @return the number of typed genomes in the collection
     */
    public int getPriorTotal() {
        return priorTotal;
    }

    /**
     * Returns the types the model can predict.
     *
     * @return the sequence types of the collection
     */
    public List<String> getTypes() {
        return types;
    }
}
