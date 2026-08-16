package org.metagene.ftdbexp.stquality;

import org.metagene.genestrip.tax.SmallTaxTree;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Every node of a database as a predictor of sequence type, and the naive Bayes classifier that
 * follows from taking the reads of one isolate as independent observations of such nodes.
 * <p>
 * The node a read is classified to is what is observed; the isolate's type is what is to be inferred.
 * Writing {@code P(s)} for the share of the database's genomes that are of type {@code s} and
 * {@code P(s | v)} for that share among the genomes below node {@code v}, Bayes on each read and the
 * naive assumption that the reads are independent given the type give
 * <pre>
 *   log P(s | v_1..v_n)  =  log P(s)  +  SUM_i log ( P(s | v_i) / P(s) )  +  const,
 * </pre>
 * the constant collecting the terms that do not depend on {@code s} and so do not affect the
 * decision. The classifier is therefore a sum of log ratios over the reads, and that form is worth
 * reading twice: a read classified to a node whose type distribution is the collection's own
 * contributes {@code log 1 = 0}. It is not evidence, and it is not counted as any. Only a node that
 * shifts the distribution away from the prior moves the score, and it moves it the further the rarer
 * the type it favours -- which is what a majority vote over reads gets wrong when most reads rest on
 * a node that has narrowed nothing.
 * <p>
 * Both probabilities are Laplace-smoothed over the types of the whole collection. Without it a
 * single read classified to a node whose subtree happens to hold no genome of the true type would
 * send that type's score to minus infinity and no number of later reads could recover it, which is
 * not a statement the data supports: a read may come from sequence the isolate shares with another
 * lineage.
 * <p>
 * The independence assumption is false here and is worth saying so. Reads of one isolate come from
 * one genome, so they are heavily correlated, and the resulting scores are far too confident to be
 * read as probabilities. What survives the assumption is the ordering, which is all the classifier
 * uses: the type with the highest score is the one the reads point at, and that is the decision
 * being scored.
 */
public class NaiveBayesSTModel {
    /** The types of the collection, in the order the ground truth assigns them. */
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
        this.types = genomeTypes.getSTNames();

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
                int index = genomeTypes.getSTIndex(node.getName());
                if (index >= 0) {
                    own.put(genomeTypes.getSTName(index), 1);
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
            double score = Math.log((priorCounts.getOrDefault(type, 0) + 1.0) / (priorTotal + types.size()));
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
