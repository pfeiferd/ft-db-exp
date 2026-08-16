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
 * inferred. The classes are the sequence types of the collection together with {@link #STAR}, for
 * the genomes the scheme does not cover. Taking the reads as independent observations given the
 * class gives
 * <pre>
 *   P(s | v_1..v_n)  ~  P(s) * PRODUCT_i  P(v_i | s)
 * </pre>
 * up to a factor common to every class -- {@code P(v_1..v_n)} above all -- which cannot affect which
 * one wins, and the isolate is called for the class that maximises the right-hand side. Evaluated as
 * a sum of logarithms, since a product over millions of reads underflows long before it is compared.
 * <p>
 * {@code P(s)} is the prior over the ISOLATE's class and is taken uniform; see {@link #classify},
 * where the reasons sit next to the line that drops it.
 * <p>
 * {@code P(v | s)} is what the database supplies. Writing {@code g(s|v)} for the number of genomes
 * of class {@code s} below node {@code v} and {@code V} for the nodes of the database, a read of an
 * isolate of class {@code s} is taken to land at {@code v} in proportion to how many of that class's
 * genomes the node covers:
 * <pre>
 *   P(v | s)  =  ( g(s|v) + 1 )  /  SUM_v' ( g(s|v') + 1 )
 * </pre>
 * The added one is Laplace smoothing, which keeps a class alive that a node does not happen to hold:
 * without it a single read classified to a node whose subtree holds no genome of the true class
 * would send that class to zero and no number of later reads could recover it, and a read may come
 * from sequence the isolate shares with another lineage. The denominator is then no more than what
 * makes the numerators sum to one over the nodes.
 * <p>
 * That normalisation is not a formality. Divide instead by {@code g(s) + 1}, the number of genomes
 * the class has in the collection, giving the intuitive "share of them below v", and a class with no
 * genome at all scores {@code 1} at every node, which is the largest value any class can reach
 * anywhere: it then wins or ties against every real class on every read, and every isolate is called
 * untyped. The equivalent ratio form {@code P(s|v)/P(s)} inherits the defect unchanged, since the two
 * differ only by a factor common to all classes. Normalised over the nodes such a class is spread
 * uniformly at {@code 1/|V|} instead, which says what it should: a class nothing has been seen of
 * explains no node in particular. In the cdiff database every genome carries a type, so {@link #STAR}
 * is exactly such a class and the distinction decides the whole experiment rather than an edge case.
 * <p>
 * The form of the product is worth reading twice. Take a node whose genomes carry the collection's
 * own class distribution, so that {@code g(s|v) = g(s) g(v) / g} for every class: up to the smoothing
 * the factor it contributes is {@code g(v)/g} divided by the number of nodes a genome of that class
 * lies under on average. It does not depend on how frequent the class is in the collection at all,
 * only on how deep its genomes sit, so for classes at a comparable depth such a node cannot move the
 * decision however many reads land there. It is not evidence, and it is not counted as any. Only a
 * node that concentrates one class more than the collection does moves the estimate, and it moves it
 * the further the rarer the class it favours -- which is what a majority vote over reads gets wrong
 * when most reads rest on a node that has narrowed nothing.
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
     * inflate the confidence of whatever types remain -- a node holding ten typed and ten untyped
     * genomes would speak for the typed ones as though nothing else sat below it. Carried as a class,
     * the untyped genomes compete for the node like any other, and the classifier may call an isolate
     * untyped, which for this organism is the informative answer that it looks like nothing the
     * scheme covers.
     * <p>
     * The two reasons for being untyped need not be told apart beforehand. Where untyped genomes
     * scatter, as a broken assembly does, the class is spread as thinly as the collection is and
     * cannot move a decision; where they cluster, as a cryptic clade does, it moves one exactly as
     * much as the clustering warrants.
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
    /** Genomes of each type in the whole collection, i.e. {@code g(s)}. */
    private final Map<String, Integer> collectionCounts = new LinkedHashMap<>();
    /** Genomes in the whole collection, i.e. {@code g}. */
    private int collectionTotal;
    /**
     * Per type, {@code SUM_v g(s|v)}: how often a genome of that type is counted over all nodes,
     * which is the sum over its genomes of the number of nodes each lies under. Together with the
     * node count it is the denominator of {@link #getLogLikelihood} -- the two add up to
     * {@code SUM_v ( g(s|v) + 1 )} -- and hence what normalises it over the nodes.
     */
    private final Map<String, Integer> incidencesByType = new HashMap<>();
    /** The number of nodes, i.e. {@code |V|}, the observation space a read is drawn from. */
    private int nodeCount;

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
            nodeCount++;
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
            for (Map.Entry<String, Integer> e : own.entrySet()) {
                incidencesByType.merge(e.getKey(), e.getValue(), Integer::sum);
            }
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
     * Returns the smoothed {@code log P(v | s)} that one read classified to the given node
     * contributes to the type's score.
     *
     * @param nodePos the node's dense position
     * @param type    the sequence type
     * @return the log likelihood, always negative since it is a probability over the nodes
     */
    public double getLogLikelihood(int nodePos, String type) {
        Map<String, Integer> counts = countsByPos.get(nodePos);
        // A node with no genome below it is not a special case: it gets the smoothed numerator like
        // any node the class does not reach. What it must not get is a factor of one, which is what
        // an unnormalised estimate hands the classes nothing has been seen of.
        int atNode = counts == null ? 0 : counts.getOrDefault(type, 0);
        return Math.log((atNode + 1.0) / (incidencesByType.getOrDefault(type, 0) + (double) nodeCount));
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
            // The leading P(s), taken uniform over the classes and hence a constant that is omitted:
            // it cannot change which class wins. The prior belonging here is the one over the
            // ISOLATE's class, and the reference collection is the wrong population to read it from
            // -- it records what has been deposited, not what walks into a ward. In this collection
            // ST 1 holds 71 genomes and none of the 37 isolates is ST 1, while ST 11 is the most
            // frequent among them. It is the model's only prior; the likelihood's denominator is not
            // a second one but the normaliser of P(v | s) over the nodes.
            double score = 0;
            for (Map.Entry<Integer, Long> e : readsPerNode.entrySet()) {
                score += e.getValue() * getLogLikelihood(e.getKey(), type);
            }
            if (score > bestScore) {
                bestScore = score;
                best = type;
            }
        }
        return best;
    }

    /**
     * Returns the number of genomes the model was built over, i.e. {@code g}.
     *
     * @return the number of genomes in the collection
     */
    public int getCollectionTotal() {
        return collectionTotal;
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
