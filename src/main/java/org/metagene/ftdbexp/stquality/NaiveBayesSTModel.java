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
 * inferred. The classes are the sequence types the collection holds a genome of, together with
 * {@link #STAR} where some genome carries none. Taking the reads as independent observations given
 * the class gives
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
 * Each factor is evaluated as the node's LIFT for the class, {@code P(s|v) / P(s)}, which Bayes makes
 * equal to {@code P(v|s) / P(v)} and hence to {@code P(v|s)} up to the class-independent {@code P(v)}
 * that the comparison drops anyway. Writing {@code g(s|v)} for the genomes of class {@code s} below
 * node {@code v} and {@code g(v)} for all of them,
 * <pre>
 *   P(s|v) = ( g(s|v) + a ) / ( g(v) + 1 )        P(s) = ( g(s) + a ) / ( g + 1 )
 * </pre>
 * with {@code a = SMOOTHING / |S'|}, the collection standing in for the all-covering node.
 * <p>
 * That arrangement is not a matter of taste, and getting it wrong is what made an earlier version of
 * this class call every isolate ST 1. Nearly every read of a sample lands at the node holding the
 * shared {@code k}-mers of the species -- in cdiff between 97 and 99 percent of them -- and that node
 * covers every genome of every class. Its lift is then exactly one: numerator and denominator are the
 * same expression over the same counts, so the quotient is the bit pattern {@code 1.0} and its
 * logarithm is {@code 0.0}, not merely something small. Millions of uninformative reads therefore
 * contribute nothing at all, and the decision is left to the few that landed lower. Estimate the same
 * quantity as a distribution over the nodes instead, and the all-covering node keeps a residue that
 * differs between classes -- by class size if the smoothing is coarse, by the mean depth of the
 * class's genomes even if it is not. Multiplied by a few million reads, a residue of a tenth of a nat
 * is half a million nats, and whichever class it favours wins every isolate regardless of the
 * evidence. A property that holds only approximately is worth nothing at this sample size.
 * <p>
 * A class with no genome in the collection is left out rather than smoothed, see the constructor: it
 * has no distribution to estimate, and any convention chosen for it decides instead of describing.
 * <p>
 * The form of the product is worth reading twice. A node whose genomes carry the collection's own
 * class distribution contributes a factor of one to every class: it is not evidence, and it is not
 * counted as any. Only a node that concentrates one class more than the collection does moves the
 * estimate, and it moves it the further the rarer the class it favours -- which is what a majority
 * vote over reads gets wrong when most reads rest on a node that has narrowed nothing.
 * <p>
 * The independence assumption is false here and is worth saying so. Reads of one isolate come from
 * one genome, so they are heavily correlated, and the resulting scores are far too confident to be
 * read as probabilities. What survives the assumption is the ordering, which is all the classifier
 * uses: the type with the highest score is the one the reads point at, and that is the decision
 * being scored.
 */
public class NaiveBayesSTModel extends STClassifier {
    /**
     * The class of a genome the typing scheme does not cover: the paper's star in
     * {@code S' = S union {star}}. The same constant as {@link STComposition#STAR}, named here as
     * well because the model's own documentation and the paper both speak of it.
     */
    public static final String STAR = STComposition.STAR;

    /**
     * One pseudo-genome, spread evenly over the classes, added to every count before a share is taken
     * of it -- Lidstone smoothing with {@code SMOOTHING / |S'|} per class rather than the textbook one
     * per class.
     * <p>
     * The size matters and one per class is far too much. A node holding 35 genomes would then receive
     * 58 imaginary ones, which outweigh what is actually there: a class the node does not hold at all
     * comes out above its own baseline, so a read landing there is counted as evidence FOR a lineage
     * that is absent. With one pseudo-genome in total the same read tells against it by 1.84 nats, and
     * with a hundredth of one by 4.61 -- the latter being more confidence than this data supports,
     * since a single read misplaced by a false positive of the k-mer index would then all but rule the
     * true lineage out. What the value cannot affect is the neutrality of an all-covering node, which
     * holds for every {@code SMOOTHING} because the baseline is that node's own estimate.
     */
    private static final double SMOOTHING = 1.0;

    /**
     * Builds the model from a database's taxonomy and the typing of the genomes it was filled from.
     *
     * @param tree        the database's taxonomy
     * @param genomeTypes the type of every genome the database was filled from
     */
    public NaiveBayesSTModel(SmallTaxTree tree, STGroundTruth genomeTypes) {
        this(new STComposition(tree, genomeTypes));
    }

    /**
     * Builds the model over a composition that has already been counted, so that several classifiers
     * scored against each other weigh the same database rather than two readings of it.
     *
     * @param composition the type composition of the database's nodes
     */
    public NaiveBayesSTModel(STComposition composition) {
        super(composition);
    }

    @Override
    public String getName() {
        return "bayes";
    }

    /**
     * Returns the type most genomes below the given node carry, i.e. what that node predicts on its
     * own.
     *
     * @param nodePos the node's dense position
     * @return the majority type, or {@code null} if no typed genome sits below the node
     */
    public String getMajorityST(int nodePos) {
        return composition.getMajorityST(nodePos);
    }

    /**
     * Returns how many types are still in question at the given node.
     *
     * @param nodePos the node's dense position
     * @return the number of candidate types, or 0 if none
     */
    public int getCandidates(int nodePos) {
        return composition.getCandidates(nodePos);
    }

    /**
     * Returns the smoothed {@code log P(v | s)} that one read classified to the given node
     * contributes to the type's score.
     *
     * @param nodePos the node's dense position
     * @param type    the sequence type
     * @return the log likelihood, always negative since it is a probability over the nodes
     */
    public double getLogLift(int nodePos, String type) {
        int below = composition.getCount(nodePos);
        if (below == 0) {
            return 0;
        }
        double alpha = SMOOTHING / composition.getTypes().size();
        double atNode = (composition.getCount(nodePos, type) + alpha) / (below + SMOOTHING);
        double inCollection = (composition.getCollectionCount(type) + alpha)
                / (composition.getCollectionTotal() + SMOOTHING);
        return Math.log(atNode / inCollection);
    }

    /**
     * Classifies an isolate from the nodes its reads were classified to.
     *
     * @param readsPerNode how many reads landed at each node, by dense node position
     * @return the type with the highest score, or {@code null} if no read carried information
     */
    @Override
    public String classify(Map<Integer, Long> readsPerNode) {
        if (readsPerNode.isEmpty()) {
            return null;
        }
        String best = null;
        double bestScore = Double.NEGATIVE_INFINITY;
        for (String type : composition.getTypes()) {
            // The leading P(s), taken uniform over the classes and hence a constant that is omitted:
            // it cannot change which class wins. The prior belonging here is the one over the
            // ISOLATE's class, and the reference collection is the wrong population to read it from
            // -- it records what has been deposited, not what walks into a ward. In this collection
            // ST 1 holds 71 genomes and none of the 37 isolates is ST 1, while ST 11 is the most
            // frequent among them. It is the model's only prior; the likelihood's denominator is not
            // a second one but the normaliser of P(v | s) over the nodes.
            double score = 0;
            for (Map.Entry<Integer, Long> e : readsPerNode.entrySet()) {
                score += e.getValue() * getLogLift(e.getKey(), type);
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
        return composition.getCollectionTotal();
    }

    /**
     * Returns the types the model can predict.
     *
     * @return the sequence types of the collection
     */
    public List<String> getTypes() {
        return composition.getTypes();
    }
}
