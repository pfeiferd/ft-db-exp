package org.metagene.ftdbexp.stquality;

import java.util.Collections;
import java.util.Map;
import java.util.Set;

/**
 * A rule that names the lineage of one isolate from what its reads matched.
 * <p>
 * Every classifier here is handed the same thing -- how much the isolate contributed at each node of
 * the database, by dense node position -- and answers with a sequence type or with nothing. What it
 * is handed may be counted in reads, in matched {@code k}-mers or in distinct matched {@code k}-mers;
 * the currency is the caller's choice and not the rule's, so the same rule is scored in all three and
 * the columns stay comparable.
 * <p>
 * The type exists because the experiment is a comparison of rules, and a comparison is only worth
 * something if the things compared are interchangeable. Two consequences follow and both are the
 * point. Every rule weighs the same {@link STComposition}, counted once, so what separates two
 * answers is the rule and never two readings of the database. And a new rule is a new subclass rather
 * than an edit to whatever writes the results, which is what the last few variants each cost.
 *
 * @see NaiveBayesSTModel
 * @see PathVoteSTModel
 */
public abstract class STClassifier {
    /** How many genomes of each type sit below each node, shared by every classifier. */
    protected final STComposition composition;

    /**
     * Creates a classifier over an already counted composition.
     *
     * @param composition the type composition of the database's nodes
     */
    protected STClassifier(STComposition composition) {
        this.composition = composition;
    }

    /**
     * Names the lineage of one isolate.
     *
     * @param countsPerNode how much the isolate contributed at each node, by dense node position
     * @return the sequence type, or {@code null} where the rule makes no call
     */
    public abstract String classify(Map<Integer, Long> countsPerNode);

    /**
     * Returns a short name for this rule, used to label its column in the results.
     *
     * @return the name
     */
    public abstract String getName();

    /**
     * Names every lineage the rule leaves in question for one isolate.
     * <p>
     * An answer need not be a single type. A rule that places the isolate at a node of the taxonomy
     * answers with the types below that node, and where several remain the isolate is classified to
     * all of them at once -- which is a weaker statement than naming one, and is scored as one:
     * {@link #getPrecision} gives such an answer the credit {@code 1 / n} that the paper's
     * {@code p_st} gives a k-mer resting at a node leaving {@code n} types open. By default a rule
     * that names one type answers with just that one.
     *
     * @param countsPerNode how much the isolate contributed at each node, by dense node position
     * @return the types left in question, empty where the rule makes no call
     */
    public Set<String> classifyAll(Map<Integer, Long> countsPerNode) {
        String one = classify(countsPerNode);
        return one == null ? Collections.<String>emptySet() : Collections.singleton(one);
    }

    /**
     * Scores one isolate: the credit the rule earns for it.
     * <p>
     * One where the rule names the isolate's lineage and nothing else, {@code 1 / n} where it leaves
     * {@code n} types in question and the true one is among them, and zero where it is not or where
     * no call is made. Naming every type of the database therefore earns almost nothing rather than
     * counting as right, and naming one wrongly earns nothing rather than being excused.
     *
     * @param countsPerNode how much the isolate contributed at each node, by dense node position
     * @param trueType      the isolate's actual sequence type
     * @return the credit, between 0 and 1
     */
    public double getPrecision(Map<Integer, Long> countsPerNode, String trueType) {
        Set<String> answer = classifyAll(countsPerNode);
        return answer.contains(trueType) ? 1.0 / answer.size() : 0.0;
    }

    /**
     * Returns how many types are still in question where this rule placed the isolate.
     *
     * @param countsPerNode how much the isolate contributed at each node, by dense node position
     * @return the number of candidate types, or 0 where no call is made
     */
    public int getCandidates(Map<Integer, Long> countsPerNode) {
        return classifyAll(countsPerNode).size();
    }

    /**
     * Returns the composition this classifier weighs.
     *
     * @return the composition
     */
    public STComposition getComposition() {
        return composition;
    }
}
