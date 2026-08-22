package org.metagene.ftdbexp.taxquality;

import java.util.Collections;
import java.util.Map;
import java.util.Set;

/**
 * A rule that names the data taxon a sample's reads point at, from what they matched.
 * <p>
 * This is the C.~difficile study's {@code STClassifier} with its class concept replaced. There the
 * unit was an isolate -- one genome, one sequence type -- and the question was which lineage it
 * belonged to. A clinical metagenome is not one genome, so the unit here is the reads a sample
 * contributes <em>below one node</em>, and the question is which species they point at. For the
 * genus Streptococcus that is exactly the question the case study asks: whether what a respiratory
 * specimen carries is the pneumococcus or one of its commensal neighbours.
 * <p>
 * Every rule is handed the same thing -- how much the sample contributed at each node, by dense node
 * position -- and answers with a tax id or with nothing. What it is handed may be counted in reads,
 * in matched {@code k}-mers or in distinct matched {@code k}-mers; the currency is the caller's
 * choice and not the rule's, so the same rule is scored in all three and the columns stay
 * comparable.
 * <p>
 * The type exists because the experiment is a comparison of rules, and a comparison is only worth
 * something if the things compared are interchangeable. Two consequences follow and both are the
 * point. Every rule weighs the same {@link TaxonComposition}, counted once, so what separates two
 * answers is the rule and never two readings of the database. And a new rule is a new subclass
 * rather than an edit to whatever writes the results.
 *
 * @see NaiveBayesTaxonModel
 * @see PathVoteTaxonModel
 */
public abstract class TaxonClassifier {
    /** How many genomes of each data taxon sit below each node, shared by every classifier. */
    protected final TaxonComposition composition;

    /**
     * Creates a classifier over an already counted composition.
     *
     * @param composition the data-taxon composition of the database's nodes
     */
    protected TaxonClassifier(TaxonComposition composition) {
        this.composition = composition;
    }

    /**
     * Names the data taxon one sample's reads point at.
     *
     * @param countsPerNode how much the sample contributed at each node, by dense node position
     * @return the tax id, or {@code null} where the rule makes no call
     */
    public abstract String classify(Map<Integer, Long> countsPerNode);

    /**
     * @return the name of this rule, as it appears in the result columns
     */
    public abstract String getName();

    /**
     * The rule's answer as a set, so that a rule which leaves several taxa open can be scored beside
     * one that always commits. A rule answering with $n$ taxa is credited $1/n$ by
     * {@link #getPrecision}, which is the candidate-weighting of Section
     * \enquote{Measuring classification quality} applied to whole samples instead of reads.
     *
     * @param countsPerNode how much the sample contributed at each node
     * @return the taxa the rule leaves standing, empty where it makes no call
     */
    public Set<String> classifyAll(Map<Integer, Long> countsPerNode) {
        String one = classify(countsPerNode);
        return one == null ? Collections.<String>emptySet() : Collections.singleton(one);
    }

    /**
     * @param countsPerNode how much the sample contributed at each node
     * @param trueTaxon     the taxon the reference standard names for this sample
     * @return $1/n$ if the rule's $n$ answers include the true one, else zero
     */
    public double getPrecision(Map<Integer, Long> countsPerNode, String trueTaxon) {
        Set<String> answer = classifyAll(countsPerNode);
        return answer.contains(trueTaxon) ? 1.0 / answer.size() : 0.0;
    }

    /**
     * @param countsPerNode how much the sample contributed at each node
     * @return how many taxa the rule leaves standing
     */
    public int getCandidates(Map<Integer, Long> countsPerNode) {
        return classifyAll(countsPerNode).size();
    }

    /** @return the composition every rule of this comparison weighs */
    public TaxonComposition getComposition() {
        return composition;
    }
}
