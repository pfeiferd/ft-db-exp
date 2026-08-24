package org.metagene.ftdbexp.taxquality;

import java.util.Map;

/**
 * Names a sample's species by the log-lift its reads give each candidate.
 * <p>
 * For every candidate taxon $t$ the rule sums, over the nodes the sample contributed at, the
 * contribution times $\log P(t \mid v) / P(t)$: how much more of node $v$'s genomes belong to $t$
 * than of the database as a whole. A read landing at a node most of whose genomes are $t$ speaks
 * for $t$; a read landing where $t$ is no better represented than in the collection at large speaks
 * for nothing, and contributes zero. The taxon with the highest sum wins.
 * <p>
 * The leading $P(t)$ is taken uniform over the classes and hence omitted -- it cannot change which
 * class wins. That is deliberate and it was deliberate in the C.~difficile version too: the prior
 * belonging there is the one over the <em>sample's</em> class, and a reference collection is the
 * wrong population to read it from, since it records what has been deposited rather than what walks
 * into a ward. RefSeq holds 9,263 \emph{S.~pneumoniae} genomes against 242 of \emph{S.~mitis}, which
 * says a great deal about sequencing effort and nothing about what is in a given sputum sample.
 * <p>
 * <b>Not reported in the paper.</b> The C.~difficile study carried this rule for a while and dropped
 * it from its results before they were written up -- the surviving
 * {@code cdiff_isolate_st_db.csv} holds path-vote and majority columns and no naive-Bayes column at
 * all. It is kept here because the comparison of rules is what {@link TaxonClassifier} exists for and
 * a rule that is absent cannot be compared against, but {@link PathVoteTaxonModel} is what the case
 * study reports: it performed better there, and it takes two sentences to explain in a paper where
 * this one takes a paragraph and a prior nobody can justify from a reference collection.
 */
public class NaiveBayesTaxonModel extends TaxonClassifier {
    /** Laplace mass spread over the classes, so that an unseen class is unlikely and not impossible. */
    private static final double SMOOTHING = 1.0;

    /**
     * @param composition the data-taxon composition of the database's nodes
     */
    public NaiveBayesTaxonModel(TaxonComposition composition) {
        super(composition);
    }

    @Override
    public String classify(Map<Integer, Long> countsPerNode) {
        if (countsPerNode.isEmpty()) {
            return null;
        }
        String best = null;
        double bestScore = Double.NEGATIVE_INFINITY;
        for (String taxon : composition.getClasses()) {
            double score = 0;
            for (Map.Entry<Integer, Long> e : countsPerNode.entrySet()) {
                score += e.getValue() * getLogLift(e.getKey(), taxon);
            }
            if (score > bestScore) {
                bestScore = score;
                best = taxon;
            }
        }
        return best;
    }

    /**
     * How much more of a node's genomes belong to a taxon than of the collection at large, in logs.
     *
     * @param nodePos dense node position
     * @param taxon   the candidate species
     * @return the log lift, zero for a node with no genome below it
     */
    public double getLogLift(int nodePos, String taxon) {
        int below = composition.getCount(nodePos);
        if (below == 0) {
            return 0;
        }
        double alpha = SMOOTHING / composition.getClasses().size();
        double atNode = (composition.getCount(nodePos, taxon) + alpha) / (below + SMOOTHING);
        double inCollection = (composition.getCollectionCount(taxon) + alpha)
                / (composition.getCollectionTotal() + SMOOTHING);
        return Math.log(atNode / inCollection);
    }

    @Override
    public String getName() {
        return "naive bayes";
    }
}
