package org.metagene.ftdbexp.eval;

import org.metagene.genestrip.tax.Rank;

/**
 * Read classification counts for one fastq file and one database variant, together with the
 * precision, recall and F1 derived from them.
 * <p>
 * A read contributes to the tally as follows. If its ground truth cannot be resolved, it counts as
 * {@link #getUnresolved() unresolved} and is ignored otherwise -- it can neither be right nor wrong
 * if we do not know what it is. Every other read counts towards {@link #getTotal() total}, and if
 * the classifier assigned it a taxon, towards {@link #getClassified() classified}. It counts as
 * correct for a rank if the lowest common ancestor of its true and its assigned taxon lies at that
 * rank or below, which is the usual reading of "classified correctly down to rank r".
 * <p>
 * Consequently, {@code correct <= classified <= total} does <em>not</em> hold in general:
 * {@code classified} may exceed {@code total} once reads are counted that the classifier assigned a
 * taxon to although their ground truth places them outside the scope under consideration.
 */
public final class AccuracyTally {
    private long classified;
    private long correctGenus;
    private long correctSpecies;
    private long correctStrain;
    private long unresolved;
    private long total;
    private double speciesCandidateScore;
    // The genus-only subset and the two sums taken over it. The subset is the same for both: what
    // separates the gated from the ungated measure is only which count is summed, q or q'. Keeping
    // one denominator is what makes the two comparable, and the difference between the numerators
    // is the whole of the distinction.
    // Two subsets, because the two measures answer questions of different kinds. The gated one is a
    // precision and is therefore taken over the positives -- the reads the unrefined database placed
    // correctly at a genus, which is governed by sigma(r). The ungated one has to be formable on a
    // sample that supplies no sigma(r) at all, so it is taken over the observable approximation of
    // that subset: the reads the database placed no further than a genus, right or wrong.
    private long genusOnlyTotal;
    private double genusOnlyScore;
    private long genusOnlyCorrectSpecies;
    private long genusOnlyZeroScoring;
    private long obsGenusOnlyTotal;
    private double obsGenusOnlyUngatedScore;

    /**
     * Records a read whose ground truth could not be resolved.
     */
    public void recordUnresolved() {
        unresolved++;
    }

    /**
     * Records a read that the classifier assigned a taxon to although its ground truth lies outside
     * the scope under consideration. Such a read is a false positive for that scope, so it counts as
     * classified but not towards the total.
     */
    public void recordOutOfScopeClassification() {
        classified++;
    }

    /**
     * Records a read whose ground truth is known and in scope.
     *
     * @param classified   whether the analysis assigned a taxon to the read at all
     * @param lcaRank      the rank of the lowest common ancestor of the read's true and assigned
     *                     taxon, or {@code null} if there is none with a rank of its own -- which is
     *                     not the same as the read being unclassified
     * @param speciesScore the reciprocal of the number of species the classification leaves in
     *                     question, or {@code 0} if the classification says nothing about the read's
     *                     true species; see {@link SpeciesCandidates}
     */
    public void record(boolean classified, Rank lcaRank, double speciesScore) {
        record(classified, lcaRank, speciesScore, speciesScore, false, false);
    }

    /**
     * Records a read, additionally noting whether it belongs to the subset of reads the unrefined
     * database left at their genus -- the only reads a refinement can improve on at all.
     * <p>
     * Two of the parameters exist for the ground-truth-free estimate. {@code ungatedSpeciesScore} is
     * the same reciprocal candidate count as {@code speciesScore} but <em>without</em> the test that
     * the read's true species is still in question at the assigned node: it states how far the
     * classification narrows the species down, not how far it narrows them down correctly, and it is
     * therefore computable on a fastq file whose ground truth is unknown. {@code obsGenusOnly} is
     * the corresponding substitute for the genus-only subset -- the unrefined database assigned the
     * read to a node at genus rank, whether or not that genus is the right one.
     *
     * @param classified          whether the analysis assigned a taxon to the read
     * @param lcaRank             the rank the read's true and assigned taxon agree at, may be {@code null}
     * @param speciesScore        the read's candidate-weighted species score
     * @param ungatedSpeciesScore the same score without the correctness test
     * @param genusOnly           whether the read is in the genus-only subset $R_g$
     * @param obsGenusOnly        whether it is in the observable subset the ungated measure uses
     */
    public void record(boolean classified, Rank lcaRank, double speciesScore, double ungatedSpeciesScore,
                       boolean genusOnly, boolean obsGenusOnly) {
        total++;
        if (obsGenusOnly) {
            obsGenusOnlyTotal++;
            obsGenusOnlyUngatedScore += ungatedSpeciesScore;
        }
        if (genusOnly) {
            genusOnlyTotal++;
            genusOnlyScore += speciesScore;
            if (speciesScore == 0) {
                genusOnlyZeroScoring++;
            }
            if (lcaRank != null && (Rank.SPECIES.equals(lcaRank) || lcaRank.isBelow(Rank.SPECIES))) {
                genusOnlyCorrectSpecies++;
            }
        }
        if (!classified) {
            return;
        }
        this.classified++;
        speciesCandidateScore += speciesScore;
        if (lcaRank == null) {
            return;
        }
        if (Rank.GENUS.equals(lcaRank) || lcaRank.isBelow(Rank.GENUS)) {
            correctGenus++;
        }
        if (Rank.SPECIES.equals(lcaRank) || lcaRank.isBelow(Rank.SPECIES)) {
            correctSpecies++;
        }
        if (Rank.STRAIN.equals(lcaRank) || lcaRank.isBelow(Rank.STRAIN)) {
            correctStrain++;
        }
    }

    /**
     * Records a read of a fastq file whose ground truth is unknown.
     * <p>
     * Only the ungated measures are defined in that case: how far the classification narrows the
     * species down can be read off the assignment alone, whether it does so <em>correctly</em>
     * cannot. The gated counters are therefore left untouched rather than filled with zeros, which
     * would read as "every read wrong" instead of "not determined". {@link #getTotal()} counts every
     * read here, since there is no ground truth by which a read could be out of scope.
     *
     * @param classified          whether the analysis assigned a taxon to the read
     * @param ungatedSpeciesScore the reciprocal number of species the assignment leaves in question
     * @param genusOnly           whether the unrefined database left the read at a genus, i.e.
     *                            whether it belongs to the genus-only subset
     */
    public void recordWithoutGroundTruth(boolean classified, double ungatedSpeciesScore, boolean obsGenusOnly) {
        total++;
        if (obsGenusOnly) {
            obsGenusOnlyTotal++;
            obsGenusOnlyUngatedScore += ungatedSpeciesScore;
        }
        if (classified) {
            this.classified++;
        }
    }

    /**
     * Returns the number of correctly classified reads for the given rank.
     *
     * @param rank {@link Rank#GENUS}, {@link Rank#SPECIES} or {@link Rank#STRAIN}
     * @return the number of reads whose assigned taxon agrees with their true taxon down to
     * {@code rank}
     * @throws IllegalArgumentException if the rank is not one of the three tracked ranks
     */
    public long getCorrect(Rank rank) {
        if (Rank.GENUS.equals(rank)) {
            return correctGenus;
        }
        if (Rank.SPECIES.equals(rank)) {
            return correctSpecies;
        }
        if (Rank.STRAIN.equals(rank)) {
            return correctStrain;
        }
        throw new IllegalArgumentException("No counts tracked for rank " + rank);
    }

    /**
     * Returns the share of classified reads that are correct down to the given rank.
     *
     * @param rank the rank to report for
     * @return the precision, or {@link Double#NaN} if no read was classified at all
     */
    public double getPrecision(Rank rank) {
        return classified == 0 ? Double.NaN : ((double) getCorrect(rank)) / classified;
    }

    /**
     * Returns the share of reads with known ground truth that are correct down to the given rank.
     *
     * @param rank the rank to report for
     * @return the recall, or {@link Double#NaN} if no read has a known ground truth
     */
    public double getRecall(Rank rank) {
        return total == 0 ? Double.NaN : ((double) getCorrect(rank)) / total;
    }

    /**
     * Returns the harmonic mean of {@link #getPrecision(Rank)} and {@link #getRecall(Rank)}.
     *
     * @param rank the rank to report for
     * @return the F1 score, or {@link Double#NaN} if either component is undefined or both are zero
     */
    public double getF1(Rank rank) {
        double precision = getPrecision(rank);
        double recall = getRecall(rank);
        if (Double.isNaN(precision) || Double.isNaN(recall) || precision + recall == 0) {
            return Double.NaN;
        }
        return 2 * precision * recall / (precision + recall);
    }

    /**
     * Returns the share of classified reads that the classification pins down to a species, counting
     * a classification that leaves n species in question as 1/n of a hit. Unlike the plain species
     * precision it credits a classification that narrows the species down without reaching a single
     * one, and unlike the plain genus precision it distinguishes a genus of three species from one
     * of forty. This is the measure the refinement is meant to improve.
     *
     * @return the candidate-weighted species precision, or {@link Double#NaN} if no read was
     * classified at all
     */
    public double getSpeciesCandidatePrecision() {
        return classified == 0 ? Double.NaN : speciesCandidateScore / classified;
    }

    /**
     * Returns the candidate-weighted counterpart of {@link #getRecall(Rank)} for the species rank.
     *
     * @return the candidate-weighted species recall, or {@link Double#NaN} if no read has a known
     * ground truth
     */
    public double getSpeciesCandidateRecall() {
        return total == 0 ? Double.NaN : speciesCandidateScore / total;
    }

    /**
     * Returns the harmonic mean of {@link #getSpeciesCandidatePrecision()} and
     * {@link #getSpeciesCandidateRecall()}.
     *
     * @return the candidate-weighted species F1, or {@link Double#NaN} if it is undefined
     */
    public double getSpeciesCandidateF1() {
        double precision = getSpeciesCandidatePrecision();
        double recall = getSpeciesCandidateRecall();
        if (Double.isNaN(precision) || Double.isNaN(recall) || precision + recall == 0) {
            return Double.NaN;
        }
        return 2 * precision * recall / (precision + recall);
    }

    /**
     * Returns the accumulated candidate weights, i.e. the numerator of the candidate-weighted
     * measures.
     *
     * @return the sum of the per-read species scores
     */
    public double getSpeciesCandidateScore() {
        return speciesCandidateScore;
    }

    /**
     * Returns the number of reads that the unrefined database classified no further than to their
     * genus. This is the subset a refinement can improve on: everything else was either already
     * pinned to a species or was wrong at the genus rank to begin with.
     *
     * @return the size of the genus-only subset
     */
    public long getGenusOnlyTotal() {
        return genusOnlyTotal;
    }

    /**
     * Returns the candidate-weighted species precision restricted to the genus-only subset. Compared
     * between the unrefined and the refined variant, the difference states the gain where a gain was
     * possible, undiluted by the reads that were already at their species.
     *
     * @return the restricted precision, or {@link Double#NaN} if the subset is empty
     */
    public double getGenusOnlyPrecision() {
        return genusOnlyTotal == 0 ? Double.NaN : genusOnlyScore / genusOnlyTotal;
    }

    /**
     * Returns the share of the genus-only subset that this variant pins down to the species.
     *
     * @return the share in {@code [0, 1]}, or {@link Double#NaN} if the subset is empty
     */
    public double getGenusOnlySpeciesShare() {
        return genusOnlyTotal == 0 ? Double.NaN : ((double) genusOnlyCorrectSpecies) / genusOnlyTotal;
    }

    /**
     * Returns the accumulated candidate weights over the genus-only subset.
     *
     * @return the sum of the per-read species scores within the subset
     */
    public double getGenusOnlyScore() {
        return genusOnlyScore;
    }

    /**
     * Returns the ungated counterpart of {@link #getGenusOnlyPrecision()}: the same average over the
     * same subset, but of {@code q'} rather than {@code q}, i.e. of the reciprocal number of species
     * a classification leaves in question regardless of whether the read's true species is among
     * them. The denominator is identical, so the two differ only in what is summed.
     *
     * @return the mean ungated score over the subset, or {@code NaN} if the subset is empty
     */
    public double getObsGenusOnlyUngatedPrecision() {
        return obsGenusOnlyTotal == 0 ? Double.NaN : obsGenusOnlyUngatedScore / obsGenusOnlyTotal;
    }

    /**
     * Returns the size of the observable subset the ungated measure is taken over.
     *
     * @return the number of reads the unrefined database placed no further than a genus
     */
    public long getObsGenusOnlyTotal() {
        return obsGenusOnlyTotal;
    }

    /**
     * Returns how many reads of the genus-only subset score zero under {@code q}, i.e. how many the
     * database placed where their true species is not in question at all. These are the reads on
     * which the gated and the ungated measure disagree, and hence the whole of the difference
     * between them.
     *
     * @return the number of such reads
     */
    public long getGenusOnlyZeroScoring() {
        return genusOnlyZeroScoring;
    }

    /**
     * Returns the number of reads the classifier assigned a taxon to.
     *
     * @return the number of classified reads
     */
    public long getClassified() {
        return classified;
    }

    /**
     * Returns the number of reads whose ground truth is known and in scope.
     *
     * @return the number of reads the recall relates to
     */
    public long getTotal() {
        return total;
    }

    /**
     * Returns the number of reads whose ground truth could not be resolved. A number substantially
     * above zero points at a mismatch between the simulated reads and the accession map, so it is
     * worth reporting rather than hiding.
     *
     * @return the number of reads without resolvable ground truth
     */
    public long getUnresolved() {
        return unresolved;
    }

    /**
     * Returns an independent copy of this tally, used to snapshot the counts once a fastq file is
     * done while the accumulating instance is reset for the next one.
     *
     * @return a copy holding the current counts
     */
    public AccuracyTally copy() {
        AccuracyTally copy = new AccuracyTally();
        copy.classified = classified;
        copy.correctGenus = correctGenus;
        copy.correctSpecies = correctSpecies;
        copy.correctStrain = correctStrain;
        copy.unresolved = unresolved;
        copy.total = total;
        copy.speciesCandidateScore = speciesCandidateScore;
        copy.genusOnlyTotal = genusOnlyTotal;
        copy.genusOnlyScore = genusOnlyScore;
        copy.genusOnlyCorrectSpecies = genusOnlyCorrectSpecies;
        copy.obsGenusOnlyUngatedScore = obsGenusOnlyUngatedScore;
        copy.obsGenusOnlyTotal = obsGenusOnlyTotal;
        copy.genusOnlyZeroScoring = genusOnlyZeroScoring;
        return copy;
    }

    /**
     * Adds another tally's counts to this one.
     * <p>
     * Used to merge the per-thread tallies of one fastq file into a single result. Every counter
     * here is a plain sum or a plain count, so merging is exact and independent of the order the
     * reads were distributed over the threads -- the merged tally is identical to what a single
     * thread recording all reads in sequence would have produced. The two derived averages,
     * {@link #getGenusOnlyPrecision()} and its relatives, divide sums by counts and are therefore
     * correct on the merged tally as well, which they would not be if they were averaged per thread
     * and then averaged again.
     *
     * @param other the tally to add
     */
    public void add(AccuracyTally other) {
        classified += other.classified;
        correctGenus += other.correctGenus;
        correctSpecies += other.correctSpecies;
        correctStrain += other.correctStrain;
        unresolved += other.unresolved;
        total += other.total;
        speciesCandidateScore += other.speciesCandidateScore;
        genusOnlyTotal += other.genusOnlyTotal;
        genusOnlyScore += other.genusOnlyScore;
        genusOnlyCorrectSpecies += other.genusOnlyCorrectSpecies;
        obsGenusOnlyUngatedScore += other.obsGenusOnlyUngatedScore;
        obsGenusOnlyTotal += other.obsGenusOnlyTotal;
        genusOnlyZeroScoring += other.genusOnlyZeroScoring;
    }

    /**
     * Resets all counts to zero.
     */
    public void reset() {
        classified = 0;
        correctGenus = 0;
        correctSpecies = 0;
        correctStrain = 0;
        unresolved = 0;
        total = 0;
        speciesCandidateScore = 0;
        genusOnlyTotal = 0;
        genusOnlyScore = 0;
        genusOnlyCorrectSpecies = 0;
        obsGenusOnlyUngatedScore = 0;
        obsGenusOnlyTotal = 0;
        genusOnlyZeroScoring = 0;
    }

    @Override
    public String toString() {
        return "classified=" + classified + ", genus=" + correctGenus + ", species=" + correctSpecies
                + ", strain=" + correctStrain + ", total=" + total + ", unresolved=" + unresolved
                + ", speciesScore=" + speciesCandidateScore;
    }
}
