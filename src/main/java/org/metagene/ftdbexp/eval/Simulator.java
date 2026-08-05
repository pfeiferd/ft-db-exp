package org.metagene.ftdbexp.eval;

import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.tax.TaxTree;

import java.util.Map;

/**
 * The read simulators the experiments use, each with the way it encodes a read's origin in the read
 * identifier. Which one produced a fastq file determines how its ground truth has to be read.
 */
public enum Simulator {
    /** InSilicoSeq, used for the Illumina reads of the viral experiments. */
    ISS {
        @Override
        public boolean needsAccessionMap() {
            return false;
        }

        @Override
        public ReadGroundTruth groundTruth(TaxTree taxTree, AccessionMap accessionMap,
                                          Map<String, String> extractedTaxIds) {
            return new IssReadGroundTruth(taxTree, extractedTaxIds);
        }
    },
    /** NanoSim, used for the Nanopore reads of the tick-borne experiments. */
    NANOSIM {
        @Override
        public ReadGroundTruth groundTruth(TaxTree taxTree, AccessionMap accessionMap,
                                          Map<String, String> extractedTaxIds) {
            return new NanoSimReadGroundTruth(accessionMap);
        }
    };

    /**
     * Returns the ground truth resolver for this simulator's read identifiers.
     *
     * @param taxTree      the taxonomy, used where a read names its taxon directly
     * @param accessionMap    the accession-to-taxon map of the database the reads are matched against
     * @param extractedTaxIds the sequence-to-taxon table written alongside the extracted genomes
     * @return the matching resolver
     */
    public abstract ReadGroundTruth groundTruth(TaxTree taxTree, AccessionMap accessionMap,
                                                Map<String, String> extractedTaxIds);

    /**
     * Parses a simulator name as given on the command line, case-insensitively.
     *
     * @param name the simulator name, e.g. {@code iss} or {@code nanosim}
     * @return the matching simulator
     * @throws IllegalArgumentException if the name denotes no known simulator
     */
    /**
     * Whether this simulator's resolver consults Genestrip's accession map. Building that map is
     * expensive, so an evaluation skips it where the ground truth comes from elsewhere.
     *
     * @return whether the accession map has to be supplied
     */
    public boolean needsAccessionMap() {
        return true;
    }

    public static Simulator parse(String name) {
        for (Simulator simulator : values()) {
            if (simulator.name().equalsIgnoreCase(name)) {
                return simulator;
            }
        }
        throw new IllegalArgumentException("Unknown simulator '" + name + "', expected one of ISS, NANOSIM");
    }
}
