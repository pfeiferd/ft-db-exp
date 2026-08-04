package org.metagene.ftdbexp.eval;

import org.metagene.genestrip.refseq.AccessionMap;

/**
 * The read simulators the experiments use, each with the way it encodes a read's origin in the read
 * identifier. Which one produced a fastq file determines how its ground truth has to be read.
 */
public enum Simulator {
    /** InSilicoSeq, used for the Illumina reads of the viral experiments. */
    ISS {
        @Override
        public ReadGroundTruth groundTruth(AccessionMap accessionMap) {
            return new IssReadGroundTruth(accessionMap);
        }
    },
    /** NanoSim, used for the Nanopore reads of the tick-borne experiments. */
    NANOSIM {
        @Override
        public ReadGroundTruth groundTruth(AccessionMap accessionMap) {
            return new NanoSimReadGroundTruth(accessionMap);
        }
    };

    /**
     * Returns the ground truth resolver for this simulator's read identifiers.
     *
     * @param accessionMap the accession-to-taxon map of the database the reads are matched against
     * @return the matching resolver
     */
    public abstract ReadGroundTruth groundTruth(AccessionMap accessionMap);

    /**
     * Parses a simulator name as given on the command line, case-insensitively.
     *
     * @param name the simulator name, e.g. {@code iss} or {@code nanosim}
     * @return the matching simulator
     * @throws IllegalArgumentException if the name denotes no known simulator
     */
    public static Simulator parse(String name) {
        for (Simulator simulator : values()) {
            if (simulator.name().equalsIgnoreCase(name)) {
                return simulator;
            }
        }
        throw new IllegalArgumentException("Unknown simulator '" + name + "', expected one of ISS, NANOSIM");
    }
}
