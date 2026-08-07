package org.metagene.ftdbexp.eval;

/**
 * Turns a fastq key into the name the paper's tables use.
 * <p>
 * The tick samples appear twice in the paper and must not be confused: Table \ref{genusonly} scores
 * the NanoSim simulation trained on a sample, Table \ref{realgain} scores the sample itself. Both
 * derive from the same fastq key, so the distinction has to be made here -- {@code Sim. Tick 3}
 * against {@code Tick 3}.
 * <p>
 * These are labels, never join keys. {@link SpecificityReport} matches its rows to the calibration
 * by the raw fastq key precisely so that renaming a label cannot silently break the join.
 */
final class SampleNames {

    private SampleNames() {
    }

    /**
     * @param fastqKey  the key of the fastq file
     * @param simulated whether the reads were simulated from the sample rather than being it
     * @return {@code tick3} as {@code Sim. Tick 3} or {@code Tick 3}; anything else unchanged,
     * which leaves the SRA accessions and the InSilicoSeq model names as they are
     */
    static String display(String fastqKey, boolean simulated) {
        if (fastqKey.length() > 4 && fastqKey.startsWith("tick")) {
            String suffix = fastqKey.substring(4);
            for (int i = 0; i < suffix.length(); i++) {
                if (!Character.isDigit(suffix.charAt(i))) {
                    return fastqKey;
                }
            }
            return (simulated ? "Sim. Tick " : "Tick ") + suffix;
        }
        return fastqKey;
    }
}
