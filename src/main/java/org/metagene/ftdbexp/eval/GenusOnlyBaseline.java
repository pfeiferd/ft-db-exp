package org.metagene.ftdbexp.eval;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Remembers which reads the unrefined database classified no further than to their genus, so that
 * the refined run can be measured on exactly those reads.
 * <p>
 * The overall gain of a refinement is easy to misread: most reads are already assigned to their
 * species before any refinement, and they dilute every average. A read can only profit from the
 * refinement if the unrefined database left it at the genus -- correct at that rank, but not
 * specific enough to name the species. Relating the gain to this subset states what the refinement
 * achieves where it can achieve anything at all.
 * <p>
 * An external classifier uses the same class for its own subset, collected from its own output and
 * never consulted afterwards, since it has no refined second pass to score.
 * <p>
 * The subset is kept per fastq file. Read identifiers are only unique within a file: a simulator
 * numbers the reads of each genome from zero, so the same identifier may well occur in two
 * different fastq files of the same collection. Both evaluation runs walk the same mapping file, so
 * the files come in the same order, and this class follows that order with a cursor. It verifies
 * the file names as it goes and fails loudly on a mismatch rather than silently comparing the reads
 * of one file against the subset of another.
 */
public class GenusOnlyBaseline {
    private final List<String> fastqKeys = new ArrayList<String>();
    private final List<Set<String>> subsets = new ArrayList<Set<String>>();

    /**
     * Reads collected for the file currently being processed. Concurrent, because the matcher
     * reports its reads from every worker thread at once; the finished per-file sets below are only
     * ever read afterwards and need no synchronisation of their own.
     */
    private Set<String> collecting = Collections.newSetFromMap(new ConcurrentHashMap<String, Boolean>());
    /** Index of the file the consulting run is currently processing. */
    private int cursor = 0;

    /**
     * Records a read of the file currently being processed as one the unrefined database left at
     * its genus.
     *
     * @param readDescriptor the read's identifier
     */
    public void collect(String readDescriptor) {
        collecting.add(readDescriptor);
    }

    /**
     * Puts the cursor back to the first file, so that a further run can be scored on the same
     * subsets. The subsets themselves are untouched; only the position in them is reset. This is
     * what lets an external classifier be measured on the reads of the unrefined Genestrip run
     * after the refined one has already walked them.
     */
    public void rewind() {
        cursor = 0;
    }

    /**
     * Closes the subset of the file just finished.
     *
     * @param fastqKey the key of that file
     */
    public void endCollecting(String fastqKey) {
        fastqKeys.add(fastqKey);
        subsets.add(collecting);
        collecting = Collections.newSetFromMap(new ConcurrentHashMap<String, Boolean>());
    }

    /**
     * Returns whether the given read of the file currently being processed belongs to the subset.
     *
     * @param readDescriptor the read's identifier
     * @return whether the unrefined database left this read at its genus
     */
    public boolean contains(String readDescriptor) {
        return cursor < subsets.size() && subsets.get(cursor).contains(readDescriptor);
    }

    /**
     * Advances to the next file.
     *
     * @param fastqKey the key of the file just finished, used to verify that both runs walk the
     *                 files in the same order
     * @throws IllegalStateException if the file does not match the one recorded at this position
     */
    public void endConsulting(String fastqKey) {
        if (cursor >= fastqKeys.size() || !fastqKeys.get(cursor).equals(fastqKey)) {
            throw new IllegalStateException("Fastq files are processed in a different order than when"
                    + " the baseline was collected: expected "
                    + (cursor < fastqKeys.size() ? fastqKeys.get(cursor) : "no further file")
                    + " but got " + fastqKey);
        }
        cursor++;
    }

    /**
     * Returns the number of reads in the subset of the file currently being processed.
     *
     * @return the size of the current subset
     */
    public int currentSize() {
        return cursor < subsets.size() ? subsets.get(cursor).size() : 0;
    }
}
