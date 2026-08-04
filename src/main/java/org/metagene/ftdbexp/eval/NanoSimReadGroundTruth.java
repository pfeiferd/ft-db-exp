package org.metagene.ftdbexp.eval;

import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.ByteArrayUtil;

/**
 * Ground truth for reads simulated by <a href="https://github.com/bcgsc/NanoSim">NanoSim</a>.
 * <p>
 * NanoSim names a read after the entry of its genome list followed by the source sequence, e.g.
 * <pre>
 *   &#64;119857x137145-NZ-PKBC01000002_8514_aligned_11841_R_45_6436_373
 * </pre>
 * The part before the first dash is the label from the genome list, which the experiments of the
 * first Genestrip paper chose to be {@code <taxid>x<index>}. What follows up to the next underscore
 * is the accession of the source sequence -- but mangled twice: NanoSim replaces the accession's own
 * underscore by a dash ({@code NZ_PKBC01000002} becomes {@code NZ-PKBC01000002}) and it drops the
 * version suffix.
 * <p>
 * Both mutilations have to be undone before the {@link AccessionMap} can be consulted. The dash is
 * simply turned back into an underscore; the version, being lost for good, is guessed by trying
 * {@code .1} through {@code .9} and taking the first one the map knows. That is what the original
 * experiment did as well, and in practice the first or second attempt hits.
 * <p>
 * Note that the read identifier also carries the taxon directly, in the label before the dash. That
 * would avoid the guessing entirely, but it only works for reads generated from a genome list that
 * follows the {@code <taxid>x<index>} convention, whereas the accession is always present. The
 * accession route is therefore the more robust one and the one implemented here.
 */
public class NanoSimReadGroundTruth implements ReadGroundTruth {
    /** Highest version suffix tried when reconstructing the accession. */
    private static final int MAX_VERSION = 9;

    private final AccessionMap accessionMap;

    /**
     * Creates the ground truth resolver.
     *
     * @param accessionMap the accession-to-taxon map of the database the reads are matched against
     */
    public NanoSimReadGroundTruth(AccessionMap accessionMap) {
        this.accessionMap = accessionMap;
    }

    @Override
    public TaxTree.TaxIdNode resolve(byte[] descriptor, int length) {
        int labelEnd = ByteArrayUtil.indexOf(descriptor, 0, length, '-');
        if (labelEnd < 0) {
            return null;
        }
        // The underscore ending the accession; the search starts behind the label so that a dash
        // inside it cannot be mistaken for the separator.
        int accessionEnd = ByteArrayUtil.indexOf(descriptor, labelEnd + 1, length, '_');
        if (accessionEnd < 0 || accessionEnd + 1 >= length) {
            return null;
        }

        // Repair the accession in place: the dash NanoSim introduced becomes an underscore again,
        // and the trailing underscore becomes the dot of the version suffix.
        int innerDash = ByteArrayUtil.indexOf(descriptor, labelEnd + 1, accessionEnd, '-');
        byte savedDash = 0;
        if (innerDash >= 0) {
            savedDash = descriptor[innerDash];
            descriptor[innerDash] = '_';
        }
        byte savedEnd = descriptor[accessionEnd];
        byte savedVersion = descriptor[accessionEnd + 1];
        descriptor[accessionEnd] = '.';
        try {
            for (int version = 1; version <= MAX_VERSION; version++) {
                descriptor[accessionEnd + 1] = (byte) ('0' + version);
                TaxTree.TaxIdNode node =
                        accessionMap.get(descriptor, labelEnd + 1, accessionEnd + 2, false);
                if (node != null) {
                    return node;
                }
            }
            return null;
        } finally {
            // Leave the caller's buffer as we found it - it is reused for the next read.
            descriptor[accessionEnd] = savedEnd;
            descriptor[accessionEnd + 1] = savedVersion;
            if (innerDash >= 0) {
                descriptor[innerDash] = savedDash;
            }
        }
    }
}
