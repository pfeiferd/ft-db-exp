package org.metagene.ftdbexp.eval;

import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.ByteArrayUtil;

import java.nio.charset.StandardCharsets;

/**
 * Ground truth for reads simulated by <a href="https://github.com/bcgsc/NanoSim">NanoSim</a>.
 * <p>
 * NanoSim names a read after the entry of its genome list followed by the source sequence, e.g.
 * <pre>
 *   &#64;119857x137145-NZ-PKBC01000002_8514_aligned_11841_R_45_6436_373
 * </pre>
 * The part before the first dash is the label from the genome list, {@code <taxid>x<index>}, and
 * the taxon is read straight out of it. What follows is the accession of the source sequence, but
 * mangled twice: NanoSim replaces the accession's own underscore by a dash ({@code NZ_PKBC01000002}
 * becomes {@code NZ-PKBC01000002}) and drops the version suffix.
 * <p>
 * Resolving via that accession would mean undoing both mutilations and then guessing the lost
 * version by trying {@code .1} through {@code .9} against Genestrip's accession map -- and it would
 * fail outright for every genome taken from Genbank, which the map does not know because Genestrip
 * attaches such a file to its taxon by file rather than by accession. For a database drawing on
 * Genbank that is most of the reads. The label carries the answer already, exactly and without
 * guessing, and {@code NanoSimGenomeList} in this project is what writes it -- so the convention is
 * guaranteed here rather than assumed.
 */
public class NanoSimReadGroundTruth implements ReadGroundTruth {
    private final TaxTree taxTree;
    /** Guards the warning below, so that a systematic mismatch is reported once and not per read. */
    private boolean warned;

    /**
     * Creates the ground truth resolver.
     *
     * @param taxTree the taxonomy the label's tax id is looked up in
     */
    public NanoSimReadGroundTruth(TaxTree taxTree) {
        this.taxTree = taxTree;
    }

    @Override
    public TaxTree.TaxIdNode resolve(byte[] descriptor, int length) {
        int start = length > 0 && descriptor[0] == '@' ? 1 : 0;
        int labelEnd = ByteArrayUtil.indexOf(descriptor, start, length, '-');
        if (labelEnd < 0) {
            warnOnce("no label before a dash", descriptor, start, length);
            return null;
        }
        // The label is "<taxid>x<index>"; everything before the 'x' is the taxon.
        int x = ByteArrayUtil.indexOf(descriptor, start, labelEnd, 'x');
        if (x <= start) {
            warnOnce("label is not of the form <taxid>x<index>", descriptor, start, length);
            return null;
        }
        String taxId = new String(descriptor, start, x - start, StandardCharsets.UTF_8);
        TaxTree.TaxIdNode node = taxTree.getNodeByTaxId(taxId);
        if (node == null) {
            warnOnce("tax id " + taxId + " is not in the taxonomy", descriptor, start, length);
        }
        return node;
    }

    /**
     * Reports the first read whose label yields no taxon, then counts the rest silently.
     *
     * @param reason     what went wrong, for the message
     * @param descriptor the raw read descriptor
     * @param start      the index the label starts at
     * @param length     the number of valid bytes in {@code descriptor}
     */
    private void warnOnce(String reason, byte[] descriptor, int start, int length) {
        if (warned) {
            return;
        }
        warned = true;
        System.err.println("WARNING: cannot resolve the ground truth of read '"
                + new String(descriptor, start, length - start, StandardCharsets.UTF_8).trim() + "': " + reason
                + ". Such reads are counted as unresolved and excluded from every measure."
                + " Was the genome list written by NanoSimGenomeList?"
                + " Further occurrences are counted but not reported.");
    }
}
