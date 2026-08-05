package org.metagene.ftdbexp.eval;

import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.ByteArrayUtil;

import java.nio.charset.StandardCharsets;
import java.util.Map;

/**
 * Ground truth for reads simulated by <a href="https://insilicoseq.readthedocs.io">InSilicoSeq</a>.
 * <p>
 * InSilicoSeq names a read after the fasta record it was drawn from, appending an underscore, a
 * running number and the mate suffix, e.g.
 * <pre>
 *   &#64;NC_001422.1_1234/1
 * </pre>
 * The leading part up to the <em>last</em> underscore before the running number is the accession of
 * the source sequence, which {@link ExtractedTaxIds} maps to a taxon.
 * <p>
 * Note that RefSeq accessions contain an underscore themselves ({@code NC_001422.1}), so the
 * separator cannot be found by simply looking for the first one. Since every accession prefix is at
 * most four characters long, the search starts behind that prefix instead.
 * <p>
 * Reads generated from the output of the goal {@code extractrefseqfasta} carry a longer name,
 * because that goal appends the taxon to every fasta header to please Kraken 2 during library
 * building:
 * <pre>
 *   &#64;NC_001422.1|kraken:taxid|10847_1234/1
 * </pre>
 * The taxon comes from {@link ExtractedTaxIds}, the table the same goal writes for exactly this
 * purpose, and from nowhere else. Genestrip's accession map is deliberately not consulted: for a
 * genome taken from Genbank it attaches the file to its taxon rather than resolving the accession
 * -- {@code ignoreAccessionMap} in {@code FastaReaderGoal} -- so it does not know the sequence at
 * all. Falling back to it would leave every such read unresolved, which for a database drawing
 * substantially on Genbank is most of them, while looking like an ordinary result.
 */
public class IssReadGroundTruth implements ReadGroundTruth {
    /**
     * Offset the separator search starts at, chosen to skip the underscore inside an accession
     * prefix such as {@code NC_}.
     */
    private static final int ACCESSION_PREFIX_END = 5;

    private final TaxTree taxTree;
    private final Map<String, String> extractedTaxIds;
    /** Guards the warning below, so that a systematic mismatch is reported once and not per read. */
    private boolean warned;

    /**
     * Creates the ground truth resolver.
     *
     * @param taxTree         the taxonomy, used to look up the taxa named by the extraction table
     * @param extractedTaxIds the sequence-to-taxon table of {@link ExtractedTaxIds}
     * @throws IllegalArgumentException if the table is empty, since every read would then be
     *                                  unresolved -- a missing table is a setup error worth failing
     *                                  on rather than a result worth reporting
     */
    public IssReadGroundTruth(TaxTree taxTree, Map<String, String> extractedTaxIds) {
        if (extractedTaxIds.isEmpty()) {
            throw new IllegalArgumentException(
                    "No extracted genomes to resolve the ground truth against. Run the goal "
                            + "'extractrefseqcsv' for this project and regenerate the reads from its output.");
        }
        this.taxTree = taxTree;
        this.extractedTaxIds = extractedTaxIds;
    }

    @Override
    public TaxTree.TaxIdNode resolve(byte[] descriptor, int length) {
        int start = accessionStart(descriptor, length);
        if (start < 0) {
            return null;
        }
        int end = accessionEnd(descriptor, start, length);
        if (end < 0) {
            return null;
        }
        String descr = new String(descriptor, start, end - start, StandardCharsets.UTF_8);
        String taxId = extractedTaxIds.get(descr);
        if (taxId == null) {
            warnOnce("no entry for '" + descr + "'", descriptor, start, length);
            return null;
        }
        TaxTree.TaxIdNode node = taxTree.getNodeByTaxId(taxId);
        if (node == null) {
            // The table names a taxon the taxonomy does not have - a stale extraction against a
            // newer NCBI dump, for instance. Worth distinguishing from an unknown sequence.
            warnOnce("tax id " + taxId + " of '" + descr + "' is not in the taxonomy", descriptor, start, length);
        }
        return node;
    }

    /**
     * Reports the first read that could not be resolved, with the reason and the offending name.
     * <p>
     * Every further one is only counted, by {@link AccuracyTally#recordUnresolved()}, whose total is
     * reported per fastq file and written to the result CSV. One read failing may be an oddity; a
     * substantial count means the reads and the extraction table do not belong together, and the
     * figures then rest on whatever fraction did resolve.
     *
     * @param reason     what went wrong, for the message
     * @param descriptor the raw read descriptor
     * @param start      the index the accession starts at
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
                + " Were these reads generated from this project's extracted genomes?"
                + " Further occurrences are counted but not reported.");
    }



    /**
     * Determines where the accession ends, i.e. the first delimiter behind it.
     * <p>
     * Two delimiters are possible. A plain fasta header leaves only the underscore that InSilicoSeq
     * puts before the running read number; a header written by {@code extractrefseqfasta} carries
     * the taxon behind a {@code '|'}, which then comes first. Taking whichever occurs earlier
     * handles both, and it stays correct if the suffix is ever extended, since the {@code '|'}
     * remains the first character that cannot belong to an accession.
     *
     * @param descriptor the raw read descriptor
     * @param start      the index the accession starts at
     * @param length     the number of valid bytes in {@code descriptor}
     * @return the index behind the accession, or {@code -1} if no delimiter follows it
     */
    private static int accessionEnd(byte[] descriptor, int start, int length) {
        // The underscore search skips the accession prefix, which contains one itself ("NC_").
        int underscore = ByteArrayUtil.indexOf(descriptor, ACCESSION_PREFIX_END, length, '_');
        int bar = ByteArrayUtil.indexOf(descriptor, start, length, '|');
        if (bar < 0) {
            return underscore;
        }
        return underscore < 0 || bar < underscore ? bar : underscore;
    }

    /**
     * Determines where the accession begins, skipping the fastq/fasta record markers that may
     * precede it.
     *
     * @param descriptor the raw read descriptor
     * @param length     the number of valid bytes in {@code descriptor}
     * @return the index the accession starts at, or {@code -1} if the descriptor is too short
     */
    private static int accessionStart(byte[] descriptor, int length) {
        if (length <= ACCESSION_PREFIX_END) {
            return -1;
        }
        int start = 0;
        if (descriptor[0] == '@' || descriptor[0] == '\t') {
            start = 1;
        }
        if (descriptor[start] == '>') {
            start++;
        }
        return start;
    }
}
