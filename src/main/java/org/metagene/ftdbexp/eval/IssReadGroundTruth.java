package org.metagene.ftdbexp.eval;

import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.ByteArrayUtil;

/**
 * Ground truth for reads simulated by <a href="https://insilicoseq.readthedocs.io">InSilicoSeq</a>.
 * <p>
 * InSilicoSeq names a read after the fasta record it was drawn from, appending an underscore, a
 * running number and the mate suffix, e.g.
 * <pre>
 *   &#64;NC_001422.1_1234/1
 * </pre>
 * The leading part up to the <em>last</em> underscore before the running number is the accession of
 * the source sequence, which the {@link AccessionMap} maps to a taxon.
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
 * Here the accession ends at the first {@code '|'} rather than at an underscore -- searching for the
 * underscore would run past the taxon and yield {@code NC_001422.1|kraken:taxid|10847}, which no
 * accession map resolves. Both forms are therefore accepted: whichever of the two delimiters comes
 * first ends the accession.
 */
public class IssReadGroundTruth implements ReadGroundTruth {
    /**
     * Offset the separator search starts at, chosen to skip the underscore inside an accession
     * prefix such as {@code NC_}.
     */
    private static final int ACCESSION_PREFIX_END = 5;

    private final AccessionMap accessionMap;

    /**
     * Creates the ground truth resolver.
     *
     * @param accessionMap the accession-to-taxon map of the database the reads are matched against
     */
    public IssReadGroundTruth(AccessionMap accessionMap) {
        this.accessionMap = accessionMap;
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
        return accessionMap.get(descriptor, start, end, false);
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
