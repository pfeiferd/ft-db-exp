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
        int end = ByteArrayUtil.indexOf(descriptor, ACCESSION_PREFIX_END, length, '_');
        if (end < 0) {
            return null;
        }
        return accessionMap.get(descriptor, start, end, false);
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
