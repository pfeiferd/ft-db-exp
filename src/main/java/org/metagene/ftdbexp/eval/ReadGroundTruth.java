package org.metagene.ftdbexp.eval;

import org.metagene.genestrip.tax.TaxTree;

/**
 * Resolves the taxon a simulated read was actually generated from, i.e. its ground truth.
 * <p>
 * Read simulators encode the source sequence in the read descriptor, but each of them does so in its
 * own dialect. An implementation of this interface knows one such dialect and maps a descriptor to
 * the corresponding node of the taxonomy.
 *
 * @see IssReadGroundTruth
 */
public interface ReadGroundTruth {
    /**
     * Resolves the taxon a read was generated from.
     *
     * @param descriptor the raw read descriptor as read from the fastq file, possibly including a
     *                   leading {@code @}
     * @param length     the number of valid bytes in {@code descriptor}
     * @return the taxonomy node the read originates from, or {@code null} if it cannot be resolved,
     * e.g. because the source sequence is unknown to the accession map
     */
    TaxTree.TaxIdNode resolve(byte[] descriptor, int length);
}
