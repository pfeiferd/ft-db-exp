package org.metagene.ftdbexp.eval;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * The sequence-to-taxon table that the goal {@code extractrefseqcsv} writes next to the extracted
 * FASTA files.
 * <p>
 * The goal records, for every region it extracts, which taxon the database associates it with. That
 * is precisely the ground truth for reads simulated from those files, and it is the only source
 * that covers <em>all</em> of them: for a genome taken from Genbank, Genestrip attaches the file to
 * its taxon rather than resolving the accession, so the accession map does not know such a sequence
 * at all.
 * <p>
 * The file is semicolon-separated with one header line:
 * <pre>
 *   refseq descr; taxid;
 *   JBHNWG010001660.1;944036;
 * </pre>
 */
public final class ExtractedTaxIds {
    private ExtractedTaxIds() {
    }

    /**
     * Loads the table, keyed by sequence description.
     *
     * @param csvFile the CSV written by {@code extractrefseqcsv}
     * @return the description-to-taxid map, or an empty map if the file does not exist, in which
     * case the caller is expected to fall back to the accession map
     * @throws IOException if the file exists but cannot be read
     */
    public static Map<String, String> load(File csvFile) throws IOException {
        if (csvFile == null || !csvFile.exists()) {
            return Collections.emptyMap();
        }
        Map<String, String> result = new HashMap<String, String>();
        BufferedReader reader = new BufferedReader(
                new InputStreamReader(new FileInputStream(csvFile), StandardCharsets.UTF_8));
        try {
            String line = reader.readLine(); // Header.
            while ((line = reader.readLine()) != null) {
                int sep = line.indexOf(';');
                if (sep <= 0) {
                    continue;
                }
                int end = line.indexOf(';', sep + 1);
                if (end < 0) {
                    end = line.length();
                }
                String descr = line.substring(0, sep).trim();
                String taxId = line.substring(sep + 1, end).trim();
                if (!descr.isEmpty() && !taxId.isEmpty()) {
                    result.put(descr, taxId);
                }
            }
        } finally {
            reader.close();
        }
        return result;
    }
}
