package org.metagene.ftdbexp.stquality;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

/**
 * The sequence type of every isolate whose reads the classification experiment scores, keyed by the
 * fastq map key the reads arrive under.
 * <p>
 * This is a second, independent ground truth beside {@link STGroundTruth}, and the two must not be
 * confused. {@code STGroundTruth} types the genomes the <em>database was filled from</em> and is
 * what the intrinsic measures compare a k-mer's placement against. This one types the isolates the
 * database has <em>never seen</em>, and is what a classification is compared against. An isolate
 * that also contributed a genome to the database would make the experiment worthless, so the two
 * files are expected to be disjoint and {@code IsolateSTAccuracyMain} says so if they are not.
 * <p>
 * The format is one row per fastq key, semicolon separated, with a header line and {@code #}
 * comments:
 * <pre>
 *   isolate;st;
 *   B11;11;
 *   B11np;11;
 * </pre>
 * The Illumina and Nanopore runs of one isolate are separate keys and each needs its own row; that
 * they carry the same type is what makes the pair a platform comparison rather than two samples.
 * <p>
 * For the isolates of {@code data/fastq/cdiff_isolates.txt} the types come from Table 2 of the
 * source study (Bejaoui et al., BMC Genomics 2025), which keys them by the same identifiers the
 * archive carries as {@code sample_alias} -- so the join is by name and needs no accession lookup.
 */
public class IsolateSTTruth {
    private final Map<String, String> stByKey = new LinkedHashMap<>();

    /**
     * Reads the typing.
     *
     * @param csv the file to read
     * @throws IOException if it cannot be read
     * @throws IllegalStateException if a key appears twice with differing types, which would make
     *                               the isolate's verdict depend on the row order
     */
    public IsolateSTTruth(File csv) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(csv))) {
            String line;
            boolean first = true;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) {
                    continue;
                }
                String[] parts = line.split(";");
                if (first && parts.length > 1 && "isolate".equalsIgnoreCase(parts[0].trim())) {
                    first = false;
                    continue;
                }
                first = false;
                if (parts.length < 2) {
                    continue;
                }
                String key = parts[0].trim();
                String st = parts[1].trim();
                if (key.isEmpty() || st.isEmpty()) {
                    continue;
                }
                String previous = stByKey.put(key, st);
                if (previous != null && !previous.equals(st)) {
                    throw new IllegalStateException("Isolate '" + key + "' is typed as both '"
                            + previous + "' and '" + st + "' in " + csv + ". Which of the two a run"
                            + " is scored against would then depend on the order of the rows.");
                }
            }
        }
    }

    /**
     * Returns the sequence type of the isolate the given fastq key belongs to.
     *
     * @param key the fastq map key
     * @return the sequence type, or {@code null} if the key is not typed
     */
    public String getST(String key) {
        return stByKey.get(key);
    }

    /**
     * Returns the fastq keys this typing covers.
     *
     * @return the keys, in the order they were read
     */
    public Set<String> getKeys() {
        return Collections.unmodifiableSet(stByKey.keySet());
    }

    /**
     * Returns how many isolates are typed.
     *
     * @return the number of typed fastq keys
     */
    public int size() {
        return stByKey.size();
    }
}
