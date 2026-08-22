package org.metagene.ftdbexp;

import org.metagene.ftdbexp.taxquality.PathVoteTaxonModel;
import org.metagene.ftdbexp.taxquality.TaxonCallReport;

import java.io.File;

/**
 * Runs the path vote of the case study: which taxon below a given node each sample's reads point at,
 * for the unrefined and for the refined database.
 * <p>
 * The counterpart of {@link FTSpecificityMain}, which measures how much more specific the refined
 * answers are without naming any of them. Here they are named, one per sample, so that they can be
 * held against a reference standard the source study states per sample rather than per read.
 */
public class FTTaxonCallMain {
    private static final File BASE_DIR = new File("./data");
    private static final File RESULTS_DIR = new File("./results");

    /**
     * @param args the database project name, the fastq mapping file, optionally a report key, the
     *             tax id the vote is restricted to (default {@code 1301}, Streptococcus) and the
     *             minimum the winning path must gather (default none)
     * @throws Exception if a database cannot be read or the report cannot be written
     */
    public static void main(String[] args) throws Exception {
        if (args.length < 2) {
            System.err.println("Usage: FTTaxonCallMain <db> <fastq map file> [<report key>] [<root taxid>] [<minimum>]");
            System.err.println("  <db>             name of the database project under data/projects");
            System.err.println("  <fastq map file> mapping file resolved as usual against data/fastq");
            System.err.println("  <report key>     short name used in the result file name");
            System.err.println("  <root taxid>     the node the vote is restricted to, default 1301");
            System.err.println("  <minimum>        how much the winning path must gather, default none");
            System.exit(1);
        }
        String db = args[0];
        String fqMapFile = args[1];
        String reportKey = args.length > 2 && !args[2].isEmpty() ? args[2] : stripExtension(fqMapFile);
        String rootTaxId = args.length > 3 && !args[3].isEmpty() ? args[3] : "1301";
        long minimum = args.length > 4 && !args[4].isEmpty()
                ? Long.parseLong(args[4]) : PathVoteTaxonModel.NO_MINIMUM;

        new TaxonCallReport(BASE_DIR, RESULTS_DIR).write(db, fqMapFile, reportKey, rootTaxId, minimum);
    }

    /**
     * @param fileName the file name
     * @return the file name without its extension
     */
    private static String stripExtension(String fileName) {
        int dot = fileName.lastIndexOf('.');
        return dot > 0 ? fileName.substring(0, dot) : fileName;
    }
}
