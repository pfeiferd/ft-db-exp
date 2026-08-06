package org.metagene.ftdbexp;

import org.metagene.ftdbexp.eval.SpecificityReport;

import java.io.File;

/**
 * Command line entry point for the experiments on real reads, i.e. on fastq files whose ground truth
 * is not known.
 * <p>
 * It is invoked from Maven, e.g.
 * <pre>
 *   mvn exec:exec@specificity -Dname=viral -Dfqmap=saliva_real.txt -Dreportkey=saliva
 * </pre>
 * and writes {@code results/<db>_<report key>_specificity.csv}.
 * <p>
 * Unlike {@link FTExpMain} this reports no precision and no recall. Neither is defined without the
 * species a read stems from. What it does report is how far each database variant narrows the
 * species down on the reads the unrefined one left at a genus, and the difference between the two --
 * the specificity gain that bounds the precision gain from above.
 */
public class FTSpecificityMain {
    private static final File BASE_DIR = new File("./data");
    private static final File RESULTS_DIR = new File("./results");

    /**
     * Runs the ground-truth-free comparison between the unrefined and the refined database.
     *
     * @param args the database project name, the fastq mapping file and, optionally, a report key
     *             used in the output file name; the report key defaults to the mapping file name
     *             without its extension
     * @throws Exception if the databases cannot be read or the report cannot be written
     */
    public static void main(String[] args) throws Exception {
        if (args.length < 2) {
            System.err.println("Usage: FTSpecificityMain <db> <fastq map file> [<report key>]");
            System.err.println("  <db>             name of the database project under data/projects");
            System.err.println("  <fastq map file> mapping file resolved as usual against data/fastq");
            System.err.println("  <report key>     short name used in the result file name");
            System.exit(1);
        }
        String db = args[0];
        String fqMapFile = args[1];
        String reportKey = args.length > 2 && !args[2].isEmpty() ? args[2] : stripExtension(fqMapFile);

        new SpecificityReport(BASE_DIR, RESULTS_DIR, db).write(db, fqMapFile, reportKey);
    }

    /**
     * Strips a trailing file extension.
     *
     * @param fileName the file name
     * @return the file name without its extension
     */
    private static String stripExtension(String fileName) {
        int dot = fileName.lastIndexOf('.');
        return dot > 0 ? fileName.substring(0, dot) : fileName;
    }
}
