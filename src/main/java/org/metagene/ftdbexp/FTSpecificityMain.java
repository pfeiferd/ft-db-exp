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
 * and writes {@code results/<db>_<report key>_specificity.csv}. With
 * {@code -Dgs.externalonly=true} it rewrites the rows of the external classifiers alone, taking the
 * read counts from the CSV of a previous full run.
 * <p>
 * Unlike {@link FTExpMain} this reports no precision and no recall. Neither is defined without the
 * species a read stems from. What it does report is how far each database variant narrows the
 * species down on the reads the unrefined one left at a genus. That difference is not a bound on the
 * precision gain in either direction; it becomes an estimate of the gated precisions only through
 * the calibration factors of a simulated run, which {@link org.metagene.ftdbexp.eval.SpecificityReport}
 * applies.
 */
public class FTSpecificityMain {
    private static final File BASE_DIR = new File("./data");
    private static final File RESULTS_DIR = new File("./results");

    /**
     * Runs the ground-truth-free comparison between the unrefined and the refined database.
     *
     * @param args the database project name, the fastq mapping file and, optionally, a report key
     *             and the report key of the simulated run supplying the calibration rho_{d,s}
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

        String calibrationKey = args.length > 3 && !args[3].isEmpty() ? args[3] : null;
        SpecificityReport report = new SpecificityReport(BASE_DIR, RESULTS_DIR, db);
        // The external rows need no Genestrip classification, so they can be redone on their own
        // once the full report has run. See SpecificityReport.writeExternalOnly.
        if (Boolean.getBoolean("gs.externalonly")) {
            report.writeExternalOnly(db, fqMapFile, reportKey, calibrationKey);
        } else {
            report.write(db, fqMapFile, reportKey, calibrationKey);
        }
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
