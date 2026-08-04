package org.metagene.ftdbexp;

import org.metagene.ftdbexp.eval.RefinementAccuracyReport;
import org.metagene.ftdbexp.eval.Simulator;

import java.io.File;

/**
 * Command line entry point for the classification-quality experiments of the Genestrip-FT paper.
 * <p>
 * It is invoked from Maven, e.g.
 * <pre>
 *   mvn exec:exec@accuracy -Dname=viral -Dfqmap=viral_sim.txt -Dreportkey=iss
 * </pre>
 * and writes {@code results/&lt;db&gt;_&lt;report key&gt;_accuracy.csv}.
 */
public class FTExpMain {
    private static final File BASE_DIR = new File("./data");
    private static final File RESULTS_DIR = new File("./results");

    /**
     * Runs the accuracy comparison between the unrefined and the refined database.
     *
     * @param args the database project name, the fastq mapping file and, optionally, a report key
     *             used in the output file name and the simulator that produced the reads; the report
     *             key defaults to the mapping file name without its extension, the simulator to ISS
     * @throws Exception if the databases cannot be read or the report cannot be written
     */
    public static void main(String[] args) throws Exception {
        if (args.length < 2) {
            System.err.println("Usage: FTExpMain <db> <fastq map file> [<report key>] [<simulator>]");
            System.err.println("  <db>             name of the database project under data/projects");
            System.err.println("  <fastq map file> mapping file under data/projects/<db>/txt");
            System.err.println("  <report key>     short name used in the result file name");
            System.err.println("  <simulator>      ISS (default) or NANOSIM - how the reads encode their origin");
            System.exit(1);
        }
        String db = args[0];
        String fqMapFile = args[1];
        String reportKey = args.length > 2 && !args[2].isEmpty() ? args[2] : stripExtension(fqMapFile);
        Simulator simulator = args.length > 3 && !args[3].isEmpty()
                ? Simulator.parse(args[3]) : Simulator.ISS;

        RefinementAccuracyReport report =
                new RefinementAccuracyReport(BASE_DIR, RESULTS_DIR, db, simulator);
        // A null scope means: count exactly the reads whose organism the database covers.
        report.write(db, fqMapFile, reportKey, null);
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
