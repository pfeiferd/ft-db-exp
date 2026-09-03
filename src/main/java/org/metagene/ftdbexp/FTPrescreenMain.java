package org.metagene.ftdbexp;

import org.metagene.ftdbexp.prescreen.PrescreenReport;

import java.io.File;

/**
 * Runs the branching-degree pre-screen of {@link PrescreenReport} for one database.
 * <p>
 * It reads the histogram the goal {@code branchhistocsv} writes and answers one question: is there
 * enough k-mer mass at a low branching degree for a refinement to be worth building? Costing that
 * answer an unrefined database rather than a refined one is the whole point.
 */
public class FTPrescreenMain {
    private static final File RESULTS_DIR = new File("./results");

    /**
     * @param args the database project name
     * @throws Exception if the histogram cannot be read or the report cannot be written
     */
    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.err.println("Usage: FTPrescreenMain <db>");
            System.err.println("  <db>  name of the database project under data/projects");
            System.err.println();
            System.err.println("Reads results/<db>_branchhistocsv.csv and writes results/<db>_prescreen.csv.");
            System.exit(1);
        }
        new PrescreenReport(RESULTS_DIR).write(args[0]);
    }
}
