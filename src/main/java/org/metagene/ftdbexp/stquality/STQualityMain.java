package org.metagene.ftdbexp.stquality;

import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.finertree.FTProject;

import java.io.File;

/**
 * Command line entry point for the sequence-type level database quality of the Genestrip-FT paper.
 * <p>
 * It measures a database refined below the species rank against an external typing of the genomes it
 * was filled from, which is the only ground truth available there: the taxonomy has no rank beneath
 * the species to count candidates at, so the lineage takes its place. See {@link STQualityCountsGoal}
 * for the measure and {@link STGroundTruth} for what the typing has to look like.
 * <p>
 * Invoked from Maven:
 * <pre>
 *   mvn exec:exec@stquality -Dname=cdiff -Dstcsv=results/cdiff_mlst.csv
 * </pre>
 * or directly:
 * <pre>
 *   java ... org.metagene.ftdbexp.stquality.STQualityMain cdiff results/cdiff_mlst.csv [both|db|ftdb]
 * </pre>
 * Both databases are measured by default, since the point of the measure is the difference between
 * them. The results are written to {@code data/projects/<db>/csv/<db>_stquality.csv} and
 * {@code <db>_ftstquality.csv}; {@code run_exps.sh} copies whatever lands there into {@code results/}.
 * <p>
 * This reads the genomes again, exactly as {@code dbquality} does, so it takes about as long - and it
 * resolves a Genbank assembly the way Genestrip does, through {@code data/common/genbank}, not through
 * a search path.
 */
public class STQualityMain {
    private static final File BASE_DIR = new File("./data");

    /**
     * Runs the measure.
     *
     * @param args {@code <db> <sequence type csv> [both|db|ftdb]}
     * @throws Exception if the database or the typing cannot be read, or a goal fails
     */
    public static void main(String[] args) throws Exception {
        if (args.length < 2) {
            System.err.println("Usage: STQualityMain <db> <sequence type csv> [both|db|ftdb]");
            System.err.println("  <db>                  name of the database project under data/projects");
            System.err.println("  <sequence type csv>   the CSV written by bin/mlst_assemblies.sh, keyed by");
            System.err.println("                        leaf name (the fasta file a genome was read from)");
            System.err.println("  both|db|ftdb          which database(s) to measure; default 'both'");
            System.exit(1);
        }
        String db = args[0];
        File stCsv = new File(args[1]);
        String which = args.length > 2 && !args[2].isEmpty() ? args[2] : "both";
        boolean doDB = "both".equals(which) || "db".equals(which);
        boolean doFTDB = "both".equals(which) || "ftdb".equals(which);
        if (!doDB && !doFTDB) {
            System.err.println("Expected 'both', 'db' or 'ftdb' but got '" + which + "'.");
            System.exit(1);
        }
        if (!stCsv.exists()) {
            System.err.println("No such sequence type file: " + stCsv.getAbsolutePath());
            System.err.println("Produce it with bin/cdiff_eval.sh, which types one genome per leaf of the");
            System.err.println("database and keys its rows by the leaf's name.");
            System.exit(1);
        }

        STGroundTruth groundTruth = new STGroundTruth(stCsv);
        System.out.println("Read " + groundTruth.getTypedLeaves() + " typed genome(s) in "
                + groundTruth.getSTCount() + " distinct sequence type(s) from " + stCsv
                + "; " + groundTruth.getUntypedLeaves() + " listed without a type.");
        if (groundTruth.getTypedLeaves() == 0) {
            System.err.println("Nothing to measure against.");
            System.exit(1);
        }

        FTProject project = new FTProject(new GSCommon(BASE_DIR), db, null, null, null,
                null, null, null, null, null, null, false);
        project.initConfigParam(GSConfigKey.THREADS, -1);

        STQualityMaker maker = new STQualityMaker(project, groundTruth);
        try {
            // One at a time and in this order: each holds a database and the filter of its own pass,
            // and the unrefined one is the cheaper of the two to find a broken join with.
            if (doDB) {
                System.out.println("=== unrefined database ===");
                run(maker, false);
            }
            if (doFTDB) {
                System.out.println("=== refined database ===");
                run(maker, true);
            }
        } finally {
            maker.dumpAll();
        }
    }

    /**
     * Makes one of the two CSV goals and says where its file went, so that a run reports its own
     * output rather than leaving it to be looked for.
     *
     * @param maker   the maker holding the goals
     * @param refined whether the refined database is meant
     */
    private static void run(STQualityMaker maker, boolean refined) {
        STQualityCSVGoal goal = (STQualityCSVGoal) maker.getCSVGoal(refined);
        goal.make();
        for (File file : goal.getFiles()) {
            System.out.println("Wrote " + file.getAbsolutePath());
        }
    }
}
