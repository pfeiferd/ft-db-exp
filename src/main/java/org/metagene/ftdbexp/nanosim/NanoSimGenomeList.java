package org.metagene.ftdbexp.nanosim;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.goals.refseq.ExtractRefSeqCSVGoal;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;

/**
 * Writes the genome list that NanoSim reads via its {@code -gl} option.
 * <p>
 * The list has one tab-separated line per reference genome, {@code <label>\t<path>}. The label
 * becomes a prefix of every simulated read's identifier, which is what later allows the read to be
 * traced back to the taxon it was generated from. Following the experiments of the first Genestrip
 * paper the label is {@code <taxid>x<index>}, so that the taxon is visible in the read name even
 * without consulting the accession map.
 * <p>
 * The genomes themselves are the per-accession fasta files that Genestrip's goal
 * {@code extractrefseqfasta} extracts from the RefSeq archives; this class derives their names from
 * the accompanying CSV of the goal {@code extractrefseqcsv}. Note that NanoSim does not accept
 * gzipped fasta files.
 * <p>
 * Very small fasta files are left out: NanoSim cannot draw reads of a realistic length from them and
 * they slow the simulation down disproportionately. Unlike the original experiment, which dropped
 * them silently, this class reports how many entries it skipped -- a genome missing from the list
 * simply never contributes reads, which quietly changes what the simulation covers.
 */
public class NanoSimGenomeList {
    /** CSV dialect of the {@code extractrefseqcsv} goal. */
    private static final CSVFormat CSV_FORMAT = CSVFormat.DEFAULT.builder()
            .setQuote(null).setCommentMarker('#').setDelimiter(';').setRecordSeparator('\n').build();

    /** Smallest fasta file still included in the genome list, in bytes. */
    private static final long MIN_FASTA_SIZE = 1024;

    private final File baseDir;

    /**
     * Creates the writer.
     *
     * @param baseDir the Genestrip base directory, i.e. the one holding {@code common} and
     *                {@code projects}
     */
    public NanoSimGenomeList(File baseDir) {
        this.baseDir = baseDir;
    }

    /**
     * Runs the {@code extractrefseqcsv} goal for the given database and turns its result into a
     * NanoSim genome list.
     *
     * @param db        the name of the database project
     * @param fastaDir  the directory holding the extracted per-accession {@code .fa} files
     * @param outFile   the genome list to write
     * @return the number of genomes written to the list
     * @throws IOException if the goal cannot be run or the list cannot be written
     */
    public int write(String db, File fastaDir, File outFile) throws IOException {
        GSProject project = new GSProject(new GSCommon(baseDir), db, null, null, null, null, null,
                null, null, null, null, false);
        project.initConfigParam(GSConfigKey.THREADS, -1);

        GSMaker maker = new GSMaker(project);
        File csvFile;
        try {
            ExtractRefSeqCSVGoal<?> goal =
                    (ExtractRefSeqCSVGoal<?>) maker.getGoal(GSGoalKey.EXTRACT_REFSEQ_CSV);
            goal.make();
            csvFile = goal.getFile();
        } finally {
            maker.dumpAll();
        }

        File parent = outFile.getParentFile();
        if (parent != null && !parent.exists() && !parent.mkdirs()) {
            throw new IOException("Cannot create directory " + parent);
        }

        int written = 0;
        int missing = 0;
        int tooSmall = 0;
        try (CSVParser parser = CSV_FORMAT.parse(
                new InputStreamReader(new FileInputStream(csvFile), StandardCharsets.UTF_8));
             PrintStream out = new PrintStream(new FileOutputStream(outFile), false,
                     StandardCharsets.UTF_8.name())) {
            int index = 0;
            for (CSVRecord record : parser) {
                index++;
                if (index == 1 || record.size() < 2) {
                    continue; // header line
                }
                String accession = record.get(0).trim();
                String taxId = record.get(1).trim();
                File fasta = new File(fastaDir, accession + ".fa");
                if (!fasta.isFile()) {
                    missing++;
                } else if (fasta.length() <= MIN_FASTA_SIZE) {
                    tooSmall++;
                } else {
                    out.print(taxId);
                    out.print('x');
                    out.print(index);
                    out.print('\t');
                    out.print(fasta.getPath());
                    out.println();
                    written++;
                }
            }
        }

        System.out.println("Wrote " + outFile + " with " + written + " genomes.");
        if (tooSmall > 0) {
            System.out.println("  skipped " + tooSmall + " fasta files of at most "
                    + MIN_FASTA_SIZE + " bytes - too short for realistic reads.");
        }
        if (missing > 0) {
            System.out.println("  skipped " + missing + " accessions without a fasta file in "
                    + fastaDir + " - run the goal 'extractrefseqfasta' if that is unexpected.");
        }
        if (written == 0) {
            System.out.println("  WARNING: the list is empty, NanoSim has nothing to simulate from.");
        }
        return written;
    }

    /**
     * Writes a NanoSim genome list from the command line.
     *
     * @param args the database project name and, optionally, the fasta directory and the output
     *             file; both default to the project's own directories
     * @throws IOException if the list cannot be written
     */
    public static void main(String[] args) throws IOException {
        if (args.length < 1) {
            System.err.println("Usage: NanoSimGenomeList <db> [<fasta dir>] [<out file>]");
            System.exit(1);
        }
        String db = args[0];
        File baseDir = new File("./data");
        File fastaDir = args.length > 1 && !args[1].isEmpty()
                ? new File(args[1]) : new File(baseDir, "projects/" + db + "/fasta");
        File outFile = args.length > 2 && !args[2].isEmpty()
                ? new File(args[2]) : new File(baseDir, "projects/" + db + "/csv/" + db + "_nanosim.tsv");
        new NanoSimGenomeList(baseDir).write(db, fastaDir, outFile);
    }
}
