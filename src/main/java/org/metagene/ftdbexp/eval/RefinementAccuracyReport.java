package org.metagene.ftdbexp.eval;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Compares the read classification quality of an unrefined Genestrip database against its refined
 * counterpart and writes the result as a CSV file for inclusion in the paper.
 * <p>
 * Both variants are evaluated over the very same simulated fastq files and against the very same
 * ground truth, so the difference between the two rows of a fastq key is attributable to the
 * refinement alone.
 */
public class RefinementAccuracyReport {
    /** The database variants compared, in the order they appear in the report. */
    public enum Variant {
        /** The database as produced by the LCA update, i.e. before any refinement. */
        UNREFINED("unrefined", GSGoalKey.MATCHRES, GSGoalKey.LOAD_DB),
        /** The database after the refinement of the goal {@code ftupdatedb}. */
        REFINED("refined", FTGoalKey.FTMATCHRES, FTGoalKey.LOAD_FTDB);

        private final String label;
        private final GoalKey matchGoalKey;
        private final GoalKey loadDbGoalKey;

        Variant(String label, GoalKey matchGoalKey, GoalKey loadDbGoalKey) {
            this.label = label;
            this.matchGoalKey = matchGoalKey;
            this.loadDbGoalKey = loadDbGoalKey;
        }

        /**
         * Returns the name used for this variant in the CSV file.
         *
         * @return the CSV label
         */
        public String getLabel() {
            return label;
        }

        /**
         * Returns the matching goal that evaluates this variant.
         *
         * @return the goal key of the variant's {@code matchres} goal
         */
        public GoalKey getMatchGoalKey() {
            return matchGoalKey;
        }

        /**
         * Returns the goal loading this variant's database, whose taxonomy determines how many
         * species a classification leaves in question.
         *
         * @return the goal key of the variant's database loading goal
         */
        public GoalKey getLoadDbGoalKey() {
            return loadDbGoalKey;
        }
    }

    /** The ranks reported, in the order their columns appear. */
    private static final Rank[] REPORTED_RANKS = {Rank.GENUS, Rank.SPECIES};

    private final AccuracyEvaluator evaluator;
    private final File resultsDir;

    /**
     * Creates the report writer.
     *
     * @param baseDir    the Genestrip base directory
     * @param resultsDir the directory the CSV file is written to
     * @param db         the name of the database project
     * @param simulator  the simulator that produced the reads
     * @throws IOException if the taxonomy or the accession map cannot be read
     */
    public RefinementAccuracyReport(File baseDir, File resultsDir, String db, Simulator simulator)
            throws IOException {
        this.evaluator = new AccuracyEvaluator(baseDir, db, simulator);
        this.resultsDir = resultsDir;
    }

    /**
     * Evaluates both database variants and writes {@code <db>_<report key>_accuracy.csv} to the
     * results directory.
     *
     * @param db        the name of the database project
     * @param fqMapFile the fastq mapping file, relative to the project's {@code txt} directory
     * @param reportKey a short name identifying the fastq collection, used in the file name
     * @param scope     restricts the reads counted towards recall, may be {@code null}
     * @return the file that was written
     * @throws IOException if the file cannot be written
     */
    public File write(String db, String fqMapFile, String reportKey, SmallTaxTree scope) throws IOException {
        Map<Variant, Map<String, AccuracyTally>> byVariant =
                new LinkedHashMap<Variant, Map<String, AccuracyTally>>();
        // The unrefined run fills the baseline of reads it left at their genus, the refined one is
        // then measured on exactly those. Hence the order of Variant.values() matters here.
        GenusOnlyBaseline baseline = new GenusOnlyBaseline();
        // The ground-truth-free substitute for that subset, fixed by the same run for the same
        // reason: both variants must be scored on one and the same set of reads.
        GenusOnlyBaseline obsBaseline = new GenusOnlyBaseline();
        for (Variant variant : Variant.values()) {
            System.out.println("Evaluating " + variant.getLabel() + " database " + db + " on " + fqMapFile);
            byVariant.put(variant, evaluator.evaluate(db, fqMapFile, variant.getMatchGoalKey(),
                    variant.getLoadDbGoalKey(), scope, baseline, obsBaseline,
                    variant == Variant.UNREFINED, false));
        }

        if (!resultsDir.exists() && !resultsDir.mkdirs()) {
            throw new IOException("Cannot create results directory " + resultsDir);
        }
        File file = new File(resultsDir, db + "_" + reportKey + "_accuracy.csv");
        try (PrintStream ps = new PrintStream(new FileOutputStream(file), false, StandardCharsets.UTF_8.name())) {
            writeHeader(ps);
            for (String fastqKey : byVariant.get(Variant.UNREFINED).keySet()) {
                for (Variant variant : Variant.values()) {
                    AccuracyTally tally = byVariant.get(variant).get(fastqKey);
                    if (tally != null) {
                        writeRow(ps, fastqKey, variant, tally);
                    }
                }
            }
        }
        writeSummary(resultsDir, db, reportKey, byVariant);
        System.out.println("Wrote " + file);
        return file;
    }

    /**
     * Writes {@code <db>_<report key>_summary.csv}: one row per fastq key, shaped so that the
     * paper's table can include it with {@code \csvreader} and nothing has to be copied by hand.
     * <p>
     * The detailed CSV beside it carries one row per database variant, which is the right shape for
     * the measurements but the wrong one for a table: the quantities the paper reports --- the gain
     * {@code delta}, its ground-truth-free counterpart {@code delta'} and their ratio {@code rho} ---
     * are differences and a quotient <em>between</em> those rows. LaTeX is poor at arithmetic across
     * rows, so they are computed here instead.
     *
     * @param db        the name of the database project
     * @param reportKey the report key, used in the file name
     * @param byVariant the tallies of both variants, keyed by fastq key
     * @throws IOException if the file cannot be written
     */
    private static void writeSummary(File resultsDir, String db, String reportKey,
                                     Map<Variant, Map<String, AccuracyTally>> byVariant) throws IOException {
        File file = new File(resultsDir, db + "_" + reportKey + "_summary.csv");
        try (PrintStream ps = new PrintStream(new FileOutputStream(file), false, StandardCharsets.UTF_8.name())) {
            ps.println("db;model;reads;classified;genus only;genus only share"
                    + ";prec g u;prec g f;delta"
                    + ";obs genus only;prec g ungated u;prec g ungated f;delta ungated;rho;");
            for (String fastqKey : byVariant.get(Variant.UNREFINED).keySet()) {
                AccuracyTally u = byVariant.get(Variant.UNREFINED).get(fastqKey);
                AccuracyTally f = byVariant.get(Variant.REFINED).get(fastqKey);
                if (u == null || f == null) {
                    continue;
                }
                double pu = u.getGenusOnlyPrecision();
                double pf = f.getGenusOnlyPrecision();
                double gu = u.getObsGenusOnlyUngatedPrecision();
                double gf = f.getObsGenusOnlyUngatedPrecision();
                double delta = pf - pu;
                double deltaUngated = gf - gu;
                ps.print(db);
                ps.print(';');
                ps.print(displayModel(fastqKey));
                ps.print(';');
                ps.print(u.getTotal());
                ps.print(';');
                ps.print(u.getClassified());
                ps.print(';');
                ps.print(u.getGenusOnlyTotal());
                ps.print(';');
                ps.print(format(u.getTotal() == 0 ? Double.NaN
                        : 100.0 * u.getGenusOnlyTotal() / u.getTotal()));
                ps.print(';');
                ps.print(format(pu));
                ps.print(';');
                ps.print(format(pf));
                ps.print(';');
                ps.print(format(delta));
                ps.print(';');
                ps.print(u.getObsGenusOnlyTotal());
                ps.print(';');
                ps.print(format(gu));
                ps.print(';');
                ps.print(format(gf));
                ps.print(';');
                ps.print(format(deltaUngated));
                ps.print(';');
                ps.print(format(deltaUngated == 0 ? Double.NaN : delta / deltaUngated));
                ps.println(';');
            }
        }
        System.out.println("Wrote " + file);
    }

    /**
     * Turns a fastq key into the label the paper prints for it, so that the CSV can be included
     * without a mapping on the LaTeX side. The label names the <em>read set</em>: for the
     * InSilicoSeq projects one set per error model, for the tick-borne data one per sample.
     *
     * @param fastqKey the key as it appears in the fastq mapping file
     * @return the label to print
     */
    private static String displayModel(String fastqKey) {
        if ("iss_miseq".equals(fastqKey)) {
            return "MiSeq";
        }
        if ("iss_hiseq".equals(fastqKey)) {
            return "HiSeq";
        }
        if ("iss_perfect".equals(fastqKey)) {
            return "error-free";
        }
        if ("iss_nanopore".equals(fastqKey)) {
            return "Nanopore";
        }
        if ("nanosim".equals(fastqKey)) {
            return "NanoSim";
        }
        // Anything else keeps its key, which is what names the read set. For the tick-borne data
        // that is the sample -- tick1, tick2 and so on -- and mapping those to the simulator would
        // collapse eight distinct read sets into eight identical labels. Which simulator produced
        // them is a property of the database's block in the table and belongs in its caption.
        return fastqKey;
    }

    /**
     * Writes the CSV header.
     *
     * @param ps the stream to write to
     */
    private static void writeHeader(PrintStream ps) {
        ps.print("fastq key;variant;classified;total;unresolved");
        for (Rank rank : REPORTED_RANKS) {
            String name = rank.getName();
            ps.print(";correct " + name + ";precision " + name + ";recall " + name + ";f1 " + name);
        }
        // The candidate-weighted species measures: a classification leaving n species in question
        // counts as 1/n of a hit, so narrowing the species down pays off even short of pinning it.
        ps.print(";species score;precision species cand;recall species cand;f1 species cand");
        // Restricted to the reads the unrefined database left at their genus: the only ones a
        // refinement can improve on. The delta between the two variants of "genus only precision
        // species cand" is the gain where a gain was possible.
        ps.print(";genus only;genus only score;genus only precision species cand"
                + ";genus only species share;genus only zero scoring");
        // The ungated measure, prec'_g of the paper. It differs from the gated one in *both* of its
        // ingredients: the summand drops the test that the read's true species is still in question
        // at the assigned node, and the subset is the observable one -- the reads the unrefined
        // database left at a genus, right or wrong -- because R_g itself is defined through the
        // ground truth and a real sample supplies none. Neither ingredient consults sigma(r), which
        // is what makes these two columns obtainable from a real fastq file. They are reported
        // alongside the gated ones so that the ratio rho_d between the two gains can be calibrated
        // where the truth *is* known.
        ps.println(";obs genus only;obs genus only ungated precision;");
    }

    /**
     * Writes one result row.
     *
     * @param ps       the stream to write to
     * @param fastqKey the key of the fastq file the row refers to
     * @param variant  the database variant the row refers to
     * @param tally    the counts to report
     */
    private static void writeRow(PrintStream ps, String fastqKey, Variant variant, AccuracyTally tally) {
        ps.print(fastqKey);
        ps.print(';');
        ps.print(variant.getLabel());
        ps.print(';');
        ps.print(tally.getClassified());
        ps.print(';');
        ps.print(tally.getTotal());
        ps.print(';');
        ps.print(tally.getUnresolved());
        for (Rank rank : REPORTED_RANKS) {
            ps.print(';');
            ps.print(tally.getCorrect(rank));
            ps.print(';');
            ps.print(format(tally.getPrecision(rank)));
            ps.print(';');
            ps.print(format(tally.getRecall(rank)));
            ps.print(';');
            ps.print(format(tally.getF1(rank)));
        }
        ps.print(';');
        ps.print(format(tally.getSpeciesCandidateScore()));
        ps.print(';');
        ps.print(format(tally.getSpeciesCandidatePrecision()));
        ps.print(';');
        ps.print(format(tally.getSpeciesCandidateRecall()));
        ps.print(';');
        ps.print(format(tally.getSpeciesCandidateF1()));
        ps.print(';');
        ps.print(tally.getGenusOnlyTotal());
        ps.print(';');
        ps.print(format(tally.getGenusOnlyScore()));
        ps.print(';');
        ps.print(format(tally.getGenusOnlyPrecision()));
        ps.print(';');
        ps.print(format(tally.getGenusOnlySpeciesShare()));
        ps.print(';');
        ps.print(tally.getGenusOnlyZeroScoring());
        ps.print(';');
        ps.print(tally.getObsGenusOnlyTotal());
        ps.print(';');
        ps.print(format(tally.getObsGenusOnlyUngatedPrecision()));
        ps.println(';');
    }

    /**
     * Formats a metric, writing an empty field for an undefined value so that the CSV does not carry
     * a locale-dependent NaN symbol into the paper.
     *
     * @param value the metric to format
     * @return the formatted value, or the empty string if it is undefined
     */
    private static String format(double value) {
        return Double.isNaN(value) ? "" : Double.toString(value);
    }
}
