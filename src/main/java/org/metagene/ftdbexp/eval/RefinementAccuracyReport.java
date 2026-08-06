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
                    variant == Variant.UNREFINED));
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
        System.out.println("Wrote " + file);
        return file;
    }

    /**
     * Writes the CSV header.
     *
     * @param ps the stream to write to
     */
    private void writeHeader(PrintStream ps) {
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
        ps.print(";genus only;genus only score;genus only precision species cand;genus only species share");
        // The ungated counterpart of the same measure -- the reciprocal candidate count without the
        // test that the read's true species is still in question -- together with the observable
        // substitute for the genus-only subset. Neither looks at the ground truth, so both can be
        // obtained from a real fastq file; the columns are reported here so that the ratio between
        // the gated and the ungated gain can be calibrated on data where the truth *is* known.
        ps.println(";genus only ungated precision;genus only gate missed"
                + ";obs genus only;obs genus only precision;obs genus only ungated precision"
                + ";obs genus only also true;");
    }

    /**
     * Writes one result row.
     *
     * @param ps       the stream to write to
     * @param fastqKey the key of the fastq file the row refers to
     * @param variant  the database variant the row refers to
     * @param tally    the counts to report
     */
    private void writeRow(PrintStream ps, String fastqKey, Variant variant, AccuracyTally tally) {
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
        ps.print(format(tally.getGenusOnlyUngatedPrecision()));
        ps.print(';');
        ps.print(tally.getGenusOnlyGateMissed());
        ps.print(';');
        ps.print(tally.getObsGenusOnlyTotal());
        ps.print(';');
        ps.print(format(tally.getObsGenusOnlyPrecision()));
        ps.print(';');
        ps.print(format(tally.getObsGenusOnlyUngatedPrecision()));
        ps.print(';');
        ps.print(tally.getObsGenusOnlyAlsoTrue());
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
