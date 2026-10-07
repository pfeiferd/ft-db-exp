package org.metagene.ftdbexp.eval;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
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
    /**
     * The external classifiers whose per-read output is picked up from {@code results/kraken} if it
     * is there. The names are the ones bin/kraken_classify.sh writes into the file names.
     */
    static final String[] EXTERNAL_TOOLS = { "k2", "ku" };

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
        writeQuality(resultsDir, db, reportKey, byVariant);
        writeSimdata(resultsDir, db, reportKey, byVariant);
        writeExternal(db, fqMapFile, reportKey, scope, byVariant.get(Variant.UNREFINED));
        System.out.println("Wrote " + file);
        return file;
    }

    /**
     * Scores whatever external classifier left its per-read output under {@code results/kraken} and
     * writes it in the layout of this report's own CSVs, with the columns of the refined variant
     * left empty: an external tool has no refined counterpart, and pretending otherwise would put a
     * number where there is none.
     * <p>
     * Runs inside {@link #write} because the fastq keys of the runs above are what the output files
     * are looked up by. Each classifier is scored on its own genus-only subset, collected from its
     * own output by {@link AccuracyEvaluator#evaluateExternal}, which is what the columns of the
     * paper's table state. Until 2026-10-03 the subsets of the unrefined Genestrip run were handed
     * down here instead.
     *
     * @param db        the name of the database project
     * @param fqMapFile the fastq mapping file the runs above used
     * @param reportKey the key the result files are named after
     * @param scope     the scope passed to the runs above
     * @param unrefined the tallies of the unrefined run, whose fastq keys name the output files
     * @throws IOException if an output file cannot be read or a CSV cannot be written
     */
    private void writeExternal(String db, String fqMapFile, String reportKey, SmallTaxTree scope,
                               Map<String, AccuracyTally> unrefined) throws IOException {
        File krakenDir = new File(resultsDir, "kraken");
        if (!krakenDir.isDirectory()) {
            return;
        }
        for (String tool : EXTERNAL_TOOLS) {
            Map<String, File> outputs = new LinkedHashMap<String, File>();
            for (String fastqKey : unrefined.keySet()) {
                File out = new File(krakenDir, db + "_" + tool + "_" + fastqKey + ".tsv");
                if (out.isFile()) {
                    outputs.put(fastqKey, out);
                }
            }
            if (outputs.isEmpty()) {
                continue;
            }
            if (outputs.size() != unrefined.size()) {
                // A partial set would write a CSV with rows missing where the paper's table expects
                // one per read set, and the rows that are there would be read as the whole run.
                // Saying so beats writing a CSV nobody can trust.
                throw new IOException("Only " + outputs.size() + " of " + unrefined.size()
                        + " read sets have a " + tool + " output under " + krakenDir
                        + ". Classify them all or none.");
            }
            System.out.println("Evaluating " + tool + " on " + fqMapFile);
            // No rewind and no baseline handed over: evaluateExternal collects the classifier's own
            // subsets. A read-name mismatch between the tool's output and the simulation is caught
            // there too, by AccuracyEvaluator.failIfGroundTruthLost -- it used to be caught here, by
            // a subset that no read of the output fell into.
            Map<String, AccuracyTally> tallies = evaluator.evaluateExternal(db, fqMapFile,
                    Variant.UNREFINED.getLoadDbGoalKey(), scope, outputs);
            writeQualityOfExternal(resultsDir, db, tool, reportKey, tallies);
            writeSummaryOfExternal(resultsDir, db, tool, reportKey, tallies);
        }
    }

    /**
     * Writes {@code <db>_<report key>_summary.csv}: one row per fastq key, shaped so that the
     * paper's table can include it with {@code \csvreader} and nothing has to be copied by hand.
     * <p>
     * The detailed CSV beside it carries one row per database variant, which is the right shape for
     * the measurements but the wrong one for a table: the paper puts the unrefined and the refined
     * figure of a measure side by side, and {@code rho} is a quotient of differences <em>between</em>
     * those rows. LaTeX is poor at arithmetic across rows, so it is computed here instead.
     * <p>
     * The two gains themselves get no column: each is the difference of the two precisions printed
     * beside it, so a column would restate what the row already says.
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
            // The raw fastq key travels with the row so that a report keyed off this one can join
            // to it without depending on how the model happens to be spelled.
            ps.println("db;fastq key;model;reads;classified;genus only;genus only share"
                    + ";prec g u;prec g f"
                    + ";obs genus only;prec g ungated u;prec g ungated f;rho u;rho f;");
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
                // The calibration is a pair of factors on the levels, not one factor on the gain.
                // Each turns an ungated precision into the gated one it stands for, and the two are
                // kept apart because they answer to different effects. rho_u is the subset
                // discrepancy |R'_g| / |R_g| times a factor that is one when no read was stopped at
                // a foreign genus -- 0.9991 to 1.0000 on the viral sets, but 0.985 to 0.993 on the
                // tick ones, where the ungated average is small enough that the few dozen such reads
                // still shift it. rho_f carries that and the sibling reads a refinement rescues,
                // which is why rho_f exceeds rho_u on every read set here.
                //
                // Ratios of levels rather than of gains, for two reasons. They stay within [0, 1]
                // and so read as the shrinkage factors they are, where a ratio of two small
                // differences has no such bound and crosses one on "MiSeq". And they let a real
                // sample be given an estimated precision before and after, rather than only an
                // estimated gain -- the difference of the two estimates is that gain anyway.
                double rhoU = gu == 0 ? Double.NaN : pu / gu;
                double rhoF = gf == 0 ? Double.NaN : pf / gf;
                ps.print(db);
                ps.print(';');
                ps.print(fastqKey);
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
                ps.print(u.getObsGenusOnlyTotal());
                ps.print(';');
                ps.print(format(gu));
                ps.print(';');
                ps.print(format(gf));
                ps.print(';');
                ps.print(format(rhoU));
                ps.print(';');
                ps.print(format(rhoF));
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
        if ("iss_saliva".equals(fastqKey)) {
            return "saliva-like";
        }
        if ("iss_mngs".equals(fastqKey)) {
            return "mNGS-like";
        }
        // Anything else keeps its identity, which is what names the read set. For the tick-borne
        // data that is the sample -- tick1, tick2 and so on -- and mapping those to the simulator
        // would collapse eight distinct read sets into eight identical labels. Which simulator
        // produced them is a property of the database's block and belongs in the caption. Only the
        // spelling is normalised, through SampleNames, which SpecificityReport shares: it joins to
        // this report's rows by exactly this name, so the two must agree on it.
        return SampleNames.display(fastqKey, true);
    }

    /**
     * Writes {@code <db>_<report key>_quality.csv}: one row per read set rather than per variant, so
     * that a table can put the unrefined and the refined figure of a measure side by side.
     * <p>
     * Every measure is given for both variants, so that a table can show what the boolean counts do
     * -- namely almost nothing: a refinement inserts nodes between a genus and its species and moves
     * k-mers onto them, which changes how many species an answer leaves open, not whether the
     * answer is right. The candidate-weighted measures do see it, which is the contrast the paper's
     * two blocks are meant to make visible.
     *
     * @param resultsDir the directory to write to
     * @param db         the name of the database project
     * @param reportKey  short name used in the result file name
     * @param byVariant  the tallies of both variants, keyed by fastq key
     * @throws IOException if the file cannot be written
     */
    /**
     * Writes the quality CSV of an external classifier in the layout of {@link #writeQuality}, with
     * every column of the refined variant empty. The database name carries the tool, so that a row
     * of {@code cv_k2} cannot be mistaken for one of {@code cv}.
     *
     * @param resultsDir the directory the CSV goes into
     * @param db         the name of the database project
     * @param tool       the external tool, as in the output file names
     * @param reportKey  the key the result files are named after
     * @param tallies    the tool's tallies, keyed by fastq key
     * @throws IOException if the CSV cannot be written
     */
    private static void writeQualityOfExternal(File resultsDir, String db, String tool, String reportKey,
                                               Map<String, AccuracyTally> tallies) throws IOException {
        File file = new File(resultsDir, db + "_" + tool + "_" + reportKey + "_quality.csv");
        try (PrintStream ps = new PrintStream(new FileOutputStream(file), false, StandardCharsets.UTF_8.name())) {
            ps.println("db;fastq key;read set"
                    + ";prec genus u;prec genus f;recall genus u;recall genus f;f1 genus u;f1 genus f"
                    + ";prec species u;prec species f;recall species u;recall species f;f1 species u;f1 species f"
                    + ";prec cand u;prec cand f;recall cand u;recall cand f;f1 cand u;f1 cand f;");
            for (Map.Entry<String, AccuracyTally> e : tallies.entrySet()) {
                AccuracyTally t = e.getValue();
                ps.print(db + "_" + tool);
                ps.print(';');
                ps.print(e.getKey());
                ps.print(';');
                ps.print(displayModel(e.getKey()));
                ps.print(';');
                for (Rank rank : new Rank[] { Rank.GENUS, Rank.SPECIES }) {
                    ps.print(format(t.getPrecision(rank)));
                    ps.print(";;");
                    ps.print(format(t.getRecall(rank)));
                    ps.print(";;");
                    ps.print(format(t.getF1(rank)));
                    ps.print(";;");
                }
                ps.print(format(t.getSpeciesCandidatePrecision()));
                ps.print(";;");
                ps.print(format(t.getSpeciesCandidateRecall()));
                ps.print(";;");
                ps.print(format(t.getSpeciesCandidateF1()));
                ps.println(";;");
            }
        }
        System.out.println("Wrote " + file);
    }

    /**
     * Writes the summary CSV of an external classifier in the layout of {@link #writeSummary}, again
     * with the refined columns empty. {@code rho} is the tool's own ratio of its gated to its
     * ungated restricted precision, measured on the same reads as Genestrip's.
     *
     * @param resultsDir the directory the CSV goes into
     * @param db         the name of the database project
     * @param tool       the external tool, as in the output file names
     * @param reportKey  the key the result files are named after
     * @param tallies    the tool's tallies, keyed by fastq key
     * @throws IOException if the CSV cannot be written
     */
    private static void writeSummaryOfExternal(File resultsDir, String db, String tool, String reportKey,
                                               Map<String, AccuracyTally> tallies) throws IOException {
        File file = new File(resultsDir, db + "_" + tool + "_" + reportKey + "_summary.csv");
        try (PrintStream ps = new PrintStream(new FileOutputStream(file), false, StandardCharsets.UTF_8.name())) {
            ps.println("db;fastq key;model;reads;classified;genus only;genus only share"
                    + ";prec g u;prec g f"
                    + ";obs genus only;prec g ungated u;prec g ungated f;rho u;rho f;");
            for (Map.Entry<String, AccuracyTally> e : tallies.entrySet()) {
                AccuracyTally t = e.getValue();
                double p = t.getGenusOnlyPrecision();
                double g = t.getObsGenusOnlyUngatedPrecision();
                ps.print(db + "_" + tool);
                ps.print(';');
                ps.print(e.getKey());
                ps.print(';');
                ps.print(displayModel(e.getKey()));
                ps.print(';');
                ps.print(t.getTotal());
                ps.print(';');
                ps.print(t.getClassified());
                ps.print(';');
                ps.print(t.getGenusOnlyTotal());
                ps.print(';');
                // Over every read of the set, as writeSummary does it above and as the paper's
                // column "Share of |R|" says. It was over the classified reads here until
                // 2026-10-03, which went unnoticed while these rows carried the counts of the
                // unrefined Genestrip run: the denominators differed but the printed share was
                // simply a different number in a column nobody could compare.
                ps.print(format(t.getTotal() == 0 ? Double.NaN
                        : 100.0 * t.getGenusOnlyTotal() / t.getTotal()));
                ps.print(';');
                ps.print(format(p));
                ps.print(";;");
                ps.print(t.getObsGenusOnlyTotal());
                ps.print(';');
                ps.print(format(g));
                ps.print(";;");
                ps.print(format(g == 0 ? Double.NaN : p / g));
                ps.println(";;");
            }
        }
        System.out.println("Wrote " + file);
    }

    private static void writeQuality(File resultsDir, String db, String reportKey,
                                     Map<Variant, Map<String, AccuracyTally>> byVariant) throws IOException {
        File file = new File(resultsDir, db + "_" + reportKey + "_quality.csv");
        try (PrintStream ps = new PrintStream(new FileOutputStream(file), false, StandardCharsets.UTF_8.name())) {
            ps.println("db;fastq key;read set"
                    + ";prec genus u;prec genus f;recall genus u;recall genus f;f1 genus u;f1 genus f"
                    + ";prec species u;prec species f;recall species u;recall species f;f1 species u;f1 species f"
                    + ";prec cand u;prec cand f;recall cand u;recall cand f;f1 cand u;f1 cand f;");
            for (String fastqKey : byVariant.get(Variant.UNREFINED).keySet()) {
                AccuracyTally u = byVariant.get(Variant.UNREFINED).get(fastqKey);
                AccuracyTally f = byVariant.get(Variant.REFINED).get(fastqKey);
                if (u == null || f == null) {
                    continue;
                }
                ps.print(db);
                ps.print(';');
                ps.print(fastqKey);
                ps.print(';');
                ps.print(displayModel(fastqKey));
                ps.print(';');
                for (Rank rank : new Rank[] { Rank.GENUS, Rank.SPECIES }) {
                    ps.print(format(u.getPrecision(rank)));
                    ps.print(';');
                    ps.print(format(f.getPrecision(rank)));
                    ps.print(';');
                    ps.print(format(u.getRecall(rank)));
                    ps.print(';');
                    ps.print(format(f.getRecall(rank)));
                    ps.print(';');
                    ps.print(format(u.getF1(rank)));
                    ps.print(';');
                    ps.print(format(f.getF1(rank)));
                    ps.print(';');
                }
                ps.print(format(u.getSpeciesCandidatePrecision()));
                ps.print(';');
                ps.print(format(f.getSpeciesCandidatePrecision()));
                ps.print(';');
                ps.print(format(u.getSpeciesCandidateRecall()));
                ps.print(';');
                ps.print(format(f.getSpeciesCandidateRecall()));
                ps.print(';');
                ps.print(format(u.getSpeciesCandidateF1()));
                ps.print(';');
                ps.print(format(f.getSpeciesCandidateF1()));
                ps.println(';');
            }
        }
        System.out.println("Wrote " + file);
    }

    /**
     * Writes {@code <db>_<report key>_simdata.csv}: what the read sets are, rather than how well
     * they were classified.
     * <p>
     * The three columns describing the simulation --- read length, per-base error and the number of
     * reads generated --- cannot be recovered here. They are measured by {@code make_fastqs.sh} as
     * it produces each set and left in {@code <db>_simparams.csv}; by the time this runs, the
     * settings and the quality strings that produced them are either gone or, for NanoSim, known to
     * be meaningless. This joins them to the two counts that only the evaluation has: how many of
     * the generated reads fall within the database's scope, and how many were unresolved.
     * <p>
     * A read set with no row in {@code <db>_simparams.csv} still gets a row here, with the
     * simulation columns empty. That is deliberate: a set generated before the script recorded its
     * parameters should appear in the table as an incomplete row rather than vanish from it.
     *
     * @param resultsDir the directory to write to
     * @param db         the name of the database project
     * @param reportKey  short name used in the result file name
     * @param byVariant  the tallies of both variants, keyed by fastq key
     * @throws IOException if the file cannot be written
     */
    private static void writeSimdata(File resultsDir, String db, String reportKey,
                                     Map<Variant, Map<String, AccuracyTally>> byVariant) throws IOException {
        Map<String, String[]> params = readSimParams(resultsDir, db);
        File file = new File(resultsDir, db + "_" + reportKey + "_simdata.csv");
        try (PrintStream ps = new PrintStream(new FileOutputStream(file), false, StandardCharsets.UTF_8.name())) {
            ps.println("db;fastq key;read set;read length;per-base error;simulated"
                    + ";in scope;in scope share;unresolved;");
            for (String fastqKey : byVariant.get(Variant.UNREFINED).keySet()) {
                AccuracyTally u = byVariant.get(Variant.UNREFINED).get(fastqKey);
                if (u == null) {
                    continue;
                }
                String[] p = params.get(fastqKey);
                String length = p == null ? "" : p[0];
                String error = p == null ? "" : p[1];
                String simulated = p == null ? "" : p[2];
                ps.print(db);
                ps.print(';');
                ps.print(fastqKey);
                ps.print(';');
                ps.print(displayModel(fastqKey));
                ps.print(';');
                ps.print(length);
                ps.print(';');
                ps.print(error);
                ps.print(';');
                ps.print(simulated);
                ps.print(';');
                ps.print(u.getTotal());
                ps.print(';');
                // The share the scope keeps of what was generated. For a database covering the whole
                // category the reads were drawn from this is near 100 %; for one requesting a few
                // genera of a large category it is the fraction the table exists to show.
                double generated = parseOrNaN(simulated);
                // The scope is a subset of what was generated, so the share cannot exceed 100 %.
                // A value above it is never a measurement: it means the two numbers come from
                // different states of the pipeline, because they come from different places -- the
                // count from this run, `simulated' from <db>_simparams.csv, which make_fastqs.sh
                // writes one row per read set and never revisits. A set whose generation is skipped
                // keeps whatever row it once got, so resuming a batch past the fastq steps (for
                // instance run_all_exps.sh --from 9) leaves an older row in place. That is how
                // `saliva-like' came to report 200 % in the paper's Table "simdata": its row still
                // held InSilicoSeq's per-file count of 500,000, from before record_simparams began
                // scaling by the number of mates, while the evaluation counted both mates of every
                // pair. Say so here rather than let the quotient travel into a table.
                if (generated > 0 && u.getTotal() > generated) {
                    System.err.printf(
                            "WARNING: %s/%s reports %d reads in scope but only %s simulated (%.1f %%).%n"
                            + "         `in scope' cannot exceed `simulated'; the row for this read set in%n"
                            + "         %s_simparams.csv is stale. Delete it and re-run the read generation%n"
                            + "         for this set -- it skips the simulation but re-measures the file.%n",
                            db, fastqKey, u.getTotal(), simulated, 100.0 * u.getTotal() / generated, db);
                }
                ps.print(format(generated > 0 ? 100.0 * u.getTotal() / generated : Double.NaN));
                ps.print(';');
                ps.print(u.getUnresolved());
                ps.println(';');
            }
        }
        System.out.println("Wrote " + file);
    }

    /**
     * Reads {@code <db>_simparams.csv}, written by {@code make_fastqs.sh}.
     *
     * @param resultsDir the directory holding it
     * @param db         the name of the database project
     * @return fastq key to {read length, per-base error, reads generated}; empty if absent
     */
    private static Map<String, String[]> readSimParams(File resultsDir, String db) {
        Map<String, String[]> byKey = new HashMap<String, String[]>();
        File file = new File(resultsDir, db + "_simparams.csv");
        if (!file.exists()) {
            System.out.println("No " + file + " - the simulation columns stay empty. "
                    + "Run bin/make_fastqs.sh to record them.");
            return byKey;
        }
        try (BufferedReader in = new BufferedReader(new InputStreamReader(
                new FileInputStream(file), StandardCharsets.UTF_8))) {
            String line = in.readLine();
            while ((line = in.readLine()) != null) {
                String[] c = line.split(";", -1);
                if (c.length >= 4 && !c[0].trim().isEmpty()) {
                    byKey.put(c[0].trim(), new String[]{c[1].trim(), c[2].trim(), c[3].trim()});
                }
            }
        } catch (IOException e) {
            System.out.println("Cannot read " + file + " (" + e.getMessage() + ") - columns stay empty.");
        }
        return byKey;
    }

    private static double parseOrNaN(String value) {
        try {
            return Double.parseDouble(value);
        } catch (RuntimeException e) {
            return Double.NaN;
        }
    }

    /**
     * Writes the CSV header.
     *
     * @param ps the stream to write to
     */
    private static void writeHeader(PrintStream ps) {
        ps.print("fastq key;read set;variant;classified;total;unresolved");
        for (Rank rank : REPORTED_RANKS) {
            String name = rank.getName();
            ps.print(";correct " + name + ";precision " + name + ";recall " + name + ";f1 " + name);
        }
        // The candidate-weighted species measures: a classification leaving n species in question
        // counts as 1/n of a hit, so narrowing the species down pays off even short of pinning it.
        ps.print(";species score;precision species cand;recall species cand;f1 species cand");
        // Restricted to the reads the unrefined database left at their genus: the only ones a
        // refinement can improve on. The increase in "genus only precision species cand" from the
        // unrefined to the refined variant is the gain where a gain was possible.
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
        ps.print(displayModel(fastqKey));
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
