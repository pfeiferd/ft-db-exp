package org.metagene.ftdbexp.eval;

import org.metagene.ftdbexp.eval.RefinementAccuracyReport.Variant;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Measures what a refinement achieves on fastq files whose ground truth is <em>not</em> known, i.e.
 * on real sequencing data.
 * <p>
 * Everything the {@link RefinementAccuracyReport} reports rests on the species a read actually stems
 * from. A real sample does not supply it, so precision and recall are simply not available. What
 * remains available is how far a classification narrows the species down: the reciprocal of the
 * number of species still in question at the assigned node, averaged over the reads the unrefined
 * database left at a genus. The difference between the two database variants is the
 * increase in the ungated precision of the paper's Section "Estimating the gain without ground
 * truth". It is not a bound on the increase in the gated precision in either direction --- the two are averages of different
 * summands over different subsets --- but each ungated precision, carried by the calibration factor
 * {@code rho_u} or {@code rho_f} determined on simulated reads of the same database, estimates the
 * gated precision it stands for. Those two estimates, and the gain between them, are what this
 * report exists for.
 * <p>
 * The subset is fixed by the unrefined run and looked up by the refined one, exactly as for the
 * simulated data, so that both variants are scored on one and the same set of reads.
 */
public class SpecificityReport {
    private final AccuracyEvaluator evaluator;
    private final File resultsDir;

    /**
     * Creates the report writer.
     *
     * @param baseDir    the Genestrip base directory holding {@code common} and {@code projects}
     * @param resultsDir the directory the CSV is written to
     * @param db         the name of the database project
     * @throws IOException if the taxonomy cannot be read
     */
    public SpecificityReport(File baseDir, File resultsDir, String db) throws IOException {
        // The simulator only determines how a read's origin would be read from its identifier, and
        // nothing does that here. ISS is passed as an arbitrary placeholder; the ground-truth-free
        // evaluation never calls the resolver.
        this.evaluator = new AccuracyEvaluator(baseDir, db, Simulator.ISS);
        this.resultsDir = resultsDir;
    }

    /**
     * Runs both database variants over the fastq files of the mapping file and writes one row per
     * fastq key to {@code <db>_<report key>_specificity.csv}.
     *
     * @param db        the name of the database project
     * @param fqMapFile the fastq mapping file
     * @param reportKey      short name used in the result file name
     * @param calibrationKey the report key of the simulated run whose {@code rho_{d,s}} calibrates
     *                       this one, or {@code null} for none. The join is by fastq key, which is
     *                       exactly right for the ticks: the simulation named {@code tickN} was
     *                       trained on the real sample named {@code tickN}, so each estimate uses
     *                       the calibration derived from its own sample.
     * @return the file that was written
     * @throws IOException if the databases cannot be read or the file cannot be written
     */
    public File write(String db, String fqMapFile, String reportKey, String calibrationKey)
            throws IOException {
        Map<Variant, Map<String, AccuracyTally>> byVariant =
                new LinkedHashMap<Variant, Map<String, AccuracyTally>>();
        GenusOnlyBaseline unused = new GenusOnlyBaseline();
        GenusOnlyBaseline obsBaseline = new GenusOnlyBaseline();
        for (Variant variant : Variant.values()) {
            System.out.println("Evaluating " + variant.getLabel() + " database " + db + " on " + fqMapFile);
            byVariant.put(variant, evaluator.evaluate(db, fqMapFile, variant.getMatchGoalKey(),
                    variant.getLoadDbGoalKey(), (SmallTaxTree) null, unused, obsBaseline,
                    variant == Variant.UNREFINED, true));
        }

        if (!resultsDir.exists() && !resultsDir.mkdirs()) {
            throw new IOException("Cannot create results directory " + resultsDir);
        }
        Map<String, double[]> calibration = readCalibration(db, calibrationKey);
        File file = new File(resultsDir, db + "_" + reportKey + "_specificity.csv");
        try (PrintStream ps = new PrintStream(new FileOutputStream(file), false, StandardCharsets.UTF_8.name())) {
            ps.println("fastq key;sample;reads;classified unrefined;classified refined"
                    + ";obs genus only;obs genus only share"
                    + ";ungated precision unrefined;ungated precision refined"
                    + ";rho u;rho f;est prec g u;est prec g f;obs genus only reported;");
            for (String fastqKey : byVariant.get(Variant.UNREFINED).keySet()) {
                AccuracyTally u = byVariant.get(Variant.UNREFINED).get(fastqKey);
                AccuracyTally f = byVariant.get(Variant.REFINED).get(fastqKey);
                if (u == null || f == null) {
                    continue;
                }
                writeRow(ps, fastqKey, u, f, calibration);
            }
        }
        writeExternal(db, fqMapFile, reportKey, calibrationKey, obsBaseline,
                byVariant.get(Variant.UNREFINED));
        System.out.println("Wrote " + file);
        return file;
    }

    /**
     * Writes one row of a specificity CSV.
     * <p>
     * Both tallies of a Genestrip run go in here. An external classifier has no refined variant
     * and three cells it cannot fill from its own tally, so it has {@link #writeExternalRow} of its
     * own; the {@code null} guards on {@code f} below are what the two shapes still share.
     *
     * @param ps          the stream to write to
     * @param fastqKey    the key of the sample
     * @param u           the tally of the unrefined run
     * @param f           the tally of the refined run, or {@code null} if there is none
     * @param calibration the rho values, keyed by fastq key
     */
    private static void writeRow(PrintStream ps, String fastqKey, AccuracyTally u,
                                 AccuracyTally f, Map<String, double[]> calibration) {
        double pu = u.getObsGenusOnlyUngatedPrecision();
        double pf = f == null ? Double.NaN : f.getObsGenusOnlyUngatedPrecision();
        ps.print(fastqKey);
        ps.print(';');
        ps.print(SampleNames.display(fastqKey, false));
        ps.print(';');
        ps.print(u.getTotal());
        ps.print(';');
        ps.print(u.getClassified());
        ps.print(';');
        ps.print(f == null ? "" : String.valueOf(f.getClassified()));
        ps.print(';');
        // The observable subset, never R_g: without sigma(r) there is no telling whether
        // the genus a read was left at was the right one, so getGenusOnlyTotal() stays zero
        // on this path by construction.
        ps.print(u.getObsGenusOnlyTotal());
        ps.print(';');
        ps.print(format(u.getClassified() == 0 ? Double.NaN
                : 100.0 * u.getObsGenusOnlyTotal() / u.getClassified()));
        ps.print(';');
        ps.print(format(pu));
        ps.print(';');
        ps.print(f == null ? "" : format(pf));
        ps.print(';');
        // Two ways a calibration applies. Per sample, when the simulated run was trained
        // on the very sample being scored -- the ticks, where both are keyed tickN. Or
        // wholesale, when the calibration run models the sample's parameters rather than
        // one sample: the saliva-like read set is a single row standing for all three
        // saliva runs, whose keys are SRA accessions and match nothing. A calibration of
        // exactly one row is therefore taken to apply to every sample; more than one row
        // means the rows are per sample and only a key match will do.
        double[] rho = calibration.get(fastqKey);
        if (rho == null && calibration.size() == 1) {
            rho = calibration.values().iterator().next();
        }
        ps.print(rho == null ? "" : format(rho[0]));
        ps.print(';');
        ps.print(rho == null || Double.isNaN(rho[1]) ? "" : format(rho[1]));
        ps.print(';');
        // The estimated gated precision before and after, each level carrying its own
        // factor. Their difference is the estimated precision gain; it gets no column, for
        // the same reason the measured gains get none.
        ps.print(rho == null ? "" : format(rho[0] * pu));
        ps.print(';');
        ps.print(rho == null || Double.isNaN(rho[1]) || f == null ? ""
                : format(rho[1] * pf));
        ps.print(';');
        ps.println(';');
    }

    /**
     * Scores whatever external classifier left per-read output for these samples under
     * {@code results/kraken} and writes it in the layout of this report's own CSV.
     * <p>
     * The subset is the one the unrefined Genestrip run collected, and that is the point of this
     * pass. The question it answers is where another tool places the very reads the unrefined
     * database could not resolve, so both must be scored on one and the same set of reads. This is
     * the opposite choice from {@link RefinementAccuracyReport}, whose external rows carry subsets
     * of their own, and the two differ because the questions differ: there the column states what a
     * classifier leaves at a genus, here it states how far another tool gets on a fixed set of hard
     * reads.
     * <p>
     * The calibration is the tool's own, read from {@code <db>_<tool>_<calibration key>_summary.csv}.
     * Its {@code rho f} is empty, since an external tool has no refined variant, and the refined
     * estimate therefore stays empty as well.
     *
     * @param db             the name of the database project
     * @param fqMapFile      the fastq mapping file the runs above used
     * @param reportKey      the key the result file is named after
     * @param calibrationKey the report key of the simulated run supplying rho, may be {@code null}
     * @param obsBaseline    the observable genus-only subset of the unrefined run
     * @param unrefined      the tallies of the unrefined run, whose keys name the output files
     * @throws IOException if an output file cannot be read or a CSV cannot be written
     */
    private void writeExternal(String db, String fqMapFile, String reportKey, String calibrationKey,
                               GenusOnlyBaseline obsBaseline, Map<String, AccuracyTally> unrefined)
            throws IOException {
        File krakenDir = new File(resultsDir, "kraken");
        if (!krakenDir.isDirectory()) {
            return;
        }
        for (String tool : RefinementAccuracyReport.EXTERNAL_TOOLS) {
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
                throw new IOException("Only " + outputs.size() + " of " + unrefined.size()
                        + " samples have a " + tool + " output under " + krakenDir
                        + ". Classify them all or none.");
            }
            System.out.println("Evaluating " + tool + " on " + fqMapFile);
            obsBaseline.rewind();
            Map<String, AccuracyTally> tallies = evaluator.evaluateExternalWithoutGroundTruth(
                    db, fqMapFile, Variant.UNREFINED.getLoadDbGoalKey(), outputs, obsBaseline);
            Map<String, double[]> calibration = readCalibration(db + "_" + tool, calibrationKey);
            File file = new File(resultsDir, db + "_" + tool + "_" + reportKey + "_specificity.csv");
            try (PrintStream ps = new PrintStream(new FileOutputStream(file), false,
                    StandardCharsets.UTF_8.name())) {
                ps.println("fastq key;sample;reads;classified unrefined;classified refined"
                        + ";obs genus only;obs genus only share"
                        + ";ungated precision unrefined;ungated precision refined"
                        + ";rho u;rho f;est prec g u;est prec g f;obs genus only reported;");
                for (Map.Entry<String, AccuracyTally> e : tallies.entrySet()) {
                    writeExternalRow(ps, e.getKey(), e.getValue(), unrefined.get(e.getKey()),
                            calibration, tool);
                }
            }
            System.out.println("Wrote " + file);
        }
    }

    /**
     * Writes one row for an external classifier, in the columns of this report's own CSV.
     * <p>
     * Three cells cannot be taken from the classifier's tally, because its input is filtered to the
     * reads it classified. {@code bin/kraken_classify.sh} drops the unclassified lines as the tool
     * produces them, which is what keeps the per-read output of a real run from reaching hundreds of
     * gigabytes, and the tally therefore counts the classified reads and nothing else.
     * <ul>
     * <li>{@code reads} comes from the Genestrip run, which saw every read of the sample.</li>
     * <li>{@code classified} is the tally's total, which under the filter <em>is</em> the number of
     * reads the tool classified.</li>
     * <li>The ungated precision divides the sum of scores by the size of the subset the unrefined
     * Genestrip run collected, not by the part of it the tool reported. A read of the subset the
     * tool left unclassified never reaches the tally, and averaging over what the tally saw would be
     * an average over the reads the tool happened to resolve.</li>
     * </ul>
     * The refined cells stay empty, since an external tool has no refined variant, and so does the
     * share of the subset in the classified reads, which would mix one tool's subset with another's
     * denominator. {@code obs genus only reported} says how many of the subset the tool did report.
     *
     * @param ps          the stream to write to
     * @param fastqKey    the key of the sample
     * @param t           the tally of the external classifier
     * @param ref         the tally of the unrefined Genestrip run on the same sample
     * @param calibration the tool's own rho values, keyed by fastq key
     * @param tool        the tool's name, for the report on the terminal
     */
    private static void writeExternalRow(PrintStream ps, String fastqKey, AccuracyTally t,
                                         AccuracyTally ref, Map<String, double[]> calibration,
                                         String tool) {
        long reported = t.getObsGenusOnlyTotal();
        long subset = ref == null ? reported : ref.getObsGenusOnlyTotal();
        double sum = reported == 0 ? 0 : t.getObsGenusOnlyUngatedPrecision() * reported;
        double pu = subset == 0 ? Double.NaN : sum / subset;
        if (reported < subset) {
            System.out.printf("  %s/%s: reported %,d of the %,d reads of the subset;"
                    + " the rest counts as zero.%n", tool, fastqKey, reported, subset);
        }
        double[] rho = calibration.get(fastqKey);
        if (rho == null && calibration.size() == 1) {
            rho = calibration.values().iterator().next();
        }
        ps.print(fastqKey);
        ps.print(';');
        ps.print(SampleNames.display(fastqKey, false));
        ps.print(';');
        ps.print(ref == null ? "" : String.valueOf(ref.getTotal()));
        ps.print(';');
        ps.print(t.getTotal());
        ps.print(";;");
        ps.print(subset);
        ps.print(";;");
        ps.print(format(pu));
        ps.print(';');
        ps.print(';');
        ps.print(rho == null ? "" : format(rho[0]));
        ps.print(";;");
        ps.print(rho == null ? "" : format(rho[0] * pu));
        ps.print(";;");
        ps.print(reported);
        ps.println(';');
    }

    /**
     * Reads the {@code rho u} and {@code rho f} columns of a simulated run's summary, keyed by
     * its raw fastq key.
     * <p>
     * The calibration cannot be computed here: {@code rho_{d,s}} is the ratio of the gated to the
     * ungated gain, and the gated one needs the ground truth this report does not have. It is
     * therefore taken from the simulated run, which is exactly the transfer the paper describes ---
     * and the reason the estimate is only as good as the resemblance between that simulation and
     * this sample.
     *
     * @param db             the name of the database project
     * @param calibrationKey the report key of the simulated run, or {@code null} for none
     * @return fastq key to {rho_u, rho_f}; empty if none was named or its summary is absent. A map
     * of exactly one entry is applied to every sample regardless of its key --- see the call site.
     */
    private Map<String, double[]> readCalibration(String db, String calibrationKey) {
        Map<String, double[]> byKey = new HashMap<String, double[]>();
        if (calibrationKey == null || calibrationKey.isEmpty()) {
            return byKey;
        }
        File summary = new File(resultsDir, db + "_" + calibrationKey + "_summary.csv");
        if (!summary.exists()) {
            System.out.println("No calibration at " + summary + " - the estimate columns stay empty.");
            return byKey;
        }
        try (BufferedReader in = new BufferedReader(new InputStreamReader(
                new FileInputStream(summary), StandardCharsets.UTF_8))) {
            String header = in.readLine();
            if (header == null) {
                return byKey;
            }
            String[] names = header.split(";", -1);
            int keyAt = indexOf(names, "fastq key");
            int rhoUAt = indexOf(names, "rho u");
            int rhoFAt = indexOf(names, "rho f");
            if (keyAt < 0 || rhoUAt < 0) {
                System.out.println("No 'fastq key'/'rho u' columns in " + summary
                        + " - estimate columns stay empty.");
                return byKey;
            }
            String line;
            while ((line = in.readLine()) != null) {
                String[] cells = line.split(";", -1);
                if (cells.length <= Math.max(keyAt, rhoUAt) || cells[rhoUAt].trim().isEmpty()) {
                    continue;
                }
                // rho f may be absent: the summary of an external classifier has no refined
                // variant to form it from. NaN carries that through to the row writer, which
                // leaves the refined cells empty rather than printing a number.
                String rhoF = rhoFAt >= 0 && cells.length > rhoFAt ? cells[rhoFAt].trim() : "";
                try {
                    byKey.put(cells[keyAt].trim(), new double[]{
                            Double.parseDouble(cells[rhoUAt].trim()),
                            rhoF.isEmpty() ? Double.NaN : Double.parseDouble(rhoF)});
                } catch (NumberFormatException e) {
                    // A row whose rho is undefined carries no calibration; skip it rather than fail.
                }
            }
        } catch (IOException e) {
            System.out.println("Cannot read " + summary + " (" + e.getMessage() + ") - estimate columns stay empty.");
        }
        return byKey;
    }

    private static int indexOf(String[] names, String name) {
        for (int i = 0; i < names.length; i++) {
            if (name.equals(names[i].trim())) {
                return i;
            }
        }
        return -1;
    }


    /**
     * Formats a metric, writing an empty field for an undefined value so that no locale-dependent
     * NaN symbol reaches the paper.
     *
     * @param value the value to format
     * @return the formatted value, or the empty string if it is undefined
     */
    private static String format(double value) {
        return Double.isNaN(value) || Double.isInfinite(value) ? "" : Double.toString(value);
    }

}
