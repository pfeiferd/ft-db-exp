package org.metagene.ftdbexp.eval;

import org.metagene.ftdbexp.eval.RefinementAccuracyReport.Variant;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
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
 * <em>specificity gain</em> {@code delta'} of the paper's Section "Estimating the gain without
 * ground truth", which bounds the precision gain from above and, multiplied by the calibration
 * factor {@code rho_d} determined on simulated reads of the same database, estimates it.
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
     * @param reportKey short name used in the result file name
     * @return the file that was written
     * @throws IOException if the databases cannot be read or the file cannot be written
     */
    public File write(String db, String fqMapFile, String reportKey) throws IOException {
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
        File file = new File(resultsDir, db + "_" + reportKey + "_specificity.csv");
        try (PrintStream ps = new PrintStream(new FileOutputStream(file), false, StandardCharsets.UTF_8.name())) {
            ps.println("fastq key;reads;classified unrefined;classified refined"
                    + ";obs genus only;obs genus only share"
                    + ";ungated precision unrefined;ungated precision refined;delta;");
            for (String fastqKey : byVariant.get(Variant.UNREFINED).keySet()) {
                AccuracyTally u = byVariant.get(Variant.UNREFINED).get(fastqKey);
                AccuracyTally f = byVariant.get(Variant.REFINED).get(fastqKey);
                if (u == null || f == null) {
                    continue;
                }
                double pu = u.getObsGenusOnlyUngatedPrecision();
                double pf = f.getObsGenusOnlyUngatedPrecision();
                ps.print(fastqKey);
                ps.print(';');
                ps.print(u.getTotal());
                ps.print(';');
                ps.print(u.getClassified());
                ps.print(';');
                ps.print(f.getClassified());
                ps.print(';');
                ps.print(u.getObsGenusOnlyTotal());
                ps.print(';');
                ps.print(format(u.getTotal() == 0 ? Double.NaN
                        : 100.0 * u.getObsGenusOnlyTotal() / u.getTotal()));
                ps.print(';');
                ps.print(format(pu));
                ps.print(';');
                ps.print(format(pf));
                ps.print(';');
                ps.print(format(pf - pu));
                ps.println(';');
            }
        }
        System.out.println("Wrote " + file);
        return file;
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
