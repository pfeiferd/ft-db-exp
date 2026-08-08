package org.metagene.ftdbexp.eval;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Collections;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Test;
import org.metagene.genestrip.tax.Rank;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * Checks that what the evaluation reports is what the paper defines.
 * <p>
 * The measures are two averages over two different subsets, and the two are easy to cross by
 * accident because they are computed side by side from the same reads:
 * <ul>
 * <li>the restricted precision, an average of the gated count $q$ over the genus-only subset $R_g$,
 * whose membership follows from the read's true species and is therefore unavailable on a real
 * sample;</li>
 * <li>the specificity, an average of the ungated count $q'$ over the observable subset
 * $R_g^\circ$, whose membership follows from the assignment alone.</li>
 * </ul>
 * Each test below states one of those definitions arithmetically and compares it against the
 * accessor the report writes into its CSV. The last one checks the CSV itself: a header with fewer
 * fields than its rows silently shifts every column a reader picks out by name.
 */
public class MeasureConsistencyTest {

    private static final double EPS = 1e-12;

    /**
     * The gated precision is the average of $q$ over $R_g$ -- summed over the reads of that subset
     * and divided by its size, with reads outside it contributing to neither.
     */
    @Test
    public void gatedPrecisionAveragesOverItsOwnSubset() {
        AccuracyTally tally = new AccuracyTally();
        // Three reads in R_g scoring 1/2, 1/4 and 0, and one read outside it scoring 1 -- which
        // must not reach the restricted average even though it is the best-scoring read of all.
        tally.record(true, Rank.GENUS, 0.5, 0.5, true, true);
        tally.record(true, Rank.GENUS, 0.25, 0.25, true, true);
        tally.record(true, Rank.GENUS, 0.0, 0.5, true, false);
        tally.record(true, Rank.SPECIES, 1.0, 1.0, false, false);

        assertEquals(3, tally.getGenusOnlyTotal());
        assertEquals((0.5 + 0.25 + 0.0) / 3, tally.getGenusOnlyPrecision(), EPS);
        assertEquals(4, tally.getTotal());
    }

    /**
     * The ungated precision is the average of $q'$ over $R_g^\circ$. It must take its denominator
     * from that subset and not from $R_g$, and must sum the ungated score and not the gated one.
     */
    @Test
    public void ungatedPrecisionAveragesOverTheObservableSubset() {
        AccuracyTally tally = new AccuracyTally();
        // A read in R_g but not in R_g' -- a sibling species was named, so q = 0 while q' = 1/2.
        tally.record(true, Rank.GENUS, 0.0, 0.5, true, false);
        // A read in both.
        tally.record(true, Rank.GENUS, 0.25, 0.25, true, true);
        // A read in R_g' but not in R_g: stopped at a genus that is not the read's own.
        tally.record(true, Rank.FAMILY, 0.0, 0.5, false, true);

        assertEquals(2, tally.getGenusOnlyTotal());
        assertEquals(2, tally.getObsGenusOnlyTotal());
        // Neither subset contains the other, and the two averages share no read but the middle one.
        assertEquals((0.0 + 0.25) / 2, tally.getGenusOnlyPrecision(), EPS);
        assertEquals((0.25 + 0.5) / 2, tally.getObsGenusOnlyUngatedPrecision(), EPS);
    }

    /**
     * On a sample without ground truth only the observable subset can be formed, so the gated
     * measure must stay undefined rather than quietly reporting an average over no reads.
     */
    @Test
    public void groundTruthFreePathFillsOnlyTheObservableSubset() {
        AccuracyTally tally = new AccuracyTally();
        tally.recordWithoutGroundTruth(true, 0.5, true);
        tally.recordWithoutGroundTruth(true, 0.25, true);
        tally.recordWithoutGroundTruth(true, 1.0, false);

        assertEquals(2, tally.getObsGenusOnlyTotal());
        assertEquals((0.5 + 0.25) / 2, tally.getObsGenusOnlyUngatedPrecision(), EPS);
        assertEquals(0, tally.getGenusOnlyTotal());
        assertTrue("the gated precision has no meaning without sigma(r)",
                Double.isNaN(tally.getGenusOnlyPrecision()));
    }

    /**
     * The evaluator shards its tallies per thread and merges them per fastq file, so both subsets
     * have to survive a merge, and a reset has to clear both -- a field forgotten in either place
     * leaks reads of one file into the next.
     */
    @Test
    public void mergingAndResettingKeepBothSubsets() {
        AccuracyTally a = new AccuracyTally();
        a.record(true, Rank.GENUS, 0.5, 0.5, true, true);
        AccuracyTally b = new AccuracyTally();
        b.record(true, Rank.GENUS, 0.25, 0.75, true, true);
        b.record(true, Rank.FAMILY, 0.0, 0.5, false, true);

        AccuracyTally merged = new AccuracyTally();
        merged.add(a);
        merged.add(b);
        assertEquals(2, merged.getGenusOnlyTotal());
        assertEquals(3, merged.getObsGenusOnlyTotal());
        assertEquals((0.5 + 0.25) / 2, merged.getGenusOnlyPrecision(), EPS);
        assertEquals((0.5 + 0.75 + 0.5) / 3, merged.getObsGenusOnlyUngatedPrecision(), EPS);

        b.reset();
        assertEquals(0, b.getGenusOnlyTotal());
        assertEquals(0, b.getObsGenusOnlyTotal());
        assertTrue(Double.isNaN(b.getObsGenusOnlyUngatedPrecision()));
    }

    /**
     * Zero-scoring reads are those of $R_g$ whose true species is not in question at the node they
     * were placed on. They are counted over $R_g$, and a read outside it must not be counted even
     * when it scores zero.
     */
    @Test
    public void zeroScoringIsCountedOverTheGatedSubset() {
        AccuracyTally tally = new AccuracyTally();
        tally.record(true, Rank.GENUS, 0.0, 0.5, true, false);
        tally.record(true, Rank.GENUS, 0.5, 0.5, true, true);
        tally.record(true, Rank.FAMILY, 0.0, 0.5, false, true);

        assertEquals(1, tally.getGenusOnlyZeroScoring());
    }

    /**
     * A header with fewer fields than its rows shifts every column read out by name. The two are
     * written by separate methods, so nothing but a comparison keeps them in step.
     */
    @Test
    public void csvHeaderAndRowHaveTheSameFieldCount() throws Exception {
        Method writeHeader = RefinementAccuracyReport.class
                .getDeclaredMethod("writeHeader", PrintStream.class);
        Method writeRow = RefinementAccuracyReport.class
                .getDeclaredMethod("writeRow", PrintStream.class, String.class,
                        RefinementAccuracyReport.Variant.class, AccuracyTally.class);
        writeHeader.setAccessible(true);
        writeRow.setAccessible(true);

        AccuracyTally tally = new AccuracyTally();
        tally.record(true, Rank.GENUS, 0.5, 0.5, true, true);

        assertEquals("header and row must describe the same number of columns",
                fieldCount(writeHeader, new Object[]{null}),
                fieldCount(writeRow, new Object[]{null, "key",
                        RefinementAccuracyReport.Variant.UNREFINED, tally}));
    }

    /**
     * The summary CSV is the file the paper reads, so its columns have to carry what their names
     * claim: the gated pair from $R_g$, the ungated pair from $R_g^\circ$, the two gains as the
     * differences of those pairs, and $\rho_d$ as their ratio. Crossing any of them would leave a
     * plausible-looking table stating the wrong thing.
     */
    @Test
    public void summaryColumnsCarryTheMeasuresTheyAreNamedAfter() throws Exception {
        // Unrefined: R_g = 2 reads averaging (0.10 + 0.20)/2 = 0.15, R_g' = 4 reads averaging
        // (0.40 + 0.40 + 0.40 + 0.40)/4 = 0.40. The two subsets deliberately differ in size, so a
        // denominator taken from the wrong one cannot pass unnoticed.
        AccuracyTally u = new AccuracyTally();
        u.record(true, Rank.GENUS, 0.10, 0.40, true, true);
        u.record(true, Rank.GENUS, 0.20, 0.40, true, true);
        u.record(true, Rank.FAMILY, 0.0, 0.40, false, true);
        u.record(true, Rank.FAMILY, 0.0, 0.40, false, true);
        // Refined: same subsets (they are fixed by the unrefined run), gated mean 0.35, ungated 0.60.
        AccuracyTally f = new AccuracyTally();
        f.record(true, Rank.GENUS, 0.30, 0.60, true, true);
        f.record(true, Rank.GENUS, 0.40, 0.60, true, true);
        f.record(true, Rank.FAMILY, 0.0, 0.60, false, true);
        f.record(true, Rank.FAMILY, 0.0, 0.60, false, true);

        Map<String, String> row = summaryRow(u, f);
        assertEquals(4L, Long.parseLong(row.get("reads")));
        assertEquals(2L, Long.parseLong(row.get("genus only")));
        assertEquals(4L, Long.parseLong(row.get("obs genus only")));
        assertEquals(0.15, Double.parseDouble(row.get("prec g u")), 1e-6);
        assertEquals(0.35, Double.parseDouble(row.get("prec g f")), 1e-6);
        assertEquals(0.40, Double.parseDouble(row.get("prec g ungated u")), 1e-6);
        assertEquals(0.60, Double.parseDouble(row.get("prec g ungated f")), 1e-6);
        // The calibration is a pair of factors on the levels, not one on the gain: rho_u turns the
        // unrefined ungated precision into the gated one (0.15 / 0.40) and rho_f does the same after
        // the refinement (0.35 / 0.60). Neither gain gets a column -- each is a difference of two
        // columns already present.
        assertEquals(0.15 / 0.40, Double.parseDouble(row.get("rho u")), 1e-6);
        assertEquals(0.35 / 0.60, Double.parseDouble(row.get("rho f")), 1e-6);
        assertFalse("the gain columns were dropped", row.containsKey("delta"));
        assertFalse("the gain columns were dropped", row.containsKey("delta ungated"));
        assertFalse("a single rho was replaced by the pair", row.containsKey("rho"));
    }

    /**
     * The label columns are read straight into a table cell by {@code \csvreader}, in text mode, so
     * a fastq key carrying an underscore -- {@code iss_saliva} did -- stops the paper's build with
     * <em>Missing $ inserted</em>. A key with a name of its own is printed under that name; one
     * without keeps its identity, escaped. The raw key column is neither, since the reports join on
     * it.
     */
    @Test
    public void labelColumnsAreSafeToTypesetAndRawKeysAreNot() throws Exception {
        AccuracyTally u = new AccuracyTally();
        u.record(true, Rank.GENUS, 0.10, 0.40, true, true);
        AccuracyTally f = new AccuracyTally();
        f.record(true, Rank.GENUS, 0.30, 0.60, true, true);

        Map<String, String> named = summaryRow(u, f, "iss_saliva");
        assertEquals("saliva-like", named.get("model"));
        assertEquals("the raw key travels unescaped, since it is what a report joins on",
                "iss_saliva", named.get("fastq key"));

        Map<String, String> fallback = summaryRow(u, f, "iss_novaseq");
        assertEquals("iss\\_novaseq", fallback.get("model"));
        assertEquals("iss_novaseq", fallback.get("fastq key"));
    }

    /**
     * The same escaping guards {@link SpecificityReport}'s sample column, which takes its labels
     * from {@link SampleNames} directly rather than through a table of known keys.
     */
    @Test
    public void sampleLabelsEscapeWhatTexWouldReadAsMarkup() {
        assertEquals("Sim. Tick 3", SampleNames.display("tick3", true));
        assertEquals("Tick 3", SampleNames.display("tick3", false));
        assertEquals("SRR5571991", SampleNames.display("SRR5571991", false));
        assertEquals("a\\_b\\%c\\&d\\#e\\$f", SampleNames.display("a_b%c&d#e$f", false));
        // A semicolon would end the cell rather than mis-typeset it, so it goes too.
        assertFalse(SampleNames.display("a;b", false).contains(";"));
    }

    /** Runs writeSummary into a temporary directory and zips its header against its single row. */
    private Map<String, String> summaryRow(AccuracyTally u, AccuracyTally f) throws Exception {
        return summaryRow(u, f, "set");
    }

    private Map<String, String> summaryRow(AccuracyTally u, AccuracyTally f, String fastqKey)
            throws Exception {
        Map<RefinementAccuracyReport.Variant, Map<String, AccuracyTally>> byVariant =
                new EnumMap<>(RefinementAccuracyReport.Variant.class);
        byVariant.put(RefinementAccuracyReport.Variant.UNREFINED,
                Collections.singletonMap(fastqKey, u));
        byVariant.put(RefinementAccuracyReport.Variant.REFINED,
                Collections.singletonMap(fastqKey, f));

        File dir = Files.createTempDirectory("summary").toFile();
        dir.deleteOnExit();
        Method writeSummary = RefinementAccuracyReport.class.getDeclaredMethod(
                "writeSummary", File.class, String.class, String.class, Map.class);
        writeSummary.setAccessible(true);
        writeSummary.invoke(null, dir, "db", "key", byVariant);

        List<String> lines = Files.readAllLines(new File(dir, "db_key_summary.csv").toPath(),
                StandardCharsets.UTF_8);
        assertEquals("one header and one data row expected", 2, lines.size());
        String[] header = lines.get(0).split(";", -1);
        String[] values = lines.get(1).split(";", -1);
        assertEquals("summary header and row must describe the same number of columns",
                header.length, values.length);
        Map<String, String> row = new HashMap<>();
        for (int i = 0; i < header.length; i++) {
            row.put(header[i], values[i]);
        }
        return row;
    }

    private int fieldCount(Method method, Object[] args) throws Exception {
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        try (PrintStream ps = new PrintStream(out, true, StandardCharsets.UTF_8.name())) {
            args[0] = ps;
            method.invoke(null, args);
        }
        String line = new String(out.toByteArray(), StandardCharsets.UTF_8).trim();
        return line.split(";", -1).length;
    }
}
