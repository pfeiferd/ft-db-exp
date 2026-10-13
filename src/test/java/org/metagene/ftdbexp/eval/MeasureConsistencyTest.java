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
 * The measures are two averages over one subset, and the two are easy to cross by accident because
 * they are computed side by side from the same reads:
 * <ul>
 * <li>the restricted precision, an average of the gated count $q$ over the genus-only subset
 * $R_g$;</li>
 * <li>the specificity, an average of the ungated count $q'$ over that same subset.</li>
 * </ul>
 * $R_g$ is read off the assignment alone -- the reads the unrefined database left at or below a
 * genus and above any species -- so it needs no ground truth and the two averages share their
 * denominator. Only $q$ consults the read's true species, which is why the gated average is the one
 * a real sample cannot supply.
 * <p>
 * Each test below states one of those definitions arithmetically and compares it against the
 * accessor the report writes into its CSV. The last ones check the CSV itself: a header with fewer
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
        // Two reads in R_g scoring 1/2 and 1/4, and one read outside it scoring 1 -- which must not
        // reach the restricted average even though it is the best-scoring read of all.
        tally.record(true, Rank.GENUS, 0.5, 0.5, true);
        tally.record(true, Rank.GENUS, 0.25, 0.25, true);
        tally.record(true, Rank.SPECIES, 1.0, 1.0, false);

        assertEquals(2, tally.getGenusOnlyTotal());
        assertEquals((0.5 + 0.25) / 2, tally.getGenusOnlyPrecision(), EPS);
        assertEquals(3, tally.getTotal());
    }

    /**
     * The ungated precision is the average of $q'$ over the very same subset. It must sum the
     * ungated score and divide by the subset's size, so that the two precisions differ in their
     * numerators and in nothing else -- which is what makes their ratio $\rho_d$ a statement about
     * the gate alone.
     */
    @Test
    public void ungatedPrecisionSharesTheDenominatorOfTheGatedOne() {
        AccuracyTally tally = new AccuracyTally();
        // A read of the subset whose assignment is an ancestor of its species: the gate lets it
        // through and both scores count it.
        tally.record(true, Rank.GENUS, 0.25, 0.25, true);
        // A read of the subset resting off its own path: the ungated score counts it, the gated one
        // scores it zero. It stays in the denominator of both.
        tally.record(true, Rank.FAMILY, 0.0, 0.5, true);
        // A read outside the subset: already at its species.
        tally.record(true, Rank.SPECIES, 1.0, 1.0, false);

        assertEquals(2, tally.getGenusOnlyTotal());
        assertEquals((0.25 + 0.0) / 2, tally.getGenusOnlyPrecision(), EPS);
        assertEquals((0.25 + 0.5) / 2, tally.getGenusOnlyUngatedPrecision(), EPS);
    }

    /**
     * On a sample without ground truth the subset is formed exactly as it is elsewhere, since
     * nothing in it needs $\sigma(r)$. Only $q$ is out of reach there, so the gated measures must
     * stay undefined rather than reading as zero, i.e. as "every read wrong".
     */
    @Test
    public void groundTruthFreePathFillsTheSubsetButNotTheGatedScore() {
        AccuracyTally tally = new AccuracyTally();
        tally.recordWithoutGroundTruth(true, 0.5, true);
        tally.recordWithoutGroundTruth(true, 0.25, true);
        tally.recordWithoutGroundTruth(true, 1.0, false);

        assertEquals(2, tally.getGenusOnlyTotal());
        assertEquals((0.5 + 0.25) / 2, tally.getGenusOnlyUngatedPrecision(), EPS);
        assertTrue("the gated precision has no meaning without sigma(r)",
                Double.isNaN(tally.getGenusOnlyPrecision()));
        assertTrue("nor has the share of the subset pinned down to a species",
                Double.isNaN(tally.getGenusOnlySpeciesShare()));
    }

    /**
     * The evaluator shards its tallies per thread and merges them per fastq file, so both sums have
     * to survive a merge, and a reset has to clear both -- a field forgotten in either place leaks
     * reads of one file into the next.
     */
    @Test
    public void mergingAndResettingKeepBothSums() {
        AccuracyTally a = new AccuracyTally();
        a.record(true, Rank.GENUS, 0.5, 0.5, true);
        AccuracyTally b = new AccuracyTally();
        b.record(true, Rank.GENUS, 0.25, 0.75, true);
        b.record(true, Rank.FAMILY, 0.0, 0.5, true);

        AccuracyTally merged = new AccuracyTally();
        merged.add(a);
        merged.add(b);
        assertEquals(3, merged.getGenusOnlyTotal());
        assertEquals((0.5 + 0.25 + 0.0) / 3, merged.getGenusOnlyPrecision(), EPS);
        assertEquals((0.5 + 0.75 + 0.5) / 3, merged.getGenusOnlyUngatedPrecision(), EPS);

        b.reset();
        assertEquals(0, b.getGenusOnlyTotal());
        assertTrue(Double.isNaN(b.getGenusOnlyPrecision()));
        assertTrue(Double.isNaN(b.getGenusOnlyUngatedPrecision()));
    }

    /**
     * Zero-scoring reads are those of $R_g$ whose true species is not in question at the node they
     * were placed on. They are counted over $R_g$, and a read outside it must not be counted even
     * when it scores zero.
     */
    @Test
    public void zeroScoringIsCountedOverTheGatedSubset() {
        AccuracyTally tally = new AccuracyTally();
        tally.record(true, Rank.GENUS, 0.0, 0.5, true);
        tally.record(true, Rank.GENUS, 0.5, 0.5, true);
        tally.record(true, Rank.FAMILY, 0.0, 0.5, false);

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
        tally.record(true, Rank.GENUS, 0.5, 0.5, true);

        assertEquals("header and row must describe the same number of columns",
                fieldCount(writeHeader, new Object[]{null}),
                fieldCount(writeRow, new Object[]{null, "key",
                        RefinementAccuracyReport.Variant.UNREFINED, tally}));
    }

    /**
     * The summary CSV is the file the paper reads, so its columns have to carry what their names
     * claim: the gated pair and the ungated pair over $R_g$, and $\rho_d$ as the ratio of the two.
     * Crossing any of them would leave a plausible-looking table stating the wrong thing.
     */
    @Test
    public void summaryColumnsCarryTheMeasuresTheyAreNamedAfter() throws Exception {
        // Unrefined: R_g = 3 of the 4 reads, gated mean (0.10 + 0.20 + 0.0)/3 = 0.10, ungated mean
        // (0.40 + 0.40 + 0.40)/3 = 0.40. The subset is deliberately smaller than the read count, so
        // a denominator taken from the wrong place cannot pass unnoticed.
        AccuracyTally u = new AccuracyTally();
        u.record(true, Rank.GENUS, 0.10, 0.40, true);
        u.record(true, Rank.GENUS, 0.20, 0.40, true);
        u.record(true, Rank.FAMILY, 0.0, 0.40, true);
        u.record(true, Rank.SPECIES, 1.0, 1.0, false);
        // Refined: same subset (it is fixed by the unrefined run), gated mean 0.30, ungated 0.60.
        AccuracyTally f = new AccuracyTally();
        f.record(true, Rank.GENUS, 0.30, 0.60, true);
        f.record(true, Rank.GENUS, 0.60, 0.60, true);
        f.record(true, Rank.FAMILY, 0.0, 0.60, true);
        f.record(true, Rank.SPECIES, 1.0, 1.0, false);

        Map<String, String> row = summaryRow(u, f);
        assertEquals(4L, Long.parseLong(row.get("reads")));
        assertEquals(3L, Long.parseLong(row.get("genus only")));
        assertFalse("the duplicate subset column was dropped", row.containsKey("obs genus only"));
        assertEquals(0.10, Double.parseDouble(row.get("prec g u")), 1e-6);
        assertEquals(0.30, Double.parseDouble(row.get("prec g f")), 1e-6);
        assertEquals(0.40, Double.parseDouble(row.get("prec g ungated u")), 1e-6);
        assertEquals(0.60, Double.parseDouble(row.get("prec g ungated f")), 1e-6);
        // The calibration is a pair of factors on the levels, not one on the gain: rho_u turns the
        // unrefined ungated precision into the gated one (0.10 / 0.40) and rho_f does the same after
        // the refinement (0.30 / 0.60). Neither gain gets a column -- each is a difference of two
        // columns already present.
        assertEquals(0.10 / 0.40, Double.parseDouble(row.get("rho u")), 1e-6);
        assertEquals(0.30 / 0.60, Double.parseDouble(row.get("rho f")), 1e-6);
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
        u.record(true, Rank.GENUS, 0.10, 0.40, true);
        AccuracyTally f = new AccuracyTally();
        f.record(true, Rank.GENUS, 0.30, 0.60, true);

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
