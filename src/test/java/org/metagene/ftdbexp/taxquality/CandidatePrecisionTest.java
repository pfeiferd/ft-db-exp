package org.metagene.ftdbexp.taxquality;

import org.junit.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

/**
 * Pins the candidate precision of {@link TaxonCallReport} to the definition the paper gives:
 * {@code q = 1/|Sigma(k)|} where the reference organism lies in {@code Sigma(k)}, and {@code 0}
 * where it does not.
 * <p>
 * The measure exists because the all-or-nothing rule that stood here before could not see what a
 * refinement does. That rule named a taxon only when the winning node had a single candidate and
 * reported nothing otherwise, so a call narrowed from forty species to three scored the same as one
 * left at the genus -- exactly the blindness the paper holds against the boolean counts. The cases
 * below are therefore mostly about the middle ground, which is the whole point of the measure.
 * <p>
 * The three sets are plain sets rather than a taxonomy: {@code candidatePrecision} takes the species
 * of the winning node and the reference organisms, and a {@link org.metagene.genestrip.tax.SmallTaxTree}
 * can only be built from the NCBI dumps, which a unit test has no business reading.
 */
public class CandidatePrecisionTest {
    private static Set<String> set(String... taxIds) {
        return new LinkedHashSet<>(Arrays.asList(taxIds));
    }

    /** A call that names the species outright: one candidate, so the full score. */
    @Test
    public void testASpeciesCallScoresOne() {
        assertEquals(1.0, TaxonCallReport.candidatePrecision(set("1313"), set("1313")), 1e-12);
    }

    /** The middle ground the boolean counts collapse: correct, but three species still in question. */
    @Test
    public void testAPartialNarrowingScoresItsShare() {
        assertEquals(1.0 / 3, TaxonCallReport.candidatePrecision(
                set("1313", "1303", "1318"), set("1313")), 1e-12);
    }

    /** The same reference organism, the same call, but a node covering forty species. */
    @Test
    public void testAWiderNodeScoresLess() {
        Set<String> forty = new LinkedHashSet<>();
        for (int i = 0; i < 40; i++) {
            forty.add(Integer.toString(1300 + i));
        }
        assertEquals(1.0 / 40, TaxonCallReport.candidatePrecision(forty, set("1313")), 1e-12);
    }

    /** A call that does not cover the organism at all is simply wrong, however narrow it is. */
    @Test
    public void testAMissScoresZero() {
        assertEquals(0.0, TaxonCallReport.candidatePrecision(set("1314"), set("1313")), 1e-12);
    }

    /** No call: nothing was narrowed down. A score of zero, and not an absent value. */
    @Test
    public void testNoCallScoresZeroRatherThanNothing() {
        assertEquals(0.0, TaxonCallReport.candidatePrecision(null, set("1313")), 1e-12);
        assertEquals(0.0, TaxonCallReport.candidatePrecision(Collections.<String>emptySet(), set("1313")), 1e-12);
    }

    /**
     * A sample whose reference names no species of this taxonomy -- NRF, or an organism from another
     * genus. The measure is undefined there and must stay empty: scoring it zero would count every
     * such sample as a failure of the database.
     */
    @Test
    public void testAnAbsentReferenceIsUndefinedAndNotZero() {
        assertNull(TaxonCallReport.candidatePrecision(set("1313"), Collections.<String>emptySet()));
        assertNull(TaxonCallReport.candidatePrecision(set("1313"), null));
    }

    /** Several organisms in one sample: covering any one of them is a correct narrowing. */
    @Test
    public void testAnyOfSeveralReferenceOrganismsCounts() {
        assertEquals(1.0 / 2, TaxonCallReport.candidatePrecision(
                set("1313", "1314"), set("727", "1314")), 1e-12);
        assertEquals(0.0, TaxonCallReport.candidatePrecision(
                set("1313", "1314"), set("727", "1280")), 1e-12);
    }
}
