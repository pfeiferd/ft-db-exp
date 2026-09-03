package org.metagene.ftdbexp.taxquality;

import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

/**
 * Pins the bucketing and the medians of the nearest-node summary, the file a table of the paper is
 * typeset from.
 * <p>
 * The coarsest bucket takes the root of the walk whatever its candidate count. Deciding it by the
 * node and not by a number is what keeps it right when the number moves, and it did move: empty
 * OTHER buckets were counted as a species until 2026-09-03, so the genus of {@code strepto} reported
 * 158 where it holds 157.
 */
public class SummaryBucketTest {
    @Test
    public void testTheRootIsAlwaysTheCoarsestBucket() {
        assertEquals(4, TaxonCallReport.bucketOf(157, true));
        assertEquals("even when the count is wrong", 4, TaxonCallReport.bucketOf(158, true));
        assertEquals("and even when it is small", 4, TaxonCallReport.bucketOf(3, true));
    }

    /** A node leaving more than a hundred species open has narrowed nothing worth reporting apart. */
    @Test
    public void testANodeLeavingMoreThanAHundredJoinsTheRoot() {
        assertEquals(3, TaxonCallReport.bucketOf(100, false));
        assertEquals(4, TaxonCallReport.bucketOf(101, false));
        assertEquals(4, TaxonCallReport.bucketOf(120, false));
    }

    @Test
    public void testTheRemainingBoundaries() {
        assertEquals(0, TaxonCallReport.bucketOf(1, false));
        assertEquals(1, TaxonCallReport.bucketOf(2, false));
        assertEquals(1, TaxonCallReport.bucketOf(10, false));
        assertEquals(2, TaxonCallReport.bucketOf(11, false));
        assertEquals(2, TaxonCallReport.bucketOf(50, false));
        assertEquals(3, TaxonCallReport.bucketOf(51, false));
    }

    /** A bucket no sample reached has no median, and says so rather than saying zero. */
    @Test
    public void testAnEmptyBucketHasNoMedian() {
        assertNull(TaxonCallReport.median(new ArrayList<Long>()));
    }

    @Test
    public void testTheMedianOfAnOddAndAnEvenCount() {
        assertEquals(Long.valueOf(5), TaxonCallReport.median(Arrays.asList(9L, 1L, 5L)));
        assertEquals(Long.valueOf(3), TaxonCallReport.median(Arrays.asList(1L, 5L)));
    }

    /** The input is not reordered: the caller's list is used again for other columns. */
    @Test
    public void testTheMedianLeavesItsInputAlone() {
        List<Long> values = new ArrayList<>(Arrays.asList(9L, 1L, 5L));
        TaxonCallReport.median(values);
        assertEquals(Arrays.asList(9L, 1L, 5L), values);
    }
}
