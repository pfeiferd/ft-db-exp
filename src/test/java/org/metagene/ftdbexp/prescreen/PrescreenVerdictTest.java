package org.metagene.ftdbexp.prescreen;

import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

/**
 * Pins how a branching-degree histogram is read into the four shares the pre-screen decides on.
 * <p>
 * The row format is {@code taxid;name;children;d1;...;dC;OTHER}: one column per branching degree up
 * to the child count, then the OTHER bucket. Mistaking the last column for a degree would count
 * k-mers pushed up from outside the database as headroom, which is the one thing no refinement of
 * these children can claim.
 */
public class PrescreenVerdictTest {
    private static String[] row(String... cells) {
        return cells;
    }

    /**
     * The real {@code Streptococcus} genus row, abbreviated: the paper reports 99.52 per cent of its
     * k-mers at degree ten or below, and this is where that number comes from.
     */
    @Test
    public void testTheHeadroomIsSummedOverDegreesTwoToTen() {
        // children=4, degrees 1..4 then OTHER: nothing single, 60 at degree 2, 30 at 3, 10 at 4.
        PrescreenReport.Verdict v = PrescreenReport.verdictOf(
                row("1301", "Streptococcus", "4", "0", "60", "30", "10", "0"));
        assertEquals(100, v.total);
        assertEquals(0.0, v.singleShare, 1e-9);
        assertEquals("degrees 2 and 3 and 4 all lie at or below ten", 100.0, v.headroomShare, 1e-9);
    }

    /** The last column is the OTHER bucket and never a branching degree. */
    @Test
    public void testTheOtherBucketIsNotCountedAsHeadroom() {
        PrescreenReport.Verdict v = PrescreenReport.verdictOf(
                row("1", "n", "2", "0", "40", "60"));
        assertEquals(100, v.total);
        assertEquals("40 at degree two", 40.0, v.headroomShare, 1e-9);
        assertEquals("60 pushed up from outside", 60.0, v.otherShare, 1e-9);
    }

    /** Above ten children the headroom stops at ten, however many degrees the row carries. */
    @Test
    public void testDegreesAboveTenAreNotHeadroom() {
        String[] cells = new String[3 + 12 + 1];
        cells[0] = "1"; cells[1] = "n"; cells[2] = "12";
        for (int d = 1; d <= 12; d++) {
            cells[2 + d] = d <= 10 ? "5" : "25";
        }
        cells[cells.length - 1] = "0";
        PrescreenReport.Verdict v = PrescreenReport.verdictOf(cells);
        assertEquals(100, v.total);
        assertEquals("degrees 2..10, nine of them at five each", 45.0, v.headroomShare, 1e-9);
    }

    /** A k-mer carried by one child alone is not there: the LCA update would have pushed it down. */
    @Test
    public void testTheSingleDegreeIsReportedSeparately() {
        PrescreenReport.Verdict v = PrescreenReport.verdictOf(row("1", "n", "2", "20", "80", "0"));
        assertEquals(20.0, v.singleShare, 1e-9);
        assertEquals(80.0, v.headroomShare, 1e-9);
    }

    /** A node holding nothing is not a verdict either way and is left out. */
    @Test
    public void testANodeWithoutKmersYieldsNoVerdict() {
        assertNull(PrescreenReport.verdictOf(row("1", "n", "2", "0", "0", "0")));
    }

    /** A malformed row is skipped rather than guessed at. */
    @Test
    public void testAMalformedRowYieldsNoVerdict() {
        assertNull(PrescreenReport.verdictOf(row("1", "n")));
        assertNull(PrescreenReport.verdictOf(row("1", "n", "notanumber", "5", "5")));
    }
}
