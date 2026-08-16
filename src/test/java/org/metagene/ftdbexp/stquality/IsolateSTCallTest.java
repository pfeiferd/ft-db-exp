package org.metagene.ftdbexp.stquality;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import org.junit.Test;
import org.metagene.ftdbexp.stquality.IsolateSTCall.Verdict;

/**
 * Unit tests for the held-out classification result: what an isolate is called for when every node
 * predicts the majority type of the genomes below it, and what the two rates mean.
 */
public class IsolateSTCallTest {
    /** Records {@code n} reads landing on a node that predicts {@code st} among {@code candidates} types. */
    private static void reads(IsolateSTCall call, int n, String st, int candidates) {
        reads(call, n, st, candidates, st == null ? 0 : st.hashCode());
    }

    /** The same, at an explicit node position, which the naive Bayes classifier scores from. */
    private static void reads(IsolateSTCall call, int n, String st, int candidates, int nodePos) {
        for (int i = 0; i < n; i++) {
            call.recordClassified(nodePos, st, candidates);
        }
    }

    @Test
    public void theTypeMostReadsPredictIsCalled() {
        IsolateSTCall c = new IsolateSTCall("B11", "42");
        reads(c, 80, "42", 1);
        reads(c, 20, "1", 58);
        assertEquals("42", c.getCalledST());
        assertEquals(80, c.getCalledReads());
        assertEquals(0.8, c.getCalledShare(), 1e-9);
        assertEquals(Verdict.CORRECT, c.getMajorityVerdict());
    }

    @Test
    public void noThresholdIsNeededForACall() {
        // A single classified read decides, and that is the point of a per-node predictor: the read
        // did not "nearly" name a type, the node it landed on names one.
        IsolateSTCall c = new IsolateSTCall("B11", "42");
        reads(c, 1, "42", 1);
        assertEquals(Verdict.CORRECT, c.getMajorityVerdict());
    }

    @Test
    public void thePriorCanCarryTheCallAway() {
        // The failure mode worth having a test for. Most reads rest on a broad node whose majority
        // is the collection's most frequent type; only a few reach a node specific to this isolate.
        // Counted one read one vote, the prior wins and the call is wrong -- which is exactly what
        // the read accuracy below reports, and why it is reported.
        IsolateSTCall c = new IsolateSTCall("B11", "42");
        reads(c, 950, "1", 58);
        reads(c, 50, "42", 1);
        assertEquals("1", c.getCalledST());
        assertEquals(Verdict.WRONG, c.getMajorityVerdict());
        assertEquals(0.05, c.getReadAccuracy(), 1e-9);
    }

    @Test
    public void readAccuracyCountsPredictionsAndNotVotes() {
        IsolateSTCall c = new IsolateSTCall("B11", "42");
        reads(c, 30, "42", 2);
        reads(c, 70, "1", 58);
        assertEquals(100, c.getClassified());
        assertEquals(30, c.getCorrectReads());
        assertEquals(0.3, c.getReadAccuracy(), 1e-9);
    }

    @Test
    public void meanCandidatesTracksWhatTheRefinementNarrows() {
        IsolateSTCall c = new IsolateSTCall("B11", "42");
        reads(c, 50, "42", 1);
        reads(c, 50, "1", 59);
        assertEquals(30.0, c.getMeanCandidates(), 1e-9);
    }

    @Test
    public void aNodeWithoutTypedGenomesPredictsNothing() {
        IsolateSTCall c = new IsolateSTCall("B11", "42");
        reads(c, 10, null, 0);
        assertEquals(10, c.getClassified());
        assertEquals(0, c.getCorrectReads());
        assertNull(c.getCalledST());
        assertEquals(Verdict.NO_CALL, c.getMajorityVerdict());
    }

    @Test
    public void unclassifiedReadsCountOnlyTowardsTheTotal() {
        IsolateSTCall c = new IsolateSTCall("B11", "42");
        for (int i = 0; i < 1000; i++) {
            c.recordUnclassified();
        }
        assertEquals(1000, c.getReads());
        assertEquals(0, c.getClassified());
        assertNull(c.getCalledST());
        assertEquals(Verdict.NO_CALL, c.getMajorityVerdict());
    }

    @Test
    public void anUntypedIsolateIsNotScored() {
        IsolateSTCall c = new IsolateSTCall("B11", null);
        reads(c, 80, "1", 1);
        assertEquals(Verdict.UNTYPED, c.getMajorityVerdict());
    }

    @Test
    public void mergingIsExact() {
        IsolateSTCall a = new IsolateSTCall("B11", "42");
        reads(a, 30, "42", 1);
        reads(a, 5, "1", 58);
        IsolateSTCall b = new IsolateSTCall("B11", "42");
        reads(b, 50, "42", 1);
        reads(b, 15, "1", 58);
        a.add(b);
        assertEquals(80, a.getPredictions().get("42").longValue());
        assertEquals(20, a.getPredictions().get("1").longValue());
        assertEquals(80, a.getCorrectReads());
        assertEquals(100, a.getReads());
    }

    @Test
    public void readsPerNodeAreKeptForTheClassifier() {
        IsolateSTCall c = new IsolateSTCall("B11", "42");
        reads(c, 30, "42", 1, 7);
        reads(c, 20, "1", 58, 3);
        reads(c, 10, "42", 1, 7);
        assertEquals(40, c.getReadsPerNode().get(7).longValue());
        assertEquals(20, c.getReadsPerNode().get(3).longValue());
    }

    @Test
    public void aVerdictIsTakenAgainstWhicheverCallIsGiven() {
        IsolateSTCall c = new IsolateSTCall("B11", "42");
        reads(c, 950, "1", 58);
        reads(c, 50, "42", 1);
        // The majority is carried away by the prior; a classifier that weighs the reads differently
        // may still land on the right type, and the verdict follows the call it is given.
        assertEquals(Verdict.WRONG, c.getVerdict(c.getCalledST()));
        assertEquals(Verdict.CORRECT, c.getVerdict("42"));
        assertEquals(Verdict.NO_CALL, c.getVerdict(null));
    }

    @Test
    public void kMersAreTalliedPerNodeAndScoredSeparately() {
        IsolateSTCall c = new IsolateSTCall("B11", "42");
        // The read level sees one broad node; the k-mer level sees that most of the matched k-mers
        // sit at a node specific to the isolate. That is the case this second predictor exists for.
        reads(c, 100, "1", 58, 3);
        c.recordKMers(3, 200, 50, "1");
        c.recordKMers(9, 800, 150, "42");
        assertEquals(1000, c.getMatchedKMers());
        assertEquals(0.8, c.getKMerAccuracy(), 1e-9);
        assertEquals(0.0, c.getReadAccuracy(), 1e-9);
        assertEquals(800, c.getKMersPerNode().get(9).longValue());
        assertEquals(200, c.getKMersPerNode().get(3).longValue());
        // Dieselbe Evidenz einmal je distinktem k-mer gezaehlt faellt hier anders aus, weil die
        // Multiplizitaet ungleich verteilt ist - genau der Unterschied, den beide Spalten zeigen.
        assertEquals(200, c.getMatchedUniqueKMers());
        assertEquals(0.75, c.getUniqueKMerAccuracy(), 1e-9);
    }

    @Test
    public void kMerCountsOfZeroAreIgnored() {
        IsolateSTCall c = new IsolateSTCall("B11", "42");
        c.recordKMers(5, 0, 0, "42");
        assertEquals(0, c.getMatchedKMers());
        assertTrue(c.getKMersPerNode().isEmpty());
    }

    @Test
    public void mergingCoversTheKMerLevelToo() {
        IsolateSTCall a = new IsolateSTCall("B11", "42");
        a.recordKMers(9, 300, 30, "42");
        a.recordKMers(3, 100, 10, "1");
        IsolateSTCall b = new IsolateSTCall("B11", "42");
        b.recordKMers(9, 500, 50, "42");
        a.add(b);
        assertEquals(900, a.getMatchedKMers());
        assertEquals(800, a.getKMersPerNode().get(9).longValue());
        assertEquals(800.0 / 900, a.getKMerAccuracy(), 1e-9);
        assertEquals(90, a.getMatchedUniqueKMers());
        assertEquals(80.0 / 90, a.getUniqueKMerAccuracy(), 1e-9);
    }

    @Test
    public void resetEmptiesTheTallyForReuse() {
        IsolateSTCall c = new IsolateSTCall("B11", "42");
        reads(c, 30, "42", 1);
        c.recordKMers(9, 100, 20, "42");
        c.reset();
        assertTrue(c.getReadsPerNode().isEmpty());
        assertTrue(c.getKMersPerNode().isEmpty());
        assertTrue(c.getUniqueKMersPerNode().isEmpty());
        assertEquals(0, c.getMatchedKMers());
        assertEquals(0, c.getMatchedUniqueKMers());
        assertEquals(0, c.getReads());
        assertEquals(0, c.getClassified());
        assertEquals(0, c.getCorrectReads());
        assertNull(c.getCalledST());
    }

    @Test
    public void truthIsReadByKeyAndRejectsContradictions() throws IOException {
        File dir = Files.createTempDirectory("isolatest").toFile();
        dir.deleteOnExit();
        File ok = new File(dir, "ok.csv");
        Files.write(ok.toPath(), ("# comment\nisolate;st;\nB11;1;\nB11np;1;\nB12;11;\n")
                .getBytes(StandardCharsets.UTF_8));
        IsolateSTTruth truth = new IsolateSTTruth(ok);
        assertEquals(3, truth.size());
        assertEquals("1", truth.getST("B11"));
        assertEquals("1", truth.getST("B11np"));
        assertEquals("11", truth.getST("B12"));
        assertNull(truth.getST("B99"));

        File clash = new File(dir, "clash.csv");
        Files.write(clash.toPath(), ("isolate;st;\nB11;1;\nB11;2;\n").getBytes(StandardCharsets.UTF_8));
        try {
            new IsolateSTTruth(clash);
            org.junit.Assert.fail("a key typed twice with differing types must be refused");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("B11"));
        }
    }
}
