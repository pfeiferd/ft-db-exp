package org.metagene.ftdbexp.stquality;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;

import org.junit.Test;

/**
 * Unit tests for the read classifier's own rule applied to a whole isolate.
 * <p>
 * The two properties worth pinning are the ones that made the rule worth trying here at all: the node
 * every candidate hangs under cannot decide anything however much of the sample rests on it, and a
 * tie between two reference genomes answers with the node grouping them rather than with nothing.
 */
public class PathVoteSTModelTest {
    private final STModelFixture fixture = new STModelFixture();

    private PathVoteSTModel model(STModelFixture.Built b) {
        return new PathVoteSTModel(b.composition());
    }

    private Map<Integer, Long> counts(Object... pairs) {
        Map<Integer, Long> m = new LinkedHashMap<>();
        for (int i = 0; i < pairs.length; i += 2) {
            m.put((Integer) pairs[i], ((Number) pairs[i + 1]).longValue());
        }
        return m;
    }

    @Test
    public void theHeaviestPathWins() throws IOException {
        STModelFixture.Built b = fixture.build(fixture.leaves("a.fna", "1", "b.fna", "1"),
                                               fixture.leaves("c.fna", "2", "d.fna", "2"));
        // Everything rests on the node above all the genomes except for a little at one group, and
        // that little decides: it is the only count that is not on both candidates' paths.
        Map<Integer, Long> c = counts(b.dataNode(), 1_000_000L, b.group(1), 5L);
        assertEquals("2", model(b).classify(c));
    }

    @Test
    public void theNodeAboveEverythingCannotDecide() throws IOException {
        // The same counts with the shared node a thousand times heavier still answer the same way.
        // It is an ancestor of every candidate, so it enters every total identically -- which is why
        // this rule needs no arrangement to keep the uninformative mass of a sample out of the
        // decision, where the weighted classifier had to be built for it.
        STModelFixture.Built b = fixture.build(fixture.leaves("a.fna", "1", "b.fna", "1"),
                                               fixture.leaves("c.fna", "2", "d.fna", "2"));
        PathVoteSTModel m = model(b);
        assertEquals("2", m.classify(counts(b.dataNode(), 1_000L, b.group(1), 5L)));
        assertEquals("2", m.classify(counts(b.dataNode(), 1_000_000_000L, b.group(1), 5L)));
    }

    @Test
    public void aTieBetweenTwoGenomesIsAnsweredWithTheirLineage() throws IOException {
        // Two genomes of one type, hit equally hard. Neither identifies the isolate, but the node
        // that groups them does, and that node is what a refinement adds: in an unrefined database
        // the same tie would resolve to the node above all the genomes and say nothing.
        STModelFixture.Built b = fixture.build(fixture.leaves("a.fna", "11", "b.fna", "11"),
                                               fixture.leaves("c.fna", "2", "d.fna", "2"));
        Map<Integer, Long> c = counts(b.leaf(0), 500L, b.leaf(1), 500L);
        PathVoteSTModel m = model(b);
        assertEquals(b.group(0), m.classifyNode(c, PathVoteSTModel.NO_MINIMUM).getPosition());
        assertEquals("11", m.classify(c));
        assertEquals("the lineage node leaves one type in question", 1,
                m.getCandidates(c));
    }

    @Test
    public void oneGenomeAheadOfTheOtherIsAnsweredWithThatGenome() throws IOException {
        // Not a tie, so the rule goes as deep as the counts take it -- nearest neighbour, which is
        // the baseline the weighted classifier has to beat and all an unrefined database can offer.
        STModelFixture.Built b = fixture.build(fixture.leaves("a.fna", "11", "b.fna", "2"));
        Map<Integer, Long> c = counts(b.leaf(0), 500L, b.leaf(1), 499L);
        assertEquals(b.leaf(0), model(b).classifyNode(c, PathVoteSTModel.NO_MINIMUM).getPosition());
    }

    @Test
    public void aMinimumRaisesTheAnswerTowardsTheRoot() throws IOException {
        // Asking for more evidence than one genome carries moves the answer up to the node that has
        // it, which is how Genestrip's minKMersForClass trades depth for support.
        STModelFixture.Built b = fixture.build(fixture.leaves("a.fna", "11", "b.fna", "11"),
                                               fixture.leaves("c.fna", "2", "d.fna", "2"));
        Map<Integer, Long> c = counts(b.leaf(0), 100L, b.group(0), 50L, b.dataNode(), 10L);
        PathVoteSTModel m = model(b);
        assertEquals(b.leaf(0), m.classifyNode(c, PathVoteSTModel.NO_MINIMUM).getPosition());
        assertEquals(b.group(0), m.classifyNode(c, 120L).getPosition());
        assertEquals(b.dataNode(), m.classifyNode(c, 155L).getPosition());
        // And beyond what the whole sample carries there is no answer rather than a bad one.
        assertNull(m.classifyNode(c, 1000L));
    }

    @Test
    public void anAnswerLeavingSeveralTypesOpenEarnsAFraction() throws IOException {
        // The winning node holds two types, so the isolate is classified to both at once. That is a
        // weaker statement than naming one and is credited as one: a half where the true type is
        // among them, nothing where it is not.
        STModelFixture.Built b = fixture.build(fixture.leaves("a.fna", "11", "b.fna", "2"),
                                               fixture.leaves("c.fna", "3", "d.fna", "3"));
        Map<Integer, Long> c = counts(b.leaf(0), 500L, b.leaf(1), 500L);
        PathVoteSTModel m = model(b);
        assertEquals(b.group(0), m.classifyNode(c, PathVoteSTModel.NO_MINIMUM).getPosition());
        assertEquals(2, m.getCandidates(c));
        assertEquals(0.5, m.getPrecision(c, "11"), 1e-9);
        assertEquals(0.5, m.getPrecision(c, "2"), 1e-9);
        assertEquals(0.0, m.getPrecision(c, "3"), 1e-9);
    }

    @Test
    public void namingOneTypeRightlyEarnsAll() throws IOException {
        STModelFixture.Built b = fixture.build(fixture.leaves("a.fna", "11", "b.fna", "11"),
                                               fixture.leaves("c.fna", "2", "d.fna", "2"));
        Map<Integer, Long> c = counts(b.leaf(0), 500L, b.leaf(1), 500L);
        PathVoteSTModel m = model(b);
        assertEquals(1.0, m.getPrecision(c, "11"), 1e-9);
        assertEquals(0.0, m.getPrecision(c, "2"), 1e-9);
    }

    @Test
    public void anAnswerCoveringEverythingIsWorthAlmostNothing() throws IOException {
        // Placing the isolate at the node above all the genomes names every type there is. It is
        // never wrong and it is nearly worthless, and the score has to say so rather than count it
        // as a hit.
        STModelFixture.Built b = fixture.build(fixture.leaves("a.fna", "11", "b.fna", "2"),
                                               fixture.leaves("c.fna", "3", "d.fna", "5"));
        Map<Integer, Long> c = counts(b.dataNode(), 1000L);
        PathVoteSTModel m = model(b);
        assertEquals(b.dataNode(), m.classifyNode(c, PathVoteSTModel.NO_MINIMUM).getPosition());
        assertEquals(4, m.getCandidates(c));
        assertEquals(0.25, m.getPrecision(c, "11"), 1e-9);
    }

    @Test
    public void theInformationGainIsZeroWhereNothingWasNarrowed() throws IOException {
        // Placing the isolate at the node above every genome tells it nothing it did not know, and
        // that is what the gain has to say -- whatever the node's type count is and however frequent
        // the isolate's own type happens to be in the collection.
        STModelFixture.Built b = fixture.build(
                fixture.leaves("a.fna", "11", "b.fna", "11", "c.fna", "11"),
                fixture.leaves("d.fna", "2", "e.fna", "3"));
        Map<Integer, Long> c = counts(b.dataNode(), 1000L);
        PathVoteSTModel m = model(b);
        assertEquals(0.0, m.getInformationGain(c, "11"), 1e-9);
        assertEquals(0.0, m.getInformationGain(c, "3"), 1e-9);
    }

    /**
     * A collection of 26 genomes over six types, so that no single lineage is most of it. The gain is
     * measured against what the collection already made likely, so a fixture in which the true type
     * is already the majority has almost nothing to give and would say nothing about the measure.
     */
    private STModelFixture.Built wideCollection() throws IOException {
        return fixture.build(
                fixture.leaves("a.fna", "11", "b.fna", "11", "c.fna", "11", "d.fna", "11",
                               "e.fna", "11", "f.fna", "2"),
                fixture.leaves("g.fna", "3", "h.fna", "3", "i.fna", "3", "j.fna", "3", "k.fna", "3"),
                fixture.leaves("l.fna", "5", "m.fna", "5", "n.fna", "5", "o.fna", "5", "p.fna", "5"),
                fixture.leaves("q.fna", "7", "r.fna", "7", "s.fna", "7", "t.fna", "7", "u.fna", "7"),
                fixture.leaves("v.fna", "8", "w.fna", "8", "x.fna", "8", "y.fna", "8", "z.fna", "8"));
    }

    @Test
    public void pinningTheLineageEarnsNearlyAll() throws IOException {
        // A node holding one lineage and nothing else is close to a complete answer but not quite
        // one: the smoothing keeps five genomes from proving a lineage outright, and the gain says so
        // rather than rounding it up. It approaches one as the node's genomes outweigh the smoothing.
        STModelFixture.Built b = wideCollection();
        PathVoteSTModel m = model(b);
        double gain = m.getInformationGain(counts(b.group(1), 10L), "3");
        assertTrue("a pure node must score near a complete answer: " + gain, gain > 0.9);
        assertEquals(1.0, m.getPrecision(counts(b.group(1), 10L), "3"), 1e-9);
    }

    @Test
    public void aLopsidedNodeIsNotReadAsACoinToss() throws IOException {
        // Five genomes of one lineage and one of another. The credit 1/n calls that a coin toss at
        // 0.5; the gain reads what the node actually says and puts it near a complete answer. That
        // difference is the reason the gain exists.
        STModelFixture.Built b = wideCollection();
        Map<Integer, Long> c = counts(b.group(0), 10L);
        PathVoteSTModel m = model(b);
        assertEquals(2, m.getCandidates(c));
        assertEquals(0.5, m.getPrecision(c, "11"), 1e-9);
        double gain = m.getInformationGain(c, "11");
        assertTrue("the lopsided node must score far above a coin toss: " + gain, gain > 0.8);
        // And the one genome of the other lineage is not thereby ruled out, only made a poor answer.
        assertTrue(m.getInformationGain(c, "2") < 0.5);
    }

    @Test
    public void narrowingTowardsTheWrongLineageEarnsNothing() throws IOException {
        STModelFixture.Built b = fixture.build(fixture.leaves("a.fna", "11", "b.fna", "11"),
                                               fixture.leaves("c.fna", "2", "d.fna", "2"));
        Map<Integer, Long> c = counts(b.group(1), 10L);
        PathVoteSTModel m = model(b);
        assertEquals(0.0, m.getInformationGain(c, "11"), 1e-9);
        assertEquals(false, m.isHit(c, "11"));
        assertEquals(true, m.isHit(c, "2"));
    }

    @Test
    public void nothingCountedIsNoCall() throws IOException {
        STModelFixture.Built b = fixture.build(fixture.leaves("a.fna", "1", "b.fna", "2"));
        assertNull(model(b).classify(new LinkedHashMap<Integer, Long>()));
    }
}
