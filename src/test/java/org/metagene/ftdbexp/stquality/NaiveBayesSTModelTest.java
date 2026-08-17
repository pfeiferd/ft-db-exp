package org.metagene.ftdbexp.stquality;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.LinkedHashMap;
import java.util.Map;

import org.junit.Test;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;

/**
 * Unit tests for the per-node predictor, and in particular for what a genome the typing scheme does
 * not cover does to it. Such a genome is evidence that the node is not what its typed genomes alone
 * suggest, so leaving it out would not remove evidence but silently strengthen what remains.
 */
public class NaiveBayesSTModelTest {
    private final STModelFixture fixture = new STModelFixture();

    /** A model together with the fixture's tree, so a test can address nodes by their place in it. */
    private static final class Built {
        private final NaiveBayesSTModel model;
        private final STModelFixture.Built tree;

        Built(STModelFixture.Built tree) {
            this.tree = tree;
            this.model = new NaiveBayesSTModel(tree.composition());
        }

        int dataNode() {
            return tree.dataNode();
        }

        int leaf(int n) {
            return tree.leaf(n);
        }

        int group(int n) {
            return tree.group(n);
        }
    }

    @SafeVarargs
    private final Built build(Map<String, String>... groups) throws IOException {
        return new Built(fixture.build(groups));
    }

    private Map<String, String> leaves(String... pairs) {
        return fixture.leaves(pairs);
    }

    @Test
    public void anUntypedGenomeIsAClassOfItsOwn() throws IOException {
        Built b = build(leaves("a.fna", "1", "b.fna", "1", "c.fna", "-"));
        // Two classes are in question below the data node, not one: ST 1 and the untyped genome.
        // Counting only the typed ones would report a single candidate and a certainty the
        // collection does not have.
        assertEquals(2, b.model.getCandidates(b.dataNode()));
        assertEquals("1", b.model.getMajorityST(b.dataNode()));
        assertTrue(b.model.getTypes().contains(NaiveBayesSTModel.STAR));
    }

    @Test
    public void untypedCanBeTheMajority() throws IOException {
        Built b = build(leaves("a.fna", "1", "b.fna", "-", "c.fna", "-"));
        // A node whose genomes are mostly outside the scheme predicts exactly that, which for this
        // organism is the informative answer rather than a missing one.
        assertEquals(NaiveBayesSTModel.STAR, b.model.getMajorityST(b.dataNode()));
    }

    @Test
    public void anUntypedGenomeDilutesTheEstimateOfTheTypedOnes() throws IOException {
        // Both nodes hold two ST 1 genomes; the second holds an untyped one beside them. Compared
        // within one model, so the baseline and the smoothing are identical and only the dilution
        // separates the two.
        Built b = build(leaves("a.fna", "1", "b.fna", "1"),
                        leaves("c.fna", "1", "d.fna", "1", "e.fna", "-"),
                        leaves("f.fna", "2", "g.fna", "2"));
        double pure = b.model.getLogLift(b.group(0), "1");
        double diluted = b.model.getLogLift(b.group(1), "1");
        assertTrue("a node diluted by an untyped genome must not speak more strongly for ST 1: "
                + diluted + " vs " + pure, diluted < pure);
        // And it speaks for the untyped class where the pure node does not.
        assertTrue(b.model.getLogLift(b.group(1), NaiveBayesSTModel.STAR)
                > b.model.getLogLift(b.group(0), NaiveBayesSTModel.STAR));
    }

    @Test
    public void aNodeCoveringEveryGenomeIsExactlyNeutral() throws IOException {
        // The property the whole classifier rests on, and the one an earlier version lost. Nearly
        // every read of a sample lands at the node holding the shared k-mers of the species, which
        // covers every genome of every class. If its contribution were merely small rather than zero,
        // a few million reads would multiply it into a verdict of their own: at a tenth of a nat per
        // read, half a million nats for whichever class the residue happened to favour, and every
        // isolate would be called that. Exactly zero is therefore not a nicety.
        Built b = build(leaves("a.fna", "1", "b.fna", "1", "c.fna", "1"),
                        leaves("d.fna", "2", "e.fna", "-"));
        for (String type : b.model.getTypes()) {
            assertEquals("the node above every genome must not move class " + type,
                    0.0, b.model.getLogLift(b.dataNode(), type), 0.0);
        }
        // And it stays neutral however many reads rest there: a majority over them would answer with
        // the collection's most frequent type, the classifier answers with nothing at all.
        Map<Integer, Long> reads = new LinkedHashMap<>();
        reads.put(b.dataNode(), 5_000_000L);
        reads.put(b.group(1), 1L);
        assertEquals("2", b.model.classify(reads));
    }

    @Test
    public void aClassWithNoGenomeIsNotACandidate() throws IOException {
        // Every genome typed, so nothing is untyped and the untyped class has no distribution to
        // estimate. Left in and smoothed, it is the one class no node can speak against and it takes
        // every isolate; left out, the question does not arise. The same holds for the sequence types
        // the typing file knows of but this database holds no genome of.
        Built b = build(leaves("a.fna", "1", "b.fna", "1", "c.fna", "1"), leaves("d.fna", "2"));
        assertTrue("an empty class must not be a candidate",
                !b.model.getTypes().contains(NaiveBayesSTModel.STAR));
        assertEquals(2, b.model.getTypes().size());
        Map<Integer, Long> reads = new LinkedHashMap<>();
        reads.put(b.group(1), 1L);
        assertEquals("2", b.model.classify(reads));
    }

    @Test
    public void aLeafMissingFromTheTypingCountsAsUntyped() throws IOException {
        // Not the same thing as a row saying "-", but the same consequence: the genome is there and
        // the scheme says nothing about it, so it must not vanish from the denominator.
        Built b = build(leaves("a.fna", "1", "b.fna", ""));
        assertEquals(2, b.model.getCandidates(b.dataNode()));
    }

    @Test
    public void theCollectionsCompositionDoesNotDecideOnItsOwn() throws IOException {
        // Nine ST 1 genomes against one ST 2, and a single read on the node holding only the ST 2
        // one. Weighted by how many genomes the collection happens to hold, ST 1 would win on the
        // strength of being frequent in a reference set; the classifier must not do that, because
        // the population the isolates come from is not the population the database records.
        Built b = build(leaves("a.fna", "1", "b.fna", "1", "c.fna", "1", "d.fna", "1", "e.fna", "1",
                               "f.fna", "1", "g.fna", "1", "h.fna", "1", "i.fna", "1"),
                        leaves("j.fna", "2"));
        Map<Integer, Long> reads = new LinkedHashMap<>();
        reads.put(b.group(1), 1L);
        assertEquals("2", b.model.classify(reads));
    }

    @Test
    public void classificationCanSettleOnUntyped() throws IOException {
        Built b = build(leaves("a.fna", "1", "b.fna", "1", "c.fna", "1",
                "d.fna", "-", "e.fna", "-"));
        Map<Integer, Long> reads = new LinkedHashMap<>();
        reads.put(b.leaf(3), 100L);   // d.fna, the fourth leaf created
        assertEquals(NaiveBayesSTModel.STAR, b.model.classify(reads));
    }

}
