package org.metagene.ftdbexp.stquality;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import org.junit.Test;

/**
 * Tests the sequence-type level measure of the {@code cdiff} case study: the arithmetic of
 * {@link STCounts} and the reading of the typing by {@link STGroundTruth}.
 * <p>
 * What is worth pinning down is the one place the measure differs structurally from the data-taxon
 * level it is modelled on. {@code |D_n|} counts leaves and aggregates upwards by addition;
 * {@code |S_n|} is the size of a union, because two genomes below a node may be the same lineage, and
 * adding there would report more candidates than exist and make every precision below the node look
 * worse than it is. The worked example below is built so that the two answers differ.
 */
public class STQualityMeasureTest {
    private static final double EPS = 1e-9;

    /**
     * Three genomes below a node, two of them the same sequence type. The set has two members, not
     * three, and that is what every precision at the node divides by.
     */
    @Test
    public void testSequenceTypesUnionRatherThanAdd() {
        STCounts leafA = new STCounts(true, 0);
        STCounts leafB = new STCounts(true, 0);
        STCounts leafC = new STCounts(true, 0);
        leafA.addST(0);
        leafB.addST(0);   // same lineage as A
        leafC.addST(1);

        STCounts node = new STCounts(false, 10);
        node.unionSTs(leafA);
        node.unionSTs(leafB);
        node.unionSTs(leafC);

        assertEquals("three genomes, two lineages", 2, node.getSTCount());
        assertEquals(1, leafA.getSTCount());
    }

    /**
     * The node precision is {@code c_st / (|S_n| * kmers)}, averaged over the node's own k-mers. With
     * two lineages below and both carrying every one of the four k-mers, it is one: the k-mers sit
     * exactly where they belong, because nothing below the node distinguishes them.
     */
    @Test
    public void testNodePrecisionIsOneWhenEveryLineageCarriesEveryKMer() {
        STCounts node = new STCounts(false, 4);
        node.addST(0);
        node.addST(1);
        for (int i = 0; i < 8; i++) {          // 4 k-mers x 2 lineages
            node.incTpForNodePrecision();
        }
        assertEquals(1.0, node.getNodePrecision(), EPS);
    }

    /**
     * The other end: each of the four k-mers is carried by one of the two lineages only, so half of
     * the candidates the node offers are wrong for any given k-mer.
     */
    @Test
    public void testNodePrecisionHalvesWhenOnlyOneOfTwoLineagesCarriesEachKMer() {
        STCounts node = new STCounts(false, 4);
        node.addST(0);
        node.addST(1);
        for (int i = 0; i < 4; i++) {          // 4 k-mers x 1 lineage each
            node.incTpForNodePrecision();
        }
        assertEquals(0.5, node.getNodePrecision(), EPS);
    }

    /**
     * A k-mer sitting at a single genome scores one whatever the database looks like, since there is
     * only one lineage it could belong to. That is the case {@code sp*} exists to exclude.
     */
    @Test
    public void testALeafScoresOneByConstructionAndIsExcludedFromTheRestrictedAverage() {
        STCounts leaf = new STCounts(true, 100);
        leaf.addST(0);
        for (int i = 0; i < 100; i++) {
            leaf.incTpForNodePrecision();
        }
        assertEquals(1.0, leaf.getNodePrecision(), EPS);

        STCounts root = new STCounts(false, 0);
        root.unionSTs(leaf);
        root.aggregateSubtree(leaf);
        assertEquals("the leaf's k-mers count towards sp", 100, root.getSubtreeKmerSum());
        assertEquals(1.0, root.getSubtreePrecision(), EPS);
        assertEquals("but not towards sp*", 0, root.getSubtreeKmersAboveData());
        assertTrue("sp* is undefined with nothing above the data taxa",
                Double.isNaN(root.getRestrictedSubtreePrecision()));
    }

    /**
     * The subtree precision is the k-mer-weighted mean over the subtree, so a node holding many poor
     * k-mers outweighs one holding few good ones - and a node holding none does not enter at all,
     * which is what keeps the average insensitive to how many nodes a refinement inserts.
     */
    @Test
    public void testSubtreePrecisionIsWeightedByKMersAndIgnoresEmptyNodes() {
        STCounts inner = new STCounts(false, 30);   // 30 k-mers at precision 1/3
        inner.addST(0);
        inner.addST(1);
        inner.addST(2);
        for (int i = 0; i < 30; i++) {
            inner.incTpForNodePrecision();
        }
        STCounts other = new STCounts(false, 10);   // 10 k-mers at precision 1
        other.addST(0);
        other.addST(1);
        for (int i = 0; i < 20; i++) {
            other.incTpForNodePrecision();
        }
        STCounts empty = new STCounts(false, 0);    // holds nothing, must not count
        empty.addST(0);

        STCounts root = new STCounts(false, 0);
        root.aggregateSubtree(inner);
        root.aggregateSubtree(other);
        root.aggregateSubtree(empty);

        assertEquals(40, root.getSubtreeKmerSum());
        // (30 * 1/3 + 10 * 1) / 40
        assertEquals(20.0 / 40.0, root.getSubtreePrecision(), EPS);
        assertEquals("nothing here is a leaf, so sp* equals sp",
                root.getSubtreePrecision(), root.getRestrictedSubtreePrecision(), EPS);
    }

    /**
     * The two recalls answer different questions and are both worth having: pooled, the node's k-mers
     * decide; unweighted, each lineage gets one vote. A collection where ten of 196 types hold well
     * over half the genomes can make those disagree sharply.
     */
    @Test
    public void testRecallIsPooledAndAveragedOverTheLineagesBelow() {
        STCounts node = new STCounts(false, 0);
        node.addTypeRecall(10, 10);     // a lineage the database claims completely
        node.addTypeRecall(1, 2);       // and one it claims by half

        assertEquals("pooled by k-mers", 11.0 / 12.0, node.getRecall(), EPS);
        assertEquals("one vote per lineage", (1.0 + 0.5) / 2, node.getAvgRecall(), EPS);
        assertEquals(11, node.getTp());
        assertEquals(12, node.getTpPlusFn());
    }

    /**
     * A lineage that contributed no pair at all - none of its genomes' k-mers is in the database -
     * has no recall, and must not be averaged in as a zero. Its absence is a fact about coverage, not
     * a failure of placement, and counting it would confuse the two.
     */
    @Test
    public void testALineageWithoutPairsDoesNotDragTheAverageDown() {
        STCounts node = new STCounts(false, 0);
        node.addTypeRecall(4, 4);
        node.addTypeRecall(0, 0);

        assertEquals(1.0, node.getRecall(), EPS);
        assertEquals(1.0, node.getAvgRecall(), EPS);
    }

    /** With nothing found below it, a node has no recall rather than a recall of zero. */
    @Test
    public void testNoPairsMeansNoRecall() {
        STCounts node = new STCounts(false, 10);
        assertTrue(Double.isNaN(node.getRecall()));
        assertTrue(Double.isNaN(node.getAvgRecall()));
    }

    /**
     * The set is walkable, which is how each node collects the read counts of exactly the lineages
     * below it without the genomes being traversed a second time.
     */
    @Test
    public void testTheSequenceTypeSetCanBeWalked() {
        STCounts node = new STCounts(false, 0);
        node.addST(2);
        node.addST(7);
        node.addST(2);
        StringBuilder seen = new StringBuilder();
        for (int st = node.nextST(0); st >= 0; st = node.nextST(st + 1)) {
            seen.append(st).append(' ');
        }
        assertEquals("2 7 ", seen.toString());
        assertEquals(-1, node.nextST(8));
    }

    /** A node with no typed genome under it has no measure, rather than a bad one. */
    @Test
    public void testNoSequenceTypeMeansNoMeasure() {
        STCounts node = new STCounts(false, 50);
        assertEquals(0, node.getSTCount());
        assertTrue(Double.isNaN(node.getNodePrecision()));
        assertTrue(Double.isNaN(node.getSubtreePrecision()));
    }

    /** The typing is joined on the leaf's name, and an untyped genome contributes no unit. */
    @Test
    public void testGroundTruthReadsTheLeafKeyedSchema() throws IOException {
        File csv = write("leaf;scheme;st;alleles;",
                "GCA_000001.1_A_genomic.fna.gz;cdifficile;11;adk(1);",
                "GCA_000002.1_B_genomic.fna.gz;cdifficile;11;adk(1);",
                "GCA_000003.1_C_genomic.fna.gz;cdifficile;37;adk(2);",
                "GCA_000004.1_D_genomic.fna.gz;cdifficile;-;adk(?);");

        STGroundTruth gt = new STGroundTruth(csv);
        assertEquals("two lineages over three typed genomes", 2, gt.getSTCount());
        assertEquals(3, gt.getTypedLeaves());
        assertEquals(1, gt.getUntypedLeaves());
        assertEquals("same ST, same index",
                gt.getSTIndex("GCA_000001.1_A_genomic.fna.gz"), gt.getSTIndex("GCA_000002.1_B_genomic.fna.gz"));
        assertFalse(gt.getSTIndex("GCA_000001.1_A_genomic.fna.gz") == gt.getSTIndex("GCA_000003.1_C_genomic.fna.gz"));
        assertEquals("untyped genomes join to nothing", -1, gt.getSTIndex("GCA_000004.1_D_genomic.fna.gz"));
        assertEquals("and so does a name the typing never saw", -1, gt.getSTIndex("nowhere.fna.gz"));
        assertEquals("11", gt.getSTName(gt.getSTIndex("GCA_000001.1_A_genomic.fna.gz")));
    }

    /**
     * The older schema of the typing script is still read, keyed by its `file' column. It is worth
     * accepting: such a file has the wrong keys in it rather than the wrong columns, and it should be
     * rejected for joining to nothing - which the counting goal says plainly - and not for a header
     * this class failed to recognise.
     */
    @Test
    public void testGroundTruthReadsTheOlderFileKeyedSchema() throws IOException {
        File csv = write("file;assembly;scheme;st;alleles;",
                "AP025558.1.fa;AP025558;cdifficile;81;adk(3);",
                "AP031492.1.fa;AP031492;cdifficile;2;adk(1);");

        STGroundTruth gt = new STGroundTruth(csv);
        assertEquals(2, gt.getSTCount());
        assertEquals("81", gt.getSTName(gt.getSTIndex("AP025558.1.fa")));
    }

    /** A file with neither key column is refused, rather than read as an empty typing. */
    @Test
    public void testGroundTruthRefusesAnUnrecognisableHeader() throws IOException {
        File csv = write("genome;type;", "x;1;");
        try {
            new STGroundTruth(csv);
            throw new AssertionError("expected an IllegalStateException");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage().contains("'st' column"));
        }
    }

    private static File write(String... lines) throws IOException {
        File f = Files.createTempFile("stquality", ".csv").toFile();
        f.deleteOnExit();
        try (PrintWriter w = new PrintWriter(f, StandardCharsets.UTF_8)) {
            for (String line : lines) {
                w.println(line);
            }
        }
        return f;
    }
}
