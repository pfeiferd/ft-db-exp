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
    /**
     * Builds a tree of the shape a refined database has below the species, with one file node per
     * genome, and the typing to go with it.
     *
     * @param leafTypes leaf name to sequence type, an empty value marking a genome with none
     */
    /** A model together with the tree it was built from, so a test can address nodes by tax id. */
    private static final class Built {
        private final NaiveBayesSTModel model;
        private final SmallTaxTree tree;

        Built(NaiveBayesSTModel model, SmallTaxTree tree) {
            this.model = model;
            this.tree = tree;
        }

        /** Position of the data node, i.e. the one all the genomes hang under. */
        int dataNode() {
            return tree.getNodeByTaxId("3").getPosition();
        }

        /** Position of the n-th file node, numbered as {@code build} created them. */
        int leaf(int n) {
            return tree.getNodeByTaxId(String.valueOf(10 + n)).getPosition();
        }

        /** Position of the n-th refined node, i.e. of the n-th group given to {@code build}. */
        int group(int n) {
            return tree.getNodeByTaxId(String.valueOf(100 + n)).getPosition();
        }
    }

    /**
     * Builds a tree of the shape a refined database has below the species: one data node carrying
     * every genome, and under it one refined node per group with that group's file nodes below it.
     * The data node is therefore what the prior is taken over, and a refined node is where a log
     * ratio is anything other than zero -- at the node the prior itself comes from it is one by
     * construction, which is the point of the ratio and not a defect to test around.
     *
     * @param groups one map per refined node, leaf name to sequence type, an empty value or
     *               {@code -} marking a genome the scheme does not cover
     */
    @SafeVarargs
    private final Built build(Map<String, String>... groups) throws IOException {
        File dir = Files.createTempDirectory("nbmodel").toFile();
        dir.deleteOnExit();
        StringBuilder nodes = new StringBuilder();
        StringBuilder names = new StringBuilder();
        // 1 root -> 2 species -> 3 data -> one file node per genome
        nodes.append("1\t|\t1\t|\tno rank\t|\t\t|\n");
        names.append("1\t|\troot\t|\t\t|\tscientific name\t|\n");
        nodes.append("2\t|\t1\t|\tspecies\t|\t\t|\n");
        names.append("2\t|\tspecies\t|\t\t|\tscientific name\t|\n");
        nodes.append("3\t|\t2\t|\t").append(Rank.DATA.getName()).append("\t|\t\t|\n");
        names.append("3\t|\tdata\t|\t\t|\tscientific name\t|\n");
        int id = 10;
        int refinedId = 100;
        StringBuilder mlst = new StringBuilder("leaf;scheme;st;alleles;\n");
        for (Map<String, String> group : groups) {
            nodes.append(refinedId).append("\t|\t3\t|\t").append(Rank.REFINED.getName()).append("\t|\t\t|\n");
            names.append(refinedId).append("\t|\tr").append(refinedId).append("\t|\t\t|\tscientific name\t|\n");
            for (Map.Entry<String, String> e : group.entrySet()) {
                nodes.append(id).append("\t|\t").append(refinedId).append("\t|\t")
                        .append(Rank.FILE.getName()).append("\t|\t\t|\n");
                names.append(id).append("\t|\t").append(e.getKey()).append("\t|\t\t|\tscientific name\t|\n");
                mlst.append(e.getKey()).append(";cdifficile;").append(e.getValue()).append(";;\n");
                id++;
            }
            refinedId++;
        }
        Files.write(new File(dir, TaxTree.NODES_DMP).toPath(), nodes.toString().getBytes(StandardCharsets.UTF_8));
        Files.write(new File(dir, TaxTree.NAMES_DMP).toPath(), names.toString().getBytes(StandardCharsets.UTF_8));
        File mlstFile = new File(dir, "mlst.csv");
        Files.write(mlstFile.toPath(), mlst.toString().getBytes(StandardCharsets.UTF_8));

        TaxTree full = new TaxTree(dir, false);
        for (String taxId : new String[] { "1", "2", "3" }) {
            full.getNodeByTaxId(taxId).markRequired();
        }
        for (int i = 10; i < id; i++) {
            full.getNodeByTaxId(String.valueOf(i)).markRequired();
        }
        for (int i = 100; i < refinedId; i++) {
            full.getNodeByTaxId(String.valueOf(i)).markRequired();
        }
        SmallTaxTree tree = full.toSmallTaxTree();
        return new Built(new NaiveBayesSTModel(tree, new STGroundTruth(mlstFile)), tree);
    }

    private Map<String, String> leaves(String... pairs) {
        Map<String, String> m = new LinkedHashMap<>();
        for (int i = 0; i < pairs.length; i += 2) {
            m.put(pairs[i], pairs[i + 1]);
        }
        return m;
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
        // within one model, so the prior and the smoothing are identical and only the dilution
        // separates the two.
        Built b = build(leaves("a.fna", "1", "b.fna", "1"),
                        leaves("c.fna", "1", "d.fna", "1", "e.fna", "-"),
                        leaves("f.fna", "2", "g.fna", "2"));
        double pure = b.model.getLogRatio(b.group(0), "1");
        double diluted = b.model.getLogRatio(b.group(1), "1");
        assertTrue("a node diluted by an untyped genome must not speak more strongly for ST 1: "
                + diluted + " vs " + pure, diluted < pure);
        // And it speaks for the untyped class where the pure node does not.
        assertTrue(b.model.getLogRatio(b.group(1), NaiveBayesSTModel.STAR)
                > b.model.getLogRatio(b.group(0), NaiveBayesSTModel.STAR));
    }

    @Test
    public void aLeafMissingFromTheTypingCountsAsUntyped() throws IOException {
        // Not the same thing as a row saying "-", but the same consequence: the genome is there and
        // the scheme says nothing about it, so it must not vanish from the denominator.
        Built b = build(leaves("a.fna", "1", "b.fna", ""));
        assertEquals(2, b.model.getCandidates(b.dataNode()));
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
