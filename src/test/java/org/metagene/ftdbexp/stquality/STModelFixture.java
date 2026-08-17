package org.metagene.ftdbexp.stquality;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.LinkedHashMap;
import java.util.Map;

import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;

/**
 * Builds the little database the classifier tests weigh their rules against, so that the two of them
 * are compared on the same shape rather than on two hand-built approximations of it.
 * <p>
 * The shape is the one a refined database has below the species: one data node carrying every genome,
 * and under it one refined node per group with that group's file nodes below it. That is the
 * arrangement both rules turn on -- the data node is what a sample's uninformative mass rests on, and
 * the refined node is what a refinement adds between a genome and the species.
 */
public class STModelFixture {
    /** A counted database together with the tree it came from, so a test can address nodes. */
    public static final class Built {
        private final STComposition composition;
        private final SmallTaxTree tree;

        Built(STComposition composition, SmallTaxTree tree) {
            this.composition = composition;
            this.tree = tree;
        }

        /** The composition both classifiers are built on. */
        public STComposition composition() {
            return composition;
        }

        /** Position of the data node, i.e. the one all the genomes hang under. */
        public int dataNode() {
            return tree.getNodeByTaxId("3").getPosition();
        }

        /** Position of the n-th file node, numbered as {@link STModelFixture#build} created them. */
        public int leaf(int n) {
            return tree.getNodeByTaxId(String.valueOf(10 + n)).getPosition();
        }

        /** Position of the n-th refined node, i.e. of the n-th group given to the fixture. */
        public int group(int n) {
            return tree.getNodeByTaxId(String.valueOf(100 + n)).getPosition();
        }
    }

    /**
     * Builds the database.
     *
     * @param groups one map per refined node, leaf name to sequence type, an empty value or
     *               {@code -} marking a genome the scheme does not cover
     * @return the counted composition and its tree
     * @throws IOException if the temporary taxonomy cannot be written
     */
    @SafeVarargs
    public final Built build(Map<String, String>... groups) throws IOException {
        File dir = Files.createTempDirectory("stmodel").toFile();
        dir.deleteOnExit();
        StringBuilder nodes = new StringBuilder();
        StringBuilder names = new StringBuilder();
        // 1 root -> 2 species -> 3 data -> one refined node per group -> one file node per genome
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
        return new Built(new STComposition(tree, new STGroundTruth(mlstFile)), tree);
    }

    /**
     * Convenience for writing a group as alternating name and type.
     *
     * @param pairs leaf name, sequence type, leaf name, sequence type, ...
     * @return the group
     */
    public Map<String, String> leaves(String... pairs) {
        Map<String, String> m = new LinkedHashMap<>();
        for (int i = 0; i < pairs.length; i += 2) {
            m.put(pairs[i], pairs[i + 1]);
        }
        return m;
    }
}
