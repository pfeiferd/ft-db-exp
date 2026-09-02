package org.metagene.ftdbexp.taxquality;

import org.junit.Test;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertSame;

/**
 * Pins the node the nearest-node columns of {@code <db>_<key>_taxoncall.csv} report.
 * <p>
 * The question they answer is how close to a named reference organism a sample's evidence got, and
 * how many species are still open where it stopped. The ranked votes beside them cannot answer it:
 * a refined node bears no name of the reference taxonomy, so it never wins a vote that has to name a
 * taxon, and every read reaching one is counted as unplaced instead.
 * <p>
 * The lineage below is the one {@code strepto} actually has. In the unrefined database
 * {@code S. pneumoniae} is a direct child of the genus, so the answer can only ever be the species
 * or the genus. A refinement inserts a ladder between them, and it is the rungs of that ladder that
 * these columns exist to see.
 */
public class NearestNodeTest {
    private final SmallTaxIdNode pneumoStrain = new SmallTaxIdNode("170187", "170187", Rank.STRAIN);
    /** The refined strain cluster the refinement moves most of the species' own k-mers into. */
    private final SmallTaxIdNode strainCluster =
            new SmallTaxIdNode("000578", "000578", null, new SmallTaxIdNode[] { pneumoStrain });
    private final SmallTaxIdNode pneumo =
            new SmallTaxIdNode("1313", "1313", Rank.SPECIES, new SmallTaxIdNode[] { strainCluster });
    /** The refined node holding the mitis complex: eight species, the pneumococcus among them. */
    private final SmallTaxIdNode mitisComplex =
            new SmallTaxIdNode("0001210", "0001210", null, new SmallTaxIdNode[] { pneumo });
    private final SmallTaxIdNode genus =
            new SmallTaxIdNode("1301", "1301", Rank.GENUS, new SmallTaxIdNode[] { mitisComplex });

    private static List<Long> counts(long... values) {
        List<Long> out = new ArrayList<>();
        for (long v : values) {
            out.add(v);
        }
        return out;
    }

    /** The plain case: the sample's evidence sits at the species itself. */
    @Test
    public void testSignalAtTheSpeciesAnswersWithTheSpecies() {
        TaxonCallReport.Near near = TaxonCallReport.nearest(
                pneumo, Arrays.asList(pneumo, genus), counts(500, 3000));
        assertSame(pneumo, near.node);
        assertEquals(500, near.at);
        assertEquals(500, near.below);
    }

    /**
     * The case that makes the subtree search necessary, and it is the common one after a refinement
     * rather than an edge case. The refinement moves 1,174,677 of the pneumococcus's 1,273,715 own
     * k-mers down into a strain cluster, leaving 99,038 at the species node -- and it can leave
     * nothing there at all. Walking upwards from an emptied species node would climb past it to the
     * mitis complex and report eight species open for a sample that had in fact reached the
     * organism: the refinement scored worst exactly where it works best.
     */
    @Test
    public void testSignalOnlyBelowTheSpeciesStillAnswersWithTheSpecies() {
        TaxonCallReport.Near near = TaxonCallReport.nearest(
                pneumo, Arrays.asList(strainCluster, mitisComplex, genus), counts(700, 400, 3000));
        assertSame("an emptied species node must not send the walk upwards", pneumo, near.node);
        assertEquals("nothing sits at the species itself", 0, near.at);
        assertEquals("but the strain cluster below it did", 700, near.below);
    }

    /** Reaching a strain below the cluster is still reaching the species. */
    @Test
    public void testSignalAtAStrainAnswersWithTheSpecies() {
        TaxonCallReport.Near near = TaxonCallReport.nearest(
                pneumo, Arrays.asList(pneumoStrain, genus), counts(120, 3000));
        assertSame(pneumo, near.node);
        assertEquals(0, near.at);
        assertEquals(120, near.below);
    }

    /**
     * Nothing reached the organism, so the answer is the deepest ancestor that carries anything.
     * The two counts differ here, and that difference is the whole reason both are reported: what
     * sits <em>at</em> the mitis complex is the material no refinement can place, while what sits at
     * or below it is everything that got that far.
     */
    @Test
    public void testWithoutSignalAtTheOrganismTheDeepestAncestorAnswers() {
        SmallTaxIdNode oralis = new SmallTaxIdNode("1303", "1303", Rank.SPECIES);
        List<SmallTaxIdNode> nodes = Arrays.asList(mitisComplex, genus);
        TaxonCallReport.Near near = TaxonCallReport.nearest(pneumo, nodes, counts(400, 3000));
        assertSame(mitisComplex, near.node);
        assertEquals("the undecidable residue", 400, near.at);
        assertEquals("and nothing below it", 400, near.below);
        assertEquals("a sibling is not on the lineage and contributes nothing", 0,
                TaxonCallReport.countAtOrBelow(pneumo, Arrays.asList(oralis), counts(9000)));
    }

    /** The unrefined tree has no ladder: the genus is the only other answer available. */
    @Test
    public void testWithoutARefinementTheAnswerFallsToTheGenus() {
        SmallTaxIdNode flatPneumo = new SmallTaxIdNode("1313", "1313", Rank.SPECIES);
        SmallTaxIdNode flatGenus =
                new SmallTaxIdNode("1301", "1301", Rank.GENUS, new SmallTaxIdNode[] { flatPneumo });
        TaxonCallReport.Near near = TaxonCallReport.nearest(
                flatPneumo, Arrays.asList(flatGenus), counts(3000));
        assertSame(flatGenus, near.node);
        assertEquals(3000, near.at);
    }

    /**
     * No signal anywhere on the lineage. The columns are then left empty rather than zeroed: a
     * sample that put nothing there did not put none there, it was never asked.
     */
    @Test
    public void testNoSignalOnTheLineageAnswersWithNothing() {
        SmallTaxIdNode elsewhere = new SmallTaxIdNode("1314", "1314", Rank.SPECIES);
        assertNull(TaxonCallReport.nearest(pneumo, Arrays.asList(elsewhere), counts(4000)));
    }

    /** Without a target the columns are not measured at all, which is the default for every other database. */
    @Test
    public void testNoTargetAnswersWithNothing() {
        assertNull(TaxonCallReport.nearest(null, Arrays.asList(genus), counts(3000)));
    }
}
