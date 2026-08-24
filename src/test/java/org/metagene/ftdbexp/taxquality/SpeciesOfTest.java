package org.metagene.ftdbexp.taxquality;

import org.junit.Test;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

import static org.junit.Assert.assertEquals;

/**
 * Pins which node of a lineage counts as a genome's class.
 * <p>
 * It must be the species, and the rule that finds it cannot stop at the first rank it recognises:
 * {@link Rank#STRAIN} sits between {@link Rank#SPECIES} and {@link Rank#FORMA_SPECIALIS}, so a walk
 * upwards meets the strain first. Stopping there made each of the 153 strains of
 * {@code Streptococcus pneumoniae} a class of its own, which inflated every candidate count and left
 * a correct call on {@code Streptococcus pyogenes} scoring zero, because the reference names a
 * species and the candidate set held only strains.
 * <p>
 * The lineages below are the ones the {@code strepto} database actually contains, read off its
 * quality CSV.
 */
public class SpeciesOfTest {
    /** Builds a lineage from the root downwards and returns its deepest node. */
    private static SmallTaxIdNode lineage(Object... rankThenTaxId) {
        SmallTaxIdNode child = null;
        for (int i = rankThenTaxId.length - 2; i >= 0; i -= 2) {
            Rank rank = (Rank) rankThenTaxId[i];
            String taxId = (String) rankThenTaxId[i + 1];
            child = child == null
                    ? new SmallTaxIdNode(taxId, taxId, rank)
                    : new SmallTaxIdNode(taxId, taxId, rank, new SmallTaxIdNode[] { child });
        }
        SmallTaxIdNode deepest = child;
        while (deepest.getSubNodes() != null && deepest.getSubNodes().length > 0) {
            deepest = deepest.getSubNodes()[0];
        }
        return deepest;
    }

    /** The common case: the fill files a genome under a strain of a species. */
    @Test
    public void testAStrainResolvesToItsSpecies() {
        assertEquals("1313", TaxonComposition.speciesOf(lineage(
                Rank.GENUS, "1301", Rank.SPECIES, "1313", Rank.STRAIN, "170187", Rank.DATA, "d1")));
    }

    /** No strain in between: the species is met directly. */
    @Test
    public void testASpeciesIsItsOwnClass() {
        assertEquals("1313", TaxonComposition.speciesOf(lineage(
                Rank.GENUS, "1301", Rank.SPECIES, "1313", Rank.DATA, "d1")));
    }

    /**
     * The {@code unclassified Streptococcus} branch: a `no rank' bucket sits between the genus and
     * species that carry no proper name. The bucket must be stepped over, not mistaken for an end.
     */
    @Test
    public void testANoRankBucketDoesNotHideTheSpeciesBelowIt() {
        assertEquals("1306", TaxonComposition.speciesOf(lineage(
                Rank.GENUS, "1301", Rank.NO_RANK, "2608887", Rank.SPECIES, "1306", Rank.DATA, "d1")));
    }

    /** An isolate sits below the species but outside the rank window, and must not end the walk. */
    @Test
    public void testAnIsolateIsSteppedOver() {
        assertEquals("1313", TaxonComposition.speciesOf(lineage(
                Rank.GENUS, "1301", Rank.SPECIES, "1313", Rank.ISOLATE, "i1", Rank.DATA, "d1")));
    }

    /** All three artificial nodes of the fill, nested as ReworkingStoreFastaReader nests them. */
    @Test
    public void testTheArtificialNodesAreSteppedOver() {
        assertEquals("1314", TaxonComposition.speciesOf(lineage(
                Rank.GENUS, "1301", Rank.SPECIES, "1314", Rank.DATA, "d1", Rank.FILE, "f1", Rank.ID, "i1")));
    }

    /** A refined node is inserted above the species and must not be taken for one. */
    @Test
    public void testARefinedNodeIsNotASpecies() {
        assertEquals("1313", TaxonComposition.speciesOf(lineage(
                Rank.GENUS, "1301", Rank.REFINED, "r1", Rank.SPECIES, "1313", Rank.DATA, "d1")));
    }

    /** Nothing at species rank anywhere: the genome is unplaced rather than filed under its genus. */
    @Test
    public void testALineageWithoutASpeciesIsUnplaced() {
        assertEquals(TaxonComposition.UNPLACED, TaxonComposition.speciesOf(lineage(
                Rank.GENUS, "1301", Rank.NO_RANK, "2608887", Rank.DATA, "d1")));
    }

    /** The walk stops at the genus and does not wander on to a species elsewhere in the tree. */
    @Test
    public void testTheWalkStopsAboveTheSpecies() {
        assertEquals(TaxonComposition.UNPLACED, TaxonComposition.speciesOf(lineage(
                Rank.SPECIES, "9999", Rank.GENUS, "1301", Rank.DATA, "d1")));
    }
}
