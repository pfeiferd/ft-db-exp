package org.metagene.ftdbexp.taxquality;

import org.junit.Test;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

import static org.junit.Assert.assertEquals;

/**
 * Pins the vote to what it is meant to be: read counts carried up to the species they belong to, and
 * the species ranked by what they gathered.
 * <p>
 * Two things are easy to get wrong here and both were wrong at some point. Counts must roll
 * <em>up</em> into the species, or a species entering the database with 153 genomes fields 153
 * candidates with a fraction of the reads each and loses to a draft genome with one. And counts that
 * reach no species must not be handed to one: a read left at a genus, a refined node or the
 * {@code unclassified Streptococcus} bucket has not been placed at a species, and crediting it to
 * every species below would let a single read outweigh fifty.
 */
public class SpeciesVoteTest {
    // genus 1301
    //   +- species 1313 -- strain 170187 -- DATA d1
    //   |                  strain 171101 -- DATA d2
    //   +- species 1314 -- DATA d3
    //   +- no rank 2608887 (unclassified) -- species 1306 -- DATA d4
    private final SmallTaxIdNode d1 = leaf("d1");
    private final SmallTaxIdNode d2 = leaf("d2");
    private final SmallTaxIdNode d3 = leaf("d3");
    private final SmallTaxIdNode d4 = leaf("d4");
    private final SmallTaxIdNode strainA = inner("170187", Rank.STRAIN, d1);
    private final SmallTaxIdNode strainB = inner("171101", Rank.STRAIN, d2);
    private final SmallTaxIdNode pneumo = inner("1313", Rank.SPECIES, strainA, strainB);
    private final SmallTaxIdNode pyo = inner("1314", Rank.SPECIES, d3);
    private final SmallTaxIdNode sp = inner("1306", Rank.SPECIES, d4);
    private final SmallTaxIdNode bucket = inner("2608887", Rank.NO_RANK, sp);
    private final SmallTaxIdNode genus = inner("1301", Rank.GENUS, pneumo, pyo, bucket);

    private static SmallTaxIdNode leaf(String taxId) {
        return new SmallTaxIdNode(taxId, taxId, Rank.DATA);
    }

    private static SmallTaxIdNode inner(String taxId, Rank rank, SmallTaxIdNode... subNodes) {
        return new SmallTaxIdNode(taxId, taxId, rank, subNodes);
    }

    /** The two strains of one species are one candidate, not two. */
    @Test
    public void testStrainsRollUpIntoTheirSpecies() {
        assertEquals("1313", TaxonComposition.speciesOf(d1));
        assertEquals("1313", TaxonComposition.speciesOf(d2));
    }

    /** A species below the unclassified bucket is still a species and keeps its own reads. */
    @Test
    public void testASpeciesBelowTheBucketKeepsItsReads() {
        assertEquals("1306", TaxonComposition.speciesOf(d4));
    }

    /** The bucket's own reads reached no species and must stay where they are. */
    @Test
    public void testTheBucketItselfReachesNoSpecies() {
        assertEquals(TaxonComposition.UNPLACED, TaxonComposition.speciesOf(bucket));
    }

    /** Neither do the genus's own reads. */
    @Test
    public void testTheGenusReachesNoSpecies() {
        assertEquals(TaxonComposition.UNPLACED, TaxonComposition.speciesOf(genus));
    }

    /**
     * What a refinement inserts sits above the species, so its reads count as unplaced until the
     * refinement pushes them further down. That is why the unplaced figure is the one to watch.
     */
    @Test
    public void testARefinedNodeReachesNoSpecies() {
        SmallTaxIdNode refined = inner("r1", Rank.REFINED, pneumo);
        inner("1301b", Rank.GENUS, refined);
        assertEquals(TaxonComposition.UNPLACED, TaxonComposition.speciesOf(refined));
    }
}
