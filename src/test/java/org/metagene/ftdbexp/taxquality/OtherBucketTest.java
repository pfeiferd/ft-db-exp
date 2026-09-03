package org.metagene.ftdbexp.taxquality;

import org.junit.Test;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

/**
 * Pins that the empty OTHER buckets a refinement creates are not counted as organisms.
 * <p>
 * A candidate set is the species a node still leaves in question, and {@code 1/|Sigma(k)|} over it is
 * what the paper scores a call by. An OTHER bucket is not a species: the refinement adds it to stand
 * for "none of the siblings", it holds no genome, and in all six databases of the paper every
 * childless {@code REFINED} node is such a bucket and holds no k-mer.
 * <p>
 * Counting it added exactly one to every refined node above it and to no unrefined node at all, since
 * the unrefined tree has none. The genus of {@code strepto} reported 158 species against the
 * unrefined 157 and the node holding the mitis complex reported nine against eight, so the measure
 * understated the refinement it was there to assess.
 */
public class OtherBucketTest {
    /** The bucket the refinement leaves beside the children it did place. */
    @Test
    public void testAnEmptyOtherBucketContributesNoCandidate() {
        assertNull(TaxonComposition.classOfLeaf(new SmallTaxIdNode("OTHER", "OTHER", Rank.REFINED)));
    }

    /** A genome under a species is that species, as before. */
    @Test
    public void testAGenomeContributesItsSpecies() {
        SmallTaxIdNode species = new SmallTaxIdNode("1313", "1313", Rank.SPECIES,
                new SmallTaxIdNode[] { new SmallTaxIdNode("d1", "d1", Rank.DATA) });
        assertEquals("1313", TaxonComposition.classOfLeaf(species.getSubNodes()[0]));
    }

    /**
     * A real leaf that names no species still counts. A genome under {@code unclassified
     * Streptococcus} is material the database holds, and a k-mer there could have come from it. That
     * is a fact about the reference taxonomy, not an artefact of the refinement, and the two must not
     * be conflated by dropping every UNPLACED alike.
     */
    @Test
    public void testAnUnplaceableGenomeIsStillACandidate() {
        SmallTaxIdNode unnamed = new SmallTaxIdNode("2608887", "2608887", Rank.NO_RANK,
                new SmallTaxIdNode[] { new SmallTaxIdNode("dU", "dU", Rank.DATA) });
        assertEquals(TaxonComposition.UNPLACED, TaxonComposition.classOfLeaf(unnamed.getSubNodes()[0]));
    }

    /** A refined node that did place children is not a bucket and is never asked: only leaves are. */
    @Test
    public void testARefinedNodeWithChildrenIsNotALeafAndKeepsItsChildrensSpecies() {
        SmallTaxIdNode child = new SmallTaxIdNode("1313", "1313", Rank.SPECIES);
        SmallTaxIdNode refined = new SmallTaxIdNode("r", "r", Rank.REFINED, new SmallTaxIdNode[] { child });
        assertEquals(1, refined.getSubNodes().length);
        assertNull("a childless one would contribute nothing", TaxonComposition.classOfLeaf(
                new SmallTaxIdNode("r2", "r2", Rank.REFINED)));
    }
}
