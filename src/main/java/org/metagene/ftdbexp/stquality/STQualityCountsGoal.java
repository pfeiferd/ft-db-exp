package org.metagene.ftdbexp.stquality;

import it.unimi.dsi.fastutil.objects.Object2LongMap;
import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.probfilter.KMerIndexFilterHelper;
import org.metagene.genestrip.finertree.refseq.AbstractUpdateFastaReader;
import org.metagene.genestrip.goals.refseq.FastaReaderGoal;
import org.metagene.genestrip.goals.refseq.RefSeqFnaFilesDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.probfilter.BlockedBloomFilter;
import org.metagene.genestrip.probfilter.ProbFilter;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AbstractStoreFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Computes the intrinsic quality of a database at the <em>sequence type</em> level, by re-reading the
 * genomes it was filled from and asking, for each stored k-mer, how many of the lineages that could
 * carry it actually do.
 * <p>
 * It is the counterpart of Genestrip-FT's {@code DBQualityCountsGoal} for a database refined
 * <em>below</em> the species rank, where the unit that measure counts - the data taxon - is a single
 * genome and says almost nothing. Sections "Case study" and "Node precision" of the paper give the
 * definition; {@link STCounts} restates it. The whole difference to the data-taxon level is which
 * unit a k-mer is credited to:
 * <ul>
 * <li>the duplicate filter is keyed by {@code (k-mer, sequence type)} rather than by
 *     {@code (k-mer, leaf)}, so a k-mer carried by fifty genomes of one lineage counts once - which is
 *     precisely {@code c_st(a)} rather than {@code c(a)};</li>
 * <li>the denominator is {@code |S_n|}, the number of distinct types below the node, in place of the
 *     number of genomes below it.</li>
 * </ul>
 * Since several genomes now share a unit, {@code |S_n| <= |D_n|} throughout and the measure is the
 * more forgiving of the two - which is the point: a k-mer shared by every genome of one lineage
 * scores one here and {@code 1/|D_n|} there, and it is the former that states what the refinement is
 * for.
 * <p>
 * This goal lives in ft-db-exp and not in genestrip-ft because the ground truth it rests on is not
 * something a Genestrip database knows or could derive: a sequence type comes from an external typing
 * of the genomes, and only this one experiment has one.
 */
public class STQualityCountsGoal extends FastaReaderGoal<Map<String, STCounts>, FTProject>
        implements Goal.LogHeapInfo {
    private final ObjectGoal<AccessionMap, FTProject> accessionMapGoal;
    private final ObjectGoal<Database, FTProject> storeGoal;
    private final STGroundTruth groundTruth;

    private SmallTaxTree tree;
    private KMerStore<SmallTaxTree.SmallTaxIdNode> kMerStore;
    private ProbFilter filter;
    private Map<String, STCounts> map;
    private List<MyFastaReader> readersList;
    /**
     * Tax ids of the nodes the measure applies to, i.e. those at or below a requested tax id. See
     * {@link #inScope}.
     */
    private Set<String> scope;
    /** Per sequence type, the true positives and the pairs found at all; backs the recalls. */
    private long[] typeTp;
    private long[] typeTpPlusFn;

    /**
     * Creates the goal.
     *
     * @param project          the FT project
     * @param key              the goal key
     * @param bundle           the execution context supplying the reader threads
     * @param categoriesGoal   the goal providing the RefSeq categories to read
     * @param taxNodesGoal     the goal providing the tax nodes to be included
     * @param fnaFilesGoal     the goal providing the downloaded genomic fasta files
     * @param additionalGoal   the goal providing additional fasta files mapped to tax nodes
     * @param accessionMapGoal the goal providing the accession-to-tax-node map
     * @param storeGoal        the goal providing the loaded database to measure
     * @param groundTruth      the sequence type of every genome the database was filled from
     * @param deps             further goals this goal depends on
     */
    @SafeVarargs
    public STQualityCountsGoal(FTProject project, GoalKey key, ExecutionContext bundle,
                               ObjectGoal<Set<RefSeqCategory>, FTProject> categoriesGoal,
                               ObjectGoal<Set<TaxTree.TaxIdNode>, FTProject> taxNodesGoal,
                               RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                               ObjectGoal<Map<File, TaxTree.TaxIdNode>, FTProject> additionalGoal,
                               ObjectGoal<AccessionMap, FTProject> accessionMapGoal,
                               ObjectGoal<Database, FTProject> storeGoal,
                               STGroundTruth groundTruth,
                               Goal<FTProject>... deps) {
        super(project, key, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal,
                Goal.append(deps, accessionMapGoal, storeGoal));
        this.storeGoal = storeGoal;
        this.accessionMapGoal = accessionMapGoal;
        this.groundTruth = groundTruth;
    }

    @Override
    protected void doMakeThis() {
        // Without artificial nodes the k-mers sit directly on the tax ids, every tax id is its own
        // leaf, and there is no genome to join a sequence type to.
        if (!booleanConfigValue(GSConfigKey.DATA_NODES)
                && !booleanConfigValue(GSConfigKey.FILE_NODES)
                && !booleanConfigValue(GSConfigKey.ID_NODES)) {
            throw new IllegalStateException("This goal requires data, file or id nodes");
        }
        try {
            tree = storeGoal.get().getTaxTree();
            Object2LongMap<String> stats = storeGoal.get().getStats();
            initScope();

            map = new HashMap<>();
            int typedLeavesInTree = 0;
            int leavesInTree = 0;
            for (SmallTaxTree.SmallTaxIdNode node : tree) {
                boolean leaf = isLeafNode(node);
                STCounts counts = new STCounts(leaf, stats.getOrDefault(node.getTaxId(), 0L));
                if (leaf && inScope(node)) {
                    leavesInTree++;
                    int stIndex = groundTruth.getSTIndex(node.getName());
                    if (stIndex >= 0) {
                        typedLeavesInTree++;
                        counts.addST(stIndex);
                    }
                }
                map.put(node.getTaxId(), counts);
            }
            if (typedLeavesInTree == 0) {
                throw new IllegalStateException("Not one of the " + leavesInTree + " leaves of this database"
                        + " carries a sequence type. The typing keys its rows by the name of the fasta file a"
                        + " genome was read from, which is what a FILE node is named after - so a typing made"
                        + " from extracted per-accession files, or from a different build of the database,"
                        + " joins to nothing.");
            }
            if (getLogger().isInfoEnabled()) {
                getLogger().info("Sequence types: " + groundTruth.getSTCount() + " over " + typedLeavesInTree
                        + " of " + leavesInTree + " leaves ("
                        + String.format("%.1f", 100.0 * typedLeavesInTree / leavesInTree)
                        + " %); the untyped ones contribute no unit and make the measure more generous.");
            }

            // Unioned upwards before the k-mers are read, because the filter is sized from |S_n| -
            // and because the CSV reports it per node, so a node whose subtree holds no typed genome
            // has to be recognisable as having no measure at all rather than a bad one.
            unionSTsUpwards();

            // One entry per (k-mer, sequence type) pair, and a k-mer stored at n can pair with no
            // type that is not in S_n, so this sum is an upper bound and a tight one.
            //
            // The obvious estimate - every k-mer on the path from a typed genome up to the root,
            // summed over the genomes, as DBQualityCountsGoal makes it - is not usable here. It
            // counts the k-mers of the shared ancestors once per genome, which is harmless while a
            // leaf is one of a handful under its taxon and ruinous when it is one of 3,500: on
            // `cdiff' it comes to 227 billion entries against the 12 billion below, and asks for a
            // filter of 265 GB.
            long size = 0;
            for (SmallTaxTree.SmallTaxIdNode node : tree) {
                STCounts counts = map.get(node.getTaxId());
                size += counts.getKmerSumForNode() * counts.getSTCount();
            }
            filter = new BlockedBloomFilter(size);
            if (getLogger().isInfoEnabled()) {
                getLogger().info("At most " + size + " (k-mer, sequence type) pair(s); filter size in MB: "
                        + (filter.getBitSize() / 8 / 1024 / 1024));
            }
            typeTp = new long[groundTruth.getSTCount()];
            typeTpPlusFn = new long[groundTruth.getSTCount()];
            kMerStore = storeGoal.get().convertKMerStore();

            readersList = new ArrayList<>();
            readFastas();

            long entries = 0;
            for (MyFastaReader reader : readersList) {
                entries += reader.entries;
            }
            if (getLogger().isInfoEnabled()) {
                getLogger().info("Filter entries: " + entries);
            }
            if (entries > 2 * size && getLogger().isErrorEnabled()) {
                getLogger().error("Entries exceed filter size by over factor 2. Something went wrong!");
            }

            aggregateRecalls();
            aggregateSubtrees();
            set(map);
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            map = null;
            tree = null;
            kMerStore = null;
            filter = null;
            readersList = null;
            scope = null;
            typeTp = null;
            typeTpPlusFn = null;
            cleanUpThreads();
        }
    }

    /**
     * Works out which nodes the measure applies to: those at or below one of the tax ids the project
     * requests, which for {@code cdiff} is the subtree of 1496.
     * <p>
     * Above them the measure would be defined and meaningless. {@code S_n} is built from an external
     * typing that covers the requested organism and nothing else, so at, say, the phylum every genome
     * of every other organism below it is untyped and contributes no unit - and the paper's caveat,
     * that untyped genomes make {@code |S_v(a)|} smaller and {@code p_st} correspondingly more
     * generous, stops being a caveat and becomes the whole of the number. A k-mer stored at the phylum
     * would be scored against the 196 lineages of one species as though nothing else could carry it.
     * <p>
     * It also keeps the pass affordable: those nodes hold the k-mers this species shares with its
     * relatives, and pairing each of them with every lineage is most of the {@code |S_n|}-weighted
     * bound the filter is built to.
     * <p>
     * An empty request means the project asked for everything, and then nothing is restricted.
     */
    private void initScope() {
        Set<TaxTree.TaxIdNode> requested = taxNodesGoal.get();
        if (requested == null || requested.isEmpty()) {
            scope = null;
            return;
        }
        Set<String> requestedIds = new HashSet<>();
        for (TaxTree.TaxIdNode node : requested) {
            requestedIds.add(node.getTaxId());
        }
        scope = new HashSet<>();
        for (SmallTaxTree.SmallTaxIdNode node : tree) {
            for (SmallTaxTree.SmallTaxIdNode a = node; a != null; a = a.getParent()) {
                if (requestedIds.contains(a.getTaxId())) {
                    scope.add(node.getTaxId());
                    break;
                }
            }
        }
    }

    /**
     * Returns whether the measure applies to the given node.
     *
     * @param node the node to test
     * @return whether it is at or below one of the requested tax ids
     */
    private boolean inScope(SmallTaxTree.SmallTaxIdNode node) {
        return scope == null || scope.contains(node.getTaxId());
    }

    /**
     * Gives every node the read counts of the lineages below it, which is what its two recalls are
     * taken over.
     * <p>
     * A lineage is counted once per node it occurs below, however many of its genomes sit there: it is
     * the unit, and {@code S_n} is a union for the same reason. Membership is read off the set the
     * union already built, so the genomes need not be walked again.
     */
    private void aggregateRecalls() {
        for (SmallTaxTree.SmallTaxIdNode node : tree) {
            STCounts counts = map.get(node.getTaxId());
            for (int st = counts.nextST(0); st >= 0; st = counts.nextST(st + 1)) {
                counts.addTypeRecall(typeTp[st], typeTpPlusFn[st]);
            }
        }
    }

    /**
     * Propagates every leaf's sequence type into all of its ancestors, so that each node's set is
     * {@code S_n}, the types represented by the genomes below it.
     * <p>
     * A union and not a sum, which is the one structural difference to the data-taxon level: two
     * children may carry the same type, and counting per child would report more lineages under a
     * node than there are.
     */
    private void unionSTsUpwards() {
        for (SmallTaxTree.SmallTaxIdNode node : tree) {
            STCounts counts = map.get(node.getTaxId());
            if (counts.isForLeaf() && counts.getSTCount() > 0) {
                // Stops at the scope boundary rather than running to the root: a node above it gets
                // no set, hence no measure and no share of the filter. See initScope().
                for (SmallTaxTree.SmallTaxIdNode a = node.getParent(); a != null && inScope(a);
                     a = a.getParent()) {
                    map.get(a.getTaxId()).unionSTs(counts);
                }
            }
        }
    }

    /**
     * Pools every node's own k-mers into each of its ancestors, which is what the two subtree averages
     * are taken over.
     *
     * @see STCounts#aggregateSubtree
     */
    private void aggregateSubtrees() {
        for (SmallTaxTree.SmallTaxIdNode node : tree) {
            STCounts counts = map.get(node.getTaxId());
            if (counts.getSTCount() > 0) {
                for (SmallTaxTree.SmallTaxIdNode a = node; a != null && inScope(a); a = a.getParent()) {
                    map.get(a.getTaxId()).aggregateSubtree(counts);
                }
            }
        }
    }

    /**
     * Whether the given node is where a genomic file's k-mers come to rest, and therefore a genome
     * that a sequence type can be attached to.
     * <p>
     * The fill nests the artificial nodes - {@code ReworkingStoreFastaReader.reworkNode()} descends a
     * tax id into its DATA child, that into a FILE child, that into an ID child, as far as the three
     * flags are enabled - and reading back, {@code AbstractUpdateFastaReader.updateLeafNode()} walks
     * the same chain from the other end. Both land on the deepest of them, which is what this
     * identifies: an origin-rank node with no origin-rank child. REFINED is deliberately not an origin
     * rank: a refined node is inserted <em>above</em> the origin nodes and is internal in exactly the
     * way a taxonomy node is.
     * <p>
     * This mirrors {@code DBQualityCountsGoal.isLeafNode}, which is protected and cannot be reached
     * from here. It is repeated rather than exposed because the definition belongs to the database
     * layout and this project may not change genestrip-ft; should the layout ever change, both have to
     * move together.
     *
     * @param node the node to test
     * @return whether the node is the deepest artificial node on its branch
     */
    static boolean isLeafNode(SmallTaxTree.SmallTaxIdNode node) {
        if (!isOriginRank(node.getRank())) {
            return false;
        }
        SmallTaxTree.SmallTaxIdNode[] subNodes = node.getSubNodes();
        if (subNodes != null) {
            for (SmallTaxTree.SmallTaxIdNode subNode : subNodes) {
                if (isOriginRank(subNode.getRank())) {
                    return false;
                }
            }
        }
        return true;
    }

    private static boolean isOriginRank(Rank rank) {
        return Rank.DATA.equals(rank) || Rank.FILE.equals(rank) || Rank.ID.equals(rank);
    }

    @Override
    protected AbstractStoreFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
        MyFastaReader reader = new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES),
                taxNodesGoal.get(),
                isIncludeRefSeqFna() ? accessionMapGoal.get() : null,
                intConfigValue(GSConfigKey.KMER_SIZE),
                intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID),
                (Rank) configValue(GSConfigKey.MAX_GENOMES_PER_TAXID_RANK),
                longConfigValue(GSConfigKey.MAX_KMERS_PER_TAXID),
                intConfigValue(GSConfigKey.MAX_DUST),
                intConfigValue(GSConfigKey.KMER_SAMPLING),
                booleanConfigValue(GSConfigKey.ASSEMBLY_ACCESSIONS_ONLY),
                regionsPerTaxid,
                booleanConfigValue(GSConfigKey.ENABLE_LOWERCASE_BASES));
        readersList.add(reader);
        return reader;
    }

    /**
     * Reader that credits each k-mer it finds in the database to the sequence type of the genome it
     * was read from, once per type.
     */
    protected class MyFastaReader extends AbstractUpdateFastaReader {
        /** Number of distinct (k-mer, sequence type) pairs this reader added to the filter. */
        protected long entries;

        public MyFastaReader(int bufferSize, Set<TaxTree.TaxIdNode> taxNodes, AccessionMap accessionMap,
                             int k, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId,
                             int maxDust, int kMerSampling, boolean assemblyAccessionsOnly,
                             StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases) {
            super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank,
                    maxKmersPerTaxId, maxDust, kMerSampling, assemblyAccessionsOnly, regionsPerTaxid,
                    enableLowerCaseBases, booleanConfigValue(GSConfigKey.ID_NODES),
                    booleanConfigValue(GSConfigKey.FILE_NODES), booleanConfigValue(GSConfigKey.DATA_NODES));
        }

        @Override
        protected SmallTaxTree getTree() {
            return tree;
        }

        @Override
        protected boolean handleStore(long kmer) {
            if (leafNode == null) {
                return false;
            }
            int stIndex = groundTruth.getSTIndex(leafNode.getName());
            if (stIndex < 0) {
                // An untyped genome contributes no unit, so it can neither raise nor lower c_st.
                return false;
            }
            SmallTaxTree.SmallTaxIdNode storedNode = kMerStore.getLong(kmer, null);
            if (storedNode == null || !inScope(storedNode)) {
                // Out of scope is not a miss: the measure simply does not reach above the requested
                // tax ids, so such a k-mer is left out of both the precision and the recall rather
                // than counted as a k-mer this lineage failed to claim.
                return false;
            }
            // Keyed by the type and not by the genome: that is the whole of what separates c_st from
            // c. Fifty genomes of one lineage carrying this k-mer add one to the tally, not fifty.
            if (!filter.putLong(KMerIndexFilterHelper.combine(kmer, stIndex))) {
                return false;
            }
            entries++;
            // Path correctness says the genomes carrying a k-mer all lie below the node it is stored
            // at, so this walk succeeds for every genome of the type and it does not matter which of
            // them won the filter. Where it fails, the database claims the k-mer for a branch this
            // genome is not on - the false negative the recall is there to count.
            SmallTaxTree.SmallTaxIdNode pathNode = leafNode;
            while (pathNode != null && pathNode != storedNode) {
                pathNode = pathNode.getParent();
            }
            boolean onPath = pathNode == storedNode;
            synchronized (typeTpPlusFn) {
                typeTpPlusFn[stIndex]++;
                if (onPath) {
                    typeTp[stIndex]++;
                }
            }
            if (onPath) {
                STCounts nodeCounts = map.get(storedNode.getTaxId());
                synchronized (nodeCounts) {
                    nodeCounts.incTpForNodePrecision();
                }
                return true;
            }
            return false;
        }
    }
}
