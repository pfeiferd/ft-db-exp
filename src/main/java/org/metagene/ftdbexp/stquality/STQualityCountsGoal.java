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
import org.metagene.genestrip.store.RadixKMerStore;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;

import java.io.File;
import java.util.ArrayDeque;
import java.util.Deque;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
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
    /** K-mers a reader buffers before one batched store lookup; the value AbstractKMerIndexGoal uses. */
    private static final int BATCH_SIZE = 128;

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
     * {@link #inScope}. Used while the per-node arrays below are being built and not on the reading
     * path, which addresses a node by its position.
     */
    private Set<String> scope;
    /** Per sequence type, the true positives and the pairs found at all; backs the recalls. */
    private long[] typeTp;
    private long[] typeTpPlusFn;

    // Everything the reading path needs about a node, indexed by its dense position rather than
    // looked up by tax id. A HashMap keyed by a String costs a string hash per access, and these are
    // accessed a few billion times; SmallTaxTree numbers its nodes densely for exactly this purpose
    // (see SmallTaxTree#getNodeCount).
    private int nodeCount;
    private SmallTaxTree.SmallTaxIdNode[] nodeByPos;
    private STCounts[] countsByPos;
    private boolean[] inScopeByPos;
    /** The sequence type of the leaf at this position, or {@code -1} where there is none. */
    private int[] stByPos;
    /** Number of distinct sequence types, i.e. the width of the per-type tallies. */
    private int typeCount;

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
        // Data nodes, and not merely one of the three kinds of artificial node. What a leaf is has to
        // agree between the fill and this measure, and `dataNodes' is what makes that agreement
        // simple: with it on, ReworkingStoreFastaReader.reworkNode() files every genome at a DATA
        // node or deeper, so no taxonomy node ever holds a genome's k-mers and a node without
        // children is exactly a node the fill filed a genome at. isLeafNode() is then one test and
        // needs to know nothing about ranks.
        //
        // Requiring it is not what used to shut this goal out of a database refined below the
        // species -- that was the second half of the old guard, which REJECTED `fileNodes' while
        // only a DATA node could be recognised as a leaf. Both halves went at once; only the first
        // is coming back. `fileNodes' and `idNodes' stay free, and the `cdiff' project of
        // ft-db-exp2, which needs file nodes because the taxonomy supplies no children below the
        // species, has data nodes on as every project here does.
        if (!booleanConfigValue(GSConfigKey.DATA_NODES)) {
            throw new IllegalStateException("This goal requires data nodes (dataNodes=true)");
        }
        try {
            tree = storeGoal.get().getTaxTree();
            Object2LongMap<String> stats = storeGoal.get().getStats();
            initScope();

            map = new HashMap<>();
            nodeCount = tree.getNodeCount();
            nodeByPos = new SmallTaxTree.SmallTaxIdNode[nodeCount];
            countsByPos = new STCounts[nodeCount];
            inScopeByPos = new boolean[nodeCount];
            stByPos = new int[nodeCount];
            Arrays.fill(stByPos, -1);
            int typedLeavesInTree = 0;
            int leavesInTree = 0;
            for (SmallTaxTree.SmallTaxIdNode node : tree) {
                boolean leaf = isLeafNode(node);
                STCounts counts = new STCounts(leaf, stats.getOrDefault(node.getTaxId(), 0L));
                int pos = node.getPosition();
                nodeByPos[pos] = node;
                countsByPos[pos] = counts;
                inScopeByPos[pos] = inScope(node);
                if (leaf && inScopeByPos[pos]) {
                    leavesInTree++;
                    // The one join to the external typing, done once per leaf here so that the
                    // reading path never hashes a file name again.
                    int stIndex = groundTruth.getSTIndex(node.getName());
                    if (stIndex >= 0) {
                        typedLeavesInTree++;
                        stByPos[pos] = stIndex;
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
            // Not the estimate DBQualityCountsGoal makes - every k-mer on the path from a genome up
            // to the root, summed over the genomes. That one is exactly right for *its* unit: summed
            // over leaves it comes to the same thing as this sum with |D_n| in place of |S_n|, and a
            // (k-mer, genome) pair is what it counts. Here the unit is the lineage, several genomes
            // share one, and |S_n| <= |D_n| throughout - on `cdiff' 196 types against 3,638 typed
            // genomes, which is 227 billion entries and a 265 GB filter against the 3.9 billion and
            // 4.5 GB this sum asks for.
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
            typeCount = groundTruth.getSTCount();
            kMerStore = storeGoal.get().convertKMerStore();

            readersList = new ArrayList<>();
            readFastas();

            long entries = 0;
            long pathViolations = 0;
            long[] probe = new long[4];
            for (MyFastaReader reader : readersList) {
                // The readers are done, so whatever the last (partial) batch still holds is looked up
                // and counted here, single-threaded, before their tallies are merged.
                reader.flushBatch();
                reader.mergeInto(typeTp, typeTpPlusFn, countsByPos);
                entries += reader.entries;
                pathViolations += reader.getPathViolations();
                long[] p = reader.getProbeCounters();
                for (int i = 0; i < probe.length; i++) {
                    probe[i] += p[i];
                }
            }
            if (getLogger().isInfoEnabled()) {
                getLogger().info("Filter entries: " + entries);
                getLogger().info("Where the k-mers went: " + probe[0] + " offered by the readers, "
                        + probe[1] + " found in the store (" + (probe[0] - probe[1]) + " missed), "
                        + probe[2] + " out of scope, " + probe[3] + " already seen for their type, "
                        + entries + " counted.");
            }
            // A recall of one alongside a non-zero count here means the recall is right about the
            // lineages and silent about the genomes: some genome carries a k-mer the database stores
            // off its path, and the per-type dedup handed the check to a sibling that was on it. The
            // placement in `ftupdatedb' is what puts a k-mer off a path - see the OTHER slot in
            // `UpdateStoreGoal.coversAll' - so a database that reports this is one to distrust before
            // any precision figure is read from it.
            if (pathViolations > 0 && getLogger().isWarnEnabled()) {
                getLogger().warn(pathViolations + " genome occurrence(s) whose k-mer is not stored on"
                        + " their path, found among the (k-mer, sequence type) pairs the filter had"
                        + " already seen and the recall therefore never checked. The recall below is"
                        + " per sequence type and stays 1 in this case; path correctness per genome"
                        + " does NOT hold for this database.");
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
            nodeByPos = null;
            countsByPos = null;
            inScopeByPos = null;
            stByPos = null;
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
        for (int pos = 0; pos < nodeCount; pos++) {
            STCounts counts = countsByPos[pos];
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
        for (int pos = 0; pos < nodeCount; pos++) {
            STCounts counts = countsByPos[pos];
            if (counts.isForLeaf() && counts.getSTCount() > 0) {
                // Stops at the scope boundary rather than running to the root: a node above it gets
                // no set, hence no measure and no share of the filter. See initScope().
                for (SmallTaxTree.SmallTaxIdNode a = nodeByPos[pos].getParent();
                     a != null && inScopeByPos[a.getPosition()]; a = a.getParent()) {
                    countsByPos[a.getPosition()].unionSTs(counts);
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
        for (int pos = 0; pos < nodeCount; pos++) {
            STCounts counts = countsByPos[pos];
            if (counts.getSTCount() > 0) {
                for (SmallTaxTree.SmallTaxIdNode a = nodeByPos[pos];
                     a != null && inScopeByPos[a.getPosition()]; a = a.getParent()) {
                    countsByPos[a.getPosition()].aggregateSubtree(counts);
                }
            }
        }
    }

    /**
     * Whether the given node is a leaf, i.e. one the database fill filed a genome at.
     * <p>
     * With data nodes on -- which {@link #doMakeThis()} requires -- that is exactly a node without
     * children: {@code ReworkingStoreFastaReader.reworkNode()} files every genome at a DATA node or
     * at something below it, so no taxonomy node ever holds a genome's k-mers, and whatever the
     * refinement inserts between a node and its original children gives that node children and so
     * keeps it internal.
     * <p>
     * Asking about the children's <em>ranks</em> instead is what got this wrong: after a refinement
     * a data node's file nodes are no longer its children, the data node passed for a leaf, and its
     * k-mers -- the ones a refinement has the most to gain on -- dropped out of every average
     * restricted to what sits above the data taxa. On cdiff that alone lifted the reported sp* from
     * 0.236 to 0.338 with no k-mer moving.
     * <p>
     * This mirrors {@code DBQualityCountsGoal.isLeafNode}, which is protected and cannot be reached
     * from here. It is repeated rather than exposed because the definition belongs to the database
     * layout and this project may not change genestrip-ft; should the layout ever change, both have
     * to move together.
     *
     * @param node the node to test
     * @return whether the fill files genomes at this node
     */
    static boolean isLeafNode(SmallTaxTree.SmallTaxIdNode node) {
        SmallTaxTree.SmallTaxIdNode[] subNodes = node.getSubNodes();
        return subNodes == null || subNodes.length == 0;
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
    protected class MyFastaReader extends AbstractUpdateFastaReader
            implements RadixKMerStore.BatchValueConsumer<SmallTaxTree.SmallTaxIdNode> {
        /** Number of distinct (k-mer, sequence type) pairs this reader added to the filter. */
        protected long entries;

        /**
         * Buffers for the batched store lookup, or {@code null} when it cannot be used. The lookup is
         * memory-latency bound and a batch lets many of its cache misses overlap - see
         * {@link RadixKMerStore#getBatch}.
         */
        private final RadixKMerStore.BatchBuffers batch;
        /** This reader's own tallies, merged by {@link #mergeInto} once every reader has finished. */
        private final long[] readerTypeTp;
        private final long[] readerTypeTpPlusFn;
        private final long[] readerTpForNode;
        /**
         * Genome occurrences whose k-mer the database does not store on their path, counted only for
         * the pairs the per-type filter had already seen - see {@link #count}. Never enters the
         * recall; it exists to say whether the recall's per-type unit is hiding anything.
         */
        private long readerPathViolations;
        // Diagnostic: where a k-mer is lost between being read and being counted. The reading pass
        // has been seen to end with `Filter entries: 0' while the store answered every lookup a probe
        // put to it, so the loss is at one of these four points and guessing which has not worked.
        private long probeOffered;
        private long probeAccepted;
        private long probeOutOfScope;
        private long probeDuplicate;
        /** The current region's leaf, resolved once per region rather than once per k-mer. */
        private SmallTaxTree.SmallTaxIdNode cachedLeaf;
        private int cachedLeafPos = -1;

        public MyFastaReader(int bufferSize, Set<TaxTree.TaxIdNode> taxNodes, AccessionMap accessionMap,
                             int k, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId,
                             int maxDust, int kMerSampling, boolean assemblyAccessionsOnly,
                             StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases) {
            super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank,
                    maxKmersPerTaxId, maxDust, kMerSampling, assemblyAccessionsOnly, regionsPerTaxid,
                    enableLowerCaseBases, booleanConfigValue(GSConfigKey.ID_NODES),
                    booleanConfigValue(GSConfigKey.FILE_NODES), booleanConfigValue(GSConfigKey.DATA_NODES));
            readerTypeTp = new long[typeCount];
            readerTypeTpPlusFn = new long[typeCount];
            readerTpForNode = new long[nodeCount];
            // Batched only while no per-taxon limit binds: a batched k-mer is counted after
            // handleStore() has returned, so its return value can no longer report the k-mer, and that
            // value feeds the per-region counters those limits are enforced from.
            boolean unlimited = maxGenomesPerTaxId == Integer.MAX_VALUE && maxKmersPerTaxId == Long.MAX_VALUE;
            batch = (kMerStore instanceof RadixKMerStore && unlimited)
                    ? new RadixKMerStore.BatchBuffers(BATCH_SIZE) : null;
        }

        @Override
        protected SmallTaxTree getTree() {
            return tree;
        }

        /**
         * Adds this reader's tallies to the shared ones. Called once, after every reader has finished,
         * so nothing here needs a lock - which is the point: locking per pair put every thread on the
         * same two monitors for the nodes they all touch.
         *
         * @param tpTarget       per-type true positives to add to
         * @param tpPlusFnTarget per-type pairs to add to
         * @param nodeTarget     per-node tallies to add to
         */
        long getPathViolations() {
            return readerPathViolations;
        }

        long[] getProbeCounters() {
            return new long[] { probeOffered, probeAccepted, probeOutOfScope, probeDuplicate };
        }

        void mergeInto(long[] tpTarget, long[] tpPlusFnTarget, STCounts[] nodeTarget) {
            for (int i = 0; i < typeCount; i++) {
                tpTarget[i] += readerTypeTp[i];
                tpPlusFnTarget[i] += readerTypeTpPlusFn[i];
            }
            for (int pos = 0; pos < nodeCount; pos++) {
                if (readerTpForNode[pos] != 0) {
                    nodeTarget[pos].addTpForNodePrecision(readerTpForNode[pos]);
                }
            }
        }

        /**
         * Looks the buffered k-mers up in one batch and counts those the database holds.
         */
        protected void flushBatch() {
            if (batch != null && !batch.isEmpty()) {
                ((RadixKMerStore<SmallTaxTree.SmallTaxIdNode>) kMerStore).getBatch(batch, this);
            }
        }

        /**
         * Resolves the region's leaf as usual and remembers its position, so that the reading path
         * costs an indexed read per k-mer instead of hashing the leaf's file name in the typing.
         */
        @Override
        protected void updateLeafNode() {
            super.updateLeafNode();
            if (leafNode != cachedLeaf) {
                cachedLeaf = leafNode;
                cachedLeafPos = leafNode == null ? -1 : leafNode.getPosition();
            }
        }

        @Override
        protected boolean handleStore(long kmer) {
            // An untyped genome contributes no unit, so it can neither raise nor lower c_st - and
            // stByPos says so without a lookup.
            if (cachedLeafPos < 0 || stByPos[cachedLeafPos] < 0) {
                return false;
            }
            probeOffered++;
            if (batch != null) {
                if (batch.add(kmer, cachedLeafPos)) {
                    flushBatch();
                }
                return false;
            }
            SmallTaxTree.SmallTaxIdNode storedNode = kMerStore.getLong(kmer, null);
            if (storedNode == null) {
                return false;
            }
            return count(kmer, cachedLeafPos, storedNode);
        }

        /**
         * Counts one k-mer of a flushed batch, which by construction the database holds.
         *
         * @param kmer       the k-mer that was looked up
         * @param leafPos    the position of the leaf it was read in, as buffered with it
         * @param storedNode the node the database stores it at
         */
        @Override
        public void accept(long kmer, int leafPos, SmallTaxTree.SmallTaxIdNode storedNode) {
            count(kmer, leafPos, storedNode);
        }

        private boolean count(long kmer, int leafPos, SmallTaxTree.SmallTaxIdNode storedNode) {
            probeAccepted++;
            int storedPos = storedNode.getPosition();
            if (!inScopeByPos[storedPos]) {
                probeOutOfScope++;
                // Out of scope is not a miss: the measure simply does not reach above the requested
                // tax ids, so such a k-mer is left out of both the precision and the recall rather
                // than counted as a k-mer this lineage failed to claim.
                return false;
            }
            int stIndex = stByPos[leafPos];
            // Keyed by the type and not by the genome: that is the whole of what separates c_st from
            // c. Fifty genomes of one lineage carrying this k-mer add one to the tally, not fifty.
            if (!filter.putLong(KMerIndexFilterHelper.combine(kmer, stIndex))) {
                // Another genome of this type already brought the pair, and only that one had its path
                // checked below. The recall is therefore blind here, and blind in a way that matters:
                // path correctness holds per GENOME, and the assumption that it then holds for every
                // genome of the type is the very thing the recall is meant to establish. Whichever
                // genome won the filter is a race between the reader threads, so a violation confined
                // to one genome of a well-represented type would show up only by chance.
                //
                // Checking it here and not before the dedup keeps this free on the common path - the
                // dedup exists to avoid exactly this work - and still reaches every occurrence the
                // recall does not see. Counted separately and never folded into tp/tp+fn: those are
                // the paper's numbers and this must not move them.
                probeDuplicate++;
                if (!tree.isAncestorOf(nodeByPos[leafPos], storedNode)) {
                    readerPathViolations++;
                }
                return false;
            }
            entries++;
            // Path correctness says the genomes carrying a k-mer all lie below the node it is stored
            // at, so this holds for every genome of the type and it does not matter which of them won
            // the filter. Where it fails, the database claims the k-mer for a branch this genome is
            // not on - the false negative the recall is there to count. The tree answers from the two
            // depths, so an unrelated node costs no walk at all.
            boolean onPath = tree.isAncestorOf(nodeByPos[leafPos], storedNode);
            readerTypeTpPlusFn[stIndex]++;
            if (onPath) {
                readerTypeTp[stIndex]++;
                readerTpForNode[storedPos]++;
                return true;
            }
            return false;
        }
    }
}
