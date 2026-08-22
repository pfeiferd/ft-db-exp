package org.metagene.ftdbexp.eval;

import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.FinerTreeMaker;
import org.metagene.genestrip.goals.MatchResultGoal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.FastqKMerMatcher;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

/**
 * Runs a Genestrip read-matching goal over a set of simulated fastq files and tallies how well the
 * resulting classifications agree with the reads' ground truth.
 * <p>
 * The same evaluator serves the unrefined and the refined database: both are reached through the
 * same {@link MatchResultGoal}, only under a different {@link GoalKey}. The NCBI taxonomy and the
 * accession map are loaded once and reused across runs, whereas the candidate species counts are
 * derived per run from the taxonomy of the database variant at hand -- that tree is precisely what
 * the refinement changes.
 */
public class AccuracyEvaluator {
    private final File baseDir;
    private final TaxTree taxTree;
    private final Simulator simulator;
    private final Map<String, String> extractedTaxIds;
    /**
     * The resolver for a read's true taxon. Built in {@link #evaluate} rather than in the
     * constructor, and only for a run that has a ground truth at all - see the note there.
     */
    private ReadGroundTruth groundTruth;

    /**
     * Creates the evaluator, loading the taxonomy and the accession map of the given database.
     *
     * @param baseDir the Genestrip base directory, i.e. the one holding {@code common} and
     *                {@code projects}
     * @param db        the name of the database project the reads are matched against
     * @param simulator the simulator that produced the reads, which determines how their ground
     *                  truth is read from the read identifiers
     * @throws IOException if the taxonomy or the accession map cannot be read
     */
    public AccuracyEvaluator(File baseDir, String db, Simulator simulator) throws IOException {
        this.baseDir = baseDir;

        FTProject project = newProject(db, null);
        FinerTreeMaker<FTProject> maker = new FinerTreeMaker<FTProject>(project);
        try {
            @SuppressWarnings("unchecked")
            ObjectGoal<TaxTree, FTProject> taxTreeGoal =
                    (ObjectGoal<TaxTree, FTProject>) maker.getGoal(GSGoalKey.TAXTREE);
            taxTree = taxTreeGoal.get();

            // The table written next to the extracted genomes is the ground truth for reads
            // simulated from them, and the only source that covers the ones taken from Genbank. It
            // is loaded here but not turned into a resolver yet: see evaluate().
            this.simulator = simulator;
            this.extractedTaxIds = ExtractedTaxIds.load(
                    project.getOutputFile(GSGoalKey.EXTRACT_REFSEQ_CSV.getName(), GSProject.GSFileType.CSV, false));
        } finally {
            maker.dumpAll();
        }
    }

    /**
     * Matches the fastq files listed in the given mapping file against one variant of the database
     * and tallies the classification accuracy per fastq key.
     *
     * @param db            the name of the database project
     * @param fqMapFile     the fastq mapping file, relative to the project's {@code txt} directory,
     *                      with one {@code <key> <path>} line per fastq file
     * @param matchGoalKey  the matching goal to run, i.e. the database variant to evaluate
     * @param loadDbGoalKey the goal loading that variant's database, whose taxonomy supplies the
     *                      candidate species counts
     * @param scope         restricts the reads that count towards recall to those whose true taxon
     *                      lies at or below a requested node of this tree; pass {@code null} to use
     *                      the database's own taxonomy, i.e. to count exactly the reads the database
     *                      was built to cover
     * @param baseline        records which reads were left at their genus, so that the gain can be
     *                        related to exactly those; may be {@code null} to skip that measurement
     * @param obsBaseline     records the observable substitute for that subset -- the reads assigned
     *                        to a node at genus rank, which needs no ground truth
     * @param collectBaseline whether this run fills the baselines (the unrefined one) or consults them
     * @param groundTruthFree whether the reads have no known ground truth, in which case only the
     *                        ungated measures are tallied and no accession map is consulted
     * @return the tallies keyed by fastq key, in the order the files were processed
     * @throws IOException if the database or the fastq files cannot be read
     */
    public Map<String, AccuracyTally> evaluate(String db, String fqMapFile, GoalKey matchGoalKey,
                                               GoalKey loadDbGoalKey, SmallTaxTree scope,
                                               GenusOnlyBaseline baseline, GenusOnlyBaseline obsBaseline,
                                               boolean collectBaseline, boolean groundTruthFree)
            throws IOException {
        // Built here rather than in the constructor, because only here is it known whether the run
        // has a ground truth at all. IssReadGroundTruth refuses to exist without an extraction
        // table -- rightly, since every read would else go unresolved while looking like a result --
        // and SpecificityReport passes ISS purely as a placeholder for a ground-truth-free run that
        // never resolves anything. Constructing eagerly therefore made that run fail on a file it
        // does not read:
        //   IllegalArgumentException: No extracted genomes to resolve the ground truth against.
        // The assignment happens before the matcher starts its worker threads, so the callback
        // below sees it without further synchronisation; a lazy initialisation inside the callback
        // would be a data race, since it runs on every consumer thread at once.
        if (!groundTruthFree && groundTruth == null) {
            groundTruth = simulator.groundTruth(taxTree, extractedTaxIds);
        }
        FTProject project = newProject(db, fqMapFile);
        FinerTreeMaker<FTProject> maker = new FinerTreeMaker<FTProject>(project);
        Map<String, AccuracyTally> result = new LinkedHashMap<String, AccuracyTally>();
        // One tally per matcher thread instead of one shared tally behind a lock. The callback below
        // runs for every single read, and guarding it globally serialised the whole evaluation: the
        // matcher's worker threads spent their time queueing for that monitor rather than matching,
        // so the JVM sat at roughly one busy core no matter how many threads were configured. Each
        // thread now accumulates on its own and the results are summed once per fastq file, which is
        // exact -- every counter is a sum or a count, so the merged tally equals what a single
        // thread would have produced.
        final Queue<AccuracyTally> threadTallies = new ConcurrentLinkedQueue<AccuracyTally>();
        final ThreadLocal<AccuracyTally> localTally = new ThreadLocal<AccuracyTally>() {
            @Override
            protected AccuracyTally initialValue() {
                AccuracyTally fresh = new AccuracyTally();
                threadTallies.add(fresh);
                return fresh;
            }
        };

        try {
            @SuppressWarnings("unchecked")
            ObjectGoal<Database, FTProject> dbGoal =
                    (ObjectGoal<Database, FTProject>) maker.getGoal(loadDbGoalKey);
            SmallTaxTree dbTree = dbGoal.get().getTaxTree();
            SpeciesCandidates candidates = new SpeciesCandidates();
            // Without an explicit scope the database's own taxonomy is the right one. A project that
            // requests only some genera of a RefSeq category - as the protozoa ones do - is fed reads
            // from the whole category, and a read whose organism the database never covered must not
            // count against its recall. Where a project requests an entire category, as the viral one
            // does, every read is in scope anyway and this changes nothing.
            SmallTaxTree effectiveScope = scope != null ? scope : dbTree;

            MatchResultGoal<?> matchResGoal = (MatchResultGoal<?>) maker.getGoal(matchGoalKey);
            matchResGoal.setAfterMatchCallback(new MatchResultGoal.AfterMatchCallback() {
                @Override
                public void afterMatch(FastqKMerMatcher.MatcherReadEntry entry, boolean found) {
                    // Deliberately unsynchronised. The tally belongs to this thread alone, and the
                    // two structures shared with the other threads -- the candidate counter and the
                    // baseline sets -- are concurrent by construction.
                    //
                    // That the merge in afterKey() sees these writes is guaranteed by Genestrip
                    // rather than by a lock here: a consumer thread of AbstractFastqReader calls
                    // this callback from nextEntry() and only then writes the volatile
                    // ReadEntry.pooled, and the producer polls every entry of the pool for that flag
                    // before it finishes the file. The volatile write and the matching read
                    // establish a happens-before edge covering everything the consumer did first,
                    // and afterKey() runs on that same producer thread. Should the reader ever stop
                    // handing entries back through `pooled', this reasoning has to be redone.
                    record(localTally.get(), entry, dbTree, candidates, effectiveScope, baseline,
                            obsBaseline, collectBaseline, groundTruthFree);
                }

                @Override
                public void afterKey(String key, MatchingResult res) {
                    // Called once the file is done, i.e. after the worker threads have finished with
                    // it, so the per-thread tallies are complete and can be summed.
                    AccuracyTally merged = new AccuracyTally();
                    for (AccuracyTally threadTally : threadTallies) {
                        merged.add(threadTally);
                        threadTally.reset();
                    }
                    // Both subsets are advanced in lockstep: each is filled by the unrefined run
                    // and consulted by the refined one, so either variant is scored on the reads the
                    // unrefined run put in that subset and not on its own.
                    if (baseline != null) {
                        if (collectBaseline) {
                            baseline.endCollecting(key);
                            obsBaseline.endCollecting(key);
                        } else {
                            baseline.endConsulting(key);
                            obsBaseline.endConsulting(key);
                        }
                    }
                    System.out.println(key + ": " + merged);
                    warnIfUnresolved(key, merged);
                    result.put(key, merged);
                }
            });
            matchResGoal.make();
        } finally {
            maker.dumpAll();
        }
        return result;
    }

    /**
     * Reports how many reads of a fastq file had no resolvable ground truth.
     * <p>
     * Such reads are excluded from every count, so a substantial number means the figures rest on
     * only part of the data. The count also lands in the result CSV, but a run is long enough that
     * nobody should have to go looking for it afterwards.
     *
     * @param key   the fastq key just finished
     * @param tally its counts
     */
    private static void warnIfUnresolved(String key, AccuracyTally tally) {
        long unresolved = tally.getUnresolved();
        if (unresolved == 0) {
            return;
        }
        long considered = tally.getTotal() + unresolved;
        System.err.printf("WARNING: %s: %,d of %,d reads (%.1f %%) have no resolvable ground truth"
                        + " and are excluded from every measure.%n",
                key, unresolved, considered, considered == 0 ? 0.0 : 100.0 * unresolved / considered);
    }

    /**
     * Adds one read to the tally by comparing its assigned taxon against its ground truth.
     *
     * @param tally      the tally to update
     * @param entry      the matcher's result for the read
     * @param dbTree     the taxonomy of the database variant being evaluated
     * @param candidates the candidate species counter for that taxonomy
     * @param scope      the scope restricting which reads count towards recall, may be {@code null}
     */
    private void record(AccuracyTally tally, FastqKMerMatcher.MatcherReadEntry entry, SmallTaxTree dbTree,
                        SpeciesCandidates candidates, SmallTaxTree scope, GenusOnlyBaseline baseline,
                        GenusOnlyBaseline obsBaseline, boolean collectBaseline, boolean groundTruthFree) {
        if (groundTruthFree) {
            recordWithoutGroundTruth(tally, entry, candidates, obsBaseline, collectBaseline);
            return;
        }
        TaxTree.TaxIdNode trueNode = groundTruth.resolve(entry.readDescriptor, entry.readDescriptorSize);
        if (trueNode == null) {
            tally.recordUnresolved();
            return;
        }
        // Everything below is compared inside the database's own taxonomy, never against the NCBI
        // taxonomy: a database contains synthetic nodes - the data nodes and, after a refinement,
        // the refined ones - whose tax ids do not exist in NCBI at all. Looking those up in the NCBI
        // tree yields nothing and would silently drop every read classified to one of them.
        SmallTaxTree.SmallTaxIdNode classNode = entry.classNode;
        if (inScope(trueNode, scope)) {
            SmallTaxTree.SmallTaxIdNode trueInDb = inDbTree(trueNode, dbTree);
            Rank lcaRank = classNode == null || trueInDb == null
                    ? null : lowestRankedCommonAncestor(dbTree, trueInDb, classNode);
            // The ungated score drops the test that the read's true species is still in question at
            // the assigned node. It therefore states how far the classification narrows the species
            // down rather than how far it narrows them down *correctly* -- and, since it never looks
            // at the ground truth, it is the one score obtainable from a real fastq file.
            double ungatedScore = classNode != null ? candidates.weightFor(classNode) : 0;
            double score = classNode != null && trueInDb != null
                    && SpeciesCandidates.areComparable(classNode, trueInDb)
                    ? ungatedScore : 0;
            boolean genusOnly = false;
            boolean obsGenusOnly = false;
            if (baseline != null) {
                String descriptor = new String(entry.readDescriptor, 0, entry.readDescriptorSize,
                        StandardCharsets.UTF_8);
                if (collectBaseline) {
                    // R_g: correct down to the genus but no further, the refinement's only
                    // opportunity. Governed by the positives, i.e. by sigma(r), because prec_g is a
                    // precision and a read placed in the wrong genus is not one a refinement can
                    // put right.
                    genusOnly = isGenusOnlyRank(lcaRank);
                    if (genusOnly) {
                        baseline.collect(descriptor);
                    }
                    // R_g': the same rank window read off the assignment alone. It is not a
                    // subset of R_g and does not contain it either -- it drops the reads whose
                    // assignment named a sibling species (their lca with sigma(r) is the genus, but
                    // the node itself sits at a species) and admits those left at a foreign genus.
                    obsGenusOnly = isGenusOnlyNode(classNode);
                    if (obsGenusOnly) {
                        obsBaseline.collect(descriptor);
                    }
                } else {
                    genusOnly = baseline.contains(descriptor);
                    obsGenusOnly = obsBaseline.contains(descriptor);
                }
            }
            // The subset is fixed by the unrefined run; on the refined pass the flag is looked up
            // rather than recomputed, so both variants are scored on identical reads.
            tally.record(classNode != null, lcaRank, score, ungatedScore, genusOnly, obsGenusOnly);
        } else if (classNode != null && inScope(taxTree.getNodeByTaxId(classNode.getTaxId()), scope)) {
            // The read does not belong to the scope but was classified into it: a false positive.
            tally.recordOutOfScopeClassification();
        }
    }

    /**
     * Adds one read of a fastq file without known ground truth to the tally.
     * <p>
     * Nothing here consults {@code groundTruth}, the scope or the read's true taxon, because none of
     * them exists for a real sample. What remains observable is the assigned node and, through the
     * database's taxonomy, how many species it leaves in question -- which is the ungated score of
     * Section \"Estimating the gain without ground truth\" in the paper.
     *
     * @param tally           the tally to update
     * @param entry           the matcher's result for the read
     * @param candidates      the candidate species counter for the variant's taxonomy
     * @param obsBaseline     the observable genus-only subset, filled by the unrefined run
     * @param collectBaseline whether this run fills that subset or consults it
     */
    private void recordWithoutGroundTruth(AccuracyTally tally, FastqKMerMatcher.MatcherReadEntry entry,
                                          SpeciesCandidates candidates, GenusOnlyBaseline obsBaseline,
                                          boolean collectBaseline) {
        SmallTaxTree.SmallTaxIdNode classNode = entry.classNode;
        double ungatedScore = classNode != null ? candidates.weightFor(classNode) : 0;
        boolean obsGenusOnly;
        String descriptor = new String(entry.readDescriptor, 0, entry.readDescriptorSize,
                StandardCharsets.UTF_8);
        if (collectBaseline) {
            obsGenusOnly = isGenusOnlyNode(classNode);
            if (obsGenusOnly) {
                obsBaseline.collect(descriptor);
            }
        } else {
            obsGenusOnly = obsBaseline.contains(descriptor);
        }
        tally.recordWithoutGroundTruth(classNode != null, ungatedScore, obsGenusOnly);
    }

    /**
     * Maps a taxon of the NCBI taxonomy onto the database's taxonomy, walking upwards until a node
     * is found that the database actually contains. Databases keep an extract of the taxonomy, so a
     * read's exact taxon may be absent even though its species or genus is present.
     *
     * @param node   the taxon to map
     * @param dbTree the taxonomy of the database variant
     * @return the corresponding node of the database's taxonomy, or {@code null} if none of the
     * taxon's ancestors is contained either
     */
    private SmallTaxTree.SmallTaxIdNode inDbTree(TaxTree.TaxIdNode node, SmallTaxTree dbTree) {
        for (TaxTree.TaxIdNode current = node; current != null; current = current.getParent()) {
            SmallTaxTree.SmallTaxIdNode dbNode = dbTree.getNodeByTaxId(current.getTaxId());
            if (dbNode != null) {
                return dbNode;
            }
        }
        return null;
    }

    /**
     * Returns the rank at which the two taxa agree, i.e. the rank of their lowest common ancestor.
     * Ancestors without a rank of their own are skipped upwards, since they carry no information
     * about how specific the agreement is.
     *
     * @param dbTree    the taxonomy of the database variant being evaluated
     * @param trueNode  the taxon the read was generated from, mapped into that taxonomy
     * @param classNode the taxon the read was classified to
     * @return the rank the two taxa agree at, or {@code null} if their common ancestor has no rank
     * up to the root
     */
    private Rank lowestRankedCommonAncestor(SmallTaxTree dbTree, SmallTaxTree.SmallTaxIdNode trueNode,
                                            SmallTaxTree.SmallTaxIdNode classNode) {
        SmallTaxTree.SmallTaxIdNode lca = dbTree.getLowestCommonAncestor(trueNode, classNode);
        while (lca != null && !isTaxonomicRank(lca.getRank())) {
            lca = lca.getParent();
        }
        return lca == null ? null : lca.getRank();
    }

    /**
     * Returns whether the rank two taxa agree at leaves the read in the window a refinement can act
     * on: at or below a genus, but still above the species. This is the membership test of the
     * genus-only subset $R_g$, and it is applied to the rank of the lowest common ancestor of the
     * read's true taxon and its assignment, so the subset is governed by the ground truth.
     * <p>
     * The window is a range of ranks rather than the genus alone because a taxonomy places further
     * ranks inside it -- see {@link #isGenusOnlyNode} for why that matters and by how much. Unranked
     * ancestors need no handling here: {@link #lowestRankedCommonAncestor} has already resolved them
     * upwards, so the two subsets are formed by one and the same rule.
     *
     * @param lcaRank the rank the read's true taxon and its assignment agree at, may be
     *                {@code null} if they share no ranked ancestor at all
     * @return whether the agreement stops inside the genus-only window
     */
    private static boolean isGenusOnlyRank(Rank lcaRank) {
        return isAtLeast(lcaRank, Rank.GENUS) && !isAtLeast(lcaRank, Rank.SPECIES);
    }

    /**
     * Returns whether the node a read was assigned to leaves it in the window a refinement can act
     * on: at or below a genus, but still above the species.
     * <p>
     * A node carrying no taxonomic rank of its own is resolved to its nearest ranked ancestor first,
     * exactly as {@link #lowestRankedCommonAncestor} does for the genus-only subset, so that the two
     * subsets are formed by the same rule. That is what makes unranked nodes count: a read placed on
     * an unranked node beneath a genus is a genus-level answer and a refinement moves it, while one
     * placed beneath a species resolves to that species and is already as specific as it can be.
     * <p>
     * Testing the node's own rank against the genus alone would be too narrow twice over.
     * {@link Rank} places {@link Rank#SUBGENUS} and {@link Rank#SPECIES_GROUP} between genus and
     * species, and a read left at either has several species in question just as one left at the
     * genus has; and {@link Rank#NO_RANK} nodes abound -- 13,458 of them in the viral database
     * alone. Excluding the ranked ones cost 13 % of the improvable reads of the tick-borne database,
     * whose six species-group nodes carry the Rickettsia spotted fever, typhus, canis, belli and
     * phagocytophilum groups, which is precisely where close relatives accumulate.
     * <p>
     * The artificial ranks Genestrip marks its own nodes with are not taxonomic ranks either, so a
     * data node resolves to the species above it and drops out, which is right. A refined node would
     * resolve to the genus above it and count -- but cannot occur here, because this subset is only
     * ever determined from the unrefined run.
     *
     * @param node the node a read was assigned to, may be {@code null} for an unclassified read
     * @return whether the read was placed no further than a genus
     */
    private static boolean isGenusOnlyNode(SmallTaxTree.SmallTaxIdNode node) {
        SmallTaxTree.SmallTaxIdNode ranked = node;
        while (ranked != null && !isTaxonomicRank(ranked.getRank())) {
            ranked = ranked.getParent();
        }
        if (ranked == null) {
            return false;
        }
        Rank rank = ranked.getRank();
        return isAtLeast(rank, Rank.GENUS) && !isAtLeast(rank, Rank.SPECIES);
    }

    /**
     * Returns whether a rank is at or below the given one.
     *
     * @param rank      the rank to test, may be {@code null}
     * @param threshold the rank to compare against
     * @return whether {@code rank} is as specific as {@code threshold} or more so
     */
    private static boolean isAtLeast(Rank rank, Rank threshold) {
        return rank != null && (threshold.equals(rank) || rank.isBelow(threshold));
    }

    /**
     * Returns whether a rank says something about the taxonomic specificity of a node.
     * <p>
     * Besides the nodes without a rank of their own, this excludes the artificial ranks Genestrip
     * uses to mark where a k-mer originates from. Those are ordered at the very bottom of
     * {@link Rank}, so {@code Rank.REFINED.isBelow(Rank.SPECIES)} holds -- which is true of a data
     * node, sitting underneath its species, but decidedly false of a refined node: a refined node
     * lies <em>between</em> a genus and its species and groups several of them. Comparing ranks
     * naively would therefore count a read assigned to a refined node as a species-level hit even
     * though that assignment is precisely what leaves the species open. Skipping the artificial
     * ranks upwards resolves both cases correctly: a data node yields the species above it, a
     * refined node the genus.
     *
     * @param rank the rank to test, may be {@code null}
     * @return whether the rank is a genuine taxonomic rank
     */
    private static boolean isTaxonomicRank(Rank rank) {
        return rank != null
                && !Rank.NO_RANK.equals(rank)
                && !Rank.REFINED.equals(rank)
                && !Rank.DATA.equals(rank)
                && !Rank.FILE.equals(rank)
                && !Rank.ID.equals(rank);
    }

    /**
     * Returns whether a taxon lies at or below a node that the given scope explicitly requests.
     *
     * @param node  the taxon to test
     * @param scope the scope, or {@code null} for no restriction
     * @return whether the taxon is in scope
     */
    private boolean inScope(TaxTree.TaxIdNode node, SmallTaxTree scope) {
        if (scope == null) {
            return true;
        }
        if (node == null) {
            return false;
        }
        SmallTaxTree.SmallTaxIdNode scopeNode = scope.getNodeByTaxId(node.getTaxId());
        while (scopeNode != null) {
            if (scopeNode.isRequested()) {
                return true;
            }
            scopeNode = scopeNode.getParent();
        }
        return false;
    }

    /**
     * Creates a project for the given database, optionally bound to a fastq mapping file.
     *
     * @param db        the name of the database project
     * @param fqMapFile the fastq mapping file, or {@code null} if no fastq files are to be processed
     * @return the project, configured to use all available threads
     * @throws IOException if the project's configuration cannot be read
     */
    private FTProject newProject(String db, String fqMapFile) throws IOException {
        FTProject project = new FTProject(new GSCommon(baseDir), db, null, null, fqMapFile,
                null, null, null, null, null, null, false);
        project.initConfigParam(GSConfigKey.THREADS, -1);
        // Downloaded once into data/fastq rather than streamed anew on every run. Without this the
        // fastq map keeps its URLs: FastqMapTransformGoal passes them through untouched, whereby
        // FastqDownloadsGoal.getFiles() finds no file to fetch and does nothing, and every pass over
        // a URL-based map re-reads it over the network. `data/fastq' and not the project's own folder,
        // so that the file lands where `mvn exec:exec@fastqdl' (which passes -ll) puts it and where
        // run_classification_exps.sh looks for it -- named by the map key, e.g. data/fastq/P1.fastq.gz.
        // Maps naming local files are unaffected: the transform only rewrites URL resources.
        project.setDownloadFastqsToCommon(true);
        return project;
    }
}
