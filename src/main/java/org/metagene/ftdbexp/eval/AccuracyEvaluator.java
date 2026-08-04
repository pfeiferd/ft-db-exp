package org.metagene.ftdbexp.eval;

import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.FinerTreeMaker;
import org.metagene.genestrip.goals.MatchResultGoal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.FastqKMerMatcher;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;

import java.io.File;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;

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
    private final ReadGroundTruth groundTruth;

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

            @SuppressWarnings("unchecked")
            ObjectGoal<AccessionMap, FTProject> accessionMapGoal =
                    (ObjectGoal<AccessionMap, FTProject>) maker.getGoal(GSGoalKey.ACCMAP);
            groundTruth = simulator.groundTruth(accessionMapGoal.get());
            accessionMapGoal.cleanThis();
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
     *                      lies at or below a requested node of this tree; may be {@code null} to
     *                      count every read whose ground truth can be resolved
     * @return the tallies keyed by fastq key, in the order the files were processed
     * @throws IOException if the database or the fastq files cannot be read
     */
    public Map<String, AccuracyTally> evaluate(String db, String fqMapFile, GoalKey matchGoalKey,
                                               GoalKey loadDbGoalKey, SmallTaxTree scope) throws IOException {
        FTProject project = newProject(db, fqMapFile);
        FinerTreeMaker<FTProject> maker = new FinerTreeMaker<FTProject>(project);
        Map<String, AccuracyTally> result = new LinkedHashMap<String, AccuracyTally>();
        AccuracyTally tally = new AccuracyTally();

        try {
            @SuppressWarnings("unchecked")
            ObjectGoal<Database, FTProject> dbGoal =
                    (ObjectGoal<Database, FTProject>) maker.getGoal(loadDbGoalKey);
            SmallTaxTree dbTree = dbGoal.get().getTaxTree();
            SpeciesCandidates candidates = new SpeciesCandidates();

            MatchResultGoal<?> matchResGoal = (MatchResultGoal<?>) maker.getGoal(matchGoalKey);
            matchResGoal.setAfterMatchCallback(new MatchResultGoal.AfterMatchCallback() {
                @Override
                public void afterMatch(FastqKMerMatcher.MatcherReadEntry entry, boolean found) {
                    // The matcher calls this from several threads, so the tally needs guarding.
                    synchronized (tally) {
                        record(tally, entry, dbTree, candidates, scope);
                    }
                }

                @Override
                public void afterKey(String key, MatchingResult res) {
                    synchronized (tally) {
                        System.out.println(key + ": " + tally);
                        result.put(key, tally.copy());
                        tally.reset();
                    }
                }
            });
            matchResGoal.make();
        } finally {
            maker.dumpAll();
        }
        return result;
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
                        SpeciesCandidates candidates, SmallTaxTree scope) {
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
            double score = classNode != null && trueInDb != null
                    && SpeciesCandidates.areComparable(classNode, trueInDb)
                    ? candidates.weightFor(classNode) : 0;
            tally.record(classNode != null, lcaRank, score);
        } else if (classNode != null && inScope(taxTree.getNodeByTaxId(classNode.getTaxId()), scope)) {
            // The read does not belong to the scope but was classified into it: a false positive.
            tally.recordOutOfScopeClassification();
        }
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
        while (lca != null && Rank.NO_RANK.equals(lca.getRank())) {
            lca = lca.getParent();
        }
        return lca == null ? null : lca.getRank();
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
        return project;
    }
}
