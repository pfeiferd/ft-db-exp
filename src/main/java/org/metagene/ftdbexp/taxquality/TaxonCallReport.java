package org.metagene.ftdbexp.taxquality;

import org.metagene.ftdbexp.eval.RefinementAccuracyReport.Variant;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.FinerTreeMaker;
import org.metagene.genestrip.goals.MatchResultGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.CountsPerTaxid;
import org.metagene.genestrip.match.FastqKMerMatcher;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;

/**
 * NOT REPORTED IN THE PAPER, and kept for the reason the naive Bayes rule beside it is kept: it was
 * needed to find out that it should not be.
 * <p>
 * A sputum sample carries a community of streptococci rather than one organism, and the reference
 * standard of the case study names several organisms for nine of the thirteen samples it can decide.
 * A rule that crowns a single taxon per sample therefore answers a question the material does not
 * pose, and what it answers is decided by the rule as much as by the database. The paper reports the
 * per-read measure of Section "Estimating the gain without ground truth" instead, which needs no
 * such rule and rests on millions of reads rather than on thirteen samples.
 * <p>
 * What this did establish, and what the case study now states, came out of running it: that the
 * unnamed `Streptococcus sp.' genomes of the database hold the pneumococcus's k-mers hostage, and
 * that a candidate set counted in strains rather than species scores a correct call as a miss.
 */
public class TaxonCallReport {
    private final File baseDir;
    private final File resultsDir;

    /** The currencies a sample's contribution is counted in. */
    /**
     * What the vote is taken over. The three are not variants of one measurement but three questions,
     * and a sample can answer them differently.
     * <p>
     * {@link #READS} counts the reads Genestrip classified to a taxon, so each read carries one vote
     * and a long read counts no more than a short one. {@link #KMERS} counts every matched k-mer
     * specific to a taxon's genome, whether or not it sits in a read that was classified there -- so
     * a read whose k-mers are split across several species still contributes all of them, and the
     * evidence of the sample is weighed rather than its reads counted. {@link #UNIQUE} is the same
     * over distinct k-mers, which keeps a single deeply covered region from outvoting the rest.
     * <p>
     * Which one to believe is not decided here. The report writes a row per currency and leaves the
     * comparison to whoever reads it, since the currencies can disagree and their disagreement is
     * itself worth seeing: a call that holds in all three rests on something other than coverage.
     */
    /** How many of the heaviest votes are reported per sample, variant and currency. */
    private static final int TOP_VOTES = 3;

    /**
     * Where the second-coarsest bucket of the summary stops.
     * <p>
     * Above this a node leaves so much of the tree open that reporting it apart from the root would
     * suggest a narrowing that did not happen: for {@code strepto} the root holds 157 species, and a
     * node still leaving 120 has narrowed nothing worth a column of its own.
     */
    private static final int COARSE_MAX = 100;

    /** The buckets of the summary, coarsest last. The ranges are on the candidate count. */
    private static final String[] BUCKETS = {"species", "2-10", "11-50", "51-100", "above 100"};

    private enum Currency {
        READS("reads"), KMERS("kmers"), UNIQUE("unique kmers");

        private final String label;

        Currency(String label) {
            this.label = label;
        }

        long of(CountsPerTaxid c) {
            switch (this) {
                case READS:  return c.getReads();
                case KMERS:  return c.getKMers();
                default:     return c.getUniqueKMers();
            }
        }
    }

    /**
     * @param baseDir    the Genestrip base directory holding {@code common} and {@code projects}
     * @param resultsDir the directory the CSV is written to
     */
    public TaxonCallReport(File baseDir, File resultsDir) {
        this.baseDir = baseDir;
        this.resultsDir = resultsDir;
    }

    /**
     * Runs both database variants over the fastq map and writes one row per sample, variant and
     * currency to {@code <db>_<report key>_taxoncall.csv}.
     *
     * @param db         the name of the database project
     * @param fqMapFile  the fastq mapping file, resolved as usual against {@code data/fastq}
     * @param reportKey  short name used in the result file name
     * @param rootTaxId  the node the vote is restricted to, e.g. {@code 1301} for Streptococcus
     * @param targetTaxId the reference organism the nearest-node columns are measured towards, e.g.
     *                    {@code 1313} for S. pneumoniae; empty or null leaves those columns blank
     * @param minimum    how much the winning path must gather, and how much a node must hold before
     *                   the nearest-node columns will stop at it, or
     *                   {@link PathVoteTaxonModel#NO_MINIMUM}
     * @return the file that was written
     * @throws IOException if a database or a fastq file cannot be read, or the file cannot be written
     */
    public File write(String db, String fqMapFile, String reportKey, String rootTaxId,
                      String targetTaxId, long minimum) throws IOException {
        // Read before the runs, not after: each run scores its calls against it as they are made.
        Map<String, List<String[]>> truth = GroundTruth.read(new File(baseDir, "projects/" + db + "/ground_truth.csv"));
        Map<Variant, Map<String, Map<Currency, List<Call>>>> byVariant = new LinkedHashMap<>();
        for (Variant variant : Variant.values()) {
            System.out.println("=== " + db + " / " + variant.getLabel() + " ===");
            byVariant.put(variant, run(db, fqMapFile, variant, rootTaxId, targetTaxId, minimum, truth));
        }

        File file = new File(resultsDir, db + "_" + reportKey + "_taxoncall.csv");
        try (PrintStream ps = new PrintStream(new FileOutputStream(file), false, StandardCharsets.UTF_8.name())) {
            // One row per sample, variant and currency; the votes are side by side in it rather than
            // under one another, so that a reader compares the three at a glance and a spreadsheet
            // sorts on any of them without regrouping. Everything describing the sample stands once,
            // at the front; a vote that does not exist leaves its block empty.
            StringBuilder header = new StringBuilder("db;sample;variant;currency;contributing nodes;total;"
                    + "unplaced;culture;wimp;culture taxid;wimp taxid;"
                    + "nearest node;nearest node name;nearest node rank;nearest candidates;"
                    + "nearest at;nearest at or below;nearest share at;nearest share at or below;");
            for (int i = 1; i <= TOP_VOTES; i++) {
                header.append("node ").append(i).append(";node name ").append(i)
                        .append(";node rank ").append(i).append(";taxon ").append(i)
                        .append(";candidates ").append(i).append(";q culture ").append(i)
                        .append(";q wimp ").append(i).append(';');
            }
            ps.println(header);
            Set<String> samples = new LinkedHashSet<>();
            for (Map<String, Map<Currency, List<Call>>> m : byVariant.values()) {
                samples.addAll(m.keySet());
            }
            for (String sample : samples) {
                for (Variant variant : Variant.values()) {
                    Map<Currency, List<Call>> perCurrency = byVariant.get(variant).get(sample);
                    if (perCurrency == null) {
                        continue;
                    }
                    for (Currency currency : Currency.values()) {
                        List<Call> calls = perCurrency.get(currency);
                        if (calls == null || calls.isEmpty()) {
                            continue;
                        }
                        // The sample's own figures, taken from the first vote because every vote
                        // carries the same ones.
                        Call first = calls.get(0);
                        ps.print(db); ps.print(';');
                        ps.print(sample); ps.print(';');
                        ps.print(variant.getLabel()); ps.print(';');
                        ps.print(currency.label); ps.print(';');
                        ps.print(first.contributingNodes); ps.print(';');
                        ps.print(first.total); ps.print(';');
                        ps.print(first.unplaced); ps.print(';');
                        ps.print(join(truth.get(sample), 0)); ps.print(';');
                        ps.print(join(truth.get(sample), 1)); ps.print(';');
                        ps.print(first.sigmaCulture == null ? "" : first.sigmaCulture); ps.print(';');
                        ps.print(first.sigmaWimp == null ? "" : first.sigmaWimp); ps.print(';');
                        // Empty and not zero where no signal reaches the organism's lineage: a
                        // sample that put nothing there did not put none there, it was never asked.
                        boolean hasNear = first.nearTaxId != null;
                        ps.print(hasNear ? first.nearTaxId : ""); ps.print(';');
                        ps.print(hasNear && first.nearName != null ? first.nearName : ""); ps.print(';');
                        ps.print(hasNear && first.nearRank != null ? first.nearRank : ""); ps.print(';');
                        ps.print(hasNear ? Integer.toString(first.nearCandidates) : ""); ps.print(';');
                        ps.print(hasNear ? Long.toString(first.nearAt) : ""); ps.print(';');
                        ps.print(hasNear ? Long.toString(first.nearBelow) : ""); ps.print(';');
                        // Both shares are of the sample's total below the root, so that a row says
                        // what part of the sample the answer rests on and not merely where it landed.
                        ps.print(hasNear && first.total > 0
                                ? format((double) first.nearAt / first.total) : ""); ps.print(';');
                        ps.print(hasNear && first.total > 0
                                ? format((double) first.nearBelow / first.total) : ""); ps.print(';');
                        for (int i = 0; i < TOP_VOTES; i++) {
                            Call c = i < calls.size() ? calls.get(i) : null;
                            ps.print(c == null || c.nodeTaxId == null ? "" : c.nodeTaxId); ps.print(';');
                            ps.print(c == null || c.nodeName == null ? "" : c.nodeName); ps.print(';');
                            ps.print(c == null || c.nodeRank == null ? "" : c.nodeRank); ps.print(';');
                            ps.print(c == null || c.taxon == null ? "" : c.taxon); ps.print(';');
                            ps.print(c == null || c.nodeTaxId == null ? "" : Integer.toString(c.candidates));
                            ps.print(';');
                            ps.print(c == null ? "" : format(c.qCulture)); ps.print(';');
                            ps.print(c == null ? "" : format(c.qWimp)); ps.print(';');
                        }
                        ps.println();
                    }
                }
            }
        }
        System.out.println("Wrote " + file);
        writeSummary(db, reportKey, rootTaxId, targetTaxId, byVariant, truth);
        return file;
    }

    /**
     * Writes the per-currency, per-variant summary of the nearest-node columns, one row for each
     * pair, to {@code <db>_<report key>_taxoncallsummary.csv}.
     * <p>
     * It exists so that a table can be typeset from a file rather than from figures somebody counted:
     * the bucketing, the medians and the sums are made here, next to the data they describe and in
     * the same pass that produced it, and the paper reads the result. A count that is wrong is then
     * wrong in one place.
     * <p>
     * A bucket no sample reached carries an empty median rather than a zero: there is nothing to take
     * a median of, and a zero would read as evidence of nothing where there is no evidence either
     * way. The coarsest bucket takes the root of the walk whatever its candidate count, so that it
     * stays right when that count moves.
     */
    private void writeSummary(String db, String reportKey, String rootTaxId, String targetTaxId,
                              Map<Variant, Map<String, Map<Currency, List<Call>>>> byVariant,
                              Map<String, List<String[]>> truth) throws IOException {
        File file = new File(resultsDir, db + "_" + reportKey + "_taxoncallsummary.csv");
        Set<String> samples = new LinkedHashSet<>();
        for (Map<String, Map<Currency, List<Call>>> m : byVariant.values()) {
            samples.addAll(m.keySet());
        }
        try (PrintStream ps = new PrintStream(new FileOutputStream(file), false, StandardCharsets.UTF_8.name())) {
            StringBuilder header = new StringBuilder("db;currency;variant;samples;target;root species;");
            for (String b : BUCKETS) {
                header.append(b).append(';').append(b).append(" median;");
            }
            header.append("off lineage;median in-genus total;median at or below species;"
                    + "max ladder candidates;unplaced;contributing nodes;top calls unchanged;"
                    + "off lineage named target;median share at;median share at or below;"
                    + "names target;names target list;misses target;misses target list;");
            ps.println(header);
            for (Currency currency : Currency.values()) {
                for (Variant variant : Variant.values()) {
                    Map<String, Map<Currency, List<Call>>> perSample = byVariant.get(variant);
                    if (perSample == null) {
                        continue;
                    }
                    int[] counts = new int[BUCKETS.length];
                    List<List<Long>> atByBucket = new ArrayList<>();
                    for (int i = 0; i < BUCKETS.length; i++) {
                        atByBucket.add(new ArrayList<Long>());
                    }
                    List<Long> totals = new ArrayList<>();
                    List<Long> belowSpecies = new ArrayList<>();
                    // The shares are of the sample's own in-genus total, so they are per sample and
                    // cannot be recovered from the sums above: a ratio of medians is not a median of
                    // ratios. Kept only for the ladder nodes, the ones whose candidate count needs
                    // guarding -- a species-level answer rests on the species and needs no warning.
                    List<Double> shareAt = new ArrayList<>();
                    List<Double> shareBelow = new ArrayList<>();
                    List<String> names = new ArrayList<>();
                    List<String> misses = new ArrayList<>();
                    long off = 0, offNamed = 0, unplaced = 0, nodes = 0, maxLadder = 0;
                    int rootSpecies = 0;
                    for (String sample : samples) {
                        Map<Currency, List<Call>> perCurrency = perSample.get(sample);
                        List<Call> calls = perCurrency == null ? null : perCurrency.get(currency);
                        if (calls == null || calls.isEmpty()) {
                            continue;
                        }
                        Call c = calls.get(0);
                        // What the top call says about the reference organism. The paper names these
                        // samples, so a run that changes which they are must change the sentence: the
                        // list is written out rather than left to be counted off the per-sample rows.
                        if (targetTaxId != null) {
                            if (targetTaxId.equals(c.nodeTaxId)) {
                                names.add(sample);
                            } else if (targetTaxId.equals(c.sigmaCulture)) {
                                misses.add(sample + " (\\emph{" + shortSpecies(c.nodeName) + "})");
                            }
                        }
                        totals.add(c.total);
                        unplaced += c.unplaced;
                        nodes += c.contributingNodes;
                        if (c.nearTaxId == null) {
                            off++;
                            // A sample with no signal on the organism's lineage has not thereby been
                            // shown free of it, and the paper calls this column `correctly off
                            // lineage' on the strength of the reference standard naming something
                            // else in every such sample. That is a property of the samples and not of
                            // the column, so it is counted rather than assumed: anything but zero
                            // here and the column has to be renamed.
                            if (targetTaxId != null
                                    && (targetTaxId.equals(c.sigmaCulture) || targetTaxId.equals(c.sigmaWimp))) {
                                offNamed++;
                            }
                            continue;
                        }
                        boolean isRoot = rootTaxId.equals(c.nearTaxId);
                        if (isRoot) {
                            rootSpecies = Math.max(rootSpecies, c.nearCandidates);
                        }
                        int bucket = bucketOf(c.nearCandidates, isRoot);
                        counts[bucket]++;
                        if (!isRoot && c.nearCandidates > 1 && c.total > 0) {
                            shareAt.add((double) c.nearAt / c.total);
                            shareBelow.add((double) c.nearBelow / c.total);
                        }
                        atByBucket.get(bucket).add(c.nearAt);
                        if (c.nearCandidates <= 1) {
                            belowSpecies.add(c.nearBelow);
                        } else if (!isRoot && c.nearCandidates <= COARSE_MAX) {
                            maxLadder = Math.max(maxLadder, c.nearCandidates);
                        }
                    }
                    ps.print(db); ps.print(';');
                    ps.print(currency.label); ps.print(';');
                    ps.print(variant.getLabel()); ps.print(';');
                    ps.print(samples.size()); ps.print(';');
                    ps.print(targetTaxId == null ? "" : targetTaxId); ps.print(';');
                    ps.print(rootSpecies == 0 ? "" : Integer.toString(rootSpecies)); ps.print(';');
                    for (int i = 0; i < BUCKETS.length; i++) {
                        ps.print(counts[i]); ps.print(';');
                        Long m = median(atByBucket.get(i));
                        ps.print(m == null ? "" : Long.toString(m)); ps.print(';');
                    }
                    ps.print(off); ps.print(';');
                    Long mt = median(totals);
                    ps.print(mt == null ? "" : Long.toString(mt)); ps.print(';');
                    Long mb = median(belowSpecies);
                    ps.print(mb == null ? "" : Long.toString(mb)); ps.print(';');
                    ps.print(maxLadder == 0 ? "" : Long.toString(maxLadder)); ps.print(';');
                    ps.print(unplaced); ps.print(';');
                    ps.print(nodes); ps.print(';');
                    ps.print(unchangedTopCalls(byVariant, samples)); ps.print(';');
                    ps.print(offNamed); ps.print(';');
                    Double msa = medianOf(shareAt), msb = medianOf(shareBelow);
                    ps.print(msa == null ? "" : String.format(Locale.ROOT, "%.4f", msa)); ps.print(';');
                    ps.print(msb == null ? "" : String.format(Locale.ROOT, "%.4f", msb)); ps.print(';');
                    ps.print(names.size()); ps.print(';');
                    ps.print(joinNames(names)); ps.print(';');
                    ps.print(misses.size()); ps.print(';');
                    ps.print(joinNames(misses)); ps.println(';');
                }
            }
        }
        System.out.println("Wrote " + file);
    }

    /** @return the index in {@link #BUCKETS} a node with this many candidates belongs to */
    static int bucketOf(int candidates, boolean isRoot) {
        if (isRoot || candidates > COARSE_MAX) {
            return BUCKETS.length - 1;
        }
        if (candidates <= 1) {
            return 0;
        }
        return candidates <= 10 ? 1 : candidates <= 50 ? 2 : 3;
    }

    /** {@link #median} for shares, which do not round to whole numbers. */
    static Double medianOf(List<Double> values) {
        if (values.isEmpty()) {
            return null;
        }
        List<Double> sorted = new ArrayList<>(values);
        Collections.sort(sorted);
        int n = sorted.size();
        return n % 2 == 1 ? sorted.get(n / 2) : (sorted.get(n / 2 - 1) + sorted.get(n / 2)) / 2;
    }

    /** {@code Streptococcus mitis} as {@code S.~mitis}, for a name the paper typesets in italics. */
    static String shortSpecies(String name) {
        if (name == null) {
            return "";
        }
        int space = name.indexOf(' ');
        return space <= 0 ? name : name.charAt(0) + ".~" + name.substring(space + 1);
    }

    /** {@code a, b and c}; the separator avoids the semicolon this file is delimited with. */
    static String joinNames(List<String> names) {
        if (names.isEmpty()) {
            return "";
        }
        if (names.size() == 1) {
            return names.get(0);
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < names.size() - 1; i++) {
            if (i > 0) {
                sb.append(", ");
            }
            sb.append(names.get(i));
        }
        return sb.append(" and ").append(names.get(names.size() - 1)).toString();
    }

    /** @return the middle value, averaged over the two middle ones for an even count, or null */
    static Long median(List<Long> values) {
        if (values.isEmpty()) {
            return null;
        }
        List<Long> sorted = new ArrayList<>(values);
        Collections.sort(sorted);
        int n = sorted.size();
        return n % 2 == 1 ? sorted.get(n / 2) : (sorted.get(n / 2 - 1) + sorted.get(n / 2)) / 2;
    }

    /**
     * How many samples the two variants give the same top call for, in every currency at once. It is
     * the figure behind the paper's claim that a classifier stated per sample cannot see a
     * refinement, so it is counted here rather than inferred from the per-sample rows afterwards.
     */
    private static int unchangedTopCalls(Map<Variant, Map<String, Map<Currency, List<Call>>>> byVariant,
                                         Set<String> samples) {
        int same = 0;
        for (String sample : samples) {
            boolean allEqual = true;
            for (Currency currency : Currency.values()) {
                String first = null;
                boolean seen = false;
                for (Variant variant : Variant.values()) {
                    Map<String, Map<Currency, List<Call>>> perSample = byVariant.get(variant);
                    Map<Currency, List<Call>> perCurrency = perSample == null ? null : perSample.get(sample);
                    List<Call> calls = perCurrency == null ? null : perCurrency.get(currency);
                    String name = calls == null || calls.isEmpty() || calls.get(0).nodeName == null
                            ? "" : calls.get(0).nodeName;
                    if (!seen) {
                        first = name;
                        seen = true;
                    } else if (!first.equals(name)) {
                        allEqual = false;
                    }
                }
            }
            if (allEqual) {
                same++;
            }
        }
        return same;
    }

    /** What one rule said about one sample, in one currency. */
    private static final class Call {
        String nodeTaxId;
        /** Dense position of the winning node, or -1 without a call; needed to ask for its candidates. */
        int nodePos = -1;
        String nodeName;
        String nodeRank;
        String taxon;
        int candidates;
        int contributingNodes;
        long total;
        /**
         * How much of {@link #total} reached no species: the counts at the genus, at refined nodes
         * and in the unranked buckets between them. The figure a refinement should lower.
         */
        long unplaced;
        /** The reference organism resolved to a taxon of this tree, per column; null where none is named. */
        String sigmaCulture;
        String sigmaWimp;
        /** The candidate precision against either column; null where that column names no organism of the subtree. */
        Double qCulture;
        Double qWimp;
        /**
         * The node nearest the reference organism that carries any signal at all, and what it says.
         * <p>
         * The ranked votes above answer "what does this sample look like"; these answer "how close to
         * the reference organism did its evidence get, and how many species are still open there".
         * They are the only place a refined node appears in this report: a refined node bears no name
         * of the reference taxonomy, so it can never win a vote that has to name a taxon, and every
         * read reaching one falls into {@link #unplaced}. {@link #nearCandidates} names it the only
         * way that means anything -- by how many species it still leaves in question.
         */
        String nearTaxId;
        String nearName;
        String nearRank;
        /** Species below the nearest node, the reference organism included; 1 at or below its species. */
        int nearCandidates;
        /** Signal whose lowest common ancestor is exactly the nearest node. */
        long nearAt;
        /** Signal at the nearest node or anywhere below it, i.e. what reached that depth. */
        long nearBelow;
    }

    /**
     * Matches the map's fastq files against one variant and votes on each.
     */
    private Map<String, Map<Currency, List<Call>>> run(String db, String fqMapFile, Variant variant,
                                                 String rootTaxId, String targetTaxId, long minimum,
                                                 final Map<String, List<String[]>> truth) throws IOException {
        final Map<String, Map<Currency, List<Call>>> result = new LinkedHashMap<>();
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
        FinerTreeMaker<FTProject> maker = new FinerTreeMaker<>(project);
        try {
            @SuppressWarnings("unchecked")
            ObjectGoal<Database, FTProject> dbGoal =
                    (ObjectGoal<Database, FTProject>) maker.getGoal(variant.getLoadDbGoalKey());
            final SmallTaxTree tree = dbGoal.get().getTaxTree();
            final TaxonComposition composition = new TaxonComposition(tree);
            final PathVoteTaxonModel rule = new PathVoteTaxonModel(composition, minimum);
            final SmallTaxTree.SmallTaxIdNode root = tree.getNodeByTaxId(rootTaxId);
            if (root == null) {
                throw new IllegalArgumentException("The database " + db + " has no node " + rootTaxId
                        + ", so there is no subtree to restrict the vote to.");
            }
            // Absent by choice leaves the columns empty; named but missing is a mistake worth saying,
            // since the whole point of the columns is that they are measured towards that organism.
            final SmallTaxTree.SmallTaxIdNode target =
                    targetTaxId == null || targetTaxId.isEmpty() ? null : tree.getNodeByTaxId(targetTaxId);
            if (targetTaxId != null && !targetTaxId.isEmpty() && target == null) {
                throw new IllegalArgumentException("The database " + db + " has no node " + targetTaxId
                        + ", so the nearest-node columns cannot be measured towards it.");
            }
            if (target != null && !isBelow(target, root)) {
                throw new IllegalArgumentException("Node " + targetTaxId + " does not lie below "
                        + rootTaxId + ", so no signal counted by this report can ever reach it.");
            }
            MatchResultGoal<?> matchResGoal = (MatchResultGoal<?>) maker.getGoal(variant.getMatchGoalKey());
            matchResGoal.setAfterMatchCallback(new MatchResultGoal.AfterMatchCallback() {
                @Override
                public void afterMatch(FastqKMerMatcher.MatcherReadEntry entry, boolean found) {
                    // Nothing per read: the vote is taken from the per-taxon totals below, which the
                    // matcher has already accumulated by the time afterKey() runs.
                }

                @Override
                public void afterKey(String key, MatchingResult res) {
                    Map<Currency, List<Call>> perCurrency = new LinkedHashMap<>();
                    for (Currency currency : Currency.values()) {
                        perCurrency.put(currency, votes(res, currency, tree, root, composition, rule,
                                minimum, target));
                    }
                    List<String[]> rows = truth.get(key);
                    Set<String> cultureSigmas = sigmasOf(rows, 2, tree, root);
                    Set<String> wimpSigmas = sigmasOf(rows, 3, tree, root);
                    for (List<Call> calls : perCurrency.values()) {
                        for (Call c : calls) {
                        c.sigmaCulture = String.join(",", cultureSigmas);
                        c.sigmaWimp = String.join(",", wimpSigmas);
                        c.qCulture = candidatePrecision(c, cultureSigmas, composition);
                        c.qWimp = candidatePrecision(c, wimpSigmas, composition);
                        }
                    }
                    result.put(key, perCurrency);
                    List<Call> readCalls = perCurrency.get(Currency.READS);
                    Call reads = readCalls.isEmpty() ? new Call() : readCalls.get(0);
                    // Only the read vote is echoed; the k-mer votes are in the CSV beside it.
                    System.out.println("  " + key + " [" + Currency.READS.label + "]: "
                            + (reads.taxon == null ? "no call" : reads.taxon)
                            + " at " + reads.nodeTaxId + " (" + reads.nodeName + "), "
                            + reads.total + " reads over " + reads.contributingNodes + " nodes");
                }
            });
            matchResGoal.make();
        } finally {
            maker.dumpAll();
        }
        return result;
    }

    /**
     * Turns one fastq file's per-taxon totals into the {@value #TOP_VOTES} heaviest votes, counting
     * only what lies below the root.
     * <p>
     * More than one, because a sputum sample carries a community rather than an organism and the
     * runner-up may be what the reference asks about -- see {@link PathVoteTaxonModel#classifyNodes}.
     * A sample the vote cannot answer at all still yields one row, with an empty node, so that it
     * appears in the report rather than silently missing from it.
     *
     * @param res         the matching result of one sample
     * @param currency    what the vote is taken over
     * @param tree        the taxonomy of the database being scored
     * @param root        the node the vote is restricted to
     * @param composition the species composition of that tree
     * @param rule        the vote
     * @param minimum     how much a winner's path must gather
     * @param target      the reference organism the nearest-node columns are measured towards, or null
     * @return one call per vote, heaviest first, never empty
     */
    private static List<Call> votes(MatchingResult res, Currency currency, SmallTaxTree tree,
                                    SmallTaxTree.SmallTaxIdNode root, TaxonComposition composition,
                                    PathVoteTaxonModel rule, long minimum,
                                    SmallTaxTree.SmallTaxIdNode target) {
        Map<Integer, Long> countsPerNode = new HashMap<>();
        // Kept beside the counts because the nearest-node walk below needs the nodes themselves and
        // not only their dense positions: there is no position-to-node lookup on the tree.
        List<SmallTaxTree.SmallTaxIdNode> signalNodes = new ArrayList<>();
        List<Long> signalCounts = new ArrayList<>();
        long total = 0;
        for (Map.Entry<String, CountsPerTaxid> e : res.getTaxid2Stats().entrySet()) {
            SmallTaxTree.SmallTaxIdNode node = tree.getNodeByTaxId(e.getKey());
            if (node == null || !isBelow(node, root)) {
                continue;
            }
            long value = currency.of(e.getValue());
            if (value <= 0) {
                continue;
            }
            countsPerNode.merge(node.getPosition(), value, Long::sum);
            signalNodes.add(node);
            signalCounts.add(value);
            total += value;
        }
        // Heaviest first; the position in this list is the vote's rank, and the report writes it
        // into the columns of that rank.
        List<Call> calls = new ArrayList<>();
        for (SmallTaxTree.SmallTaxIdNode winner : rule.classifyNodes(countsPerNode, minimum, TOP_VOTES)) {
            Call call = new Call();
            call.nodeTaxId = winner.getTaxId();
            call.nodePos = winner.getPosition();
            call.nodeName = winner.getName();
            call.nodeRank = winner.getRank() == null ? "" : winner.getRank().getName();
            call.candidates = composition.getCandidates(winner.getPosition());
            call.taxon = call.candidates == 1 ? composition.getMajorityClass(winner.getPosition()) : null;
            calls.add(call);
        }
        if (calls.isEmpty()) {
            calls.add(new Call());
        }
        Near near = nearest(target, signalNodes, signalCounts, minimum);
        long unplaced = rule.unplaced(countsPerNode);
        for (Call call : calls) {
            if (near != null) {
                call.nearTaxId = near.node.getTaxId();
                call.nearName = near.node.getName();
                call.nearRank = near.node.getRank() == null ? "" : near.node.getRank().getName();
                call.nearCandidates = composition.getCandidates(near.node.getPosition());
                call.nearAt = near.at;
                call.nearBelow = near.below;
            }
            // The same for every vote of a sample: they describe the sample, not the answer.
            call.total = total;
            call.contributingNodes = countsPerNode.size();
            call.unplaced = unplaced;
        }
        return calls;
    }

    /**
     * @param node the node to test
     * @param root the node it must lie below, itself included
     * @return whether it does
     */
    /**
     * The candidate precision of one call against one reference organism, i.e. the measure q of the
     * paper's section on classification quality, taken per sample instead of per read.
     * <p>
     * It is {@code 1/|Sigma(k)|} where the reference organism lies in {@code Sigma(k)}, the species
     * of the winning node's subtree, and {@code 0} where it does not. A call that names the species
     * outright has {@code |Sigma(k)| = 1} and scores one; a call left at a genus of forty species
     * scores a fortieth. This is what the earlier all-or-nothing rule threw away: it named a taxon
     * only for {@code |Sigma(k)| = 1} and reported nothing at all for every partial narrowing, which
     * is precisely the improvement a refinement makes.
     * <p>
     * A sample may name more than one organism. The call is credited if it covers any of them, since
     * each is present and a call that narrows down to one of them has narrowed down correctly.
     *
     * @param call        the call to score
     * @param sigmas      the reference organisms as taxids, empty where the column names no species
     * @param composition the species composition of the tree the call was made in
     * @return the candidate precision, or {@code null} where the column names no organism at all and
     * the measure is therefore undefined rather than zero
     */
    private static Double candidatePrecision(Call call, Set<String> sigmas, TaxonComposition composition) {
        return candidatePrecision(call.nodePos < 0 ? null : composition.getClassesAt(call.nodePos), sigmas);
    }

    /**
     * The measure itself, as a function of the two sets it is defined over, so that it can be checked
     * without a taxonomy.
     *
     * @param candidates the species of the winning node's subtree, or {@code null} where no call was made
     * @param sigmas     the reference organisms as taxids
     * @return the candidate precision, or {@code null} where {@code sigmas} is empty
     */
    static Double candidatePrecision(Set<String> candidates, Set<String> sigmas) {
        if (sigmas == null || sigmas.isEmpty()) {
            return null;
        }
        if (candidates == null || candidates.isEmpty()) {
            // No call was made, or the node carries no species at all: nothing was narrowed down, which
            // is a score of zero and not an absent value.
            return 0.0;
        }
        for (String sigma : sigmas) {
            if (candidates.contains(sigma)) {
                return 1.0 / candidates.size();
            }
        }
        return 0.0;
    }

    /**
     * Renders a candidate precision, leaving the field empty where the measure is undefined -- which
     * is not the same as zero and must not read as it.
     *
     * @param q the value, or {@code null} if undefined
     * @return the field's text
     */
    private static String format(Double q) {
        return q == null ? "" : String.format(Locale.ROOT, "%.6f", q);
    }

    /**
     * The reference organisms of one column that this vote could name at all: the taxids of the
     * sample's rows, kept only where the taxon lies in the subtree the vote is restricted to.
     * <p>
     * Without that restriction a sample whose reference reads {@code H. influenzae} would score zero
     * rather than count as undefined, and a database of streptococci would be marked wrong for not
     * naming an organism it does not contain. Over the 83 samples of {@code strepto} that is 52 of 58
     * samples scoring a structural zero, which no classifier could lift: the average would be unable
     * to show an improvement however good the calls became. Restricting the reference is the same
     * move the paper makes with its genus-only subset, and for the same reason.
     * <p>
     * A taxon the tree does not know is dropped as well. The database cannot name what it does not
     * hold, so such a sample is undecidable here rather than failed.
     *
     * @param rows       the sample's reference rows, or {@code null} if it has none
     * @param taxidField the field holding the resolved taxid of the column in question
     * @param tree       the taxonomy of the database being scored
     * @param root       the node the vote is restricted to
     * @return the reference organisms that lie below {@code root}, possibly empty
     */
    private static Set<String> sigmasOf(List<String[]> rows, int taxidField, SmallTaxTree tree,
                                        SmallTaxTree.SmallTaxIdNode root) {
        Set<String> res = new LinkedHashSet<>();
        if (rows != null) {
            for (String[] r : rows) {
                if (taxidField < r.length && !r[taxidField].isEmpty()) {
                    SmallTaxTree.SmallTaxIdNode node = tree.getNodeByTaxId(r[taxidField]);
                    if (node != null && isBelow(node, root)) {
                        res.add(r[taxidField]);
                    }
                }
            }
        }
        return res;
    }

    /** What {@link #nearest} found: a node and the two counts belonging to it. */
    static final class Near {
        SmallTaxTree.SmallTaxIdNode node;
        long at;
        long below;
    }

    /**
     * The node nearest the reference organism that carries any signal in this currency.
     * <p>
     * The subtree is searched first, and that order is not cosmetic. A refinement moves a species'
     * own $k$-mers down into strain clusters beneath it -- for {@code S. pneumoniae} it moves
     * 1,174,677 of 1,273,715, leaving 99,038 at the species node -- so the species node of a refined
     * database is routinely empty while everything below it is not. Walking upwards from the species
     * without looking below it first would climb past an emptied node to an ancestor spanning eight
     * species, and report a sample that had resolved to the organism as one that had not. The
     * refinement would be scored worst exactly where it works best.
     * <p>
     * At or below the species there is only ever one species in question, so the answer needs no
     * special case: {@link TaxonComposition#speciesOf} steps over the artificial DATA, FILE, ID and
     * REFINED nodes and resolves every genome under the species to the species, which makes
     * {@code getCandidates} one for the species node and for anything beneath it alike.
     *
     * @return the node and its two counts, or {@code null} if no signal reaches the organism's
     *         lineage at all -- which is not the same as a count of zero and is left empty, not zeroed
     */
    static Near nearest(SmallTaxTree.SmallTaxIdNode target,
                        List<SmallTaxTree.SmallTaxIdNode> signalNodes, List<Long> signalCounts,
                        long minimum) {
        if (target == null) {
            return null;
        }
        // Anything at or below the organism answers with the organism itself. Which of its strain
        // clusters carried the signal is a distinction this measure does not make, since all of them
        // leave the same one species open.
        long below = countAtOrBelow(target, signalNodes, signalCounts);
        if (meets(below, minimum)) {
            return near(target, signalNodes, signalCounts, below);
        }
        // Upwards from the organism, stopping at the first ancestor that holds enough itself. What
        // sits at or below that ancestor is asked for only once it is the answer.
        for (SmallTaxTree.SmallTaxIdNode n = target.getParent(); n != null; n = n.getParent()) {
            if (meets(countAt(n, signalNodes, signalCounts), minimum)) {
                return near(n, signalNodes, signalCounts, countAtOrBelow(n, signalNodes, signalCounts));
            }
        }
        return null;
    }

    /**
     * Whether a count is worth stopping the walk at.
     * <p>
     * Without a minimum a single read decides how deep the answer goes, and on the clinical runs it
     * did: the median count at the reported node was one read, so a sample could be said to have
     * reached the mitis complex on 0.009 per cent of its evidence. The candidate count is then
     * impressive and means nothing. The minimum is the same one the vote uses -- what counts as
     * enough evidence does not change between the two questions -- and it is a count in whatever
     * currency the row is written in, so one value has to serve reads and k-mers alike.
     *
     * @param count   what the node holds
     * @param minimum the floor, or {@link PathVoteTaxonModel#NO_MINIMUM} for none
     * @return whether the node holds anything at all, and enough of it
     */
    private static boolean meets(long count, long minimum) {
        return count > 0 && (minimum <= PathVoteTaxonModel.NO_MINIMUM || count >= minimum);
    }

    private static Near near(SmallTaxTree.SmallTaxIdNode node,
                             List<SmallTaxTree.SmallTaxIdNode> signalNodes, List<Long> signalCounts,
                             long below) {
        Near result = new Near();
        result.node = node;
        result.at = countAt(node, signalNodes, signalCounts);
        result.below = below;
        return result;
    }

    /**
     * @return the signal whose lowest common ancestor is exactly this node
     */
    static long countAt(SmallTaxTree.SmallTaxIdNode node,
                        List<SmallTaxTree.SmallTaxIdNode> signalNodes, List<Long> signalCounts) {
        long sum = 0;
        for (int i = 0; i < signalNodes.size(); i++) {
            if (signalNodes.get(i) == node) {
                sum += signalCounts.get(i);
            }
        }
        return sum;
    }

    /**
     * @return the signal at the node or anywhere below it. Node identity and not the dense position
     *         is what decides, so that this works on a hand-built lineage as well as on a loaded
     *         tree: positions are assigned by the tree's own traversal and a detached node has none.
     */
    static long countAtOrBelow(SmallTaxTree.SmallTaxIdNode node,
                               List<SmallTaxTree.SmallTaxIdNode> signalNodes, List<Long> signalCounts) {
        long sum = 0;
        for (int i = 0; i < signalNodes.size(); i++) {
            if (isBelow(signalNodes.get(i), node)) {
                sum += signalCounts.get(i);
            }
        }
        return sum;
    }

    private static boolean isBelow(SmallTaxTree.SmallTaxIdNode node, SmallTaxTree.SmallTaxIdNode root) {
        for (SmallTaxTree.SmallTaxIdNode n = node; n != null; n = n.getParent()) {
            if (n == root) {
                return true;
            }
        }
        return false;
    }

    /**
     * @param rows  the ground truth rows of one sample, or {@code null}
     * @param field 0 for the culture column, 1 for the source study's own call
     * @return the field of every row, comma separated, or the empty string
     */
    private static String join(List<String[]> rows, int field) {
        if (rows == null) {
            return "";
        }
        StringBuilder sb = new StringBuilder();
        for (String[] r : rows) {
            if (field < r.length && !r[field].isEmpty()) {
                if (sb.length() > 0) {
                    sb.append(", ");
                }
                sb.append(r[field]);
            }
        }
        return sb.toString();
    }

    /** Reads the transcribed per-sample reference standard, if the project has one. */
    private static final class GroundTruth {
        static Map<String, List<String[]>> read(File file) throws IOException {
            Map<String, List<String[]>> bySample = new LinkedHashMap<>();
            if (!file.exists()) {
                System.out.println("No " + file + " - the culture columns stay empty.");
                return bySample;
            }
            try (BufferedReader r = new BufferedReader(new InputStreamReader(
                    new FileInputStream(file), StandardCharsets.UTF_8))) {
                String line;
                while ((line = r.readLine()) != null) {
                    if (line.startsWith("#") || line.trim().isEmpty()) {
                        continue;
                    }
                    String[] p = line.split(";", -1);
                    if (p.length < 5) {
                        continue;
                    }
                    // culture, wimp, and the taxids they were resolved to (empty where the entry
                    // names no species). The names are for reading, the taxids for measuring: `S.
                    // aureus' is a Staphylococcus and must never be taken for a Streptococcus.
                    bySample.computeIfAbsent(p[0], k -> new ArrayList<>())
                            .add(new String[] { p[3], p[4],
                                    p.length > 5 ? p[5] : "", p.length > 6 ? p[6] : "" });
                }
            }
            return bySample;
        }
    }
}
