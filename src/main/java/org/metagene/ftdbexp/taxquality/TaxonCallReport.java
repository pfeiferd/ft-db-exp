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
 * Names, per sample and per database variant, the taxon the path vote of {@link PathVoteTaxonModel}
 * settles on below one node of the tree.
 * <p>
 * This is what the C.~difficile study's {@code IsolateSTAccuracyMain} did for isolates and sequence
 * types, applied to samples and taxa. The unit is deliberately not a read: for the Streptococcus case
 * study the question is what a respiratory specimen carries, and the reads bearing on it are the ones
 * the database placed anywhere below the genus. Restricting to that subtree is what makes a single
 * answer meaningful for a specimen that is otherwise polymicrobial.
 * <p>
 * Every sample is answered three times over, in reads, in matched {@code k}-mers and in distinct
 * matched {@code k}-mers. The three are not variants of the rule but three readings of the same
 * evidence, and they disagree in a way worth seeing: a single repetitive region can lend a species
 * many matched {@code k}-mers and one read, while a sample of low depth may have several reads and
 * few distinct {@code k}-mers. The C.~difficile run had them at 12, 22 and 5 hits of 74, so which
 * currency is used decided more there than the rule did.
 */
public class TaxonCallReport {
    private final File baseDir;
    private final File resultsDir;

    /** The currencies a sample's contribution is counted in. */
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
     * @param minimum    how much the winning path must gather, or {@link PathVoteTaxonModel#NO_MINIMUM}
     * @return the file that was written
     * @throws IOException if a database or a fastq file cannot be read, or the file cannot be written
     */
    public File write(String db, String fqMapFile, String reportKey, String rootTaxId, long minimum)
            throws IOException {
        Map<Variant, Map<String, Map<Currency, Call>>> byVariant = new LinkedHashMap<>();
        for (Variant variant : Variant.values()) {
            System.out.println("=== " + db + " / " + variant.getLabel() + " ===");
            byVariant.put(variant, run(db, fqMapFile, variant, rootTaxId, minimum));
        }
        Map<String, List<String[]>> truth = GroundTruth.read(new File(baseDir, "projects/" + db + "/ground_truth.csv"));

        File file = new File(resultsDir, db + "_" + reportKey + "_taxoncall.csv");
        try (PrintStream ps = new PrintStream(new FileOutputStream(file), false, StandardCharsets.UTF_8.name())) {
            ps.println("db;sample;variant;currency;node;node name;node rank;taxon;candidates;"
                    + "contributing nodes;total;culture;wimp;");
            Set<String> samples = new LinkedHashSet<>();
            for (Map<String, Map<Currency, Call>> m : byVariant.values()) {
                samples.addAll(m.keySet());
            }
            for (String sample : samples) {
                for (Variant variant : Variant.values()) {
                    Map<Currency, Call> perCurrency = byVariant.get(variant).get(sample);
                    if (perCurrency == null) {
                        continue;
                    }
                    for (Currency currency : Currency.values()) {
                        Call c = perCurrency.get(currency);
                        if (c == null) {
                            continue;
                        }
                        ps.print(db); ps.print(';');
                        ps.print(sample); ps.print(';');
                        ps.print(variant.getLabel()); ps.print(';');
                        ps.print(currency.label); ps.print(';');
                        ps.print(c.nodeTaxId == null ? "" : c.nodeTaxId); ps.print(';');
                        ps.print(c.nodeName == null ? "" : c.nodeName); ps.print(';');
                        ps.print(c.nodeRank == null ? "" : c.nodeRank); ps.print(';');
                        ps.print(c.taxon == null ? "" : c.taxon); ps.print(';');
                        ps.print(c.candidates); ps.print(';');
                        ps.print(c.contributingNodes); ps.print(';');
                        ps.print(c.total); ps.print(';');
                        ps.print(join(truth.get(sample), 0)); ps.print(';');
                        ps.print(join(truth.get(sample), 1)); ps.println(';');
                    }
                }
            }
        }
        System.out.println("Wrote " + file);
        return file;
    }

    /** What one rule said about one sample, in one currency. */
    private static final class Call {
        String nodeTaxId;
        String nodeName;
        String nodeRank;
        String taxon;
        int candidates;
        int contributingNodes;
        long total;
    }

    /**
     * Matches the map's fastq files against one variant and votes on each.
     */
    private Map<String, Map<Currency, Call>> run(String db, String fqMapFile, Variant variant,
                                                 String rootTaxId, long minimum) throws IOException {
        final Map<String, Map<Currency, Call>> result = new LinkedHashMap<>();
        FTProject project = new FTProject(new GSCommon(baseDir), db, null, null, fqMapFile,
                null, null, null, null, null, null, false);
        project.initConfigParam(GSConfigKey.THREADS, -1);
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
            MatchResultGoal<?> matchResGoal = (MatchResultGoal<?>) maker.getGoal(variant.getMatchGoalKey());
            matchResGoal.setAfterMatchCallback(new MatchResultGoal.AfterMatchCallback() {
                @Override
                public void afterMatch(FastqKMerMatcher.MatcherReadEntry entry, boolean found) {
                    // Nothing per read: the vote is taken from the per-taxon totals below, which the
                    // matcher has already accumulated by the time afterKey() runs.
                }

                @Override
                public void afterKey(String key, MatchingResult res) {
                    Map<Currency, Call> perCurrency = new LinkedHashMap<>();
                    for (Currency currency : Currency.values()) {
                        perCurrency.put(currency, vote(res, currency, tree, root, composition, rule, minimum));
                    }
                    result.put(key, perCurrency);
                    Call reads = perCurrency.get(Currency.READS);
                    System.out.println("  " + key + ": " + (reads.taxon == null ? "no call" : reads.taxon)
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
     * Turns one fastq file's per-taxon totals into a vote, counting only what lies below the root.
     */
    private static Call vote(MatchingResult res, Currency currency, SmallTaxTree tree,
                             SmallTaxTree.SmallTaxIdNode root, TaxonComposition composition,
                             PathVoteTaxonModel rule, long minimum) {
        Map<Integer, Long> countsPerNode = new HashMap<>();
        Call call = new Call();
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
            call.total += value;
        }
        call.contributingNodes = countsPerNode.size();
        SmallTaxTree.SmallTaxIdNode winner = rule.classifyNode(countsPerNode, minimum);
        if (winner != null) {
            call.nodeTaxId = winner.getTaxId();
            call.nodeName = winner.getName();
            call.nodeRank = winner.getRank() == null ? "" : winner.getRank().getName();
            call.candidates = composition.getCandidates(winner.getPosition());
            call.taxon = call.candidates == 1 ? composition.getMajorityClass(winner.getPosition()) : null;
        }
        return call;
    }

    /**
     * @param node the node to test
     * @param root the node it must lie below, itself included
     * @return whether it does
     */
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
                    bySample.computeIfAbsent(p[0], k -> new ArrayList<>())
                            .add(new String[] { p[3], p[4] });
                }
            }
            return bySample;
        }
    }
}
