package org.metagene.ftdbexp.stquality;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Collections;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Writes the per-tax-id sequence-type level precisions of {@link STQualityCountsGoal} to a CSV, in the
 * order of the taxonomy tree.
 * <p>
 * The columns carry the same names as the corresponding ones of Genestrip-FT's
 * {@code DBQualityCSVGoal}, so that the two files can be read by the same tooling and compared column
 * for column; what differs is the unit behind them, which is a sequence type here and a data taxon
 * there. Semicolon-separated, as every other CSV of this project.
 */
public class STQualityCSVGoal extends FileGoal<FTProject> {
    private static final DecimalFormat DF = new DecimalFormat("0.00000000", new DecimalFormatSymbols(Locale.US));

    private final ObjectGoal<Database, FTProject> storeGoal;
    private final ObjectGoal<Map<String, STCounts>, FTProject> countsGoal;

    /**
     * Creates the goal.
     *
     * @param project    the FT project
     * @param key        the goal key, which also names the output file
     * @param storeGoal  the goal providing the database, for its taxonomy tree
     * @param countsGoal the goal providing the per-tax-id tallies
     * @param deps       further goals this goal depends on
     */
    @SafeVarargs
    public STQualityCSVGoal(FTProject project, GoalKey key, ObjectGoal<Database, FTProject> storeGoal,
                            ObjectGoal<Map<String, STCounts>, FTProject> countsGoal, Goal<FTProject>... deps) {
        super(project, key, Goal.append(deps, storeGoal, countsGoal));
        this.storeGoal = storeGoal;
        this.countsGoal = countsGoal;
    }

    @Override
    public List<File> getFiles() {
        return Collections.singletonList(
                getProject().getOutputFile(getKey().getName(), GSProject.GSFileType.CSV, false));
    }

    @Override
    protected void makeFile(File file) throws IOException {
        Map<String, STCounts> counts = countsGoal.get();
        SmallTaxTree tree = storeGoal.get().getTaxTree();

        try (PrintStream ps = new PrintStream(file, StandardCharsets.UTF_8)) {
            // "sequence types" is |S_n|, the denominator of every precision in the row. "kmers at
            // node" and "st pairs" are the two halves of the node precision: the latter is the sum of
            // c_st over the former. The subtree columns average p_st over the k-mers of the whole
            // subtree, and the last two restrict that to the k-mers stored above the data taxa, which
            // is where a refinement can act at all.
            // The recall columns pool over the lineages below a node - "avg recall" giving each of
            // them one vote whatever its size, "recall" weighting by the k-mers they contribute.
            ps.println("taxid;name;rank;parent taxid;sequence types;kmers at node;st pairs;node precision;"
                    + "subtree kmers;subtree precision;subtree kmers above data;restricted subtree precision;"
                    + "tp;tp+fn;recall;avg recall");
            for (SmallTaxTree.SmallTaxIdNode node : tree) {
                STCounts c = counts.get(node.getTaxId());
                // A node with no typed genome underneath carries no measure - not a bad one. That
                // includes the OTHER nodes the refinement introduces and every branch of the taxonomy
                // the typing does not reach, and they are left out rather than reported as NaN.
                if (c == null || c.getSTCount() == 0) {
                    continue;
                }
                SmallTaxTree.SmallTaxIdNode parent = node.getParent();
                ps.print(node.getTaxId());
                ps.print(";");
                ps.print(node.getName());
                ps.print(";");
                ps.print(node.getRank() == null ? "null" : node.getRank().getName());
                ps.print(";");
                ps.print(parent == null ? "null" : parent.getTaxId());
                ps.print(";");
                ps.print(c.getSTCount());
                ps.print(";");
                ps.print(c.getKmerSumForNode());
                ps.print(";");
                ps.print(c.getTpForNodePrecision());
                ps.print(";");
                print(ps, c.getNodePrecision());
                ps.print(";");
                ps.print(c.getSubtreeKmerSum());
                ps.print(";");
                print(ps, c.getSubtreePrecision());
                ps.print(";");
                ps.print(c.getSubtreeKmersAboveData());
                ps.print(";");
                print(ps, c.getRestrictedSubtreePrecision());
                ps.print(";");
                ps.print(c.getTp());
                ps.print(";");
                ps.print(c.getTpPlusFn());
                ps.print(";");
                print(ps, c.getRecall());
                ps.print(";");
                print(ps, c.getAvgRecall());
                ps.println(";");
            }
        }
    }

    /**
     * Prints a precision, leaving the field empty where the measure is undefined. An empty field says
     * "not applicable" to every reader of a CSV; the string {@code NaN} would have to be special-cased
     * by each of them, and read as a number by the ones that forget.
     *
     * @param ps    the stream to print to
     * @param value the value, possibly {@link Double#NaN}
     */
    private static void print(PrintStream ps, double value) {
        if (!Double.isNaN(value)) {
            ps.print(DF.format(value));
        }
    }
}
