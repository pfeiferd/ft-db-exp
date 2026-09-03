package org.metagene.ftdbexp.prescreen;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Decides whether a genus is worth refining at all, from the branching-degree histogram alone.
 * <p>
 * The refinement can only move a $k$-mer that few of a node's children carry. One carried by a single
 * child is not there to begin with: the LCA update would have pushed it down to that child, which is
 * why the degree-one column of every histogram is zero. One carried by nearly all of them is core
 * genome and cannot move whatever the clustering does. What is left in between is the headroom, and
 * measuring it costs an unrefined database and a $k$-mer index rather than the refinement itself.
 * <p>
 * This exists because the <em>C. difficile</em> study was built three times before anyone measured
 * that its headroom was a few per cent. The rule of thumb it is used with: a genus whose degree
 * {@code 2-10} share is in the low single digits will not repay a refinement, and the right response
 * is to walk away rather than to build again and hope.
 * <p>
 * Reads {@code <db>_branchhistocsv.csv}, whose rows are {@code taxid;name;children;d1;...;dC;OTHER}:
 * one column per branching degree from one up to the node's child count, and a last column for the
 * OTHER bucket, which holds what was pushed up by a taxon outside the database entirely.
 */
public class PrescreenReport {
    /** Where the headroom is taken to stop. A node reached by more than ten children names too many. */
    private static final int HEADROOM_MAX = 10;

    private final File resultsDir;

    public PrescreenReport(File resultsDir) {
        this.resultsDir = resultsDir;
    }

    /** One node of the histogram, reduced to the four shares the decision rests on. */
    public static final class Verdict {
        public String taxId;
        public String name;
        public int children;
        public long total;
        public double singleShare;
        public double headroomShare;
        public double coreShare;
        public double otherShare;
    }

    /**
     * @param row the histogram row, already split on the delimiter
     * @return the verdict for it, or {@code null} if the row carries no k-mers or cannot be read
     */
    static Verdict verdictOf(String[] row) {
        if (row.length < 4) {
            return null;
        }
        Verdict v = new Verdict();
        v.taxId = row[0].trim();
        v.name = row[1].trim();
        try {
            v.children = Integer.parseInt(row[2].trim());
        } catch (NumberFormatException e) {
            return null;
        }
        List<Long> degrees = new ArrayList<>();
        for (int i = 3; i < row.length; i++) {
            String cell = row[i].trim();
            if (cell.isEmpty()) {
                continue;
            }
            try {
                degrees.add(Long.parseLong(cell));
            } catch (NumberFormatException e) {
                return null;
            }
        }
        if (degrees.isEmpty()) {
            return null;
        }
        // The last column is the OTHER bucket, not a branching degree, so it is taken off before the
        // degrees are read and added back only to the total: a k-mer there was pushed up by something
        // the database does not contain, and no refinement of these children can place it.
        long other = degrees.remove(degrees.size() - 1);
        long total = other;
        for (long d : degrees) {
            total += d;
        }
        if (total == 0) {
            return null;
        }
        long single = degrees.isEmpty() ? 0 : degrees.get(0);
        long headroom = 0;
        for (int d = 2; d <= Math.min(HEADROOM_MAX, degrees.size()); d++) {
            headroom += degrees.get(d - 1);
        }
        // Core: carried by more than half the children. Nothing above this can be moved down, so it
        // is reported to show how much of the node is fixed before the clustering starts.
        long core = 0;
        for (int d = degrees.size(); d > v.children / 2 && d >= 1; d--) {
            core += degrees.get(d - 1);
        }
        v.total = total;
        v.singleShare = 100.0 * single / total;
        v.headroomShare = 100.0 * headroom / total;
        v.coreShare = 100.0 * core / total;
        v.otherShare = 100.0 * other / total;
        return v;
    }

    /**
     * Writes {@code <db>_prescreen.csv}, one row per node of the histogram, and prints the verdict
     * for the heaviest node so that a decision can be taken from the console.
     *
     * @param db the database project name
     * @return the file that was written
     * @throws IOException if the histogram cannot be read or the report cannot be written
     */
    public File write(String db) throws IOException {
        File in = new File(resultsDir, db + "_branchhistocsv.csv");
        if (!in.exists()) {
            throw new IOException("No " + in + ". Run the goal `branchhistocsv' for " + db + " first;"
                    + " it needs the unrefined database and the k-mer index, and not the refinement.");
        }
        List<Verdict> verdicts = new ArrayList<>();
        try (BufferedReader r = new BufferedReader(
                new InputStreamReader(new FileInputStream(in), StandardCharsets.UTF_8))) {
            String line = r.readLine();
            while ((line = r.readLine()) != null) {
                if (line.trim().isEmpty()) {
                    continue;
                }
                Verdict v = verdictOf(line.split(";", -1));
                if (v != null) {
                    verdicts.add(v);
                }
            }
        }
        File out = new File(resultsDir, db + "_prescreen.csv");
        try (PrintStream ps = new PrintStream(new FileOutputStream(out), false, StandardCharsets.UTF_8.name())) {
            ps.println("db;tax id;name;children;kmers;degree 1 %;degree 2-" + HEADROOM_MAX
                    + " %;core %;other %;");
            for (Verdict v : verdicts) {
                ps.print(db); ps.print(';');
                ps.print(v.taxId); ps.print(';');
                ps.print(v.name); ps.print(';');
                ps.print(v.children); ps.print(';');
                ps.print(v.total); ps.print(';');
                ps.print(pct(v.singleShare)); ps.print(';');
                ps.print(pct(v.headroomShare)); ps.print(';');
                ps.print(pct(v.coreShare)); ps.print(';');
                ps.print(pct(v.otherShare)); ps.println(';');
            }
        }
        System.out.println("Wrote " + out);
        Verdict heaviest = null;
        for (Verdict v : verdicts) {
            if (heaviest == null || v.total > heaviest.total) {
                heaviest = v;
            }
        }
        if (heaviest != null) {
            System.out.println();
            System.out.println("############ " + db + ": prescreen verdict ############");
            System.out.printf(Locale.ROOT, "  heaviest node   %s (%s), %d children, %,d k-mers%n",
                    heaviest.name, heaviest.taxId, heaviest.children, heaviest.total);
            System.out.printf(Locale.ROOT, "  degree 1        %6s %%  (empty by construction)%n",
                    pct(heaviest.singleShare));
            System.out.printf(Locale.ROOT, "  degree 2-%-2d     %6s %%  <- the headroom a refinement can claim%n",
                    HEADROOM_MAX, pct(heaviest.headroomShare));
            System.out.printf(Locale.ROOT, "  core (> C/2)    %6s %%  (cannot move)%n", pct(heaviest.coreShare));
            System.out.printf(Locale.ROOT, "  OTHER           %6s %%  (pushed up from outside the database)%n",
                    pct(heaviest.otherShare));
            System.out.println();
            System.out.println(heaviest.headroomShare < 10
                    ? "  VERDICT: walk away. The headroom is in the low single digits, which is what"
                      + " C. difficile\n           looked like before it was built three times."
                    : "  VERDICT: worth refining. Build the refinement and measure it properly.");
        }
        return out;
    }

    private static String pct(double v) {
        return String.format(Locale.ROOT, "%.2f", v);
    }
}
