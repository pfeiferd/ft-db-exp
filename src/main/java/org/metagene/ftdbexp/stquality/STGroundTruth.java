package org.metagene.ftdbexp.stquality;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * The sequence type of each genome the database was filled from, read from the CSV that
 * {@code bin/mlst_assemblies.sh} writes.
 * <p>
 * This is the unit the measures of {@link STQualityCountsGoal} count, and it is external to the
 * taxonomy: a sequence type is not a node, so it cannot be read off the tree and has to be joined in
 * from a typing. The join key is the name of the leaf a genome became - which is the name of the
 * fasta file it was read from, since {@code fileNodes=true} names an artificial FILE node after it
 * ({@code TaxTree.fileNode(node, file.getName(), ...)}). The typing script keys its rows by exactly
 * that string, so the join is by equality and needs nothing reconciled.
 * <p>
 * A genome that types to no profile carries no unit. The paper says what follows: it is left out of
 * {@code S_n} rather than counted as a type of its own, which makes {@code |S_n|} smaller and the
 * precision correspondingly more generous, so the share of untyped genomes belongs beside any figure
 * reported from this. {@link #getUntypedLeaves()} is that share's numerator.
 * <p>
 * Sequence types are handed out dense indices here, because the counting pass needs to address them
 * as bits of a set per node and as a key half in a duplicate filter.
 */
public class STGroundTruth {
    /** The value {@code mlst} writes when a profile is novel or incomplete. */
    private static final String NO_TYPE = "-";

    private final Map<String, Integer> leafToST = new HashMap<>();
    private final List<String> stNames = new ArrayList<>();
    private final Map<String, Integer> stToIndex = new HashMap<>();
    private int untypedLeaves;

    /**
     * Reads the typing from the given CSV.
     * <p>
     * The columns are located by their header names rather than by position, and the key column may
     * be called either {@code leaf} or {@code file}: the script used to write one row per sequence
     * accession under the latter name, and a file from that era should be rejected for having the
     * wrong keys in it, not silently misread because a column moved.
     *
     * @param csv the CSV written by {@code bin/mlst_assemblies.sh}
     * @throws IOException           if the file cannot be read
     * @throws IllegalStateException if it carries neither of the expected key columns, or no ST column
     */
    public STGroundTruth(File csv) throws IOException {
        try (BufferedReader r = new BufferedReader(
                new InputStreamReader(new FileInputStream(csv), StandardCharsets.UTF_8))) {
            String header = r.readLine();
            if (header == null) {
                throw new IllegalStateException("Empty sequence type file: " + csv);
            }
            String[] cols = header.split(";", -1);
            int keyCol = indexOf(cols, "leaf");
            if (keyCol < 0) {
                keyCol = indexOf(cols, "file");
            }
            int stCol = indexOf(cols, "st");
            if (keyCol < 0 || stCol < 0) {
                throw new IllegalStateException("Expected a 'leaf' (or 'file') and an 'st' column in "
                        + csv + ", but its header is: " + header);
            }
            String line;
            while ((line = r.readLine()) != null) {
                if (line.isEmpty()) {
                    continue;
                }
                String[] f = line.split(";", -1);
                if (f.length <= keyCol || f.length <= stCol) {
                    continue;
                }
                String leaf = f[keyCol].trim();
                String st = f[stCol].trim();
                if (leaf.isEmpty()) {
                    continue;
                }
                if (st.isEmpty() || NO_TYPE.equals(st)) {
                    untypedLeaves++;
                    continue;
                }
                Integer index = stToIndex.get(st);
                if (index == null) {
                    index = stNames.size();
                    stToIndex.put(st, index);
                    stNames.add(st);
                }
                // A leaf named twice would mean the typing disagrees with itself; the first row wins
                // and the second is dropped rather than silently replacing it.
                leafToST.putIfAbsent(leaf, index);
            }
        }
    }

    private static int indexOf(String[] cols, String name) {
        for (int i = 0; i < cols.length; i++) {
            if (name.equalsIgnoreCase(cols[i].trim())) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Returns the index of the sequence type of the genome that became the given leaf.
     *
     * @param leafName the name of the leaf, i.e. of the fasta file the genome was read from
     * @return the dense index of its sequence type, or {@code -1} if the genome is untyped or absent
     * from the typing altogether
     */
    public int getSTIndex(String leafName) {
        Integer index = leafToST.get(leafName);
        return index == null ? -1 : index;
    }

    /**
     * Returns the name of a sequence type by its index.
     *
     * @param stIndex the dense index
     * @return the sequence type as the typing scheme names it
     */
    public String getSTName(int stIndex) {
        return stNames.get(stIndex);
    }

    /**
     * Returns the number of distinct sequence types in the typing.
     *
     * @return the number of distinct sequence types, i.e. the size of the index space
     */
    public int getSTCount() {
        return stNames.size();
    }

    /**
     * Returns the number of typed genomes.
     *
     * @return the number of leaves that carry a sequence type
     */
    public int getTypedLeaves() {
        return leafToST.size();
    }

    /**
     * Returns the number of genomes the typing lists without a sequence type, i.e. those whose profile
     * is novel or incomplete. They contribute no unit to any {@code S_n}.
     *
     * @return the number of untyped genomes in the typing
     */
    public int getUntypedLeaves() {
        return untypedLeaves;
    }

    /**
     * Returns the sequence types, indexed by the index this class assigns them.
     *
     * @return an unmodifiable view of the sequence type names in index order
     */
    public List<String> getSTNames() {
        return Collections.unmodifiableList(stNames);
    }
}
