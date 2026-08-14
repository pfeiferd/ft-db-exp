package org.metagene.ftdbexp.stquality;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.FinerTreeMaker;
import org.metagene.genestrip.goals.refseq.RefSeqFnaFilesDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.TaxTree;

import java.io.File;
import java.util.Map;
import java.util.Set;

/**
 * The Genestrip-FT goal graph plus the sequence-type quality goals of this project.
 * <p>
 * Four goals are added, a counting goal and a CSV writer for each of the two databases, so that what
 * the refinement does to the measure can be read off the pair. They are registered here rather than in
 * genestrip-ft because they rest on an external typing of the genomes, which no Genestrip database
 * knows about; a subclass of the maker is all it takes, since the goals they depend on - the tax
 * nodes, the accession map, the two loaded databases - are the ordinary ones and are simply looked up.
 * <p>
 * The maker is also the only way to reach the execution context that supplies the reader threads,
 * {@code GSMaker.getExecutionContext} being protected.
 */
public class STQualityMaker extends FinerTreeMaker<FTProject> {
    /**
     * A key whose goal is dropped when something it depends on is cleaned, as Genestrip's own report
     * keys are ({@code ftdbinfo}, {@code ftsvgtaxtree}, ...).
     * <p>
     * {@code GoalKey.DefaultGoalKey} answers {@code false} to this, and a goal that answers
     * {@code false} is reached by neither {@code -t clean} nor {@code -t cleanall}: the former applies
     * only to the internal goal the maker aggregates the request into and never touches a real one,
     * and the latter descends only into dependencies that permit it. Such a CSV would then survive a
     * rebuilt database and describe the previous one.
     */
    private static final class ReportGoalKey implements GoalKey {
        private final String name;

        ReportGoalKey(String name) {
            this.name = name;
        }

        @Override
        public String getName() {
            return name;
        }

        @Override
        public boolean isTransClean() {
            return true;
        }

        @Override
        public String toString() {
            return name;
        }
    }

    /** Key of the counting goal over the unrefined database. */
    public static final GoalKey ST_QUALITY_COUNTS = new ReportGoalKey("stqualcounts");
    /** Key of the CSV goal over the unrefined database; also names {@code <db>_stquality.csv}. */
    public static final GoalKey ST_QUALITY = new ReportGoalKey("stquality");
    /** Key of the counting goal over the refined database. */
    public static final GoalKey FT_ST_QUALITY_COUNTS = new ReportGoalKey("ftstqualcounts");
    /** Key of the CSV goal over the refined database; also names {@code <db>_ftstquality.csv}. */
    public static final GoalKey FT_ST_QUALITY = new ReportGoalKey("ftstquality");

    private final STGroundTruth groundTruth;

    /**
     * Creates the maker.
     *
     * @param project     the FT project
     * @param groundTruth the sequence type of every genome the database was filled from
     */
    public STQualityMaker(FTProject project, STGroundTruth groundTruth) {
        super(project);
        this.groundTruth = groundTruth;
        // createGoals() has already run inside super(), before this field was set, so the goals below
        // could not have been registered there. They are added here instead, which is also why
        // registerSTGoals() is separate and not an override of createGoals().
        registerSTGoals();
    }

    @SuppressWarnings("unchecked")
    private void registerSTGoals() {
        FTProject project = getProject();

        ObjectGoal<Set<RefSeqCategory>, FTProject> categoriesGoal =
                (ObjectGoal<Set<RefSeqCategory>, FTProject>) getGoal(GSGoalKey.CATEGORIES);
        ObjectGoal<Set<TaxTree.TaxIdNode>, FTProject> taxNodesGoal =
                (ObjectGoal<Set<TaxTree.TaxIdNode>, FTProject>) getGoal(GSGoalKey.TAXNODES);
        ObjectGoal<TaxTree, FTProject> taxTreeGoal =
                (ObjectGoal<TaxTree, FTProject>) getGoal(GSGoalKey.TAXTREE);
        RefSeqFnaFilesDownloadGoal fnaFilesGoal = (RefSeqFnaFilesDownloadGoal) getGoal(GSGoalKey.REFSEQFNA);
        ObjectGoal<Map<File, TaxTree.TaxIdNode>, FTProject> additionalGoal =
                (ObjectGoal<Map<File, TaxTree.TaxIdNode>, FTProject>) getGoal(GSGoalKey.ADD_FASTAS);
        ObjectGoal<AccessionMap, FTProject> accessionMapGoal =
                (ObjectGoal<AccessionMap, FTProject>) getGoal(GSGoalKey.ACCMAP);
        ObjectGoal<Database, FTProject> dbGoal =
                (ObjectGoal<Database, FTProject>) getGoal(GSGoalKey.LOAD_DB);
        ObjectGoal<Database, FTProject> ftdbGoal =
                (ObjectGoal<Database, FTProject>) getGoal(FTGoalKey.LOAD_FTDB);

        // taxTreeGoal is a dependency only so that the large tree is not dropped before the readers
        // are done with it - the same reason DBQualityCountsGoal is wired to it in FinerTreeMaker.
        STQualityCountsGoal counts = new STQualityCountsGoal(project, ST_QUALITY_COUNTS,
                getExecutionContext(project), categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal,
                accessionMapGoal, dbGoal, groundTruth, taxTreeGoal);
        registerGoal(counts);
        registerGoal(new STQualityCSVGoal(project, ST_QUALITY, dbGoal, counts));

        STQualityCountsGoal ftCounts = new STQualityCountsGoal(project, FT_ST_QUALITY_COUNTS,
                getExecutionContext(project), categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal,
                accessionMapGoal, ftdbGoal, groundTruth, taxTreeGoal);
        registerGoal(ftCounts);
        registerGoal(new STQualityCSVGoal(project, FT_ST_QUALITY, ftdbGoal, ftCounts));
    }

    /**
     * Returns the goal writing the CSV for the given database variant.
     *
     * @param refined whether the refined database is meant
     * @return the CSV goal
     */
    public Goal<FTProject> getCSVGoal(boolean refined) {
        return getGoal(refined ? FT_ST_QUALITY : ST_QUALITY);
    }
}
