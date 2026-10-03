#!/bin/sh
#
# Runs the classification experiments, in two independent parts:
#
#   quality      every simulated fastq file is matched against the unrefined and against the
#                refined database, and the classifications are scored against the reads' known
#                ground truth. The measures are the ones of Section "Classification quality" in
#                the Genestrip-FT paper: the boolean positive counts at the genus and the species
#                rank, plus the candidate-weighted species count, which credits a classification in
#                proportion to how far it narrows the species down. The latter is the one that can
#                tell a refined database from an unrefined one at all.
#
#   performance  wall time and maximum RAM of classifying the same files, per database variant,
#                for Section "Performance" of the paper. Measured through the plain goals `match'
#                and `ftmatch' rather than derived from the quality runs -- see run_perf() below.
#
# Usage:
#   sh ./bin/run_classification_exps.sh [viral|protozoa|strepto|nocardia|tick-borne|saliva|mngs|kraken|k24|accuracy|perf|real|all]
#
#   ERROR_SALIVA=1 sh ./bin/run_classification_exps.sh viral    # the saliva-matched read set
#
#   A database name runs both parts for it, following ERROR_FREE; `accuracy' and `perf' run one
#   part for every database, and `accuracy' and `all' evaluate the two Illumina read sets and the
#   error-free one of every InSilicoSeq project. The Nanopore regimes of make_fastqs.sh are not
#   among them -- NanoSim covers that ground -- but ERROR_NANOPORE=1 or ERROR_NANOPORE_LONG=1 still
#   scores them if the reads exist.
#
# Prerequisites:
#   sh ./bin/install_tools.sh     installs InSilicoSeq, NanoSim and cgmemtime
#   sh ./bin/run_exps.sh          builds the databases (or at least the goals `db` and `ftdb`)
#   sh ./bin/make_fastqs.sh       generates the simulated reads
#
# Results land in ./results as <db>_<report key>_accuracy.csv and as match_<db>.log and
# ftmatch_<db>.log; they are meant to be copied into the paper's own results folder, from where
# they are included directly.
#
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)

what=${1:-all}
res_path="${basedir}/results"
mkdir -p "$res_path"

# $1 = database project name. ERROR_FREE selects the reads generated without any sequencing error,
# which are kept apart from the regular ones down to the name of the result file.
run_iss() {
  db=$1
  if [ -n "${ERROR_FREE:-}" ]; then
    mapname="${db}_sim_perfect.txt"
    reportkey="iss_perfect"
    what="error-free InSilicoSeq reads"
  elif [ -n "${ERROR_NANOPORE:-}" ]; then
    mapname="${db}_sim_nanopore.txt"
    reportkey="iss_nanopore"
    what="InSilicoSeq reads at a Nanopore-level per-base error"
  elif [ -n "${ERROR_NANOPORE_LONG:-}" ]; then
    mapname="${db}_sim_nanoporelong.txt"
    reportkey="iss_nanopore_long"
    what="long InSilicoSeq reads at a Nanopore-level per-base error"
  elif [ -n "${ERROR_SALIVA:-}" ]; then
    mapname="${db}_sim_saliva.txt"
    reportkey="iss_saliva"
    what="InSilicoSeq reads matched to the human saliva runs"
  elif [ -n "${ERROR_MNGS:-}" ]; then
    mapname="${db}_sim_mngs.txt"
    reportkey="iss_mngs"
    what="InSilicoSeq reads matched to the BGISEQ mNGS runs"
  else
    mapname="${db}_sim.txt"
    reportkey="iss"
    what="InSilicoSeq reads"
  fi
  map="${basedir}/data/fastq/${mapname}"
  if [ ! -f "$map" ]; then
    echo "Missing ${map} - run 'sh ./bin/make_fastqs.sh ${db}' first." >&2
    return 1
  fi
  echo "############ ${db}: ${what}, unrefined vs. refined ############"
  mvn exec:exec@accuracy -Dname="$db" -Dfqmap="$mapname" -Dreportkey="$reportkey" -Dsimulator=ISS
}

# The read sets the paper reports: the two Illumina models and the error-free one.
# See make_fastqs.sh for what each of them is for, and why the Nanopore ones are not included.
#
# Every regime variable run_iss reads is cleared here, not just the ones this function sets, and for
# the same reason make_fastqs.sh clears them: run_iss takes the first non-empty one, so one left over
# in the environment redirects these calls without saying so. The saliva case is the one that bites.
# `ERROR_SALIVA=1 run_classification_exps.sh all' would score the saliva-matched reads under the
# report key iss_saliva where the plain Illumina run belongs, leaving no viral_iss_summary.csv at
# all -- and `all' goes on to call run_real for the saliva samples, which reads exactly that report
# key as its calibration, so the damage would surface as a calibration silently taken from the wrong
# read set rather than as a failure. Run the saliva regime on its own:
#
#   ERROR_SALIVA=1 sh ./bin/run_classification_exps.sh viral
#
# and run it *before* the `real' part, since that is what produces the calibration `real' consumes.
run_iss_all_regimes() {
  saved_error_free=${ERROR_FREE:-}
  saved_nanopore=${ERROR_NANOPORE:-}
  saved_nanopore_long=${ERROR_NANOPORE_LONG:-}
  saved_saliva=${ERROR_SALIVA:-}
  ERROR_FREE=""; ERROR_NANOPORE=""; ERROR_NANOPORE_LONG=""; ERROR_SALIVA=""
  run_iss "$1"
  ERROR_FREE=1; ERROR_NANOPORE=""; ERROR_NANOPORE_LONG=""; ERROR_SALIVA=""
  run_iss "$1"
  # No Nanopore regimes here either, for the reason given in make_fastqs.sh: NanoSim replaces them.
  ERROR_FREE=$saved_error_free
  ERROR_NANOPORE=$saved_nanopore
  ERROR_NANOPORE_LONG=$saved_nanopore_long
  ERROR_SALIVA=$saved_saliva
}

# The experiments on real reads, which have no ground truth: `cv' against the five human saliva runs
# of the first paper and `tb' against its tick samples. Neither precision nor recall is defined here,
# so this runs the `specificity' entry point instead of `accuracy': it reports how far each database
# variant narrows the species down on the reads the unrefined one left at a genus. The difference
# between the two is not a bound on the precision gain; it becomes an estimate of it only through the
# calibration factors rho_u and rho_f of a simulated run -- see SpecificityReport.
#
# $1 = database project name, $2 = fastq mapping file, $3 = report key,
# $4 = report key of the simulated run supplying the calibration rho (optional). For the ticks that
# is `nanosim': its simulation named tickN was trained on the real sample named tickN, so joining
# the two by fastq key gives each estimate the calibration derived from its own sample.
run_real() {
  map="${basedir}/data/fastq/$2"
  if [ ! -f "$map" ]; then
    echo "Missing ${map}." >&2
    return 1
  fi
  # A map that declares its runs by URL -- ticks_real.txt and the two PRJEB30781 maps of `strepto' --
  # can fetch what it names, through Genestrip's own goal `fastqdownload' rather than a script of
  # ours: the goal resolves, downloads and names the files itself, and skips whatever is in place.
  # Maps that name local files instead (the generated simulation maps, saliva_real.txt) have no URL
  # in the second column and are left alone; saliva is fetched by bin/fetch_saliva.sh, which exists
  # because those runs are paired and need their URLs resolved per accession.
  if awk '$1 !~ /^#/ && $2 ~ /^https?:/ {found=1} END {exit !found}' "$map"; then
    missing=""
    for key in $(awk '$1 !~ /^#/ && NF >= 2 {print $1}' "$map" | sort -u); do
      [ -s "${basedir}/data/fastq/${key}.fastq.gz" ] || missing="yes"
    done
    if [ -n "$missing" ]; then
      echo "=== $1: fetching the runs of $2 ==="
      ( cd "$basedir" && mvn exec:exec@fastqdl -Dname="$1" -Dfqmap="$2" )
    fi
  fi
  # Every file the map names must be on disk. Genestrip drops a missing one silently and classifies
  # what is left, which for a paired run means one mate: SRR5571991 was analysed that way on
  # 2026-09-19 because its second mate had stopped downloading at 426 MB of 69 GB, and the result --
  # half the reads of a run whose siblings were counted whole -- looked like a measurement. A key
  # that resolves to fewer files than its siblings is reported for the same reason.
  _rr_absent=""
  for _rr_f in $(awk '!/^#/ && NF >= 2 && $2 !~ /^https?:/ { print $2 }' "$map"); do
    [ -s "${basedir}/data/fastq/${_rr_f}" ] || _rr_absent="${_rr_absent} ${_rr_f}"
  done
  if [ -n "$_rr_absent" ]; then
    echo "WARNING: $2 names files that are not on disk:${_rr_absent}" >&2
    echo "         Those runs are classified without them -- a paired run then yields one mate." >&2
  fi
  awk '!/^#/ && NF >= 2 { n[$1]++ }
       END { for (k in n) if (n[k] > most) most = n[k]
             for (k in n) if (n[k] < most)
               printf "WARNING: %s resolves to %d file(s) where others have %d.\n", k, n[k], most > "/dev/stderr" }' "$map"

  echo "############ $1: real reads (${3}), unrefined vs. refined ############"
  mvn exec:exec@specificity -Dname="$1" -Dfqmap="$2" -Dreportkey="$3" -Dgs.project.calibration="${4:-}"
}

# Wall time and maximum RAM of classifying the *real* reads, measured the same way as for the
# simulated ones and for the same reason kept apart from the quality run above -- see the comment on
# run_perf(). The log key distinguishes these from the simulated runs of the same database.
#
# $1 = database project name, $2 = fastq mapping file, $3 = report key
run_real_perf() {
  run_perf "$1" "$2" "$3"
}

run_ticks() {
  map="${basedir}/data/fastq/ticks_sim.txt"
  if [ ! -f "$map" ]; then
    echo "Missing ${map} - run 'sh ./bin/make_fastqs.sh tick-borne' first." >&2
    return 1
  fi
  echo "############ tick-borne: NanoSim reads, unrefined vs. refined ############"
  mvn exec:exec@accuracy -Dname=tick-borne -Dfqmap=ticks_sim.txt -Dreportkey=nanosim -Dsimulator=NANOSIM
}

# Wall time and maximum RAM of classifying the same fastq files, once against the unrefined and
# once against the refined database. This is measured separately from the accuracy runs above and
# not derived from them, for two reasons: those run both variants inside a single JVM, so their
# cost cannot be separated, and they score every read against its ground truth on top, which is
# evaluation work rather than classification work.
#
# Each variant therefore runs on its own through the plain goals `match' and `ftmatch', in its own
# JVM, under cgmemtime. The goal is cleaned first: Genestrip skips a goal whose result files are
# already in place, which would otherwise be measured as a runtime of nearly zero.
#
# $1 = database project name, $2 = fastq mapping file, $3 = key for the log file name (optional,
# defaults to the database name). The key exists because a database is timed on more than one read
# collection -- the simulated one and the real one -- and the logs must not overwrite each other.
run_perf() {
  db=$1
  map=$2
  logkey=${3:-$db}
  # The error-free reads exist to bound what the refinement can achieve, not to time it: they are
  # the same volume of data through the same code path, so timing them again adds nothing.
  if [ -n "${ERROR_FREE:-}" ]; then
    echo "SKIP  classification performance for ${db}: measured on the regular reads only."
    return 0
  fi
  if [ ! -f "${basedir}/data/fastq/${map}" ]; then
    echo "Missing ${basedir}/data/fastq/${map} - run 'sh ./bin/make_fastqs.sh ${db}' first." >&2
    return 1
  fi
  # A warning rather than an error in both cases below: the quality results above are complete and
  # worth keeping, and only the performance part of Section "Performance" would be missing.
  if [ ! -x ./tools/cgmemtime/cgmemtime ]; then
    echo "WARNING: cgmemtime is missing - skipping the classification performance of ${db}." >&2
    echo "         Run ./bin/install_tools.sh to enable it." >&2
    return 0
  fi
  # Being present is not enough: cgmemtime needs a cgroup it may create, which it cannot do without
  # a systemd user session. Probing it once here beats discovering it after the measured run, where
  # the failure would be indistinguishable from a failure of the goal itself.
  if ! ./tools/cgmemtime/cgmemtime true >/dev/null 2>&1; then
    echo "WARNING: cgmemtime is installed but cannot run here - skipping the classification" >&2
    echo "         performance of ${db}. It needs a systemd user scope to create its cgroup." >&2
    return 0
  fi
  # goal name -> log file prefix; `match' uses the unrefined database, `ftmatch' the refined one.
  for goal in match ftmatch; do
    echo "############ ${logkey}: classification performance, goal ${goal} ############"
    mvn exec:exec@match -Dname="$db" -Dgoal="$goal" -Dfqmap="$map" -Dgs.target=clean
    ./tools/cgmemtime/cgmemtime mvn exec:exec@match -Dname="$db" -Dgoal="$goal" -Dfqmap="$map" \
        > "${res_path}/${goal}_${logkey}.log"
    publish_match_results "$db" "$goal" "$logkey"
  done
}

# Copies the per-taxon match results of the run just measured into <results>.
#
# These are the goal's own output -- one CSV per fastq key, giving per taxon how many reads were
# assigned to it -- and they are the only per-taxon view of a run that survives it. The reports
# beside them aggregate: <db>_<key>_summary.csv states how large the genus-only subset |R'_g| is,
# these say *which* taxa it consists of, which is what an analysis of why that subset differs
# between a simulation and the real sample it models has to work from.
#
# They must be copied here, inside the loop, rather than by a sweep at the end. Genestrip names them
# after the fastq key, and the simulated and the real tick runs use the same keys tick1 .. tick8 from
# two different maps (ticks_sim.txt and eightticks.txt), so the second run overwrites the first's
# files in the project folder -- and the `clean' target above deletes them outright. The log key,
# which already keeps the two runs' performance logs apart, keeps their CSVs apart the same way.
#
# Note that this rides on the performance measurement: run_perf returns early when cgmemtime is
# unavailable, and then the match goals never run and there is nothing to copy.
#
# $1 = database project name, $2 = goal name, $3 = log key distinguishing this run
publish_match_results() {
  _pm_db=$1; _pm_goal=$2; _pm_logkey=$3
  _pm_dir="${basedir}/data/projects/${_pm_db}/csv"
  _pm_n=0
  for _pm_src in "${_pm_dir}/${_pm_db}_${_pm_goal}_"*.csv; do
    [ -e "$_pm_src" ] || continue
    _pm_key=$(basename "$_pm_src" .csv)
    _pm_key=${_pm_key#"${_pm_db}_${_pm_goal}_"}
    cp "$_pm_src" "${res_path}/${_pm_goal}_${_pm_logkey}_${_pm_key}.csv"
    _pm_n=$((_pm_n + 1))
  done
  if [ "$_pm_n" -gt 0 ]; then
    echo "  kept ${_pm_n} per-taxon result(s) as ${res_path}/${_pm_goal}_${_pm_logkey}_*.csv"
  else
    echo "  WARNING: no ${_pm_db}_${_pm_goal}_*.csv in ${_pm_dir} - no per-taxon results kept." >&2
  fi
}

# The first of the two runs reads its database from a cold page cache while the second may find
# parts of the file system cache still warm. Both databases are of nearly the same size, so the
# effect is comparable, but for a close comparison it is worth repeating the pair in the opposite
# order and keeping the slower figure of each variant.

case "$what" in
  viral)        run_iss viral; run_perf viral viral_sim.txt ;;
  protozoa)     run_iss protozoa; run_perf protozoa protozoa_sim.txt ;;
  strepto)      run_iss strepto; run_perf strepto strepto_sim.txt ;;
  # `nocardia' has real reads but they are few and very unevenly spread -- 29,175 Nocardia reads over
  # fourteen samples, of which two hold 82 per cent and one holds two. The simulated set is what
  # makes a calibration possible at all; the real one is read case by case, not averaged.
  nocardia)     run_iss nocardia; run_perf nocardia nocardia_sim.txt ;;
  tick-borne)   run_ticks; run_perf tick-borne ticks_sim.txt ;;
  # The saliva-matched simulations, one per database the real saliva runs are classified against.
  # Accuracy only: this set exists to calibrate those runs, its own performance is not reported, and
  # a per-database target would drag the plain simulated runs along with it. ERROR_SALIVA is set and
  # restored here rather than expected in the environment, for the reason run_iss_all_regimes gives:
  # one left standing sends a later call to the wrong read set without saying so.
  saliva)       _cl_saved=${ERROR_SALIVA:-}
                ERROR_SALIVA=1
                run_iss viral; run_iss strepto
                ERROR_SALIVA=$_cl_saved ;;
  # Kraken 2 and KrakenUniq on the same four simulated sets of `viral', which is the one database
  # another tool can be given the same scope. Their per-read output lands in results/kraken, and the
  # next `accuracy' run for `viral' picks it up by itself: RefinementAccuracyReport scores whatever
  # it finds there on the very reads of the unrefined Genestrip run. So the order matters -- build
  # and classify first, evaluate after.
  kraken)       sh ./bin/kraken_build.sh all
                sh ./bin/kraken_classify.sh all ;;
  # The k = 24 twin of `cv', on the same four read sets and with the same measures. It is a control
  # for what governs the refinement's gain: a shorter k-mer is shared by more species, so more of
  # the database's k-mers end up above the data taxa, which is the mass a refinement can move. The
  # twin is created by bin/k24_projects.sh, which also links the read maps under its name, so the
  # reads here are the very files `viral' was scored on. ERROR_SALIVA is set and restored around the
  # saliva-matched set for the reason given above.
  k24)          _cl_db=${K24_DB:-viral-k24}
                _cl_saved=${ERROR_SALIVA:-}
                ERROR_SALIVA=""
                run_iss_all_regimes "$_cl_db"
                ERROR_SALIVA=1
                run_iss "$_cl_db"
                ERROR_SALIVA=$_cl_saved
                run_perf "$_cl_db" "${_cl_db}_sim.txt" ;;
  # The same for `nocardia', whose real runs are BGISEQ-500 rather than saliva: one simulated set at
  # their length and error, so that the estimate columns of the real-read table have a factor.
  mngs)         _cl_saved=${ERROR_MNGS:-}
                ERROR_MNGS=1
                run_iss nocardia
                ERROR_MNGS=$_cl_saved ;;
  accuracy)     run_iss_all_regimes viral; run_iss_all_regimes protozoa; run_iss_all_regimes strepto
                run_iss_all_regimes nocardia; run_ticks ;;
  # `strepto' runs on the saliva files. Saliva is a streptococcal habitat -- the oral species of
  # the mitis and salivarius groups are what a healthy mouth carries -- so the same five runs put a
  # bacterial database of one genus beside the viral one on identical material, with the calibration
  # taken from the saliva-matched simulation of that database. Its clinical runs are downloaded and
  # still train the Nanopore simulation, but they are no longer classified for the paper.
  real)         run_real viral "${SALIVA_MAP:-saliva_real.txt}" saliva iss_saliva
                run_real strepto "${SALIVA_MAP:-saliva_real.txt}" saliva iss_saliva
                run_real tick-borne eightticks.txt ticks nanosim
                run_real nocardia nocardia_mngs.txt mngs iss_mngs
                run_real_perf viral "${SALIVA_MAP:-saliva_real.txt}" saliva
                run_real_perf strepto "${SALIVA_MAP:-saliva_real.txt}" saliva
                run_real_perf tick-borne eightticks.txt ticks ;;
  perf)         run_perf viral viral_sim.txt; run_perf protozoa protozoa_sim.txt
                run_perf strepto strepto_sim.txt; run_perf nocardia nocardia_sim.txt
                run_perf tick-borne ticks_sim.txt
                run_real_perf viral "${SALIVA_MAP:-saliva_real.txt}" saliva
                run_real_perf strepto "${SALIVA_MAP:-saliva_real.txt}" saliva
                run_real_perf tick-borne eightticks.txt ticks ;;
  all)          run_iss_all_regimes viral; run_iss_all_regimes protozoa; run_iss_all_regimes strepto
                run_iss_all_regimes nocardia; run_ticks
                run_perf viral viral_sim.txt; run_perf protozoa protozoa_sim.txt
                run_perf strepto strepto_sim.txt; run_perf nocardia nocardia_sim.txt
                run_perf tick-borne ticks_sim.txt
                run_real viral "${SALIVA_MAP:-saliva_real.txt}" saliva iss_saliva
                run_real strepto "${SALIVA_MAP:-saliva_real.txt}" saliva iss_saliva
                run_real tick-borne eightticks.txt ticks nanosim
                run_real nocardia nocardia_mngs.txt mngs iss_mngs ;;
  *)          echo "Usage: $0 [viral|protozoa|strepto|nocardia|tick-borne|saliva|mngs|kraken|k24|accuracy|perf|real|all]" >&2; exit 1 ;;
esac

echo
echo "=== results ==="
ls -la "$res_path"/*_accuracy.csv 2>/dev/null || echo "no accuracy results yet"
ls -la "$res_path"/match_*.log "$res_path"/ftmatch_*.log 2>/dev/null || echo "no performance logs yet"
ls -la "$res_path"/match_*_*.csv "$res_path"/ftmatch_*_*.csv 2>/dev/null \
    || echo "no per-taxon match results yet"
