#!/bin/sh
#
# Runs every experiment behind the Genestrip-FT paper, end to end, in the one order that satisfies
# the dependencies between them. It drives the existing scripts rather than replacing them: each
# step below is a command that can also be run on its own, and this file's contribution is the order
# and the handful of preconditions that are easy to get wrong when the steps are run by hand.
#
#   sh ./bin/run_all_exps.sh              # everything, from the build to the paper's results folder
#   sh ./bin/run_all_exps.sh --from 6     # resume at step 6 after a failure
#   sh ./bin/run_all_exps.sh --only 7     # a single step
#   sh ./bin/run_all_exps.sh --list       # what the steps are, without running any
#   sh ./bin/run_all_exps.sh --dry-run    # print the commands instead of running them
#
# Re-running is safe. Every underlying step skips work it has already done -- Genestrip treats a
# goal whose result files exist as made, make_fastqs.sh skips a fastq file that is there, and
# run_exps.sh refuses to re-time a database it would not actually rebuild -- so a run interrupted
# after two days can simply be started again.
#
# Environment:
#   FRESH_DBS=1     delete the databases first, via clean_all.sh, so that step 4 measures their
#                   generation. Without it an existing database is kept and its db_gen_*.log and
#                   ftdb_gen_*.log are left as they are: run_exps.sh will not overwrite a real
#                   measurement with the timing of a goal that did nothing. Set this whenever the
#                   paper's performance table is to be recomputed from one consistent batch.
#   FETCH_SALIVA=1  fetch the human saliva runs (about 400 GB) before the real-read step. Off by
#                   default, since they are usually already on the machine or deliberately absent.
#   KEEP_TICK_SIMS=1  never delete a simulated tick fastq -- see step 5, which by default deletes
#                   exactly those whose NanoSim abundance table was not preserved.
#   PAPER_RESULTS   a folder to copy the results to when the run is done, in addition to leaving
#                   them in ./results. Unset by default: this project stands on its own and knows
#                   nothing about where a paper or any other consumer keeps its inputs. The folder
#                   must exist; step 14 copies into it and regenerates its LaTeX macros there.
#   SKIP_BUILD=1    do not run `mvn install' in step 1.
#
# Wall time is days, not hours, and the disk needs well over a terabyte -- see README.md.
#
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)

# Empty unless asked for -- see PAPER_RESULTS above. No default path here: a consumer's layout is
# not this project's to guess, and a guess that misses produces a warning after a run of days.
paper_results=${PAPER_RESULTS:-}
paper_results_said=""

# States once what becomes of the results at the end of the run. The preflight calls it so that a
# PAPER_RESULTS pointing nowhere is known in the first minute rather than after the last step, and
# step 14 calls it again -- where it stays silent, having already said its piece. Without that the
# same warning is printed once per invocation, three times over a run resumed twice.
check_paper_results() {
  if [ -n "$paper_results_said" ]; then
    return 0
  fi
  paper_results_said=1
  if [ -z "$paper_results" ]; then
    echo "OK    results stay in ${basedir}/results (PAPER_RESULTS is unset)"
  elif [ -d "$paper_results" ]; then
    echo "OK    step 14 copies them to ${paper_results}"
  else
    echo "WARNING: PAPER_RESULTS=${paper_results} does not exist - the results will stay in" >&2
    echo "         ${basedir}/results. Create the folder, or copy them by hand afterwards." >&2
  fi
}
# The samples make_fastqs.sh simulates, mirrored here for step 5. Kept in step with its own default.
tick_samples=${SAMPLES:-"tick1 tick2 tick3 tick4 tick5 tick6 tick7 tick8"}

from=1
only=""
dry=""
while [ $# -gt 0 ]; do
  case "$1" in
    --from) from=$2; shift 2 ;;
    --only) only=$2; from=$2; shift 2 ;;
    --dry-run) dry=1; shift ;;
    --list) list=1; shift ;;
    -h|--help) sed -n '2,40p' "$0"; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

step_no=0
started=$(date '+%Y-%m-%d %H:%M:%S')

# Announces a step and decides whether it runs. Returns 1 when the step is to be skipped, so each
# step below reads `step "..." || true' and guards its body with a plain if.
step() {
  step_no=$((step_no + 1))
  if [ -n "${list:-}" ]; then
    echo "  ${step_no}. $1"
    return 1
  fi
  if [ "$step_no" -lt "$from" ] || { [ -n "$only" ] && [ "$step_no" != "$only" ]; }; then
    echo "SKIP  step ${step_no}: $1"
    return 1
  fi
  echo
  echo "==================================================================================="
  echo "=== step ${step_no}: $1"
  echo "=== $(date '+%Y-%m-%d %H:%M:%S')  (run started ${started})"
  echo "==================================================================================="
  return 0
}

# Runs a command, or prints it under --dry-run.
run() {
  if [ -n "$dry" ]; then
    echo "    would run: $*"
  else
    "$@"
  fi
}

[ -n "${list:-}" ] && echo "Steps:"

########################################################################################
step "build the Java" || skip=1
if [ -z "${skip:-}" ]; then
  # First, and not optional in practice. The evaluation reports -- the labels their CSV files carry
  # and the columns they are made of -- are Java in this project, so a run against a stale build
  # reproduces whatever those classes did last time. That is not hypothetical: the saliva read set
  # was printed under its raw fastq key until RefinementAccuracyReport learned its name, and every
  # paper table including that CSV failed to typeset.
  if [ -z "${SKIP_BUILD:-}" ]; then
    run mvn install
  else
    echo "SKIP_BUILD is set - using whatever is in target/."
  fi
fi
skip=""

########################################################################################
step "preflight" || skip=1
if [ -z "${skip:-}" ]; then
  # cgmemtime is checked here rather than discovered in step 4 or 10, because it gates more than the
  # performance figures: run_perf returns early without it, and the per-taxon match results that
  # step 10 preserves are produced by the goals it measures. No cgmemtime, no per-taxon CSVs.
  if [ ! -x ./tools/cgmemtime/cgmemtime ]; then
    echo "WARNING: ./tools/cgmemtime/cgmemtime is missing - run ./bin/install_tools.sh." >&2
    echo "         Without it the performance figures and the per-taxon match results are lost." >&2
  elif ! ./tools/cgmemtime/cgmemtime true >/dev/null 2>&1; then
    echo "WARNING: cgmemtime is installed but cannot run here -- it needs a systemd user scope to" >&2
    echo "         create its cgroup. The performance figures and the per-taxon match results" >&2
    echo "         will be missing from this run." >&2
  else
    echo "OK    cgmemtime"
  fi
  # The real-read step fails on a missing map rather than skipping it, so say so now instead of
  # after however many hours step 10 takes to get there.
  for map in "${SALIVA_MAP:-saliva_real.txt}" eightticks.txt nocardia_mngs.txt; do
    if [ -f "data/fastq/${map}" ]; then
      echo "OK    data/fastq/${map}"
    else
      echo "WARNING: data/fastq/${map} is missing - the real-read runs of step 10 will fail." >&2
    fi
  done
  echo "OK    results go to ${basedir}/results"
  check_paper_results
fi
skip=""

########################################################################################
step "delete the databases (FRESH_DBS)" || skip=1
if [ -z "${skip:-}" ]; then
  # Only on request, and this is the one destructive step in the script. run_exps.sh will not
  # measure the generation of a database that already exists, because re-running the goal would do
  # nothing and the timing would capture Maven's startup -- so it keeps the old log instead. That is
  # the right default for a re-run, and the wrong one when the performance table is to be rebuilt
  # from a single consistent batch, which is what FRESH_DBS is for. The RefSeq downloads are not
  # touched; only the generated databases are.
  if [ -n "${FRESH_DBS:-}" ]; then
    run sh ./bin/clean_all.sh
  else
    echo "FRESH_DBS is not set - keeping the existing databases and their generation logs."
    echo "The intrinsic quality CSVs are recomputed either way; only the timings are not."
  fi
fi
skip=""

########################################################################################
step "databases, intrinsic quality, tree figures" || skip=1
if [ -z "${skip:-}" ]; then
  # Has to precede the fastq generation: the InSilicoSeq reads are drawn from the goal
  # `extractrefseqcsv' and NanoSim's genome list from `extractrefseqfasta', and both need a built
  # database. Also collects every CSV, SVG and LaTeX file under data/ into results/ and runs
  # paper_stats.sh over them.
  run sh ./bin/run_exps.sh
fi
skip=""

########################################################################################
step "clear tick simulations whose abundance was not preserved" || skip=1
if [ -z "${skip:-}" ]; then
  # NanoSim's abundance table is what says how the simulated composition relates to the real sample
  # it was trained on, and make_fastqs.sh preserves it only when it actually simulates. A sample
  # whose _sim.fastq is already there is skipped, and no amount of re-running produces the table for
  # it -- the simulation has to happen again.
  #
  # So exactly those samples are cleared, and only those: a sample whose abundance is already beside
  # its reads is left alone. That keeps this from being a blanket "delete the simulations" step. It
  # costs a NanoSim training run per deleted sample, which is hours, and after one complete run it
  # is a no-op. KEEP_TICK_SIMS=1 suppresses it entirely, at the price of an incomplete analysis.
  if [ -n "${KEEP_TICK_SIMS:-}" ]; then
    echo "KEEP_TICK_SIMS is set - leaving every simulated tick fastq in place."
  else
    cleared=""
    for s in $tick_samples; do
      if [ -s "data/fastq/${s}_sim.fastq" ] && [ ! -f "data/fastq/${s}_sim_quantification.tsv" ]; then
        run rm -f "data/fastq/${s}_sim.fastq"
        cleared="${cleared} ${s}"
      fi
    done
    if [ -n "$cleared" ]; then
      echo "Cleared for re-simulation (no abundance table):${cleared}"
    else
      echo "Nothing to clear - every simulated tick sample has its abundance table."
    fi
  fi
fi
skip=""

########################################################################################
step "simulated reads: the reported regimes" || skip=1
if [ -z "${skip:-}" ]; then
  # The two Illumina models and the error-free set for every InSilicoSeq project, plus the eight
  # NanoSim tick simulations and the single NanoSim set of `strepto', which calibrates its 83
  # nanopore runs wholesale rather than one by one. ERROR_SALIVA is explicitly cleared: it is not one of the reported
  # regimes, and make_iss picks the first regime variable that is non-empty, so one inherited from
  # the caller's environment would silently redirect this call.
  run env ERROR_SALIVA= sh ./bin/make_fastqs.sh all
fi
skip=""

########################################################################################
step "simulated reads: the calibration sets" || skip=1
if [ -z "${skip:-}" ]; then
  # Its own invocation, and it has to be. This set is not one of the three the paper reports side by
  # side; it exists so that the calibration factors of the real saliva runs can be measured at the
  # parameters of the data they are applied to -- 101 bp at 2.07 % per-base error, which no stock
  # InSilicoSeq model reproduces.
  #
  # For both databases the saliva runs are classified against: `cv' and `strepto'. A calibration is
  # a property of the pair (database, read set), so each needs its own -- `strepto' is not calibrated
  # by viral reads.
  run env ERROR_SALIVA=1 sh ./bin/make_fastqs.sh viral
  run env ERROR_SALIVA=1 sh ./bin/make_fastqs.sh strepto
  # The same for `nocardia' and its BGISEQ runs, at 70 bp and 2.05 % rather than 101 bp and 2.07 %.
  run env ERROR_MNGS=1 sh ./bin/make_fastqs.sh nocardia
fi
skip=""

########################################################################################
step "fetch the human saliva runs (FETCH_SALIVA)" || skip=1
if [ -z "${skip:-}" ]; then
  # About 400 GB of gzipped fastq for the three runs, and transiently as much again for the .sra
  # files and fasterq-dump's scratch space -- hence off unless asked for. Anything already fetched
  # is reused rather than downloaded again.
  if [ -n "${FETCH_SALIVA:-}" ]; then
    run sh ./bin/fetch_saliva.sh
  else
    echo "FETCH_SALIVA is not set - assuming the saliva runs are already in data/fastq."
  fi
fi
skip=""

########################################################################################
step "classification: the calibration sets" || skip=1
if [ -z "${skip:-}" ]; then
  # Before step 10, not after. This run writes results/<db>_iss_saliva_summary.csv, and the real
  # saliva runs of step 10 read exactly that file for their calibration -- `run_real <db> ... saliva
  # iss_saliva'. Run in the other order, the estimate columns of the real-read table come from
  # whatever stale copy happened to be lying in results/, or from nothing at all.
  run sh ./bin/run_classification_exps.sh saliva
  # And the mNGS-matched set of `nocardia', for the same reason and before the same step.
  run sh ./bin/run_classification_exps.sh mngs
fi
skip=""

########################################################################################
step "classification: everything else, simulated and real" || skip=1
if [ -z "${skip:-}" ]; then
  # Accuracy for every simulated read set, the performance measurements, and the real-read runs for
  # both collections. This is also where the per-taxon match results are preserved, one CSV per
  # fastq key per database variant, keyed by the run they belong to.
  run env ERROR_SALIVA= sh ./bin/run_classification_exps.sh all
fi
skip=""

########################################################################################
step "classification performance: the four scenarios" || skip=1
if [ -z "${skip:-}" ]; then
  # Table "matchperf" of the paper reads results/matchperf.csv, and nothing else writes it: the
  # measurement lives in perf_scenarios.sh, which until now had to be remembered by hand -- so the
  # table showed whatever the last manual run had left, or its markers.
  #
  # After the classification steps, because it deletes the `match' and `ftmatch' results of the keys
  # it measures: a goal whose output is still there is not remade, and the run then times a JVM
  # start instead of a classification. Step 10 has copied what it needs into results/ by then.
  run sh ./bin/perf_scenarios.sh all
fi
skip=""

########################################################################################
step "record the machine and regenerate the macros" || skip=1
if [ -z "${skip:-}" ]; then
  # sysinfo.sh writes results/sysinfo.txt and results/sysinfo.tex: the hardware and JVM facts a
  # paper states about the machine its measurements were taken on. Here rather than in the preflight
  # so that it describes the run that has just happened. Run under sudo it also records DMI and
  # SMART detail; nothing in the LaTeX file needs it.
  run sh ./bin/sysinfo.sh -q
  # paper_stats.sh runs a second time here -- run_exps.sh already ran it at the end of step 4, which
  # is before any classification result exists, so the dbstats.tex it left describes half a run.
  # Regenerating it now means ./results is complete and consistent on its own, which is what a
  # results folder carried off this machine by hand consists of. Step 14 runs it once more against
  # PAPER_RESULTS, where the folder may hold files from other batches as well.
  run sh ./bin/paper_stats.sh
fi
skip=""

########################################################################################
step "copy the results to PAPER_RESULTS" || skip=1
if [ -z "${skip:-}" ]; then
  # A consumer that reads the CSV files from a folder of its own sees nothing until they are copied
  # there. run_exps.sh runs paper_stats.sh at the end of step 4, i.e. before any classification
  # result exists, which is why the copy happens here and the macros are regenerated afterwards from
  # what the target folder actually holds.
  if [ -z "$paper_results" ] || [ ! -d "$paper_results" ]; then
    check_paper_results
    echo "Nothing copied - the results are in ${basedir}/results."
  else
    # -R because the results hold the `logs' folder as well, which a plain cp would refuse to copy
    # and so would fail the step.
    run sh -c "cp -R \"${basedir}\"/results/* \"${paper_results}\"/"
    run sh ./bin/paper_stats.sh "$paper_results"
    echo "Copied to ${paper_results} and regenerated its LaTeX macros."
  fi
fi
skip=""

if [ -n "${list:-}" ]; then
  exit 0
fi

echo
echo "==================================================================================="
echo "=== done.  started ${started}, finished $(date '+%Y-%m-%d %H:%M:%S')"
echo "==================================================================================="
ls -la "${basedir}/results" | tail -n +2 | wc -l | sed 's/^/files in results: /'
