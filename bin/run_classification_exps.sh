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
#   sh ./bin/run_classification_exps.sh [viral|protozoa|gut-protozoa|tick-borne|accuracy|perf|all]
#
#   A database name runs both parts for it, following ERROR_FREE; `accuracy' and `perf' run one
#   part for every database, and `accuracy' and `all' evaluate the error-containing *and* the
#   error-free reads of every InSilicoSeq project.
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

# Both evaluations of a database: on the error-containing reads and on the error-free ones. See
# make_fastqs.sh for why a comprehensive run needs both.
run_iss_both() {
  saved_error_free=${ERROR_FREE:-}
  ERROR_FREE=""
  run_iss "$1"
  ERROR_FREE=1
  run_iss "$1"
  ERROR_FREE=$saved_error_free
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
# $1 = database project name, $2 = fastq mapping file
run_perf() {
  db=$1
  map=$2
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
  if [ ! -x ./tools/cgmemtime/cgmemtime ]; then
    # A warning rather than an error: the quality results above are complete and worth keeping,
    # and only the performance part of Section "Performance" is missing.
    echo "WARNING: cgmemtime is missing - skipping the classification performance of ${db}." >&2
    echo "         Run ./bin/install_tools.sh to enable it." >&2
    return 0
  fi
  # goal name -> log file prefix; `match' uses the unrefined database, `ftmatch' the refined one.
  for goal in match ftmatch; do
    echo "############ ${db}: classification performance, goal ${goal} ############"
    mvn exec:exec@match -Dname="$db" -Dgoal="$goal" -Dfqmap="$map" -Dgs.target=clean
    ./tools/cgmemtime/cgmemtime mvn exec:exec@match -Dname="$db" -Dgoal="$goal" -Dfqmap="$map" \
        > "${res_path}/${goal}_${db}.log"
  done
}

# The first of the two runs reads its database from a cold page cache while the second may find
# parts of the file system cache still warm. Both databases are of nearly the same size, so the
# effect is comparable, but for a close comparison it is worth repeating the pair in the opposite
# order and keeping the slower figure of each variant.

case "$what" in
  viral)        run_iss viral; run_perf viral viral_sim.txt ;;
  protozoa)     run_iss protozoa; run_perf protozoa protozoa_sim.txt ;;
  gut-protozoa) run_iss gut-protozoa; run_perf gut-protozoa gut-protozoa_sim.txt ;;
  tick-borne)   run_ticks; run_perf tick-borne ticks_sim.txt ;;
  accuracy)     run_iss_both viral; run_iss_both protozoa; run_iss_both gut-protozoa; run_ticks ;;
  perf)         run_perf viral viral_sim.txt; run_perf protozoa protozoa_sim.txt
                run_perf gut-protozoa gut-protozoa_sim.txt; run_perf tick-borne ticks_sim.txt ;;
  all)          run_iss_both viral; run_iss_both protozoa; run_iss_both gut-protozoa; run_ticks
                run_perf viral viral_sim.txt; run_perf protozoa protozoa_sim.txt
                run_perf gut-protozoa gut-protozoa_sim.txt; run_perf tick-borne ticks_sim.txt ;;
  *)          echo "Usage: $0 [viral|protozoa|gut-protozoa|tick-borne|accuracy|perf|all]" >&2; exit 1 ;;
esac

echo
echo "=== results ==="
ls -la "$res_path"/*_accuracy.csv 2>/dev/null || echo "no accuracy results yet"
ls -la "$res_path"/match_*.log "$res_path"/ftmatch_*.log 2>/dev/null || echo "no performance logs yet"
