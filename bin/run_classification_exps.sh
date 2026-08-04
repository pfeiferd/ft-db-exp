#!/bin/sh
#
# Runs the classification quality experiments: every simulated fastq file is matched against the
# unrefined and against the refined database, and the resulting classifications are scored against
# the reads' known ground truth.
#
# The measures are the ones of Section "Classification quality" in the Genestrip-FT paper: the
# boolean positive counts at the genus and the species rank, plus the candidate-weighted species
# count, which credits a classification in proportion to how far it narrows the species down. The
# latter is the one that can tell a refined database from an unrefined one at all.
#
# Prerequisites:
#   sh ./bin/install_tools.sh     installs InSilicoSeq and NanoSim
#   sh ./bin/run_exps.sh          builds the databases (or at least the goals `db` and `ftdb`)
#   sh ./bin/make_fastqs.sh       generates the simulated reads
#
# Results land in ./results as <db>_<report key>_accuracy.csv and are meant to be copied into the
# paper's own results folder, from where they are included directly.
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

run_ticks() {
  map="${basedir}/data/fastq/ticks_sim.txt"
  if [ ! -f "$map" ]; then
    echo "Missing ${map} - run 'sh ./bin/make_fastqs.sh tick-borne' first." >&2
    return 1
  fi
  echo "############ tick-borne: NanoSim reads, unrefined vs. refined ############"
  mvn exec:exec@accuracy -Dname=tick-borne -Dfqmap=ticks_sim.txt -Dreportkey=nanosim -Dsimulator=NANOSIM
}

case "$what" in
  viral)      run_iss viral ;;
  protozoa)   run_iss protozoa ;;
  gut-protozoa) run_iss gut-protozoa ;;
  tick-borne) run_ticks ;;
  all)        run_iss viral; run_iss protozoa; run_ticks ;;
  *)          echo "Usage: $0 [viral|protozoa|gut-protozoa|tick-borne|all]" >&2; exit 1 ;;
esac

echo
echo "=== results ==="
ls -la "$res_path"/*_accuracy.csv 2>/dev/null || echo "no accuracy results yet"
