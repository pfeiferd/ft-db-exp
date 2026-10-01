#!/bin/sh
#
# The k-mer store comparison of the paper's appendix: the same database, once in the radix store it
# was built into and once in a sorted array copied from it, measured for memory and for lookup
# throughput. Everything runs through this project's dependencies, so no Genestrip checkout is
# needed.
#
# Usage:
#   sh ./bin/store_bench.sh                            # viral, the smallest file of its scenario map
#   sh ./bin/store_bench.sh strepto                    # the same reads against a denser database
#   sh ./bin/store_bench.sh viral <file.fastq.gz>      # an explicit input, resolved against data/fastq
#   ALL=1 sh ./bin/store_bench.sh strepto              # every file of the map, which takes hours
#   THREADS=-1 sh ./bin/store_bench.sh strepto         # one consumer per processor less one
#
# DB picks the database file, by default the project's refined one. MAP overrides the fastq map the
# input is taken from. GS_XMX caps the heap of the measured JVM; both stores are held at once,
# with one filter each, so it needs about 18 bytes per k-mer of the database plus some 2.5 bytes
# for the two filters. That is 10G for `viral' with its 386 million k-mers and 14G for `strepto'
# with 555 million; the default of 56G from the pom is plenty where the machine has it.
#
# `strepto' says more about the lookups than `viral' does: on the saliva runs it classifies seven
# to ten percent of the reads where `viral' classifies two, so more of the time is spent in the
# store. The memory figures do not care which database is measured.
#
# THREADS sets the consumer threads and is two by default, which is what the paper's appendix uses
# throughout. Reading, inflating and parsing a fastq file all happen on one producer thread, so with
# one consumer per processor that thread is the limit of the run and a faster lookup hardly shows.
# Two puts the limit back on the lookups. It also makes a run several times longer than the pipeline
# default would, which is the price of measuring the store rather than the reader.
#
# Writes results/storebench_<project>.csv with one row per store and results/storebench_<project>.log.
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)
data="${basedir}/data"
res="${basedir}/results"

project=${1:-viral}
input=$2
db=${DB:-${data}/projects/${project}/db/${project}_ftdb.zip}
log="${res}/storebench_${project}.log"

[ -f "$db" ] || { echo "No database at ${db}; build it or set DB=<file>." >&2; exit 1; }

# Without an explicit input the files come from the map that perf_scenarios.sh uses for the same
# database, so that the stores are compared on what the paper's scenarios are measured on. ALL=1
# takes every file of that map; by default the smallest one is taken, which is enough for a
# comparison and costs the least.
if [ -z "$input" ]; then
  case "$project" in
    viral|strepto) map="${MAP:-saliva_real.txt}" ;;
    tick-borne) map="${MAP:-eightticks.txt}" ;;
    nocardia) map="${MAP:-nocardia_mngs.txt}" ;;
    *) map="${MAP:-}" ;;
  esac
  [ -n "$map" ] || { echo "No default map for ${project}; name a file or set MAP=<map>." >&2; exit 1; }
  [ -f "${data}/fastq/${map}" ] || { echo "No map at ${data}/fastq/${map}." >&2; exit 1; }
  files=""
  for name in $(awk '!/^#/ && NF >= 2 {print $2}' "${data}/fastq/${map}"); do
    file="${data}/fastq/${name}"
    [ -f "$file" ] || continue
    files="${files}${files:+ }${file}"
  done
  [ -n "$files" ] || { echo "None of the files of ${map} is on disk." >&2; exit 1; }
  if [ -n "${ALL:-}" ]; then
    input=$(echo "$files" | tr ' ' ',')
  else
    # shellcheck disable=SC2086
    input=$(ls -S $files | tail -1)
  fi
else
  case "$input" in
    /*) ;;
    *) input="${data}/fastq/${input}" ;;
  esac
  [ -f "$input" ] || { echo "No input at ${input}." >&2; exit 1; }
fi

mkdir -p "$res"
echo "############ store comparison: ${project} ############"
echo "database: ${db}"
echo "input:    ${input}"
echo "log:      ${log}"

xmx_opt=""
if [ -n "${GS_XMX:-}" ]; then
  xmx_opt="-Dgs.xmx=${GS_XMX}"
fi
threads_opt=""
if [ -n "${THREADS:-}" ]; then
  threads_opt="-Dgs.store.threads=${THREADS}"
fi

mvn exec:exec@storebench -Dname="$project" -Dfqfile="$input" -Dgs.store.db="$db" $xmx_opt $threads_opt 2>&1 | tee "$log"
