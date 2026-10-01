#!/bin/sh
#
# The k-mer store comparison of the paper's appendix: the same database, once in the radix store it
# was built into and once in a sorted array copied from it, measured for memory and for lookup
# throughput. Everything runs through this project's dependencies, so no Genestrip checkout is
# needed.
#
# Usage:
#   ALL=1 sh ./bin/store_bench.sh strepto              # every simulated saliva-like file
#   ALL=1 sh ./bin/store_bench.sh tick-borne           # every simulated tick file
#   sh ./bin/store_bench.sh viral <file.fastq.gz>      # an explicit input, resolved against data/fastq
#   REAL=1 ALL=1 sh ./bin/store_bench.sh tick-borne    # the real tick runs instead
#   REAL=1 ALL=1 FILTER=both sh ./bin/store_bench.sh tick-borne   # real runs, with and without it
#   THREADS=-1 sh ./bin/store_bench.sh strepto         # one consumer per processor less one
#   THREADS=2 sh ./bin/store_bench.sh strepto          # two consumers, three threads in all
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
# THREADS sets the CONSUMER threads and is one by default, so that a run has two threads in all:
# the producer that reads, inflates and parses the file, and one consumer that looks the k-mers up.
# That is what the paper's appendix uses throughout. With one consumer per processor the producer is
# the limit of the run and a faster lookup hardly shows, while one consumer puts the limit on the
# lookups. It also makes a run several times longer than the pipeline default would, which is the
# price of measuring the store rather than the reader.
#
# Each file is classified twice per store and only the second pass counts. The first warms the JIT
# and the page cache, and it does so for either store alike: without it the store that is measured
# first reads a cold file and compiles the matching loop, while the one measured second gets both for
# free. WARMUP=0 turns that off and halves the time a run takes.
#
# FILTER takes `on' (the default), `off' or `both'. The Bloom filter in front of a store answers most
# k-mers of a real read set, which then never reach the store, so a run with it on measures the
# pipeline and one with it off measures the layout. `both' does the two in one go and so also says
# what the filter itself is worth, at twice the time.
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

# Without an explicit input the files come from the project's simulated read sets. Those are drawn
# from the genomes the database was built from, so nearly every k-mer of a read is in the store and
# the Bloom filter in front of it rejects almost nothing. That is what makes the comparison about the
# two layouts: on real data most k-mers never reach a store at all, and then the layout barely shows.
# REAL=1 takes the real read sets instead, MAP=<map> any other map. ALL=1 takes every file of the
# map, which is what the appendix uses for the simulated sets; by default the smallest one is taken.
if [ -z "$input" ]; then
  if [ -n "${REAL:-}" ]; then
    case "$project" in
      viral|strepto) map="${MAP:-saliva_real.txt}" ;;
      tick-borne) map="${MAP:-eightticks.txt}" ;;
      nocardia) map="${MAP:-nocardia_mngs.txt}" ;;
      *) map="${MAP:-}" ;;
    esac
  else
    case "$project" in
      strepto) map="${MAP:-strepto_sim_saliva.txt}" ;;
      tick-borne) map="${MAP:-ticks_sim.txt}" ;;
      viral) map="${MAP:-viral_sim.txt}" ;;
      nocardia) map="${MAP:-nocardia_sim_mngs.txt}" ;;
      *) map="${MAP:-}" ;;
    esac
  fi
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
filter_opt=""
if [ -n "${FILTER:-}" ]; then
  filter_opt="-Dgs.store.filter=${FILTER}"
fi
warmup_opt=""
if [ -n "${WARMUP:-}" ]; then
  case "$WARMUP" in
    0|no|false) warmup_opt="-Dgs.store.warmup=false" ;;
    *) warmup_opt="-Dgs.store.warmup=true" ;;
  esac
fi

mvn exec:exec@storebench -Dname="$project" -Dfqfile="$input" -Dgs.store.db="$db" $xmx_opt $threads_opt $warmup_opt $filter_opt 2>&1 | tee "$log"
