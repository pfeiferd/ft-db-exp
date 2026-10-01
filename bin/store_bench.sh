#!/bin/sh
#
# The k-mer store comparison of the paper's appendix: the same database, once in the radix store it
# was built into and once in a sorted array copied from it, measured for memory and for lookup
# throughput. Everything runs through this project's dependencies, so no Genestrip checkout is
# needed.
#
# Usage:
#   sh ./bin/store_bench.sh                            # viral, the smallest fastq of the project
#   sh ./bin/store_bench.sh tick-borne                 # another project
#   sh ./bin/store_bench.sh viral <file.fastq.gz>      # an explicit input
#
# DB picks the database file, by default the project's refined one. GS_XMX caps the heap of the
# measured JVM; both stores are held at once, so it needs about 18 bytes per k-mer of the database
# plus its filters, e.g. 12G for `viral'.
#
# Writes results/storebench.csv with one row per store and results/storebench_<project>.log.
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

if [ -z "$input" ]; then
  # The smallest file the project offers, so that a first run is the cheapest one. The project's own
  # fastq folder first, then the shared one, which is where the downloads land.
  input=$(ls -S "${data}/projects/${project}/fastq"/*.fastq.gz 2>/dev/null | tail -1)
  if [ -z "$input" ]; then
    input=$(ls -S "${data}/fastq"/*.fastq.gz 2>/dev/null | tail -1)
  fi
  [ -n "$input" ] || { echo "No fastq.gz for ${project}; name one explicitly." >&2; exit 1; }
fi
case "$input" in
  /*) ;;
  *) input="${data}/projects/${project}/fastq/${input}" ;;
esac
[ -f "$input" ] || { echo "No input at ${input}." >&2; exit 1; }

mkdir -p "$res"
echo "############ store comparison: ${project} ############"
echo "database: ${db}"
echo "input:    ${input}"
echo "log:      ${log}"

xmx_opt=""
if [ -n "${GS_XMX:-}" ]; then
  xmx_opt="-Dgs.xmx=${GS_XMX}"
fi

mvn exec:exec@storebench -Dname="$project" -Dfqfile="$input" -Dgs.store.db="$db" $xmx_opt 2>&1 | tee "$log"
