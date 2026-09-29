#!/bin/sh
#
# Classifies the four simulated read sets of `viral' with Kraken 2 and KrakenUniq, so that their
# answers can be scored with the measures of the paper and put beside Genestrip's own.
#
# The read sets are the ones Table "classquality" of the paper reports for `cv': the error-free set,
# the two stock Illumina models "MiSeq" and "HiSeq", and the saliva-matched one. They carry their
# ground truth in the read name, `>ACCESSION|kraken:taxid|TAXID_...', which is what makes a per-read
# comparison possible at all -- the same names Genestrip is scored against.
#
# Only the first mate of each pair is classified, as everywhere else in this project: Genestrip has
# no notion of a pair and scores every record on its own, so reading the second mate would double
# the reads without adding an independent observation.
#
# What comes out is one tab-separated file per tool and read set under results/kraken, in each tool's
# own output format. Both start their lines with C or U, the read name and the assigned tax id, which
# is all the scoring needs.
#
# Usage:
#   sh ./bin/kraken_classify.sh            # both tools, all four read sets
#   sh ./bin/kraken_classify.sh k2         # only Kraken 2
#   sh ./bin/kraken_classify.sh ku         # only KrakenUniq
#
# Environment:
#   KRAKEN_PROJECT   the Genestrip project whose databases and read sets are used, default `viral'
#   THREADS          classification threads, default the number of CPUs
#   READ_SETS        the read set keys to classify, default the four of the paper
#
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)

what=${1:-all}
project=${KRAKEN_PROJECT:-viral}
threads=${THREADS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}
sets=${READ_SETS:-iss_perfect iss_miseq iss_hiseq iss_saliva}

fastqdir="${basedir}/data/fastq"
outdir="${basedir}/results/kraken"
k2db="${basedir}/data/kraken/${project}_k2"
kudb="${basedir}/data/kraken/${project}_ku"
k2bin="${basedir}/tools/kraken2/bin"
kubin="${basedir}/tools/krakenuniq/bin"

mkdir -p "$outdir"

# The fastq file of a read set key, as make_fastqs.sh named it.
fastq_of() {
  echo "${fastqdir}/${project}_${1}_reads_R1.fastq.gz"
}

run_tool() {
  _rt_tool=$1; _rt_db=$2; _rt_bin=$3
  if [ ! -d "$_rt_db" ]; then
    echo "No database at ${_rt_db} - run ./bin/kraken_build.sh first." >&2
    return 1
  fi
  for key in $sets; do
    _rt_fq=$(fastq_of "$key")
    _rt_out="${outdir}/${project}_${_rt_tool}_${key}.tsv"
    if [ ! -s "$_rt_fq" ]; then
      echo "  SKIP ${key} -- ${_rt_fq} is not on disk" >&2
      continue
    fi
    if [ -s "$_rt_out" ]; then
      echo "SKIP  ${_rt_out} exists"
      continue
    fi
    echo "=== ${_rt_tool}: ${key} ==="
    case "$_rt_tool" in
      k2) "${_rt_bin}/kraken2" --threads "$threads" --db "$_rt_db" --gzip-compressed \
              --output "$_rt_out" --report "${_rt_out%.tsv}.report" "$_rt_fq" ;;
      ku) "${_rt_bin}/krakenuniq" --threads "$threads" --db "$_rt_db" --gzip-compressed \
              --output "$_rt_out" --report-file "${_rt_out%.tsv}.report" "$_rt_fq" ;;
    esac
  done
}

case "$what" in
  k2|all)
    [ -x "${k2bin}/kraken2" ] || { echo "Kraken 2 is missing - run ./bin/install_tools.sh first." >&2; exit 1; }
    run_tool k2 "$k2db" "$k2bin"
    ;;
esac

case "$what" in
  ku|all)
    [ -x "${kubin}/krakenuniq" ] || { echo "KrakenUniq is missing - run ./bin/install_tools.sh first." >&2; exit 1; }
    run_tool ku "$kudb" "$kubin"
    ;;
esac

case "$what" in
  k2|ku|all) ;;
  *) echo "Usage: $0 [k2|ku|all]" >&2; exit 1 ;;
esac

echo
ls -la "$outdir"/*.tsv 2>/dev/null || echo "no classifications yet"
