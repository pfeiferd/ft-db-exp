#!/bin/sh
#
# The sampled control of `cv': a complete evaluation strand for Genestrip, unrefined and refined, on
# the same genomes, the same taxonomy and the same four read sets as `viral', entering only one
# k-mer in n instead of every one of them.
#
# Why. Kraken 2's database is smaller than Genestrip's for two independent reasons, and only one of
# them is the seed weight that bin/k24_exps.sh isolates. The other is the store: Kraken 2 keeps a
# 4-byte cell per minimizer, 17 bits of truncated hash and 15 bits of taxon index, and enters only
# the minimizers, which at its defaults k = 35 and l = 31 is a window of five and so a density of
# 2/(5+1) = 1/3. Genestrip enters every k-mer in a 64-bit word that holds the k-mer exactly, plus a
# 10-bit-per-entry Bloom filter in front. Measured on `cv' that is 9.25 bytes per entry in memory
# against 5.71 for Kraken 2.
#
# So a Genestrip database of Kraken 2's size needs a sampling rate of about four, and this strand
# says what that costs where it has to be paid, namely on the reads. The sampling selects by the
# k-mer rather than by its position, so a read keeps about one k-mer in n at random positions; on
# clean reads that is harmless, on error-rich ones the surviving error-free k-mers are thinned by the
# same factor. That is the number this strand produces.
#
# Three steps, as in bin/k24_exps.sh, which describes the mechanics they share:
#
#   1. bin/sampling_projects.sh writes the twin `viral-s4'.
#   2. run_exps.sh builds and times it and reports over its contents. CSV_SUFFIX keeps the disk sizes
#      and timings in CSVs of their own. The genomes are extracted afterwards, for the twin's own
#      accession-to-taxon table.
#   3. run_classification_exps.sh scores the four read sets and times the classification.
#
# The twin is smaller and faster to build than the original, since a quarter of the k-mers go in.
#
# Usage:
#   sh ./bin/sampling_exps.sh all          # all three steps
#   sh ./bin/sampling_exps.sh projects     # only the twin project, seconds
#   sh ./bin/sampling_exps.sh build        # only the database and the reports over it
#   sh ./bin/sampling_exps.sh accuracy     # only the classification
#
# Environment:
#   SAMPLING_PROJECT  the project to make the twin of, default `viral'. Only `cv' is reported.
#   S                 one k-mer in how many to enter, default 4. It names the twin as well, so S=5
#                     builds and evaluates `viral-s5' throughout.
#   K                 a k-mer length to set along with the sampling, unset by default. K=24 builds
#                     and evaluates `viral-k24-s4', which makes the database smaller both ways at
#                     once -- Kraken 2's seed weight and a comparable entry count.
#
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)

what=${1:-}
src=${SAMPLING_PROJECT:-viral}
s=${S:-4}
k=${K:-}
# The twin's name and the CSV suffix carry every key that differs from the original, so two
# configurations can never write into one another's results.
suffix="${k:+-k${k}}-s${s}"
twin="${src}${suffix}"

if [ -z "$what" ]; then
  echo "Usage: $0 [projects|build|accuracy|all]" >&2
  echo "  No default: 'all' builds a database and classifies four read sets, which takes hours." >&2
  echo "  Name the step." >&2
  exit 1
fi

case "$what" in
  projects|all)
    echo "############ 1/3: the one-in-${s}${k:+, k = ${k}} twin of ${src} ############"
    K="$k" S="$s" sh ./bin/sampling_projects.sh "$src"
    ;;
esac

case "$what" in
  build|all)
    echo "############ 2/3: building and timing ${twin} ############"
    CSV_SUFFIX="$suffix" sh ./bin/run_exps.sh "$twin"

    # The twin's own accession-to-taxon table, for the reason bin/k24_exps.sh gives: the accuracy
    # runs resolve their ground truth against it, and ExtractedTaxIds.load returns an empty map for
    # a missing file rather than failing. The reads themselves are not simulated again -- they are
    # the files `viral' was scored on -- so the accessions of this table have to be the same ones.
    # ExtractRefSeqFastasGoal does read kMerSampling, but only for the per-tax-id k-mer cap, and
    # these projects set none; a warning about unresolvable ground truth in step 3 would be the sign
    # that this assumption is wrong.
    extractcsv="${basedir}/data/projects/${twin}/csv/${twin}_extractrefseqcsv.csv"
    fastadir="${basedir}/data/projects/${twin}/fasta"
    if [ -z "$(ls -A "$fastadir" 2>/dev/null)" ] || [ ! -s "$extractcsv" ]; then
      echo "=== ${twin}: extracting the genomes the database was built from ==="
      mvn exec:exec@db -Dname="$twin" -Dgoal=extractrefseqcsv
    fi
    if [ ! -s "$extractcsv" ]; then
      echo "Goal extractrefseqcsv produced no ${extractcsv}." >&2
      echo "Without it the accuracy step cannot resolve the ground truth." >&2
      exit 1
    fi
    ;;
esac

case "$what" in
  accuracy|all)
    echo "############ 3/3: the four read sets against ${twin} ############"
    TWIN_DB="$twin" sh ./bin/run_classification_exps.sh twin
    ;;
esac

case "$what" in
  projects|build|accuracy|all) ;;
  *) echo "Usage: $0 [projects|build|accuracy|all]" >&2; exit 1 ;;
esac

echo
echo "=== results ==="
ls -la "results/${twin}_dbinfo.csv" "results/${twin}_dbquality.csv" "results/${twin}_ftquality.csv" \
       "results/db_gen_perf${suffix}.csv" "results/db_disk_sizes${suffix}.csv" \
       "results/${twin}_iss_accuracy.csv" "results/${twin}_iss_perfect_accuracy.csv" \
       "results/${twin}_iss_saliva_accuracy.csv" 2>/dev/null \
  || echo "not every file is there yet - see the steps above"
