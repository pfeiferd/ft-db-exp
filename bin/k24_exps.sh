#!/bin/sh
#
# The k = 24 control of `cv': a second, complete evaluation strand for Genestrip, unrefined and
# refined, on the same genomes, the same taxonomy and the same four read sets as `viral', differing
# in the k-mer length alone.
#
# Why k = 24. Kraken 2 classifies with a spaced seed of weight 24 -- 31 minimizer positions of which
# 7 are masked -- whereas Genestrip matches exact 31-mers. The weight of a seed is the number of
# positions that must agree, and under independently placed substitutions it is the whole story:
# a spaced seed of weight 24 and a contiguous 24-mer match with the same probability, both against a
# read carrying errors and against a related species. So a Genestrip database at k = 24 sits where
# Kraken 2 sits on that trade-off, and the strand says what the trade costs and what it buys.
#
# It also answers a question of the paper's own: what governs the gain of the refinement. k is the
# one knob that moves the share of k-mers stored above the data taxa -- the mass a refinement can
# push down -- without touching the genomes, the read sets, the measures or the code. The strand
# therefore reports the same quantities as the main run and nothing new: the share above the data
# taxa, sp and sp* for both variants, and the per-read precisions on the four read sets.
#
# Three steps, in this order:
#
#   1. bin/k24_projects.sh writes the twin `viral-k24': symlinks to the original's inputs, a
#      config.properties generated from the original's with `kMerSize=24', and the fastq maps linked
#      under the twin's name so that it is scored on the very reads `viral' was scored on.
#   2. run_exps.sh builds and times it, and computes the reports over its contents -- dbinfo,
#      branchhistocsv, dbquality, ftquality. Those are the point here, unlike for the `-sa' twins of
#      store_compare.sh, where they would recompute figures that must come out identical. CSV_SUFFIX
#      keeps the disk sizes and the timings in CSVs of their own instead of overwriting the main
#      run's. The genomes are extracted afterwards, because the accuracy runs resolve their ground
#      truth against the twin's own accession-to-taxon table.
#   3. run_classification_exps.sh scores the four read sets and times the classification.
#
# Wall time is about one to two hours on the machine the paper reports, and the twin needs some
# 8 GB of disk beside the original. The reads are not simulated again, so nothing here depends on
# InSilicoSeq.
#
# Usage:
#   sh ./bin/k24_exps.sh all          # all three steps
#   sh ./bin/k24_exps.sh projects     # only step 1, seconds
#   sh ./bin/k24_exps.sh build        # only step 2
#   sh ./bin/k24_exps.sh accuracy     # only step 3
#
# Environment:
#   K24_PROJECT   the project to make the twin of, default `viral'. Only `cv' is reported; the
#                 others would each cost a database build of hours for a control the paper does not
#                 draw on.
#   K             the k-mer length of the twin, default 24. It names the twin as well, so K=27
#                 builds and evaluates `viral-k27' throughout.
#
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)

what=${1:-}
src=${K24_PROJECT:-viral}
k=${K:-24}
twin="${src}-k${k}"

if [ -z "$what" ]; then
  echo "Usage: $0 [projects|build|accuracy|all]" >&2
  echo "  No default: 'all' builds a database and classifies four read sets, which takes hours." >&2
  echo "  Name the step." >&2
  exit 1
fi

case "$what" in
  projects|all)
    echo "############ 1/3: the k = ${k} twin of ${src} ############"
    K="$k" sh ./bin/k24_projects.sh "$src"
    ;;
esac

case "$what" in
  build|all)
    echo "############ 2/3: building and timing ${twin} ############"
    CSV_SUFFIX="-k${k}" sh ./bin/run_exps.sh "$twin"

    # The accession-to-taxon table the accuracy runs resolve their ground truth against, plus the
    # per-accession FASTAs the goal writes on the way. run_exps.sh does not make it -- in the main
    # run make_fastqs.sh does, because it also needs the FASTAs to simulate from -- and
    # ExtractedTaxIds.load returns an empty map for a missing file rather than failing, so without
    # this the accuracy step dies with "No extracted genomes to resolve the ground truth against".
    #
    # Via `extractrefseqcsv', not `extractrefseqfasta': the latter is an object goal, which
    # Genestrip treats as a weak dependency and never makes when asked for it on the command line.
    # The guard tests both outputs, for the reason make_iss() in make_fastqs.sh gives.
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
    K24_DB="$twin" sh ./bin/run_classification_exps.sh k24
    ;;
esac

case "$what" in
  projects|build|accuracy|all) ;;
  *) echo "Usage: $0 [projects|build|accuracy|all]" >&2; exit 1 ;;
esac

echo
echo "=== results ==="
ls -la "results/${twin}_dbinfo.csv" "results/${twin}_dbquality.csv" "results/${twin}_ftquality.csv" \
       "results/db_gen_perf-k${k}.csv" "results/db_disk_sizes-k${k}.csv" \
       "results/${twin}_iss_accuracy.csv" "results/${twin}_iss_perfect_accuracy.csv" \
       "results/${twin}_iss_saliva_accuracy.csv" 2>/dev/null \
  || echo "not every file is there yet - see the steps above"
