#!/bin/sh
#
# Clears the generated state of the named projects alone, so that a configuration change really takes
# effect, and keeps everything whose regeneration would cost hours or change the inputs.
#
# Why beside clean_all.sh: that script is for "results must come from one code state" and takes the
# whole experiment with it, the simulated reads included -- and re-simulating changes them, after
# which the new numbers cannot be held against the tables that were not rebuilt.
#
# WHY DELETING IS NECESSARY AT ALL. Genestrip's goals are make-like: a goal whose output file is
# there reports itself made and does nothing, and a changed config.properties invalidates nothing on
# disk. `clear'/`ftclear' removes only krakenout/, csv/ and log/ and leaves db/ alone, so a database
# built under the old configuration is silently reused.
#
# WHAT IT DELETES, per project:
#   csv/ log/ krakenout/   through the goal graph, with `ftclear'
#   db/                    explicitly -- the databases, their k-mer index and their Bloom filter
#   results/<p>_* etc.     the collected copies, so that a step which silently did not re-run shows
#                          up as a missing file instead of as a stale number in the paper
#
# WHAT IT KEEPS: the extracted per-accession genomes under data/projects/<p>/fasta and everything in
# data/fastq, because a setting like maxDust changes what is stored from them and not which genomes
# are selected, and the simulated reads were drawn from exactly these files; results/kraken, because
# re-classifying the real runs is the most expensive step of the whole experiment and no
# Genestrip-side change invalidates it; and every other project.
#
# `ftclear' takes csv/ with it, and that holds <db>_extractrefseqcsv.csv, the accession-to-taxon table
# the ground truth is resolved against. A classification run that starts without it dies hours in
# with "No extracted genomes to resolve the ground truth against", so this regenerates it right away.
#
# Usage:
#   sh ./bin/clear_projects.sh viral viral-k24 viral-k24-s4
#   DRY_RUN=1 sh ./bin/clear_projects.sh viral        list what would happen, do nothing
#
set -e
scriptdir=$(dirname "$0")
cd "$scriptdir/.."

[ $# -gt 0 ] || { echo "Usage: sh ./bin/clear_projects.sh <project> [project ...]" >&2; exit 1; }

if [ -n "${DRY_RUN:-}" ]; then
  echo "DRY_RUN: nothing will be deleted and no goal will run."
fi

drop() {
  [ -e "$1" ] || return 0
  if [ -n "${DRY_RUN:-}" ]; then
    echo "  would delete  $1"
  else
    echo "  deleting      $1"
    rm -rf "$1"
  fi
}

for p in "$@"; do
  [ -d "data/projects/${p}" ] || { echo "No such project: data/projects/${p}" >&2; exit 1; }
  echo "### ${p}"
  if [ -n "${DRY_RUN:-}" ]; then
    echo "  would run     ftclear for ${p}  (removes csv/, log/, krakenout/)"
  else
    mvn exec:exec@db -Dname="${p}" -Dgoal=ftclear
  fi
  drop "data/projects/${p}/db"
  # The collected copies. The patterns are anchored on the project name followed by `_' or `.', so
  # `viral' never catches `viral-k24'. One exception is kept: <p>_simparams.csv describes the
  # simulated *reads* -- their measured length and per-base error -- and the reads survive here, so
  # throwing it away would only force make_fastqs.sh to measure them again. It has gone missing that
  # way before, and Table 6's tick columns rendered empty for weeks as a result.
  for g in "results/${p}_"* "results/match_${p}_"* "results/ftmatch_${p}_"* \
           "results/db_gen_${p}."* "results/ftdb_gen_${p}."*; do
    case "$g" in
      *_simparams.csv) [ -e "$g" ] && echo "  keeping       $g" ;;
      *) drop "$g" ;;
    esac
  done
  if [ -n "${DRY_RUN:-}" ]; then
    echo "  would run     extractrefseqcsv for ${p}"
  else
    mvn exec:exec@db -Dname="${p}" -Dgoal=extractrefseqcsv
  fi
done

echo
echo "Done. Kept: data/projects/*/fasta, data/fastq, results/kraken and every other project."
echo "results/db_gen_perf*.csv is shared between projects and was left alone -- re-measure it with"
echo "sh ./bin/db_gen_perf.sh once the databases are back."
