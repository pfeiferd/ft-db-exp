#!/bin/sh
#
# Clears the generated state of the named projects alone, so that a configuration change really
# takes effect, and keeps everything whose regeneration would cost hours or change the inputs.
#
# Why this exists beside clean_all.sh: that script is for "results must come from one code state"
# and it takes the whole experiment with it -- every project, the simulated reads, and all of
# results/ including results/kraken. After a change to one project family that is far too much.
# Re-simulating the reads would also change them, and then the new numbers could no longer be held
# against the tables that were not rebuilt.
#
# WHY DELETING IS NECESSARY AT ALL. Genestrip's goals are make-like: a goal whose output file is
# already there reports itself made and does nothing. A changed config.properties does not
# invalidate anything on disk. In particular `clear'/`ftclear' removes only krakenout/, csv/ and
# log/ (GSMaker.java: the CLEAR goal's file list) and leaves db/ alone, so a database built under
# the old configuration is silently reused. That is what this script removes.
#
# WHAT IT DELETES, per project:
#   csv/ log/ krakenout/   through the goal graph, with `ftclear'
#   db/                    explicitly -- the databases, their k-mer index and their Bloom filter
#   results/<p>_* etc.     the collected copies, so that a step which silently did not re-run shows
#                          up as a missing file instead of as a stale number in the paper
#
# WHAT IT KEEPS, and why:
#   data/projects/<p>/fasta   the extracted per-accession genomes. A setting like maxDust changes
#                             what is *stored* from them, never which genomes are selected, and the
#                             simulated reads were drawn from exactly these files -- re-extracting
#                             would cost hours and risk a different read set.
#   data/fastq                every simulated and downloaded fastq, for the same reason.
#   results/kraken            the per-read output of KrakenUniq and Kraken 2. Re-classifying the
#                             real runs is the single most expensive step of the whole experiment,
#                             and nothing about a Genestrip-side change invalidates it.
#   every other project       untouched.
#
# `ftclear' takes csv/ with it, and that holds <db>_extractrefseqcsv.csv, the accession-to-taxon
# table the ground truth is resolved against. The goals rebuild it, but a classification run that
# starts without it dies hours in with "No extracted genomes to resolve the ground truth against",
# so this script regenerates it right away rather than leaving that to chance.
#
# Usage:
#   sh ./bin/clear_projects.sh viral viral-k24 viral-s4 viral-k24-s4
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
