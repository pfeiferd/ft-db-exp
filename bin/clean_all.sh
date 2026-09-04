#!/bin/sh
#
# Deletes everything the experiments generate, and nothing that was downloaded.
#
# Run it when results must come from one code state rather than from whatever happened to be on
# disk. It is also step 3 of run_all_exps.sh under FRESH_DBS=1.
#
# Usage:
#   sh ./bin/clean_all.sh              delete
#   DRY_RUN=1 sh ./bin/clean_all.sh    list what would be deleted, delete nothing
#
# WHAT SURVIVES, and why it is safe. `cleanall' cleans a goal and then its dependencies, but only
# those that permit it: Goal.transitiveClean() consults isAllowTransitiveClean(), which reads the
# transClean flag of the goal key. Every download goal is declared with that flag false --
# commonsetup, setup, taxdownload, refseqcat, checksummap, refseqfna, assemblydownload,
# fastasgenbankdl, adddownloads, fastqdownload, fastadownload, dbdownload -- so the RefSeq and
# GenBank sequence, the taxonomy dump, the assembly summary and every downloaded fastq are out of
# reach of this script by construction rather than by our care.
#
# `cleantotal' is the target that ignores that flag and would take the downloads with it. This
# script does not use it, and neither should anything else here: re-downloading the RefSeq bacteria
# division to save a few gigabytes is not a trade worth making.
#
# WHAT THE GOALS DO NOT REACH, which is why this script does more than call them. The simulated
# reads are written by InSilicoSeq and NanoSim, which are external tools rather than Genestrip
# goals, so nothing in the goal graph cleans them; and results/ is filled by a shell copy. Both are
# derived, both go stale when a database's genome selection changes, and both are removed here.
#
# The per-accession FASTAs under data/projects/<db>/fasta are removed explicitly for the same
# reason. They are not downloads -- `extractrefseqfasta' derives them from data/common/refseq -- and
# they encode which genomes the per-species cap admitted. make_fastqs.sh skips its extraction when
# that directory is non-empty, so leaving it behind means the next run simulates reads from a
# genome selection the rebuilt database no longer holds. That is the trap this line exists to shut.
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."

if [ -n "${DRY_RUN:-}" ]; then
  echo "DRY_RUN: nothing will be deleted."
fi

# $1 = a path to remove, with its reason for the log.
drop() {
  [ -e "$1" ] || return 0
  if [ -n "${DRY_RUN:-}" ]; then
    echo "  would delete  $1"
  else
    echo "  deleting      $1"
    rm -rf "$1"
  fi
}

echo "### databases and per-project generated files"
for p in data/projects/*/; do
  name=$(basename "$p")
  if [ -z "${DRY_RUN:-}" ]; then
    # The goal graph first: it knows what it made, and it will not touch a download.
    mvn exec:exec@cleanall -Dname="${name}" -Dgoal=ftgenall
    mvn exec:exec@db -Dname="${name}" -Dgoal=ftclear
  else
    echo "  would run     cleanall/ftclear for ${name}"
  fi
  # Then the two things the goals leave behind. `db' because a stored database that survives a
  # config change is silently reused, the stored accession map included, which is where a
  # per-species genome cap lives; `fasta' for the reason given at the head of this file.
  drop "${p}db"
  drop "${p}fasta"
  drop "${p}csv"
  drop "${p}log"
  drop "${p}krakenout"
  drop "${p}tex"
done

echo "### simulated reads"
# Named by the simulation maps themselves rather than by a pattern, so a downloaded fastq cannot be
# caught: a map entry with a URL in its second field is a download and is skipped here. Only the
# maps whose names carry `_sim' are read, and those declare local files by construction.
for m in data/fastq/*_sim*.txt; do
  [ -e "$m" ] || continue
  awk '$1 !~ /^#/ && NF >= 2 && $2 !~ /^https?:/ { print $2 }' "$m" | sort -u | while read -r f; do
    drop "data/fastq/${f}"
  done
  drop "$m"
done

echo "### collected results"
drop results

echo
echo "Done. Kept: data/common (RefSeq and GenBank sequence, taxonomy, assembly summary) and every"
echo "downloaded fastq under data/fastq. Rebuild with: sh ./bin/run_all_exps.sh"
