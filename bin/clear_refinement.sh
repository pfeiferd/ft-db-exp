#!/bin/sh
#
# Deletes exactly what the OTHER-slot refinement defect invalidated, and nothing else.
#
# THE DEFECT. UpdateStoreGoal built each refined node's bit set with one slot per child and none for
# the trailing OTHER slot, while the membership vector a k-mer arrives with has both. A refined node
# whose cluster contains OTHER could therefore never be seen to cover it, and every k-mer an unnamed
# genome touches stayed at the node being refined instead of moving down. Measured on orthopox: the
# genus kept 233,212 k-mers where it should keep 135,738, so 97,474 of 239,260 -- two in five -- were
# stranded. The refinement still produced a correct database, only a far less specific one, so what
# follows is not corrupted data but understated results.
#
# WHAT IS AFFECTED. Only the refined database and what reads it:
#
#   ftupdatedb -> ftdb -> ftdbinfo, ftquality, ftsvgtaxtree, and every classification run of the
#   goal `ftmatch' together with the accuracy, specificity and taxon-call numbers derived from it.
#
# WHAT IS NOT. Everything before the refinement, which is most of the cost:
#
#   the accession map, the downloaded RefSeq archives and the extracted fastas; the unrefined
#   database `db'; the k-mer index `storekmerindex', which is built and serialised before the
#   refinement runs; the intersection counts and the dendrogram, which were both correct -- the
#   clustering put OTHER inside a refined node exactly as it should, and only the update phase could
#   not act on it; the simulated and real fastq files, which are drawn from the extracted genomes and
#   never see a database; and every result of the goal `match', which uses the unrefined database.
#
# So this script deletes the refined databases and the results that read them, and leaves the
# expensive passes over the genomes alone. Afterwards run:
#
#   sh ./bin/run_all_exps.sh
#
# which rebuilds `ftdb' for each project and redoes the classification. Do NOT pass FRESH_DBS=1:
# that would delete the unrefined databases too and cost days for nothing.
#
# Usage:
#   sh ./bin/clear_refinement.sh            # delete
#   sh ./bin/clear_refinement.sh --dry-run  # list what would be deleted and delete nothing
#
set -e

dry=""
case "${1:-}" in
  --dry-run|-n) dry=1 ;;
  "") ;;
  *) echo "usage: $0 [--dry-run]" >&2; exit 2 ;;
esac

cd "$(dirname "$0")/.."
basedir=$(pwd)
res_path="${basedir}/results"

n=0
bytes=0

# Refuses outright to touch anything the refinement did not produce. The unrefined databases, the
# serialised k-mer indices, the extracted fastas and the fastq files are days of work and are sound;
# a mistaken pattern above must fail loudly here rather than quietly cost a week.
protected() {
  case "$1" in
    *_db.zip|*/db/*_db|*_storekmerindex.ser.gz|*/fasta/*|*/data/fastq/*|*/data/common/*) return 0 ;;
    */db_gen_*.log) return 0 ;;
  esac
  return 1
}

# Removes a path if it exists and says so. A glob that matched nothing arrives here as the literal
# pattern, which is why every caller guards with -e and why this does too.
drop() {
  [ -e "$1" ] || return 0
  if protected "$1"; then
    echo "REFUSING to delete ${1#$basedir/} - it is not a product of the refinement." >&2
    echo "This is a guard against a mistaken pattern; nothing has been deleted." >&2
    exit 3
  fi
  sz=$(du -sk "$1" 2>/dev/null | cut -f1)
  sz=${sz:-0}
  bytes=$((bytes + sz))
  n=$((n + 1))
  if [ -n "$dry" ]; then
    printf '  would delete  %8s MB  %s\n' "$((sz / 1024))" "${1#$basedir/}"
  else
    printf '  deleting      %8s MB  %s\n' "$((sz / 1024))" "${1#$basedir/}"
    rm -rf "$1"
  fi
}

# Every project that has a project directory, so a database added later is covered without editing
# this list.
projects=$(ls -1 "${basedir}/data/projects" 2>/dev/null || true)

echo "=== refined databases (goal ftdb) ==="
for db in $projects; do
  # Both forms: Genestrip keeps a database as a zip and may also leave it unpacked beside it.
  drop "${basedir}/data/projects/${db}/db/${db}_ftdb.zip"
  drop "${basedir}/data/projects/${db}/db/${db}_ftdb"
done

echo "=== everything a project derives from the refined database ==="
# By pattern rather than by a list of goal names, so a goal added later is covered and a goal I have
# mis-remembered is not missed. Every artefact of the refinement is named <db>_ft*; nothing that
# comes from before it is -- dbinfo, dbquality, branchhistocsv, branchhistorankcsv, intersectcsv,
# kmerrankstatscsv, dendrolatex, extractrefseqcsv and the nanosim list all lack the prefix, and are
# left alone deliberately.
for db in $projects; do
  dir="${basedir}/data/projects/${db}/csv"
  if [ -d "$dir" ]; then
    for f in "${dir}/${db}_ft"*; do
      if [ -e "$f" ]; then drop "$f"; fi
    done
  fi
  # Classification output of the refined database, which Genestrip writes outside csv/ as well.
  for sub in krakenout tex; do
    d2="${basedir}/data/projects/${db}/${sub}"
    if [ -d "$d2" ]; then
      for f in "${d2}/${db}_ft"*; do
        if [ -e "$f" ]; then drop "$f"; fi
      done
    fi
  done
done

echo "=== collected results ==="
# The results folder is refilled by run_exps.sh's copy loop and by paper_stats.sh, so deleting
# generously here costs nothing but a copy. Two things are kept on purpose: db_gen_*.log, because
# the unrefined databases are not being rebuilt and their timings cannot be recovered once gone,
# and every input that predates the refinement.
for f in "${res_path}"/*_ft[a-z]*.csv "${res_path}"/*_ftsvgtaxtree.svg "${res_path}"/ftmatch_*.csv \
         "${res_path}"/ftmatch_*.log "${res_path}"/ftdb_gen_*.log \
         "${res_path}"/*_accuracy.csv "${res_path}"/*_summary.csv "${res_path}"/*_specificity.csv \
         "${res_path}"/*_st_ftdb.csv "${res_path}"/matchperf.csv "${res_path}"/dbstats.tex \
         "${res_path}"/db_disk_sizes.csv "${res_path}"/db_gen_perf.csv; do
  if [ -e "$f" ]; then drop "$f"; fi
done

echo
if [ -n "$dry" ]; then
  echo "$n path(s), about $((bytes / 1024)) MB -- nothing was deleted."
  echo "Re-run without --dry-run to delete."
else
  echo "Deleted $n path(s), about $((bytes / 1024)) MB."
  echo
  echo "Left in place on purpose: the unrefined databases, the k-mer indices"
  echo "(<db>_storekmerindex.ser.gz), the RefSeq archives, the extracted fastas and every fastq."
  echo
  echo "Now run:  sh ./bin/run_all_exps.sh"
  echo "Do not pass FRESH_DBS=1 -- the unrefined databases are sound and cost days to rebuild."
fi
