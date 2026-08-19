#!/bin/sh
#
# Decides whether a genus is worth a full FT case study, before any of it is built.
#
#   sh ./bin/prescreen.sh streptococcus
#   sh ./bin/prescreen.sh streptococcus nocardia
#
# WHAT THIS ANSWERS. The refinement can only claim k-mers that a *subset* of a node's children
# carries. A k-mer in one child alone is already below the node; a k-mer in all of them is core
# genome and stays where it is. Everything between those two is the headroom, and a genus either
# has it or it does not. The C. difficile study was abandoned after three full builds because it
# did not, and the number that settled it could have been had in an afternoon. This script is that
# afternoon.
#
# The measurement is the `branchhisto' goal: per refined node, how many of its k-mers occur in
# exactly 1, 2, ... of its child subtrees, plus a trailing OTHER bucket. It needs the unrefined
# database and the k-mer index, and nothing beyond them -- in particular it does NOT need
# `dendrogram', `ftupdatedb' or the refined database, which are the expensive parts.
#
# READING THE RESULT. Per refined node the summary below reports where its k-mers sit:
#
#   deg 1        k-mers seen in one child subtree. Large here is normal and means nothing either
#                way; these are not what the refinement moves.
#   deg 2..10    the headroom. A refined node can hold these, and a tight group of two to ten
#                species is what a clinically meaningful subgroup looks like. THIS IS THE NUMBER
#                THE DECISION HANGS ON.
#   deg > C/2    carried by more than half the children: core genome of the genus. Unrefinable, and
#                a genus whose above-species k-mers are nearly all here is a genus to walk away from.
#   OTHER        the trailing bucket, i.e. carried by something outside the node's own children.
#
# A rough rule from the databases that did work: if `deg 2..10' is under a couple of per cent of the
# node's k-mers, stop. If it is a double-digit percentage, build the case study.
#
# COST. The database build dominates: it reads the RefSeq division files for the project's
# categories. Everything after it is minutes. A database that already exists is not rebuilt.
#
# For streptococcus the build is the whole expense -- 22,664 assemblies and 46.6 Gbp, against
# tick-borne's 4,135 and 5.8 Gbp -- so its config.properties currently sets maxGenomesPerTaxid=50 /
# maxPerTaxidRank=species, which brings the store down to tick-borne scale without losing a single
# species. That cap applies to the fill and to the quality goals, but NOT to the LCA update and NOT
# to the k-mer index, both of which override the capping code path on purpose; see the "DATABASE
# SIZE" section of that file for the trace. So the branching-degree shape reported below is
# computed over all 22,664 genomes either way, and the screen answers the same question the full
# build would. Remove the two lines for the final build, and delete
# data/projects/streptococcus/db first, or the screening database is silently reused.
set -e

scriptdir=$(dirname "$0")
cd $scriptdir/..

if [ $# -eq 0 ]; then
  echo "usage: sh ./bin/prescreen.sh <project> [<project> ...]" >&2
  exit 2
fi

mkdir -p results logs

for p in "$@"; do
  if [ ! -d "data/projects/$p" ]; then
    echo "prescreen: no such project: data/projects/$p" >&2
    exit 1
  fi

  echo
  echo "=== $p: building the unrefined database (skipped if it is already there) ==="
  mvn exec:exec@db -Dname=$p -Dgoal=db

  echo
  echo "=== $p: k-mers per node of the unrefined database ==="
  mvn exec:exec@db -Dname=$p -Dgoal=dbinfo

  # branchhisto pulls in kmerindexbloom -- and, under the default AUTO sizing, kmerindexsize before
  # it, which reads the genomes once more to sketch how many (k-mer, child) pairs the filter will
  # hold. That extra pass is the price of not sizing the filter against the conservative bound.
  echo
  echo "=== $p: branching-degree histograms ==="
  mvn exec:exec@db -Dname=$p -Dgoal=branchhistocsv
  mvn exec:exec@db -Dname=$p -Dgoal=branchhistorankcsv

  # The Jaccard matrix per refined node. Cheap once the index exists, and the second opinion on the
  # same question: a matrix that is near-constant off the diagonal is the index filter's noise
  # floor rather than biology, and means ftIndexBloomFilterFpp is set too high.
  echo
  echo "=== $p: pairwise intersection counts and Jaccard indices ==="
  mvn exec:exec@db -Dname=$p -Dgoal=intersectcsv

  csv=data/projects/$p/csv/${p}_branchhistocsv.csv
  if [ ! -f "$csv" ]; then
    echo "prescreen: expected $csv, which is not there -- look at the goal's output above" >&2
    exit 1
  fi
  cp "$csv" results/${p}_branchhistocsv.csv
  cp data/projects/$p/csv/${p}_branchhistorankcsv.csv results/ 2>/dev/null || true
  cp data/projects/$p/csv/${p}_dbinfo.csv results/ 2>/dev/null || true

  echo
  echo "=== $p: VERDICT ==="
  awk -F';' '
    NR == 1 { next }                       # header
    NF < 5  { next }
    {
      taxid = $1; nm = $2; C = $3 + 0
      # Buckets are $4 .. $(3 + C + 1): degree 1 .. C, then the trailing OTHER bucket.
      total = 0; d1 = 0; d2_10 = 0; core = 0
      other = $(3 + C + 1) + 0
      for (i = 1; i <= C + 1; i++) {
        v = $(3 + i) + 0
        total += v
        if (i > C) continue              # the OTHER bucket, counted separately
        if (i == 1)             d1    += v
        else if (i <= 10)       d2_10 += v
        if (i > C / 2)          core  += v
      }
      if (total == 0) next
      printf "%d\t%-9s %-40.40s C=%-5d kmers=%-12d deg1=%5.1f%%  deg2-10=%5.1f%%  deg>C/2=%5.1f%%  OTHER=%5.1f%%\n", \
             total, taxid, nm, C, total, 100*d1/total, 100*d2_10/total, 100*core/total, 100*other/total
    }
  ' "$csv" | sort -rn -k1,1 | cut -f2- | head -20

  echo
  echo "    (full table in results/${p}_branchhistocsv.csv; the rows above are the refined nodes,"
  echo "     the biggest first. deg2-10 is the headroom -- see the header of this script.)"
done
