#!/bin/sh
#
# Measures Genestrip's two k-mer stores against each other on the databases of the paper's
# classification performance scenarios: `viral', `tick-borne' and `strepto'.
#
# The default store is RadixKMerStore, which resolves the low bits of a k-mer by address arithmetic
# and searches a bucket of a few dozen entries. Before it there was KMerSortedArray, which keeps its
# k-mers in one sorted array and binary-searches it. Both are still in the code, and a database is
# filled into the second one by setting `sortedArrayStore=true'. What that costs, in fill time, in
# refinement time, in classification throughput and in database size, is what the paper's appendix
# states and what this script produces.
#
# Three steps, in this order:
#
#   1. bin/sa_projects.sh writes the `-sa' twin of each project: symlinks to the original's inputs
#      and a config.properties generated from the original's, plus the store switch. Nothing there is
#      maintained by hand, so the twins cannot drift away from the databases they are compared with.
#   2. run_exps.sh builds and times the twins. DB_ONLY keeps it to the databases: the reports over
#      their contents would recompute figures that must come out identical, since the k-mers and the
#      taxonomy are the same. CSV_SUFFIX keeps their disk sizes and timings in CSVs of their own
#      instead of overwriting the originals'.
#   3. perf_scenarios.sh runs the four scenarios against the twins. DB_SUFFIX appends `-sa' to every
#      project name it classifies against and to the names of its logs and its CSV, so the two
#      measurements sit side by side and their rows line up.
#
# Wall time is a day or more -- the fill of `strepto' alone took 431 minutes with the radix store --
# and the twins need some 20 GB of disk beside the originals. Both are why this is a script of its
# own rather than a step of run_all_exps.sh.
#
# The step has to be named. There is no default, deliberately: `all' starts a build of days, and a
# bare call of a script is too easy a way to start one.
#
# Usage:
#   sh ./bin/store_compare.sh all          # all three steps, a day or more
#   sh ./bin/store_compare.sh projects     # only step 1, seconds
#   sh ./bin/store_compare.sh build        # only step 2
#   sh ./bin/store_compare.sh perf         # only step 3
#
# Environment:
#   SA_PROJECTS   the projects to compare, default `viral tick-borne strepto'. Mind the value
#                 capacity of the binary store: it indexes its values with a short and holds 65,535
#                 of them, which is enough for these three (`viral' needs 17,193) but not for a
#                 database of more taxa. FillDBGoal refuses such a fill rather than truncating it.
#
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."

what=${1:-}
projects=${SA_PROJECTS:-viral tick-borne strepto}

if [ -z "$what" ]; then
  echo "Usage: $0 [projects|build|perf|all]" >&2
  echo "  No default: 'all' builds three databases in both variants and measures them, which takes" >&2
  echo "  a day or more. Name the step." >&2
  exit 1
fi

# The twin of every project named, in the same order.
twins=""
for p in $projects; do
  twins="${twins} ${p}-sa"
done

case "$what" in
  projects|all)
    echo "############ 1/3: the -sa twins of ${projects} ############"
    sh ./bin/sa_projects.sh $projects
    ;;
esac

case "$what" in
  build|all)
    echo "############ 2/3: building and timing${twins} ############"
    DB_ONLY=1 CSV_SUFFIX=-sa sh ./bin/run_exps.sh $twins
    ;;
esac

case "$what" in
  perf|all)
    echo "############ 3/3: the four scenarios against the twins ############"
    DB_SUFFIX=-sa sh ./bin/perf_scenarios.sh all
    ;;
esac

case "$what" in
  projects|build|perf|all) ;;
  *) echo "Usage: $0 [projects|build|perf|all]" >&2; exit 1 ;;
esac

echo
echo "=== results ==="
ls -la results/db_gen_perf-sa.csv results/db_disk_sizes-sa.csv results/matchperf-sa.csv 2>/dev/null \
  || echo "not all three CSVs are there yet - see the steps above"
