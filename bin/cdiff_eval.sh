#!/bin/sh
#
# The evaluation of the `cdiff' database, which is refined *below* the species rank and therefore
# needs a ground truth the taxonomy cannot supply. This script carries it as far as that ground
# truth: it makes sure the database's inventory of nodes exists and assigns a sequence type to every
# genome that became a leaf of it.
#
#   sh ./bin/cdiff_eval.sh              # write dbinfo if missing, then type
#   sh ./bin/cdiff_eval.sh --force      # write dbinfo again even if it is already there
#
# Run it after the database exists, i.e. after run_exps.sh has built `cdiff'.
#
# WHY dbinfo IS PART OF THIS. What has to be typed is one genome per leaf of the dendrogram, and the
# database is the only place that records which genomes those are: `genbank.maxPerTaxid' admits a
# subset of the assemblies, and `additional.txt' contributes leaves under other taxa. `dbinfo' writes
# the tree out with one row per node, so bin/mlst_assemblies.sh can read the FILE nodes below the
# requested tax ids straight off it -- see there for why the join key that comes with them is the
# whole point.
#
# THE EXTRACTION IS GONE, and it is worth saying what it was for. `cdiff' used to be filled from the
# RefSeq release as well, and a release file is a chunk of many organisms rather than a genome, so
# typing it as it stands would have been meaningless; `extractrefseqcsv' wrote one fasta per sequence
# region instead, which then had to be regrouped into assemblies by parsing accessions. With
# `refseq.filldb=false' and one Genbank fasta per assembly there is no chunk to unpack and no group
# to reconstruct, and the detour cost the join: its rows were keyed by an accession stem, which is
# not what a leaf is called. Should a project need it again, the goal to ask for is
# `extractrefseqcsv' and not `extractrefseqfasta' -- the latter is an ObjectGoal, which Maker.make()
# skips as a weak dependency, so it reports success within a second and does nothing at all.
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)

db=cdiff
force=""
[ "$1" = "--force" ] && force=1

dbinfo="${basedir}/data/projects/${db}/csv/${db}_dbinfo.csv"

echo "############ 1/2  the node inventory of ${db} ############"
if [ -f "$dbinfo" ] && [ -z "$force" ]; then
  echo "  SKIP  ${dbinfo} is already there"
  echo "        (pass --force to write it again; do that whenever the database was rebuilt,"
  echo "         since a dbinfo older than the database beside it describes another tree)"
else
  mvn exec:exec@db -Dname="$db" -Dgoal=dbinfo
  if [ ! -f "$dbinfo" ]; then
    echo "  WARNING: the goal did not write ${dbinfo}." >&2
    exit 1
  fi
  echo "  wrote ${dbinfo}"
fi

echo "############ 2/2  typing the leaves ############"
sh ./bin/mlst_assemblies.sh "$db"

echo
echo "Done. results/${db}_mlst.csv holds one sequence type per leaf of the database, keyed by the"
echo "leaf's name - which is the fasta file the genome was read from, so it joins to the dendrogram"
echo "directly."
echo "Still missing for the evaluation: the classification of the isolate reads of PRJNA1148956"
echo "against the unrefined and the refined database, and the sub-species precision scored against"
echo "these sequence types. Fetch the reads with bin/fetch_cdiff.sh."
