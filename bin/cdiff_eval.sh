#!/bin/sh
#
# The evaluation of the `cdiff' database, which is refined *below* the species rank and therefore
# needs a ground truth the taxonomy cannot supply. This script carries it as far as that ground
# truth: it extracts the genomes the database was filled from and assigns a sequence type to each.
#
#   sh ./bin/cdiff_eval.sh              # extract if needed, then type
#   sh ./bin/cdiff_eval.sh --force      # extract again even if the fasta folder is populated
#
# Run it after the database exists, i.e. after run_exps.sh has built `cdiff'.
#
# WHY THE EXTRACTION IS PART OF THIS. The genomes are not lying about ready to be typed. `cdiff'
# draws on two sources at once -- some 320 GenBank assemblies under data/common/genbank and some 190
# RefSeq release files under data/common/refseq -- and a RefSeq release file is a chunk of many
# organisms rather than a genome, so typing it as it stands would be meaningless. The extraction
# writes one fasta per sequence region that actually went into the database, from both sources, and
# that is what can be typed.
#
# THE GOAL TO ASK FOR IS `extractrefseqcsv'. The extraction itself is `extractrefseqfasta', but that
# one is an ObjectGoal, and Maker.make() aggregates the requested goals as dependencies of an
# internal goal whose make() skips every weak dependency -- which every ObjectGoal is by default.
# Asking for it directly therefore prints BUILD SUCCESS after a second and leaves the folder empty.
# `extractrefseqcsv' is a FileListGoal, so it runs, and it pulls the extraction in as its own
# dependency. It writes the map from sequence description to tax id as well, which is what tells an
# extracted file which taxon it came from.
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)

db=cdiff
force=""
[ "$1" = "--force" ] && force=1

fastadir="${basedir}/data/projects/${db}/fasta"

echo "############ 1/2  extracting the genomes of ${db} ############"
# Counting rather than testing for the folder: Genestrip creates it early and leaves it empty when
# the extraction did not run, which is exactly the state this script exists to get out of.
extracted=0
if [ -d "$fastadir" ]; then
  extracted=$(find "$fastadir" -maxdepth 1 -type f \( -name '*.fa' -o -name '*.fa.gz' \) 2>/dev/null | wc -l)
fi
if [ "$extracted" -gt 0 ] && [ -z "$force" ]; then
  echo "  SKIP  ${extracted} fasta file(s) already in ${fastadir}"
  echo "        (pass --force to extract them again)"
else
  echo "  this reads every source file of the database and takes a while"
  mvn exec:exec@db -Dname="$db" -Dgoal=extractrefseqcsv
  extracted=$(find "$fastadir" -maxdepth 1 -type f \( -name '*.fa' -o -name '*.fa.gz' \) 2>/dev/null | wc -l)
  echo "  extracted ${extracted} fasta file(s) into ${fastadir}"
  if [ "$extracted" -eq 0 ]; then
    echo "  WARNING: the extraction produced nothing. The goal reports what it read on stderr;" >&2
    echo "           a run of a second or two means it did nothing at all." >&2
    exit 1
  fi
fi

echo "############ 2/2  typing them ############"
sh ./bin/mlst_assemblies.sh "$db"

echo
echo "Done. results/${db}_mlst.csv holds one sequence type per extracted genome."
echo "Still missing for the evaluation: the classification of the isolate reads of PRJNA1148956"
echo "against the unrefined and the refined database, and the sub-species precision scored against"
echo "these sequence types. Fetch the reads with bin/fetch_cdiff.sh."
