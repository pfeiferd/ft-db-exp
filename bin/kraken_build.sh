#!/bin/sh
#
# Builds the Kraken 2 and KrakenUniq counterparts of the Genestrip database `viral', so that the
# three can be held against each other on the same reads: cv, cv_k2 and cv_ku.
#
# WHY ONLY `viral'. It is the one database of the paper another tool can be given the same scope.
# `strepto' and `nocardia' hold a single genus each, and a minimizer-based tool cannot be built to
# that scope in a comparable way -- which is itself a finding and belongs in the discussion, not in
# a script. `viral' therefore carries the comparison, even though it is where the refinement gains
# the least.
#
# WHAT MAKES THE THREE COMPARABLE. All of them are built from the very same sequences: the
# per-accession FASTA files the goal `extractrefseqfasta' wrote into data/projects/viral/fasta, whose
# headers carry the tax id in Kraken's own form, `>ACCESSION|kraken:taxid|TAXID'. No download, no
# other release, no second selection of genomes. The taxonomy comes from data/common as well, so all
# three see the same tree.
#
# Masking is off (`--no-masking'), as in the first study on Genestrip. Kraken 2 would otherwise run
# dustmasker over the library and drop low-complexity stretches that Genestrip keeps, which would
# make the comparison one of masking policies rather than of classification.
#
# Usage:
#   sh ./bin/kraken_build.sh           # both databases
#   sh ./bin/kraken_build.sh k2        # only Kraken 2
#   sh ./bin/kraken_build.sh ku        # only KrakenUniq
#
# Environment:
#   KRAKEN_PROJECT   the Genestrip project whose FASTA files are used, default `viral'
#   THREADS          build threads, default the number of CPUs
#
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)

what=${1:-all}
project=${KRAKEN_PROJECT:-viral}
threads=${THREADS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}

fastadir="${basedir}/data/projects/${project}/fasta"
common="${basedir}/data/common"
krakendir="${basedir}/data/kraken"
k2db="${krakendir}/${project}_k2"
kudb="${krakendir}/${project}_ku"
k2bin="${basedir}/tools/kraken2/bin"
kubin="${basedir}/tools/krakenuniq/bin"

if [ ! -d "$fastadir" ] || [ -z "$(ls -A "$fastadir" 2>/dev/null)" ]; then
  echo "No FASTA files in ${fastadir}." >&2
  echo "Run 'mvn exec:exec@db -Dname=${project} -Dgoal=extractrefseqfasta' first." >&2
  exit 1
fi
for f in nodes.dmp names.dmp; do
  [ -f "${common}/${f}" ] || { echo "Missing ${common}/${f}." >&2; exit 1; }
done

# Both tools want the taxonomy inside their database directory, and both read the same two files.
# Linking rather than copying: the dumps are large and neither tool writes to them.
prepare_taxonomy() {
  mkdir -p "$1/taxonomy"
  for f in nodes.dmp names.dmp; do
    [ -e "$1/taxonomy/${f}" ] || ln -s "${common}/${f}" "$1/taxonomy/${f}"
  done
}

case "$what" in
  k2|all)
    echo "############ Kraken 2: ${k2db} ############"
    [ -x "${k2bin}/kraken2-build" ] || { echo "Kraken 2 is missing - run ./bin/install_tools.sh first." >&2; exit 1; }
    if [ -f "${k2db}/hash.k2d" ]; then
      echo "SKIP  ${k2db} exists"
    else
      mkdir -p "$k2db"
      prepare_taxonomy "$k2db"
      # One call per file, as the first study did. kraken2-build reads the tax id from the header,
      # so no accession-to-taxid map is needed and nothing is fetched from the NCBI.
      count=0
      for f in "$fastadir"/*.fa; do
        "${k2bin}/kraken2-build" --no-masking --add-to-library "$f" --db "$k2db" >/dev/null
        count=$((count + 1))
      done
      echo "  added ${count} sequences"
      "${k2bin}/kraken2-build" --no-masking --build --threads "$threads" --db "$k2db"
    fi
    ;;
esac

case "$what" in
  ku|all)
    echo "############ KrakenUniq: ${kudb} ############"
    [ -x "${kubin}/krakenuniq-build" ] || { echo "KrakenUniq is missing - run ./bin/install_tools.sh first." >&2; exit 1; }
    # krakenuniq-build counts distinct k-mers with Jellyfish and takes its path from the environment.
    # It is installed beside the KrakenUniq sources, so it is looked up there unless already set.
    if [ -z "${JELLYFISH_BIN:-}" ]; then
      JELLYFISH_BIN=$(find "${basedir}/tools/krakenuniq" -type f -name jellyfish -perm -u+x 2>/dev/null | head -1)
      export JELLYFISH_BIN
    fi
    [ -n "${JELLYFISH_BIN:-}" ] && [ -x "$JELLYFISH_BIN" ] || {
      echo "No Jellyfish under ${basedir}/tools/krakenuniq; krakenuniq-build needs it in JELLYFISH_BIN." >&2
      echo "Re-run ./bin/install_tools.sh, which now keeps the sources and the jellyfish-install beside them." >&2
      exit 1; }
    echo "  Jellyfish: ${JELLYFISH_BIN}"
    if [ -f "${kudb}/database.kdb" ]; then
      echo "SKIP  ${kudb} exists"
    else
      mkdir -p "${kudb}/library"
      prepare_taxonomy "$kudb"
      # KrakenUniq takes the library as files under library/ together with a seqid-to-taxid map. Both
      # come out of the same headers, so the map is generated here rather than by a goal of its own.
      cat "$fastadir"/*.fa > "${kudb}/library/${project}.fa"
      awk '/^>/ { name = substr($1, 2); split(name, p, "|"); print name "\t" p[3] }' \
          "${kudb}/library/${project}.fa" > "${kudb}/seqid2taxid.map"
      echo "  library of $(grep -c '^>' "${kudb}/library/${project}.fa") sequences"
      "${kubin}/krakenuniq-build" --db "$kudb" --kmer-len 31 --threads "$threads" --taxids-for-genomes
    fi
    ;;
esac

case "$what" in
  k2|ku|all) ;;
  *) echo "Usage: $0 [k2|ku|all]" >&2; exit 1 ;;
esac

echo
du -sh "$k2db" "$kudb" 2>/dev/null || true
