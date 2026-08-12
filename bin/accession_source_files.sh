#
# Writes, for every sequence of a database's species, which source file it came from.
#
# The refined taxonomy of a Genestrip-FT database carries artificial FILE nodes, one per source
# fasta, and for a RefSeq-built database those are the release's `bacteria.N.1.genomic.fna.gz'
# chunks rather than single genomes. A chunk holds whatever sequences happen to sit in it, so
# "which genomes are below this node" - the question every sub-species measure asks - cannot be
# answered from the database alone. It also cannot be answered from the RefSeq catalog: its fourth
# column names the directory (`bacteria|complete') and not the file. What does answer it is the
# fasta headers themselves, which this script reads once.
#
# The result joins the two sides of a sub-species evaluation:
#
#   accession -> source file   (here)          accession -> sequence type   (mlst_assemblies.sh)
#
# and hence file node -> set of sequence types, from which every inner node follows by union along
# the tree that `ftdbinfo' prints.
#
# Only the files the database actually used are read, taken from its ftdbinfo CSV - for cdiff that
# is 465 of the release's chunks, some 143 GB compressed, which is an hour of streaming rather than
# the day a full rescan would take. Nothing is decompressed to disk; only header lines are kept.
#
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)

db=${1:-cdiff}
# The organism as it appears in the fasta headers, as an extended regular expression. A species
# travels under more than one name in the archives - C. difficile was Clostridium difficile and
# briefly Peptoclostridium difficile - and a header written years ago keeps the name of its day.
pattern=${2:-"Clostridioides difficile|Clostridium difficile|Peptoclostridium difficile"}

res_path="${basedir}/results"
info="${res_path}/${db}_ftdbinfo.csv"
out="${res_path}/${db}_accession_files.csv"
mkdir -p "$res_path"

if [ ! -f "$info" ]; then
  echo "Missing ${info}." >&2
  echo "Produce it with: mvn exec:exec@db -Dname=${db} -Dgoal=ftdbinfo" >&2
  echo "and copy the CSV from the project's csv folder into ${res_path}." >&2
  exit 1
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
mkdir -p "${work}/hdr"

# ---- which files the database used --------------------------------------------------------------
# ftdbinfo writes one row per node as `pos;level;name;rank;taxid;stored kmers;...'; the FILE rows
# carry the source file's name in the name column.
wanted="${work}/wanted.txt"
awk -F';' '$4 == "FILE" { print $3 }' "$info" | LC_ALL=C sort -u > "$wanted"
count=$(wc -l < "$wanted" | tr -d ' ')
if [ "$count" -eq 0 ]; then
  echo "No FILE nodes in ${info}, so there is nothing to resolve." >&2
  echo "That database was built without 'fileNodes=true'." >&2
  exit 1
fi
echo "=== ${count} source file(s) named by the taxonomy of ${db} ==="

# Where those files live. Searched rather than assumed: RefSeq chunks and GenBank assemblies sit in
# different folders, and a project may add fastas of its own.
index="${work}/index.txt"
# The whole of data/common rather than the folders one would name: Genestrip keeps the RefSeq
# release, the Genbank assemblies and any hand-placed fastas in separate folders below it, and
# naming them one by one is how the Genbank ones got missed. -L so that a data folder symlinked onto
# a larger disk - the normal arrangement at these sizes - is followed rather than passed over.
find -L "${basedir}/data/common" "${basedir}/data/projects/${db}/fasta" \
     -maxdepth 2 -type f 2>/dev/null > "$index" || true

# What a previous run already resolved is kept rather than read again: the scan is bound by
# streaming hundreds of gigabytes, and a run that had to be repeated - because a folder was missed,
# say - should only pay for what it missed. Set RESCAN=1 to read everything afresh.
skipped=0
if [ -f "$out" ] && [ -z "${RESCAN:-}" ]; then
  awk -F';' 'NR > 1' "$out" > "${work}/hdr/00-previous.csv"
  awk -F';' 'NR > 1 { print $2 }' "$out" | LC_ALL=C sort -u > "${work}/done.txt"
  skipped=$(wc -l < "${work}/done.txt" | tr -d ' ')
  if [ "$skipped" -gt 0 ]; then
    echo "  ${skipped} file(s) already in ${out}, keeping those (RESCAN=1 to read them again)"
    LC_ALL=C comm -23 "$wanted" "${work}/done.txt" > "${work}/wanted.new"
    mv "${work}/wanted.new" "$wanted"
  fi
fi

paths="${work}/paths.txt"
: > "$paths"
missing=0
while IFS= read -r name; do
  # A plain grep would match a file whose name merely ends with the wanted one, so the basename is
  # compared as a whole.
  path=$(awk -v n="$name" 'BEGIN { FS = "/" } $NF == n { print; exit }' "$index")
  if [ -n "$path" ]; then
    printf '%s\n' "$path" >> "$paths"
  else
    missing=$((missing + 1))
    echo "  WARNING: no file named ${name} below data/common or the project's fasta folder" >&2
  fi
done < "$wanted"
found=$(wc -l < "$paths" | tr -d ' ')
echo "OK    ${found} to read, ${missing} missing"
if [ "$found" -eq 0 ]; then
  if [ "$skipped" -gt 0 ]; then
    echo "  Nothing left to read; ${out} already covers every file that could be found."
    if [ "$missing" -gt 0 ]; then
      echo "  ${missing} file(s) remain unresolved - see the warnings above." >&2
    fi
    exit 0
  fi
  echo "  None of the ${count} file(s) named by the taxonomy could be found." >&2
  exit 1
fi

# ---- read their headers -------------------------------------------------------------------------
# One job per file, several at a time: this is bound by reading and decompressing, so it scales with
# whatever the machine will give. Each job writes a file of its own; nothing is appended to a shared
# one, which would interleave.
jobs=${SCAN_JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}
cat > "${work}/scan.sh" <<'EOF'
f=$1
name=$(basename "$f")
# `gzip -cdf' decompresses what is compressed and passes through what is not, so the same command
# serves .fna.gz and plain .fna. grep finding nothing is a result, not a failure, hence the guard.
gzip -cdf "$f" \
  | grep '^>' \
  | grep -E "$SCAN_PATTERN" \
  | awk -v file="$name" '{ acc = substr($1, 2); print acc ";" file ";" }' \
  > "${SCAN_WORK}/hdr/${name}.csv" || true
EOF
export SCAN_PATTERN="$pattern"
export SCAN_WORK="$work"
echo "=== reading the headers of ${found} file(s), ${jobs} at a time ==="
echo "  this streams the whole of them; only header lines are kept"
xargs -P "$jobs" -n 1 sh "${work}/scan.sh" < "$paths"

{
  echo "accession;file;"
  # Sorted so that the result is the same whichever order the jobs happened to finish in.
  cat "${work}/hdr/"*.csv 2>/dev/null | LC_ALL=C sort
} > "$out"

sequences=$(( $(wc -l < "$out" | tr -d ' ') - 1 ))
files=$(awk -F';' 'NR > 1 { f[$2] } END { print length(f) }' "$out")
echo
echo "Wrote ${out}"
echo "  ${sequences} sequence(s) of '${db}' across ${files} source file(s)."
if [ "$sequences" -eq 0 ]; then
  echo "  Nothing matched. The organism appears in a fasta header under a name this pattern does" >&2
  echo "  not cover; pass the right one as the second argument, e.g." >&2
  echo "    sh ./bin/accession_source_files.sh ${db} 'Clostridioides difficile'" >&2
  exit 1
fi
echo "  Join it to results/${db}_mlst.csv on the accession (the mlst file column is the accession"
echo "  followed by '.fa'), and to ${db}_ftdbinfo.csv on the file name, to obtain the sequence"
echo "  types below each node of the refined taxonomy."
