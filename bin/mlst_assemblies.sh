#!/bin/sh
#
# Assigns a multi-locus sequence type to every genome that became a leaf of the `cdiff' database,
# and writes the result to results/cdiff_mlst.csv.
#
#   sh ./bin/mlst_assemblies.sh                    # the cdiff project
#   sh ./bin/mlst_assemblies.sh <project> <scheme> # any other project and PubMLST scheme
#
# WHY THIS EXISTS. The refinement of `cdiff' clusters assemblies by k-mer similarity alone, and the
# leaves of the resulting dendrogram are fasta files rather than named lineages. Nothing in that
# construction knows what a cluster *is*. An ST per assembly supplies the missing half in two ways:
#
#   - it makes the clusters interpretable, and lets them be checked against an established typing
#     scheme instead of against the refinement itself;
#   - it is the ground truth of a sub-species precision measure. Below the species there is no rank
#     to count candidates at, and the isolates of PRJNA1148956 are not in the database anyway --
#     only their reads were deposited -- so a read has no true node to be scored against. Scoring
#     against the ST of the isolate, and counting candidates as the distinct STs still in question
#     at the assigned node, is what makes the measure definable at all. Without this file there is
#     no gated measure below the species.
#
# WHAT GETS TYPED, AND WHY IT IS READ OUT OF THE DATABASE. The unit of this table has to be the unit
# of the dendrogram, which is a leaf: `fileNodes=true' puts one artificial node of rank FILE under
# the tax id a genome's k-mers map to, named after the fasta file it was read from
# (`TaxTree.fileNode(node, file.getName(), ...)'). So the list of genomes to type is taken from
# `<db>_dbinfo.csv' -- every FILE node inside the subtree of the requested tax ids -- and not from a
# folder. Two things follow that a folder cannot give:
#
#   - the join key is right by construction. A row is keyed by the leaf's name, which *is* the file
#     name, `.gz' and all, so the table joins to the dendrogram with nothing to reconcile.
#   - nothing is typed that is not in the tree, and nothing in the tree is missed. `genbank.maxPerTaxid'
#     admits a subset of the assemblies and the database is the only place that records which subset,
#     while `additional.txt' contributes leaves of its own -- for `cdiff' the human decoy under 9606,
#     which the subtree filter drops because it is not below the requested species.
#
# One name may be several leaves. `fileNode()' creates its node under the tax id the region resolved
# to, so a genome read under two of them -- the species and one of its strains, say -- becomes a
# FILE node under each, both named after the same file. The names are therefore made unique before
# typing: it is one genome and it has one sequence type, and since everything downstream joins on
# the name, one row serves every leaf that carries it.
#
# THIS REPLACES A DETOUR, and the reason it does is worth keeping. The genomes used to be taken from
# `data/projects/<db>/fasta', which `extractrefseqfasta' fills with one file per sequence accession:
# a RefSeq release file is a chunk of many organisms, so it cannot be typed as it stands. Since a
# seven-locus profile is practically never complete on a single WGS contig, those files then had to
# be regrouped into assemblies by parsing accessions, and the resulting rows were keyed by an
# accession stem -- which is not what a leaf is called, so the join to the dendrogram silently
# matched nothing. With one Genbank fasta per assembly there is no chunk to unpack and no group to
# reconstruct: one file is one assembly is one leaf.
#
# VALIDATE THE PIPELINE BEFORE TRUSTING IT. The source study typed its own 37 isolates; assembling
# their Illumina reads and running this over the result must reproduce the STs of its Table 2. If
# it does, the same command is trustworthy over the database's genomes -- and if it does not, that
# is much better learned here than after a dendrogram has been interpreted.
#
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)

db=${1:-cdiff}
# tseemann/mlst's name for the PubMLST C. difficile scheme. Confirm with `mlst --list'; naming it
# explicitly rather than letting the tool guess keeps a genome from being typed under the scheme of
# some other species that happens to score better on a poor assembly.
scheme=${2:-cdifficile}

res_path="${basedir}/results"
out="${res_path}/${db}_mlst.csv"
mkdir -p "$res_path"

# Where a leaf's fasta may be found, in the order searched. Genbank first, since that is where the
# per-assembly genomes come from; the others cover a project that also lists fastas of its own.
# Held in the positional parameters rather than in a variable, which is the one list a POSIX shell
# has that survives a path with a space in it. The two arguments this script takes have been read
# into `db' and `scheme' above, so nothing is lost by overwriting them.
#
# data/common/refseq is deliberately NOT among them. A release file is a chunk of many organisms, so
# a leaf named after one is not a genome and typing it would produce a sequence type for nothing in
# particular. Leaving the folder out makes such a leaf turn up as unresolved, where it is reported,
# rather than as a row that looks like a result.
set -- "${basedir}/data/common/genbank" \
       "${basedir}/data/common/fasta" \
       "${basedir}/data/projects/${db}/fasta"

# install_tools.sh puts mlst below tools/bin, and its own `export PATH' lives no longer than that
# script does. Prepending the directory here is what make_fastqs.sh does with the simulators for the
# same reason, and it leaves an mlst installed system-wide in charge, since that one comes first if
# tools/bin holds none.
export PATH="${basedir}/tools/bin:${PATH}"

if ! command -v mlst >/dev/null 2>&1; then
  echo "mlst is not on the PATH, and ${basedir}/tools/bin holds none either." >&2
  echo "Run ./bin/install_tools.sh first - it clones https://github.com/tseemann/mlst there and" >&2
  echo "installs the blast+ it drives. A failure of that step is only a warning, so check its" >&2
  echo "output for '5/5  mlst' if it appeared to succeed." >&2
  exit 1
fi

# The database's own inventory of nodes. Written by `dbinfo' into the project's csv folder, and
# copied to results/ by run_exps.sh afterwards - either will do, so both are looked at.
dbinfo=""
for cand in "${basedir}/data/projects/${db}/csv/${db}_dbinfo.csv" "${res_path}/${db}_dbinfo.csv"; do
  if [ -f "$cand" ]; then
    dbinfo="$cand"
    break
  fi
done
if [ -z "$dbinfo" ]; then
  echo "No ${db}_dbinfo.csv, so there is no list of what the database's leaves are." >&2
  echo "Build it from the database (which must exist) with:" >&2
  echo "  mvn exec:exec@db -Dname=${db} -Dgoal=dbinfo" >&2
  exit 1
fi

taxidfile="${basedir}/data/projects/${db}/taxids.txt"
if [ ! -f "$taxidfile" ]; then
  echo "Missing ${taxidfile}, so the subtree to type cannot be determined." >&2
  exit 1
fi

# The scheme has to exist, or mlst silently types everything as `-' and the output looks like a
# result. `mlst --list' prints them space-separated on one line.
#
# The two ways this can go wrong need telling apart, because they are fixed in entirely different
# places: mlst may know no scheme at all, which means its installation is incomplete, or it may know
# plenty while none of them carries the name we ask for. Reporting the number it knows separates the
# two at a glance, and its stderr is passed on rather than dropped, since a broken installation
# explains itself there and nowhere else.
mlst_stderr=$(mktemp)
schemes=$(mlst --list 2>"$mlst_stderr" | tr ' ' '\n' | grep -v '^[[:space:]]*$' | sort)
scheme_count=$(printf '%s\n' "$schemes" | grep -c . || true)

if [ "$scheme_count" -eq 0 ]; then
  echo "mlst knows no schemes at all, so its installation is incomplete." >&2
  echo "It keeps them in db/pubmlst below its own directory; check that this one has them:" >&2
  echo "  ls \"\$(dirname \"\$(readlink -f \"\$(command -v mlst)\")\")/../db/pubmlst\" | head" >&2
  echo "and that it can run at all:" >&2
  echo "  mlst --check" >&2
  if [ -s "$mlst_stderr" ]; then
    echo "What 'mlst --list' reported:" >&2
    sed 's/^/  /' "$mlst_stderr" >&2
  fi
  rm -f "$mlst_stderr"
  exit 1
fi

if ! printf '%s\n' "$schemes" | grep -qx "$scheme"; then
  echo "mlst does not know a scheme called '${scheme}', though it knows ${scheme_count} others." >&2
  candidates=$(printf '%s\n' "$schemes" | grep -i "diff\|clostr" || true)
  if [ -n "$candidates" ]; then
    echo "These look related:" >&2
    printf '%s\n' "$candidates" | sed 's/^/  /' >&2
  else
    echo "None of them looks related to C. difficile. The full list starts with:" >&2
    printf '%s\n' "$schemes" | head -20 | sed 's/^/  /' >&2
    echo "  ... (see 'mlst --list' for all ${scheme_count})" >&2
  fi
  echo "Pass the right one as the second argument." >&2
  rm -f "$mlst_stderr"
  exit 1
fi
rm -f "$mlst_stderr"

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
asmdir="${work}/asm"
mkdir -p "$asmdir"

echo "=== reading the leaves of ${db} from $(basename "$dbinfo") ==="
# dbinfo lists the tree in pre-order with a depth column, so a subtree is the run of rows that
# follows its root and is deeper than it. Entering at a requested tax id and leaving at the first
# row no deeper than the one entered at therefore delimits exactly what belongs to it - which is
# what keeps the human decoy of additional.txt, a FILE node under 9606, out of a C. difficile
# typing. The TOTAL row is a summary and carries a checksum where a tax id belongs, hence skipped.
leaves="${work}/leaves.txt"
awk -F';' -v taxidfile="$taxidfile" '
    BEGIN {
      while ((getline line < taxidfile) > 0) {
        sub(/#.*$/, "", line)
        gsub(/[[:space:]]/, "", line)
        if (line != "") {
          want[line] = 1
        }
      }
    }
    NR == 1 || $3 == "TOTAL" { next }
    {
      level = $2 + 0
      if (inside && level <= rootlevel) {
        inside = 0
      }
      if (!inside) {
        if ($5 in want) {
          inside = 1
          rootlevel = level
        }
        next
      }
      if ($4 == "FILE") {
        print $3
      }
    }' "$dbinfo" | LC_ALL=C sort -u > "$leaves"

leafcount=$(wc -l < "$leaves" | tr -d ' ')
if [ "$leafcount" -eq 0 ]; then
  echo "No FILE nodes below the requested tax ids of ${db}." >&2
  echo "Either the database was built with 'fileNodes=false', in which case there is nothing to" >&2
  echo "type here, or ${dbinfo} is from an older build than the database beside it." >&2
  exit 1
fi
echo "OK    ${leafcount} leaf/leaves"

# Each leaf is named after the file it was read from, so resolving it is a lookup over the search
# path rather than a parse. A name that resolves nowhere is reported and not passed over: it means
# the database was filled from something this script cannot see, and a ground truth quietly missing
# those genomes is worse than one that says so.
echo "=== locating their fasta files ==="
resolved="${work}/resolved.tsv"
missing="${work}/missing.txt"
: > "$resolved"
: > "$missing"
while IFS= read -r leaf; do
  found=""
  for dir in "$@"; do
    if [ -f "${dir}/${leaf}" ]; then
      found="${dir}/${leaf}"
      break
    fi
  done
  if [ -n "$found" ]; then
    printf '%s\t%s\n' "$leaf" "$found" >> "$resolved"
  else
    printf '%s\n' "$leaf" >> "$missing"
  fi
done < "$leaves"

foundcount=$(wc -l < "$resolved" | tr -d ' ')
missingcount=$(wc -l < "$missing" | tr -d ' ')
echo "OK    ${foundcount} located, ${missingcount} not found"
if [ "$missingcount" -gt 0 ]; then
  echo "  These leaves have no fasta on the search path, so they cannot be typed:" >&2
  head -10 "$missing" | sed 's/^/    /' >&2
  [ "$missingcount" -gt 10 ] && echo "    ... and $((missingcount - 10)) more" >&2
  echo "  Searched:" >&2
  printf '%s\n' "$@" | sed 's/^/    /' >&2
  # The one cause worth naming, because it is not a missing file but a database built the other way.
  chunks=$(grep -c '^[a-z_]*\.[0-9][0-9]*\.[0-9][0-9]*\.genomic\.fna\(\.gz\)\?$' "$missing" || true)
  if [ "$chunks" -gt 0 ]; then
    echo "  ${chunks} of them are named after RefSeq release files, so this database was filled" >&2
    echo "  from the release ('refseq.filldb' is not false) and those leaves are chunks of many" >&2
    echo "  organisms rather than genomes. They cannot carry a sequence type at all, and a" >&2
    echo "  dendrogram over them clusters by accession order as much as by lineage - see" >&2
    echo "  data/projects/cdiff/config.properties on why cdiff switched the release off." >&2
  fi
fi
if [ "$foundcount" -eq 0 ]; then
  echo "Nothing to type." >&2
  exit 1
fi

# mlst reads plain fasta; `gzip -cdf' decompresses what is compressed and copies through what is
# not, so one invocation per genome unpacks it whatever it is stored as.
#
# The genomes are typed in batches rather than one per invocation: mlst loads the scheme's allele
# database on every start, which for a single genome costs more than the typing itself. The batch is
# also what bounds the extra disk -- only the genomes of the current batch exist unpacked.
#
# The unpacked file is named by a serial number and not after its leaf, with a map beside it. A leaf
# name is a file name of someone else's choosing, and letting it name a file again - after being
# concatenated with a suffix, passed through the shell unquoted into an mlst argument list and split
# back off the tool's output - is a chain of assumptions about characters it may contain. A serial
# number assumes nothing and the map is exact.
BATCH=${MLST_BATCH:-50}
threadopt=""
if mlst --help 2>&1 | grep -q -- "--threads"; then
  threadopt="--threads $(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"
fi

rawout="${work}/mlst.tsv"
map="${work}/map.tsv"
: > "$rawout"
: > "$map"
batch=""
nbatch=0
typedsofar=0
serial=0

run_batch() {
  [ "$nbatch" -gt 0 ] || return 0
  # A batch that dies takes no more than its own genomes with it: they are recorded without an ST,
  # exactly as a genome mlst could not type, and the run goes on.
  # Not named `out': that is the result CSV this script writes at the end.
  if batchout=$(mlst --quiet --scheme "$scheme" $threadopt $batch 2>/dev/null); then
    printf '%s\n' "$batchout" >> "$rawout"
  else
    echo "  WARNING: mlst failed on a batch of ${nbatch} genome(s)" >&2
    for a in $batch; do
      printf '%s\t%s\t-\n' "$a" "$scheme" >> "$rawout"
    done
  fi
  rm -f $batch
  batch=""
  nbatch=0
}

echo "=== typing ${foundcount} genome(s) of ${db} under scheme ${scheme} ==="
while IFS="$(printf '\t')" read -r leaf path; do
  serial=$((serial + 1))
  id=$(printf 'g%08d' "$serial")
  printf '%s\t%s\n' "$id" "$leaf" >> "$map"
  asmfile="${asmdir}/${id}.fa"
  gzip -cdf "$path" > "$asmfile"
  batch="${batch} ${asmfile}"
  nbatch=$((nbatch + 1))
  typedsofar=$((typedsofar + 1))
  if [ "$nbatch" -ge "$BATCH" ]; then
    run_batch
    echo "  ${typedsofar}/${foundcount}"
  fi
done < "$resolved"
run_batch
echo "  ${typedsofar}/${foundcount}"

# Semicolon-separated, as every other CSV this project writes, so that the paper's \csvreader and
# the reports read it the same way. The allele columns are kept whole in one field: their number
# varies with the scheme, and a fixed header cannot describe them.
#
# One row per leaf, and the leaf's name is the key: that is what the dendrogram calls the genome, so
# everything downstream joins on it directly.
{
  echo "leaf;scheme;st;alleles;"
  # mlst leaves the ST as `-' when the profile is novel or incomplete. That is a result, not a
  # failure, and it is passed through rather than dropped: a genome whose ST is unknown still
  # belongs in the table, and a cluster made only of such genomes is itself worth seeing.
  # `scheme' has to be handed to awk explicitly; it is a shell variable and would otherwise be empty
  # inside, which would leave the scheme column blank for a genome missing from mlst's output.
  awk -F'\t' -v OFS=';' -v scheme="$scheme" '
    NR == FNR {
      n = split($1, parts, "/"); key = parts[n]
      sub(/\.fa$/, "", key)
      sch[key] = $2
      st[key] = ($3 == "" ? "-" : $3)
      alleles = ""
      for (i = 4; i <= NF; i++) {
        alleles = (i == 4 ? $i : alleles "," $i)
      }
      all[key] = alleles
      next
    }
    {
      id = $1
      print $2, (id in sch ? sch[id] : scheme), (id in st ? st[id] : "-"), all[id], ""
    }' "$rawout" "$map"
} > "$out"

typed=$(awk -F';' 'NR > 1 && $3 != "-" && $3 != "" { n++ } END { print n + 0 }' "$out")
novel=$(awk -F';' 'NR > 1 && $3 == "-" { n++ } END { print n + 0 }' "$out")
sts=$(awk -F';' 'NR > 1 && $3 != "-" && $3 != "" { print $3 }' "$out" | sort -u | wc -l | tr -d ' ')
echo
echo "Wrote ${out}"
echo "  ${typed} of ${leafcount} leaf/leaves carry an ST, in ${sts} distinct ST(s)."
if [ "$novel" -gt 0 ]; then
  echo "  ${novel} do not. An ST of '-' means a novel or incomplete profile." >&2
  echo "  Those genomes are kept in the table; they cannot contribute a candidate ST to a" >&2
  echo "  sub-species measure and have to be counted as their own unit or excluded deliberately --" >&2
  echo "  not ignored by accident." >&2
fi
if [ "$missingcount" -gt 0 ]; then
  echo "  ${missingcount} further leaf/leaves had no fasta at all and are absent from the table." >&2
fi
