#!/bin/sh
#
# Assigns a multi-locus sequence type to every genome the `cdiff' database was built from, and
# writes the result to results/cdiff_mlst.csv.
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
# THE JOIN IS FREE. `fileNodes=true' names each artificial node after the fasta file it came from
# (TaxTree.fileNode(node, file.getName(), ...)), and `mlst' reports one row per fasta file. So the
# ST table joins to the dendrogram's leaves on the file name, with nothing to reconcile by hand.
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

fastadir="${basedir}/data/projects/${db}/fasta"
res_path="${basedir}/results"
out="${res_path}/${db}_mlst.csv"
mkdir -p "$res_path"

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
if [ ! -d "$fastadir" ]; then
  echo "Missing ${fastadir}." >&2
  echo "Build the database first, so that Genestrip downloads the genomes:" >&2
  echo "  mvn exec:exec@db -Dname=${db} -Dgoal=db" >&2
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

# Genestrip stores the genomes gzipped; mlst reads plain fasta, and unpacking the whole folder just
# to type it would double the disk. Each file is therefore decompressed to a temporary copy, typed,
# and removed again.
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

count=0
files=""
for f in "$fastadir"/*.fna.gz "$fastadir"/*.fa.gz "$fastadir"/*.fasta.gz "$fastadir"/*.fna "$fastadir"/*.fa "$fastadir"/*.fasta; do
  [ -e "$f" ] || continue
  count=$((count + 1))
done
if [ "$count" -eq 0 ]; then
  echo "No fasta files in ${fastadir}." >&2
  exit 1
fi
echo "=== typing ${count} genome(s) of ${db} under scheme ${scheme} ==="

# One invocation per file rather than one over the whole folder: the folder is thousands of files
# for some projects, which overruns the argument list, and a single unreadable genome would take
# the whole run down with it instead of one row.
tmpout="${work}/mlst.tsv"
: > "$tmpout"
i=0
for f in "$fastadir"/*.fna.gz "$fastadir"/*.fa.gz "$fastadir"/*.fasta.gz "$fastadir"/*.fna "$fastadir"/*.fa "$fastadir"/*.fasta; do
  [ -e "$f" ] || continue
  i=$((i + 1))
  name=$(basename "$f")
  case "$f" in
    *.gz) plain="${work}/$(basename "$f" .gz)"; gunzip -c "$f" > "$plain" ;;
    *)    plain="$f" ;;
  esac
  # mlst prints: <file> <scheme> <ST> <locus(allele)> ... . The file column is the temporary path,
  # so it is replaced by the name Genestrip knows the genome by, which is the join key.
  if line=$(mlst --quiet --scheme "$scheme" "$plain" 2>/dev/null); then
    printf '%s\t%s\n' "$name" "$(printf '%s' "$line" | cut -f2-)" >> "$tmpout"
  else
    echo "  WARNING: mlst failed on ${name}" >&2
    printf '%s\t%s\t-\n' "$name" "$scheme" >> "$tmpout"
  fi
  case "$f" in *.gz) rm -f "$plain" ;; esac
  [ $((i % 25)) -eq 0 ] && echo "  ${i}/${count}"
done

# Semicolon-separated, as every other CSV this project writes, so that the paper's \csvreader and
# the reports read it the same way. The allele columns are kept whole in one field: their number
# varies with the scheme, and a fixed header cannot describe them.
{
  echo "file;scheme;st;alleles;"
  # mlst leaves the ST as `-' when the profile is novel or incomplete. That is a result, not a
  # failure, and it is passed through rather than dropped: a genome whose ST is unknown still
  # belongs in the table, and a cluster made only of such genomes is itself worth seeing.
  while IFS="$(printf '\t')" read -r name sch st rest; do
    printf '%s;%s;%s;%s;\n' "$name" "$sch" "$st" "$(printf '%s' "$rest" | tr '\t' ',')"
  done < "$tmpout"
} > "$out"

typed=$(awk -F';' 'NR > 1 && $3 != "-" && $3 != "" { n++ } END { print n + 0 }' "$out")
novel=$(awk -F';' 'NR > 1 && $3 == "-" { n++ } END { print n + 0 }' "$out")
sts=$(awk -F';' 'NR > 1 && $3 != "-" && $3 != "" { print $3 }' "$out" | sort -u | wc -l)
echo
echo "Wrote ${out}"
echo "  ${typed} genome(s) typed, ${novel} without an ST, ${sts} distinct ST(s)."
if [ "$novel" -gt 0 ]; then
  echo "  An ST of '-' means a novel or incomplete profile. Those genomes are kept in the table;" >&2
  echo "  they cannot contribute a candidate ST to a sub-species measure and have to be counted" >&2
  echo "  as their own unit or excluded deliberately -- not ignored by accident." >&2
fi
