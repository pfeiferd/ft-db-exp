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
  echo "The genomes a database was filled from are not in that folder until they are extracted" >&2
  echo "there, which bin/cdiff_eval.sh does, or by hand:" >&2
  echo "  mvn exec:exec@db -Dname=${db} -Dgoal=extractrefseqcsv" >&2
  echo "Ask for that goal rather than for 'extractrefseqfasta': the latter is an ObjectGoal, which" >&2
  echo "Maker.make() skips as a weak dependency, so it reports success within a second and does" >&2
  echo "nothing at all." >&2
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

# RefSeq deposits a draft genome as hundreds of WGS contigs, and the extraction writes one file per
# accession. A seven-locus scheme is practically never complete on a single contig, so typing the
# files one by one leaves everything but the finished chromosomes untyped -- on the C. difficile
# database that was 640 of 463906 files. The contigs of one assembly are recognisable by their
# accession, which is <letters><2-digit assembly version><contig number>, and are concatenated and
# typed together. That is both a far better ground truth and far less work: 13184 typings instead
# of 463906.
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
lists="${work}/lists"
asmdir="${work}/asm"
mkdir -p "$lists" "$asmdir"

echo "=== grouping the files of ${db} into assemblies ==="
# find rather than a glob: some projects hold hundreds of thousands of files, which a glob expands
# into one argument list.
members="${work}/members.tsv"
find "$fastadir" -maxdepth 1 -type f \( -name '*.fna.gz' -o -name '*.fa.gz' -o -name '*.fasta.gz' \
        -o -name '*.fna' -o -name '*.fa' -o -name '*.fasta' \) \
  | LC_ALL=C sort \
  | awk -v lists="$lists" '
      {
        path = $0
        n = split(path, parts, "/"); name = parts[n]
        acc = name
        sub(/\.(fna|fa|fasta)(\.gz)?$/, "", acc)
        sub(/\.[0-9]+$/, "", acc)              # accession version
        core = acc
        sub(/^[A-Z][A-Z]_/, "", core)          # RefSeq NZ_ / NC_ prefix
        key = acc
        if (match(core, /^[A-Z]+[0-9]+$/)) {
          letters = core; sub(/[0-9].*$/, "", letters)
          digits = substr(core, length(letters) + 1)
          # Two digits of assembly version plus at least six of contig number make a WGS contig.
          # Anything shorter is an accession in its own right: a finished chromosome or a plasmid.
          if (length(digits) >= 8) {
            key = letters substr(digits, 1, 2)
          }
        }
        # The key doubles as a file name, so it must not carry a path separator. Accessions do not,
        # but a stray one would silently write outside the work directory.
        gsub(/[^A-Za-z0-9._-]/, "_", key)
        print path >> (lists "/" key)
        close(lists "/" key)
        print key "\t" name
      }' > "$members"

count=$(wc -l < "$members" | tr -d ' ')
if [ "$count" -eq 0 ]; then
  echo "No fasta files in ${fastadir}." >&2
  exit 1
fi
asmcount=$(find "$lists" -type f | wc -l | tr -d ' ')
echo "OK    ${count} file(s) in ${asmcount} assembl(y/ies)"

# mlst reads plain fasta; `gzip -cdf' decompresses what is compressed and copies through what is
# not, so one invocation per assembly concatenates its contigs whatever they are stored as.
#
# The assemblies are typed in batches rather than one per invocation: mlst loads the scheme's allele
# database on every start, which for a single genome costs more than the typing itself. The batch is
# also what bounds the extra disk -- only the assemblies of the current batch exist unpacked.
BATCH=${MLST_BATCH:-50}
threadopt=""
if mlst --help 2>&1 | grep -q -- "--threads"; then
  threadopt="--threads $(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"
fi

rawout="${work}/mlst.tsv"
: > "$rawout"
batch=""
nbatch=0
typedsofar=0

run_batch() {
  [ "$nbatch" -gt 0 ] || return 0
  # A batch that dies takes no more than its own assemblies with it: they are recorded without an
  # ST, exactly as an assembly mlst could not type, and the run goes on.
  # Not named `out': that is the result CSV this script writes at the end.
  if batchout=$(mlst --quiet --scheme "$scheme" $threadopt $batch 2>/dev/null); then
    printf '%s\n' "$batchout" >> "$rawout"
  else
    echo "  WARNING: mlst failed on a batch of ${nbatch} assembl(y/ies)" >&2
    for a in $batch; do
      printf '%s\t%s\t-\n' "$a" "$scheme" >> "$rawout"
    done
  fi
  rm -f $batch
  batch=""
  nbatch=0
}

echo "=== typing ${asmcount} assembl(y/ies) of ${db} under scheme ${scheme} ==="
for lst in "$lists"/*; do
  [ -e "$lst" ] || continue
  key=${lst##*/}
  asmfile="${asmdir}/${key}.fa"
  xargs gzip -cdf < "$lst" > "$asmfile"
  batch="${batch} ${asmfile}"
  nbatch=$((nbatch + 1))
  typedsofar=$((typedsofar + 1))
  if [ "$nbatch" -ge "$BATCH" ]; then
    run_batch
    echo "  ${typedsofar}/${asmcount}"
  fi
done
run_batch
echo "  ${typedsofar}/${asmcount}"

# Semicolon-separated, as every other CSV this project writes, so that the paper's \csvreader and
# the reports read it the same way. The allele columns are kept whole in one field: their number
# varies with the scheme, and a fixed header cannot describe them.
#
# One row per file, not per assembly: the join key of everything downstream is the accession, which
# is what Genestrip knows a genome by. Every contig of an assembly therefore carries the ST of the
# assembly it belongs to, and the assembly it was typed as is named beside it.
{
  echo "file;assembly;scheme;st;alleles;"
  # mlst leaves the ST as `-' when the profile is novel or incomplete. That is a result, not a
  # failure, and it is passed through rather than dropped: a genome whose ST is unknown still
  # belongs in the table, and a cluster made only of such genomes is itself worth seeing.
  # `scheme' has to be handed to awk explicitly; it is a shell variable and would otherwise be empty
  # inside, which would leave the scheme column blank for an assembly missing from mlst's output.
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
      key = $1
      print $2, key, (key in sch ? sch[key] : scheme), (key in st ? st[key] : "-"), all[key], ""
    }' "$rawout" "$members"
} > "$out"

typed=$(awk -F';' 'NR > 1 && $4 != "-" && $4 != "" { n++ } END { print n + 0 }' "$out")
novel=$(awk -F';' 'NR > 1 && $4 == "-" { n++ } END { print n + 0 }' "$out")
asmtyped=$(awk -F';' 'NR > 1 && $4 != "-" && $4 != "" { print $2 }' "$out" | sort -u | wc -l | tr -d ' ')
sts=$(awk -F';' 'NR > 1 && $4 != "-" && $4 != "" { print $4 }' "$out" | sort -u | wc -l | tr -d ' ')
echo
echo "Wrote ${out}"
echo "  ${asmtyped} of ${asmcount} assembl(y/ies) typed, ${sts} distinct ST(s)."
echo "  ${typed} sequence file(s) carry an ST, ${novel} do not."
if [ "$novel" -gt 0 ]; then
  echo "  An ST of '-' means a novel or incomplete profile, or a plasmid that was typed on its own." >&2
  echo "  Those genomes are kept in the table; they cannot contribute a candidate ST to a" >&2
  echo "  sub-species measure and have to be counted as their own unit or excluded deliberately --" >&2
  echo "  not ignored by accident." >&2
fi
