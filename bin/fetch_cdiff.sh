#!/bin/sh
#
# Fetches the reads of the 37 Clostridioides difficile isolates of BioProject PRJNA1148956 -- the
# Illumina/Nanopore comparison of Bejaoui et al., BMC Genomics 2025 (doi:10.1186/s12864-025-11267-9)
# -- into data/fastq, named as the archive names them, which is what data/fastq/cdiff_isolates.txt
# refers to.
#
#   sh ./bin/fetch_cdiff.sh              # everything: 74 runs, 111 files, about 29 GB
#   sh ./bin/fetch_cdiff.sh illumina     # the 37 paired Illumina runs only (about 12 GB)
#   sh ./bin/fetch_cdiff.sh nanopore     # the 37 Nanopore runs only (about 17 GB)
#   sh ./bin/fetch_cdiff.sh B11 B12np    # named isolates, by the study's own identifiers
#   DRY_RUN=1 sh ./bin/fetch_cdiff.sh    # list what would be fetched, and how much, without doing it
#
# Unlike the saliva runs this needs no sra-tools: ENA serves the fastq files directly over https,
# so a plain download does it. The file list is not hard-coded but queried from ENA's portal API at
# run time, keyed by the BioProject -- so a run added to or corrected in the archive is picked up
# rather than silently missed, and this script cannot drift out of step with the map beside it.
#
# Whatever is already in place is left alone, so an interrupted run is resumed by repeating the
# command. Every file is verified with `gzip -t' after download and removed if it fails, since a
# host that answers with an error page would otherwise leave something that looks like data.
#
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)
fastqdir="${basedir}/data/fastq"
mkdir -p "$fastqdir"

project=${PROJECT:-PRJNA1148956}
what=${1:-all}

# The portal API returns one row per run: alias, platform and the semicolon-separated file paths.
# `fastq_ftp' gives host-relative paths without a scheme; https on the same host serves them.
api="https://www.ebi.ac.uk/ena/portal/api/filereport"
query="${api}?accession=${project}&result=read_run&fields=run_accession,sample_alias,instrument_platform,fastq_ftp,fastq_bytes&format=tsv"

echo "=== querying ENA for the runs of ${project} ==="
runs=$(mktemp)
trap 'rm -f "$runs"' EXIT
if ! curl -sSf -m 120 "$query" > "$runs"; then
  echo "Could not reach ENA's portal API." >&2
  echo "If this is a network policy block, allow www.ebi.ac.uk and ftp.sra.ebi.ac.uk." >&2
  exit 1
fi
# A single header line means the project resolved but holds no runs, which is a different problem
# from a failed request and deserves its own message rather than a silent zero-file success.
if [ "$(wc -l < "$runs")" -le 1 ]; then
  echo "ENA returned no runs for ${project}." >&2
  exit 1
fi
echo "OK    $(( $(wc -l < "$runs") - 1 )) runs"

# Selects the rows to fetch. With no argument everything; with `illumina'/`nanopore' one platform;
# otherwise the arguments are taken as isolate identifiers and matched against the alias column.
select_rows() {
  case "$what" in
    all)       awk -F'\t' 'NR > 1' "$runs" ;;
    illumina)  awk -F'\t' 'NR > 1 && $3 ~ /ILLUMINA/' "$runs" ;;
    nanopore)  awk -F'\t' 'NR > 1 && $3 ~ /NANOPORE/' "$runs" ;;
    *)         for a in "$@"; do awk -F'\t' -v a="$a" 'NR > 1 && $2 == a' "$runs"; done ;;
  esac
}

rows=$(select_rows "$@")
if [ -z "$rows" ]; then
  echo "Nothing selected by '${what}'. Use all, illumina, nanopore, or isolate ids like B11 B12np." >&2
  echo "Known aliases:" >&2
  awk -F'\t' 'NR > 1 { print "  " $2 }' "$runs" | sort | tr '\n' ' ' >&2
  echo >&2
  exit 1
fi

if [ -n "${DRY_RUN:-}" ]; then
  echo "$rows" | awk -F'\t' -v dir="$fastqdir" '
    { n = split($4, files, ";"); split($5, sizes, ";")
      for (i = 1; i <= n; i++) {
        if (files[i] == "") continue
        m = split(files[i], p, "/"); name = p[m]
        cmd = "test -s \"" dir "/" name "\""
        have = (system(cmd) == 0)
        printf "  %-5s %-28s %-8s %6.1f MB  %s\n", $2, name, ($3 ~ /ILLUMINA/ ? "illumina" : "nanopore"),
               sizes[i] / 1e6, have ? "(already here)" : ""
        total += sizes[i]; if (!have) todo += sizes[i]
      } }
    END { printf "\n  %.1f GB selected, %.1f GB still to fetch\n", total / 1e9, todo / 1e9 }'
  exit 0
fi

fetched=0; skipped=0; failed=0
echo "$rows" | while IFS="$(printf '\t')" read -r run alias platform ftp bytes; do
  [ -n "$ftp" ] || continue
  # One run carries one file for Nanopore and two for paired Illumina, semicolon-separated.
  echo "$ftp" | tr ';' '\n' | while read -r path; do
    [ -n "$path" ] || continue
    name=$(basename "$path")
    out="${fastqdir}/${name}"
    if [ -s "$out" ]; then
      echo "SKIP  ${name} exists  (${alias}, ${platform})"
      skipped=$((skipped + 1))
      continue
    fi
    echo "GET   ${name}  (${alias}, ${platform})"
    # --fail so an HTML error page is not written out as if it were data, and a partial file from
    # an interrupted transfer is removed rather than left to look complete on the next run.
    if ! curl -SfL -m 3600 --retry 3 --retry-delay 5 -o "${out}.part" "https://${path}"; then
      echo "  FAILED ${name}" >&2
      rm -f "${out}.part"
      failed=$((failed + 1))
      continue
    fi
    if ! gzip -t "${out}.part" 2>/dev/null; then
      echo "  FAILED ${name}: not valid gzip - the host probably sent something else." >&2
      rm -f "${out}.part"
      failed=$((failed + 1))
      continue
    fi
    mv "${out}.part" "$out"
    fetched=$((fetched + 1))
  done
done

echo
echo "=== data/fastq ==="
# Counted from the map rather than from the loop above: the while loop runs in a subshell of the
# pipe, so its counters do not survive it, and what matters is how much of the map is satisfied.
map="${fastqdir}/cdiff_isolates.txt"
if [ -f "$map" ]; then
  want=0; have=0
  while read -r key file; do
    case "$key" in ''|\#*) continue ;; esac
    want=$((want + 1))
    [ -s "${fastqdir}/${file}" ] && have=$((have + 1))
  done < "$map"
  echo "${have} of ${want} files of $(basename "$map") are in place."
  if [ "$have" -lt "$want" ]; then
    echo "Repeat this command to fetch the rest; what is already here is skipped." >&2
  fi
else
  echo "WARNING: ${map} is missing - the classification runs read the isolates through it." >&2
fi
