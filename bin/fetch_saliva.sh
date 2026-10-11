#!/bin/sh
#
# Downloads the human saliva runs of the first Genestrip paper and writes them as gzipped fastq
# files into data/fastq, where the map data/fastq/saliva_real.txt expects them. They are the real
# data without ground truth that `cv' is applied to.
#
# Usage:
#   sh ./bin/fetch_saliva.sh                      # all five runs, first mate of each (326 GB)
#   sh ./bin/fetch_saliva.sh SRR5571991           # a single run
#   SKIP_MD5=1 sh ./bin/fetch_saliva.sh           # do not verify (saves ~10 min per 70 GB file)
#   MATES=2 sh ./bin/fetch_saliva.sh ERR1395613   # both mates of a run, which nothing reads today
#
# Files are named after their accession, exactly as bin/make_fastqs.sh of the original
# genestrip-db-exp project named them, so that runs already downloaded there are picked up rather
# than fetched again -- which at these volumes is the difference between minutes and days.
#
# BEWARE OF THE VOLUME -- though it is the *download* that is large, not the working set. All five
# runs in the first paper's order, with the sizes the ENA reports today. Only the first mate of each
# is fetched, which is what the analysis reads and what halves the total:
#
#   SRR5571985  #1   812,085,208 read pairs   164 Gbp    55.5 GB (mate 1)   112 GB (both)
#   ERR1395613  #2   900,709,176 read pairs   180 Gbp    85.7 GB            171 GB
#   SRR5571991  #3   980,879,835 read pairs   198 Gbp    69.1 GB            139 GB   <- reported in
#   ERR1395610  #4   824,479,570 read pairs   165 Gbp    77.1 GB            154 GB      detail there
#   SRR5571990  #5   605,561,636 read pairs   122 Gbp    38.5 GB             78 GB
#
#   all five                                            326 GB             654 GB
#
# All five are Illumina HiSeq 2000; the SRR runs are 101 bp per mate, the ERR runs 100 bp.
#
# Fetched from the ENA over HTTPS and not through sra-tools: `prefetch' + `fasterq-dump' write the
# fastq out uncompressed first, i.e. some 615 GB transiently to arrive at 139 GB for SRR5571991, and
# fasterq-dump refuses up front rather than filling the disk ("disk-limit exeeded!"). The ENA mirrors
# every run already gzipped and split into mates, so the download IS the final file, resumable
# (curl -C -) and checksummed against the md5 it reports. ticks_real.txt takes the same route.
#
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)

fastqdir="${basedir}/data/fastq"
mkdir -p "$fastqdir"


# All five runs of Table "errorsviral" of the first paper, in its order. That project fetched only
# the three SRR runs although it reports measurements for all five, so leaving the two ERR ones out
# would report on less data than the study this one builds upon. One mate each keeps it affordable
# -- 326 GB against 654 -- and is what the first paper analysed as well. A single accession is still
# fetched on its own: sh ./bin/fetch_saliva.sh ERR1395613 ERR1395610 .
DEFAULT_RUNS="SRR5571985 ERR1395613 SRR5571991 ERR1395610 SRR5571990"

# How many mates of each run to fetch. One, because that is what is analysed: Genestrip classifies
# mates independently, and saliva_real.txt lists the first of each run so that the read counts stay
# comparable with the first paper's. Fetching both would double a download of hundreds of gigabytes
# for files nothing opens. MATES=2 restores the old behaviour.
MATES=${MATES:-1}

runs=${*:-$DEFAULT_RUNS}

# --- md5, wherever it lives ------------------------------------------------------------------
md5_of() {
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "$1" | cut -d' ' -f1
  else
    md5 -q "$1"
  fi
}

# Verifying a 70 GB file costs minutes, so the result is remembered next to it. A run that is
# already verified is then skipped without re-reading it, which matters when this script is
# re-invoked after an interruption.
verify() {
  # $1 = file, $2 = expected md5
  if [ -z "$2" ]; then
    echo "  no md5 published for $(basename "$1") - cannot verify" >&2
    return 0
  fi
  if [ -n "${SKIP_MD5:-}" ]; then
    echo "  skipping md5 check of $(basename "$1") (SKIP_MD5 set)"
    return 0
  fi
  echo "  verifying $(basename "$1") ..."
  actual=$(md5_of "$1")
  if [ "$actual" != "$2" ]; then
    echo "md5 mismatch for $1: expected $2, got $actual" >&2
    return 1
  fi
  echo "$2" > "$1.md5"
  return 0
}

# Already present and known good?
have_it() {
  # $1 = file, $2 = expected md5
  [ -s "$1" ] || return 1
  [ -n "${SKIP_MD5:-}" ] && return 0
  [ -z "$2" ] && return 0
  [ -f "$1.md5" ] && [ "$(cat "$1.md5")" = "$2" ]
}

# --- the ENA route ---------------------------------------------------------------------------
fetch_via_ena() {
  acc=$1

  # The ENA's filereport API resolves an accession to its file URLs, md5s and sizes. Deriving the
  # path by hand is possible -- vol1/fastq/<first 6>/<0-padded last digits>/<acc>/ -- but the
  # padding rule depends on the accession's length and has changed before, so ask rather than guess.
  meta=$(curl -sS --max-time 120 \
      "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv" \
      | tail -n +2 | head -1)
  if [ -z "$meta" ]; then
    echo "the ENA knows no files for ${acc} - is the accession right?" >&2
    return 1
  fi

  urls=$(echo "$meta" | cut -f2)
  md5s=$(echo "$meta" | cut -f3)
  if [ -z "$urls" ]; then
    echo "the ENA lists no fastq files for ${acc}" >&2
    return 1
  fi

  mate=1
  while [ "$mate" -le "$MATES" ]; do
    url=$(echo "$urls" | cut -d';' -f"$mate")
    md5=$(echo "$md5s" | cut -d';' -f"$mate")
    if [ -z "$url" ]; then
      echo "${acc} has no mate ${mate} at the ENA - is it really paired?" >&2
      return 1
    fi
    target="${fastqdir}/${acc}_${mate}.fastq.gz"

    if have_it "$target" "$md5"; then
      echo "SKIP  $(basename "$target") - present and verified"
      mate=$((mate + 1))
      continue
    fi

    echo "=== ${acc} mate ${mate} ==="
    # -C - resumes the .part rather than starting the 70 GB over, and the file is moved into place
    # only once its md5 checks out, so an interrupted run never leaves one the presence check above
    # would accept. A transfer of this length rarely fails outright, it stalls -- and curl's --retry
    # fires on an error, not on a connection that goes quiet. --speed-limit/--speed-time turn a stall
    # into such an error, and the loop resumes across what retries do not cover.
    tries=0
    until curl -L --fail --retry 10 --retry-delay 15 --retry-connrefused ${retry_all_errors} \
          --speed-limit 1024 --speed-time 120 \
          -C - -o "${target}.part" "https://${url}"; do
      tries=$((tries + 1))
      if [ "$tries" -ge "${FETCH_TRIES:-20}" ]; then
        echo "gave up on ${acc} mate ${mate} after ${tries} attempts - ${target}.part is kept," >&2
        echo "so re-running resumes from where it stopped." >&2
        return 1
      fi
      echo "  attempt ${tries} ended early at $(du -h "${target}.part" 2>/dev/null | cut -f1) - resuming in 30 s" >&2
      sleep 30
    done

    if verify "${target}.part" "$md5"; then
      mv "${target}.part" "$target"
      [ -f "${target}.part.md5" ] && mv "${target}.part.md5" "${target}.md5"
      echo "OK    $target"
    else
      echo "leaving ${target}.part in place; re-run to resume or delete it to start over" >&2
      return 1
    fi
    mate=$((mate + 1))
  done
}

# --- go --------------------------------------------------------------------------------------
command -v curl >/dev/null 2>&1 || { echo "curl is missing." >&2; exit 1; }

# Retrying a transfer that ended in an HTTP error rather than a connection failure needs curl 7.71;
# older ones reject the option outright, which would fail every download rather than none.
retry_all_errors=""
if curl --help all 2>/dev/null | grep -q -- "--retry-all-errors"; then
  retry_all_errors="--retry-all-errors"
fi

for acc in $runs; do
  fetch_via_ena "$acc"
done

echo
echo "=== data/fastq ==="
# List what was asked for, not a fixed glob: passing an accession outside the default three used to
# print "nothing downloaded yet" over a directory that had just been filled.
for acc in $runs; do
  ls -la "$fastqdir"/"${acc}"_[12].fastq.gz 2>/dev/null || echo "  ${acc}: nothing present"
done
