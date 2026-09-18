#!/bin/sh
#
# Downloads the human saliva runs of the first Genestrip paper and writes them as gzipped fastq
# files into data/fastq, where the map data/fastq/saliva_real.txt expects them. They are the real
# data without ground truth that `cv' is applied to.
#
# Usage:
#   sh ./bin/fetch_saliva.sh                      # the three runs the first paper used
#   sh ./bin/fetch_saliva.sh SRR5571991           # a single run
#   SKIP_MD5=1 sh ./bin/fetch_saliva.sh           # do not verify (saves ~10 min per 70 GB file)
#
# Files are named after their accession, exactly as bin/make_fastqs.sh of the original
# genestrip-db-exp project named them, so that runs already downloaded there are picked up rather
# than fetched again -- which at these volumes is the difference between minutes and days.
#
# BEWARE OF THE VOLUME, but note that it is the *download* that is large, not the working set.
# All five runs of that table, in its own order, with the sizes the ENA reports today:
#
#   SRR5571985  #1   812,085,208 read pairs   164 Gbp   112 GB gzipped
#   ERR1395613  #2   900,709,176 read pairs   180 Gbp   171 GB gzipped
#   SRR5571991  #3   980,879,835 read pairs   198 Gbp   139 GB gzipped   <- reported in detail there
#   ERR1395610  #4   824,479,570 read pairs   165 Gbp   154 GB gzipped
#   SRR5571990  #5   605,561,636 read pairs   122 Gbp    78 GB gzipped
#
#   all five                                            654 GB gzipped
#
# The two ERR runs are the largest of the set, so fetching them roughly doubles what the other three
# cost. All five are Illumina HiSeq 2000; the SRR runs are 101 bp per mate, the ERR runs 100 bp.
#
# ---------------------------------------------------------------------------------------------
# Why this fetches from the ENA over HTTPS rather than through sra-tools
#
# The obvious route -- `prefetch' followed by `fasterq-dump' -- fails on these runs with
#
#     disk-limit exeeded!
#     fasterq-dump quit with error code 3
#
# and it fails for a good reason, even on a filesystem with hundreds of gigabytes free. The SRA
# stores reads in a compressed columnar format, so producing fastq from it means writing the data
# out *uncompressed* first. For SRR5571991 the three stages coexist on disk:
#
#     .sra archive         ~55 GB      written by prefetch
#     uncompressed fastq  ~420 GB      written by fasterq-dump   <- 198 Gbp of bases plus qualities,
#     gzipped fastq        139 GB      what we actually want        headers and newlines
#
# i.e. roughly 615 GB transiently to arrive at 139 GB. fasterq-dump estimates that up front and
# refuses rather than filling the disk halfway through a multi-hour extraction. Raising the ceiling
# with --disk-limit/--disk-limit-tmp does not make the space appear.
#
# The ENA mirrors every SRA run and serves it already gzipped, split into mates, so the download is
# the final file: 139 GB fetched, 139 GB stored, no scratch space and no conversion step. It is also
# resumable (curl -C -) and checksummed (the ENA reports an md5 per file), which were the two
# properties sra-tools was chosen for in the first place. This is the same approach ticks_real.txt
# takes for the Nanopore runs, so both real datasets now arrive the same way.
#
# That route is therefore not offered at all, and install_tools.sh no longer installs sra-toolkit.
# ---------------------------------------------------------------------------------------------
#
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)

fastqdir="${basedir}/data/fastq"
mkdir -p "$fastqdir"


# The runs of Table "errorsviral" of the first Genestrip paper -- all five of them, in that table's
# order. The original project's make_fastqs.sh fetched only the three SRR runs and left ERR1395613
# (#2) and ERR1395610 (#4) commented out, which is why they were missing here too; but that paper
# reports measurements for all five (see its results/viral_human_virus_errors_gs_ku_comp.csv and
# results/viral_saliva_rel_accuracy.csv), so leaving them out would report on less data than the
# study this one builds upon. The ENA resolves ERR accessions natively -- fetch_via_ena() asks the
# filereport API and derives nothing from the prefix -- so they present no special difficulty.
#
# The default is the three SRR runs. The two ERR ones are the largest of the five -- 171 GB and
# 154 GB, more than the other three together -- and were dropped again on 2026-08-22 for disk and
# wall time. They are still fetched on request, since the ENA resolves them like any other
# accession: sh ./bin/fetch_saliva.sh ERR1395613 ERR1395610 .
DEFAULT_RUNS="SRR5571985 SRR5571991 SRR5571990"

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
  while [ "$mate" -le 2 ]; do
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
    # -C - resumes a partial .part rather than starting the 70 GB over; --retry rides out the
    # transient failures a transfer of this length inevitably meets. The download goes to .part and
    # is moved into place only once its md5 checks out, so an interrupted run never leaves a file
    # that the presence check above would accept.
    # A transfer of this length rarely fails outright; it stalls. curl's --retry fires on an error,
    # not on a connection that goes quiet without closing, so such a transfer hangs and the run
    # appears to "break" without saying anything. --speed-limit/--speed-time turn a stall into an
    # error that --retry can act on, and the loop resumes across what retries do not cover: an
    # exhausted retry budget, a peer that drops the connection at the same byte count every time, a
    # proxy that times the request out. Every attempt resumes the same .part, so the cost of another
    # attempt is the bytes still missing, never the file.
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
