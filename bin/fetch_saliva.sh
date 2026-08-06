#!/bin/sh
#
# Downloads the human saliva runs of the first Genestrip paper and writes them as gzipped fastq
# files into data/fastq, where the map data/fastq/saliva_real.txt expects them. They are the real
# data without ground truth that `cv' is applied to.
#
# Usage:
#   sh ./bin/fetch_saliva.sh                      # the three runs the first paper used
#   sh ./bin/fetch_saliva.sh SRR5571991           # a single run
#   sh ./bin/fetch_saliva.sh SRR5571991 ERR1395613
#   KEEP_SRA=1 sh ./bin/fetch_saliva.sh           # keep the .sra files instead of deleting them
#
# Files are named after their accession, exactly as bin/make_fastqs.sh of the original
# genestrip-db-exp project named them, so that runs already downloaded there are picked up rather
# than fetched again -- which at these volumes is the difference between minutes and days.
#
# BEWARE OF THE VOLUME. These are deep metagenomic runs of 605 to 981 million read pairs, 122 to 198
# Gbp each. Expect roughly 400 GB of gzipped fastq for the three, and transiently about the same again
# for the .sra files and fasterq-dump's scratch space. There is no point starting this without a
# quarter of a terabyte free per run. Run it one at a time unless the machine has the room: the
# paper's argument needs one such file, not three, and the first paper reported SRR5571991 in detail.
#
# Why sra-tools rather than a plain URL, as ticks_real.txt uses: these runs are paired, and the two
# mates have to arrive as two files that Genestrip can read as one key. `prefetch' also resumes an
# interrupted transfer, which matters at 100+ GB per run, and verifies the download against NCBI's
# own checksum -- neither of which a bare HTTP GET does.
#
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)

fastqdir="${basedir}/data/fastq"
sracache="${basedir}/tools/sra-cache"
mkdir -p "$fastqdir" "$sracache"

cpus=${CPUS:-$(nproc 2>/dev/null || echo 4)}

for tool in prefetch fasterq-dump; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "${tool} is missing - run ./bin/install_tools.sh first." >&2
    exit 1
  fi
done

# pigz saves hours at these volumes; gzip does the same job serially.
if command -v pigz >/dev/null 2>&1; then
  compress="pigz -p ${cpus}"
else
  compress="gzip"
fi

# The runs of Table "errorsviral" of the first Genestrip paper. That paper lists five, but the
# original project's make_fastqs.sh fetched only these three through sra-tools and left the other two
# -- ERR1395613 (#2) and ERR1395610 (#4) -- commented out. They are ENA-native accessions, which
# prefetch does not resolve as readily; pass them explicitly if they are wanted.
#
#   SRR5571985  #1   812,085,208 read pairs   164 Gbp
#   SRR5571991  #3   980,879,835 read pairs   198 Gbp   <- the run that paper reports in detail
#   SRR5571990  #5   605,561,636 read pairs   122 Gbp
DEFAULT_RUNS="SRR5571991 SRR5571990 SRR5571985"

runs=${*:-$DEFAULT_RUNS}

for acc in $runs; do

  # Both mates present means done. Checking only the first would restart a sample whose second file
  # is still being written, and silently leave the pair inconsistent.
  if [ -s "${fastqdir}/${acc}_1.fastq.gz" ] && [ -s "${fastqdir}/${acc}_2.fastq.gz" ]; then
    echo "SKIP  ${acc} - both mates present"
    continue
  fi

  echo "=== ${acc}: prefetch ==="
  # The default refusal above 20 GB has to be lifted; these runs are far larger. 200g is what the
  # original project passed.
  prefetch --max-size 200g --output-directory "$sracache" "$acc"

  echo "=== ${acc}: fasterq-dump ==="
  # --split-3 is fasterq-dump's default and writes <acc>_1.fastq and <acc>_2.fastq for a paired run.
  # The scratch directory is put next to the output rather than in /tmp, which is rarely large enough
  # for a run of this size.
  fasterq-dump --threads "$cpus" \
      --temp "$sracache" --outdir "$sracache" "${sracache}/${acc}"

  for mate in 1 2; do
    raw="${sracache}/${acc}_${mate}.fastq"
    if [ ! -s "$raw" ]; then
      echo "fasterq-dump produced no ${raw} - is ${acc} really paired?" >&2
      exit 1
    fi
    echo "  compressing mate ${mate} ..."
    # Via .part and mv, so that an interrupted run leaves nothing that looks finished. Without it a
    # half-written .fastq.gz would satisfy the presence check above on the next run and the run would
    # be silently skipped -- after hours of transfer, and with a truncated file in place.
    $compress -c "$raw" > "${fastqdir}/${acc}_${mate}.fastq.gz.part"
    mv "${fastqdir}/${acc}_${mate}.fastq.gz.part" "${fastqdir}/${acc}_${mate}.fastq.gz"
    rm -f "$raw"
  done

  if [ -z "${KEEP_SRA:-}" ]; then
    rm -rf "${sracache}/${acc}"
  fi
  echo "OK    ${fastqdir}/${acc}_[12].fastq.gz"
done

echo
echo "=== data/fastq ==="
ls -la "$fastqdir"/SRR55719*.fastq.gz "$fastqdir"/ERR139*.fastq.gz 2>/dev/null || echo "nothing downloaded yet"
