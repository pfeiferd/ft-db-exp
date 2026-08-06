#!/bin/sh
#
# Generates the simulated fastq files with known ground truth that the classification experiments
# are based on, following the approach of the first Genestrip paper:
#
#   viral        InSilicoSeq applies its Illumina "MiSeq" and "HiSeq" error models to exactly
#   protozoa     those genomes the database was built from, extracted by the goal
#   gut-protozoa `extractrefseqfasta' (see make_iss below).
#   tick-borne   NanoSim trains an error model on the real Nanopore reads of a tick sample and
#                applies it to the RefSeq genomes of the twelve tick-borne genera.
#
# Every part skips whatever is already present, so the script can be re-run safely. Each can be run
# on its own:
#
#   sh ./bin/make_fastqs.sh viral
#   sh ./bin/make_fastqs.sh protozoa
#   sh ./bin/make_fastqs.sh gut-protozoa
#   sh ./bin/make_fastqs.sh tick-borne
#   N_READS=10k sh ./bin/make_fastqs.sh viral      # quick smoke test instead of a full run
#   ERROR_FREE=1 sh ./bin/make_fastqs.sh protozoa  # error-free reads, see below
#   ERROR_NANOPORE=1 sh ./bin/make_fastqs.sh protozoa      # Nanopore-level per-base error
#   NANOPORE_ERROR_PCT=15 ERROR_NANOPORE=1 sh ./bin/make_fastqs.sh protozoa   # ... at 15 %
#   NANOPORE_READ_LENGTH=1000 ERROR_NANOPORE=1 sh ./bin/make_fastqs.sh protozoa  # ... 1 kb reads
#
# A single project follows ERROR_FREE; `all' ignores it and generates both read sets of every
# InSilicoSeq project, since the paper reports them side by side.
#
# ERROR_FREE uses InSilicoSeq's "perfect" mode, which fragments the genomes into reads of realistic
# length but introduces no sequencing errors at all. The resulting figures are an upper bound: they
# show what the refinement achieves when nothing but the taxonomy limits the classification, which
# separates the effect of the refinement from the effect of read errors.
#
# Raw tick data:
#   The real Nanopore runs NanoSim trains on come from the tick surveillance study the first paper
#   builds on (SRA study PRJNA790938). They are declared by URL in the Genestrip fastq map
#   data/fastq/ticks_real.txt, and Genestrip's own goal `fastqdownload' fetches them into data/fastq
#   as tick1.fastq.gz .. tick8.fastq.gz:
#
#     mvn exec:exec@fastqdl -Dname=tick-borne -Dfqmap=ticks_real.txt
#
#   The tick-borne part below runs that for whatever is missing, so there is nothing to do by hand.
#   Set SKIP_FETCH=1 to suppress it, e.g. when the reads were copied over from another machine.
#
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)

what=${1:-all}
fastqdir="${basedir}/data/fastq"
workdir="${basedir}/tools/work"
mkdir -p "$fastqdir"

cpus=${CPUS:-$(nproc 2>/dev/null || echo 4)}

############################## viral / InSilicoSeq ##############################

# $1 = database project name
#
# The reads are drawn from exactly the genomes the database was built from, which the goal
# `extractrefseqfasta' writes out one FASTA per region. That goal applies the same selection as the
# database fill -- the same completeness setting, the same caps per taxon -- so every read has a
# counterpart in the database.
#
# Drawing them from the whole RefSeq category instead, as this script did before, meant that most
# reads came from organisms the database never covered: for a project requesting ten genera out of
# the category `protozoa', about seven of eight reads were unusable, and they had to be excluded
# from the evaluation as out of scope. Extracting first keeps the full million reads in play.
make_iss() {
  db=$1

  iss="${basedir}/tools/iss-venv/bin/iss"
  if [ ! -x "$iss" ]; then
    echo "InSilicoSeq is missing - run ./bin/install_tools.sh first." >&2
    exit 1
  fi

  fastadir="${basedir}/data/projects/${db}/fasta"
  if [ -z "$(ls -A "$fastadir" 2>/dev/null)" ]; then
    echo "=== ${db}: extracting the genomes the database was built from ==="
    # Via `extractrefseqcsv', not `extractrefseqfasta'. The latter is an object goal, and Genestrip
    # treats object goals as weak dependencies: asking for one on the command line makes the
    # internal wrapper goal but never the object goal itself, so it silently does nothing. The CSV
    # goal is an ordinary file goal that reads the object goal's value, which makes it run and write
    # the per-accession FASTA files as a side effect.
    ( cd "$basedir" && mvn exec:exec@db -Dname="$db" -Dgoal=extractrefseqcsv )
  fi
  if [ -z "$(ls -A "$fastadir" 2>/dev/null)" ]; then
    echo "Goal extractrefseqcsv produced no files in ${fastadir}." >&2
    echo "Build the database first: mvn exec:exec@db -Dname=${db} -Dgoal=db" >&2
    exit 1
  fi

  # One million reads per model, as in the first paper.
  n_reads=${N_READS:-1M}
  # The first Genestrip paper puts the per-base error of Nanopore devices at 5 % to 15 %; ten is
  # the middle of that range, and its worked example of 6 % lies within it too.
  nanopore_error=${NANOPORE_ERROR_PCT:-10}
  # The same paper reports a mean read length of 3,926 bp for its NanoSim data, and this regime was
  # first run that way. One million reads of that length are 7.9 Gbp, some thirty times the other
  # regimes, which costs hours in generation and again in classification for a result that differs
  # from the short-read one mainly in how *few* reads are left for a refinement to improve. The
  # default is therefore the short read length, which isolates the error rate at a fraction of the
  # cost; set NANOPORE_READ_LENGTH=3926 to reproduce the long-read variant.
  nanopore_read_length=${NANOPORE_READ_LENGTH:-125}
  mkdir -p "$workdir"

  # InSilicoSeq wants one uncompressed multi-FASTA. The concatenation is removed again at the end,
  # since it is large and quick to recreate from the extracted files.
  genomes="${workdir}/${db}.genomic.fna"
  if [ ! -f "$genomes" ]; then
    echo "Collecting the extracted genomes of ${db} ..."
    # The headers keep the form extractrefseqfasta writes them in, ">ACCESSION|kraken:taxid|TAXID",
    # so the extracted files stay usable for Kraken 2 library building as well. InSilicoSeq names
    # every read after the whole header; IssReadGroundTruth cuts the accession out of it.
    : > "$genomes"
    for f in "$fastadir"/*; do
      case "$f" in
        *.gz) gunzip -c "$f" >> "$genomes" ;;
        *)    cat "$f" >> "$genomes" ;;
      esac
    done
    echo "  $(grep -c '^>' "$genomes") sequence(s), $(du -h "$genomes" | cut -f1)"
  fi

  if [ -n "${ERROR_FREE:-}" ]; then
    models="perfect"
    # Own map name, so that a run does not clobber the map of another error regime.
    mapfile="${fastqdir}/${db}_sim_perfect.txt"
  elif [ -n "${ERROR_NANOPORE:-}" ]; then
    models="nanopore"
    mapfile="${fastqdir}/${db}_sim_nanopore.txt"
  else
    models="miseq hiseq"
    mapfile="${fastqdir}/${db}_sim.txt"
  fi

  for model in $models; do
    prefix="${fastqdir}/${db}_iss_${model}_reads"
    if [ -f "${prefix}_R1.fastq.gz" ]; then
      echo "SKIP  ${prefix}_R1.fastq.gz exists"
      continue
    fi
    echo "=== ${db}: generating ${n_reads} ${model} reads ==="
    if [ "$model" = perfect ]; then
      # No error model at all: the reads differ from the reference only by where they were cut.
      "$iss" generate --genomes "$genomes" --mode perfect --n_reads "$n_reads" \
        --cpus "$cpus" --compress --output "$prefix"
    elif [ "$model" = nanopore ]; then
      # InSilicoSeq's `basic' model substitutes a base with probability 10^(-q/10), so a phred
      # value of q gives a uniform per-base error rate of choice; ISS_BASIC_PHRED is honoured by
      # the patch install_tools.sh applies. Note that these are substitutions only, whereas real
      # Nanopore error is indel-heavy -- what carries over is how many k-mers survive a read,
      # which is what the classification depends on.
      phred=$(awk -v e="$nanopore_error" 'BEGIN{printf "%d", -10*log(e/100)/log(10) + 0.5}')
      echo "  per-base error ${nanopore_error} % -> phred ${phred}, read length ${nanopore_read_length} bp"
      # Records shorter than the read length are skipped by InSilicoSeq, so a long read length
      # restricts the simulation to the longer contigs of the extracted genomes.
      ISS_BASIC_PHRED="$phred" ISS_BASIC_READ_LENGTH="$nanopore_read_length" \
        "$iss" generate --genomes "$genomes" --mode basic \
        --n_reads "$n_reads" --cpus "$cpus" --compress --output "$prefix"
    else
      "$iss" generate --genomes "$genomes" --model "$model" --n_reads "$n_reads" \
        --cpus "$cpus" --compress --output "$prefix"
    fi
    # None of InSilicoSeq's scratch output is part of the ground truth: the abundance table, the
    # per-thread VCFs and, when a run is interrupted, the partial fastq files it leaves behind.
    rm -f "${prefix}_abundance.txt" "${prefix}".iss.tmp.*
  done

  rm -f "$genomes"
  rmdir "$workdir" 2>/dev/null || true

  # Reads sharing a key are reported together, so each model forms one key over both mates.
  # The map goes next to the fastq files: Genestrip resolves a map file name against the literal
  # path, the project's fastq directory and data/fastq -- but not against the project's txt folder.
  # The entries themselves are bare file names, which GSProject.fastqFilesFromPath() resolves the
  # same way and finds in data/fastq. Absolute paths would work too but would tie the map to the
  # machine that wrote it, and the maps are regenerated with the reads anyway.
  {
    for model in $models; do
      for mate in 1 2; do
        echo "iss_${model} ${db}_iss_${model}_reads_R${mate}.fastq.gz"
      done
    done
  } > "$mapfile"
  echo "Wrote ${mapfile}"
}

# All three error regimes of a database, which the paper reports side by side:
#
#   error-free   bounds what the refinement can achieve when nothing but the taxonomy limits the
#                classification
#   MiSeq/HiSeq  the Illumina error models, at 301 bp and 126 bp per mate respectively
#   Nanopore     long reads at a uniform per-base error, both in the range the first Genestrip
#                paper reports for Nanopore devices, which tests whether the effect survives a far
#                higher error rate and a far greater read length at once
make_iss_all_regimes() {
  saved_error_free=${ERROR_FREE:-}
  saved_nanopore=${ERROR_NANOPORE:-}
  ERROR_FREE=""; ERROR_NANOPORE=""
  make_iss "$1"
  ERROR_FREE=1; ERROR_NANOPORE=""
  make_iss "$1"
  ERROR_FREE=""; ERROR_NANOPORE=1
  make_iss "$1"
  ERROR_FREE=$saved_error_free
  ERROR_NANOPORE=$saved_nanopore
}

############################## tick-borne / NanoSim ##############################

make_ticks() {
  nsvenv="${basedir}/tools/nanosim-venv"
  nanosimdir="${basedir}/tools/NanoSim"
  if [ ! -x "${nsvenv}/bin/python" ] || [ ! -f "${nanosimdir}/src/simulator.py" ]; then
    echo "NanoSim is missing - run ./bin/install_tools.sh first." >&2
    exit 1
  fi
  export PATH="${basedir}/tools/bin:${PATH}"
  python="${nsvenv}/bin/python"

  # The genome list maps each reference genome to a label that ends up as a prefix of every
  # simulated read's name. Following the first paper the label is "<taxid>x<index>", which keeps
  # the read's true taxon visible. Generated by org.metagene.ftdbexp.nanosim.NanoSimGenomeList,
  # which in turn needs the per-accession fasta files of the goal `extractrefseqfasta`.
  genomes=${GENOME_LIST:-${basedir}/data/projects/tick-borne/csv/tick-borne_nanosim.tsv}
  if [ ! -f "$genomes" ]; then
    echo "Missing genome list ${genomes}." >&2
    echo "Create it with: mvn exec:exec@nanosimlist -Dname=tick-borne" >&2
    exit 1
  fi

  nswork="${basedir}/tools/work-nanosim"
  mkdir -p "$nswork"

  # NanoSim undershoots the requested read count by roughly two orders of magnitude, so the first
  # paper asked for ten million to end up at a few hundred thousand. Kept for comparability.
  reads=${READS:-10000000}
  # Eight samples, as in the first paper, which generated one simulated fastq file per real one.
  samples=${SAMPLES:-"tick1 tick2 tick3 tick4 tick5 tick6 tick7 tick8"}

  # The real reads are the training input, so a missing one silently changes what is simulated
  # rather than causing a visible failure. Fetch whatever is absent before starting, and stop if
  # any is still missing afterwards -- a run of several hours is a bad place to discover it.
  #
  # The download is Genestrip's own goal `fastqdownload' rather than a script of ours: the map
  # data/fastq/ticks_real.txt declares the eight runs by URL, and Genestrip resolves, downloads and
  # names them, skipping whatever is already there. See that file for the accessions.
  if [ -z "${SKIP_FETCH:-}" ]; then
    missing=""
    for sample in $samples; do
      [ -s "${fastqdir}/${sample}.fastq.gz" ] || missing="${missing} ${sample}"
    done
    if [ -n "$missing" ]; then
      echo "=== downloading the real tick reads:${missing} ==="
      ( cd "$basedir" && mvn exec:exec@fastqdl -Dname=tick-borne -Dfqmap=ticks_real.txt )
    fi
  fi
  absent=""
  for sample in $samples; do
    [ -s "${fastqdir}/${sample}.fastq.gz" ] || absent="${absent} ${sample}"
  done
  if [ -n "$absent" ]; then
    echo "Missing real tick reads:${absent}" >&2
    echo "Run 'sh ./bin/fetch_tick_reads.sh' first, or set SAMPLES to the ones actually present." >&2
    exit 1
  fi

  for sample in $samples; do
    out="${fastqdir}/${sample}_sim.fastq"
    if [ -s "$out" ]; then
      echo "SKIP  ${out} exists"
      continue
    fi
    input="${fastqdir}/${sample}.fastq.gz"

    echo "=== ${sample}: training the error model on the real reads ==="
    ( cd "$nswork" && "$python" "${nanosimdir}/src/read_analysis.py" metagenome \
        -q --fastq -gl "$genomes" -i "$input" -t "$cpus" )

    echo "=== ${sample}: simulating ${reads} reads ==="
    sed -i "s/Abundance/${reads}/g" "${nswork}/training_quantification.tsv"
    ( cd "$nswork" && "$python" "${nanosimdir}/src/simulator.py" metagenome \
        --seed 42 --fastq -gl "$genomes" -t "$cpus" -a training_quantification.tsv )

    mv "${nswork}/simulated_sample0_aligned_reads.fastq" "$out"
    rm -f "${nswork}"/training* "${nswork}/reference_metagenome.fasta"
    echo "OK    ${out}"
  done

  # One key per sample, as in the first paper's ticks_sim.txt, and again next to the fastq files.
  {
    for sample in $samples; do
      [ -s "${fastqdir}/${sample}_sim.fastq" ] && echo "${sample} ${sample}_sim.fastq"
    done
  } > "${fastqdir}/ticks_sim.txt"
  echo "Wrote ${fastqdir}/ticks_sim.txt"
}

case "$what" in
  viral)         make_iss viral ;;
  protozoa)      make_iss protozoa ;;
  gut-protozoa)  make_iss gut-protozoa ;;
  tick-borne)    make_ticks ;;
  all)           make_iss_all_regimes viral; make_iss_all_regimes protozoa
                 make_iss_all_regimes gut-protozoa
                 make_ticks ;;
  *)             echo "Usage: $0 [viral|protozoa|gut-protozoa|tick-borne|all]" >&2; exit 1 ;;
esac

echo
echo "=== data/fastq ==="
ls -la "$fastqdir"
