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
#   The real Nanopore runs come from the tick surveillance study the first paper builds on and are
#   available from the SRA as tick1=SRR17281117, tick2=SRR17281105, tick3=SRR17281103,
#   tick4=SRR17281101, tick5=SRR17281100, tick6=SRR17281099, tick7=SRR17281116, tick8=SRR17281115:
#
#     wget "https://www.be-md.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=<acc>" -O tickN.fastq.gz
#
#   Place them in data/fastq. That NCBI host is not reachable from every network.
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
    # Own map name, so that an error-free run does not clobber the map of the regular one.
    mapfile="${fastqdir}/${db}_sim_perfect.txt"
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
  {
    for model in $models; do
      for mate in 1 2; do
        echo "iss_${model} ${fastqdir}/${db}_iss_${model}_reads_R${mate}.fastq.gz"
      done
    done
  } > "$mapfile"
  echo "Wrote ${mapfile}"
}

# Both read sets of a database: with an Illumina error model and error-free. A comprehensive run
# needs both, since the paper reports them side by side -- the error-free figures bound what the
# refinement can achieve when nothing but the taxonomy limits the classification.
make_iss_both() {
  saved_error_free=${ERROR_FREE:-}
  ERROR_FREE=""
  make_iss "$1"
  ERROR_FREE=1
  make_iss "$1"
  ERROR_FREE=$saved_error_free
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
  samples=${SAMPLES:-"tick1 tick2 tick3 tick4 tick5 tick6 tick7"}

  for sample in $samples; do
    out="${fastqdir}/${sample}_sim.fastq"
    if [ -s "$out" ]; then
      echo "SKIP  ${out} exists"
      continue
    fi
    input="${fastqdir}/${sample}.fastq.gz"
    if [ ! -f "$input" ]; then
      echo "WARN  ${input} missing - skipping ${sample}" >&2
      continue
    fi

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
      [ -s "${fastqdir}/${sample}_sim.fastq" ] && echo "${sample} ${fastqdir}/${sample}_sim.fastq"
    done
  } > "${fastqdir}/ticks_sim.txt"
  echo "Wrote ${fastqdir}/ticks_sim.txt"
}

case "$what" in
  viral)         make_iss viral ;;
  protozoa)      make_iss protozoa ;;
  gut-protozoa)  make_iss gut-protozoa ;;
  tick-borne)    make_ticks ;;
  all)           make_iss_both viral; make_iss_both protozoa; make_iss_both gut-protozoa
                 make_ticks ;;
  *)             echo "Usage: $0 [viral|protozoa|gut-protozoa|tick-borne|all]" >&2; exit 1 ;;
esac

echo
echo "=== data/fastq ==="
ls -la "$fastqdir"
