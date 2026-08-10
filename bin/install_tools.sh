#!/bin/sh
#
# Installs every external tool the experiments need, all of them below ./tools so that nothing
# leaks into the system beyond a handful of distribution packages.
#
#   cgmemtime      measures wall time and peak RAM of the database generation
#   InSilicoSeq    simulates the Illumina reads for the viral experiments
#   NanoSim        simulates the Nanopore reads for the tick-borne experiments
#   mlst           assigns sequence types to the genomes of the `cdiff' database, which is what
#                  its sub-species refinement is interpreted and scored against
#
# Nothing is needed here for the real sequencing runs: fetch_saliva.sh and ticks_real.txt both pull
# gzipped fastq files straight over HTTPS with curl, which every machine already has. sra-toolkit
# used to be installed for that and is not any more -- see the header of fetch_saliva.sh for why
# going through prefetch/fasterq-dump cannot work at these volumes.
#
# The script is idempotent: anything already in place is skipped, so it is cheap to re-run after a
# partial failure.
#
# NanoSim deliberately deviates from the conda recipe of the first Genestrip paper. That recipe
# cannot be reproduced on every machine -- the conda hosts may be unreachable, and NanoSim's
# requirements.txt pins versions targeting Python 3.7 (numpy 1.21.5, scikit-learn 0.22.1,
# pysam 0.15.3) that do not build against a current interpreter. The external binaries therefore
# come from the distribution, the Python dependencies from PyPI without the pins. Should NanoSim
# ever fail on an API change in a newer numpy or scikit-learn, the honest fix is to run it in the
# original conda environment on another machine rather than to patch around it here.
#
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)

toolsdir="${basedir}/tools"
bindir="${toolsdir}/bin"
mkdir -p "$toolsdir" "$bindir"

echo "############ 1/5  Distribution packages ############"
# minimap2 and LAST align reads during NanoSim's training, samtools and genometools handle the
# sequences. python3-dev supplies the headers pybedtools compiles its C extension against, bedtools
# the binary it drives. All of this is what the conda recipe would have pulled from bioconda.
#
sudo apt-get install -y -q \
  minimap2 samtools last-align genometools bedtools \
  build-essential python3-dev python3-venv >/dev/null
echo "  minimap2:    $(minimap2 --version 2>&1 | head -1)"
echo "  samtools:    $(samtools --version 2>&1 | head -1)"
echo "  lastal:      $(lastal --version 2>&1 | head -1)"
echo "  genometools: $(gt --version 2>&1 | head -1)"

echo "############ 2/5  cgmemtime ############"
if [ -x "${toolsdir}/cgmemtime/cgmemtime" ]; then
  echo "  already built"
else
  tmp=$(mktemp -d)
  curl -sSL "https://codeload.github.com/gsauthof/cgmemtime/tar.gz/refs/heads/master" \
    | tar xz -C "$tmp" --strip-components=1
  # Let it print to stdout so the experiment scripts can capture it.
  sed -i 's/print_result(stderr/print_result(stdout/' "$tmp/cgmemtime.c"
  ( cd "$tmp" && make >/dev/null )
  mkdir -p "${toolsdir}/cgmemtime"
  cp "$tmp/cgmemtime" "${toolsdir}/cgmemtime/"
  rm -rf "$tmp"
  echo "  built ${toolsdir}/cgmemtime/cgmemtime"
fi

echo "############ 3/5  InSilicoSeq ############"
issvenv="${toolsdir}/iss-venv"
# A venv records absolute paths, so a moved or renamed one is broken even though its files look
# fine. Test that iss actually runs rather than that the file exists.
if "${issvenv}/bin/iss" --version >/dev/null 2>&1; then
  echo "  already installed: $("${issvenv}/bin/iss" --version 2>&1 | head -1)"
else
  python3 -m venv "$issvenv"
  "${issvenv}/bin/pip" install --quiet --upgrade pip
  "${issvenv}/bin/pip" install --quiet InSilicoSeq
  echo "  installed $("${issvenv}/bin/iss" --version 2>&1 | head -1)"
fi

# InSilicoSeq 2.0.1 cannot generate error-free reads: PerfectErrorModel does not set
# `store_mutations', unlike BasicErrorModel and KDErrorModel, so the base class raises an
# AttributeError while writing the reads. `iss' catches it, prints its own help and exits 0 -- the
# run leaves scratch files behind, no output, and a zero exit status that hides the failure. The
# patch adds the missing attribute; it is idempotent and a no-op once upstream fixes this.
perfectmodel=$(ls -d "${issvenv}"/lib/python*/site-packages/iss/error_models/perfect.py 2>/dev/null | head -1)
if [ -n "$perfectmodel" ] && ! grep -q "store_mutations" "$perfectmodel"; then
  echo "  patching InSilicoSeq: PerfectErrorModel is missing store_mutations"
  python3 - "$perfectmodel" <<'PATCH'
import sys
path = sys.argv[1]
src = open(path, encoding='utf-8').read()
src = src.replace(
    "    def __init__(self, fragment_length=None, fragment_sd=None):",
    "    def __init__(self, fragment_length=None, fragment_sd=None, store_mutations=False):", 1)
src = src.replace(
    "        self.fragment_sd = fragment_sd\n        self.quality_forward",
    "        self.fragment_sd = fragment_sd\n        self.store_mutations = store_mutations\n"
    "        self.quality_forward", 1)
open(path, 'w', encoding='utf-8').write(src)
PATCH
fi

# The `basic' error model hardcodes a mean phred score of 30 (a per-base error of 0.1 %) and a read
# length of 125. The experiments also need a far higher, Nanopore-like error rate and read length,
# so both are made configurable through ISS_BASIC_PHRED and ISS_BASIC_READ_LENGTH. mut_sequence()
# substitutes a base with probability 10^(-q/10), hence q = -10*log10(e) for a target error e.
# Without the variables the model behaves exactly as before.
basicmodel=$(ls -d "${issvenv}"/lib/python*/site-packages/iss/error_models/basic.py 2>/dev/null | head -1)
if [ -n "$basicmodel" ] && ! grep -q "ISS_BASIC_READ_LENGTH" "$basicmodel"; then
  echo "  patching InSilicoSeq: making the basic model's phred score and read length configurable"
  python3 - "$basicmodel" <<'PATCH'
import sys
path = sys.argv[1]
src = open(path, encoding='utf-8').read()
src = src.replace(
    "        self.quality_forward = self.quality_reverse = 30",
    '        self.quality_forward = self.quality_reverse = int(os.environ.get("ISS_BASIC_PHRED", "30"))', 1)
src = src.replace(
    "        self.read_length = 125",
    '        self.read_length = int(os.environ.get("ISS_BASIC_READ_LENGTH", "125"))', 1)
if "import os" not in src:
    src = src.replace("import numpy as np", "import os\n\nimport numpy as np", 1)
open(path, 'w', encoding='utf-8').write(src)
PATCH
fi

echo "############ 4/5  NanoSim ############"
nsvenv="${toolsdir}/nanosim-venv"
nanosimdir="${toolsdir}/NanoSim"

if [ -x "${bindir}/sam2pairwise" ]; then
  echo "  sam2pairwise already built"
else
  # Not packaged anywhere, so it is built from source. A single small C++ program.
  tmp=$(mktemp -d)
  curl -sSL "https://codeload.github.com/mlafave/sam2pairwise/tar.gz/refs/heads/master" \
    | tar xz -C "$tmp" --strip-components=1
  ( cd "$tmp/src" && make >/dev/null )
  cp "$tmp/src/sam2pairwise" "${bindir}/"
  rm -rf "$tmp"
  echo "  built ${bindir}/sam2pairwise"
fi

if [ -f "${nanosimdir}/src/simulator.py" ]; then
  echo "  NanoSim source already present"
else
  # ~660 MB, because the tarball ships the pre-trained error models.
  mkdir -p "$nanosimdir"
  curl -sSL "https://codeload.github.com/bcgsc/NanoSim/tar.gz/refs/heads/master" \
    | tar xz -C "$nanosimdir" --strip-components=1
  echo "  unpacked ${nanosimdir}"
fi

# The pre-trained error models are removed rather than kept, because they cannot be used with the
# unpinned dependencies installed below and failing to notice that is easy. They were pickled under
# scikit-learn 0.22, and loading them under a current one runs into four successive incompatibilities:
#
#   sklearn.neighbors.kde              renamed to sklearn.neighbors._kde        (0.22)
#   sklearn.neighbors._dist_metrics    moved to sklearn.metrics._dist_metrics   (1.3)
#   EuclideanDistance                  split into 32- and 64-bit variants       (1.3)
#   KernelDensity.bandwidth_           fitted-attribute contract changed
#
# The first three can be aliased away, and doing so is what makes this dangerous: the model then
# unpickles into a structurally invalid state, and `simulator.py' reports "Finished!" while writing
# zero-byte fastq files. A run that fails this way is indistinguishable from a successful one by
# exit status alone. Checked on 2026-08-05 against scikit-learn 1.9.0.
#
# Two ways to get Nanopore reads regardless, neither of which needs these files:
#   - train a model with `read_analysis.py', which *creates* the pickles under the scikit-learn that
#     is actually installed, so no version skew arises. This is what make_fastqs.sh does for the
#     tick-borne data, and it is the route the first Genestrip paper took.
#   - run the simulation on a machine with the original conda environment of that paper
#     (Python 3.7, scikit-learn 0.22.1), where the shipped models load as intended.
if [ -d "${nanosimdir}/pre-trained_models" ]; then
  echo "  removing NanoSim's pre-trained models - unusable with current scikit-learn, see the"
  echo "  comment in this script; train a model with read_analysis.py instead"
  rm -rf "${nanosimdir}/pre-trained_models"
fi

if [ ! -x "${nsvenv}/bin/python" ]; then
  python3 -m venv "$nsvenv"
fi
"${nsvenv}/bin/pip" install --quiet --upgrade pip
# Unpinned on purpose, see the header. `regex` is not in requirements.txt although NanoSim imports
# it (model_homopolymer_lengths.py); the conda recipe installed it separately for the same reason.
"${nsvenv}/bin/pip" install --quiet \
  numpy scipy scikit-learn pysam HTSeq joblib six pybedtools piecewise-regression regex

# NanoSim always simulates the "unaligned" reads too -- the share of a real run that maps to nothing,
# which it reproduces from the aligned/unaligned ratio it measured while training. They are of no use
# here: they are drawn from a flat error table and from no genome at all, so they carry no taxon in
# their name and no ground truth can be recovered from them, and make_fastqs.sh discards the file
# straight away. Generating them is not free either -- simulation_unaligned() walks base by base in
# Python -- and metagenome mode offers no option to skip it.
#
# The patch adds one: with NANOSIM_NO_UNALIGNED set, the block that simulates and merges them is
# skipped entirely, so neither the per-process subfiles nor the merged file are written. Without the
# variable NanoSim behaves exactly as before, and `perfect' mode already skipped the block anyway.
nosim_sim="${nanosimdir}/src/simulator.py"
if [ -f "$nosim_sim" ] && ! grep -q 'NANOSIM_NO_UNALIGNED' "$nosim_sim"; then
  echo "  patching NanoSim: making the simulation of unaligned reads skippable"
  python3 - "$nosim_sim" <<'PATCH'
import sys
path = sys.argv[1]
src = open(path, encoding='utf-8').read()
old = """    # Simulate unaligned reads, if per, number_unaligned = 0, taken care of in read_ecdf
    if not per:"""
new = """    # Simulate unaligned reads, if per, number_unaligned = 0, taken care of in read_ecdf
    # NANOSIM_NO_UNALIGNED skips them entirely; see bin/install_tools.sh of genestrip's ft-db-exp.
    if not per and not os.environ.get("NANOSIM_NO_UNALIGNED"):"""
if old not in src:
    sys.exit("simulator.py does not look as expected - not patching")
open(path, 'w', encoding='utf-8').write(src.replace(old, new, 1))
PATCH
fi

# Optional patch, off by default. With a reference larger than minimap2's default index batch size
# of 4G the index gets split, which slows NanoSim's training down considerably. The experiments of
# the first paper therefore raised it to 24G -- but minimap2 then holds the whole index in memory,
# so a value beyond the machine's RAM turns a slow run into a failing one. Set it deliberately:
#   MINIMAP2_INDEX_SIZE=24G sh ./bin/install_tools.sh
if [ -n "${MINIMAP2_INDEX_SIZE:-}" ]; then
  if grep -q 'map-ont -t " + num_threads + " -I ' "${nanosimdir}/src/read_analysis.py"; then
    echo "  minimap2 index size already patched"
  else
    sed -i "s/map-ont -t \" + num_threads + \" \"/map-ont -t \" + num_threads + \" -I ${MINIMAP2_INDEX_SIZE} \"/g" \
      "${nanosimdir}/src/read_analysis.py"
    echo "  patched read_analysis.py with -I ${MINIMAP2_INDEX_SIZE}"
  fi
fi

echo "############ 5/5  mlst ############"
# tseemann/mlst types an assembly against the PubMLST schemes it ships with. It is needed only by
# bin/mlst_assemblies.sh, which supplies the sequence types the `cdiff' database's dendrogram is
# interpreted and scored against -- see that script for why. It is a Perl program driving blast+
# and any2fasta, with no build step of its own.
#
# A failure here is a warning rather than an error. Everything else this script installs is needed
# for the read simulations, which are the bulk of the experiments; MLST typing concerns one database
# and can be caught up later.
#
# THE PERL PROBLEM. mlst declares `use 5.32.0', and its older releases still want 5.26, so on a
# distribution whose Perl predates that it does not even compile. Replacing the system Perl is no
# option -- a Debian or Ubuntu carries a great deal of its own tooling on it -- so a Perl of its own
# is built below tools/ with perlbrew when, and only when, the one on the PATH is too old. That
# build takes some twenty to forty minutes; it happens once and is skipped from then on.
mlst_perl_min_major=5
mlst_perl_min_minor=32
mlstperl=""

perl_version_ok() {
  # $1 = a perl executable. True if it is at least the version mlst asks for.
  [ -x "$1" ] || command -v "$1" >/dev/null 2>&1 || return 1
  "$1" -e "require ${mlst_perl_min_major}.${mlst_perl_min_minor}.0;" >/dev/null 2>&1
}

if perl_version_ok perl; then
  mlstperl=$(command -v perl)
  echo "  OK  system perl is recent enough: $(perl -e 'print $^V')"
else
  echo "  system perl is $(perl -e 'print $^V' 2>/dev/null || echo unknown), which is older than" \
       "v${mlst_perl_min_major}.${mlst_perl_min_minor} - building one below tools/"
  export PERLBREW_ROOT="${toolsdir}/perlbrew"
  perlbrew_perl=perl-5.36.3
  perlbrew_bin="${PERLBREW_ROOT}/perls/${perlbrew_perl}/bin/perl"
  if [ -x "$perlbrew_bin" ]; then
    echo "  SKIP  ${perlbrew_bin} exists"
    mlstperl="$perlbrew_bin"
  else
    if [ ! -x "${PERLBREW_ROOT}/bin/perlbrew" ]; then
      # perlbrew's own installer is a Perl script and runs on the old Perl; only the Perl it builds
      # has to be new.
      if ! curl -sSL https://install.perlbrew.pl | bash >/dev/null 2>&1; then
        echo "  WARNING: could not install perlbrew into ${PERLBREW_ROOT}." >&2
      fi
    fi
    if [ -x "${PERLBREW_ROOT}/bin/perlbrew" ]; then
      echo "  building ${perlbrew_perl} - this takes 20-40 minutes and happens once"
      # -n skips the test suite, which is the bulk of the time and tests Perl rather than anything
      # this project depends on. -j uses the cores the machine has.
      if "${PERLBREW_ROOT}/bin/perlbrew" install -n -j "$(nproc 2>/dev/null || echo 4)" \
           "$perlbrew_perl" >"${toolsdir}/perlbrew-install.log" 2>&1; then
        echo "  built ${perlbrew_perl}"
        mlstperl="$perlbrew_bin"
      else
        echo "  WARNING: building ${perlbrew_perl} failed - see ${toolsdir}/perlbrew-install.log" >&2
      fi
    fi
  fi
fi

# The modules mlst and the MLST::* packages beside it load. They go into whichever Perl was settled
# on above: a freshly built one starts out with core modules only, and the distribution packages an
# older run may have installed belong to the system Perl, not to it.
if [ -n "$mlstperl" ]; then
  if ! "$mlstperl" -MMoo -MList::MoreUtils -MJSON -MPath::Tiny -e1 >/dev/null 2>&1; then
    echo "  installing the Perl modules mlst needs"
    cpanm_bin="$(dirname "$mlstperl")/cpanm"
    if [ ! -x "$cpanm_bin" ]; then
      curl -sSL https://cpanmin.us -o "$cpanm_bin" 2>/dev/null && chmod +x "$cpanm_bin" || true
    fi
    if [ -x "$cpanm_bin" ]; then
      "$mlstperl" "$cpanm_bin" --quiet --notest Moo List::MoreUtils JSON Path::Tiny \
        >"${toolsdir}/cpanm-install.log" 2>&1 || true
    fi
    if ! "$mlstperl" -MMoo -MList::MoreUtils -MJSON -MPath::Tiny -e1 >/dev/null 2>&1; then
      echo "  WARNING: some of Moo, List::MoreUtils, JSON, Path::Tiny are still missing -" >&2
      echo "           see ${toolsdir}/cpanm-install.log" >&2
    fi
  else
    echo "  OK  the Perl modules mlst needs are present"
  fi
fi

# blast+ and any2fasta are the two executables mlst requires at run time (`require_exe' in its
# bin/mlst). any2fasta is a self-contained Perl script that runs on any Perl, so it is cloned rather
# than packaged, and blast+ comes from the distribution.
sudo apt-get install -y -q ncbi-blast+ >/dev/null 2>&1 || true
if ! command -v blastn >/dev/null 2>&1 && [ ! -x "${bindir}/blastn" ]; then
  echo "  WARNING: blastn is not on the PATH - mlst cannot type anything without it." >&2
  echo "           Install ncbi-blast+ (needs root) or put a blastn into ${bindir}." >&2
fi
a2fdir="${toolsdir}/any2fasta"
if [ -x "${a2fdir}/any2fasta" ]; then
  echo "  SKIP  ${a2fdir} exists"
elif git clone -q --depth 1 https://github.com/tseemann/any2fasta.git "$a2fdir" >/dev/null 2>&1; then
  echo "  cloned ${a2fdir}"
else
  echo "  WARNING: could not install any2fasta - mlst requires it beside blastn." >&2
fi
[ -x "${a2fdir}/any2fasta" ] && ln -sf "${a2fdir}/any2fasta" "${bindir}/any2fasta"

mlstdir="${toolsdir}/mlst"
if [ -d "${mlstdir}/bin" ]; then
  echo "  SKIP  ${mlstdir} exists"
elif git clone -q --depth 1 https://github.com/tseemann/mlst.git "$mlstdir"; then
  echo "  cloned ${mlstdir}"
else
  echo "  WARNING: could not clone mlst - bin/mlst_assemblies.sh will not run." >&2
  mlstdir=""
fi

if [ -n "$mlstdir" ] && [ -f "${mlstdir}/bin/mlst" ] && [ -n "$mlstperl" ]; then
  # A wrapper rather than a symlink: mlst has to run under the Perl settled on above, which is not
  # the one a `#!/usr/bin/env perl' would find, and it has to find any2fasta and blastn no matter
  # what the caller's PATH holds. FindBin still resolves to the real script, so mlst locates its own
  # db/pubmlst as it expects.
  cat > "${bindir}/mlst" <<WRAPPER
#!/bin/sh
# Generated by bin/install_tools.sh - edit that instead.
PATH="${bindir}:\${PATH}"
export PATH
exec "${mlstperl}" "${mlstdir}/bin/mlst" "\$@"
WRAPPER
  chmod +x "${bindir}/mlst"

  # The scheme name matters: bin/mlst_assemblies.sh defaults to `cdifficile', which is what
  # PubMLST's database pubmlst_cdifficile_seqdef becomes here. Checking it now beats discovering
  # after a typing run that every genome came back as `-'.
  if "${bindir}/mlst" --list 2>/dev/null | tr ' ' '\n' | grep -qx cdifficile; then
    echo "  OK  mlst runs and knows the 'cdifficile' scheme"
  else
    echo "  WARNING: mlst does not run, or has no 'cdifficile' scheme. What it says:" >&2
    "${bindir}/mlst" --list 2>&1 >/dev/null | sed 's/^/           /' >&2 || true
  fi
fi

echo
echo "############ Smoke test ############"
export PATH="${bindir}:${PATH}"
for script in simulator read_analysis; do
  version=$("${nsvenv}/bin/python" "${nanosimdir}/src/${script}.py" --version 2>&1 | grep -o 'NanoSim [0-9.]*' | head -1)
  if [ -n "$version" ]; then
    echo "  OK  ${script}.py: ${version}"
  else
    echo "  FAILED  ${script}.py did not run:" >&2
    "${nsvenv}/bin/python" "${nanosimdir}/src/${script}.py" --version 2>&1 | tail -15 >&2
    exit 1
  fi
done
echo "  OK  iss: $("${issvenv}/bin/iss" --version 2>&1 | head -1)"
echo
echo "All tools installed below ${toolsdir}. Next: sh ./bin/make_fastqs.sh"
