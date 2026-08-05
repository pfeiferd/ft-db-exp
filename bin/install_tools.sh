#!/bin/sh
#
# Installs every external tool the experiments need, all of them below ./tools so that nothing
# leaks into the system beyond a handful of distribution packages.
#
#   cgmemtime      measures wall time and peak RAM of the database generation
#   InSilicoSeq    simulates the Illumina reads for the viral experiments
#   NanoSim        simulates the Nanopore reads for the tick-borne experiments
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

echo "############ 1/4  Distribution packages ############"
# minimap2 and LAST align reads during NanoSim's training, samtools and genometools handle the
# sequences. python3-dev supplies the headers pybedtools compiles its C extension against, bedtools
# the binary it drives. All of this is what the conda recipe would have pulled from bioconda.
sudo apt-get install -y -q \
  minimap2 samtools last-align genometools bedtools \
  build-essential python3-dev python3-venv >/dev/null
echo "  minimap2:    $(minimap2 --version 2>&1 | head -1)"
echo "  samtools:    $(samtools --version 2>&1 | head -1)"
echo "  lastal:      $(lastal --version 2>&1 | head -1)"
echo "  genometools: $(gt --version 2>&1 | head -1)"

echo "############ 2/4  cgmemtime ############"
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

echo "############ 3/4  InSilicoSeq ############"
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

echo "############ 4/4  NanoSim ############"
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

if [ ! -x "${nsvenv}/bin/python" ]; then
  python3 -m venv "$nsvenv"
fi
"${nsvenv}/bin/pip" install --quiet --upgrade pip
# Unpinned on purpose, see the header. `regex` is not in requirements.txt although NanoSim imports
# it (model_homopolymer_lengths.py); the conda recipe installed it separately for the same reason.
"${nsvenv}/bin/pip" install --quiet \
  numpy scipy scikit-learn pysam HTSeq joblib six pybedtools piecewise-regression regex

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
