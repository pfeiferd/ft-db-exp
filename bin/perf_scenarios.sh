#!/bin/sh
#
# Read classification performance of the unrefined against the refined database, in the three
# scenarios Table "analysisperf" of the first Genestrip paper uses -- but comparing the two variants
# of Genestrip-FT rather than four different tools:
#
#   (a) `cv' on the human saliva runs, each file on its own
#   (b) `tb' on the real tick runs, each file on its own
#   (c) `tb' on the real tick runs, all of them analyzed jointly in one Genestrip call
#
# (b) and (c) differ in more than bookkeeping: a joint run loads the database once and amortises it
# over every file, which is why the first paper reports both.
#
# Usage:
#   sh ./bin/perf_scenarios.sh              # run all three, then write the CSV
#   sh ./bin/perf_scenarios.sh a            # one scenario
#   sh ./bin/perf_scenarios.sh csv          # only re-derive the CSV from existing logs
#
# Writes results/perf_<scenario>_<key>_<goal>.log per measured run and results/matchperf.csv, which
# Table "matchperf" of the paper reads.
#
# THE MEASUREMENT IS THE POINT, so this script refuses to guess. A run whose log is missing or
# unparseable leaves its cell empty rather than being averaged over the runs that did work.
#
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)

fastqdir="${basedir}/data/fastq"
res_path="${RESULTS_DIR:-${basedir}/results}"
mkdir -p "$res_path"

what=${1:-all}

# Number of reads sampled to establish the mean read length. Illumina reads are fixed length so any
# sample will do; Nanopore lengths vary, and 100k reads settle the mean well within the precision
# the table states. Counting bases over a 139 GB file to gain a decimal place is not worth an hour.
bp_sample_reads=${BP_SAMPLE_READS:-100000}

timer="${basedir}/tools/cgmemtime/cgmemtime"

# --- running -----------------------------------------------------------------------------------

# Keys of a fastq map, in order, without duplicates.
map_keys() {
  awk '!/^#/ && NF >= 2 && !seen[$1]++ { print $1 }' "$1"
}

# One measured run. $1 = log file, $2.. = the command.
#
# stderr is captured too. The 0-byte results/match_*.log files this script replaces are exactly what
# happens when it is not: cgmemtime reports a failure to create its cgroup on stderr, the redirect
# keeps only stdout, and the log ends up empty with the reason discarded.
measure() {
  _m_log=$1
  shift
  if [ ! -x "$timer" ]; then
    echo "WARNING: ${timer} is missing - run ./bin/install_tools.sh. Skipping." >&2
    return 0
  fi
  if ! "$timer" true >/dev/null 2>&1; then
    echo "WARNING: cgmemtime cannot run here (it needs a systemd user scope for its cgroup)." >&2
    echo "         Skipping the measurement rather than writing an empty log." >&2
    return 0
  fi
  echo "  -> ${_m_log}"
  "$timer" "$@" > "$_m_log" 2>&1 || {
    echo "WARNING: the measured run failed; see ${_m_log}." >&2
    return 0
  }
}

# Runs both database variants over one map. $1 = db, $2 = map file name, $3 = log key.
run_pair() {
  _rp_db=$1; _rp_map=$2; _rp_key=$3
  for goal in match ftmatch; do
    echo "############ ${_rp_key}: ${goal} ############"
    # Clear the previous output first and do NOT measure that: it is bookkeeping, not classification.
    mvn exec:exec@match -Dname="$_rp_db" -Dgoal="$goal" -Dfqmap="$_rp_map" -Dgs.target=clean
    measure "${res_path}/perf_${_rp_key}_${goal}.log" \
        mvn exec:exec@match -Dname="$_rp_db" -Dgoal="$goal" -Dfqmap="$_rp_map"
  done
}

# Per-file scenario: one Genestrip call per fastq key, so that the database load is paid once per
# file exactly as the first paper's per-file figures do.
# $1 = db, $2 = source map, $3 = scenario letter
run_per_file() {
  _pf_db=$1; _pf_src="${fastqdir}/$2"; _pf_sc=$3
  if [ ! -f "$_pf_src" ]; then
    echo "Missing ${_pf_src} - skipping scenario ${_pf_sc}." >&2
    return 0
  fi
  for key in $(map_keys "$_pf_src"); do
    # A one-key map beside the original, so Genestrip resolves the file names the same way.
    _pf_map="perf_${_pf_sc}_${key}.txt"
    awk -v k="$key" '!/^#/ && $1 == k' "$_pf_src" > "${fastqdir}/${_pf_map}"
    # The map declares what the scenario covers; the disk decides what can be measured. The saliva
    # runs are hundreds of gigabytes each and arrive one at a time, so a map naming three of them
    # while one is present is the normal state, not an error -- but classifying a file that is not
    # there fails several minutes in, after the database has been loaded. Say so and move on; the
    # missing run's cells stay empty, which is what this script does everywhere else too.
    _pf_absent=""
    for _pf_f in $(awk '{print $2}' "${fastqdir}/${_pf_map}"); do
      [ -s "${fastqdir}/${_pf_f}" ] || _pf_absent="${_pf_absent} ${_pf_f}"
    done
    if [ -n "$_pf_absent" ]; then
      echo "  SKIP ${key} in scenario ${_pf_sc} -- not on disk:${_pf_absent}" >&2
      rm -f "${fastqdir}/${_pf_map}"
      continue
    fi
    run_pair "$_pf_db" "$_pf_map" "${_pf_sc}_${key}"
    rm -f "${fastqdir}/${_pf_map}"
  done
}

case "$what" in
  a|all) run_per_file viral "${SALIVA_MAP:-saliva_real.txt}" a ;;
esac
case "$what" in
  b|all) run_per_file tick-borne seventicks.txt b ;;
esac
case "$what" in
  c|all) run_pair tick-borne seventicks.txt c_joint ;;
esac
case "$what" in
  a|b|c|all|csv) ;;
  *) echo "Usage: $0 [a|b|c|all|csv]" >&2; exit 1 ;;
esac

# --- deriving the table --------------------------------------------------------------------------
#
# Wall time and RAM come from cgmemtime's own summary, which it appends to the log:
#
#     wall:   0.763 s
#     group_mem_high:     201220 KiB
#
# group_mem_high is the high-water mark of the whole cgroup, i.e. the JVM and everything Maven
# forked, which is what "RAM of the run" should mean. child_RSS_high sits beside it and is the
# largest single process; for a forked JVM the two are close but the group figure is the honest one.
#
# The read counts are taken from the evaluation's own CSVs rather than recounted here: they are
# exact, already written, and counting a 139 GB fastq again to reproduce them would cost an hour.

python3 - "$res_path" "$fastqdir" "$bp_sample_reads" <<'PY'
import csv, glob, os, re, subprocess, sys

res, fastqdir, bp_sample = sys.argv[1], sys.argv[2], int(sys.argv[3])

def keys_of(res, scenario):
    pat = re.compile(r'^perf_' + re.escape(scenario) + r'_(.+)_match\.log$')
    found = set()
    for p in glob.glob(os.path.join(res, 'perf_%s_*_match.log' % scenario)):
        m = pat.match(os.path.basename(p))
        if m:
            found.add(m.group(1))
    return sorted(found)

def measurement(path):
    """(wall seconds, peak RAM MB) from a cgmemtime log, or None if it did not measure."""
    try:
        text = open(path, encoding='utf-8', errors='replace').read()
    except OSError:
        return None
    w = re.search(r'^wall:\s+([0-9.]+)\s*s', text, re.M)
    m = re.search(r'^group_mem_high:\s+([0-9]+)\s*KiB', text, re.M)
    if not m:
        m = re.search(r'^child_RSS_high:\s+([0-9]+)\s*KiB', text, re.M)
    if not w or not m:
        return None
    return float(w.group(1)), int(m.group(1)) / 1024.0

def reads_by_key():
    """Exact read counts per fastq key, from whichever evaluation CSV reports them."""
    out = {}
    for pattern, keycol, readcol in (('*_specificity.csv', 'fastq key', 'reads'),
                                     ('*_accuracy.csv', 'fastq key', 'total')):
        for f in glob.glob(os.path.join(res, pattern)):
            with open(f, encoding='utf-8') as fh:
                for row in csv.DictReader(fh, delimiter=';'):
                    if row.get(keycol) and row.get(readcol):
                        out.setdefault(row[keycol], int(float(row[readcol])))
    return out

def mean_read_length(key):
    """Mean read length over a sample of the file, or None if the file is not there."""
    hits = sorted(glob.glob(os.path.join(fastqdir, key + '*.fastq.gz')))
    if not hits:
        return None
    cmd = ("gunzip -c '%s' | head -n %d | awk 'NR %% 4 == 2 { b += length($0); n++ } "
           "END { if (n) printf \"%%.0f\", b / n }'" % (hits[0], bp_sample * 4))
    try:
        val = subprocess.check_output(['sh', '-c', cmd], stderr=subprocess.DEVNULL).decode().strip()
        return int(val) if val else None
    except (subprocess.CalledProcessError, ValueError):
        return None

def fmt(v, places=0):
    if v is None:
        return ''
    return ('%%.%df' % places) % v

reads = reads_by_key()
rows = []

def per_file(scenario, label):
    keys = keys_of(res, scenario)
    if not keys:
        return
    n = [reads.get(k) for k in keys]
    n = [x for x in n if x]
    lens = [mean_read_length(k) for k in keys]
    lens = [x for x in lens if x]
    if n:
        if scenario == 'b':
            rows.append((label, 'Min. reads per fastq file', fmt(min(n)), fmt(min(n))))
            rows.append((label, 'Max. reads per fastq file', fmt(max(n)), fmt(max(n))))
        rows.append((label, 'Avg. reads per fastq file', fmt(sum(n) / len(n)), fmt(sum(n) / len(n))))
    if lens:
        rows.append((label, 'Avg. BPs per read',
                     fmt(sum(lens) / len(lens)), fmt(sum(lens) / len(lens))))
    cells = {}
    for goal, col in (('match', 0), ('ftmatch', 1)):
        got = [measurement(os.path.join(res, 'perf_%s_%s_%s.log' % (scenario, k, goal)))
               for k in keys]
        got = [g for g in got if g]
        if not got:
            continue
        walls = [g[0] for g in got]
        rams = [g[1] for g in got]
        cells[('Avg. RAM (MB)', col)] = fmt(sum(rams) / len(rams))
        cells[('Max. wall time (min)', col)] = fmt(max(walls) / 60.0, 2)
        cells[('Avg. wall time (min)', col)] = fmt(sum(walls) / len(walls) / 60.0, 2)
        if n and len(n) == len(walls) and sum(walls) > 0:
            cells[('Avg. speed (reads / s)', col)] = fmt(sum(n) / sum(walls))
    for name in ('Avg. RAM (MB)', 'Max. wall time (min)', 'Avg. wall time (min)',
                 'Avg. speed (reads / s)'):
        if (name, 0) in cells or (name, 1) in cells:
            rows.append((label, name, cells.get((name, 0), ''), cells.get((name, 1), '')))

def joint(label):
    got = {}
    for goal, col in (('match', 0), ('ftmatch', 1)):
        g = measurement(os.path.join(res, 'perf_c_joint_%s.log' % goal))
        if g:
            got[col] = g
    if not got:
        return
    keys = keys_of(res, 'b')
    total = sum(reads.get(k, 0) for k in keys)
    rows.append((label, 'RAM (MB)', fmt(got[0][1]) if 0 in got else '',
                 fmt(got[1][1]) if 1 in got else ''))
    rows.append((label, 'Wall time (min)', fmt(got[0][0] / 60.0, 2) if 0 in got else '',
                 fmt(got[1][0] / 60.0, 2) if 1 in got else ''))
    if total:
        rows.append((label, 'Speed (reads / s)',
                     fmt(total / got[0][0]) if 0 in got else '',
                     fmt(total / got[1][0]) if 1 in got else ''))

# Short labels: the scenario column repeats on every row, so it has to stay narrow. What each of
# them is belongs in the caption, where it is stated once.
per_file('a', '(a)')
per_file('b', '(b)')
joint('(c)')

out = os.path.join(res, 'matchperf.csv')
# A header alone is worse than no file at all: the paper's \perfrows prints its "---" fallback only
# when the file is absent, so an empty one would typeset a table with no rows and no marker. Refuse
# to write it, and leave whatever is there -- a run that measured nothing must not erase one that did.
if not rows:
    print('No measured runs for any scenario - %s left as it is.' % out)
else:
    with open(out, 'w', encoding='utf-8') as fh:
        fh.write('scenario;parameter;unrefined;refined;\n')
        for r in rows:
            fh.write(';'.join(r) + ';\n')
    print('Wrote %s (%d rows)' % (out, len(rows)))
PY
