#!/bin/sh
#
# Read classification performance of the unrefined against the refined database, in the three
# scenarios Table "analysisperf" of the first Genestrip paper uses -- but comparing the two variants
# of Genestrip-FT rather than four different tools:
#
#   (a) `cv' on the human saliva runs, each file on its own
#   (b) `tb' on the real tick runs, each file on its own
#   (c) `tb' on the real tick runs, all of them analyzed jointly in one Genestrip call
#   (d) `strepto' on the same saliva runs as (a), each file on its own
#
# (b) and (c) differ in more than bookkeeping: a joint run loads the database once and amortises it
# over every file, which is why the first paper reports both. (a) and (d) differ in the database
# alone -- same files, same machine, a viral database against a bacterial one of a single genus --
# which is what makes the cost of classifying comparable across databases rather than across tools.
#
# Usage:
#   sh ./bin/perf_scenarios.sh              # run all four, then write the CSV
#   sh ./bin/perf_scenarios.sh a            # one scenario
#   sh ./bin/perf_scenarios.sh csv          # only re-derive the CSV from existing logs
#   sh ./bin/perf_scenarios.sh store        # only the k-mer store comparison of the appendix
#
# Writes results/perf_<scenario>_<key>_<goal>.log per measured run and results/matchperf.csv, which
# Table "matchperf" of the paper reads.
#
# DB_SUFFIX runs the same four scenarios against differently named databases and keeps their results
# apart: with DB_SUFFIX=-sa the scenarios classify against `viral-sa', `tick-borne-sa' and
# `strepto-sa', and write perf_<scenario>sa_... logs and matchperf-sa.csv. That is how the k-mer store
# comparison of bin/store_compare.sh measures the binary-search store without disturbing the
# measurement of the default one. The scenario letters stay the same, so the two CSVs line up row by
# row.
#
# GS_XMX caps the JVM heap of every measured run, e.g. GS_XMX=6G. It matters for the RAM column:
# cgmemtime reports the high-water resident set, which for a 56G limit is the heap G1 decided to
# commit and not what the run needs. The ftmatch log of scenario (a) shows 6208 MB committed against
# 3491 MB used, so the column says more about the garbage collector than about the database. A cap
# just above what a run needs makes the number mean the database again. It also costs: a tighter
# heap collects more often, so wall time and throughput of such a run are not comparable with the
# figures measured at the default.
#
# THE MEASUREMENT IS THE POINT, so this script refuses to guess. A run whose log is missing or
# unparseable leaves its cell empty rather than being averaged over the runs that did work.
#
set -e

# Passed on to every measured run; empty keeps the pom's gs.xmx.
xmx_opt=""
if [ -n "${GS_XMX:-}" ]; then
  xmx_opt="-Dgs.xmx=${GS_XMX}"
fi

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
basedir=$(pwd)

fastqdir="${basedir}/data/fastq"
res_path="${RESULTS_DIR:-${basedir}/results}"
mkdir -p "$res_path"

what=${1:-all}

# Appended to every project name below, empty by default. The log and CSV names take it too, with the
# dash dropped, so a suffixed run writes perf_asa_... beside perf_a_... rather than over it.
db_suffix=${DB_SUFFIX:-}
key_suffix=$(echo "$db_suffix" | tr -d '-')

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
# The output files of one goal for every key of a map, removed so that the goal is actually made.
#
# `-Dgs.target=clean' below is supposed to do this and evidently did not on 2026-09-19 -- every tick
# run reported its goals as made and classified nothing. Removing the files here as well costs a
# millisecond and does not depend on the clean target behaving as expected.
clear_goal_outputs() {
  _cg_db=$1; _cg_goal=$2; _cg_map="${fastqdir}/$3"
  [ -f "$_cg_map" ] || return 0
  for _cg_key in $(map_keys "$_cg_map"); do
    rm -f "${basedir}/data/projects/${_cg_db}/csv/${_cg_db}_${_cg_goal}_${_cg_key}.csv" \
          "${basedir}/data/projects/${_cg_db}/csv/${_cg_db}_${_cg_goal}_${_cg_key}.csv.gz"
  done
}

# Whether every fastq file a map names is on disk. A goal whose input is missing is not made either,
# and that looks exactly like a goal that was already made: no work, no error, a one-second log.
map_files_present() {
  _mf_map="${fastqdir}/$1"
  _mf_absent=""
  for _mf_f in $(awk '!/^#/ && NF >= 2 { print $2 }' "$_mf_map"); do
    [ -s "${fastqdir}/${_mf_f}" ] || _mf_absent="${_mf_absent} ${_mf_f}"
  done
  if [ -n "$_mf_absent" ]; then
    echo "  SKIP $1 -- fastq files not on disk:${_mf_absent}" >&2
    return 1
  fi
  return 0
}

run_pair() {
  _rp_db=$1; _rp_map=$2; _rp_key=$3
  map_files_present "$_rp_map" || return 0
  for goal in match ftmatch; do
    echo "############ ${_rp_key}: ${goal} ############"
    # Clear the previous output first and do NOT measure that: it is bookkeeping, not classification.
    # Both ways, because (b) and (c) run over the same map and so write the same files: whatever (b)
    # leaves behind would make (c) a no-op, which is what the batch of 2026-09-19 recorded.
    clear_goal_outputs "$_rp_db" "$goal" "$_rp_map"
    mvn exec:exec@match -Dname="$_rp_db" -Dgoal="$goal" -Dfqmap="$_rp_map" -Dgs.target=clean $xmx_opt
    measure "${res_path}/perf_${_rp_key}_${goal}.log" \
        mvn exec:exec@match -Dname="$_rp_db" -Dgoal="$goal" -Dfqmap="$_rp_map" $xmx_opt
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

# Everything a measured run must not find already there.
#
# Genestrip does not remake a file goal whose output is on disk: it reports the goal as made, exits
# after a second, and cgmemtime then times a JVM start -- 0.9 s and 180 MB, which reads like a
# measurement and is none. That is exactly what the tick runs of 2026-09-19 recorded. So the two
# goals' results go first, for every key of every map this script measures, together with this
# script's own logs, so that no figure of an earlier batch survives into the next CSV.
#
# Deliberately narrow: only `match' and `ftmatch' of the projects and keys below. The databases, the
# fastq files and every other goal's results are left alone -- they are days of work and are not
# what stands in the way here. KEEP_RESULTS=1 skips this, for a re-run that only wants the logs of
# runs it is about to do anyway.
clear_previous_results() {
  if [ -n "${KEEP_RESULTS:-}" ]; then
    echo "KEEP_RESULTS is set - leaving previous match results in place."
    echo "  A goal whose result is still there will not be remade, and its run measures nothing." >&2
    return 0
  fi
  for _cp_sc in $1; do
    rm -f "${res_path}"/perf_${_cp_sc}${key_suffix}_*.log
  done
  for _cp_spec in "a:viral${db_suffix}:${SALIVA_MAP:-saliva_real.txt}" "b:tick-borne${db_suffix}:eightticks.txt" \
                  "c:tick-borne${db_suffix}:eightticks.txt" "d:strepto${db_suffix}:${SALIVA_MAP:-saliva_real.txt}"; do
    case " $1 " in *" ${_cp_spec%%:*} "*) ;; *) continue ;; esac
    _cp_spec=${_cp_spec#*:}
    _cp_db=${_cp_spec%%:*}
    _cp_map="${fastqdir}/${_cp_spec#*:}"
    [ -f "$_cp_map" ] || continue
    for _cp_key in $(map_keys "$_cp_map"); do
      for _cp_f in "${basedir}/data/projects/${_cp_db}/csv/${_cp_db}_match_${_cp_key}.csv" \
                   "${basedir}/data/projects/${_cp_db}/csv/${_cp_db}_ftmatch_${_cp_key}.csv"; do
        for _cp_x in "$_cp_f" "${_cp_f}.gz"; do
          if [ -f "$_cp_x" ]; then
            echo "  removing $(basename "$_cp_x")"
            rm -f "$_cp_x"
          fi
        done
      done
    done
  done
}

# Only for the scenarios about to run: a measurement of the others is days of wall time and must
# survive a re-run of one of them. `all' clears everything, since it remeasures everything.
case "$what" in
  all)     clear_previous_results "a b c d" ;;
  a|b|c|d) clear_previous_results "$what" ;;
esac

case "$what" in
  a|all) run_per_file "viral${db_suffix}" "${SALIVA_MAP:-saliva_real.txt}" "a${key_suffix}" ;;
esac
case "$what" in
  b|all) run_per_file "tick-borne${db_suffix}" eightticks.txt "b${key_suffix}" ;;
esac
case "$what" in
  c|all) run_pair "tick-borne${db_suffix}" eightticks.txt "c${key_suffix}_joint" ;;
esac
case "$what" in
  d|all) run_per_file "strepto${db_suffix}" "${SALIVA_MAP:-saliva_real.txt}" "d${key_suffix}" ;;
esac
case "$what" in
  a|b|c|d|store|all|csv) ;;
  *) echo "Usage: $0 [a|b|c|d|store|all|csv]" >&2; exit 1 ;;
esac

# The k-mer store comparison of the paper's appendix, which is a measurement of the same kind and
# belongs to the same table-filling run. It takes the real read sets, `strepto' the smallest saliva run
# and `tick-borne' all eight tick runs, and classifies each of them twice per store, once with the
# Bloom filter of the match and once without it. With the filter on, most k-mers of a real read set
# never reach the store and the layout hardly shows; with it off, every k-mer does. The two together
# also say what the filter itself is worth. `viral' is left out: scenario (a) already shows what its
# lookups are worth, and `strepto' covers the same reads against a database that classifies far more
# of them.
#
# Skipped for a DB_SUFFIX run: the twins exist to measure one store at a time, while this compares
# both in one JVM and needs no twin at all.
case "$what" in
  store|all)
    if [ -n "$db_suffix" ]; then
      echo "############ store comparison skipped for the ${db_suffix} twins ############"
    else
      REAL=1 FILTER=both sh ./bin/store_bench.sh strepto
      REAL=1 ALL=1 FILTER=both sh ./bin/store_bench.sh tick-borne
    fi
    ;;
esac

# --- deriving the table --------------------------------------------------------------------------
#
# Wall time and RAM come from cgmemtime's own summary, which it appends to the log:
#
#     wall:   0.763 s
#     child_RSS_high:     181716 KiB
#     group_mem_high:     201220 KiB
#
# child_RSS_high is the high-water mark of the largest single process -- the JVM -- and that is what
# this reports. group_mem_high sits beside it and is the whole cgroup, which sounds like the honest
# figure and is not: the kernel charges the page cache of everything the group reads to it, and
# scenario (a) reads 78 to 171 GB of gzip per run. Measured on the saliva runs of 2026-09-19, where
# the database itself takes 3,485 MB of heap: child_RSS_high averages 6,756 MB unrefined and 5,977 MB
# refined, group_mem_high 38,619 and 47,861 -- a factor of 6 to 8 of cache, and the refined figures
# sit at 49 GB for all three runs, i.e. against the 56 GB heap ceiling rather than at any demand.
# Table "ftdbperf" of the paper reports child_RSS_high too, so the two tables now agree.
#
# The read counts are taken from the evaluation's own CSVs rather than recounted here: they are
# exact, already written, and counting a 139 GB fastq again to reproduce them would cost an hour.

python3 - "$res_path" "$fastqdir" "$bp_sample_reads" "$key_suffix" "$db_suffix" <<'PY'
import csv
import glob, os, re, subprocess, sys

res, fastqdir, bp_sample = sys.argv[1], sys.argv[2], int(sys.argv[3])
# The scenario letters stay a, b, c, d in the table; the logs of a suffixed run carry the suffix.
key_suffix, db_suffix = sys.argv[4], sys.argv[5]

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
    # A goal whose results are already on disk is not remade: Genestrip reports the file goals as
    # made in no time and exits, and cgmemtime then times a JVM start -- 0.9 s and 180 MB, which
    # looks like a measurement and is not one. The run that classified reads says so in its log
    # ("Making match took 1342 s"); the one that did nothing never mentions the goal at all.
    if not re.search(r'Making (ft)?match took', text):
        print('  not measured: %s classified nothing (its goal was already made)'
              % os.path.basename(path), file=sys.stderr)
        return None
    w = re.search(r'^wall:\s+([0-9.]+)\s*s', text, re.M)
    m = re.search(r'^child_RSS_high:\s+([0-9]+)\s*KiB', text, re.M)
    if not m:
        m = re.search(r'^group_mem_high:\s+([0-9]+)\s*KiB', text, re.M)
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
    missing = [k for k, x in zip(keys, n) if not x]
    n = [x for x in n if x]
    # The read counts are not reported here any more -- the paper's table of the real data carries
    # them -- but they still divide the wall time into a speed, so a run whose classification was
    # never evaluated has none and the speed is left out rather than averaged over the rest. Say which.
    if missing:
        print('  %s: no read count for %s - its speed is left out'
              % (label, ', '.join(missing)), file=sys.stderr)
        n = []
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
        g = measurement(os.path.join(res, 'perf_c%s_joint_%s.log' % (key_suffix, goal)))
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
per_file('a' + key_suffix, '(a)')
per_file('b' + key_suffix, '(b)')
joint('(c)')
per_file('d' + key_suffix, '(d)')

out = os.path.join(res, 'matchperf%s.csv' % db_suffix)
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
