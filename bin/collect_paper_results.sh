#!/bin/sh
#
# Collects the files of results/ that the paper reads and writes them to a zip beside that folder.
# results/ accumulates everything the experiments ever produced -- some 500 files -- and barely a
# fifth of it reaches the manuscript. This picks out that fifth, so the results can be handed over,
# archived or diffed without the rest.
#
# THE LIST BELOW IS FIXED ON PURPOSE. This project does not know where a paper lives and must not
# take its path as input, so the names are kept here rather than read out of a .tex file. The list
# was taken from the manuscript on 2026-10-08 and needs an edit whenever a table of the paper gains
# or drops a source file. Nothing goes wrong quietly if that is forgotten: the script names every
# entry it cannot find in results/, and an entry that is no longer read does no harm beyond a few
# kilobytes in the archive.
#
# The second group is not read by the paper at all. bin/paper_stats.sh needs it to produce
# results/dbstats.tex, and the paper does read that. The manuscript compiles from a shipped
# dbstats.tex alone, but it cannot be regenerated without these files, and a macro that silently
# keeps an old value is the one failure this project has already had. Pass --paper-only to leave
# them out.
#
# Usage:
#   sh ./bin/collect_paper_results.sh
#   sh ./bin/collect_paper_results.sh --paper-only
#   OUT=handover.zip sh ./bin/collect_paper_results.sh
#
set -e
scriptdir=$(dirname "$0")
cd "$scriptdir/.."

# ---------------------------------------------------------------------------------------------
# What the paper reads. Grouped as the manuscript uses them.
# ---------------------------------------------------------------------------------------------
# Included with \input: the generated macro files and the one dendrogram drawn in LaTeX.
PAPER_TEX="
dbstats.tex
sysinfo.tex
orthopox_dendrolatex_10242.tex
"
# The tax tree figures, included through \dbfig as SVG.
PAPER_FIGURES="
orthopox_svgtaxtree.svg
orthopox_ftsvgtaxtree.svg
orthopox2_ftsvgtaxtree.svg
"
# Read by \csvreader, one per table block. Database generation and classification performance,
# then the intrinsic quality figures, then the simulated and the real classification results.
PAPER_CSV="
db_gen_perf.csv
db_gen_perf-k24.csv
db_gen_perf-k24-s4.csv
matchperf.csv
storebench_strepto.csv
storebench_tick-borne.csv
protozoa_branchhistorankcsv.csv
protozoa_kmerrankstatscsv.csv
tick-borne_branchhistorankcsv.csv
tick-borne_kmerrankstatscsv.csv
orthopox_intersectcsv_10242.csv
viral_iss_simdata.csv
viral_iss_quality.csv
viral_iss_summary.csv
viral_iss_perfect_simdata.csv
viral_iss_perfect_quality.csv
viral_iss_perfect_summary.csv
viral_iss_saliva_simdata.csv
viral_iss_saliva_quality.csv
viral_iss_saliva_summary.csv
viral_ku_iss_quality.csv
viral_ku_iss_summary.csv
viral_ku_iss_perfect_quality.csv
viral_ku_iss_perfect_summary.csv
viral_ku_iss_saliva_quality.csv
viral_ku_iss_saliva_summary.csv
viral_k2_iss_quality.csv
viral_k2_iss_summary.csv
viral_k2_iss_perfect_quality.csv
viral_k2_iss_perfect_summary.csv
viral_k2_iss_saliva_quality.csv
viral_k2_iss_saliva_summary.csv
viral-k24_iss_quality.csv
viral-k24_iss_perfect_quality.csv
viral-k24_iss_saliva_quality.csv
viral-k24-s4_iss_quality.csv
viral-k24-s4_iss_perfect_quality.csv
viral-k24-s4_iss_saliva_quality.csv
strepto_iss_saliva_simdata.csv
strepto_iss_saliva_quality.csv
strepto_iss_saliva_summary.csv
nocardia_iss_mngs_simdata.csv
nocardia_iss_mngs_quality.csv
nocardia_iss_mngs_summary.csv
tick-borne_nanosim_simdata.csv
tick-borne_nanosim_quality.csv
tick-borne_nanosim_summary.csv
viral_saliva_specificity.csv
viral_ku_saliva_specificity.csv
viral_k2_saliva_specificity.csv
strepto_saliva_specificity.csv
nocardia_mngs_specificity.csv
tick-borne_ticks_specificity.csv
"

paperonly=
case "${1:-}" in
  --paper-only) paperonly=1; shift ;;
  "") ;;
  *) echo "Usage: sh ./bin/collect_paper_results.sh [--paper-only]" >&2; exit 1 ;;
esac
[ -d results ] || { echo "No results/ directory here. Run the experiments first." >&2; exit 1; }

out="${OUT:-paper_results.zip}"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT HUP INT TERM

: > "$tmp/want"; : > "$tmp/missing"
for f in $PAPER_TEX $PAPER_FIGURES $PAPER_CSV; do
  if [ -f "results/$f" ]; then echo "$f" >> "$tmp/want"; else echo "$f" >> "$tmp/missing"; fi
done
npaper=$(sort -u "$tmp/want" | wc -l | tr -d ' ')

# The inputs of bin/paper_stats.sh. Its DBS list is read out of that script rather than repeated
# here, so the two cannot drift apart. A database it skips for want of a dbinfo file is skipped
# here as well, without a complaint -- the list covers twins that may never have been built.
if [ -z "$paperonly" ]; then
  dbs=$(sed -n '/^DBS = \[/,/\]/p' bin/paper_stats.sh \
        | tr -d "[]'" | sed 's/^DBS = //' | tr ',' '\n' | tr -d ' ' | grep -v '^$')
  for db in $dbs; do
    for f in "${db}_dbinfo.csv" "${db}_ftdbinfo.csv" "${db}_dbquality.csv" "${db}_ftquality.csv"; do
      [ -f "results/$f" ] && echo "$f" >> "$tmp/want"
    done
    for f in results/"${db}"_*_specificity.csv; do
      [ -f "$f" ] && basename "$f" >> "$tmp/want"
    done
  done
fi

sort -u "$tmp/want" > "$tmp/list"
n=$(wc -l < "$tmp/list" | tr -d ' ')
[ "$n" -gt 0 ] || { echo "Nothing to collect -- results/ holds none of the listed files." >&2; exit 1; }

# The results/ prefix is kept, so that unpacking beside a paper puts every file where its
# \csvreader already looks for it.
rm -f "$out"
sed 's|^|results/|' "$tmp/list" > "$tmp/paths"
if command -v zip >/dev/null 2>&1; then
  zip -q -X "$out" -@ < "$tmp/paths"
else
  # No zip on this machine; python's zipfile writes the same archive and is always there.
  python3 - "$out" "$tmp/paths" <<'PY'
import sys, zipfile
out, listfile = sys.argv[1], sys.argv[2]
with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as z:
    for line in open(listfile, encoding='utf-8'):
        p = line.strip()
        if p:
            z.write(p, p)
PY
fi

listed=$(for f in $PAPER_TEX $PAPER_FIGURES $PAPER_CSV; do echo "$f"; done | wc -l | tr -d ' ')
held=$(ls -1 results | wc -l | tr -d ' ')
echo "Wrote $out: $n of the $held entries of results/, $(wc -c < "$out" | tr -d ' ') bytes."
echo "  $npaper of the $listed files the paper reads$([ -z "$paperonly" ] && echo ", plus the inputs of paper_stats.sh")."
if [ -s "$tmp/missing" ]; then
  echo
  echo "LISTED BUT NOT IN results/ -- $(wc -l < "$tmp/missing" | tr -d ' ') of them:"
  sed 's/^/  /' "$tmp/missing"
  echo
  echo "Each of these is an empty table, a dash or a stale number in the paper. Generate them"
  echo "before trusting the PDF; sh ./bin/run_all_exps.sh --from 6 covers most. If the paper no"
  echo "longer reads one of them, drop its name from the list at the head of this script."
fi
