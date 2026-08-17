#!/bin/sh
#
# Collects the resource consumption of generating each unrefined database and of the
# subsequent refinement into one CSV, and expresses the refinement as a percentage of
# the generation it builds upon.
#
# The figures come from the cgmemtime logs that run_exps.sh writes while it builds the
# databases (results/db_gen_<project>.log and results/ftdb_gen_<project>.log) and, for
# the disk columns, from results/db_disk_sizes.csv. This script only reads them, so it
# can be re-run at any time without re-measuring anything; but it cannot invent a
# measurement that was never taken, and it leaves the corresponding fields empty then.
#
# The results feed the database refinement performance table of the ft-paper, which
# reads this CSV directly instead of restating its numbers. Run it after the databases
# have been built and after db_disk_sizes.sh, whose CSV it consumes.
#
# Two memory figures are reported per run. `ram' is cgmemtime's child_RSS_high, the
# high-water resident set of the measured process itself, and is what the paper's table
# states. `group_ram' is its group_mem_high, the high-water mark of the whole process
# group including the page cache, which is much larger and is kept here so that the
# choice stays visible rather than being buried in this script.
#
# Usage: bin/db_gen_perf.sh [project ...]   (default: the projects of the paper's table)

set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

mkdir -p results

res_path=./results
csv=${res_path}/db_gen_perf.csv
sizes=${res_path}/db_disk_sizes.csv

if [ $# -gt 0 ]; then
    projects="$@"
else
    # In the order the paper's table lists them.
    projects="viral tick-borne protozoa vineyard parasites"
fi

# Echoes the short name the paper uses for project $1, or the project name itself.
label_of() {
    case "$1" in
        viral)      echo "cv" ;;
        tick-borne) echo "tb" ;;
        *)          echo "$1" ;;
    esac
}

# Echoes the value of the cgmemtime field $2 in the log $1, or nothing if either is
# missing. cgmemtime appends its fields after the output of the command it measured,
# one per line, as `<field>: <value> <unit>'.
field_of() {
    if [ -f "$1" ]; then
        grep -m1 "^$2:" "$1" 2>/dev/null | awk '{print $2}'
    fi
}

# Echoes $1 KiB as MiB, rounded, or nothing if $1 is empty.
as_mb() {
    awk -v k="$1" 'BEGIN { if (k == "") print ""; else printf "%.0f", k / 1024 }'
}

# Echoes $1 seconds as minutes with one decimal, or nothing if $1 is empty.
as_min() {
    awk -v s="$1" 'BEGIN { if (s == "") print ""; else printf "%.1f", s / 60 }'
}

# Echoes 100 * $2 / $1 with one decimal, or nothing if either is empty or $1 is zero.
as_pct() {
    awk -v a="$1" -v b="$2" \
        'BEGIN { if (a == "" || b == "" || a + 0 == 0) print ""; else printf "%.1f", b * 100 / a }'
}

# The seven columns the paper's table reads come first: csvsimple addresses columns by
# number only up to the tenth, so anything the table needs has to stay within that range.
echo "project;label;db_wall_min;db_ram_mb;db_disk_mb;ft_wall_pct;ft_ram_pct;ft_disk_pct;ft_wall_min;ft_ram_mb;db_group_ram_mb;ft_group_ram_mb" > $csv

printf '%-12s %10s %10s %10s %9s %8s %8s\n' \
       "PROJECT" "WALL(min)" "RAM(MB)" "DISK(MB)" "WALL(%)" "RAM(%)" "DISK(%)"

for p in $projects; do
    db_log=${res_path}/db_gen_${p}.log
    ft_log=${res_path}/ftdb_gen_${p}.log

    db_wall=$(field_of "$db_log" wall)
    db_rss=$(field_of "$db_log" child_RSS_high)
    db_grp=$(field_of "$db_log" group_mem_high)
    ft_wall=$(field_of "$ft_log" wall)
    ft_rss=$(field_of "$ft_log" child_RSS_high)
    ft_grp=$(field_of "$ft_log" group_mem_high)

    # The disk columns are the ones db_disk_sizes.sh has already determined.
    if [ -f "$sizes" ]; then
        db_disk=$(awk -F';' -v p="$p" '$1 == p { print $3 }' "$sizes")
        ft_disk_pct=$(awk -F';' -v p="$p" '$1 == p { print $6 }' "$sizes")
    else
        db_disk=""
        ft_disk_pct=""
    fi

    db_wall_min=$(as_min "$db_wall")
    ft_wall_min=$(as_min "$ft_wall")
    db_ram_mb=$(as_mb "$db_rss")
    ft_ram_mb=$(as_mb "$ft_rss")
    db_group_mb=$(as_mb "$db_grp")
    ft_group_mb=$(as_mb "$ft_grp")
    ft_wall_pct=$(as_pct "$db_wall" "$ft_wall")
    ft_ram_pct=$(as_pct "$db_rss" "$ft_rss")

    echo "$p;$(label_of "$p");$db_wall_min;$db_ram_mb;$db_disk;$ft_wall_pct;$ft_ram_pct;$ft_disk_pct;$ft_wall_min;$ft_ram_mb;$db_group_mb;$ft_group_mb" >> $csv

    printf '%-12s %10s %10s %10s %9s %8s %8s\n' "$p" \
           "${db_wall_min:--}" "${db_ram_mb:--}" "${db_disk:--}" \
           "${ft_wall_pct:--}" "${ft_ram_pct:--}" "${ft_disk_pct:--}"

    if [ ! -f "$db_log" ]; then
        echo "  warning: $db_log is missing - the database was never measured" >&2
    fi
    if [ ! -f "$ft_log" ]; then
        echo "  warning: $ft_log is missing - the refinement was never measured" >&2
    fi
    if [ -f "$db_log" ] && [ -z "$db_wall" ]; then
        echo "  warning: $db_log holds no cgmemtime fields - it was written without measuring" >&2
    fi
    if [ -f "$ft_log" ] && [ -z "$ft_wall" ]; then
        echo "  warning: $ft_log holds no cgmemtime fields - it was written without measuring" >&2
    fi
done

echo
echo "Memory is cgmemtime's child_RSS_high in MiB; group_mem_high is carried along as well."
echo "Wrote $csv"
