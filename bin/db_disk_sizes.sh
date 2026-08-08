#!/bin/sh
#
# Determines the disk size of the unrefined and the refined (FT) database of each
# project and expresses the refined size as a percentage of the unrefined one.
#
# The results feed the "Disk (MB)" and "Disk (%)" columns of the database refinement
# performance table in the ft-paper. Run this after the databases have been built,
# i.e. after the "DB Build Performance" loop of run_exps.sh, and before clean_all.sh
# removes them again.
#
# Usage: bin/db_disk_sizes.sh [project ...]      (default: the projects of run_exps.sh)

set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

mkdir -p results

res_path=./results
csv=${res_path}/db_disk_sizes.csv

if [ $# -gt 0 ]; then
    projects="$@"
else
    projects="viral tick-borne protozoa parasites vineyard cdiff"
fi

# Echoes the path of the database artifact of project $1 for the goal $2 ("db" or
# "ftdb"), or nothing if none exists. The artifact may be a single file or a folder,
# depending on the Genestrip version, so several candidates are probed.
find_db() {
    _dir=data/projects/$1/db
    for _cand in "$_dir/$1_$2.zip" "$_dir/$1_$2.ser" "$_dir/$1_$2"; do
        if [ -e "$_cand" ]; then
            echo "$_cand"
            return
        fi
    done
}

# Echoes the size of $1 in KiB, or nothing if $1 is empty or missing.
size_kib() {
    if [ -n "$1" ] && [ -e "$1" ]; then
        du -sk "$1" | awk '{print $1}'
    fi
}

echo "project;unrefined path;unrefined MB;refined path;refined MB;refined %" > $csv

printf '%-12s %12s %12s %10s\n' "PROJECT" "UNREF. (MB)" "REFINED (MB)" "REFINED %"

for p in $projects; do
    db_path=$(find_db "$p" db)
    ftdb_path=$(find_db "$p" ftdb)

    db_kib=$(size_kib "$db_path")
    ftdb_kib=$(size_kib "$ftdb_path")

    db_mb=$(awk -v k="$db_kib" 'BEGIN { if (k == "") print ""; else printf "%.0f", k / 1024 }')
    ftdb_mb=$(awk -v k="$ftdb_kib" 'BEGIN { if (k == "") print ""; else printf "%.0f", k / 1024 }')
    pct=$(awk -v a="$db_kib" -v b="$ftdb_kib" \
          'BEGIN { if (a == "" || b == "" || a + 0 == 0) print ""; else printf "%.1f", b * 100 / a }')

    echo "$p;$db_path;$db_mb;$ftdb_path;$ftdb_mb;$pct" >> $csv

    printf '%-12s %12s %12s %10s\n' "$p" "${db_mb:--}" "${ftdb_mb:--}" "${pct:--}"

    if [ -z "$db_path" ]; then
        echo "  warning: no unrefined database found under data/projects/$p/db" >&2
    fi
    if [ -z "$ftdb_path" ]; then
        echo "  warning: no refined database found under data/projects/$p/db" >&2
    fi
    # A folder-shaped artifact stems from an older Genestrip version and may hold the
    # filtering database next to the k-mer database, which would inflate the size.
    for _path in "$db_path" "$ftdb_path"; do
        if [ -n "$_path" ] && [ -d "$_path" ]; then
            echo "  warning: $_path is a folder - check whether it also contains the filtering database" >&2
        fi
    done
done

echo
echo "Sizes are reported as disk usage in MiB (du -sk / 1024)."
echo "Wrote $csv"
