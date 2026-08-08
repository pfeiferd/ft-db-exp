#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

mkdir -p results

res_path=./results

# For figures in "Methods" section:
mvn exec:exec@orthopox3 -Dname=orthopox -Dgoal=svgtaxtree
mvn exec:exec@orthopox3 -Dname=orthopox -Dgoal=ftsvgtaxtree
mv data/projects/orthopox/csv/orthopox_svgtaxtree.svg ${res_path}/orthopox3_svgtaxtree.svg
mv data/projects/orthopox/csv/orthopox_ftsvgtaxtree.svg ${res_path}/orthopox3_ftsvgtaxtree.svg
mvn exec:exec@db -Dname=orthopox -Dgoal=dbinfo
mvn exec:exec@db -Dname=orthopox -Dgoal=ftdbinfo
mvn exec:exec@db -Dname=orthopox -Dgoal=intersectcsv
mvn exec:exec@db -Dname=orthopox -Dgoal=dendrolatex
mvn exec:exec@db -Dname=orthopox -Dgoal=svgtaxtree
mvn exec:exec@db -Dname=orthopox -Dgoal=ftsvgtaxtree
mvn exec:exec@db -Dname=orthopox2 -Dgoal=svgtaxtree
mvn exec:exec@db -Dname=orthopox2 -Dgoal=ftsvgtaxtree
mvn exec:exec@db -Dname=borrelia -Dgoal=ftsvgtaxtree

# Times the generation of one database and records it, but only if that database is not there yet.
# Genestrip treats an existing database as already made, so re-running the goal does nothing at all:
# the measurement would then capture Maven's startup and nothing else -- around a second of wall
# time and a couple of hundred MB of RAM -- and, worse, it would overwrite a log that does hold a
# real measurement. Re-running this script after a completed run must therefore leave the logs
# alone. Delete the database to measure its generation again.
#
# $1 = project name, $2 = goal, either `db' for the unrefined database or `ftdb' for the refinement
timed_db_goal() {
  db=$1
  goal=$2
  case "$goal" in
    db)   dbfile="data/projects/${db}/db/${db}_db.zip";   log="${res_path}/db_gen_${db}.log" ;;
    ftdb) dbfile="data/projects/${db}/db/${db}_ftdb.zip"; log="${res_path}/ftdb_gen_${db}.log" ;;
    *)    echo "timed_db_goal: unknown goal '${goal}'" >&2; return 1 ;;
  esac
  if [ -f "$dbfile" ]; then
    echo "WARNING: not measuring '${goal}' for ${db}: ${dbfile} already exists, so the goal would" >&2
    echo "         do nothing and the timing would be meaningless. ${log} is left untouched." >&2
    echo "         Delete the database first if you want to measure its generation again." >&2
    return 0
  fi
  # Being present is not enough: cgmemtime needs a cgroup it may create, which it cannot do without
  # a systemd user session. Probing it here keeps a broken run from leaving a truncated log behind.
  if ! ./tools/cgmemtime/cgmemtime true >/dev/null 2>&1; then
    echo "WARNING: cgmemtime cannot run here - generating ${db}/${goal} without measuring it." >&2
    mvn exec:exec@db -Dname="$db" -Dgoal="$goal"
    return $?
  fi
  ./tools/cgmemtime/cgmemtime mvn exec:exec@db -Dname="$db" -Dgoal="$goal" > "$log"
}

# DB Build Performance
for db in viral tick-borne protozoa gut-protozoa parasites vineyard cdiff;
  do
    mvn exec:exec@db -Dname=$db -Dgoal=refseqfna
    mvn exec:exec@db -Dname=$db -Dgoal=fastasgenbankdl
    # Performance for creating unrefined db, and for the refinement built on top of it:
    timed_db_goal "$db" db
    timed_db_goal "$db" ftdb
    mvn exec:exec@db -Dname=$db -Dgoal=dbinfo
    mvn exec:exec@db -Dname=$db -Dgoal=ftdbinfo
  done

# Now that the DBs are built:
# Stats for figures in "Introduction" section
mvn exec:exec@db -Dname=viral -Dgoal=kmerrankstatscsv
mvn exec:exec@db -Dname=viral -Dgoal=branchhistorankcsv
mvn exec:exec@db -Dname=tick-borne -Dgoal=kmerrankstatscsv
mvn exec:exec@db -Dname=tick-borne -Dgoal=branchhistorankcsv
mvn exec:exec@db -Dname=protozoa -Dgoal=kmerrankstatscsv
mvn exec:exec@db -Dname=protozoa -Dgoal=branchhistorankcsv

# Intrinsic quality of DBs
for db in viral tick-borne protozoa gut-protozoa parasites vineyard cdiff;
  do
    mvn exec:exec@db -Dname=$db -Dgoal=dbquality
    mvn exec:exec@db -Dname=$db -Dgoal=ftquality
  done

find "data" -type f \( -name "*.csv" -o -name "*.svg" -o -name "*.tex" \) | while read -r file; do
	target="$res_path/$(basename "$file")"
	if [ -e "$target" ]; then
		echo "$0: overwriting $file, $target already exists" >&2
	fi
	cp "$file" "$target"
	echo "$file -> $target"
done

# Disk sizes of the unrefined and the refined databases. This must run while the
# databases still exist, i.e. before clean_all.sh removes them. The script writes
# its CSV to ${res_path} by itself, so it comes after the copy loop above.
./bin/db_disk_sizes.sh

# The statistics the paper's Tables "Genestrip databases..." and "Subtree precision..." state, as
# LaTeX macros. It runs last and reads what the copy loop above has just placed in ${res_path}, so
# dbstats.tex can never describe a different set of CSVs than the ones beside it -- which is the
# whole point of generating it rather than copying the numbers by hand.
./bin/paper_stats.sh "$res_path"

