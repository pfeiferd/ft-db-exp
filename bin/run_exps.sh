#!/bin/sh
#
# Builds the databases of the paper, measures what that costs, and gathers the figures and CSV files
# the paper reads.
#
#   sh ./bin/run_exps.sh                 # everything
#   sh ./bin/run_exps.sh vineyard          # only that project
#   sh ./bin/run_exps.sh vineyard viral     # only those
#
# Naming projects restricts every per-project step to them and skips the steps that belong to other
# projects - the Orthopox figures of the Methods section and the rank statistics of the Introduction.
# Everything after them still runs, since it reads whatever is in results/ and reports what is not
# there rather than failing.
#
# Note that a database which already exists is not rebuilt, and its generation is therefore not
# measured either: delete data/projects/<project>/db first if that is what you are after.
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

mkdir -p results

# Which projects to work on. Naming none means all of them, in the order the paper's tables list.
if [ $# -gt 0 ]; then
  projects="$*"
  restricted=1
else
  projects="viral tick-borne strepto protozoa parasites vineyard"
  restricted=""
fi

# True if project $1 is among the ones to work on.
wanted() {
  case " ${projects} " in
    *" $1 "*) return 0 ;;
    *) return 1 ;;
  esac
}

res_path=./results
# Where the unshortened logs are kept while a goal runs. Only their shortened form is copied on to
# ${res_path}/logs, so that the results hold nothing but readable files.
raw_log_path=./logs

# For figures in "Methods" section. These are about the Orthopox example rather than about any of the
# databases below, so a restricted run leaves them alone.
if [ -z "$restricted" ]; then
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
fi

# Times the generation of one database and records it, but only if that database is not there yet.
# Genestrip treats an existing database as already made, so re-running the goal does nothing at all:
# the measurement would then capture Maven's startup and nothing else -- around a second of wall
# time and a couple of hundred MB of RAM -- and, worse, it would overwrite a log that does hold a
# real measurement. Re-running this script after a completed run must therefore leave the logs
# alone. Delete the database to measure its generation again.
#
# Writes the log $1 to $2 with the intermediate renderings of the progress bar dropped, keeping the
# last one of each line. The bar is on by default and repaints every second, so a run of several
# hours would otherwise leave megabytes of overwritten frames in a file that is meant to be read.
# Only the shortened copy goes to the results, while the raw one stays behind in ${raw_log_path} in
# case a run has to be examined as it really came out.
strip_progress_bar() {
  [ -f "$1" ] || return 0
  sed 's/.*\r//' "$1" > "$2"
}

# Runs a command, keeping its output both on screen and on file: stdout goes to the file $1, and
# everything the command writes to stderr goes to the file $2 *and* to the terminal. The latter is
# what Genestrip logs through, progress bars included, so a run remains as watchable as it was
# before these logs were kept at all. Once the command is done, $2 is shortened into $3, which is
# the copy that ends up among the results.
#
# The status returned is the command's own, not the pipeline's, which would otherwise be tee's and
# would leave a failed run undetected under `set -e'.
#
# $1 = file for stdout, $2 = raw file for stderr, $3 = shortened copy of it, $4... = the command
run_logged() {
  _out=$1
  _raw=$2
  _dest=$3
  shift 3
  _statusfile=$(mktemp)
  { "$@" > "$_out" 2>&3 3>&-; echo $? > "$_statusfile"; } 3>&1 | tee "$_raw"
  _status=$(cat "$_statusfile")
  rm -f "$_statusfile"
  strip_progress_bar "$_raw" "$_dest"
  return $_status
}

# $1 = project name, $2 = goal, either `db' for the unrefined database or `ftdb' for the refinement
timed_db_goal() {
  db=$1
  goal=$2
  case "$goal" in
    db)   dbfile="data/projects/${db}/db/${db}_db.zip";   log="${res_path}/db_gen_${db}.log" ;;
    ftdb) dbfile="data/projects/${db}/db/${db}_ftdb.zip"; log="${res_path}/ftdb_gen_${db}.log" ;;
    *)    echo "timed_db_goal: unknown goal '${goal}'" >&2; return 1 ;;
  esac
  # Genestrip logs through Commons Logging's SimpleLog, which writes to stderr, so redirecting
  # stdout alone captures Maven and cgmemtime but drops everything Genestrip itself reports --
  # including the "Used heap size" that each goal implementing Goal.LogHeapInfo prints before and
  # after it runs, which is the only per-goal account of where the memory of a run goes. It goes to
  # a file of its own so that the measurement stays small and machine-readable, and it is written
  # raw to ${raw_log_path} first: only its shortened form belongs among the results.
  gsraw="${raw_log_path}/${goal}_gen_${db}.genestrip.log"
  gslog="${res_path}/logs/${goal}_gen_${db}.genestrip.log"
  if [ -f "$dbfile" ]; then
    echo "WARNING: not measuring '${goal}' for ${db}: ${dbfile} already exists, so the goal would" >&2
    echo "         do nothing and the timing would be meaningless. ${log} is left untouched." >&2
    echo "         Delete the database first if you want to measure its generation again." >&2
    return 0
  fi
  # The goal is going to run, so the logs of whatever ran here before must not survive it: a reader
  # could not tell them apart afterwards, and a run that fails early would otherwise leave the
  # previous log standing as though it described this one. Only reached once the guard above has
  # let us through, so a preserved measurement is never touched.
  mkdir -p "${raw_log_path}" "${res_path}/logs"
  rm -f "$gsraw" "$gslog"

  # Being present is not enough: cgmemtime needs a cgroup it may create, which it cannot do without
  # a systemd user session. Probing it here keeps a broken run from leaving a truncated log behind.
  if ! ./tools/cgmemtime/cgmemtime true >/dev/null 2>&1; then
    echo "WARNING: cgmemtime cannot run here - generating ${db}/${goal} without measuring it." >&2
    run_logged "$log" "$gsraw" "$gslog" mvn exec:exec@db -Dname="$db" -Dgoal="$goal"
    return $?
  fi
  run_logged "$log" "$gsraw" "$gslog" ./tools/cgmemtime/cgmemtime mvn exec:exec@db -Dname="$db" -Dgoal="$goal"
}

# DB Build Performance
echo "### building: ${projects}"
for db in $projects;
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
for db in viral tick-borne protozoa; do
  if wanted "$db"; then
    mvn exec:exec@db -Dname=$db -Dgoal=kmerrankstatscsv
    mvn exec:exec@db -Dname=$db -Dgoal=branchhistorankcsv
  fi
done

# Intrinsic quality of DBs.
#
# A database that uses file nodes is measured here like any other. DBQualityCountsGoal used to
# refuse one outright ("This goal does not support file nodes"), because it recognised a leaf by
# testing for rank DATA while the fasta reader resolved records to the FILE node below it.
# isLeafNode() now defines a leaf as the deepest artificial node on its branch, which is what both
# halves of the goal actually mean, and the two agree again. Nothing changes for the databases
# below: with no file or id nodes the deepest artificial node *is* the data node. Worth keeping in
# mind before adding a per-assembly database to `projects', which is the case that first needed it.
for db in $projects;
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
./bin/db_disk_sizes.sh $projects

# Wall time and memory of each database generation and of the refinement built on top
# of it, gathered from the cgmemtime logs written above and joined with the disk sizes
# just determined. The paper's refinement performance table reads the resulting CSV
# directly, so the numbers it prints can no longer drift apart from the runs they come
# from. This reads logs only and measures nothing itself, hence it may be re-run at any
# time; it does need db_disk_sizes.sh to have gone first.
./bin/db_gen_perf.sh $projects

# The statistics the paper's Tables "Genestrip databases..." and "Subtree precision..." state, as
# LaTeX macros. It runs last and reads what the copy loop above has just placed in ${res_path}, so
# dbstats.tex can never describe a different set of CSVs than the ones beside it -- which is the
# whole point of generating it rather than copying the numbers by hand.
./bin/paper_stats.sh "$res_path"

