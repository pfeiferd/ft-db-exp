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

# DB Build Performance
for db in viral tick-borne protozoa parasites vineyard;
  do
    mvn exec:exec@db -Dname=$db -Dgoal=refseqfna
    mvn exec:exec@db -Dname=$db -Dgoal=fastasgenbankdl
    # Performance for creating unrefined db:
    ./tools/cgmemtime/cgmemtime mvn exec:exec@db -Dname=$db -Dgoal=db  > ${res_path}/db_gen_${db}.log
    # Performance for creating refined db on top:
    ./tools/cgmemtime/cgmemtime mvn exec:exec@db -Dname=$db -Dgoal=ftdb  > ${res_path}/ftdb_gen_${db}.log
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
for db in viral tick-borne protozoa parasites vineyard;
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
