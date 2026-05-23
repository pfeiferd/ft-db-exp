#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..
basedir=$(pwd)

### Genestrip DBs ###
cd $basedir


mvn exec:exec@orthopox3 -Dname=orthopox -Dgoal=svgtaxtree
mvn exec:exec@orthopox3 -Dname=orthopox -Dgoal=ftsvgtaxtree
mv data/projects/orthopox/csv/orthopox_svgtaxtree.svg data/projects/orthopox/csv/orthopox3_svgtaxtree.svg
mv data/projects/orthopox/csv/orthopox_ftsvgtaxtree.svg data/projects/orthopox/csv/orthopox3_ftsvgtaxtree.svg

mvn exec:exec@db -Dname=orthopox -Dgoal=dbinfo
mvn exec:exec@db -Dname=orthopox -Dgoal=ftdbinfo
mvn exec:exec@db -Dname=orthopox -Dgoal=intersectcsv
mvn exec:exec@db -Dname=orthopox -Dgoal=dendrolatex
mvn exec:exec@db -Dname=orthopox -Dgoal=svgtaxtree
mvn exec:exec@db -Dname=orthopox -Dgoal=ftsvgtaxtree

mvn exec:exec@db -Dname=orthopox2 -Dgoal=svgtaxtree
mvn exec:exec@db -Dname=orthopox2 -Dgoal=ftsvgtaxtree
