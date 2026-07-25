#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

mkdir -p results
mkdir -p results/logs

res_path=./results

for db in viral tick-borne protozoa parasites vineyard;
  do
    mvn exec:exec@db -Dname=$db -Dgoal=clear
    mvn exec:exec@db -Dname=$db -Dgoal=refseqfna
    mvn exec:exec@db -Dname=$db -Dgoal=fastasgenbankdl
    # Performance for creating unrefined db:
    ./cgmemtime/cgmemtime mvn exec:exec@db -Dname=$db -Dgoal=db  > ${res_path}/logs/db_gen_${db}.log
    # Performance for creating refined db on top:
    ./cgmemtime/cgmemtime mvn exec:exec@db -Dname=$db -Dgoal=ftdb  > ${res_path}/logs/ftdb_gen_${db}.log
    mvn exec:exec@db -Dname=$db -Dgoal=dbinfo
    mvn exec:exec@db -Dname=$db -Dgoal=dbquality
    mvn exec:exec@db -Dname=$db -Dgoal=ftdbinfo
    mvn exec:exec@db -Dname=$db -Dgoal=ftquality
  done