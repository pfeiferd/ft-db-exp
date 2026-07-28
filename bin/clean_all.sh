#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

for p in data/projects/*/; do
        name=$(basename "$p")
    mvn exec:exec@cleanall -Dname=${name} -Dgoal=ftgenall
    mvn exec:exec@db -Dname=${name} -Dgoal=ftclear
  done