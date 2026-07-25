#!/bin/sh
set -e

scriptdir=$(dirname "$0")

### cgmemtime ###

cd $scriptdir/..
basedir=$(pwd)

# gm fix - only build cgmemtime if not already there
if [ ! -d "cgmemtime" ]; then
  wget https://github.com/gsauthof/cgmemtime/archive/refs/heads/master.zip
  unzip master.zip
  mv  cgmemtime-master  cgmemtime
  cd  cgmemtime
  # Let's print to stdout:
  sed -i 's/print_result(stderr/print_result(stdout/' cgmemtime.c
  make
  cd $scriptdir/..
  rm master.zip
fi
