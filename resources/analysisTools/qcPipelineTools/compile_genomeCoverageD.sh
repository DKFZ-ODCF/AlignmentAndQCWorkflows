#!/usr/bin/env bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#


​# This is a short script to compile the coverageQcD tool.
#
​# ldc2-0.12.1-linux-x86_64/bin/rdmd --compiler=ldc2-0.12.1-linux-x86_64/bin/ldmd2 --build-only -IBioD-master -O -release -version=1 -inline genomeCoverage/genomeCoverage.d
​
# Path to rdmd and ldmd2 in the LDC 0.12.1 installation (https://github.com/ldc-developers/ldc/releases/tag/v0.12.1)
RDMD="ldc2-0.12.1-linux-x86_64/bin/rdmd"
LDMD2="​ldc2-0.12.1-linux-x86_64/bin/ldmd2"

# Path to checked out BioD master branch (commit 8b633de)
BIOD="BioD-master"

​usage="Usage: $0 [options]"
​counter="0"
​while getopts ":d:l:b:o:h" options
​do
​  case $options in
​  d ) dcomp=$OPTARG
​    ;;
​  l ) ldmd=$OPTARG
​    ;;
​  b ) bioD=$OPTARG
​    ;;
​  o ) outfile=$OPTARG
​    ;;
​  h ) echo $usage
​    echo -e "OPTIONS:\n\t-d\tThe rdmd compiler, default=$RDMD\n\t-l\tThe ldmd2 compiler, default=$LDMD2\n\t-b\tPath to the BioD version, default=$BIOD\n\t-o\tOutput binary (please use a descriptive name and version)"
​    exit 0
​    ;;
​   \? ) echo $usage
​    echo -e "OPTIONS:\n\t-d\tThe rdmd compiler, default=$RDMD\n\t-l\tThe ldmd2 compiler, default=$LDMD2\n\t-b\tPath to the BioD version, default=$BIOD\n\t-o\tOutput binary (please use a descriptive name and version)"
​    exit 1
​    ;;
​  esac
​done
​
​if [[ -z "$dcomp" ]]
​then
​  dcomp="$RDMD"
​fi
​
​if [[ -z "$ldmd" ]]
​then
​  ldmd="$LDMD2"
​fi
​
​if [[ -z "$bioD" ]]
​then
​  bioD=$BIOD
​fi
​
​touch $outfile
​
​if [[ ! -f "$dcomp" || ! -x "$dcomp" || ! -f "$ldmd" || ! -x "$ldmd" || ! -d "$bioD" || ! -f "$outfile" ]]
​then
​  echo "Wrong usage, the right way is:"
​  echo $usage
​  echo -e "OPTIONS:\n\t-d\tThe rdmd compiler, default=$RDMD\n\t-l\tThe ldmd2 compiler, default=$LDMD2\n\t-b\tPath to the BioD version, default=$BIOD\n\t-o\tOutput binary (please use a descriptive name and version)"
​  rm $outfile
​  exit 1
​fi
​
​if [[ ! -f "genomeCoverage/genomeCoverage.d" ]]
​then
​  echo "The D script coverageQcD/coverageQc.d was not found. Make sure it is placed in the right subdirectory:"
​  echo $usage
​  echo -e "OPTIONS:\n\t-d\tThe rdmd compiler, default=$RDMD\n\t-l\tThe ldmd2 compiler, default=$LDMD2\n\t-b\tPath to the BioD version, default=$BIOD\n\t-o\tOutput binary (please use a descriptive name and version)"
​  rm $outfile
​  exit 2
​fi
​
​$dcomp --compiler=$ldmd --build-only -I${bioD} -O -release -inline genomeCoverage/genomeCoverage.d
​
​if [[ "$?" == 0 ]]
​then
​  mv genomeCoverage/genomeCoverage $outfile && chmod 755 $outfile && echo "CoverageQc successfully compiled and moved to $outfile" && exit 0
​else
​  rm $outfile
​  echo "An error occured while compiling coverageQc."
​  exit 3
​fi

