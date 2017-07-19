#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "usage: bash $0 <in_fastore_file> <out_fastq_prefix> <threads>"
    exit
fi

set -e
#set -v

# processing params
#
TH_PACK=$3
PAR_PE="-z"

# input/output files
#
IN="$1"
OUT="$2""_1.fastq $2""_2.fastq"

echo "decompressing: $IN --> $OUT"

./fastore_pack d "-i$IN" "-o$OUT" "-t$TH_PACK" $PAR_PE



