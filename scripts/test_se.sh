#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "usage: bash $0 <in_fastq> <threads>"
    exit
fi

set -e


# test config
#
IN=$1
TH=$2

PACK="__pack"
FQ_OUT="__out"


# run very basic tests
#
echo "--------------------------------"
echo "testing: lossless"
echo "--------------------------------"
bash fastore_compress.sh --lossless --in $IN --out $PACK --threads $TH 
bash fastore_decompress.sh --in $PACK --out $FQ_OUT --threads $TH

bash fastore_compress.sh --lossless --in $IN --out $PACK --threads $TH --fast
bash fastore_decompress.sh --in $PACK --out $FQ_OUT --threads $TH


echo "--------------------------------"
echo "testing: reduced"
echo "--------------------------------"
bash fastore_compress.sh --reduced --in $IN --out $PACK --threads $TH
bash fastore_decompress.sh --in $PACK --out $FQ_OUT --threads $TH

bash fastore_compress.sh --reduced --in $IN --out $PACK --threads $TH --fast
bash fastore_decompress.sh --in $PACK --out $FQ_OUT --threads $TH


echo "--------------------------------"
echo "testing: lossy"
echo "--------------------------------"
bash fastore_compress.sh --lossy --in $IN --out $PACK --threads $TH
bash fastore_decompress.sh --in $PACK --out $FQ_OUT --threads $TH

bash fastore_compress.sh --lossy --in $IN --out $PACK --threads $TH --fast
bash fastore_decompress.sh --in $PACK --out $FQ_OUT --threads $TH


echo "--------------------------------"
echo "testing: max"
echo "--------------------------------"
bash fastore_compress.sh --max --in $IN --out $PACK --threads $TH
bash fastore_decompress.sh --in $PACK --out $FQ_OUT --threads $TH

bash fastore_compress.sh --max --in $IN --out $PACK --threads $TH --fast
bash fastore_decompress.sh --in $PACK --out $FQ_OUT --threads $TH
