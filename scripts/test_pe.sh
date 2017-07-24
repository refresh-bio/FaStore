#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "usage: bash $0 <in_fastq_1> <in_fastq_2> <threads>"
    exit
fi

set -e


# test config
#
IN_1=$1
IN_2=$2
TH=$3

PACK="__pack"
OUT_1="__out_1.fastq"
OUT_2="__out_2.fastq"



# run very basic tests
#
echo "--------------------------------"
echo "testing: lossless"
echo "--------------------------------"
bash compress.sh --lossless --in $IN_1 --pair $IN_2 --out $PACK --threads $TH
bash decompress.sh --in $PACK --out $OUT_1 --pair $OUT_2 --threads $TH

bash compress.sh --lossless --in $IN_1 --pair $IN_2 --out $PACK --threads $TH --fast
bash decompress.sh --in $PACK --out $OUT_1 --pair $OUT_2 --threads $TH


echo "--------------------------------"
echo "testing: reduced"
echo "--------------------------------"
bash compress.sh --reduced --in $IN_1 --pair $IN_2 --out $PACK --threads $TH
bash decompress.sh --in $PACK --out $OUT_1 --pair $OUT_2 --threads $TH

bash compress.sh --reduced --in $IN_1 --pair $IN_2 --out $PACK --threads $TH --fast
bash decompress.sh --in $PACK --out $OUT_1 --pair $OUT_2 --threads $TH


echo "--------------------------------"
echo "testing: lossy"
echo "--------------------------------"
bash compress.sh --lossy --in $IN_1 --pair $IN_2 --out $PACK --threads $TH
bash decompress.sh --in $PACK --out $OUT_1 --pair $OUT_2 --threads $TH

bash compress.sh --lossy --in $IN_1 --pair $IN_2 --out $PACK --threads $TH --fast
bash decompress.sh --in $PACK --out $OUT_1 --pair $OUT_2 --threads $TH


echo "--------------------------------"
echo "testing: max"
echo "--------------------------------"
bash compress.sh --max --in $IN_1 --pair $IN_2 --out $PACK --threads $TH
bash decompress.sh --in $PACK --out $OUT_1 --pair $OUT_2 --threads $TH

bash compress.sh --max --in $IN_1 --pair $IN_2 --out $PACK --threads $TH --fast
bash decompress.sh --in $PACK --out $OUT_1 --pair $OUT_2 --threads $TH