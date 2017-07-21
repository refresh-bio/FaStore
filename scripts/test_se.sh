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
bash compress_lossless_se.sh $IN $PACK $TH
bash decompress_se.sh $PACK $FQ_OUT $TH


echo "--------------------------------"
echo "testing: reduced"
echo "--------------------------------"
bash compress_reduced_se.sh $IN $PACK $TH
bash decompress_se.sh $PACK $FQ_OUT $TH


echo "--------------------------------"
echo "testing: lossy"
echo "--------------------------------"
bash compress_lossy_se.sh $IN $PACK $TH
bash decompress_se.sh $PACK $FQ_OUT $TH


echo "--------------------------------"
echo "testing: max"
echo "--------------------------------"
bash compress_max_se.sh $IN $PACK $TH
bash decompress_se.sh $PACK $FQ_OUT $TH