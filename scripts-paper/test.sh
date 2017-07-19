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
FQ_OUT="__out"


# run very basic tests
#
echo "--------------------------------"
echo "testing: lossless"
echo "--------------------------------"
bash compress_lossless_pe.sh $IN_1 $IN_2 $PACK $TH
bash decompress_pe.sh $PACK $FQ_OUT $TH


echo "--------------------------------"
echo "testing: reduced"
echo "--------------------------------"
bash compress_reduced_pe.sh $IN_1 $IN_2 $PACK $TH
bash decompress_pe.sh $PACK $FQ_OUT $TH


echo "--------------------------------"
echo "testing: lossy"
echo "--------------------------------"
bash compress_lossy_pe.sh $IN_1 $IN_2 $PACK $TH
bash decompress_pe.sh $PACK $FQ_OUT $TH


echo "--------------------------------"
echo "testing: max"
echo "--------------------------------"
bash compress_max_pe.sh $IN_1 $IN_2 $PACK $TH
bash decompress_pe.sh $PACK $FQ_OUT $TH