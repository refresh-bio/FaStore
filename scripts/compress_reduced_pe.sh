#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "usage: bash $0 <in_fastq_1> <in_fastq_2> <out_fastore_prefix> <threads>"
    exit
fi

set -e
#set -v

# processing params
#
PAR_BIN_C1="-p8 -s0 -b256"
PAR_REBIN_C1="-r -w1024 -W1024"
PAR_PACK_C1="-r -f256 -c10 -d8 -w1024 -W1024"
#PAR_GZ="-g"

TH_BIN=$4
TH_REBIN=$4
TH_PACK=$4

PAR_ID="-H -C"
PAR_QUA="-q2"
PAR_PE="-z"


# temporary and output files configuration
#
TMP_PFX="__tmp-dna"
OUT_PFX="$3"

TMP_BIN="$TMP_PFX-bin_pe"
TMP_REBIN="$TMP_PFX-rebin_pe"

OUT_PACK="$OUT_PFX"

#IN_1=$(ls *_1.fastq.gz)
#IN_2=$(ls *_2.fastq.gz)
#IN="$@"
IN="$1 $2"

echo "processing files: $IN"

echo "- binning ..."
./fastore_bin e "-i$IN" "-o$TMP_BIN" "-t$TH_BIN" $PAR_ID $PAR_QUA $PAR_BIN_C1 $PAR_PE $PAR_GZ

echo "- rebinning: 0 -> 2 ..."
./fastore_rebin e "-i$TMP_BIN" "-o$TMP_REBIN-2" "-t$TH_REBIN" $PAR_REBIN_C1 $PAR_PE -p2
rm $TMP_BIN*

echo "- rebinning: 2 -> 4 ..."
./fastore_rebin e "-i$TMP_REBIN-2" "-o$TMP_REBIN-4" "-t$TH_REBIN" $PAR_REBIN_C1 $PAR_PE -p4
rm $TMP_REBIN-2*

echo "- rebinning: 4 -> 8 ..."
./fastore_rebin e "-i$TMP_REBIN-4" "-o$TMP_REBIN-8" "-t$TH_REBIN" $PAR_REBIN_C1 $PAR_PE -p8
rm $TMP_REBIN-4*

echo "- packing..."
./fastore_pack e "-i$TMP_REBIN-8" "-o$OUT_PACK" "-t$TH_PACK" $PAR_PACK_C1 $PAR_PE -v 2>__err.log
rm $TMP_REBIN-8*
