#!/bin/bash

set -e


# path to binaries
#
FASTORE_BIN=./fastore_bin
FASTORE_REBIN=./fastore_rebin
FASTORE_PACK=./fastore_pack
#
#


print_usage()
{
	echo "usage: bash $0 <mode> [--fast] --in <in.fq> [--pair <pair.fq>] --out <archive> [--threads <th>]"
	echo ""
	echo "where: <mode> can be one of: --lossless, --reduced, --lossy or --max"
	exit 1
}


# parse input args
#
if [ $# -eq 0 ]
  then
	print_usage
	exit 1
fi

while [[ $# -gt 0 ]]
do
	ARG="$1"

	case $ARG in
		--lossless)
			MODE=0
			shift 1;;
		--reduced)
			MODE=1
			shift 1;;
		--lossy)
			MODE=2
			shift 1;;
		--max)
			MODE=3
			shift 1;;
		--fast)
			FAST_MODE=1
			shift 1;;
		--in)
			IN="$2"
			shift 2;;
		--pair)
			PAIR="$2"
			PAR_PE="-z"
			shift 2;;
		--out)
			OUT_PFX="$2"
			shift 2;;
		--threads)
			THREADS="$2"
			shift 2;;
		*) 
			echo "Unkown argument: \"$ARG\""
			print_usage
			exit 1;;
	esac
done


# check params
#
if [ -z ${MODE+x} ]; then
	echo "ERROR: compression mode has not been specified"
	exit 1
fi

if [ -z ${IN+x} ]; then
	echo "ERROR: no input files have been specified"
	exit 1
fi

if [ -z ${THREADS+x} ]; then
	echo "WARN: number of threads mode has not been specified --> setting to 1"
	THREADS=1
fi

if [ -z ${OUT_PFX+x} ]; then
	echo "WARN: no output files prefix name has been specified --> setting to OUT"
	OUT_PFX="OUT"
fi


# parse params
#
case $MODE in
	0) PAR_ID="-H"; PAR_QUA="-q0";;
	1) PAR_ID="-H -C"; PAR_QUA="-q2";;
	2) PAR_ID="-H -C"; PAR_QUA="-q3";;
	3) PAR_ID=""; PAR_QUA="-q1";;
esac


# set processing params
#
PAR_BIN_C1="-p8 -s0 -b256"
PAR_REBIN_C1="-r -w1024 -W1024"
PAR_PACK_C1="-r -f256 -c10 -d8 -w1024 -W1024"

PAR_BIN_C0="-p8 -s10 -b256"
PAR_PACK_C0="-f256 -c10 -d8 -w256 -W256"

TH_BIN=$THREADS
TH_REBIN=$THREADS
TH_PACK=$THREADS


# temporary and output files configuration
#
TMP_PFX="__tmp-dna"
TMP_BIN="$TMP_PFX-bin_pe"
TMP_REBIN="$TMP_PFX-rebin_pe"

OUT_PACK="$OUT_PFX"


if [ -z ${PAIR+x} ]; then
	IN_FQ="$IN"
else
	IN_FQ="$IN $PAIR"
fi

#echo "$PAR_ID, $PAR_QUA, $PAR_PE, $IN, $PAIR, $OUT_PFX, $THREADS"



echo "processing files: $IN_FQ"

if [ -z ${FAST_MODE+x} ]; then

	echo "- binning ..."
	$FASTORE_BIN e "-i$IN_FQ" "-o$TMP_BIN" "-t$TH_BIN" $PAR_ID $PAR_QUA $PAR_BIN_C1 $PAR_PE

	echo "- rebinning: 0 -> 2 ..."
	$FASTORE_REBIN e "-i$TMP_BIN" "-o$TMP_REBIN-2" "-t$TH_REBIN" $PAR_REBIN_C1 $PAR_PE -p2
	rm $TMP_BIN*

	echo "- rebinning: 2 -> 4 ..."
	$FASTORE_REBIN e "-i$TMP_REBIN-2" "-o$TMP_REBIN-4" "-t$TH_REBIN" $PAR_REBIN_C1 $PAR_PE -p4
	rm $TMP_REBIN-2*

	echo "- rebinning: 4 -> 8 ..."
	$FASTORE_REBIN e "-i$TMP_REBIN-4" "-o$TMP_REBIN-8" "-t$TH_REBIN" $PAR_REBIN_C1 $PAR_PE -p8
	rm $TMP_REBIN-4*

	echo "- packing..."
	$FASTORE_PACK e "-i$TMP_REBIN-8" "-o$OUT_PACK" "-t$TH_PACK" $PAR_PACK_C1 $PAR_PE -v 2>__err.log
	rm $TMP_REBIN-8*

else

	echo "- binning ..."
	$FASTORE_BIN e "-i$IN_FQ" "-o$TMP_BIN" "-t$TH_BIN" $PAR_ID $PAR_QUA $PAR_BIN_C0 $PAR_PE

	echo "- packing..."
	$FASTORE_PACK e "-i$TMP_BIN" "-o$OUT_PACK" "-t$TH_PACK" $PAR_PACK_C0 $PAR_PE -v 2>__err.log
	rm $TMP_BIN*

fi