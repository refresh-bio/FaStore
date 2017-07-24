#!/bin/bash

set -e


# path to binaries
#
FASTORE_PACK=./fastore_pack
#
#


print_usage()
{
	echo "usage: bash $0 --in <archive> --out <out.fastq> [--pair <pair.fastq>] [--threads <th>]"
	echo ""
	echo "info: when decompressing files compressed in paired-end mode, two output files (out and pair) need to be specified"
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
		--in)
			IN_PFX="$2"
			shift 2;;
		--pair)
			PAIR="$2"
			PAR_PE="-z"
			shift 2;;
		--out)
			OUT="$2"
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
if [ -z ${IN_PFX+x} ]; then
	echo "ERROR: no input files prefix have been specified"
	exit 1
fi

if [ -z ${THREADS+x} ]; then
	echo "WARN: number of threads mode has not been specified --> setting to 1"
	THREADS=1
fi

if [ -z ${OUT+x} ]; then
	OUT="OUT.fastq"
	echo "WARN: no output files name has been specified --> setting to $OUT"
fi



# run
#
if [ -z ${PAIR+x} ]; then
	OUT_FQ="$OUT"
else
	OUT_FQ="$OUT $PAIR"
fi

echo "decompressing: $IN_PFX --> $OUT_FQ"

$FASTORE_PACK d "-i$IN_PFX" "-o$OUT_FQ" "-t$THREADS" $PAR_PE
