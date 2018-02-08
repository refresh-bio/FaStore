#!/bin/bash

set -e


# path to binaries
#
FASTORE_PACK=./fastore_pack
#
#

log()
{
	if [ ! -z ${VERBOSE+x} ]; then
		printf "$1\n"
	fi
}

print_usage()
{
	echo ""
	echo "    FaStore -- a space saving solution for raw sequencing data"
	echo ""
	echo "usage: bash $0 --in <archive> --out <out.fastq> [--pair <pair.fastq>] [--threads <th>]"
	echo ""
	echo "where:"
	echo "    --in <archive>    - archive filename"
	echo "    --out <out.fastq> - the output FASTQ filename"
	echo "    --pair <pair.fq>  - the paired output FASTQ filename (for paired-end archives only)"
	echo "    --threads <th>    - the number of processing threads"
	echo "    --verbose         - print additional information while decompressing"
	echo "    --help            - displays this message"
	echo ""
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
		--verbose)
			VERBOSE=1
			shift 1;;
		--help|*) 
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
	log "WARN: number of threads mode has not been specified --> setting to 1"
	THREADS=1
fi

if [ -z ${OUT+x} ]; then
	OUT="OUT.fastq"
	log "WARN: no output files name has been specified --> setting to $OUT"
fi



# run
#
if [ -z ${PAIR+x} ]; then
	OUT_FQ="$OUT"
else
	OUT_FQ="$OUT $PAIR"
fi

log ":: decompressing: $IN_PFX --> $OUT_FQ"

$FASTORE_PACK d "-i$IN_PFX" "-o$OUT_FQ" "-t$THREADS" $PAR_PE

log "decompressed:"
log "$(ls -s $OUT_FQ)"