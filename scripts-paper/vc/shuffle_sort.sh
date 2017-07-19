#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "usage: bash $0 <in_1.fastq> <in_2.fastq>"
    exit
fi

#set -v

# generated output files
OUT_1="shuf-$1"
OUT_2="shuf-$2"

# sorting performance parameters
THREADS=1
MAX_MEM=$(($THREADS * 5))

TMP_DIR=__tmp-data-sort
TMP_FILE=$TMP_DIR/__tmp-shuf

SORT=sort

if [ -e $OUT_1 ]; then rm $OUT_1; fi
if [ -e $OUT_2 ]; then rm $OUT_2; fi

if [ ! -d $TMP_DIR ]; then
    mkdir $TMP_DIR
fi


# concatenate the two paired records into one line and
# store them all in one temporary file, line by line
paste <(cat $1) <(cat $2) | \
	paste - - - - > $TMP_FILE

# parallel sort the records and
# parse the output into valid paired FASTQ records stored in separate files
$SORT -R -T $TMP_DIR -S "$MAX_MEM"G --parallel=$THREADS $TMP_FILE | \
	awk -v o1="$OUT_1" -v o2="$OUT_2" -F'\t' '{OFS="\n"; print $1,$3,$5,$7 >> o1; print $2,$4,$6,$8 >> o2 }'

# cleanup
rm -r $TMP_DIR
