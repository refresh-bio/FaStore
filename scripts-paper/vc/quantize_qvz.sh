#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 3 ]; then
    printf "Usage: $0 input_fastq T output_fastq\n"
    printf "T being the distortion level\n"
    exit -1
fi

input_fastq=$1
printf "Input FASTQ file: $input_fastq\n"
T=$2
printf "T: $T\n"
output_fastq=$3
printf "Output FASTQ file: $output_fastq\n"

printf "Checking input FASTQ file $input_fastq ... "
if [ ! -f $input_fastq ]; then printf "did not find input FASTQ file: $input_fastq\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                Executables                                  #
###############################################################################

#                              Modify paths                                   #
###############################################################################

# QVZ2
qvz2="/home/qvz2"

# C script
# Need to run first:
# gcc -o /home/replace_qual_fastq /home/replace_qual_fastq.c
replace_qual_fastq="/home/replace_qual_fastq"


printf "Checking executables ... "
if [ ! -x $qvz2 ]; then printf "did not find $qvz2\n"; exit -1; fi
if [ ! -e $replace_qual_fastq ]; then printf "did not find $replace_qual_fastq\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                  Compress                                   #
###############################################################################

printf "Extracting quality values from FASTQ file ... "
awk 'NR % 4 == 0' $input_fastq > $input_fastq.qual
printf "OK\n"

printf "Running QVZ2 ... "
$qvz2 -t $T -v -u $input_fastq.t$T.qual $input_fastq.qual $input_fastq.t$T.qvz
printf "OK\n"

printf "Constructing new FASTQ file with QVZ'd quality values ... "
$replace_qual_fastq $input_fastq $input_fastq.t$T.qual $output_fastq
printf "OK\n"

