#!/bin/bash

LINES=10000
OUT_1=test_1.fq
OUT_2=test_2.fq

curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR024/ERR024163/ERR024163_1.fastq.gz 2>/dev/null | gunzip | head -n $LINES > $OUT_1
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR024/ERR024163/ERR024163_2.fastq.gz 2>/dev/null | gunzip | head -n $LINES > $OUT_2