#!/bin/bash
CORE=$1
OUTDIR=$2
perl ../PLR-GEN.pl -1 SRR2822457_1.fastq -2 SRR2822457_2.fastq -ref list.txt -p $CORE -o $OUTDIR > log.txt 2>&1
