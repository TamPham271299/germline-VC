#!/bin/bash

RAW="/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/$proj_n/raw"
R1="$RAW/$SAMPLE/${SAMPLE}_1.fq.gz"
R2="$RAW/$SAMPLE/${SAMPLE}_2.fq.gz"
echo $R1 
echo $R2
LOG="$RAW/$SAMPLE/011_QC_pe.log"

echo "====== `date`: Start fastqc ======" > $LOG
CMD="fastqc --nogroup -t 16 --dir $(dirname $R1) -o $(dirname $R1) $R1 $R2"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish fastqc ======" >> $LOG
