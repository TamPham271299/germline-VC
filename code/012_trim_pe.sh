#!/bin/bash

RAW="/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/$proj_n/raw"
TRIMMED="/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/$proj_n/trimmed"
R1="$RAW/$SAMPLE/${SAMPLE}_1.fq.gz"
R2="$RAW/$SAMPLE/${SAMPLE}_2.fq.gz"
mkdir -p $TRIMMED/$SAMPLE
LOG="$TRIMMED/$SAMPLE/012_trim_pe.log"

echo "====== `date`: Start trim-galore ======" > $LOG
CMD="trim_galore --cores 4 -q 30 --illumina --length 20 --paired --fastqc_args \"--nogroup\" --output_dir $TRIMMED/$SAMPLE $R1 $R2"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish trim-galore ======" >> $LOG

echo "====== `date`: Start rename file ======" >> $LOG
CMD='for i in `find $TRIMMED/$SAMPLE -name "*_val_*.fq.gz"`; do mv "$i" "${i/_val_?.fq.gz/.fq.gz}"; done'
echo "$CMD" >> $LOG
eval "$CMD" >> $LOG 2>&1
echo "====== `date`: Finish rename file ======" >> $LOG
