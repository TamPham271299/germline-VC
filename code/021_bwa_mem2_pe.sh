#!bin/bash

export PATH=/home/mdnl/Documents/Tam_14Mar25/tools/bwa-mem2:$PATH
export PATH=/home/mdnl/Documents/Tam_14Mar25/tools/gatk-4.6.1.0:$PATH
picard="/home/mdnl/Documents/Tam_14Mar25/tools/picard-tools"

TRIMMED="/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/$proj_n/trimmed"
ALIGNED="/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/$proj_n/bwa_aligned"
ref="/home/mdnl/Documents/Tam_14Mar25/ref"
FASTA="$ref/Gencode/bwa-mem2-idx/hg38/Homo_sapiens_assembly38.fasta"
gatk_ref="$ref/Gencode/GATK_resource_bundle/hg38"
R1="$TRIMMED/$SAMPLE/${SAMPLE}_1.fq.gz"
R2="$TRIMMED/$SAMPLE/${SAMPLE}_2.fq.gz"
BAM="$ALIGNED/$SAMPLE/${SAMPLE}.bam"
LOG="$ALIGNED/$SAMPLE/021_bwa_mem2_pe.log"
mkdir -p $ALIGNED/$SAMPLE

# Alignment
echo "====== `date`: Start BWA-MEM2 mem ======" > $LOG
CMD="bwa-mem2 mem -M -t 8 \
    -R \"@RG\tID:${SAMPLE}_RG\tSM:${SAMPLE}\tLB:${SAMPLE}_LIB\tPL:ILLUMINA\" \
    $FASTA \
    $R1 \
    $R2 | \
    samtools view -b -h -o $BAM -" 
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish BWA-MEM2 mem ======" >> $LOG

# Sort BAM
echo "====== `date`: Start picard SortSam ======" >> $LOG
CMD="java -Xmx8G -jar $picard/picard.jar SortSam \
    --INPUT $BAM \
    --OUTPUT ${BAM/.bam/.sort.bam} \
    --VALIDATION_STRINGENCY LENIENT \
    --SORT_ORDER coordinate \
    --MAX_RECORDS_IN_RAM 3000000 \
    --CREATE_INDEX true"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish picard SortSam ======" >> $LOG

# Flagstat summary
echo "====== `date`: Start samtools flagstat ======" >> $LOG
CMD="samtools flagstat ${BAM/.bam/.sort.bam}"
echo $CMD >> $LOG
eval $CMD > $ALIGNED/$SAMPLE/${SAMPLE}.summary.txt
echo "====== `date`: Finish samtools flagstat ======" >> $LOG

# Post-alignment processing
## Markduplicates
echo "====== `date`: Start picard MarkDuplicates ======" >> $LOG
CMD="java -Xmx8G -jar $picard/picard.jar MarkDuplicates \
    --QUIET true \
    --INPUT ${BAM/.bam/.sort.bam} \
    --OUTPUT ${BAM/.bam/.sort.dup.bam} \
    --METRICS_FILE $ALIGNED/$SAMPLE/${SAMPLE}.marked_dup.metrics \
    --CREATE_INDEX true \
    --VALIDATION_STRINGENCY LENIENT"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1 
echo "====== `date`: Finish picard MarkDuplicates ======" >> $LOG

## Base recalibrator (Optional)
echo "====== `date`: Start gatk BaseRecalibrator ======" >> $LOG
CMD="gatk --java-options \"-Xmx8G\" BaseRecalibrator \
    -I ${BAM/.bam/.sort.dup.bam} \
    -R $FASTA \
    --known-sites $gatk_ref/Homo_sapiens_assembly38.dbsnp138.vcf \
    -O $ALIGNED/$SAMPLE/${SAMPLE}.recal_data.table"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish gatk BaseRecalibrator ======" >> $LOG

echo "====== `date`: Start gatk ApplyBQSR ======" >> $LOG
CMD="gatk --java-options \"-Xmx8G\" ApplyBQSR \
    -I ${BAM/.bam/.sort.dup.bam} \
    -R $FASTA \
    --bqsr-recal-file $ALIGNED/$SAMPLE/${SAMPLE}.recal_data.table \
    -O ${BAM/.bam/.sort.dup.bqsr.bam}"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish gatk ApplyBQSR ======" >> $LOG
