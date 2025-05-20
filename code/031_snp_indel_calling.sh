#!bin/bash

# Variant calling with GATK HaplotypeCaller
## Refer to https://hpc.nih.gov/training/gatk_tutorial/
## Refer to https://www.melbournebioinformatics.org.au/tutorials/tutorials/variant_calling_gatk1/variant_calling_gatk1/

export PATH=/home/mdnl/Documents/Tam_14Mar25/tools/gatk-4.6.1.0:$PATH

ALIGNED="/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/$proj_n/bwa_aligned"
BAM="$ALIGNED/$SAMPLE/${SAMPLE}.sort.dup.bqsr.bam"
ref="/home/mdnl/Documents/Tam_14Mar25/ref"
FASTA="$ref/Gencode/bwa-mem2-idx/hg38/Homo_sapiens_assembly38.fasta"
gatk_o="/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/$proj_n/gatk"
VCF="$gatk_o/$SAMPLE/${SAMPLE}.snp_indel.vcf.gz"
LOG="$gatk_o/$SAMPLE/031_snp_indel_calling.log"
mkdir -p $gatk_o/$SAMPLE

# Call germline SNPs and short indels using local de novo assembly of haplotype regions
echo "====== `date`: Start gatk HaplotypeCaller ======" > $LOG
CMD="gatk --java-options \"-Xmx8G\" HaplotypeCaller \
    -I $BAM \
    -R $FASTA \
    -ERC GVCF \
    -O $VCF"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish gatk HaplotypeCaller ======" >> $LOG
