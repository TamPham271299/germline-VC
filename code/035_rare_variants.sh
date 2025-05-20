#!bin/bash

export PATH=/home/mdnl/Documents/Tam_14Mar25/tools/gatk-4.6.1.0:$PATH

gatk_o="/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/$proj_n/gatk"
LOG="$gatk_o/035_rare_variants.log"

# SNP
echo "====== `date`: Start gatk VariantFiltration for rare SNPs MAF_0.05 ======" > $LOG
CMD="gatk --java-options \"-Xmx8G\" VariantFiltration \
    -V $gatk_o/cohort.snp.${snp_cutoff}.hard_filter.vcf \
    -filter \"AF >= $maf\" --filter-name \"AF_$maf\" \
    -O $gatk_o/cohort.snp.${snp_cutoff}.hard_filter.MAF.0.05.vcf"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish gatk VariantFiltration for rare SNPs MAF_0.05 ======" >> $LOG

# Extract rare SNPs
echo "====== `date`: Start extract rare SNPs MAF_0.05 ======" >> $LOG
CMD="gatk --java-options \"-Xmx8G\" SelectVariants \
    -V $gatk_o/cohort.snp.${snp_cutoff}.hard_filter.MAF.0.05.vcf \
    --exclude-filtered true \
    -O $gatk_o/cohort.snp.${snp_cutoff}.hard_filter.MAF.0.05.rare.vcf"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish extract rare SNPs MAF_0.05 ======" >> $LOG

# Indel
echo "====== `date`: Start gatk VariantFiltration for rare indels MAF_0.05 ======" >> $LOG
CMD="gatk --java-options \"-Xmx8G\" VariantFiltration \
    -V $gatk_o/cohort.indel.${indel_cutoff}.hard_filter.vcf \
    -filter \"AF >= $maf\" --filter-name \"AF_$maf\" \
    -O $gatk_o/cohort.indel.${indel_cutoff}.hard_filter.MAF.0.05.vcf"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish gatk VariantFiltration for rare indels MAF_0.05 ======" >> $LOG

# Extract rare indels
echo "====== `date`: Start extract rare indels MAF_0.05 ======" >> $LOG     
CMD="gatk --java-options \"-Xmx8G\" SelectVariants \
    -V $gatk_o/cohort.indel.${indel_cutoff}.hard_filter.MAF.0.05.vcf \
    --exclude-filtered true \
    -O $gatk_o/cohort.indel.${indel_cutoff}.hard_filter.MAF.0.05.rare.vcf"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish extract rare indels MAF_0.05 ======" >> $LOG
