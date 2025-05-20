#!bin/bash

export PATH=/home/mdnl/Documents/Tam_14Mar25/tools/gatk-4.6.1.0:$PATH

gatk_o="/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/$proj_n/gatk"
LOG="$gatk_o/034_snp_indel_hard_filter.log"

# Hard filter SNP only
echo "====== `date`: Start extract SNP callsets ======" > $LOG
CMD="gatk --java-options \"-Xmx8G\" SelectVariants \
    -V $gatk_o/cohort.snp.${snp_cutoff}.indel.${indel_cutoff}.vqsr.vcf \
    -select-type SNP \
    -O $gatk_o/cohort.snp_only.${snp_cutoff}.vqsr.vcf"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish extract SNP callsets ======" >> $LOG

echo "====== `date`: Start gatk VariantFiltration for hard filter SNP only ======" >> $LOG
CMD="gatk --java-options \"-Xmx8G\" VariantFiltration \
    -V $gatk_o/cohort.snp_only.${snp_cutoff}.vqsr.vcf \
    $filter_options_SNP \
    -O $gatk_o/cohort.snp.${snp_cutoff}.hard_filter.vcf"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish gatk VariantFiltration for hard filter SNP only ======" >> $LOG

# Hard filter indel only
echo "====== `date`: Start extract INDEL callsets ======" >> $LOG
CMD="gatk --java-options \"-Xmx8G\" SelectVariants \
    -V $gatk_o/cohort.snp.${snp_cutoff}.indel.${indel_cutoff}.vqsr.vcf \
    -select-type INDEL \
    -select-type MIXED \
    -O $gatk_o/cohort.indel_only.${indel_cutoff}.vqsr.vcf"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish extract INDEL callsets ======" >> $LOG

echo "====== `date`: Start gatk VariantFiltration for hard filter INDEL only ======" >> $LOG
CMD="gatk --java-options \"-Xmx8G\" VariantFiltration \
    -V $gatk_o/cohort.indel_only.${indel_cutoff}.vqsr.vcf \
    $filter_options_INDEL \
    -O $gatk_o/cohort.indel.${indel_cutoff}.hard_filter.vcf"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish gatk VariantFiltration for hard filter INDEL only ======" >> $LOG
