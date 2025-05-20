#!bin/bash

export PATH=/home/mdnl/Documents/Tam_14Mar25/tools/gatk-4.6.1.0:$PATH

ref="/home/mdnl/Documents/Tam_14Mar25/ref"
FASTA="$ref/Gencode/bwa-mem2-idx/hg38/Homo_sapiens_assembly38.fasta"
gatk_ref="$ref/Gencode/GATK_resource_bundle/hg38"
# -rw-rw-r-- 1 mdnl mdnl 2.1M Jul 22  2016 1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
# -rw-rw-r-- 1 mdnl mdnl 1.8G Jul 22  2016 1000G_phase1.snps.high_confidence.hg38.vcf.gz
# -rw-rw-r-- 1 mdnl mdnl 1.5M Jul 22  2016 1000G_omni2.5.hg38.vcf.gz.tbi
# -rw-rw-r-- 1 mdnl mdnl  51M Jul 22  2016 1000G_omni2.5.hg38.vcf.gz
# -rw-rw-r-- 1 mdnl mdnl 1.5M Jul 22  2016 hapmap_3.3.hg38.vcf.gz.tbi
# -rw-rw-r-- 1 mdnl mdnl  60M Jul 22  2016 hapmap_3.3.hg38.vcf.gz
# -rw-rw-r-- 1 mdnl mdnl 1.5M Jul 22  2016 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
# -rw-rw-r-- 1 mdnl mdnl  20M Jul 22  2016 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
# -rw-rw-r-- 1 mdnl mdnl  12M Jul 22  2016 Homo_sapiens_assembly38.dbsnp138.vcf.idx
# -rw-rw-r-- 1 mdnl mdnl  11G Jul 22  2016 Homo_sapiens_assembly38.dbsnp138.vcf
gatk_o="/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/$proj_n/gatk"
LOG="$gatk_o/033_snp_indel_VQSR.log"

# Filter SNP
## Step 1 - Recalibrating SNPs in exome data
### DP should not be included in annotation for recalibration when working with exome
### To achieve the best exome result, need to use an exome SNP/indel callset with at least 30 samples.
### Input VCF need to be annotated with QD, FS, DP ... (otherwise should run VariantAnnotator)
echo "====== `date`: Start gatk VariantRecalibrator (SNP mode) ======" > $LOG
CMD="gatk --java-options \"-Xmx8G\" VariantRecalibrator \
    -V $gatk_o/cohort.vcf.gz \
    --max-gaussians 4 \
    -mode SNP \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $gatk_ref/hapmap_3.3.hg38.vcf.gz \
    --resource:omni,known=false,training=true,truth=true,prior=12.0 $gatk_ref/1000G_omni2.5.hg38.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 $gatk_ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $gatk_ref/Homo_sapiens_assembly38.dbsnp138.vcf \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR \
    -O $gatk_o/cohort.snp.recal \
    --tranches-file $gatk_o/cohort.snp.tranches"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish gatk VariantRecalibrator (SNP mode) ======" >> $LOG

## Step 2 - ApplyVQSR
echo "====== `date`: Start gatk ApplyVQSR (SNP mode) ======" >> $LOG
CMD="gatk --java-options \"-Xmx8G\" ApplyVQSR \
    -R $FASTA \
    -V $gatk_o/cohort.vcf.gz \
    -O $gatk_o/cohort.snp.${snp_cutoff}.vqsr.vcf \
    --truth-sensitivity-filter-level $snp_cutoff \
    --tranches-file $gatk_o/cohort.snp.tranches \
    --recal-file $gatk_o/cohort.snp.recal \
    -mode SNP \
    --create-output-variant-index true"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish gatk ApplyVQSR (SNP mode) ======" >> $LOG

# Filter indels
## Step 1 - Recalibrating indels + mixed (auto detected as indels)
echo "====== `date`: Start gatk VariantRecalibrator (short INDEL mode) ======" >> $LOG
CMD="gatk --java-options \"-Xmx8G\" VariantRecalibrator \
    -V $gatk_o/cohort.vcf.gz \
    --max-gaussians 4 \
    -mode INDEL \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 $gatk_ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $gatk_ref/Homo_sapiens_assembly38.dbsnp138.vcf \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR \
    -O $gatk_o/cohort.indel.recal \
    --tranches-file $gatk_o/cohort.indel.tranches"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish gatk VariantRecalibrator (short INDEL mode) ======" >> $LOG

## Step 2 - ApplyVQSR
echo "====== `date`: Start gatk ApplyVQSR (short INDEL mode) ======" >> $LOG
CMD="gatk --java-options \"-Xmx8G\" ApplyVQSR \
    -R $FASTA \
    -V $gatk_o/cohort.snp.${snp_cutoff}.vqsr.vcf \
    -O $gatk_o/cohort.snp.${snp_cutoff}.indel.${indel_cutoff}.vqsr.vcf \
    --truth-sensitivity-filter-level $indel_cutoff \
    --tranches-file $gatk_o/cohort.indel.tranches \
    --recal-file $gatk_o/cohort.indel.recal \
    -mode INDEL \
    --create-output-variant-index true"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish gatk ApplyVQSR (short INDEL mode) ======" >> $LOG
