#!bin/bash

export PATH=/home/mdnl/Documents/Tam_14Mar25/tools/gatk-4.6.1.0:$PATH

ref="/home/mdnl/Documents/Tam_14Mar25/ref"
FASTA="$ref/Gencode/bwa-mem2-idx/hg38/Homo_sapiens_assembly38.fasta"
gatk_o="/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/$proj_n/gatk"
target="/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/$proj_n/config/target.tsv"
LOG="$gatk_o/032_snp_indel_jointCalling.log"

# ${SAMPLE}.snp_indel.g.vcf (genomic VCF # VCF). GVCF contain all positions (even no ALT), while VCF contains only variants.
# GVCF required for GenotypeGVCFs
# In GVCF of a sample, <NONE-REF> (ref block) means position with no ALT detected in a single sample, but kept for joint calling (GVCF: input for joint calling)
sampleID_idx=$(head -1 $target| tr "\t" "\n"| nl| grep "sampleID"| awk '{print $1}')
vcf_inputs=""
for SAMPLE in `cut -f$sampleID_idx $target| sed '1d'`; do
    vcf_inputs="$vcf_inputs -V $gatk_o/$SAMPLE/${SAMPLE}.snp_indel.g.vcf.gz"
done

# Joint calling
# Apply CombineGVCFs (combine multiple single sample GVCF files)
echo "====== `date`: Start gatk CombineGVCFs ======" > $LOG
CMD="gatk --java-options \"-Xmx8G\" CombineGVCFs \
    -R $FASTA \
    $vcf_inputs \
    -O $gatk_o/cohort.g.vcf.gz"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish gatk CombineGVCFs ======" >> $LOG
 
# Joint genotyping (0/0, 1/1 - only unphased genotypes)
# GATK uses a modified version (to include multi-allelic variants) to calculate the posterior probability of a non-reference allele
echo "====== `date`: Start gatk GenotypeGVCFs ======" >> $LOG
CMD="gatk --java-options \"-Xmx8G\" GenotypeGVCFs \
    -R $FASTA \
    -V $gatk_o/cohort.g.vcf.gz \
    -O $gatk_o/cohort.vcf.gz"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish gatk GenotypeGVCFs ======" >> $LOG
