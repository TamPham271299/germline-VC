#!bin/bash

export PATH=/home/mdnl/Documents/Tam_14Mar25/tools/gatk-

# ref="/home/mdnl/Documents/Tam_14Mar25/ref"
ref="/home/mdnl/KianaRDS/Tam_14Mar25/ref"
gatk_ref="$ref/Gencode/GATK_resource_bundle/hg38"
FASTA_DICT="$ref/Gencode/bwa-mem2-idx/hg38/Homo_sapiens_assembly38.dict"
gatk_o="/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/$proj_n/gatk"
LOG="$gatk_o/036_evaluate_callset.log"
VCF="$gatk_o/$file_n"

# Compare callset with truth set (dbSNP)
echo "====== `date`: Start gatk CollectVariantCallingMetrics for SNP ======" > $LOG
CMD="gatk --java-options \"-Xmx8G\" CollectVariantCallingMetrics \
    -I $VCF \
    -O ${VCF/.vcf/.eval_metrics} \
    --DBSNP $gatk_ref/Homo_sapiens_assembly38.dbsnp138.vcf \
    -SD $FASTA_DICT"
echo $CMD >> $LOG
eval $CMD >> $LOG 2>&1
echo "====== `date`: Finish gatk CollectVariantCallingMetrics for SNP ======" >> $LOG
