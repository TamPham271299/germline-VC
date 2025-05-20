#!bin/bash

proj_n="WES_Korean"
scripts="/home/mdnl/Documents/Tam_14Mar25/code/WGS_WES/scripts"

export proj_n="$proj_n"

SAMPLE="YUHL937-21"
cd /home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/WES_Korean
mkdir -p raw/$SAMPLE/
mv ~/Downloads/${SAMPLE}_*.fq.gz raw/$SAMPLE/

samples=("YUHL937-21")
for SAMPLE in "${samples[@]}"; do
    echo $SAMPLE

    conda activate trim_galore

    # 011. QC
    env SAMPLE="$SAMPLE" bash "$scripts/011_QC_pe.sh"

    # 012. Trim
    env SAMPLE="$SAMPLE" bash "$scripts/012_trim_pe.sh"

    conda deactivate

    conda activate variant_calling
    # 021. BWA
    env SAMPLE="$SAMPLE" bash "$scripts/021_bwa_mem2_pe.sh"

    # 031. per-sample variant calling (SNP + short indel)
    env SAMPLE="$SAMPLE" bash "$scripts/031_snp_indel_calling.sh"

    conda deactivate

    # Store processed data in KianaRDS
    cp -r raw/$SAMPLE/ ~/KianaRDS/Tam_14May25/MD/variant_calling/WES_Korean/raw/
    cp -r trimmed/$SAMPLE/ ~/KianaRDS/Tam_14May25/MD/variant_calling/WES_Korean/trimmed/
    cp -r bwa_aligned/$SAMPLE/ ~/KianaRDS/Tam_14May25/MD/variant_calling/WES_Korean/bwa_aligned/

    # # Remove processed data in WS to free up space
    # rm -r raw/$SAMPLE/ trimmed/$SAMPLE/ bwa_aligned/$SAMPLE/
done

# 032. Joint calling
bash "$scripts/032_joint_calling.sh"

# 033. Variant filtering using GATK VQSR machine learning
snp_cutoff=90.0
indel_cutoff=90.0
env proj_n="$proj_n" snp_cutoff="$snp_cutoff" indel_cutoff="$indel_cutoff" bash "$scripts/033_variant_filtering.sh"

## 034. Variant filtering using hard filters (Optional)
snp_cutoff=90.0
indel_cutoff=90.0
filter_options_SNP="-filter \"QD < 2.0\" --filter-name \"QD_2\" \
                    -filter \"MQ < 40.0\" --filter-name \"MQ_40\" \
                    -filter \"FS > 60.0\" --filter-name \"FS_60\" \
                    -filter \"SOR > 3.0\" --filter-name \"SOR_3\" \
                    -filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum_-12.5\" \
                    -filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum_-8.0\" \
                    -filter \"QUAL < 30.0\" --filter-name \"QUAL_30\""
filter_options_INDEL="-filter \"QD < 2.0\" --filter-name \"QD_2\" \
                    -filter \"FS > 200.0\" --filter-name \"FS_200\" \
                    -filter \"QUAL < 30.0\" --filter-name \"QUAL_30\" \
                    -filter \"ReadPosRankSum < -20.0\" --filter-name \"ReadPosRankSum_-20.0\""
env proj_n="$proj_n" snp_cutoff="$snp_cutoff" indel_cutoff="$indel_cutoff" filter_options_SNP="$filter_options_SNP" filter_options_INDEL="$filter_options_INDEL" bash "$scripts/034_variant_filtering_hard.sh"

## 035. Rare/Common variants
maf=0.05

env proj_n="$proj_n" maf="$maf" bash "$scripts/035_rare_variants.sh" 

## 036. Evaluate callset
inputs=("cohort.indel.90.0.hard_filter.MAF.0.05.vcf" "cohort.indel.90.0.hard_filter.MAF.0.05.rare.vcf" "cohort.snp.90.0.hard_filter.MAF.0.05.vcf" "cohort.snp.90.0.hard_filter.MAF.0.05.rare.vcf")

for input in "${inputs[@]}"; do
    echo $input
    env proj_n="$proj_n" input="$input" bash "$scripts/036_evaluate_callset.sh"
done

## 041. Annotation
conda activate vep
env proj_n="$proj_n" bash "$scripts/041_vep.sh"
conda deactivate

## 042. Filter VEP annotation
conda activate variant_calling
env proj_n="$proj_n" bash "$scripts/042_vep_filter.sh"
conda deactivate

## 043. Convert VCF into tab_delimited
env proj_n="$proj_n" bash "$scripts/043_vcf_to_tab.sh"

# 04. Data analysis
## 043. SVA (MAF < 0.05)

## 044. GBA (MAF < 0.1)





