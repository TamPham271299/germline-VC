#!bin/bash

proj_n="WES_Korean"

export proj_n="$proj_n"

# Specify samples
# samples=("YUHL937-21")
for SAMPLE in "${samples[@]}"; do
    echo $SAMPLE
    mkdir -p raw/$SAMPLE/
    
    # 011. QC
    env SAMPLE="$SAMPLE" bash "code/011_QC_pe.sh"

    # 012. Trim
    env SAMPLE="$SAMPLE" bash "code/012_trim_pe.sh"

    # 021. BWA
    env SAMPLE="$SAMPLE" bash "code/021_bwa_mem2_pe.sh"

    # 031. per-sample variant calling (SNP + short indel)
    env SAMPLE="$SAMPLE" bash "code/031_snp_indel_calling.sh"

done

# 032. Joint calling
bash "code/032_joint_calling.sh"

# 033. Variant filtering using GATK VQSR machine learning
snp_cutoff=90.0
indel_cutoff=90.0
env proj_n="$proj_n" snp_cutoff="$snp_cutoff" indel_cutoff="$indel_cutoff" bash "code/033_variant_filtering.sh"

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
env proj_n="$proj_n" snp_cutoff="$snp_cutoff" indel_cutoff="$indel_cutoff" filter_options_SNP="$filter_options_SNP" filter_options_INDEL="$filter_options_INDEL" bash "code/034_variant_filtering_hard.sh"

## 035. Rare/Common variants
maf=0.05

env proj_n="$proj_n" maf="$maf" bash "code/035_rare_variants.sh" 

# ## 036. Evaluate callset
# inputs=("cohort.indel.90.0.hard_filter.MAF.0.05.vcf" "cohort.indel.90.0.hard_filter.MAF.0.05.rare.vcf" "cohort.snp.90.0.hard_filter.MAF.0.05.vcf" "cohort.snp.90.0.hard_filter.MAF.0.05.rare.vcf")

# for input in "${inputs[@]}"; do
#     echo $input
#     env proj_n="$proj_n" input="$input" bash "code/036_evaluate_callset.sh"
# done

## 041. Annotation
env proj_n="$proj_n" bash "code/041_vep.sh"

## 042. Filter VEP annotation
env proj_n="$proj_n" bash "code/042_vep_filter.sh"
## 043. Convert VCF into tab_delimited
env proj_n="$proj_n" bash "code/043_vcf_to_tab.sh"

# 04. Data analysis
## 043. SVA (MAF < 0.05) (Run code/051_SVA_GBA.r)
## 044. GBA (MAF < 0.05) (Run code/051_SVA_GBA.r)





