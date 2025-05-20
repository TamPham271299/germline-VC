#!bin/bash

vep_o="/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/$proj_n/vep"

# # SNV
# input="$vep_o/cohort.snp.90.0.hard_filter.MAF.0.05.rare.vep.vcf.gz"
# missense="Consequence is missense_variant and gnomADe_FILTER is PASS and gnomADe_AF and gnomADe_AF_eas"
# LoF="(Consequence is stop_gained or Consequence is splice_acceptor_variant or Consequence is splice_donor_variant or Consequence is start_lost or Consequence is stop_lost) and gnomADe_FILTER is PASS and gnomADe_AF and gnomADe_AF_eas"
# syno="Consequence is synonymous_variant and gnomADe_FILTER is PASS and gnomADe_AF and gnomADe_AF_eas"
# utr="(Consequence is 5_prime_UTR_variant or Consequence is 3_prime_UTR_variant) and gnomADe_FILTER is PASS and gnomADe_AF and gnomADe_AF_eas"
# all="gnomADe_FILTER is PASS and gnomADe_AF and gnomADe_AF_eas"

# for filter_exp in missense LoF syno utr all; do
#     echo $filter_exp
#     echo ${!filter_exp}
#     filter_vep --gz --format vcf --input_file $input \
#                 --filter "${!filter_exp}" \
#                 --output_file ${input/.vep.vcf.gz}.vep_filt.${filter_exp}.vcf --force_overwrite \
#                 --only_matched
# done

# # indels
# input="$vep_o/cohort.indel.90.0.hard_filter.MAF.0.05.rare.vep.vcf.gz"
# frameshift="Consequence is frameshift_variant and gnomADe_FILTER is PASS and gnomADe_AF and gnomADe_AF_eas"

# for filter_exp in frameshift; do
#     echo $filter_exp
#     echo ${!filter_exp}
#     filter_vep --gz --format vcf --input_file $input \
#                 --filter "${!filter_exp}" \
#                 --output_file ${input/.vep.vcf.gz}.vep_filt.${filter_exp}.vcf --force_overwrite \
#                 --only_matched
# done

## Re-filter with match string not match exactly
# SNV
input="$vep_o/cohort.snp.90.0.hard_filter.MAF.0.05.rare.vep.vcf.gz"
missense="Consequence match missense_variant and gnomADe_FILTER is PASS and gnomADe_AF and gnomADe_AF_eas"
LoF="(Consequence match stop_gained or Consequence match splice_acceptor_variant or Consequence match splice_donor_variant or Consequence match start_lost or Consequence match stop_lost) and gnomADe_FILTER is PASS and gnomADe_AF and gnomADe_AF_eas"

for filter_exp in missense LoF; do
    echo $filter_exp
    echo ${!filter_exp}
    filter_vep --gz --format vcf --input_file $input \
                --filter "${!filter_exp}" \
                --output_file ${input/.vep.vcf.gz}.vep_filt.${filter_exp}.v2.vcf --force_overwrite \
                --only_matched
done

# indels
input="$vep_o/cohort.indel.90.0.hard_filter.MAF.0.05.rare.vep.vcf.gz"
frameshift="Consequence match frameshift_variant and gnomADe_FILTER is PASS and gnomADe_AF and gnomADe_AF_eas"

for filter_exp in frameshift; do
    echo $filter_exp
    echo ${!filter_exp}
    filter_vep --gz --format vcf --input_file $input \
                --filter "${!filter_exp}" \
                --output_file ${input/.vep.vcf.gz}.vep_filt.${filter_exp}.v2.vcf --force_overwrite \
                --only_matched
done
