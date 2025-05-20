#!bin/bash

vep_o="/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/$proj_n/vep"

# Extract sample names
samples=$(zcat $vep_o/cohort.snp.90.0.hard_filter.MAF.0.05.rare.vep.vcf.gz| grep "^#CHROM" | cut -f10-)

# SNV
for cons_type in LoF missense; do
    echo $cons_type
    VCF="$vep_o/cohort.snp.90.0.hard_filter.MAF.0.05.rare.vep_filt.${cons_type}.v2.vcf"
    # bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%AC\t%AN\t%CSQ\n' \
    #     -d -A tab -HH \
    #     $VCF| sed 's/#CHROM/CHROM/g' > ${VCF/.vcf/.tsv}

    bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%AC\t%AN\t%CSQ[\t%GT]\n' \
        -d -A tab -HH \
        $VCF| \
        sed -e 's/#CHROM/CHROM/g' -e "s/GT\tGT\tGT\tGT\tGT\tGT\tGT\tGT\tGT\tGT\tGT\tGT\tGT\tGT\tGT\tGT/$(echo -e "$samples")/g" > ${VCF/.vcf/.tsv}
done

# indels
for cons_type in frameshift; do
    echo $cons_type
    VCF="$vep_o/cohort.indel.90.0.hard_filter.MAF.0.05.rare.vep_filt.${cons_type}.v2.vcf"
    # bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%AC\t%AN\t%CSQ\n' \
    #                 -d -A tab -HH \
    #                 $VCF| sed 's/#CHROM/CHROM/g' > ${VCF/.vcf/.tsv}

    bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%AC\t%AN\t%CSQ[\t%GT]\n' \
        -d -A tab -HH \
        $VCF| \
        sed -e 's/#CHROM/CHROM/g' -e "s/GT\tGT\tGT\tGT\tGT\tGT\tGT\tGT\tGT\tGT\tGT\tGT\tGT\tGT\tGT\tGT/$(echo -e "$samples")/g" > ${VCF/.vcf/.tsv}
done

# SNV
VCF="$vep_o/cohort.snp.90.0.hard_filter.MAF.0.05.rare.vep.vcf.gz"
bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%AC\t%AN\t%CSQ\n' \
    -d -A tab -HH \
    $VCF| \
    sed -e 's/#CHROM/CHROM/g' > ${VCF/.vcf.gz/.tsv}

# indel
VCF="$vep_o/cohort.indel.90.0.hard_filter.MAF.0.05.rare.vep.vcf.gz"
bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%AC\t%AN\t%CSQ\n' \
    -d -A tab -HH \
    $VCF| \
    sed -e 's/#CHROM/CHROM/g' > ${VCF/.vcf.gz/.tsv}