#!bin/bash

FASTA="/home/mdnl/.vep/homo_sapiens/113_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
ref="/home/mdnl/Documents/Tam_14Mar25/ref"
gatk_o="/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/$proj_n/gatk"
vep_o="/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/$proj_n/vep"
CADD="/home/mdnl/KianaRDS/Tam_14Mar25/ref/vep_plugins/CADD"

mkdir -p $vep_o/tmp

VCFs=("cohort.snp.90.0.hard_filter.MAF.0.05.rare.vcf" "cohort.indel.90.0.hard_filter.MAF.0.05.rare.vcf")

for VCF in "${VCFs[@]}"; do
    echo $VCF
    input="$gatk_o/$VCF"

    for chr in chr{1..22} chrX chrY; do
        echo $chr
        output="$vep_o/tmp/${VCF/.vcf}.${chr}.vep.vcf.gz"
        gnomAD="$ref/gnomADv4_1/exome/hg38/gnomad.exomes.v4.1.sites.${chr}_trimmed.vcf.bgz"

        vep --format vcf --input_file $input --no_stats \
            --vcf --compress_output gzip --output_file $output --force_overwrite \
            --cache --fasta $FASTA --buffer_size 10000 \
            --variant_class --sift b --polyphen b --gene_phenotype --numbers \
            --hgvs --protein --symbol --ccds --uniprot --biotype --domains --canonical --transcript_version --gene_version \
            --check_existing --exclude_null_alleles --var_synonyms --pubmed --failed 1 \
            --per_gene --chr $chr \
            --custom file=$gnomAD,short_name=gnomADe,format=vcf,type=exact,coords=0,fields=FILTER%AF%AC%AN%AF_eas%AC_eas%AN_eas \
            --plugin CADD,snv=$CADD/whole_genome_SNVs.tsv.gz,indels=$CADD/gnomad.genomes.r4.0.indel.tsv.gz \
            --plugin GeneBe
    done

    # Concatenate VEP output
    # conda activate variant_calling
    bcftools concat -O z --threads 8 \
                    --write-index \
                    -o "$vep_o/${VCF/.vcf}.vep.vcf.gz" \
                    $vep_o/tmp/${VCF/.vcf}.chr{1..22}.vep.CADD.vcf.gz $vep_o/tmp/${VCF/.vcf}.chrX.vep.CADD.vcf.gz $vep_o/tmp/${VCF/.vcf}.chrY.vep.CADD.vcf.gz
    # conda deactivate
done

