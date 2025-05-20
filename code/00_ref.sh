#!bin/bash

ref="/home/mdnl/Documents/Tam_14Mar25/ref"
mkdir -p $ref

# Known datasets: GATK bundle for human hg38 reference
# Move GATK resources to "/home/mdnl/KianaRDS/Tam_14Mar25/ref"
ref="/home/mdnl/KianaRDS/Tam_14Mar25/ref"
mkdir -p $ref/Gencode/GATK_resource_bundle/hg38
cd $ref/Gencode/GATK_resource_bundle/hg38
fl=(
    "1000G_omni2.5.hg38.vcf.gz" "1000G_omni2.5.hg38.vcf.gz.tbi" 
    "1000G_phase1.snps.high_confidence.hg38.vcf.gz" "1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi" 
    "Homo_sapiens_assembly38.dbsnp138.vcf" "Homo_sapiens_assembly38.dbsnp138.vcf.idx" 
    "hapmap_3.3.hg38.vcf.gz" "hapmap_3.3.hg38.vcf.gz.tbi"
    "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
    )

for f in "${fl[@]}"; do
    echo $f
    wget -c "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/$f"
done

# human hg38 reference FASTA for bwa-mem2 index & GATK
mkdir -p $ref/Gencode/bwa-mem2-idx/hg38
cd $ref/Gencode/bwa-mem2-idx/hg38
wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
# samtools faidx Homo_sapiens_assembly38.fasta
wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict
# gatk-launch CreateSequenceDictionary -R Homo_sapiens_assembly38.fasta

# index genome
export PATH=/home/mdnl/Documents/Tam_14Mar25/tools/bwa-mem2:$PATH
FASTA="$ref/Gencode/bwa-mem2-idx/hg38/Homo_sapiens_assembly38.fasta"
bwa-mem2 index $FASTA
# [bwa_index] Pack FASTA... 8.72 sec
# * Entering FMI_search
# init ticks = 242505940968
# ref seq len = 6434693834
# binary seq ticks = 118351257300
# build suffix-array ticks = 2659445000712
# ref_seq_len = 6434693834
# count = 0, 1882204623, 3217346917, 4552489211, 6434693834
# BWT[2736215060] = 4
# CP_SHIFT = 6, CP_MASK = 63
# sizeof CP_OCC = 64
# pos: 804336730, ref_seq_len__: 804336729
# max_occ_ind = 100542091
# build fm-index ticks = 567971541648
# Total time taken: 1011.4222

# Download gnomAD custome databases to query allele #, AF...
## Exome
mkdir -p $ref/gnomADv4_1/exome/hg38
cd $ref/gnomADv4_1/exome/hg38

for chr in chr{1..22} chrX chrY; do
    echo $chr
    wget -c "https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad/v4.1/exomes/gnomad.exomes.v4.1.sites.${chr}_trimmed.vcf.bgz"
    wget -c "https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad/v4.1/exomes/gnomad.exomes.v4.1.sites.${chr}_trimmed.vcf.bgz.tbi"
done

## Genome
ref="~/KianaRDS/Tam_14Mar25/..."
mkdir -p $ref/gnomADv4_1/genome/hg38
cd $ref/gnomADv4_1/genome/hg38

for chr in chr{1..22} chrX chrY; do
    echo $chr
    wget -c "https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad/v4.1/genomes/gnomad.genomes.v4.1.sites.${chr}_trimmed.vcf.bgz"
    wget -c "https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad/v4.1/genomes/gnomad.genomes.v4.1.sites.${chr}_trimmed.vcf.bgz.tbi"
done

# Download CADD score
ref="/home/mdnl/KianaRDS/Tam_14Mar25/ref"
mkdir -p $ref/vep_plugins/CADD
cd $ref/vep_plugins/CADD
wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz
wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz.tbi
wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz
wget -c https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz.tbi

# Download spliceAI data

# Download SpliceVault data
ref="/home/mdnl/KianaRDS/Tam_14Mar25/ref"
mkdir -p $ref/vep_plugins/SpliceVault
cd $ref/vep_plugins/SpliceVault
wget -c https://ftp.ensembl.org/pub/current_variation/SpliceVault/SpliceVault_data_GRCh38.tsv.gz
wget -c https://ftp.ensembl.org/pub/current_variation/SpliceVault/SpliceVault_data_GRCh38.tsv.gz.tbi
