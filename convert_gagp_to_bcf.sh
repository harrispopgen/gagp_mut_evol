#!/bin/bash

# module load htslib/1.7
# module load bcftools/1.8

# Directory should point to the vcf files from GAGP (https://eichlerlab.gs.washington.edu/greatape/data/)
bcftools view /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Pan_troglodytes.vcf.gz \
	-O b \
	-o /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Pan_troglodytes.bcf

bcftools index /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Pan_troglodytes.bcf

bcftools view /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Pan_paniscus.vcf.gz \
	-O b \
	-o /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Pan_paniscus.bcf

bcftools index /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Pan_paniscus.bcf

bcftools view /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Pongo_abelii.vcf.gz \
	-O b \
	-o /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Pongo_abelii.bcf

bcftools index /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Pongo_abelii.bcf

bcftools view /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Homo.vcf.gz \
	-O b \
	-o /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Homo.bcf

bcftools index /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Homo.bcf

bcftools view /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Gorilla.vcf.gz \
	-O b \
	-o /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Gorilla.bcf

bcftools index /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Gorilla.bcf

bcftools view /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Pongo_pygmaeus.vcf.gz \
	-O b \
	-o /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Pongo_pygmaeus.bcf

bcftools index /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Pongo_pygmaeus.bcf


