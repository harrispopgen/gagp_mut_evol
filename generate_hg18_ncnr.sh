#!/bin/bash

# phastCons 44-way primate alignment from UCSC table browser
PHASTCONS_FILENAME="/net/harris/vol1/data/hg18/phastConsElements44wayPrimates.txt.gz"
# repeatMasker annotation of hg18
REPEATMASKER_FILENAME="/net/harris/vol1/data/hg18/hg18_repeatMasker.bed"
# Generated from refGene_exons_to_bed.py
REFGENE_EXONS_FILENAME="./hg18_refGene_exons.bed"
# From generate_cpgIslands_compartment.sh
CPG_ISLAND_FILENAME="~/temp/hg18_cpgIsland.bed"
# From UCSC table browser. Just a file with the names of chrs and lengths
GENOME_FILENAME="/net/harris/vol1/data/hg18/hg18.genome"

while [ "$1" != "" ]; do
    case $1 in
        --phastcons_filename)       shift
                               PHASTCONS_FILENAME=$1
                               ;;
        --repeatmasker_filename)    shift
                               REPEATMASKER_FILENAME=$1
                               ;;
        --refgene_exons_filename)    shift
                               REFGENE_EXONS_FILENAME=$1
                               ;;
        --cpg_island_filename)  shift
                               CPG_ISLAND_FILENAME=$1
                               ;;
        --genome_filename)      shift
                               GENOME_FILENAME=$1
                               ;;
    esac
    shift
done

gunzip -c $PHASTCONS_FILENAME | \
    awk '{print $2 "\t" $3 "\t" $4}' | \
    cat - $REPEATMASKER_FILENAME $REFGENE_EXONS_FILENAME $CPG_ISLAND_FILENAME | \
    awk '{print $1 "\t" $2 "\t" $3}' | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i stdin | \
    bedtools complement -i stdin -g $GENOME_FILENAME \
    > ./hg18_ncnr.bed