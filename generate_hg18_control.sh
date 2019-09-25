#!/bin/bash

# These hg18 data tables from UCSC are 0 indexed, not 1 indexed
# awk '{if ($1 !~/^#/) print $2 "\t" $3 "\t" $4}' /net/harris/vol1/data/hg18/hg18_cpgIslandExtUnmasked.txt \
# 	> ~/temp/hg18_cpgIsland.bed

gunzip -c /net/harris/vol1/data/hg18/phastConsElements44wayPrimates.txt.gz | \
    awk '{print $2 "\t" $3 "\t" $4}' | \
    cat - /net/harris/vol1/data/hg18/hg18_repeatMasker.bed ./hg18_refGene_exons.bed ~/temp/hg18_cpgIsland.bed | \
    awk '{print $1 "\t" $2 "\t" $3}' | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i stdin | \
    bedtools complement -i stdin -g /net/harris/vol1/data/hg18/hg18.genome \
    > ./hg18_control.bed