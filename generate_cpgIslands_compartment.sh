#!/bin/bash

# These hg18 data tables from UCSC are 0 indexed, not 1 indexed
# Need CpG island table (unmasked) from UCSC table browser hg18 build

FILENAME="/net/harris/vol1/data/hg18/hg18_cpgIslandExtUnmasked.txt"

if [ "$1" != "" ]
then
    FILENAME=$1
fi

awk '{if ($1 !~/^#/) print $2 "\t" $3 "\t" $4}' $FILENAME \
	> ./hg18_cpgIsland.bed
