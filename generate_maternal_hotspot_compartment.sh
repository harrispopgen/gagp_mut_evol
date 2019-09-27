#!/bin/bash

FILENAME="./data/nature24018-s3.txt"

if [ "$1" != "" ]
then
    FILENAME=$1
fi

# Need supplemental table 12 (third supplemental file) from Jonsson et al., 2017 Nature
awk '{if ($3 ~ 1) print $1 "\t" $2*1000000 "\t" ($2+1)*1000000}' $FILENAME \
    > hg38_maternal_hotspots_zero_indexed_mb.bed