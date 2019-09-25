#!/bin/bash

awk '{if ($3 ~ 1) print $1 "\t" ($2-1)*1000000 "\t" $2*1000000}' nature24018-s3.txt \
    > hg38_maternal_hotspots_one_indexed_mb.bed

awk '{if ($3 ~ 1) print $1 "\t" $2*1000000 "\t" ($2+1)*1000000}' nature24018-s3.txt \
    > hg38_maternal_hotspots_zero_indexed_mb.bed