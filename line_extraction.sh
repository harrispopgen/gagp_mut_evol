#!/bin/bash

awk '{if ($12 ~ "LINE") print $0}' /net/harris/vol1/data/hg18/hg18_repeatMasker.txt > \
     ./line_data/hg18_all_line.txt

awk '{if ($12 ~ "LINE") print $6 "\t" $7 "\t" $8 "\t" NR-1}' /net/harris/vol1/data/hg18/hg18_repeatMasker.txt > \
     ./line_data/hg18_all_line.bed
