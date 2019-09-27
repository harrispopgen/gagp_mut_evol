#!/bin/bash

REPEATMASKER_FILENAME="/net/harris/vol1/data/hg18/hg18_repeatMasker.txt"

if [ "$1" != "" ]
then
    REPEATMASKER_FILENAME=$1
fi


# Replace with pw to output from repeatMasker from hg18 UCSC table browser
awk '{if ($12 ~ "LTR" && ($13 ~ "ERV1" || $13 ~ "ERVK" || $13 ~ "ERVL" || $13 ~ "Gypsy")) print $6 "\t" $7 "\t" $8 "\t" NR-1}' $REPEATMASKER_FILENAME > \
     ./hg18_all_erv.bed
