#!/bin/bash

REPEATMASKER_FILENAME="/net/harris/vol1/data/hg18/hg18_repeatMasker.bed"

if [ "$1" != "" ]
then
    REPEATMASKER_FILENAME=$1
fi


bedtools subtract -a ./hg18_replication_timing_q0.bed \
     -b $REPEATMASKER_FILENAME \
     > ./hg18_replication_timing_q0_nonrepetitive.bed

bedtools subtract -a ./hg18_replication_timing_q3.bed \
     -b $REPEATMASKER_FILENAME \
     > ./hg18_replication_timing_q3_nonrepetitive.bed
