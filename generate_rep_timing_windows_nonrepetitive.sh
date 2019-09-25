#!/bin/bash

bedtools subtract -a ./replication_timing_compartments/hg18_replication_timing_q0.bed \
     -b /net/harris/vol1/data/hg18/hg18_repeatMasker.bed \
     > ./replication_timing_compartments/hg18_replication_timing_q0_nonrepetitive.bed

bedtools subtract -a ./replication_timing_compartments/hg18_replication_timing_q3.bed \
     -b /net/harris/vol1/data/hg18/hg18_repeatMasker.bed \
     > ./replication_timing_compartments/hg18_replication_timing_q3_nonrepetitive.bed