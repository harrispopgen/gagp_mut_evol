#!/bin/bash

# For chromHMM tracks: download the .txt.gz files from http://hgdownload.soe.ucsc.edu/goldenPath/hg18/database/
# Then convert to bed file (just slice the file for the chrom, start, & end positions
# I kept them zipped bc they were big

bedtools intersect -a <(gunzip -c ./data/chromHMM_tracks/wgEncodeBroadHmmGm12878HMM.bed.gz | awk '{if ($4 == 13) print $0}') \
                   -b <(gunzip -c ./data/chromHMM_tracks/wgEncodeBroadHmmH1hescHMM.bed.gz | awk '{if ($4 == 13) print $0}') |
    bedtools intersect -a stdin \
                       -b <(gunzip -c ./data/chromHMM_tracks/wgEncodeBroadHmmHepg2HMM.bed.gz | awk '{if ($4 == 13) print $0}') |
    bedtools intersect -a stdin \
                       -b <(gunzip -c ./data/chromHMM_tracks/wgEncodeBroadHmmHmecHMM.bed.gz | awk '{if ($4 == 13) print $0}') |
    bedtools intersect -a stdin \
                       -b <(gunzip -c ./data/chromHMM_tracks/wgEncodeBroadHmmHsmmHMM.bed.gz | awk '{if ($4 == 13) print $0}') |
    bedtools intersect -a stdin \
                       -b <(gunzip -c ./data/chromHMM_tracks/wgEncodeBroadHmmHuvecHMM.bed.gz | awk '{if ($4 == 13) print $0}') |
    bedtools intersect -a stdin \
                       -b <(gunzip -c ./data/chromHMM_tracks/wgEncodeBroadHmmK562HMM.bed.gz | awk '{if ($4 == 13) print $0}') |
    bedtools intersect -a stdin \
                       -b <(gunzip -c ./data/chromHMM_tracks/wgEncodeBroadHmmNhekHMM.bed.gz | awk '{if ($4 == 13) print $0}') |
    bedtools intersect -a stdin \
                       -b <(gunzip -c ./data/chromHMM_tracks/wgEncodeBroadHmmNhlfHMM.bed.gz | awk '{if ($4 == 13) print $0}') \
                       > ./data/chromHMM_tracks/chromHMM_heterochromatin.bed

# Need a bed file of the hg18 repeatMasker output (also from UCSC table browser, hg18 build)
# Again, slice .txt file for chrom, start, & end positions and convert to bed file
bedtools subtract -a ./data/chromHMM_heterochromatin.bed \
                  -b ./data/hg18_repeatMasker.bed \
                  > ./data/hg18_nonrepetitive_chromHMM_heterochromatin.bed