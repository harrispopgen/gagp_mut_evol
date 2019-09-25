#!/bin/bash

bedtools intersect -a <(gunzip -c ./chromHMM_tracks/wgEncodeBroadHmmGm12878HMM.bed.gz | awk '{if ($4 == 13) print $0}') \
                   -b <(gunzip -c ./chromHMM_tracks/wgEncodeBroadHmmH1hescHMM.bed.gz | awk '{if ($4 == 13) print $0}') |
    bedtools intersect -a stdin \
                       -b <(gunzip -c ./chromHMM_tracks/wgEncodeBroadHmmHepg2HMM.bed.gz | awk '{if ($4 == 13) print $0}') |
    bedtools intersect -a stdin \
                       -b <(gunzip -c ./chromHMM_tracks/wgEncodeBroadHmmHmecHMM.bed.gz | awk '{if ($4 == 13) print $0}') |
    bedtools intersect -a stdin \
                       -b <(gunzip -c ./chromHMM_tracks/wgEncodeBroadHmmHsmmHMM.bed.gz | awk '{if ($4 == 13) print $0}') |
    bedtools intersect -a stdin \
                       -b <(gunzip -c ./chromHMM_tracks/wgEncodeBroadHmmHuvecHMM.bed.gz | awk '{if ($4 == 13) print $0}') |
    bedtools intersect -a stdin \
                       -b <(gunzip -c ./chromHMM_tracks/wgEncodeBroadHmmK562HMM.bed.gz | awk '{if ($4 == 13) print $0}') |
    bedtools intersect -a stdin \
                       -b <(gunzip -c ./chromHMM_tracks/wgEncodeBroadHmmNhekHMM.bed.gz | awk '{if ($4 == 13) print $0}') |
    bedtools intersect -a stdin \
                       -b <(gunzip -c ./chromHMM_tracks/wgEncodeBroadHmmNhlfHMM.bed.gz | awk '{if ($4 == 13) print $0}') \
                       > ./chromHMM_tracks/chromHMM_heterochromatin.bed

bedtools subtract -a ./chromHMM_tracks/chromHMM_heterochromatin.bed \
                  -b /net/harris/vol1/data/hg18/hg18_repeatMasker.bed \
                  > ./chromHMM_tracks/hg18_nonrepetitive_chromHMM_heterochromatin.bed