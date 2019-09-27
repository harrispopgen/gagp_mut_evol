#!/bin/bash

# concatenate then sort variant summary files for aa lookups 
# (output from preprocess preprocess_gagp_bcfs_for_ancestral_id, input to generate_ancestral_allele_table)

for i in {1..22}
do
   cat ./preprocessed_gagp_bcfs_for_ancestral_identification/*_"chr$i"_var_summaries.txt | \
      sort -k 2,2 -n > ./preprocessed_gagp_bcfs_for_ancestral_identification/chr"$i".txt
done