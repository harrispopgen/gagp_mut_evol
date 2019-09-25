#!/bin/bash

for i in {1..22}
do
   cat ./*_"chr$i"_var_summaries.txt | \
      sort -k 2,2 -n > ./chr"$i".txt
done