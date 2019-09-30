# gagp_mut_evol
Code for generating the main figures of Goldberg &amp; Harris 2019 preprint

## Data preprocessing
Download the Great Ape Genome Project [data](https://eichlerlab.gs.washington.edu/greatape/data/) and convert VCFs to BCF using `convert_gagp_to_bcf.sh`.

## Generating ancestral allele lookup files and species-specific references
GAGP variants for all species are called to the human reference genome (hg18 build) rather than to specific species references; coordinates are also in hg18. To streamline the process of assigning ancestral states and removing recurrent variants (variants that occur in multiple lineages), I created species-specific versions of hg18 references, in which sites with an allele ancestral to a species' lineage different from that in humans were replaced to reflect the proper ancestral state.
1. Preprocess the GAGP BCFs to generate a summary file of the variants using `preprocess_gagp_bcfs_for_ancestral_id.py`. Run for each species and chromosome e.g. `python preprocess_gagp_bcfs_for_ancestral_id.py Pan_troglodytes chr12`
2. Sort and concatenate the variant summary files using `cat_sort_files.sh`.
3. Build a table summarizing the variation at every site in all species using `generate_ancestral_allele_table.py`. This table mimics the ancestral lookup table from Prado-Martinez et al. 2013. Need to do so for each chromosome separately e.g. `python generate_ancestral_allele_table.py chr12`.
4. Generate the species-specific edited hg18 reference using `generate_hg18_reference_with_gagp_ancestral_alleles.py`. Needs to be run for each species and chromosome e.g. `python generate_hg18_reference_with_gagp_ancestral_alleles.py Pan_troglodytes chr12`.

## Generating compartment .bed files
A large amount of the data preprocessing involves generating .bed files of the genomic segments that make up each of the many compartments used in this analysis.
* NCNR
* ERVs, with and without hydroxymethylation
* ERVs
* LINEs
* Maternal mutation hotspots
* CpG islands
* Early and late replication timing
* Early and late replication timing, split by repetitive content
* Nonrepetitive heterochromatin

## Generating trinucleotide content

## Generating triplet mutation counts

## Analysis, generating figures
