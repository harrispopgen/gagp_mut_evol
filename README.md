# gagp_mut_evol
Code for generating the main figures of Goldberg &amp; Harris 2020

## Data preprocessing
Download the Great Ape Genome Project [data](https://eichlerlab.gs.washington.edu/greatape/data/) and convert VCFs to BCF using `convert_gagp_to_bcf.sh`.

## Generating ancestral allele lookup files and species-specific references
GAGP variants for all species are called to the human reference genome (hg18 build) rather than to specific species references; coordinates are also in hg18. To streamline the process of assigning ancestral states and removing recurrent variants (variants that occur in multiple lineages), I created species-specific versions of hg18 references, in which sites with an allele ancestral to a species' lineage different from that in humans were replaced to reflect the proper ancestral state.
1. Preprocess the GAGP BCFs to generate a summary file of the variants using `preprocess_gagp_bcfs_for_ancestral_id.py`. Run for each species and chromosome e.g. `python preprocess_gagp_bcfs_for_ancestral_id.py Pan_troglodytes chr12`
2. Sort and concatenate the variant summary files using `cat_sort_files.sh`.
3. Build a table summarizing the variation at every site in all species using `generate_ancestral_allele_table.py`. This table mimics the ancestral lookup table from Prado-Martinez et al. 2013. Need to do so for each chromosome separately e.g. `python generate_ancestral_allele_table.py chr12`.
4. Generate the species-specific edited hg18 reference using `generate_hg18_reference_with_gagp_ancestral_alleles.py`. Needs to be run for each species and chromosome e.g. `python generate_hg18_reference_with_gagp_ancestral_alleles.py Pan_troglodytes chr12`.

## Generating compartment .bed files
A large amount of the data preprocessing involves generating .bed files of the genomic segments that make up each of the many compartments used in this analysis. Several compartments require others to be generated already (e.g., NCNR needs the CpG islands compartment).
* **Early and late replication timing**: need to download Supplemental Table 2 from Koren et al. 2012 paper with replication timing (z-normalized read counts @ particular time in S1) from http://mccarrolllab.org/resources/. Run `calculate_replication_timing_per_20kb.py [path to replication timing file]`; will generate .bed files with 20kb segments of the genome segregated into replication timing quartile. q0 is the latest replicating quartile, q3 is the earliest.
* **CpG islands**: need CpG island table (unmasked) from UCSC table browser hg18 build. Run `generate_cpgIslands_compartment.sh [path to cpg island table]`.
* **NCNR**: Need the following as input to `generate_hg18_ncnr.sh`
  * phastCons 44-way primate alignment from UCSC table browser
  * repeatMasker annotation of hg18 (also from UCSC)
  * Output from `refGene_exons_to_bed.py`, which requires refGene annotation of hg18 from UCSC
  * CpG islands compartment from `generate_cpgIslands_compartment.sh`
  * hg18.genome, from UCSC table browser. Just a file with the names of chrs and lengths
* **ERVs**: need repeatMasker annotation of hg18, run `generate_erv_compartment.sh [path to repeatMasker table]`
* **LINEs**: need repeatMasker annotation of hg18, run `generate_line_compartment.sh [path to repeatMasker table]`
* **Early and late replication timing, split by repetitive content**: Need repeatMasker table and the q0 and q3 replication timing compartments (above). Run `generate_rep_timing_windows_nonrepetitive.sh` and `generate_rep_timing_windows_repetitive.sh`.
* **Nonrepetitive heterochromatin**: Need the hg18 chromHMM annotations for nine cell types (from UCSC) and repeatMasker annotation as input for `chromhmm_heterochromatin.sh`. See script for cell types.
* **Maternal mutation hotspots**: Need supplemental table 12 (third supplemental file) from Jonsson et al., 2017 Nature then run `generate_maternal_hotspot_compartment.sh [path to table]`
* **ERVs, with and without hydroxymethylation**: Need ERV compartment bed file (above) and 5hmC data in H1 escs, from Yu et al., 2012. (GSM882245_H1.hmC_sites.FDR_0.0502.hg18.txt.gz) as input for `split_ervs_hydroxymethylation_compartments.py [--erv_file [path] --hmc_file [path]`

## Generating trinucleotide content
Comparing the mutation spectra between different compartments requires accounting for differences in triplet nucleotide content. I store triplet nucleotide content (literally counts of each triplet that occurs in a compartment) in a simple tab-delimited .txt file. These files need to be generated for each compartment as listed above using the `count_nmers.py` script. The naming convention for a nucleotide content file is `[compartment_name]_3mer_content.txt`. Required arguments are type of analysis (always use "bed_to_total_nmer"), the path to the compartment bed file, the length of the nmer being counted (always 3 for main figures), and the path to the outfile. An example command-line prompt is `count_nmers.py bed_to_total_nmer ./hg18_ncnr.bed 3 ./hg18_ncnr_3mer_content.txt`.

## Generating triplet mutation counts
1. Run the `process_bcf_to_mutational_signature.py` with for each species, chromosome, and compartment name. Also need to include nmer length (3) and pathway to GAGP BCF dir. Sample command line prompt: `python process_bcf_to_mutational_signature.py Pan_troglodytes chr12 hg18_ncnr 3 --gagp_bcf_dir ./data/`. Will generate the mutation counts at the individual and species level. Edit the hardcoded pathways to output files as necessary.
2. Sum the mutation counts by chromosome with `sum_mutation_counts_by_chr.py`. Provide pathway to the output files from step 1, compartment name, nmer length, and a True/False argument of if the separate chromosome output files should be deleted or not. E.g., `python sum_mutation_counts_by_chr.py ./nmer_mutation_counts/ hg18_ncnr 3 True'

## Analysis, generating figures
Code for each of the main figures has been separated into R scripts. Input data is generated above.
