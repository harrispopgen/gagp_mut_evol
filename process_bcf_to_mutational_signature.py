#!/usr/bin/env python

# Code to convince us that assigning each copy of a mutation to all haplotypes that share it isn't a terrible idea

import csv
import glob
import re
import gzip
import pdb
import argparse
import os
import pysam
import math
import itertools
import random
import numpy as np
from common import hwe, filter_var, rev_comp, generate_nmer_mutation_list, load_hg18_ref, complement

# random.seed(0)

def process_gagp_bcf_to_nmers_randomize_mut_assignment(species, chr, bed_file, ref_seq, nmer, nmer_mutations, gagp_bcf_dir = "/net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/"):
	
	upstream_n = int(nmer/2)
	downstream_n = int(nmer/2)
		
	indiv_dict = {}
	mut_dict = dict(zip(nmer_mutations, [0 for _ in nmer_mutations]))
	
	gagp_bcf_filename = gagp_bcf_dir + species + '.bcf'
	
	with pysam.VariantFile(gagp_bcf_filename) as bcf_in:
		
		indivs = list(bcf_in.header.samples)
		indiv_dict = dict(zip(indivs, [dict(zip(nmer_mutations, [0 for _ in nmer_mutations])) for _ in indivs]))
		
		with open(bed_file, 'rU') as open_bed:
			
			line = ''
			while not line.startswith(chr + '\t'):
				line = open_bed.readline() # Incrementing bed file to correct chromosome
				
			while line.startswith(chr + '\t'):
				
				(chr, start, stop) = line.split()[:3]
				vars = bcf_in.fetch(chr, int(start), int(stop))
				for v in vars:
					
					if not filter_var(v):
						continue
					
					curr_nmer = hg18_ref[v.pos-(upstream_n+1):v.pos+downstream_n].upper()
					
					if 'N' in curr_nmer or 'X' in curr_nmer: # Filtering N, X
						continue
					if v.ref == curr_nmer[upstream_n]:
						ancestral_allele = v.ref
						derived_allele = v.alts[0]
						derived_gt = 1
					elif v.alts[0] == curr_nmer[upstream_n]:
						ancestral_allele = v.alts[0]
						derived_allele = v.ref
						derived_gt = 0
					else: # If neither of the alleles match the hg18 reference ASK ABOUT THIS
						continue
					
					if ancestral_allele in ['G', 'T']:
						nmer_w_mut = rev_comp(curr_nmer) + '>' + complement[derived_allele]
					else:
						nmer_w_mut = curr_nmer + '>' + derived_allele
					
					# Excluding MAF >= 0.5
					if (derived_gt == 1 and v.info['AF'][0] >= 0.5) or (derived_gt == 0 and v.info['AF'][0] <= 0.5):
						continue
					
					# Here is where we randomly assign a mutation to a single haplotype
					# List of samples with 1 or 2 copies of derived allele
					samples_with_mut = [i for i in indivs if v.samples[i]['GT'].count(derived_gt) > 0]
					# Randomly give only one of the individuals the mutation (only a single time)
					indiv_dict[random.choice(samples_with_mut)][nmer_w_mut] += 1
					
					mut_dict[nmer_w_mut] += 1
					
				line = open_bed.readline()
# 	pdb.set_trace()
	return mut_dict, indiv_dict, indivs
		
def write_species_mut_dict_to_file(species, chr, snp_type, nmer, mut_dict, nmer_mutations):
	
	outfilename = './nmer_mutation_counts/' + species + '_' + chr + '_' + snp_type + '_species_' + str(nmer) + 'mer_counts.txt'
	with open(outfilename, 'w') as outfile:
		outfile.write('\t'.join(sorted(nmer_mutations)))
		outfile.write('\n')
		outfile.write('\t'.join([str(mut_dict[m]) for m in sorted(nmer_mutations)]))
		outfile.write('\n')
	
	return None

def write_indiv_mut_dict_to_file(species, chr, snp_type, nmer, indiv_dict, indivs, nmer_mutations):
	
	outfilename = './nmer_mutation_counts/' + species + '_' + chr + '_' + snp_type + '_indiv_' + str(nmer) + 'mer_counts.txt'
	with open(outfilename, 'w') as outfile:
		outfile.write('indiv\t')
		outfile.write('\t'.join(sorted(nmer_mutations)))
		outfile.write('\n')
		for i in indivs:
			outfile.write(i + '\t')
			outfile.write('\t'.join([str(indiv_dict[i][m]) for m in sorted(nmer_mutations)]))
			outfile.write('\n')
	
	return None

def argument_parse():
	
	parser = argparse.ArgumentParser()
	parser.add_argument("species", help="Name of species, as defined by GAGP")
	parser.add_argument("chr", help="Chromosome (e.g. 'chr10', 'chrY')")
	parser.add_argument("bed_file", help="Path to a bed file containing chromosomal regions (0-index) to analyze")
	parser.add_argument("nmer", help="Mutational nmer (must be odd number, mutation in the middle)", type=int)
	parser.add_argument("--gagp_bcf_dir", help="Directory holding the GAGP BCFs")
	args = parser.parse_args()
	return args

if __name__ == '__main__':
	
	args = argument_parse()
	hg18_ref = load_hg18_ref(args.species, args.chr, exclude_recurrent = True)
	if args.nmer % 2 != 1:
		raise ValueError("nmer must be an odd number")
	nmer_mutations = generate_nmer_mutation_list(args.nmer)
	if args.gagp_bcf_dir:
		mut_dict, indiv_dict, indivs = process_gagp_bcf_to_nmers_randomize_mut_assignment(args.species, args.chr, args.bed_file, hg18_ref, args.nmer, nmer_mutations, args.gagp_bcf_dir)
	else:
		mut_dict, indiv_dict, indivs = process_gagp_bcf_to_nmers_randomize_mut_assignment(args.species, args.chr, args.bed_file, hg18_ref, args.nmer, nmer_mutations)
	snp_type = args.bed_file.split('/')[-1][:-4]
	write_indiv_mut_dict_to_file(args.species, args.chr, snp_type, args.nmer, indiv_dict, indivs, nmer_mutations)
	write_species_mut_dict_to_file(args.species, args.chr, snp_type, args.nmer, mut_dict, nmer_mutations)
	
