#!/usr/bin/env python

# Generate a summary file of the variants in a single species & chromosome to build references w/ ancestral alleles

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
import numpy as np
from common import hwe, filter_var, rev_comp, generate_nmer_mutation_list, load_hg18_ref, complement

def preprocess_bcf(species, chr):
	
	bed_file = "./data/hg18.bed"
	
	# point to bcf
	gagp_bcf_filename = '/net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/' + species + '.bcf'
	
	# Need to make this directory to store preprocessing files
	with open('./preprocessed_gagp_bcfs_for_ancestral_identification/' + species + '_' + chr + '_var_summaries.txt', 'w') as open_output_file:
	
		with pysam.VariantFile(gagp_bcf_filename) as bcf_in:
		
			with open(bed_file, 'rU') as open_bed: # Slice by chromosome
			
				curr_output = []
				line = ''
				while not line.startswith(chr + '\t'):
					line = open_bed.readline() # Incrementing bed file to correct chromosome
				
				while line.startswith(chr + '\t'):
				
					(chr, start, stop) = line.split()[:3]
					vars = bcf_in.fetch(chr, int(start), int(stop))
					
					
					for v in vars:
						
						var_state = None
						if len(v.alts) > 1:
							var_state = 'multiallelic'
						elif v.info['AF'][0] == 0.0:
							var_state = 'fixed_ref'
						elif v.info['AF'][0] == 1.0:
							var_state = 'fixed_alt'
						else:
							var_state = 'biallelic'
						curr_output.append('%s\t%i\t%s\t%s\t%s\t%s' % (chr, v.pos, v.ref, '|'.join(v.alts), var_state, species))
						if len(curr_output) > 5000:
							open_output_file.write('\n'.join(curr_output) + '\n')
							curr_output = []
						
					open_output_file.write('\n'.join(curr_output) + '\n')
					curr_output = []
					
					line = open_bed.readline()
		
	return None

def argument_parse():
	
	parser = argparse.ArgumentParser()
	parser.add_argument("species", help="Name of species, as defined by GAGP (e.g. Homo, Pan_troglodytes)")
	parser.add_argument("chr", help="Chromosome (e.g. 'chr10')")
	args = parser.parse_args()
	return args

if __name__ == '__main__':
	
	args = argument_parse()
	preprocess_bcf(args.species, args.chr)
	


