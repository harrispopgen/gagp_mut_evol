#!/usr/bin/env python

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

def process_curr_alts(curr_alts_dict):
	
	curr_line = []
	for genus in ['Homo', 'Pan', 'Gorilla', 'Pongo']:
		if curr_alts_dict[genus]:
# 			pdb.set_trace()
			if len(curr_alts_dict[genus]) == 1:
				curr_alt = curr_alts_dict[genus][0]
				if curr_alt[4] == 'fixed_ref':
					curr_line.append(curr_alt[2])
				elif curr_alt[4] == 'fixed_alt':
					curr_line.append(curr_alt[3])
				elif curr_alt[4] == 'biallelic':
					curr_line.append('/'.join(curr_alt[2:4]))
				elif curr_alt[4] == 'multiallelic':
					pdb.set_trace()
					curr_line.append('N')
				else: # Something else is wrong
					pdb.set_trace()
					curr_line.append('N')
			elif len(curr_alts_dict[genus]) == 2: # two species
				curr_alts = curr_alts_dict[genus]
				curr_calls = [x[4] for x in curr_alts]
# 				pdb.set_trace()
				if 'multiallelic' in curr_calls:
					pdb.set_trace()
					curr_line.append('N')
				elif 'fixed_alt' in curr_calls:
					if curr_calls[0] == curr_calls[1]: # both fixed_alt
						if curr_alts[0][3] != curr_alts[1][3]: # Different alts
							curr_line.append('N')
						else:
							curr_line.append(curr_alts[0][3]) # Same alt
					elif 'biallelic' in curr_calls:
						if curr_alts[0][3] != curr_alts[1][3]:
							curr_line.append('N') # different alt
						else: # Same alt, we'll call it segregating in the genus
							curr_line.append('/'.join(curr_alts[0][2:4]))
					elif curr_calls[0] == 'fixed_alt': # Because one alt is fixed and other is ref, it's segregating in the genus
						curr_line.append('/'.join(curr_alts[0][2:4]))
					elif curr_calls[1] == 'fixed_alt':
						curr_line.append('/'.join(curr_alts[1][2:4]))
					else:
						pdb.set_trace()
				elif 'biallelic' in curr_calls:
					if curr_calls[0] == curr_calls[1]: # both biallelic
						if curr_alts[0][3] != curr_alts[1][3]: # Different alts
							curr_line.append('N')
						else:
							curr_line.append('/'.join(curr_alts[0][2:4])) # Same alt
					elif curr_calls[0] == 'biallelic':
						curr_line.append('/'.join(curr_alts[0][2:4]))
					elif curr_calls[1] == 'biallelic':
						curr_line.append('/'.join(curr_alts[1][2:4]))
					else: # error
						pdb.set_trace()
				elif curr_calls[0] == curr_calls[1] and curr_calls[0] == 'fixed_ref':
					curr_line.append(curr_alts[0][2])
				else:
					pdb.set_trace()
					curr_line.append('N')
			else: # problem, should only have two species max
				pdb.set_trace()
		else:
			curr_line.append('-')
	
	return curr_line

def generate_aa_table(chr):
	
	# human	bonobo/chimp	gorilla	orang/orang
	
	with open('./ancestral_allele_tables/' + chr + '.txt', 'w') as open_output_file:
		output_lines = [['chromosome', 'position', 'homo', 'pan', 'gorilla', 'pongo']]
		with open('./preprocessed_gagp_bcfs_for_ancestral_identification/' + chr + '.txt', 'rU') as open_preprocessed_file:
			curr_pos = None
			curr_alts = {'Homo': [], 'Pan': [], 'Gorilla': [], 'Pongo': []}
			counter = 0
			for l in open_preprocessed_file:
				l_split = l.rstrip('\n').split('\t')
				if curr_pos != l_split[1] and curr_pos is not None:
					counter += 1
					if counter == 5000:
						open_output_file.write('\n'.join(['\t'.join(x) for x in output_lines]))
						open_output_file.write('\n')
						output_lines = []
						counter = 0
					# process lines
					output_lines.append([chr, curr_pos] + process_curr_alts(curr_alts))
					# iterate to next position
					curr_alts = {'Homo': [], 'Pan': [], 'Gorilla': [], 'Pongo': []}
				curr_pos = l_split[1]
				curr_alts[l_split[5].split('_')[0]].append(l_split)
			output_lines.append([chr, curr_pos] + process_curr_alts(curr_alts))
		open_output_file.write('\n'.join(['\t'.join(x).upper() for x in output_lines]))
	
	return None

def argument_parse():
	
	parser = argparse.ArgumentParser()
	parser.add_argument("chr", help="Chromosome (e.g. 'chr10')")
	args = parser.parse_args()
	return args

if __name__ == '__main__':
	
	args = argument_parse()
	generate_aa_table(args.chr)
