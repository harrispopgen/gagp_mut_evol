#!/usr/bin/env python

import gzip
import argparse
import pdb

def generate_references(species, chr):
	
	hg18_ref = '' # This will be the hg18 reference genome with ancestral alleles for this species
	ancestral_allele_path = './ancestral_allele_tables/' # Need to have generated all ancestral allele tables
	genus_column = None # Which column in the table is the genus of the species?
	if species == 'Homo':
		genus_column = 2
	elif species in ['Pan_troglodytes', 'Pan_paniscus']:
		genus_column = 3
	elif species == 'Gorilla':
		genus_column = 4
	elif species in ['Pongo_abelii', 'Pongo_pygmaeus']:
		genus_column = 5
	else:
		print("Unrecognized species:", species)
		return None
	chr_filename = ancestral_allele_path + chr + '.txt'
	
	# Need hg18 as separate chr .fa files
	with open('./data/hg18/' + chr + '.fa', 'rU') as open_ref:
		open_ref.readline()
		hg18_ref = ''.join([l.rstrip('\n') for l in open_ref.readlines()])
	
# 	print(len(hg18_ref)) # This should be the same after ancestral allele edits
	
	# Need this output file within directory; will hold the species-specific references
	with open('./hg18_references_with_gagp_ancestral_alleles/' + species + '_' + chr + '.fa', 'w') as output_file_open:
		with open(chr_filename, 'rU') as ancestral_alleles_open:
			output_str = '>' + chr + '\n'
			aa_header = ancestral_alleles_open.readline().rstrip('\n').split('\t')
			curr_index = -1
			old_index = None
			for line in ancestral_alleles_open:
				line_split = line.rstrip('\n').split('\t')
				old_index = curr_index
				curr_index = int(line_split[1]) - 1 # Getting the current variant's position. Remember that these positions come from a VCF! It's 1-indexed!
				ref_allele = hg18_ref[curr_index].upper() # Getting the actual hg18 reference allele
				line_split[2:] = [ref_allele if x == '-' else x for x in line_split[2:]] # Replacing all '-' with reference allele (i.e. no variant called for genus)
				output_str += hg18_ref[old_index + 1:curr_index] # Note: old_index + 1 accounts for the -1 initial value for curr_index!
				# Filtering out problematic sites
				if 'N' in line_split[2:]: # Multiallelic sites or other issues
					output_str += 'X'
				elif len([x for x in line_split[2:] if len(x) > 1]) > 1: # Multiple genus segregating
					output_str += 'X'
				elif len(set(line_split[2:])) > 2: # > 2 different fixed sites, or fixed_alt + segregating + fixed_reference
					output_str += 'X'
				# Not a problematic site; identify ancestral allele
				# NOTE: this next if statement won't catch situations where, e.g., homo == pongo while pan == gorilla.
				elif len(line_split[genus_column]) == 1: # Fixed site; replace in species reference
					output_str += line_split[genus_column]
				else: # Segregating site; take a parsimony approach
					other_alleles = set(line_split[2:genus_column] + line_split[genus_column + 1:])
					if len(other_alleles) == 1:
						output_str += list(other_alleles)[0]
					else: # this should have been caught??
						pdb.set_trace()
				if len(output_str) > 500000:
					output_file_open.write(output_str)
					output_str = ''
			output_str += hg18_ref[curr_index + 1:]
			output_file_open.write(output_str)
	
	return None

def argument_parse():
	
	parser = argparse.ArgumentParser()
	parser.add_argument("species", help="Name of species, as defined by GAGP")
	parser.add_argument("chr", help="Chromosome (e.g. 'chr10', 'chrY')")
	args = parser.parse_args()
	return args

		
if __name__ == '__main__':
	
	args = argument_parse()
	generate_references(args.species, args.chr)