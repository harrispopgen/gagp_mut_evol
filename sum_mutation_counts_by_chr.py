#!/usr/bin/env python

import argparse
import glob
import pdb
import os
from common import str2bool, generate_nmer_mutation_list, generate_nmer_list

species = ['Homo', 'Pan_paniscus', 'Pan_troglodytes', 'Gorilla', 'Pongo_abelii', 'Pongo_pygmaeus']
chrs = ['chr' + str(n) for n in range(1, 23)]

def sum_mutation_counts(path, snp_type, nmer, delete_files=False):
	
	nmer_mutations = generate_nmer_mutation_list(nmer)
	for s in species:
		indiv_mut_dict = {}
		species_mut_dict = dict(zip(nmer_mutations, [0 for _ in nmer_mutations]))
		indivs = []
		for c in chrs:
			species_chr_filename = path + '%s_%s_%s_species_%imer_counts.txt' % (s, c, snp_type, nmer)
			indiv_chr_filename = path + '%s_%s_%s_indiv_%imer_counts.txt' % (s, c, snp_type, nmer)
			if os.path.isfile(species_chr_filename):
				with open(species_chr_filename, 'rU') as open_species_file:
					header = open_species_file.readline()[:-1].split('\t')
					counts = open_species_file.readline()[:-1].split('\t')
					for x, y in enumerate(header):
						species_mut_dict[y] += int(counts[x])
			else:
				print("Unable to find %s" % (species_chr_filename))
			if os.path.isfile(indiv_chr_filename):
				with open(indiv_chr_filename, 'rU') as open_indiv_file:
					header = open_indiv_file.readline()[:-1].split('\t')
					for l in open_indiv_file.readlines():
						l_split = l[:-1].split('\t')
						curr_indiv = l_split[0]
						if curr_indiv not in indiv_mut_dict: # First time seeing indiv: initialize dictionary
							indivs.append(curr_indiv)
							indiv_mut_dict[curr_indiv] = dict(zip(nmer_mutations, [0 for _ in nmer_mutations]))
						for x, y in enumerate(header[1:]):
							indiv_mut_dict[curr_indiv][y] += int(l_split[x + 1]) # enumerate initializes at 0 but I'm starting starting with element #1
			else:
				print("Unable to find %s" % (indiv_chr_filename))
			if delete_files:
				if os.path.isfile(species_chr_filename):
					os.remove(species_chr_filename)
				if os.path.isfile(indiv_chr_filename):
					os.remove(indiv_chr_filename)
		with open(path + '%s_%s_species_%imer_counts.txt' % (s, snp_type, nmer), 'w') as outfile_species:
			outfile_species.write('\t'.join(sorted(nmer_mutations)))
			outfile_species.write('\n')
			outfile_species.write('\t'.join([str(species_mut_dict[m]) for m in sorted(nmer_mutations)]))
			outfile_species.write('\n')
		with open(path + '%s_%s_indiv_%imer_counts.txt' % (s, snp_type, nmer), 'w') as outfile_indiv:
			outfile_indiv.write('indiv\t')
			outfile_indiv.write('\t'.join(sorted(nmer_mutations)))
			outfile_indiv.write('\n')
			for i in indivs:
				outfile_indiv.write(i + '\t')
				outfile_indiv.write('\t'.join([str(indiv_mut_dict[i][m]) for m in sorted(nmer_mutations)]))
				outfile_indiv.write('\n')

def argument_parse():
	
	parser = argparse.ArgumentParser()
	parser.add_argument("path", help="Pathway to mutation count directory")
	parser.add_argument("snp_type", help="Name of types of SNPs to sum across. Originally the name of the bed file that included the genomic regions in which the SNPs fall.")
	parser.add_argument("nmer", help="Nmer (int)", type=int)
	parser.add_argument("delete_files", help="Delete chromosome files?", type=str2bool, nargs='?', default=False)
	args = parser.parse_args()
	return args

if __name__ == '__main__':
	
	args = argument_parse()
	sum_mutation_counts(args.path, args.snp_type, args.nmer, args.delete_files)
	