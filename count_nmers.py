#!/usr/bin/env python

from Bio import SeqIO
import itertools
import argparse
import pdb
from common import filter_var, generate_nmer_mutation_list, generate_nmer_list, rev_comp, complement, autosomes, occurrences


def count_nmers(seq, nmer_list):
	
	nmer_dict = dict(zip(nmer_list, [0 for _ in nmer_list]))

	for i in nmer_list: # Pretty sure I'd always want to include reverse complements
		nmer_dict[i] = occurrences(seq.upper(), i)
		nmer_dict[i] += occurrences(seq.upper(), rev_comp(i))
	
	return nmer_dict

def fasta_to_seg_nmer(filename, outfile_name, nmer_list):
	
	output = '%s\t%s\t' % ('N', 'length')
	output += '\t'.join(nmer_list)
	
	with open(filename, 'rU') as f_open:
		for record in SeqIO.parse(f_open, 'fasta'):
			curr_nmer_count = count_nmers(record.seq)
			if int(record.id) % 1000 == 0:
				print(record.id)
			output += '\n%s\t%i\t' % (record.id, len(record.seq))
			output += '\t'.join([str(curr_nmer_count[n]) for n in nmer_list])
	
	outfile = open(outfile_name, 'w')
	outfile.write(output)
	outfile.close()

	return None

def get_hg18_chr_ref(chr):
	with open('/net/harris/vol1/data/hg18/' + chr + '.fa', 'rU') as open_file:
		return ''.join(open_file.read().split()[1:])

def bed_to_total_nmer(bedfile, nmer_list, autosomal_only = True):
	
	nmer = len(nmer_list[0])
	nmer_dict = dict(zip(nmer_list, [0 for _ in nmer_list]))
	total_len = 0
	with open(bedfile, 'rU') as open_file:
		curr_chr = ''
		curr_ref = None
		counter = 0
		for line in open_file:
			line_split = line.split('\t')
			if autosomal_only: 
				if line_split[0] not in autosomes:
					continue
			if line_split[0] != curr_chr:
				curr_ref = get_hg18_chr_ref(line_split[0])
				curr_chr = line_split[0]
				print(curr_chr)
			line_split[1:3] = [int(x) for x in line_split[1:3]]
			if line_split[2] - line_split[1] < 3:
				continue
			seq = curr_ref[line_split[1]:line_split[2]].upper() # REMEMBER! Both .bed files and python are 0-indexed!
			for n in range(len(seq) - (nmer-1)):
				try:
					if seq[n + int(nmer/2)] in 'GT':
						nmer_dict[rev_comp(seq[n:n + nmer])] += 1
					else:
						nmer_dict[seq[n:n + nmer]] += 1
				except KeyError: # Catching N et al.
					continue
	
	return nmer_dict

def bed_to_total_nmer_with_negative(bedfile, nmer_list, negative_bedfile, autosomal_only = True):
	# The negative_bedfile MUST be in the same order as the bedfile. This happens normally when using liftover
	# Also the negative_bedfile is 
	nmer = len(nmer_list[0])
	nmer_dict = dict(zip(nmer_list, [0 for _ in nmer_list]))
	total_len = 0
	with open(bedfile, 'rU') as open_file:
		with open(negative_bedfile, 'rU') as open_neg_file:
			curr_chr = ''
			curr_ref = None
			open_neg_file.readline() # incrementing past the comment
			neg_seg = open_neg_file.readline().split('\t')
			for line in open_file:
				line_split = line.split('\t')
				if autosomal_only: 
					if line_split[0] not in autosomes:
						continue
				if neg_seg is [''] or line_split[0:3] == neg_seg[0:3]:
					open_neg_file.readline()
					neg_seg = open_neg_file.readline().split('\t')
					continue
				if line_split[0] != curr_chr:
					curr_ref = get_hg18_chr_ref(line_split[0])
					curr_chr = line_split[0]
					while line_split[0] != neg_seg[0] and neg_seg[0] is not ['']:
						open_neg_file.readline()
						neg_seg = open_neg_file.readline().split('\t')
					print(curr_chr)
				line_split[1:3] = [int(x) for x in line_split[1:3]]
				if line_split[2] - line_split[1] < 3:
					continue
				seq = curr_ref[line_split[1]:line_split[2]].upper() # REMEMBER! Both .bed files and python are 0-indexed!
				for n in range(len(seq) - (nmer-1)):
					try:
						if seq[n + int(nmer/2)] in 'GT':
							nmer_dict[rev_comp(seq[n:n + nmer])] += 1
						else:
							nmer_dict[seq[n:n + nmer]] += 1
					except KeyError: # Catching N et al.
						continue
	
	return nmer_dict

def bed_to_seg_nmer(bedfile, nmer_list, outfile_name, autosomal_only = True):
	
	output = '%s\t%s\t%s\t%s\t' % ('chr', 'start', 'stop', 'length')
	output += '\t'.join(nmer_list)
	
	with open(bedfile, 'rU') as f_open:
		curr_chr = ''
		curr_ref = None
		counter = 0
		for line in f_open:
			line_split = line.split('\t')
			if autosomal_only: 
				if line_split[0] not in autosomes:
					continue
			if line_split[0] != curr_chr:
				curr_ref = get_hg18_chr_ref(line_split[0])
				curr_chr = line_split[0]
				print(curr_chr)
			line_split[1:3] = [int(x) for x in line_split[1:3]]
			if line_split[2] - line_split[1] < 3:
				continue
			seq = curr_ref[line_split[1]:line_split[2]].upper() # REMEMBER! Both .bed files and python are 0-indexed!
			curr_nmer_count = count_nmers(seq, nmer_list)
			output += '\n%s\t%i\t' % ('\t'.join([str(i) for i in line_split[0:3]]), len(seq))
			output += '\t'.join([str(curr_nmer_count[n]) for n in nmer_list])
	
	outfile = open(outfile_name, 'w')
	outfile.write(output)
	outfile.close()

def nmer_dict_to_outfile(nmer_dict, nmer_list, filename):
	
	with open(filename, 'w') as outfile:
		
		outstr = '\t'.join(nmer_list)
		outstr += '\n'
		outstr += '\t'.join([str(nmer_dict[n]) for n in nmer_list])
		outstr += '\n'
		outfile.write(outstr)
	return

def argument_parse():
	
	parser = argparse.ArgumentParser()
	parser.add_argument("analysis", help="bed_to_seg_nmer, bed_to_total_nmer, bed_to_total_nmer_with_negative, or fasta_to_seg_nmer")
	parser.add_argument("file", help="Path to a bed or fasta file containing chromosomal regions (0-index) to analyze or sequences, respectively")
	parser.add_argument("nmer", help="Nmer (must be odd number)", type=int)
	parser.add_argument("outfile", help="Path to outfile")
	parser.add_argument("--negative_bed_file", help="Path to a bed file containing segments to exclude from analysis", default=None)
	args = parser.parse_args()
	return args

if __name__ == '__main__':
	
	args = argument_parse()
	if args.nmer % 2 != 1:
		raise ValueError("nmer must be an odd number")
	nmer_list = generate_nmer_list(args.nmer)
	if args.analysis == "bed_to_seg_nmer":
		bed_to_seg_nmer(args.file, nmer_list, args.outfile, autosomal_only = True)
	elif args.analysis == "bed_to_total_nmer":
		nmer_dict_to_outfile(bed_to_total_nmer(args.file, nmer_list, autosomal_only = True), nmer_list, args.outfile)
	elif args.analysis == "fasta_to_seg_nmer":
		fasta_to_seg_nmer(args.file, args.outfile, nmer_list)
	elif args.analysis == 'bed_to_total_nmer_with_negative':
		nmer_dict_to_outfile(bed_to_total_nmer_with_negative(args.file, nmer_list, args.negative_bed_file, autosomal_only = True), nmer_list, args.outfile)
	else:
		raise ValueError("Unrecognized analysis. Supported analyses are bed_to_seg_nmer, bed_to_total_nmer, or fasta_to_seg_nmer")
	
