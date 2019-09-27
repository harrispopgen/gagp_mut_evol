#!/usr/bin/env python

import pdb

def parse_refgene_to_exon_bed(infile):
	
	exon_list = []
	
	with open(infile, 'r') as open_file:
		
		open_file.readline()
		for l in open_file:
			l_split = l.rstrip('\n').split('\t')
			curr_chr = l_split[2]
			exon_starts = l_split[9].rstrip(',').split(',')
			exon_ends = l_split[10].rstrip(',').split(',')
			exon_frames = l_split[15].rstrip(',').split(',')
			curr_exon_list = [(x[0], x[1]) for x in zip(exon_starts, exon_ends, exon_frames) if x[2] in '012']
# 			if len(curr_exon_list) != len(exon_starts): # This works okay
# 				pdb.set_trace()
			exon_list += [(curr_chr, x[0], x[1]) for x in curr_exon_list]
	
	return exon_list

def list_to_bed(exon_list, outfile):
	
	with open(outfile, 'w') as open_outfile:
		
		open_outfile.write('\n'.join(['\t'.join(x) for x in exon_list]))
		open_outfile.write('\n')
	
	return None

def argument_parse():
	
	parser = argparse.ArgumentParser()
	parser.add_argument("--refgene", help="refGene filename")
	args = parser.parse_args()
	return args


if __name__ == '__main__':
	args = argument_parse()
	# Need the refGene output from UCSC table browser, hg18 build
	if args.refgene:
		refgene_file = args.refgene
	else:
		refgene_file = './data/hg18/hg18_refGene.txt'
	exon_list = parse_refgene_to_exon_bed(refgene_file)
	list_to_bed(exon_list, './data/hg18_refGene_exons.bed')