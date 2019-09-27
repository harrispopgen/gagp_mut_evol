#!/usr/bin/env python

import pdb
import gzip

def load_bed_file(infile):
	
	seg_dict = {}
	with open(infile, 'rU') as open_file:
		for l in open_file:
			l_split = l.rstrip('\n').split('\t')
			if l_split[0] not in seg_dict:
				seg_dict[l_split[0]] = []
			seg_dict[l_split[0]].append((int(l_split[1]), int(l_split[2])))
	
	seg_dict = {chr: seg_dict[chr] for chr in seg_dict if chr in ['chr' + str(x) for x in range(1, 23)]}
	
	return seg_dict

def load_hmc_data(infile):
	
	pos_dict = {}
	
	with gzip.open(infile, 'rt') as open_file:
		open_file.readline()
		for l in open_file:
			l_split = l.rstrip('\n').split('\t')
			if l_split[0] not in pos_dict:
				pos_dict[l_split[0]] = []
			pos_dict[l_split[0]].append(int(l_split[2]) - 1) # IT IS ONE-INDEXED
	
	pos_dict = {chr: pos_dict[chr] for chr in pos_dict if chr in ['chr' + str(x) for x in range(1, 23)]}
	
# 	pos_dict = {chr:sorted(pos_dict[chr]) for chr in pos_dict}
	
	return pos_dict

# maybe later I'll subset into ervs that have high hmC to C ratios, right now i'm just ordering by #hmc sites
def fraction_hmc(compartment, hmc_pos, n_subsets = 4):
	
	fraction_hmc_list = []
	for chr in compartment:
		print(chr)
		curr_list = [] # in case it gets too big
		for curr_start, curr_end in compartment[chr]:
			frac_hmc = float(len([x for x in hmc_pos[chr] if x < curr_end and x >= curr_start]))/ float(curr_end - curr_start)
			curr_list.append((chr, curr_start, curr_end, frac_hmc))
		fraction_hmc_list += curr_list
	
	return sorted(fraction_hmc_list, key=lambda x: x[3])

def write_fraction_hmc_list_to_output(outfile, seg_list):
	
	with open(outfile, 'w') as open_file:
		open_file.write('\n'.join(['\t'.join([x[0], str(x[1]), str(x[2])]) for x in seg_list]))
		open_file.write('\n')
	
	return None

def argument_parse():
	
	parser = argparse.ArgumentParser()
	parser.add_argument("--erv_file", help="ERV compartment bed file")
	parser.add_argument("--hmc_file", help=".txt.gz 5hmC data file, from Yu et al., 2012. GSM882245_H1.hmC_sites.FDR_0.0502.hg18.txt.gz")
	args = parser.parse_args()
	return args


if __name__ == '__main__':
	
	args = argument_parse()
	if args.erv_file:
		erv_compartment = args.erv_file
	else:
		erv_compartment = load_bed_file('../erv_data/hg18_all_erv.bed')
	if args.hmc_file:
		hmc_data = load_hmc_data(args.hmc_file)
	else:
		hmc_data = load_hmc_data('./ftp.ncbi.nlm.nih.gov/geo/series/GSE36nnn/GSE36173/suppl/GSM882245_H1.hmC_sites.FDR_0.0502.hg18.txt.gz')
	frac_hmc = fraction_hmc(erv_compartment, hmc_data)
	hmc_low = [x for x in frac_hmc if x[3] == 0.0]
	hmc_high = [x for x in frac_hmc if x[3] > 0.0]
	hmc_low = sorted(hmc_low, key=lambda x: (x[0], x[1]))
	hmc_high = sorted(hmc_high, key=lambda x: (x[0], x[1]))
	write_fraction_hmc_list_to_output('./hg18_all_erv_hmc_low.bed', hmc_low)
	write_fraction_hmc_list_to_output('./hg18_all_erv_hmc_high.bed', hmc_high)
	