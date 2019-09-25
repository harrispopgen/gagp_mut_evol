#!/usr/bin/env python

import math
import pdb
import numpy as np
import itertools
import argparse
import pandas

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
autosomes = ['chr' + str(n) for n in range(1, 23)]
species_list = ['Homo', 'Gorilla', 'Pongo_abelii', 'Pongo_pygmaeus', 'Pan_paniscus', 'Pan_troglodytes']

def rev_comp(seq):
	return "".join(complement.get(base, base) for base in reversed(seq))

def occurrences(seq, sub):
    count = start = 0
    while True:
        start = seq.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count

def generate_nmer_list(nmer):
	
	upstream_n = int(nmer/2)
	downstream_n = int(nmer/2)	
	
	upstream_nmers = [''.join(x) for x in itertools.product('ACGT', repeat=upstream_n)]
	downstream_nmers = upstream_nmers.copy()
	central_nuc = 'AC'
	
	nmer_list = [x + y + z for x in upstream_nmers for y in central_nuc for z in downstream_nmers]
	
	return nmer_list

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# Implements algorithm for HWE from Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005)
# See http://csg.sph.umich.edu/abecasis/Exact/snp_hwe.r for original R code
def hwe(obs_hets, obs_hom1, obs_hom2, option = 'two-sided'):
	
	if (obs_hom1 + obs_hom2 + obs_hets) == 0:
		return 1.
	if obs_hom1 < 0 or obs_hom2 < 0 or obs_hets < 0:
		raise ValueError("Cannot calculate HWE: negative allele values")
	
	# Total number of genotypes
	n = obs_hom1 + obs_hom2 + obs_hets
	
	# Rare, common homozygotes
	obs_homr = min(obs_hom1, obs_hom2)
	obs_homc = max(obs_hom1, obs_hom2)
	
	# Number of rare allele copies
	rare = obs_homr * 2 + obs_hets
	
	# Initialize probability array
	probs = np.zeros(rare + 1)
	
	# Find midpoint of distribution
	mid = math.floor(rare * (2 * n - rare) / (2 * n))
	if mid % 2 != rare % 2:
		mid += 1
	
	probs[mid] = 1
	mysum = 1
	
	# Calculate probabilities from midpoint down
	curr_hets = mid
	curr_homr = (rare - mid) / 2
	curr_homc = n - curr_hets - curr_homr
	
	while curr_hets >= 2:
		probs[curr_hets - 2] = probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
		mysum += probs[curr_hets - 2]
		
		# 2 fewer heterozygotes -> add 1 rare homozygote, 1 common homozygote
		curr_hets -= 2
		curr_homr += 1
		curr_homc += 1
	
	# Now calculate probabilities from midpoint up
	curr_hets = mid
	curr_homr = (rare - mid) / 2
	curr_homc = n - curr_hets - curr_homr
	
	while curr_hets <= (rare - 2):
		probs[curr_hets + 2] = probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
		mysum += probs[curr_hets + 2]
		
		# Add 2 hets and subtract one of each homozygote
		curr_hets += 2
		curr_homr -= 1
		curr_homc -= 1
	
	# P-val calculation
	target = probs[obs_hets]
	
	if option == 'two-sided':
		p = min(1.0, sum(probs[probs <= target]) / mysum)
		return p
	elif option == 'excess':
		p_hi = min(1.0, sum(probs[obs_hets: rare + 1]) / mysum)
		return p_hi
	elif option == 'depletion':
		p_lo = min(1.0, sum(probs[0:obs_hets + 1]) / mysum)
		return p_lo
	else:
		raise ValueError("Unrecognized option!")
	
	return None

def filter_var(var):
	
	gts = [i['GT'] for i in var.samples.values()]
	obs_hom1 = gts.count((0,0))
	obs_hom2 = gts.count((1,1))
	obs_hets = gts.count((0,1)) # GAGP is unphased data; all heterozygotes should be (0,1) rather than (1,0)
# 	pdb.set_trace()
	return(len(var.alleles) == 2 and # Biallelic
		   var.info['AC'][0] > 1 and # No singletons
		   (var.info['AN'] - var.info['AC'][0]) > 1 and # No singletons (other allele)
		   hwe(obs_hets, obs_hom1, obs_hom2, option='excess') > 0.05) # HW excess het filter

def generate_nmer_mutation_list(nmer):
	
	upstream_n = int(nmer/2)
	downstream_n = int(nmer/2)
	
	# Setting up nmer lists
	nmer_mutations = []
	upstream_nmers = [''.join(x) for x in itertools.product('ACGT', repeat=upstream_n)]
	downstream_nmers = upstream_nmers.copy()
	central_nuc = ['A', 'C']
	nmer_mutations = ([x + y + z + _ + a for x in upstream_nmers for y in 'A' for z in downstream_nmers for _ in '>' for a in 'CGT'] +
					  [x + y + z + _ + a for x in upstream_nmers for y in 'C' for z in downstream_nmers for _ in '>' for a in 'AGT'])
	
	return nmer_mutations
	
def load_hg18_ref(species, chr, exclude_recurrent = False):
	# Loads up the reference sequence for a specific chromosome. Each species in the GAGP has a special hg18 reference with their ancestral alleles
	hg18_ref = ''
	if exclude_recurrent:
		hg18_ref_filename = './hg18_references_with_gagp_ancestral_alleles_exclude_recurrent/' + species + '_' + chr + '.fa'
	else:
		hg18_ref_filename = './hg18_references_with_gagp_ancestral_alleles/' + species + '_' + chr + '.fa'
	with open(hg18_ref_filename) as open_hg18_ref:
		open_hg18_ref.readline()
		hg18_ref = open_hg18_ref.read()
	
	return hg18_ref

def normalize_indiv_mutation_matrix(mutation_matrix, nmer_content):
	return np.divide(
			(
			 (mutation_matrix.transpose() / 
				mutation_matrix.sum(1))
			 ).transpose(), 
			(
			 (nmer_content[[x[:3] for x in mutation_matrix.columns]] / 
			 	nmer_content.sum(1).iloc[0])
			 	)
			 )

if __name__ == '__main__':
	pdb.set_trace()
	print(hwe(11, 5, 84)) # Should be .999952, .000919