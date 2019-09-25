#!/usr/bin/env python

import pdb
from numpy import mean, isnan, percentile, array

def calculate_replication_timing(window_size):
	
	with open('./hg18_replication_timing_%i_bp_windows.bed' % (window_size), 'w') as output_file:
	
		output_str = ''
		
		with open('./Koren_et_al_Table_S2.txt', 'rU') as open_rep_file:
		
			curr_chr = None
			curr_start_pos = None
			curr_scores = []
			line_counter = 0
		
			for line in open_rep_file:
			
				chr, pos, score = line[:-1].split('\t')
				if chr == '23':
					chr = 'X'
				elif chr == '24':
					chr = 'Y'
				chr = 'chr' + chr
				pos = int(float(pos))
			
				if chr != curr_chr: # new chromosome: increment segment
					if curr_chr is not None:
						output_str += '%s\t%i\t%i\t%f\n' % (curr_chr, curr_start_pos, curr_start_pos + window_size, mean(curr_scores))
					curr_chr = chr
					curr_start_pos = pos
					curr_scores = [float(score)]
				elif pos - curr_start_pos > window_size: # close window, increment segments
					output_str += '%s\t%i\t%i\t%f\n' % (curr_chr, curr_start_pos, curr_start_pos + window_size, mean(curr_scores))
					line_counter += 1
					if line_counter == 5000:
						line_counter = 0
						output_file.write(output_str)
						output_str = ''
					curr_chr = chr
					curr_start_pos = pos
					curr_scores = [float(score)]
				else:
					curr_scores.append(float(score))
			
			output_str += '%s\t%i\t%i\t%f\n' % (curr_chr, curr_start_pos, curr_start_pos + window_size, mean(curr_scores))
			output_file.write(output_str)
	
	return None

def calculate_replication_timing_breakpoints(window_size, n_subsets):
	
	score_array = []
	
	with open('./Koren_et_al_Table_S2.txt', 'rU') as open_rep_file:
	
		curr_chr = None
		curr_start_pos = None
		curr_scores = []
		line_counter = 0
	
		for line in open_rep_file:
		
			chr, pos, score = line[:-1].split('\t')
			if chr == '23':
				chr = 'X'
			elif chr == '24':
				chr = 'Y'
			chr = 'chr' + chr
			pos = int(float(pos))
		
			if chr != curr_chr: # new chromosome: increment segment
				if curr_chr is not None:
					if not isnan(mean(curr_scores)):
						score_array.append((curr_chr, curr_start_pos, mean(curr_scores)))
				curr_chr = chr
				curr_start_pos = pos
				curr_scores = [float(score)]
			elif pos - curr_start_pos > window_size: # close window, increment segments
				if not isnan(mean(curr_scores)):
					score_array.append((curr_chr, curr_start_pos, mean(curr_scores)))
				curr_chr = chr
				curr_start_pos = pos
				curr_scores = [float(score)]
			else:
				curr_scores.append(float(score))
		score_array.append((curr_chr, curr_start_pos, mean(curr_scores)))
	
# 	score_array.sort(key=lambda x: x[2])
	breakpoints = percentile(array([x[2] for x in score_array]), q=[100 * i / n_subsets for i in range(1, n_subsets)])
	
	score_array_breakpoints = [[] for _ in range(n_subsets)]
	for seg in score_array:
		for n, i in enumerate(breakpoints):
			if seg[2] < i:
				score_index = n
				break
		else:
			score_index = len(breakpoints)
		score_array_breakpoints[score_index].append(seg)
		
	def sort_seg_function(seg):
		if seg[0] == 'chrX':
			return [23, seg[1]]
		elif seg[0] == 'chrY':
			return [24, seg[1]]
		else:
			return [int(seg[0][3:]), seg[1]]
	
	for n, curr_score_array in enumerate(score_array_breakpoints):
		with open('./replication_timing_compartments/hg18_replication_timing_q%i.bed' % n, 'w') as output_file:
			curr_score_array.sort(key=sort_seg_function)
			output = '\n'.join(['\t'.join([x[0], str(x[1]), str(x[1] + window_size)]) for x in curr_score_array])
			output_file.write(output)
	
	return None


if __name__ == '__main__':
	
	calculate_replication_timing_breakpoints(20000, 4)