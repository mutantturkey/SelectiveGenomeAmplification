#!/usr/bin/env python2.7
import sys
import os

import argparse

from multiprocessing import Pool
from multiprocessing import cpu_count

from subprocess import Popen
from subprocess import PIPE

from itertools  import combinations

import numpy as np

fg_mers = {}
bg_mers = {}

heterodimer_dic = {}

seq_ends = []

fg_genome_length = 0
bg_genome_length = 0

output_file = ""

# import our variables
cpus             = int(os.environ.get("cpus", cpu_count()))
debug            = os.environ.get("debug", False)
max_select       = int(os.environ.get("max_select", 15))
max_check        = int(os.environ.get("max_check", 35))
max_mer_distance = int(os.environ.get("max_mer_distance", 5000))
max_consecutive_binding = int(os.environ.get("max_consecutive_binding", 4))
primer_weight = float(os.environ.get("primer_weight", 0))
score_str = os.environ.get("score_func", None)

if score_str is not None:
	score_func = compile('mer_score = ' + score_str, '<string>', 'exec')
else:
	score_func = compile('mer_score = (nb_primers**primer_weight) * (fg_mean_dist * fg_std_dist) / bg_ratio', '<string>', 'exec')
	score_str = '(nb_primers**primer_weight) * (fg_mean_dist * fg_std_dist) / bg_ratio'

def get_max_consecutive_binding(mer1, mer2):
	'''
	Return the maximum number of consecutively binding mers
	when comparing two different mers, using the reverse compliment.
	'''

	binding = { 'A': 'T', 
							'T': 'A',
							'C': 'G', 
							'G': 'C',	
							'_':  False
						}

  # Swap variables if the second is longer than the first
	if len(mer2) > len(mer1):
		mer1, mer2 = mer2, mer1
	
	# save the len because it'll change when we do a ljust
	mer1_len = len(mer1)
	# reverse mer2,
	mer2 = mer2[::-1]
	# pad mer one to avoid errors
	mer1 = mer1.ljust(mer1_len + len(mer2), "_")

	max_bind = 0
	for offset in range(mer1_len):
		consecutive = 0
		for x in range(len(mer2)):
			if binding[mer1[offset+x]] == mer2[x]:
				consecutive += 1
				if consecutive > max_bind:
					max_bind = consecutive
			else:
				consecutive = 0

	return max_bind

def size_iterator(mer_arr):
	for mer_len in set([ len(it) for it in mer_arr]):
		mers = filter(lambda x: len(x) == mer_len, mer_arr)
		yield mers

def populate_locations(mers, mer_dic, input_fn, length):
	''' 
	Run the kmer_locations command, and parse in the integers that are output
	by the command, and add it to mers[mer] 

	We also do the reverse compliment, so we subtract the length of the genome for it
	'''
	import tempfile

	# write our mers out to a fifi
	for n_mers in size_iterator(mers):
		_, merlist_fn = tempfile.mkstemp()
		merlist_fh = open(merlist_fn, 'w')
		for mer in n_mers:
			merlist_fh.write(mer + '\n')

		merlist_fh.flush()
		# add our merlist fn to our command
		cmd = "kmer_locations -k " + str(len(n_mers[0])) +  " -c  -m \"" + merlist_fn + "\" -l -i \"" + input_fn + "\""

		if(debug):
			print "Executing", cmd

		cmd_pipe = Popen(cmd, stdout=PIPE, shell=True)
		for line in cmd_pipe.stdout:
			line_s = line[:-1].split(" ")
			if len(line_s) == 2:
				mer_dic[line_s[0]].append(int(line_s[1]))
			elif len(line_s) == 3:
				mer_dic[line_s[0]].append(length - int(line_s[1]))

		if cmd_pipe.wait() is not 0:
			print "Executing", cmd, "failed"
			exit()

		merlist_fh.close()


def load_end_points(fn):
	''' get all the points of the end of each sequence in a sample '''

	end_points = [0]

	cmd = "sequence_end_points < " + fn

	if debug:
		print "loading sequence end points"
		print "executing: " + cmd

	points_fh = Popen(cmd, stdout=PIPE, shell=True)

	for line in points_fh.stdout:
		end_points.append(int(line))
	

	if points_fh.wait() is not 0:
		print "executing", cmd, "failed"

	return end_points

def get_length(fn):
	''' get length of a genome ( number of base pairs )'''

	cmd = 'grep "^>" ' + fn + " -v | tr -d '\\n' | wc -c"

	if debug:
		print "loading sequence end points"
		print "executing: " + cmd
	grep = Popen(cmd, stdout=PIPE, shell=True)

	length = grep.stdout.readline()

	length = int(length)

	if grep.wait() is not 0:
		print "executing", cmd, "failed"
		exit()


	return length

def load_heterodimer_dic(selected_mers):
	'''
	Generate a heterodimer dict which contains every possible combination of
	selected mers, so later we can check each combination without re-running the
	max_consecutive_binding function. 

	The stored values are Booleans, True if the result is larger than acceptable.

	'''
	for (mer1, mer2) in combinations(selected_mers, 2):
		res = get_max_consecutive_binding(mer1, mer2)
		heterodimer_dic[(mer1, mer2)] = res > max_consecutive_binding
		heterodimer_dic[(mer2, mer1)] = res > max_consecutive_binding
		# print res, heterodimer_dic[(mer1, mer2)]


def check_feasible(selected):
	total = 0
	for mer in selected:
		total += len(fg_mers[mer])

	if total is 0:
		print "something went wrong, no mers found in the foreground. Consider this a bug!"
		print fg_mers
		exit(1)

	if (fg_genome_length / (total + 1 )) > max_mer_distance:
		print "even if we select all top", max_select, "of", total
		print "mers disregarding any critera, and they were perfectly evenly spaced we would ",
		print "still not meet the right max mer distance < ", max_mer_distance, "requirement."
	
		print total, " / ", fg_genome_length, " = ", total / fg_genome_length 
		exit(1)

def percentage(part, whole, precision=2):

	part = float(part)
	whole = float(whole)
	if whole == 0:
		return '0.0%'

	percent = round(part / whole * 100, precision)
	if(percent < 10):
		percent = " " + str(percent)
	
	return str(percent) + "%"

def write_header(fh):
	fh.write("# variables used: max_select=" + str(max_select) + " max_check=" + str(max_check) + " max_mer_distance=" + str(max_mer_distance) + " max_consecutive_binding=" + str(max_consecutive_binding) + " primer_weight=" + str(primer_weight) + " ")
	fh.write("fg_genome_length=" + str(fg_genome_length) + " bg_genome_length=" + str(bg_genome_length) + "\n")
	fh.write("# scoring function: " + str(score_str) + "\n")
	fh.write("#nb_primers\tCombination\tScore\tFG_mean_dist\tFG_stdev_dist\tBG_ratio\n")

def write_result(fh, score_res):
	combination, score_val, fg_mean_dist, fg_stddev_dist, bg_ratio = score_res
	fh.write(str(len(combination)) + "\t")
	fh.write(' '.join(combination) + "\t")
	fh.write(str(score_val) + "\t")
	fh.write(str(fg_mean_dist) + "\t")
	fh.write(str(fg_stddev_dist) + "\t")
	fh.write(str(bg_ratio) + "\n")

def print_rejected(total_reject, total_checked, total_scored, excluded):
	print ""
	print "Reasons mers were excluded:\n"
	print "  max distance: " + percentage(excluded[0], total_reject) + " (" + str(excluded[0]) + ")"
	print "  mers overlap: " + percentage(excluded[1], total_reject) + " (" + str(excluded[1]) + ")"
	print "  heterodimers: " + percentage(excluded[2], total_reject) + " (" + str(excluded[2]) + ")"
	print ""
	print "  total combinations checked: ", total_checked
	print "  total combinations scored:  ", total_scored
	print "  percent rejected:  " + percentage(total_reject, total_checked) 
	print ""

def score_specific_combinations(mers):

	total_scored = 0
	total_checked = 0
	excluded = [0, 0, 0]

	p = Pool(cpus)

	fh = open(output_file, 'wb')
	write_header(fh)
	
	score_it = p.map(score_no_check, mers)
	for score_res in score_it:
		if type(score_res) is list:
			total_scored += 1
			write_result(fh, score_res)
		else:
			excluded[score_res] += 1;
	
	total_reject = len(mers) - total_scored
	print_rejected(total_reject, len(mers), total_scored, excluded)

def score_all_combinations(mers):
	import time

	total_scored = 0
	total_checked = 0
	excluded = [0, 0, 0]

	check_feasible(mers)

	p = Pool(cpus)

	fh = open(output_file, 'wb')
	write_header(fh)

	max_size = max_select 
	if len(mers) < max_select:
		max_size = len(mers) + 1

	for select_n in range(1, max_size + 1 ):
		print "scoring size ", select_n,
		t = time.time()
		scores_it = p.imap_unordered(score, combinations(mers, select_n), chunksize=8192)
		for score_res in scores_it:
			total_checked += 1
			if type(score_res) is list:
				total_scored += 1
				write_result(fh, score_res)
			else:
				excluded[score_res] += 1;

		print "size ", select_n, "took:", time.time()   - t

	total_reject = total_checked - total_scored
	
	print_rejected(total_reject, total_checked, total_scored, excluded)

	if(total_scored == 0):
		print "NO RESULTS FOUND"
		fh.write("NO RESULTS FOUND\n")
	

def check_max_distance(dist_arr):
	return any(dist > max_mer_distance for dist in dist_arr)

def score(combination):
	# input is a string of mers like 
	# ['ACCAA', 'ACCCGA', 'ACGTATA']

	# check if the combination passes our filters
	for combo in combinations(combination, 2):
		if heterodimer_dic[combo]:
			return 2

	for mer in combination:
		for other_mer in combination:
			if not mer == other_mer:
				if mer in other_mer:
					return 1

	# fg points
	fg_pts = []
	fg_dist = []

	for mer in combination:
		fg_pts = fg_pts + fg_mers[mer]

	fg_pts = fg_pts + seq_ends 

	fg_pts.sort()

	if fg_pts[0] is not 0:
		fg_pts = [0] + fg_pts

	# fg distances
	fg_dist = np.diff(fg_pts)

  # return without calculating scores if any objects are higher than our max distance
	if check_max_distance(fg_dist) is True:
		return 0

	# bg counts 
	bg_counts = 0

	for mer in combination:
		bg_counts += bg_mers[mer]

	if bg_counts <= 1:
		bg_counts = 1 

	bg_ratio = (bg_genome_length / bg_counts)


	nb_primers = len(combination)
	fg_mean_dist = np.mean(fg_dist)
	fg_std_dist = np.std(fg_dist)

	# this is our equation
	exec score_func

	return [combination, mer_score, fg_mean_dist, fg_std_dist, bg_ratio] 

def score_no_check(combination):
	# input is a string of mers like 
	# ['ACCAA', 'ACCCGA', 'ACGTATA']

	# fg points
	fg_pts = []
	fg_dist = []

	for mer in combination:
		fg_pts = fg_pts + fg_mers[mer]

	fg_pts = fg_pts + seq_ends 

	fg_pts.sort()

	if fg_pts[0] is not 0:
		fg_pts = [0] + fg_pts

	# fg distances
	fg_dist = np.diff(fg_pts)

	# bg counts 
	bg_counts = 0

	for mer in combination:
		bg_counts += bg_mers[mer]

	if bg_counts <= 1:
		bg_counts = 1 

	bg_ratio = (bg_genome_length / bg_counts)

	nb_primers = len(combination)
	fg_mean_dist = np.mean(fg_dist)
	fg_std_dist = np.std(fg_dist)

	# this is our equation
	exec score_func

	return [combination, mer_score, fg_mean_dist, fg_std_dist, bg_ratio] 



def initialize_mers(foreground, background, load_background=True):
	print "Calculating heterodimer distances"
	load_heterodimer_dic(fg_mers.keys())

	print "Populating foreground locations"
	populate_locations(fg_mers.keys(), fg_mers, foreground, fg_genome_length)

	if load_background:
		print "Populating background locations"
		populate_locations(fg_mers.keys(), bg_mers, background, bg_genome_length)

		for mer in bg_mers:
			bg_mers[mer] = len(bg_mers[mer])


def main():
	'''
	Basic worflow:

	Load Top X Selective Primers
	Populate Locations of Primers
	Score Combinations For All Sizes 

	'''
	global fg_genome_length
	global bg_genome_length
	global seq_ends
	global output_file
	global score_func

	parser = argparse.ArgumentParser(description="score mers")
	parser.add_argument("-f", "--foreground", help="foreground fasta file", required=True)
	parser.add_argument("-b", "--background", help="background fasta file", required=True)
	parser.add_argument("-o", "--output", help="output fasta with UIDs in the file", required=True)
	parser.add_argument("-s", "--selectivity-file", help="mer selectivity file generated by select_mers.py", required=False)
	parser.add_argument("-c", "--combination-file", help="a set of combinations you want to score", required=False)
	parser.add_argument("-m", "--mer-file", help="a set of you want to score all combinations of", required=False)
	parser.add_argument("-r", "--rescore-file", help="rescore an already scored output file", required=False)

	args = parser.parse_args()

	nb_flags = len(filter(lambda x: x is None, [args.combination_file, args.selectivity_file,args.mer_file, args.rescore_file]))
	if nb_flags != 3:
		if nb_flags == 4:
			parser.error("you must have at least one input file to score from [-s -c -m -r]")
		else:
			parser.error("you can only have one input file to score from" )
		exit(1)

	if not os.path.isfile(args.foreground):
		parser.error(args.foreground + " not found")
	if not os.path.isfile(args.background):
		parser.error(args.background + " not found")

	output_file =	args.output

	print "Getting genome length"
	fg_genome_length = get_length(args.foreground)
	bg_genome_length = get_length(args.background)

	print "fg_genome_length:", fg_genome_length
	print "bg_genome_length:", bg_genome_length

	print "Populating sequence end points"
	seq_ends = load_end_points(args.foreground)


	if args.selectivity_file is not None:
	  
		print "Scoring all mer combinations"

		selectivity_fh = open(args.selectivity_file, "r")
	
		# load our mer list into python
		mer_selectivity = selectivity_fh.readlines()
		mer_selectivity = [ x for x in mer_selectivity if not x.startswith('#')]

		# get the last max_check (it's sorted)
		if len(mer_selectivity) > max_check:
			selected_mers = mer_selectivity[-max_check:]
		else:
			selected_mers = mer_selectivity

		# load it into our fg and bg counts into their dictionaries
		for mer in selected_mers:
			split_mer = mer.split()
			fg_mers[split_mer[0]] = []
			bg_mers[split_mer[0]] = int(split_mer[2])

		selected_mers = [x.split()[0] for x in selected_mers]

		if len(selected_mers) is 0:
			print "no mers found."
			exit(1)

		# we already have our background counts
		initialize_mers(args.foreground, args.background, load_background=False)

		print "Scoring mer combinations"
		score_all_combinations(selected_mers)

	elif args.combination_file is not None:

		print "Scoring specific mer combinations"

		combinations = []

		combination_fh = open(args.combination_file, "r")
		for line in combination_fh:
			if line.startswith("#"):
				continue
			mers = line.split()
			combinations.append(mers)
			for mer in mers:
				fg_mers[mer] = []
				bg_mers[mer] = []

		if len(combinations) is 0:
			print "no combinations found."
			exit(1)
		initialize_mers(args.foreground, args.background)

		score_specific_combinations(combinations)

	elif args.mer_file is not None:
		print "Scoring all possible mer combinations from ", args.mer_file

		mer_fh = open(args.mer_file, "r")
		for mer in mer_fh:
			if mer.startswith("#"):
				continue
			mer = mer.strip()
			if(len(mer.split()) > 1):
				print "skipping line:", mer, "each line should contain only one mer"
				continue

			fg_mers[mer] = []
			bg_mers[mer] = []

		if len(fg_mers.keys()) is 0:
			print "no mers found."
			exit(1)

		initialize_mers(args.foreground, args.background)
		score_all_combinations(fg_mers.keys())

	elif args.rescore_file is not None:
		print "Scoring all mer combinations from ", args.rescore_file

		combinations = []

		score_fh = open(args.rescore_file, "r")
		for line in score_fh:
			if line.startswith("#"):
				continue
			split_line = line.split('\t')
			combination = split_line[1].split()
			combinations.append(combination)
			for mer in combination:
				fg_mers[mer] = []
				bg_mers[mer] = []
		
		if len(combinations) is 0:
			print "no combinations found."
			exit(1)

		initialize_mers(args.foreground, args.background)

		print "re-scoring scores file"

		score_specific_combinations(combinations)


	print "output file:", output_file

if __name__ == "__main__":
	sys.exit(main())
