#!/usr/bin/env python
import sys
import os

import argparse

from multiprocessing import Pool
from multiprocessing import cpu_count

from subprocess import Popen
from subprocess import PIPE

from itertools  import combinations

import numpy as np
import pdb

fg_mers = {}
bg_mers = {}

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


def populate_locations(selected_mers, mer_dic, input_fn):
	''' 
	Run the strstreamone command, and parse in the integers that are output
	by the command, and add it to mers[mer] 

	strstreamone just prints the location of a string argv[1] in stdout.

	We also do the reverse compliment, using tac and tr piped together.
	'''
	import tempfile


	cmds = []
	# strip file of header and delete newlines
	cmds.append("grep -v '^>' " + input_fn  +  " | tr -d '\\n' | strstream ")
	# reverse file, strip and delete newlines
	cmds.append("tac " + input_fn + \
							"| rev " \
							"| grep -v '>$' " \
							"| tr -d '\\n' " \
							"| tr [ACGT] [TGCA] | strstream ")
	
	for cmd in cmds:
		_, merlist_fn = tempfile.mkstemp()

		# write our mers out to a fifi
		merlist_fh = open(merlist_fn, 'w')
		for mer in selected_mers:
			merlist_fh.write(mer + '\n')

		merlist_fh.flush()
		# add our merlist fn to our command
		cmd = cmd + " " + merlist_fn

		strstream = Popen(cmd, stdout=PIPE, shell=True)
		for line in strstream.stdout:
			(mer, pos) = line.strip().split(" ")
			mer_dic[selected_mers[int(mer)]].append(int(pos))

		merlist_fh.close()


def check_feasible(selected):
	total = 0
	for mer in selected:
		total += len(fg_mers[mer])
	if (fg_genome_length / (total + 1 )) > max_mer_distance:
		print "even if we select all top ", max_select, 
		print "mers disregarding any critera, and they were perfectly evenly spaced we would ",
		print "still not meet the right max mer distance < ", max_mer_distance, "requirement."
	
		print total, " / ", fg_genome_length, " = ", total / fg_genome_length 
		exit()

def percentage(part, whole, precision=2):

	part = float(part)
	whole = float(whole)

	percent = round(part / whole * 100, precision)
	if(percent < 10):
		percent = " " + str(percent)
	
	return str(percent) + "%"

def write_header(fh):
	fh.write("Combination\tScore\tFG_mean_dist\tFG_stdev_dist\tBG_ratio\n")

def write_result(fh, score_res):
	combination, score_val, fg_mean_dist, fg_stddev_dist, bg_ratio = score_res
	fh.write(str(combination) + "\t")
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

	score_it = p.map(score, mers)
	for score_res in score_it:
		if type(score_res) is list:
			total_scored += 1
			write_result(fh, score_res)
		else:
			excluded[score_res] += 1;
	
	total_reject = len(mers) - total_scored
	print_rejected(total_reject, len(mers), total_scored, excluded)

def score_all_combinations(selected):
	import time

	total_scored = 0
	total_checked = 0
	excluded = [0, 0, 0]

	check_feasible(selected)

	p = Pool(cpus)

	fh = open(output_file, 'wb')
	write_header(fh)

	for select_n in range(1, max_select+1):
		print "scoring size ", select_n,
		t = time.time()
		scores_it = p.imap_unordered(score, combinations(selected, select_n), chunksize=8192)
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
	

heterodimer_dic = {}
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
	if any(dist > max_mer_distance for dist in fg_dist):
		#return [combination, "max", max(fg_dist)]
		return 0

	# bg counts 
	bg_counts = 0

	for mer in combination:
		bg_counts += bg_mers[mer]

	if bg_counts <= 1:
		bg_counts = 1 

	bg_ratio = (bg_genome_length / bg_counts)


	nb_primers = len(combination)
	fg_mean_dist =  np.mean(fg_dist)
	fg_std_dist = np.std(fg_dist)

	# this is our equation
	mer_score = (nb_primers**primer_weight) * (fg_mean_dist * fg_std_dist) / bg_ratio

	return [combination, mer_score, fg_mean_dist, fg_std_dist, bg_ratio] 

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
	
	return end_points

def get_length(fn):
	''' get length of a genome ( number of base pairs )'''

	cmd = 'grep "^>" ' + fn + " -v | tr -d '\\n' | wc -c"

	if debug:
		print "loading sequence end points"
		print "executing: " + cmd
	points_fh = Popen(cmd, stdout=PIPE, shell=True)

	length = points_fh.stdout.readline()

	length = int(length)

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

	parser = argparse.ArgumentParser(description="score mers")
	parser.add_argument("-f", "--foreground", help="foreground fasta file", required=True)
	parser.add_argument("-b", "--background", help="background fasta file", required=True)
	parser.add_argument("-o", "--output", help="output fasta with UIDs in the file", required=True)
	parser.add_argument("-s", "--selectivity-file", help="mer selectivity file generated by select_mers.py", required=False)
	parser.add_argument("-c", "--combination-file", help="a set of combinations you want to score", required=False)

	args = parser.parse_args()

	if args.selectivity_file is None and args.combination_file is None:
		print "you must either have a selectivity file or a combination file to score from"
		exit()
	if args.selectivity_file is not None and args.combination_file is not None:
		print "you can only select either a selectivity file or a combination file to score from"
		exit()

	output_file =	args.output

	print "Getting genome length"
	fg_genome_length = get_length(args.foreground)
	bg_genome_length = get_length(args.background)

	print "Populating sequence end points"
	seq_ends = load_end_points(args.foreground)


	if(args.selectivity_file is not None):
	  
		print "Scoring all mer combinations"

		selectivity_fh = open(args.selectivity_file, "r")
	
		# load our mer list into python
		mer_selectivity = selectivity_fh.readlines()

		# get the last max_check (it's sorted)
		selected_mers = mer_selectivity[-max_check:]

		# load it into our fg and bg counts into their dictionaries
		for mer in selected_mers:
			split_mer = mer.split()
			fg_mers[split_mer[0]] = []
			bg_mers[split_mer[0]] = int(split_mer[2])

		selected_mers = [x.split()[0] for x in selected_mers]

		print "Populating foreground locations"
		populate_locations(selected_mers, fg_mers, args.foreground)

		print "calculating heterodimer distances"
		load_heterodimer_dic(selected_mers)

		print "scoring mer combinations"
		score_all_combinations(selected_mers)
	
	else:
		print "Scoring specific mer combinations"

		combinations = []

		combination_fh = open(args.combination_file, "r")
		for line in combination_fh:
			mers = line.split()
			combinations.append(mers)
			for mer in mers:
				fg_mers[mer] = []
				bg_mers[mer] = []

		print "calculating heterodimer distances"
		load_heterodimer_dic(fg_mers.keys())

		print "Populating foreground locations"
		populate_locations(fg_mers.keys(), fg_mers, args.foreground)

		print "Populating background locations"
		populate_locations(fg_mers.keys(), bg_mers, args.background)

		for mer in bg_mers:
			bg_mers[mer] = len(bg_mers[mer])

		score_specific_combinations(combinations)

	print "output_file:", output_file

if __name__ == "__main__":
	sys.exit(main())
