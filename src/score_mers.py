#!/usr/bin/env python
import sys
import os

from multiprocessing import Pool
from multiprocessing import cpu_count
from subprocess import *
from itertools  import combinations
from itertools  import ifilter
from itertools  import imap

import numpy as np
import pdb

fg_mers = {}
bg_mers = {}

seq_ends = []

if len(sys.argv) == 5:
	selectivity_fn =  sys.argv[1]
	fg_fasta_fn    =  sys.argv[2]
	bg_fasta_fn    =  sys.argv[3]
	output_file    =  sys.argv[4]

	fg_genome_length = os.path.getsize(fg_fasta_fn)
	bg_genome_length = os.path.getsize(bg_fasta_fn)
else:
	print "please specify your inputs"
	print "ex: score_mers.py selectivity_file fg_fasta bg_fasta output_file"
	exit()

# import our variables
cpus             = int(os.environ.get("cpus", cpu_count()))
debug            = int(os.environ.get("debug", False))
min_mer_range    = int(os.environ.get("min_mer_range", 6))
max_mer_range    = int(os.environ.get("max_mer_range", 12))
min_mer_count    = int(os.environ.get("min_mer_count", 0))
max_select       = int(os.environ.get("max_select", 15))
max_check        = int(os.environ.get("max_check", 35))
max_mer_distance = int(os.environ.get("max_mer_distance", 5000))
max_consecutive_binding = int(os.environ.get("max_consecutive_binding", 4))


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
	cmds.append("tac " + input_fn + " | rev | grep -v '>$' | tr -d '\\n' | tr [ACGT] [TGCA] | strstream ")
	
	for cmd in cmds:
		fid, merlist_fn = tempfile.mkstemp()

		# write our mers out to a fifi
		merlist_fh = open(merlist_fn, 'w')
		for mer in selected_mers:
			print mer
			merlist_fh.write(mer + '\n')

		merlist_fh.flush()
		# add our merlist fn to our command
		cmd = cmd + " " + merlist_fn

		strstream = Popen(cmd, stdout=PIPE, shell=True)
		for line in strstream.stdout:
			(mer, pos) = line.strip().split(" ")
			mer_dic[selected_mers[int(mer)]].append(int(pos))

		merlist_fh.close()
	
	

def filter_mers(combination):
	for combo in combinations(combination, 2):
		if heterodimer_dic[combo]:
			return True

	for mer in combination:
		for other_mer in combination:
			if not mer == other_mer:
				if mer in other_mer:
					return True

	return False

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

def score_mers(selected):
	import time
	total_scored = 0

	check_feasible(selected)

	p = Pool(cpus)

	fh = open(output_file, 'wb')
	fh.write("Combination\tScore\tFG_mean_dist\tFG_stdev_dist\tBG_mean_dist\tBG_var_dist\n")

	for select_n in range(1, max_select+1):
		print "scoring size ", select_n,
		t = time.time()
		scores_it = p.imap_unordered(score, combinations(selected, select_n), chunksize=8192)
		for score_res in scores_it:
			if score_res is not None:
				total_scored += 1
				combination, score_val, fg_mean_dist, fg_stddev_dist, bg_ratio = score_res
				fh.write(str(combination) + "\t")
				fh.write(str(score_val) + "\t")
				fh.write(str(fg_mean_dist) + "\t")
				fh.write(str(fg_stddev_dist) + "\t")
				fh.write(str(bg_ratio) + "\n")
		print "size ", select_n, "took:", time.time()   - t

	if(total_scored == 0):
		print "NO RESULTS FOUND"
		fh.write("NO RESULTS FOUND\n")

heterodimer_dic = {}
def score(combination):
# input is a string of mers like 
# ['ACCAA', 'ACCCGA', 'ACGTATA']

	# check if the combination passes our filters
	if filter_mers(combination):
		return None

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
		 return None

	# bg counts 
	bg_counts = 0

	for mer in combination:
		bg_counts += len(bg_mers[mer])

	if bg_counts <= 1:
		bg_counts = 1 

	bg_sum = len(bg_counts)
	bg_ratio = (bg_genome_length / bg_sum)


	nb_primers = len(combination)
	fg_mean_dist =  np.mean(fg_dist)
	fg_std_dist = np.std(fg_dist)

	# this is our equation
	mer_score = (nb_primers * fg_mean_dist * fg_std_dist) / bg_ratio

	return [combination, mer_score, fg_mean_dist, fg_std_dist, bg_ratio] 

def load_end_points(fn):

	end_points = [0]

	cmd = "sequence_end_points < " + fn

	if debug:
		print "loading sequence end points"
		print "executing: " + cmd

	points_fh = Popen(cmd, stdout=PIPE, shell=True)

	for line in points_fh.stdout:
		end_points.append(int(line))
	
	return end_points

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
	import time
	selected = []
	selectivity_fh = open(selectivity_fn, "r")
	
	# load our mer list into python
	mer_selectivity = selectivity_fh.readlines()

	# get the last max_check (it's sorted)
	selected_mers = mer_selectivity[-max_check:]
	selected_mers = [x.split()[0] for x in selected_mers]

	# load it into our fg and bg dictionary points
	for mer in selected_mers:
		fg_mers[mer] = []
		bg_mers[mer] = []

	print "Populating sequence end points"
	seq_ends = load_end_points(fg_fasta_fn)

	print "Populating foreground locations"
	populate_locations(selected_mers, fg_mers, fg_fasta_fn)

	print "Populating background locations"
	populate_locations(selected_mers, bg_mers, bg_fasta_fn)

	print "calculating heterodimer distances"
	load_heterodimer_dic(selected_mers)

	print "scoring mer combinations"
	score_mers(selected_mers)

	print "output_file:", output_file
	

if __name__ == "__main__":
	sys.exit(main())
