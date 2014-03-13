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

if(len(sys.argv) == 5):
	selectivity_fn =  sys.argv[1]
	fg_fasta_fn    =  sys.argv[2]
	bg_fasta_fn    =  sys.argv[3]
	output_file    =  sys.argv[4]

	fg_genome_length = os.path.getsize(fg_fasta_fn)
	bg_genome_length = os.path.getsize(bg_fasta_fn)
else:
	print "please specify your inputs"
	print "ex: score_mers.py selectivity_file fg_fasta_file bg_fasta_file"
	exit()

# empty class to fill up mer information with
class Mer:
	pass

# import our variables
cpus = int(os.environ.get("cpus", cpu_count()));
min_mer_range    = int(os.environ.get("min_mer_range", 6));
max_mer_range    = int(os.environ.get("max_mer_range", 12));
min_mer_count    = int(os.environ.get("min_mer_count", 0));
max_select       = int(os.environ.get("max_select", 15));
max_check        = int(os.environ.get("max_check", 35));
max_mer_distance = int(os.environ.get("max_mer_distance", 5000));
max_consecutive_binding = int(os.environ.get("max_consecutive_binding", 4));


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

	max_bind = 0;
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

def pop_fg(mer):
	''' helper for map function '''
	populate_locations(fg_fasta_fn, fg_mers, mer)

def pop_bg(mer):
	''' helper for map function '''
	populate_locations(bg_fasta_fn, bg_mers, mer)


def populate_locations(input_fn, mers, mer):
	''' 
	Run the strstreamone command, and parse in the integers that are output
	by the command, and add it to mers[mer].pts 

	strstreamone just prints the location of a string argv[1] in stdout.

	We also do the reverse compliment, using tac and tr piped together.
	'''

	cmd = 'strstreamone ' + mer + " < " + input_fn
	
	strstream = Popen(cmd, stdout=PIPE, shell=True)
	for line in strstream.stdout:
		mers[mer].pts.append(int(line))

	cmd = 'tac ' + input_fn + " | tr '[ACGT]' '[TGCA]' | strstreamone " + mer + " " + input_fn
	strstream = Popen(cmd, stdout=PIPE, shell=True)
	for line in strstream.stdout:
		mers[mer].pts.append(int(line))


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
	total = 0;
	for mer in selected:
		total += len(fg_mers[mer].pts)
	if (fg_genome_length / total) > max_mer_distance:
		print "even if we select all top ", max_select, 
		print "mers disregarding any critera, and they were perfectly evenly spaced we would ",
		print "still not meet the right max mer distance < ", max_mer_distance, "requirement."
	
		print total, " / ", fg_genome_length, " = ", total / fg_genome_length 
		exit()

def score_mers(selected):
	import time
	total_scored = 0;

	check_feasible(selected)

	p = Pool(cpus)

	fh = open(output_file, 'wb');
	fh.write("Combination\tScore\tFG_mean_dist\tFG_var_dist\tBG_mean_dist\tBG_var_dist\n");
	for select_n in range(1, max_select+1):
		print "scoring size ", select_n,
		t = time.time()
		scores_it = p.imap_unordered(score, combinations(selected, select_n), chunksize=8192)
		for score_res in scores_it:
			if score_res is not None:
				total_scored += 1;
				combination, scores, fg_mean_dist, fg_variance_dist, bg_mean_dist, bg_variance_dist = score_res
				fh.write(str(combination) + "\t");
				fh.write(str(scores) + "\t");
				fh.write(str(fg_mean_dist) + "\t");
				fh.write(str(fg_variance_dist) + "\t");
				fh.write(str(bg_mean_dist) + "\t");
				fh.write(str(bg_variance_dist) + "\n");
		print "size ", select_n, "took:", time.time()   - t

	if(total_scored == 0):
		print "NO RESULTS FOUND"

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
		fg_pts = fg_pts + fg_mers[mer].pts

	fg_pts.sort()

	# fg distances
	fg_dist = np.diff(fg_pts)

  # return without calculating scores if any objects are higher than our max distance
	if any(dist > max_mer_distance for dist in fg_dist):
		#return [combination, "max", max(fg_dist)]
		 return None

	min_mer_distance = max(len(i) for i in combination)
	# return without calculating scores if any mers are closer than the length of
	# our longest mer in the combination
	if any(dist < min_mer_distance for dist in fg_dist):
		#return [combintaion, 'max']
		return None


	# bg points
	bg_pts = []
	bg_dist = []

	for mer in combination:
		bg_pts = bg_pts + bg_mers[mer].pts

	if len(bg_pts()) <= 0:
		bg_pts.append(0, 1, fg_genome_length)

	bg_pts.sort()

	# bg distances
	bg_dist = np.diff(bg_pts)

	nb_primers = len(combination)
	fg_mean_dist =  np.mean(fg_dist)
	fg_variance_dist = np.var(fg_dist)
	bg_mean_dist = np.mean(bg_dist)
	bg_variance_dist = np.var(bg_dist)

	# this is our equation
	score = (nb_primers * fg_mean_dist * fg_variance_dist) / ((bg_mean_dist * bg_variance_dist) + .000001)

	return [combination, score, fg_mean_dist, fg_variance_dist, bg_mean_dist, bg_variance_dist]

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
	
	# get our genome length

	for row in selectivity_fh:
		(mer, fg_count, bg_count, selectivity) = row.split()
		fg_mers[mer] = Mer()
		fg_mers[mer].pts = []
		fg_mers[mer].count = fg_count
		bg_mers[mer] = Mer()
		bg_mers[mer].pts = []
		bg_mers[mer].count = bg_count
		selected.append([mer, selectivity])
		
	selected =	selected[-max_check:] 
	selected_mers = [row[0] for row in selected]
	# print selected_mers

	print "Populating foreground locations"
	map(pop_fg, selected_mers)


	print "Populating background locations"
	map(pop_bg, selected_mers)

	print "calculating heterodimer distances"
	load_heterodimer_dic(selected_mers)

	print "scoring mer combinations"
	score_mers(selected_mers)

	print "output_file:", output_file
	

if __name__ == "__main__":
	sys.exit(main())
