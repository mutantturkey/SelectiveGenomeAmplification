#!/usr/bin/env python
import sys
import os

from multiprocessing import Pool
from subprocess import *
import numpy as np
import pdb

fg_mers = {}
bg_mers = {}

if(len(sys.argv) != 5):
	fg_count_fn =  sys.argv[1]
	fg_fasta_fn =  sys.argv[2]
	bg_count_fn =  sys.argv[3]
	bg_fasta_fn =  sys.argv[4]
	output_file =  sys.argv[5]
else:
	print "please specify your inputs"
	print "ex: select_mers.py fg_counts_file fg_fasta_file bg_counts_file bg_fasta_file output_file"
	exit()

# empty class to fill up mer information with
class Mer:
	pass

# import our variables
min_mer_range    = int(os.environ.get("min_mer_range", 6));
max_mer_range    = int(os.environ.get("max_mer_range", 10));
min_mer_count    = int(os.environ.get("min_mer_count", 0));
max_select       = int(os.environ.get("max_select", 15));
max_mer_distance = int(os.environ.get("max_mer_distance", 5000));


def populate_locations(input_fn, mers, mer):
	''' Run the strstreamone command, and parse in the integers that are output
			by the command, and add it to mers[mer].pts 
	'''

	cmd = 'strstreamone ' + mer + " < " + input_fn
	
	strstream = Popen(cmd, stdout=PIPE, shell=True)
	for line in strstream.stdout:
		mers[mer].pts.append(int(line))


def score_mers(selected):
	from itertools  import combinations
	import time

	scores = []

	p = Pool()

	fh = open(output_file, 'w');
	fh.write("scores:\n");
	for select_n in range(1, max_select+1):
		print "scoring size ", select_n,
		t = time.time()
		scores_it = p.imap_unordered(score, combinations(selected, select_n))
		for score_res in scores_it:
			fh.write(str(score_res) + "\n");
	 print "size ", select_n, "took:", t - time.time()	
	return scores


def score(combination):
# input is a string of mers like 
# ['ACCAA', 'ACCCGA', 'ACGTATA']

	ret = [combination]

	for mer in combination:
		for other_mer in combination:
			if not mer == other_mer:
				if mer in other_mer:
					ret.append("duplicates")
					return ret 

	fg_pts = []
	fg_dist = []

	bg_pts = []
	bg_dist = []

	for mer in combination:
		fg_pts = fg_pts + fg_mers[mer].pts
		bg_pts = bg_pts + bg_mers[mer].pts

	fg_pts.sort()
	bg_pts.sort()

	# remove any exact duplicates
	# fg_pts = list(set(fg_pts))
	# bg_pts = list(set(bg_pts))

	#	distances
	min_mer_distance = max(len(i) for i in combination)
	fg_dist = np.array([abs(fg_pts[i] - fg_pts[i-1]) for i in range(1, len(fg_pts))])
	bg_dist = np.array([abs(bg_pts[i] - bg_pts[i-1]) for i in range(1, len(bg_pts))])

  # return without calculating scores if any objects are higher than our max distance
	if any(dist > max_mer_distance for dist in fg_dist):
		ret.append("max")
		ret.append(max(fg_dist))
		return ret 

	# return without calculating scores if any mers are closer than the length of our longest mer in the combination
	if any(dist < min_mer_distance for dist in fg_dist):
		ret.append("min")
		ret.append(min(fg_dist))
		return ret 

	nb_primers = len(combination)
	fg_mean_dist =  np.mean(fg_dist)
	fg_variance_dist = np.var(fg_dist)
	bg_mean_dist = np.mean(bg_dist)
	bg_variance_dist = np.var(bg_dist)

	# this is our equation
	score = (nb_primers * fg_mean_dist * fg_variance_dist) / ((bg_mean_dist * bg_variance_dist) + .000001)

	ret.append(score)
	ret.append(fg_mean_dist)
	ret.append(fg_variance_dist)
	ret.append(bg_mean_dist)
	ret.append(bg_variance_dist)

	return ret


# select mers based on our 'selectivity' measure. (count in fg) / (count in bg)
def select_mers(fg_mers, bg_mers, select_nb):

	mers   = [] # contains mer strings
	fg_arr = [] # contains fg counts
	bg_arr = [] # contains bg counts

	# populate our bg_arr and fg_arr as well as our mer arr.
	for mer in fg_mers.keys():
		mers.append(mer);
		bg_arr.append(bg_mers[mer].count);
		fg_arr.append(fg_mers[mer].count);

	fg_arr = np.array(fg_arr, dtype='f');
	bg_arr = np.array(bg_arr, dtype='f');

	selectivity = fg_arr / bg_arr

	arr = [(mers[i], fg_arr[i], bg_arr[i], selectivity[i]) for i in range(len(mers))]

	# filter results less than 1 ( indicates that the bg is more present than the fg)
	# arr = filter(lambda i: i[3] > 1, arr)

	# sort by the selectivity 
	arr = sorted(arr, key = lambda row: row[3])

	# return only our mers, without our selectivity scores
	arr = list(row[0] for row in arr)
	return arr

def pop_fg(mer):
	''' helper for map function '''
	populate_locations(fg_fasta_fn, fg_mers, mer)

def pop_bg(mer):
	''' helper for map function '''
	populate_locations(bg_fasta_fn, bg_mers, mer)

def main():
	import time
	fg_count_fh = open(fg_count_fn, "r")
	bg_count_fh = open(bg_count_fn, "r")
	
	# get our genome length
	fg_genome_length = os.path.getsize(fg_fasta_fn)
	bg_genome_length = os.path.getsize(bg_fasta_fn)

	# copy in our fg_mers and counts
	for mers,fh in [(fg_mers, fg_count_fh), (bg_mers, bg_count_fh)]:

		t = time.time()
		for line in fh:
			(mer, count) = line.split()
			mers[mer] = Mer() 
			mers[mer].count = int(count)
			mers[mer].pts = []
		print time.time() - t
	
	if min_mer_count >= 1:
		print "removing that are less frequent than: ", min_mer_count
		for mer in fg_mers.keys():
			if(fg_mers[mer].count < min_mer_count):
				del fg_mers[mer]
				if mer in bg_mers:
					del bg_mers[mer]
	
	print "removing useless mers from the background"
	for mer in bg_mers.keys():
		if mer not in fg_mers:
			del bg_mers[mer]

	print "adding empty mers to the background"
	for mer in fg_mers:
		if mer not in bg_mers:
			bg_mers[mer] = Mer()
			bg_mers[mer].count = 2
			bg_mers[mer].pts = [0, bg_genome_length]

		
	exhaustive = False

	if exhaustive:
		selected = fg_mers.keys()
	else:
		selected = select_mers(fg_mers, bg_mers, max_select)
		selected =	selected[-100:] 
		print "searching through combinations of"
		print selected

#	pdb.set_trace()

	# print "selected the top ", max_select, " mers"
	# print "selected:", ", ".join(selected)

	print "Populating foreground locations"

	# p = Pool()

	map(pop_fg, selected)
	map(pop_bg, selected)

	scores = score_mers(selected)

	print "fg_genome_length", fg_genome_length
	print "bg_genome_length", bg_genome_length
	print "output_file:", output_file
	



if __name__ == "__main__":
	sys.exit(main())
