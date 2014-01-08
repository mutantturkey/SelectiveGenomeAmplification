import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import os

from multiprocessing import Pool
from subprocess import *
import numpy as np
import pdb

fg_mers = {}
bg_mers = {}
# empty class to fill up mer information with
class Mer:
	pass

class Score:
	pass

# import our variables
min_mer_range = int(os.getenv("min_mer_range"));
max_mer_range = int(os.getenv("max_mer_range"));
min_mer_count = int(os.getenv("min_mer_count"));
max_select    = int(os.getenv("max_select"));

def populate_locations(input_fn, mers, selected):

	mer_str = ' '.join(selected)
	cmd = './strstream ' + mer_str + " < " + input_fn
	
	strstream = Popen(cmd, stdout=PIPE, shell=True)
	for line in strstream.stdout:
		(index, pos) = line.split(' ')
		mers[selected[int(index)]].pts.append(int(pos))


	return None

def select_mers(selected):
	from itertools  import combinations

	scores = []

	p = Pool()

	for select_n in range(1, max_select):
		print "scoring size ", select_n
		scores_it = []
	 	scores_it = p.map(score, combinations(selected, select_n))
	 	scores = scores + scores_it


	p.close()


	scores = sorted(scores, key=lambda score: score.score)
	return scores

def score(combination):
	
  #	 print "scoring", combination

	fg_pts_arr = []
	fg_dist_arr = []

	bg_pts_arr = []
	bg_dist_arr = []

	for mer in combination:
		fg_pts_arr = fg_pts_arr + fg_mers[mer].pts
		bg_pts_arr = bg_pts_arr + bg_mers[mer].pts

	fg_pts_arr.sort()
	bg_pts_arr.sort()

	# remove any exact duplicates
	fg_pts_arr = list(set(fg_pts_arr))
	bg_pts_arr = list(set(bg_pts_arr))

	#	distances
	for i in range(1, len(fg_pts_arr)):
		fg_dist_arr.append(abs(fg_pts_arr[i-1] - fg_pts_arr[i]));

	for i in range(1, len(bg_pts_arr)):
		bg_dist_arr.append(abs(bg_pts_arr[i-1] - bg_pts_arr[i]));

	fg_dist_arr = np.array(fg_dist_arr)
	bg_dist_arr = np.array(bg_dist_arr)

	nb_primers = len(combination)
	fg_mean_dist =  np.mean(fg_dist_arr)
	fg_variance_dist = np.var(fg_dist_arr)
	bg_mean_dist = np.mean(bg_dist_arr)
	bg_variance_dist = np.var(bg_dist_arr)

	# this is our equation
	
	score = Score()
	score.mers = combination
	score.score = (nb_primers * (fg_mean_dist * fg_variance_dist)) / (bg_mean_dist * bg_variance_dist)
	# score.bg_dist_arr = bg_dist_arr
	# score.fg_dist_arr = fg_dist_arr
	return score
				


def most_frequent(mers, select_nb):
		print "selecting the top ", select_nb, " mers"
		selected = []
		for mer in mers.keys():
			selected.append([mer, mers[mer].count])
			
		selected = sorted(selected, key=lambda mer: mer[1])
		selected = selected[-select_nb:]
		selected = list(row[0] for row in selected)

		return selected


def main():
	fg_count_fn =  sys.argv[1]
	fg_fasta_fn =  sys.argv[2]

	bg_count_fn =  sys.argv[3]
	bg_fasta_fn =  sys.argv[4]

	fg_count_fh = open(fg_count_fn, "r")
	bg_count_fh = open(bg_count_fn, "r")
	
	# get our genome length
	fg_genome_length = os.path.getsize(fg_fasta_fn)
	bg_genome_length = os.path.getsize(bg_fasta_fn)

	# copy in our fg_mers and counts
	for mers,fh in [(fg_mers, fg_count_fh), (bg_mers, bg_count_fh)]:

		for line in fh:
			(mer, count) = line.split()
			count = int(count)
			mers[mer] = Mer() 
			mers[mer].count = count
			mers[mer].pts = []
	
	# if min_mer_count is less than one, then we are looking at a bottom cutoff threshold
	if(min_mer_count < 1 and min_mer_count > 0):
		pass
		#	counts = []
		#	for mer in fg_mers.keys():
		#		counts.append(fg_mers[mer]['count'])
		#
		#	counts = counts.sort()
		#	
		#	count_set = set(counts):
		#	len(count_set)
		#	count_set(
	# if it is 1 or great then use it as a numerical cutoff
	print "removing mers that don't meet our minimum count"
	if min_mer_count >= 1:
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

		
	exhaustive = True

	if exhaustive:
		selected = fg_mers.keys()
	if not exhaustive:
		selected = most_frequent(fg_mers, max_select)

	print "selected the top ", max_select, " mers"
	print "selected:", ", ".join(selected)

	print "Populating foreground locations"
	populate_locations(fg_fasta_fn, fg_mers, selected)
	print "Populating background locations"
	populate_locations(bg_fasta_fn, bg_mers, selected)

	scores = select_mers(selected)

	print "fg_genome_length", fg_genome_length
	print "bg_genome_length", bg_genome_length
	

	for score in scores:
		sys.stdout.write(str(score.score) + "\t")
		sys.stdout.write(", ".join(list(score.mers)) + "\t")
		# sys.stdout.write("\t".join(score.dist_arr));
		sys.stdout.write("\n");

		#matplotlib.pyplot.hist(score.fg_dist_arr)	
		#matplotlib.pyplot.savefig('graph/' + str(x) + ".png")
		#matplotlib.pyplot.clf()
	import pdb
	pdb.set_trace()

if __name__ == "__main__":
	sys.exit(main())
