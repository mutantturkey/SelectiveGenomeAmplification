#!/usr/bin/env python2.7
import sys
import os

fg_mers = {}
bg_mers = {}

fg_weight = float(os.environ.get("fg_weight", 0))
max_check = int(os.environ.get("max_check", 0))

if(len(sys.argv) == 3):
	fg_count_fn =  sys.argv[1]
	bg_count_fn =  sys.argv[2]
else:
	print len(sys.argv)
	sys.stderr.write("please specify your inputs\n")
	sys.stderr.write("ex: select_mers.py fg_counts bg_count\n")
	exit(1)


# select mers based on our 'selectivity' measure. (count in fg) / (count in bg)
def select_mers(fg_mers, bg_mers):

	# populate our bg_arr and fg_arr as well as our mer arr.

	score = {}

	for mer in fg_mers.keys():
		score[mer] = (fg_mers[mer] / bg_mers[mer]) * (fg_mers[mer]**fg_weight)

	sorted_scored_mers =	sorted(score, key=score.get)

	for mer in sorted_scored_mers: 
		print mer, int(fg_mers[mer]), int(bg_mers[mer]), (fg_mers[mer] / bg_mers[mer]) * (fg_mers[mer]**fg_weight)


def main():

	fg_count_fh = open(fg_count_fn, "r")
	bg_count_fh = open(bg_count_fn, "r")
	
	# copy in our fg_mers and counts
	for mers,fh in [(fg_mers, fg_count_fh), (bg_mers, bg_count_fh)]:
		for line in fh:
			(mer, count) = line.split()
			mers[mer] = float(count)
	
	for mer in fg_mers.keys():
		if mer not in bg_mers:
			bg_mers[mer] = 1

	for mer in bg_mers.keys():
		if mer not in fg_mers:
			del bg_mers[mer]

	selected = select_mers(fg_mers, bg_mers)

if __name__ == "__main__":
	sys.exit(main())
