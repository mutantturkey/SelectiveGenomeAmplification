#!/usr/bin/env python2.7
import sys
import os

fg_weight = float(os.environ.get("fg_weight", 0))
max_check = int(os.environ.get("max_check", 0))


def main():

	if(len(sys.argv) == 3):
		fg_count_fn =  sys.argv[1]
		bg_count_fn =  sys.argv[2]
	else:
		sys.stderr.write("please specify your inputs\n")
		sys.stderr.write("ex: select_mers.py fg_counts bg_count\n")
		exit(1)

	# mers dictionary:
	# 
	# Key: mer name, eg AAAACT 
	# Value: fg_mer_count, bg_mer_count
	mers = {}

	fg_count_fh = open(fg_count_fn, "r")
	bg_count_fh = open(bg_count_fn, "r")
	
	# copy in our foreground mers and counts into mers dictionary
	for line in fg_count_fh:
		(mer, count) = line.split()
		mers[mer] = [float(count), 1]
	
	
	for line in bg_count_fh:
		(mer, count) = line.split()
		if mer in mers:
			mers[mer][1] = float(count)

	score = []

	for mer in mers:
		score.append([mer, (mers[mer][0] / mers[mer][1]) * (mers[mer][0]**fg_weight)])

	sorted_scored_mers = sorted(score, key=lambda x: x[1])

	sys.stdout.write('#MERS\tFG_COUNT\tBG_COUNT\tSCORE\n')
	for scores in sorted_scored_mers: 
		sys.stdout.write(scores[0] + '\t' + str(int(mers[scores[0]][0])) + '\t' +  str(int(mers[scores[0]][1])) + '\t' +  str(scores[1])+ '\n')

if __name__ == "__main__":
	sys.exit(main())
