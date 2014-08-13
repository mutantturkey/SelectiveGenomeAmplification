#!/usr/bin/env python2.7
import sys, os
	
def main():

	if(len(sys.argv) < 2):
		print "cutoff and bg_counts is expected as an argument"
		exit()
	else:
		cutoff = int(sys.argv[1])
		bg_count_fn = sys.argv[2]
	
	# if  cutoff, is less than zero, we ignore, aka so we can do -1 by default,
	# we can't do 0, because that might have a valid use case
	if cutoff < 0:
		for line in sys.stdin:
			sys.stdout.write(line)
	else:

		mers = {}

		bg_count_fh = open(bg_count_fn, "r")
		
		# copy in our foreground mers and counts into mers dictionary, then process it
		for line in sys.stdin:
			(mer, count) = line.split()
			mers[mer] = [int(count), -1]
		
		for line in bg_count_fh:
			(mer, count) = line.split()
			if mer in mers:
				mers[mer][1] = int(count)

		for mer in mers:
			if mers[mer][1] == -1 or mers[mer][1] <= cutoff:
				sys.stdout.write(mer + '\t' + str(mers[mer][0]) + '\n')

if __name__ == "__main__":
	sys.exit(main())
