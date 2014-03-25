#!/usr/bin/env python
import sys
import os

fg_mers = {}
bg_mers = {}

if(len(sys.argv) == 3):
	fg_count_fn =  sys.argv[1]
	bg_count_fn =  sys.argv[2]
else:
	print len(sys.argv)
	print "please specify your inputs"
	print "ex: select_mers.py fg_counts bg_counts"
	exit()


# select mers based on our 'selectivity' measure. (count in fg) / (count in bg)
def select_mers(fg_mers, bg_mers):
	import numpy as np

	mers   = [] # contains mer strings
	fg_arr = [] # contains fg counts
	bg_arr = [] # contains bg counts
	
	# populate our bg_arr and fg_arr as well as our mer arr.
	for mer in fg_mers.keys():
		mers.append(mer);
		bg_arr.append(bg_mers.get(mer, 1));
		fg_arr.append(fg_mers[mer]);

	fg_arr = np.array(fg_arr, dtype='f');
	bg_arr = np.array(bg_arr, dtype='f');

	selectivity = (fg_arr / bg_arr)

	arr = [(mers[i], fg_arr[i], bg_arr[i], selectivity[i]) for i in range(len(mers))]

	# filter results less than 1 ( indicates that the bg is more present than the fg)
	# arr = filter(lambda i: i[3] > 1, arr)

	# sort by the selectivity 
	arr = sorted(arr, key = lambda row: row[3])

	# return only our mers, without our selectivity scores
	return arr


def main():

	fg_count_fh = open(fg_count_fn, "r")
	bg_count_fh = open(bg_count_fn, "r")
	
	# copy in our fg_mers and counts
	for mers,fh in [(fg_mers, fg_count_fh), (bg_mers, bg_count_fh)]:
		for line in fh:
			(mer, count) = line.split()
			mers[mer] = int(count)
	
	for mer in bg_mers.keys():
		if mer not in fg_mers:
			del bg_mers[mer]

	selected = select_mers(fg_mers, bg_mers)
	for row in selected:
		print row[0] +"\t"+str("%d" % row[1]) + "\t" + str("%d" % row[2]) + "\t" + str("%.5f" % row[3])

if __name__ == "__main__":
	sys.exit(main())
