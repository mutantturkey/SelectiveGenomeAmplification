#!/usr/bin/env python2.7
import sys
import os
import argparse

mers = {} 
debug = os.environ.get("debug", False)

from subprocess import Popen
from subprocess import PIPE


def get_sequence(pt):
	
	for it, seq in enumerate(seq_ends, start=1):
		if pt <= seq:
			return it

def load_end_points(fn):
	''' get all the points of the end of each sequence in a sample '''

	end_points = []

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


	return length


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
			print mer
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
				mer_dic[line_s[0]].append([int(line_s[1]), 3])
			elif len(line_s) == 3:
				mer_dic[line_s[0]].append([length - int(line_s[1]), 5])

		if cmd_pipe.wait() is not 0:
			print "Executing", cmd, "failed"
			exit()

		merlist_fh.close()
def main():

	parser = argparse.ArgumentParser(description="take a top-scores file, and \
		create a set of files that contain the positions (and other information) in \
		the foreground fasta file of each mer in that combination. ") 
	parser.add_argument("-f", "--fasta", help="input fasta file", required=True)
	parser.add_argument("-s", "--scores", help="input scores file", required=True)
	parser.add_argument("-o", "--output_directory", help="directory to output sets to.", required=True)
	parser.add_argument("-n", "--number", help="output the first n sets", required=False, type=int, default=20)

	args = parser.parse_args()

	if os.path.isfile(args.output_directory):
		parser.error(args.output_directory + "must point to a directory")
	elif not os.path.isdir(args.output_directory):
		os.mkdir(args.output_directory)
	score_fh = open(args.scores, "r")

	global seq_ends
	seq_ends = load_end_points(args.fasta)

	length = get_length(args.fasta)

	nb_done = 0;
	for line in score_fh:
		# skip headers
		if line.startswith("#"):
			continue

		# get mers
		combination = line.split('\t')[1].split()
		
		# open output
		fh = open(args.output_directory + "/" + '-'.join(combination), 'w')

		# populate only un-populated mers
		new_populate = []
		for mer in combination:
			if mer not in mers:
				mers[mer] = []
				new_populate.append(mer)
		
		if len(new_populate) is not 0:
			populate_locations(new_populate, mers, args.fasta, length)

		pts = []
		for mer in combination:
			for pt in mers[mer]:
				pts.append([pt[0], pt[1], mer, get_sequence(pt[0])])
		
		pts = sorted(pts, key = lambda row: row[0])
			
		fh.write("pt\tstrand\tmer\tsequence\n")
		for pt in pts:
			fh.write('\t'.join(str(x) for x in pt) + '\n')

		nb_done += 1 
		if nb_done == args.number:
			exit()

if __name__ == "__main__":
	sys.exit(main())
