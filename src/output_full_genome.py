#!/usr/bin/env python
import sys
import os
import argparse

mers = {} 
debug = False

from subprocess import Popen
from subprocess import PIPE


def get_sequence(pt):
	
	for it, seq in enumerate(seq_ends):
		if pt <= seq:
			return it

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

def populate_locations(selected_mers, mer_dic, input_fn, length):
	''' 
	Run the strstreamone command, and parse in the integers that are output
	by the command, and add it to mers[mer] 

	strstreamone just prints the location of a string argv[1] in stdout.

	We also do the reverse compliment, using tac and tr piped together.
	'''
	import tempfile


	cmds = []
	# strip file of header and delete newlines
	cmds.append(["grep -v '^>' " + input_fn  +  " | tr -d '\\n' | strstream ", False])
	# reverse file, strip and delete newlines
	cmds.append(["tac " + input_fn + \
							"| rev " \
							"| grep -v '>$' " \
							"| tr -d '\\n' " \
							"| tr [ACGT] [TGCA] | strstream ", True])
	
	for (cmd, reverse) in cmds:
		_, merlist_fn = tempfile.mkstemp()

		# write our mers out to a fifi
		merlist_fh = open(merlist_fn, 'w')
		for mer in selected_mers:
			merlist_fh.write(mer + '\n')

		merlist_fh.flush()
		# add our merlist fn to our command
		cmd = cmd + " " + merlist_fn

		strstream = Popen(cmd, stdout=PIPE, shell=True)
		if reverse:
			for line in strstream.stdout:
				(mer, pos) = line.split(" ")
				pos = length - int(pos)
				mer_dic[selected_mers[int(mer)]].append(pos)
		else:
			for line in strstream.stdout:
				(mer, pos) = line.split(" ")
				mer_dic[selected_mers[int(mer)]].append(int(pos))

		if strstream.wait() is not 0:
			print "executing", cmd, "failed"

		merlist_fh.close()


def main():

	parser = argparse.ArgumentParser(description="score mers")
	parser.add_argument("-f", "--fasta", help="input fasta file", required=True)
	parser.add_argument("-s", "--scores", help="scores file", required=True)
	parser.add_argument("-o", "--output_directory", help="scores file", required=True)
	parser.add_argument("-n", "--number", help="output top n", required=False, type=int, default=20)

	args = parser.parse_args()

	if os.path.isfile(args.output_directory):
		parser.error(args.output_directory + "must point to a directory")
	elif not os.path.isdir(args.output_directory):
		os.mkdir(args.output_directory)

	score_fh = open(args.scores, "r")

	global seq_ends
	seq_ends = load_end_points(args.fasta)

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
			populate_locations(new_populate, mers, args.fasta, get_length(args.fasta))

		pts = []
		for mer in combination:
			for pt in mers[mer]:
				pts.append([pt, mer, get_sequence(pt)])
		
		pts = sorted(pts, key = lambda row: row[0])
			
		for pt in pts:
			fh.write(str(pt))
			fh.write("\n")

		nb_done += 1 
		if nb_done == args.number:
			exit()

if __name__ == "__main__":
	sys.exit(main())
