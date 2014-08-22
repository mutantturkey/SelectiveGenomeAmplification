#!/usr/bin/env python2.7
import sys
import argparse
import os

debug = os.environ.get("debug", False)

def get_length(fn):
	'''
		get the length of a fasta file by piping it through several unix 
		programs.

		1) remove headers by grepping for any ">" at the start of a line
		2) delete all occurances of a new line, to join sequences together
		3) sum the number of characters.
	'''
	from subprocess import Popen
	from subprocess import PIPE

	cmd = 'grep "^>" ' + fn + " -v | tr -d '\\n' | wc -c"

	if debug:
		print "loading sequence end points"
		print "executing: " + cmd

	points_fh = Popen(cmd, stdout=PIPE, shell=True)

	return int(points_fh.stdout.readline())

def main(): 
	'''
	This filter removes mers where the count / (length of the genome) is below a
	certain threshold as specified by -m
	'''

	parser = argparse.ArgumentParser(description="Filter mers where (k-mer count / length of the genome) < minimum") 
	parser.add_argument("-f", "--fasta", help="foreground fasta file", required=True )
	parser.add_argument("-c", "--counts", help="kmer counts of the foreground fasta file", required=True )
	parser.add_argument("-m", "--minimum", help="the minium average foreground binding distance", required=True, type=float)

	args = parser.parse_args()

	if not os.path.exists(args.fasta):
		exit("foreground fasta file " + args.fasta + " not found.")

	if not os.path.exists(args.counts):
		exit("count file " + args.counts + " not found.")
			

	# get genome length
	genome_length = float(get_length(args.fasta))

	count_fh = open(args.counts, "rU")

	for line in count_fh:
		(_, count) = line.split()
		if (genome_length / float(count)) < args.minimum:
			sys.stdout.write(line)

if __name__ == "__main__":
	sys.exit(main())

