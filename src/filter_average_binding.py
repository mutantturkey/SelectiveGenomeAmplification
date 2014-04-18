#!/usr/bin/env python2.7
import sys
import os

debug = os.environ.get("debug", False)

from subprocess import Popen
from subprocess import PIPE

def get_length(fn):

	cmd = 'grep "^>" ' + fn + " -v | tr -d '\\n' | wc -c"

	if debug:
		print "loading sequence end points"
		print "executing: " + cmd

	points_fh = Popen(cmd, stdout=PIPE, shell=True)

	return int(points_fh.stdout.readline())

if len(sys.argv) < 2:
	print "filter_average_binding.py foreground.fa binding_distance"
	exit()
	
foreground = sys.argv[1]
distance = int(sys.argv[2])

genome_length = get_length(foreground)

for line in sys.stdin:
	(mer, count) = line.split()
	if (genome_length / int(count)) < distance:
		sys.stdout.write(line)

