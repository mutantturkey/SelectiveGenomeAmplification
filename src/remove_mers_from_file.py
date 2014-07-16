#!/usr/bin/env python2.7
import sys

if __name__ == "__main__":

	if len(sys.argv) is 1:
		exit("Filter mers, input is stdin, output is stdout, mers are argv")
	
	mers_to_delete = set()

	for mer in open(sys.argv[1], 'r'):
		mer = mer.strip().split()[0]
		mers_to_delete.add(mer)

	for line in sys.stdin:
		if line.split()[0] not in mers_to_delete:
			sys.stdout.write(line)

