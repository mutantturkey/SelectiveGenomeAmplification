import sys

def in_temp_range(kmer):

	A = kmer.count('A')
	C = kmer.count('C')
	G = kmer.count('G')
	T = kmer.count('T')

	melt_temp = 0.0;

	if len(kmer) < 13:
		melt_temp = ((A+T) * 2) + ((C+G) * 4)
	else:
		melt_temp = 64.9 + 41*(G+C-16.4)/(A+T+G+C)

	return min_melting_temp < melt_temp < max_melting_temp

min_melting_temp = float(sys.argv[1])
max_melting_temp = float(sys.argv[2])


output = []
for line in sys.stdin:
	if in_temp_range(line.split("\t")[0]):
		sys.stdout.write(line)