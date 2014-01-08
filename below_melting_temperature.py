import sys
from multiprocessing import Pool

def above_melting_temperature(kmer_with_count):
	kmer = kmer_with_count.split("\t")[0]

	A = kmer.count('A')
	C = kmer.count('C')
	G = kmer.count('G')
	T = kmer.count('T')

	melt_temp = 0.0;

	if len(kmer) < 13:
		melt_temp = ((A+T) * 2) + ((C+G) * 4)
	else:
		melt_temp = 64.9 + 41*(G+C-16.4)/(A+T+G+C)

	if melt_temp < max_melting_temp:
		return kmer_with_count
	else:
		return 0


lines = sys.stdin.readlines()
max_melting_temp = float(sys.argv[1])

p = Pool()
output = p.map(above_melting_temperature, lines)

for line in output:
	if line != 0:
		sys.stdout.write(line)
