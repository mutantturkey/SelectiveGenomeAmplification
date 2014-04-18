#!/usr/bin/env python2.7
import sys, os

binding = { 'A': 'T', 'T': 'A',	'C': 'G', 'G': 'C',	'_': False }


def max_consecutive_binding(mer1, mer2):
	if len(mer2) > len(mer1):
		mer1, mer2 = mer2, mer1
	
	# save the len because it'll change when we do a ljust
	mer1_len = len(mer1)
	# reverse mer2,
	mer2 = mer2[::-1]
	# pad mer one to avoid errors
	mer1 = mer1.ljust(mer1_len + len(mer2), "_")

	max_bind = 0;
	for offset in range(mer1_len):
		consecutive = 0
		for x in range(len(mer2)):
			if binding[mer1[offset+x]] == mer2[x]:
				consecutive += 1
				if consecutive > max_bind:
					max_bind = consecutive
			else:
				consecutive = 0

	return max_bind

def test():
	# mer 1  				 mer 2					# correct ans
	arr = [
	("ATATAT",			 "TATATA", 					5),
	("ACAGGGAT",		 "ATATATAT", 				2),
	("CATATATAT", 	 "ATATATATATAT",	  8),
	("ATATATATATAT", "ATATATATAT", 		 10),
	("ATATAT",			 "TATAT",					  5),
	("AACGATACCATG", "GGATCATACGTA", 		3),
	("CGT",          "ACG",			 				3),
	("ACG",          "CGT",			 				3),
	("CACC", 				 "GGTGT", 					4),
	("GGTGT", 			 "CACC", 					4),
	("CCCCCCCCATATAT", "TATA", 4),
	]
	
	print 'pass\tmer1\tmer2\tres\tcorr'
	for mer_combination in arr:
		response = []
		ans = max_consecutive_binding(mer_combination[0], mer_combination[1])
	
		response.append(str(ans == mer_combination[2]))
		response.append(mer_combination[0])
		response.append(mer_combination[1])
		response.append(str(ans))
		response.append(str(mer_combination[2]))
	
		print '\t'.join(response)
	
def main():

	if(len(sys.argv) < 2):
		print "cutoff is expected as an argument"
		exit()
	else:
		cutoff = int(sys.argv[1])
	
	for line in sys.stdin:
		mer = line.split()[0]
		if max_consecutive_binding(mer, mer) < cutoff:
			sys.stdout.write(line)


if __name__ == "__main__":
	sys.exit(main())
