import sys

binding = { 'A': 'T', 'T': 'A',	'C': 'G', 'G': 'C',	'_': False }


def max_consecutive_bindings(mer1, mer2):
	if len(mer2) > len(mer1):
		mer1, mer2 = mer2, mer1
	
	# reverse mer2,
	mer2 = mer2[::-1]
	mer1 = mer1.ljust(len(mer1) + len(mer1), "_")

	max_bind = 0;
	for offset in range(len(mer2)):
		consecutive = 0
		for x in range(len(mer2)):
			if binding[mer1[offset+x]] == mer2[x]:
					consecutive = consecutive + 1
			else:
				consecutive = 0

			max_bind = max(consecutive,max_bind)

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
	]
	
	print 'pass\tmer1\tmer2\tres\tcorr'
	for mer_combination in arr:
		response = []
		ans = max_consecutive_bindings(mer_combination[0], mer_combination[1])
	
		response.append(str(ans == mer_combination[2]))
		response.append(mer_combination[0])
		response.append(mer_combination[1])
		response.append(str(ans))
		response.append(str(mer_combination[2]))
	
		print '\t'.join(response)
	
