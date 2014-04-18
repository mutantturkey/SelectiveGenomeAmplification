#!/usr/bin/env python2.7
import argparse
import sys
import math

"""Calculate the thermodynamic melting temperatures of nucleotide sequences."""

def Tm_staluc(s, dna=5000, na=10, mg=20, dntps=10):
		"""Returns DNA/DNA tm using nearest neighbor thermodynamics.

		dnac is DNA concentration [nM]
		na is salt concentration [mM].
		mg, the concentration of magnesium in mM
		dNTPs, the concentration of dNTPs in mM
		
		Sebastian Bassi <sbassi@genesdigitales.com>"""
		
		#Credits: 
		#Main author: Sebastian Bassi <sbassi@genesdigitales.com>
		#Overcount function: Greg Singer <singerg@tcd.ie>
		#Based on the work of Nicolas Le Novere <lenov@ebi.ac.uk> Bioinformatics.
		#17:1226-1227(2001)

		#This function returns better results than EMBOSS DAN because it uses
		#updated thermodynamics values and takes into account inicialization
		#parameters from the work of SantaLucia (1998).
		
		#Things to do:
		#+Detect complementary sequences. Change K according to result.
		#+Add support for heteroduplex (see Sugimoto et al. 1995).
		#+Correction for Mg2+. Now supports only monovalent ions.
		#+Put thermodinamics table in a external file for users to change at will
		#+Add support for danglings ends (see Le Novele. 2001) and mismatches.
		
		dh = 0 #DeltaH. Enthalpy
		ds = 0 #deltaS Entropy

		def tercorr(stri):
				deltah = 0
				deltas = 0
				#DNA/DNA
				#Allawi and SantaLucia (1997). Biochemistry 36 : 10581-10594
				if stri.startswith('G') or stri.startswith('C'):
						deltah -= 0.1
						deltas += 2.8
				elif stri.startswith('A') or stri.startswith('T'):
						deltah -= 2.3
						deltas -= 4.1
				if stri.endswith('G') or stri.endswith('C'):
						deltah -= 0.1
						deltas += 2.8
				elif stri.endswith('A') or stri.endswith('T'):
						deltah -= 2.3
						deltas -= 4.1
				dhL = dh + deltah
				dsL = ds + deltas
				return dsL,dhL

		def overcount(st,p):
				"""Returns how many p are on st, works even for overlapping"""
				ocu = 0
				x = 0
				while 1:
						try:
								i = st.index(p,x)
						except ValueError:
								break
						ocu += 1
						x = i + 1
				return ocu

		R = 1.987 # universal gas constant in Cal/degrees C*Mol
		sup = s.upper()
		vsTC,vh = tercorr(sup)
		vs = vsTC
		
		k = (dna/4.0)*1e-9
		#With complementary check on, the 4.0 should be changed to a variable.
		
		#DNA/DNA
		#Allawi and SantaLucia (1997). Biochemistry 36 : 10581-10594
		vh = vh + (overcount(sup,"AA"))*7.9 + (overcount(sup,"TT"))*\
		7.9 + (overcount(sup,"AT"))*7.2 + (overcount(sup,"TA"))*7.2 \
		+ (overcount(sup,"CA"))*8.5 + (overcount(sup,"TG"))*8.5 + \
		(overcount(sup,"GT"))*8.4 + (overcount(sup,"AC"))*8.4
		vh = vh + (overcount(sup,"CT"))*7.8+(overcount(sup,"AG"))*\
		7.8 + (overcount(sup,"GA"))*8.2 + (overcount(sup,"TC"))*8.2
		vh = vh + (overcount(sup,"CG"))*10.6+(overcount(sup,"GC"))*\
		9.8 + (overcount(sup,"GG"))*8 + (overcount(sup,"CC"))*8
		vs = vs + (overcount(sup,"AA"))*22.2+(overcount(sup,"TT"))*\
		22.2 + (overcount(sup,"AT"))*20.4 + (overcount(sup,"TA"))*21.3
		vs = vs + (overcount(sup,"CA"))*22.7+(overcount(sup,"TG"))*\
		22.7 + (overcount(sup,"GT"))*22.4 + (overcount(sup,"AC"))*22.4
		vs = vs + (overcount(sup,"CT"))*21.0+(overcount(sup,"AG"))*\
		21.0 + (overcount(sup,"GA"))*22.2 + (overcount(sup,"TC"))*22.2
		vs = vs + (overcount(sup,"CG"))*27.2+(overcount(sup,"GC"))*\
		24.4 + (overcount(sup,"GG"))*19.9 + (overcount(sup,"CC"))*19.9
		ds = vs
		dh = vh
		
		fgc = lambda s: len(filter(lambda x: x =='G' or x =='C', s)) / float(len(s))

		tm = ((1000* (-dh))/(-ds+(R * (math.log(k)))))-273.15
		Mmg = mg * 1E-3
		Mna = na * 1E-3
		Mdntp = dntps * 1E-3

		Fmg = (-(3E4 * Mdntp - 3E4 * Mmg + 1) + ((3E4 * Mdntp - 3E4 * Mmg + 1)**2 + 4 * 3E4 * Mmg )**0.5 )/ (2 * 3E4)
		cationratio = Fmg**0.5 / Mna
		
		if cationratio < 0.22:
			SaltCorrectedTm = 1 / ( (1 / tm) + ((4.29 * fgc(s) - 3.95) * math.log(Mna) + 0.940 * (math.log(Mna))**2) * 1E-5 )
		else:
			a = 0; d = 0; g = 0;
			if cationratio < 6.0:
				a = 3.92 * (0.843 - 0.352 * (Mna)^0.5 * math.log(Mna) )
				d = 1.42 * (1.279 - 4.03 * math.log(Mna) * 1E-3 - 8.03 * (math.log(Mna))^2 * 1E-3 )
				g = 8.31 * (0.486 - 0.258 * math.log(Mna) + 5.25 * (math.log(Mna))^3 * 1E-3)
			elif cationratio > 6.0:
				a = 3.92
				d = 1.42
				g = 8.31

			SaltcorrectedTm = 1 / ( (1 / tm) + (a - 0.91 * math.log(Fmg) + fgc(s) * (6.26 + d * math.log(Fmg)) + 1/(2 * (len(s) - 1)) * (-48.2 + 52.5 * math.log(Fmg) + g * (math.log(Fmg))**2)) * 1E-5 )

		return tm

# naiive
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



if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="Filter mers, input is stdin, output is stdout") 
	parser.add_argument("-m", "--min", help="min melting temp [c]", required=True, type=float)
	parser.add_argument("-x", "--max", help="max melting temp [c]", required=True, type=float)
	parser.add_argument("-d", "--dna", help="DNA concentation [nM]", required=False, type=int, default=5000)
	parser.add_argument("-n", "--na", help="sodium concentation [mM]", required=False, type=int, default=10)
	parser.add_argument("-g", "--mg", help="magnesium concentation [mM]", required=False, type=int, default=20)
	parser.add_argument("-p", "--dntp", help="dNTPs concentation [mM]", required=False, type=int, default=10)

	ar = parser.parse_args()

	for line in sys.stdin:
		if ar.min < Tm_staluc(line.split("\t")[0], dna=ar.dna, na=ar.na, mg=ar.mg, dntps=ar.dntp) < ar.max:
			sys.stdout.write(line)
