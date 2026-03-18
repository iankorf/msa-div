import os
import sys
import msalib

types = ('Coiled-coil', 'Disordered', 'Domain', 'Family', 'Motif', 'Repeat')
for t in types:
	os.system(f'mkdir -p build/{t}')

aas = 'ACDEFGHIKLMNPQRSTVWY'
comp = {}
for t in types:
	comp[t] = {}
	for aa in aas: comp[t][aa] = 1

for i, msa in enumerate(msalib.read_stockholm(sys.argv[1])):
	#if i == 50: break
	t = msa.type
	a = msa.accession
	with open(f'build/{t}/{a}', 'w') as output: msa.write(output)
	for seq in msa.seqs:
		for aa in seq:
			if aa not in comp[t]: continue
			comp[t][aa] += 1

totals = {}
for t in types:
	totals[t] = sum(comp[t].values())

for aa in aas:
	print(aa, end='\t')
	for t in types:
		print(round(comp[t][aa]/totals[t], ndigits=4), end='\t')
	print()