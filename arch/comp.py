import sys
import msalib

count = {}
for msa in msalib.read_stockholm(sys.argv[1]):
	for seq in msa.seqs:
		for aa in seq:
			if aa not in count: count[aa] = 0
			count[aa] += 1

del count['.'] # there are many
del count['X'] # there are 3390
del count['Z'] # there are 19
del count['B'] # there are 35
del count['O'] # there are 19
del count['U'] # there are 126

total = sum(count.values())
for aa in sorted(count):
	print(aa, count[aa], count[aa]/total)
