import sys

def read_blosum(filename):
	"""Reads blosum scoring matrix into 2D dictionary"""
	alphabet = []
	matrix = {}
	with open(filename) as fp:
		for line in fp:
			if line.startswith('#'): continue
			if line.startswith(' '):
				f = line.split()
				for c in f: alphabet.append(c)
			elif line:
				f = line.split()
				c1 = f[0]
				if c1 not in matrix: matrix[c1] = {}
				for c2, v in zip(alphabet, f[1:]):
					matrix[c1][c2] = int(v)
	return matrix

blosum = read_blosum(sys.argv[1])

for letter, row in blosum.items():
	print("'", letter, "':{", sep='', end='')
	for letter, value in blosum[letter].items():
		print("'", letter, "':", value, ',', sep='', end='')
	print('},')
