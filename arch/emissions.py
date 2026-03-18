import argparse
import random
import sys

from statistics import mean, stdev

import msalib

parser = argparse.ArgumentParser()
parser.add_argument('seeds', metavar='<file>',
	help='Pfam seed alignments in STOCKHOLM format')
parser.add_argument('--min-depth', type=int, metavar='<int>', default=40,
	help='minimum alignment depth [%(default)i]')
parser.add_argument('--max-depth', type=int, metavar='<int>', default=80,
	help='maximum alignment depth [%(default)i]')
parser.add_argument('--min-length', type=int, metavar='<int>', default=50,
	help='minimum alignment length [%(default)i]')
parser.add_argument('--max-length', type=int, metavar='<int>', default=200,
	help='maximum alignment length [%(default)i]')
parser.add_argument('--testing', action='store_true')
arg = parser.parse_args()

# expected pairings
alph = list(msalib.AA.keys())
w8s = list(msalib.AA.values())
exp_scores = {}
limit = 1_000_000
for _ in range(limit):
	scores = []
	seq = random.choices(alph, weights=w8s, k=5)
	for i, a in enumerate(seq):
		for b in seq[i+1:]:
			s = msalib.B62[a][b]
			scores.append(s)
	rs = round(mean(scores))
	if rs not in exp_scores: exp_scores[rs] = 0
	exp_scores[rs] += 1
print('random expectation with k=5')
for score in sorted(exp_scores):
	print(score, exp_scores[score], exp_scores[score]/limit)

# observed pairings
print('observed alignment columns')
obs_scores = {}
for i, msa in enumerate(msalib.read_stockholm(arg.seeds)):
	if msa.type != 'Domain': continue
	if msa.depth < arg.min_depth: continue
	if msa.depth > arg.max_depth: continue
	if msa.length < arg.min_length: continue
	if msa.length > arg.max_length: continue

	# observed pairings
	for n in range(msa.length):
		col = msa.column(n)
		scores = []
		for i, a in enumerate(col):
			if a not in msalib.B62: continue
			for b in col[i+1:]:
				if b not in msalib.B62: continue
				s = msalib.B62[a][b]
				if s < -3: s = -3
				if s > 7: s = 7
				scores.append(s)
		if len(scores) == 0: continue
		rs = round(mean(scores))
		if rs not in obs_scores: obs_scores[rs] = 0
		obs_scores[rs] += 1

	if arg.testing and i > 70: break

total = sum(obs_scores.values())
for score in sorted(obs_scores):
	print(score, obs_scores[score], obs_scores[score]/total)
