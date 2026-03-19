import argparse
import json
import math
import statistics
import sys

import msalib

def file_type(filename):
	with open(filename) as fp:
		firstline = next(fp)
		if firstline.startswith('>'): return 'fasta'
		if firstline.startswith('# STOCKHOLM'): return 'stockholm'
		sys.exit(f'unknown file type: {filename}')

def display_matrix(hmm, score, trace, beg, end):
	for j in range(len(hmm['states'])):
		print(hmm['states'][j], end='\t')
		for i in range(beg, end+1):
			if score[i][j] is None: print('N', end='|')
			else:                   print(f'{score[i][j]:.4f}', end='|')
			if trace[i][j] is None: print('N', end='\t')
			else:                   print(trace[i][j], end='\t')
		print()

##########
## CLI ##
########

parser = argparse.ArgumentParser()
parser.add_argument('hmm', metavar='<*.json>', help='hmm in json')
parser.add_argument('maf', metavar='<file>',
	help='multiple alignment file in stockholm or fasta')
parser.add_argument('--verbose', action='store_true')
arg = parser.parse_args()

###############
## Sequence ##
#############
seqs = []
uids = []
match file_type(arg.maf):
	case 'fasta':
		for defline, seq in msalib.read_fasta(arg.maf):
			seqs.append(seq)
			uids.append(defline.split()[0])
	case 'stockholm':
		with open(arg.maf) as fp:
			for line in fp:
				if line.startswith('#'): continue
				if line.startswith('//'): break
				uid, seq = line.split()
				seqs.append(seq.replace('.', '-'))
				uids.append(uid)

seq = []
for i in range(len(seqs[0])):
	col = []
	for j in range(len(seqs)):
		col.append(seqs[j][i])
	x = msalib.column_discretizer(col)
	seq.append(x)
if arg.verbose: print(seq, file=sys.stderr)

##########
## HMM ##
########

with open(arg.hmm) as fp: hmm = json.load(fp)
if arg.verbose: print('hmm:', hmm, file=sys.stderr)

# create emission lut in log2
emit = []
for probs in hmm['emissions']:
	scores = []
	for p in probs: scores.append(math.log2(p))
	scores.append(0) # gap score
	emit.append(scores)
if arg.verbose: print('emissions:', emit, file=sys.stderr)

# create transmission lut in log2
trans = []
for probs in hmm['transitions']:
	trans.append( [math.log2(x) for x in probs] )
if arg.verbose: print('transitions:', trans, file=sys.stderr)

###############
## Decoding ##
#############

## shortcuts
states = len(hmm['states'])

## init
score = [] # indexed as score[i:position][j:state]
trace = []
for i in range(len(seq)+1):
	score.append([None] * states)
	trace.append([None] * states)
for j in range(states):
	score[0][j] = math.log2(hmm['inits'][j])
	trace[0][j] = '*'

#if arg.verbose: display_matrix(hmm, score, trace, 0, 3)

## induction
for i in range(1, len(seq)+1):
	for j in range(states):
		# find maximum prev state
		max_score = -math.inf
		max_state = None
		for k in range(states):
			s = score[i-1][k] + trans[k][j] + emit[j][seq[i-1]]
			if s > max_score:
				max_score = s
				max_state = k
		# set maximum prev state
		score[i][j] = max_score
		trace[i][j] = max_state
	#if arg.verbose: display_matrix(hmm, score, trace, i, i)

## finalization
max_score = -math.inf
max_state = None
for j in range(states):
	if score[-1][j] > max_score:
		max_score = score[-1][j]
		max_state = trace[-1][j]

## trace back
path = []
i = len(seq)
j = max_state
while i > 1:
	path.append(j)
	j = trace[i][j]
	i -= 1
path.append(j)
path.reverse()

spans = []
beg = 0
state = path[0]
while i < len(path):
	if path[i] != state:
		spans.append( (state, beg, i-1) )
		state = path[i]
		beg = i
	i += 1
spans.append( (state, beg, i-1) )

## output in 1-based coordinates with state labels
for state, beg, end in spans:
	print(hmm['states'][state], beg+1, end+1, sep='\t')
