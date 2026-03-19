msa-div
=======

Goal: Divide Multiple Sequence Alignments into regions of "good" and "bad".

## Preamble ##

What does it mean for a region to be highly or poorly conserved? Given some
categorization of good and bad, an HMM would be a good choice for decoding.

## Parameter Estimation ##

A "domain" is by definition a region of function and/or conservation.
Therefore, a reasonable way to start defining "good" is to examine protein
domains previously identified by experts.

Pfam-A 38.1 contains curated multiple alignment seeds organized in the
following types:

```
    232  Coiled-coil
    155  Disordered
  13585  Domain
  12046  Family
    136  Motif
   1327  Repeat
```

The amino acid composition of each type is shown below. They are all
pretty similar. Disoredered domains have fewer hydrophobics. Coiled have
a lot of Glutamate. There may be other patterns.

```
AA	Coiled	Disord	Domain	Family	Motif	Repeat
A	0.0839	0.0837	0.0832	0.0876	0.0957	0.078
C	0.0074	0.0094	0.0149	0.0127	0.0142	0.016
D	0.0534	0.0499	0.0583	0.0534	0.0475	0.0584
E	0.1317	0.0773	0.0665	0.0623	0.0598	0.0613
F	0.0178	0.0243	0.0409	0.0443	0.0352	0.0407
G	0.0282	0.0795	0.069	0.0684	0.0731	0.0619
H	0.018	0.0231	0.0227	0.0212	0.0211	0.0235
I	0.0454	0.0251	0.058	0.0596	0.0517	0.0572
K	0.0913	0.0607	0.0532	0.0515	0.0599	0.052
L	0.1163	0.0629	0.0988	0.1048	0.1112	0.1117
M	0.023	0.0147	0.0201	0.022	0.0239	0.02
N	0.0412	0.0385	0.0395	0.0381	0.0358	0.0446
P	0.0165	0.095	0.0442	0.0451	0.0403	0.0401
Q	0.0789	0.0569	0.0377	0.0368	0.0401	0.0396
R	0.0718	0.0555	0.058	0.0564	0.0683	0.0503
S	0.0608	0.1194	0.0635	0.0641	0.0681	0.0756
T	0.0449	0.0609	0.0544	0.0536	0.0565	0.055
V	0.0457	0.0414	0.071	0.0705	0.0641	0.066
W	0.0057	0.0063	0.0137	0.0148	0.0097	0.0148
Y	0.018	0.0152	0.0324	0.0327	0.024	0.0331
```

Gap characters are least common in Coiled and Motif (ignored in the previous
table).

- Coiled 0.2017
- Disorderd 0.3026
- Domain 0.4382
- Family 0.4135
- Motif 0.2238
- Repeat 0.4472

A domain is a useful point of interest because it has a lot of alignments and
is a unit of function. Domains range from 16 to 1843 amino acids with multiple
alignment depths ranging from 1 to 4028 sequences. Typical domains are around
100 amino acids. The pieces between domains tend to be smaller. A simple
transition adjaceny matrix would be as follows:

| From | To   | Prob | Notes
|:-----|:-----|:-----|:------
| good | good | 0.99 | 100 aa long on average
| good | bad  | 0.01 |
| bad  | good | 0.02 |
| bad  | bad  | 0.98 | 50 aa long on average

What exactly are the emission probabilities of good and bad states? On can
imagine various ways to characterize a column of a multiple alignment.

- Entropy
- Percent identity
- Percent similariy
- Pairwise percent identity (average of all pairwise combinations)
- Pairwise percent similarity
- Pairwise average score (average BLOSUM62 score over all combinations)

After trying entropy, percent identity, pairwise percent identity, and pairwise
score, I chose pairwiwse score to represent the "goodness" of a multiple
alignment column. Entropy did not correlate with the other methods. While
percent identity and pairwise percent identity were highly correlated with each
other and with pairwise score, I feel that score is better than percent
identity because it captures amino acid chemistry and because protein
alignments have preferred scoring matrices to match-mismatch scoring for
decades.

Given the amino acid frequencies in Pfam 38.1, the average expected BLOSUM62
score is about -0.95. An HMM needs discrete emission values, not just one
value. I decided to generate these values by creating 1 million random
alignments of depth 5. I chose 5 to spread the range of scores and also because
I imagined that 5 is about the minimal useful depth of a multiple alignment.
Here is the discretized scoring space of 1 million random alignments where
scores are represented as integers.

```
-3 846 0.000846
-2 256632 0.256632
-1 490276 0.490276
0 215207 0.215207
1 31646 0.031646
2 4677 0.004677
3 620 0.00062
4 87 8.7e-05
5 6 6e-06
6 3 3e-06
```

To compute the observed scores of multiple alignments, I collected alignment
columns from "typical" domains in Pfam 38.1 with the following criteria:

- Alignment depth 40-80
- Alignment length 50-200

The observed columns of multiple alignments have the following average score
distributions. Scores < -3 were converted to 3 and scores > 7 were converted to
7.

```
-3 956 0.005449839811193834
-2 2287 0.013037430594351777
-1 12879 0.07341891938113534
0 50283 0.2866467523287234
1 37532 0.21395751861268514
2 25677 0.14637608455232645
3 19316 0.1101141273985566
4 11338 0.0646341880536775
5 6801 0.03877025162754107
6 5170 0.029472460066811843
7 3179 0.01812242757299707
```

Idealized emission probabilities for bad (Random) and good (Observed) are shown
in the table below (the Code is how the value is implemented in the program).

| Score | Random | Obs  | Code |
|:------|:-------|:-----|:-----|
|  -2   | 0.250  | 0.02 |   0  |
|  -1   | 0.500  | 0.07 |   1  |
|   0   | 0.200  | 0.30 |   2  |
|   1   | 0.035  | 0.20 |   3  |
|   2   | 0.005  | 0.15 |   4  |
|   3   | 0.004  | 0.11 |   5  |
|   4   | 0.003  | 0.06 |   6  |
|   5   | 0.002  | 0.05 |   7  |
|   6   | 0.001  | 0.04 |   8  |
|  gap  | x      | x    |   9  |

Gaps cannot be considered to be paired with amino acids and must be handled
heuristically in decoding.

## Decoding ##

Here is an HMM in JSON where "hi" is observed and "lo" is random. Note that
there are 9 emission values representing scores from -2 to 6. Gaps are not
described by the HMM and are "hacked" into the algorithm, having the same
emission probability in all states.

```
{
  "states": ["hi", "lo"],
  "inits": [0.5, 0.5],
  "transitions": [
    [0.99, 0.01],
    [0.02, 0.98]
  ],
  "emissions": [
    [0.020, 0.070, 0.300, 0.200, 0.150, 0.110, 0.060, 0.050, 0.040],
    [0.250, 0.500, 0.200, 0.035, 0.005, 0.004, 0.003, 0.002, 0.001]
  ]
}
```

Given the sequences below, the HMM should split this MSA into 3 parts:

1) Good from 1-10 (high conservation, all A)
2) Bad from 11-20 (very little conservation)
3) Good from 21-50 (the gap region should be ignored)

```
S1_TEST/1-50	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
S2_TEST/1-50	AAAAAAAAAACCCCCCCCCCAAAAAAAAAA----------AAAAAAAAAA
S3_TEST/1-50	AAAAAAAAAADDDDDDDDDDAAAAAAAAAA----------AAAAAAAAAA
S4_TEST/1-50	AAAAAAAAAAEEEEEEEEEEAAAAAAAAAA----------AAAAAAAAAA
S5_TEST/1-50	AAAAAAAAAAFFFFFFFFFFAAAAAAAAAA----------AAAAAAAAAA
S6_TEST/1-50	AAAAAAAAAAGGGGGGGGGGAAAAAAAAAA----------AAAAAAAAAA
S7_TEST/1-50	AAAAAAAAAAHHHHHHHHHHAAAAAAAAAA----------AAAAAAAAAA
S8_TEST/1-50	AAAAAAAAAAIIIIIIIIIIAAAAAAAAAA----------AAAAAAAAAA
S9_TEST/1-50	AAAAAAAAAAKKKKKKKKKKAAAAAAAAAA----------AAAAAAAAAA
```

You can test the code as follows:

```
python3 msa-div.py hilo.json test.msa
```

Which produces the following output:

```
hi	1	10
lo	11	20
hi	21	50
```
