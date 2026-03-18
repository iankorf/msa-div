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

A domain is based on a region of sequence whereas a family is the entire
protein. Domains range from 16 to 1843 amino acids with multiple alignment
depths ranging from 1 to 4028 sequences. Typical domains are around 100 amino
acids. The pieces between domains tend to be smaller. A simple transition
adjaceny matrix would be as follows:

| From | To   | Prob |
|:-----|:-----|:------
| good | good | 0.99 |
| good | bad  | 0.01 |
| bad  | good | 0.02 |
| bad  | bad  | 0.98 |

What exactly are the emission probabilities of good and bad states? Some
measure of how uniform the alignment is. After trying entropy, percent
identity, pairwise percent identity, and pairwise score, I chose pairwiwse
score (entropy did not correlate with the others and percent identity of any
measure doesn't capture amino acid similarity as well as a scoring matrix).

Given the amino acid frequencies in Pfam 38.1, the average expected BLOSUM62
score is about -0.95. An HMM needs discrete emission values, not just one
value. I decided to generate these values by creating 1 million random
alignments of depth 5. I chose 5 to spread the range of stochasticity and also
because I imagined that 5 is about the minimal useful depth of a multiple
alignment. Here is the discretized scoring space of 1 million random
alignments.

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
in the table below. These will be looked up in a emission table via the Code.

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

A column that is mostly gaps is due to a couple sequences with different
lengths rather than different sequences. Such gaps should not affect decoding.

## Decoding ##

HMM in JSON where "hi" is observed and "lo" is random. Note that the properties
of gaps are not in the HMM, as the HMM describes the conservation of columns
and gaps represent different lengths of sequences (when large). The decoding of
gaps will be handled heuristically.

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
