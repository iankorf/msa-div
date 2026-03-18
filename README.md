msa-div
=======

Goal: Divide Multiple Sequence Alignments into regions of "good" and "bad".

## Preamble ##

What does it mean for a region to be highly or poorly conserved?

## Data ##

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
depths ranging from 1 to 4028 sequences.

## Emission Probabilities ##

Given the amino acid frequencies in Pfam 38.1, the average expected BLOSUM62
score is about -0.95. Given alignment depths of 5, a random distribution of 1
million columns have the following score distribution.

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

The observed columns of multiple alignments have the following average score
distributions (scores < -3 = 3 and > 7 = 7).

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

An idealized version looks something like this:

| Score | Random | Obs  |
|:------|:-------|:-----|
|  -2   | 0.250  | 0.02 |
|  -1   | 0.500  | 0.07 |
|   0   | 0.200  | 0.30 |
|   1   | 0.035  | 0.20 |
|   2   | 0.005  | 0.15 |
|   3   | 0.004  | 0.11 |
|   4   | 0.003  | 0.06 |
|   5   | 0.002  | 0.05 |
|   6   | 0.001  | 0.04 |

## Categorizing ##

