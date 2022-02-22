# sumt: command line program for computing consensus tree and other phylogenetic tree summaries

![](https://img.shields.io/badge/version-2.1.3-blue)

The command-line program `sumt` takes as input one or more files containing samples of phylogenetic trees (e.g., from a Bayesian MCMC analysis or a bootstrap procedure), and produces as output files containing (1) the majority-rule consensus tree with clade support values, (2) a summary of observed bipartitions along with branch length means and variances, and optionally (3) a list of tree topologies and how frequently they were observed. The name is taken from the `sumt` command in [MrBayes](https://nbisweden.github.io/MrBayes/index.html), whose functionality it also is meant to resemble. Clade support values and topology frequencies can be interpreted as posterior probabilities if the input trees are from a Bayesian MCMC analysis.


## Availability

The `sumt` source code is available on GitHub: https://github.com/agormp/sumt. The executable can be installed from PyPI: https://pypi.org/project/sumt/

## Installation

```
python3 -m pip install sumt
```

## Dependencies

`sumt` relies on the [phylotreelib library](https://github.com/agormp/phylotreelib), which is automatically included when using pip to install.

## Overview

* Input:
	* One or more files containing phylogenetic trees (all trees must have same leaf names), in NEXUS or Newick format.
	* Typically trees are from a bootstrap procedure or a Bayesian MCMC analysis
* Output:
	* File containing consensus tree with clade support values (= frequency of bipartition in input trees)
	* File containing list of bipartitions present in input trees, along with mean and variance of corresponding branch lengths
	* (Optionally) File containing list of observed tree topologies with posterior and cumulated probabilities
	* (Optionally) Progress indication is written to screen
* Reasonably fast and light on memory usage: 100,000 trees with 41 leaves processed in 2.5 minutes and using 13 GB memory on 2019 MacBook.
* Option to discard fraction of trees as burn-in (for Bayesian analyses)
* Option to compute average standard deviation of split frequencies (ASDSF) when multiple input files are given. This can be used as a measure of convergence of Bayesian analyses, assuming that different files represent independent MCMC samples.
* Option to include all compatible bipartitions in consensus tree (in addition to those that are present in more than 50% of input trees).
* Option to root consensus tree on outgroup or using midpoint rooting.
* Option to assign specific weights to different input files.
* Option to automatically assign weights so all files have equal impact regardless of number of trees in them.
* Option to get verbose progress indication:
	* Running counts of number of distinct bipartitions and topologies seen, along with theoretical maximum possible
	* Count of number of bipartitions in consensus tree, along with theoretical maximum
	* Report of maximum memory usage during processing
* Option to include zero length terms when computing branch lengths and ASDF (length of branch is taken to be zero in those input trees where the corresponding bipartition is not present). Not sure whether this is a good idea...

## Usage

```
Usage: sumt [options] FILE [FILE ...]
       sumt [options] -w WEIGHT FILE -w WEIGHT FILE ...

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -I FORM               format of input: nexus or newick [default: nexus]
  -O FORM               format of output: nexus or newick [default: nexus]
  -q                    quiet: don't print progress indication to terminal
                        window. NOTE: also turns on the -n option
  -v                    verbose: more information, longer error messages
  -n                    no warning when overwriting files
  --basename=NAME       base name of output files (default: derived from input
                        file)
  -b NUM                burnin: fraction of trees to discard [0 - 1; default:
                        0.25]
  -t NUM                compute tree probabilities, report NUM percent
                        credible interval [0 - 1; default: 1.00]
  -s                    compute average standard deviation of split
                        frequencies (ASDSF)
  -f NUM                Min. frequency for including bipartitions in report
                        and in computation of ASDSF [default: 0.1]
  -a                    add all compatible bipartitions to consensus tree
  -z                    include zero length terms when computing branch length
                        and average standard deviation of split frequencies
  -w WEIGHT FILE -w WEIGHT FILE ...
                        put different weights on different FILEs
  --autow               automatically assign file weights based on tree
                        counts, so all files have equal impact
  -m                    perform midpoint rooting of tree
  -r TAX [-r TAX ...]   root consensus tree on specified outgroup taxon/taxa
  --rootfile=FILE       root consensus tree on outgroup taxa listed in file
                        (one name per line)
```

## Usage examples

### Example 1: Two smallish files from MCMC analysis, with tree probabilities

The command below causes `sumt` to do the following:

* Summarise the tree samples in the files `rpoB.nexus.run1.t` and `rpoB.nexus.run2.t`
* Discard 25% of tree samples as burn-in
* Report 80% credible set of topologies (i.e., all the most frequently seen topologies to a cumulated probability of 80%)
* Compute average standard deviation of split frequencies as a measure of MCMC convergence
* Report mean, variance, and standard error of the mean of branch lengths for all bipartitions seen in more than 10% of input trees

```
sumt -b 0.25 -t 0.8 -s -f 0.1 rpoB.nexus.run1.t rpoB.nexus.run2.t
```

#### Screen output

This is printed to screen during run:

```
   Counting trees in file 'rpoB.nexus.run1.t':          1,001
   Counting trees in file 'rpoB.nexus.run2.t':          1,001

   Analyzing file: rpoB.nexus.run1.t (Weight: 0.500)
   Discarded 250 of 1,001 trees (burnin fraction=0.25)
   Processing trees ('.' signifies 100 trees):

   .......


   Analyzing file: rpoB.nexus.run2.t (Weight: 0.500)
   Discarded 250 of 1,001 trees (burnin fraction=0.25)
   Processing trees ('.' signifies 100 trees):

   .......

   Average standard deviation of split frequencies: 0.032065

   Bipartition list written to rpoB.nexus.parts
   Consensus tree written to rpoB.nexus.con
   Tree probabilities written to rpoB.nexus.trprobs

   Done. 1,502 trees analyzed.
   Time spent: 0:00:00 (h:m:s)

```

#### Bipartition overview

This is the contents of the file `rpoB.nexus.parts`. Bipartitions are indicated using the "asterisk and dots" notation also used by e.g. MrBayes: Columns correspond to taxa (with column 1 = taxon 1). All taxa with "\*" (or ".") are in the same half of the bipartition:

```
List of taxa in bipartitions:

 1 -- s100
 2 -- s501
 3 -- s502
 4 -- s503
 5 -- s504
 6 -- s505
 7 -- s506
 8 -- s507
 9 -- s508
10 -- s509
11 -- s510
12 -- s511
13 -- s512
14 -- s513


List of bipartitions:

PART = Description of partition in .* format
PROB = Posterior probability of the partition
BLEN = Mean branch length
VAR  = Branch length variance
SEM  = Standard error of the mean for branch length

PART             PROB      BLEN      VAR         SEM
*.............   1.000000  0.022560  (0.000038)  (0.000160)
.*............   1.000000  0.010498  (0.000022)  (0.000121)
..*...........   1.000000  0.002119  (0.000004)  (0.000049)
...*..........   1.000000  0.016401  (0.000031)  (0.000144)
....*.........   1.000000  0.001863  (0.000004)  (0.000053)
.....*........   1.000000  0.005484  (0.000011)  (0.000085)
......*.......   1.000000  0.008589  (0.000013)  (0.000094)
.......*......   1.000000  0.008552  (0.000017)  (0.000106)
........*.....   1.000000  0.010439  (0.000019)  (0.000112)
.........*....   1.000000  0.002088  (0.000004)  (0.000052)
..........*...   1.000000  0.002368  (0.000004)  (0.000051)
...........*..   1.000000  0.009009  (0.000016)  (0.000103)
............*.   1.000000  0.003896  (0.000008)  (0.000073)
.............*   1.000000  0.005444  (0.000010)  (0.000082)
....*........*   1.000000  0.007544  (0.000014)  (0.000095)
**.*...*......   1.000000  0.013263  (0.000023)  (0.000123)
....*......***   1.000000  0.009635  (0.000022)  (0.000120)
......*.*.....   0.976032  0.005354  (0.000010)  (0.000082)
**.*.*.*.**...   0.974035  0.003716  (0.000007)  (0.000069)
..*.*......***   0.957390  0.003515  (0.000007)  (0.000072)
.........**...   0.946738  0.003772  (0.000007)  (0.000071)
...*...*......   0.633822  0.006356  (0.000014)  (0.000123)
**............   0.483356  0.005468  (0.000013)  (0.000136)
**.*...*.**...   0.372836  0.002258  (0.000004)  (0.000086)
...........**.   0.362184  0.002125  (0.000004)  (0.000088)
**.*.*.*......   0.344208  0.003060  (0.000007)  (0.000118)
....*.......**   0.320905  0.001740  (0.000003)  (0.000084)
....*......*.*   0.316911  0.001702  (0.000003)  (0.000084)
**.*..........   0.316245  0.006348  (0.000014)  (0.000174)
.....*...**...   0.276298  0.002065  (0.000003)  (0.000090)
*..*...*......   0.250999  0.003611  (0.000006)  (0.000123)
.*.*..........   0.246338  0.006094  (0.000011)  (0.000175)
```

#### Tree probabilities

This is the content of the file `rpoB.nexus.trprobs`. Note: For data sets with more than about 15 taxa each sampled tree will typically be unique and all topologies therefore have the same probability, meaning the credible set is less useful. (Bipartitions on those trees will, however, not be unique, and clade probabilities carry useful information).

```
#NEXUS

[This file contains the 80% most probable trees found during the
MCMC search, sorted by posterior probability (the 80% HPD interval).
Lower case 'p' indicates the posterior probability of a tree.
Upper case 'P' indicates the cumulative posterior probability.]

begin trees;
   tree tree_1    [p = 0.049268] [P = 0.049268] = ((((s509,s510),(((s506,s508),(s502,((s511,s512),(s513,s504)))),s505)),(s503,s507)),s501,s100);
   tree tree_2    [p = 0.043276] [P = 0.092543] = (((s503,s507),(((s509,s510),((s506,s508),(s502,(((s513,s504),s511),s512)))),s505)),s501,s100);
   tree tree_3    [p = 0.042610] [P = 0.135153] = ((((s509,s510),(((s506,s508),(((s513,s504),(s511,s512)),s502)),s505)),s507),(s503,s501),s100);
   tree tree_4    [p = 0.039947] [P = 0.175100] = (((s503,s507),(((s509,s510),s505),((s506,s508),(((s513,s504),(s511,s512)),s502)))),s501,s100);
   tree tree_5    [p = 0.039947] [P = 0.215047] = (((((s509,s510),s505),((s506,s508),(s502,(s512,((s513,s504),s511))))),s501),(s503,s507),s100);
   tree tree_6    [p = 0.037949] [P = 0.252996] = (((s503,s507),(((s509,s510),s505),((s506,s508),(s502,(s511,((s513,s504),s512)))))),s501,s100);
   tree tree_7    [p = 0.035952] [P = 0.288948] = (((s503,s507),(((s509,s510),((s506,s508),(((s513,s504),(s511,s512)),s502))),s505)),s501,s100);
   tree tree_8    [p = 0.034621] [P = 0.323569] = (((((s509,s510),((s506,s508),((((s513,s504),s512),s511),s502))),s505),s501),(s503,s507),s100);
   tree tree_9    [p = 0.033955] [P = 0.357523] = ((((s509,s510),(((s506,s508),(s502,(((s513,s504),s512),s511))),s505)),(s503,s507)),s501,s100);
   tree tree_10   [p = 0.033289] [P = 0.390812] = (((s503,s507),((s509,s510),(s505,((s506,s508),((((s513,s504),s511),s512),s502))))),s501,s100);
   tree tree_11   [p = 0.029294] [P = 0.420107] = (((((s509,s510),s505),((s506,s508),((((s513,s504),s511),s512),s502))),(s503,s507)),s501,s100);
   tree tree_12   [p = 0.027963] [P = 0.448069] = ((s503,s507),(s501,(((s509,s510),s505),((s506,s508),(s502,(s511,((s513,s504),s512)))))),s100);
   tree tree_13   [p = 0.027963] [P = 0.476032] = ((((s509,s510),(((s506,s508),((((s513,s504),s512),s511),s502)),s505)),s501),(s503,s507),s100);
   tree tree_14   [p = 0.027297] [P = 0.503329] = ((s503,((((s509,s510),((s506,s508),(((s513,s504),(s511,s512)),s502))),s505),s507)),s501,s100);
   tree tree_15   [p = 0.026631] [P = 0.529960] = ((((s509,s510),(((s506,s508),((((s513,s504),s511),s512),s502)),s505)),s507),(s503,s501),s100);
   tree tree_16   [p = 0.025300] [P = 0.555260] = ((s503,s507),(s501,(((s509,s510),s505),((s506,s508),(s502,((s513,s504),(s511,s512)))))),s100);
   tree tree_17   [p = 0.024634] [P = 0.579893] = ((s503,((((s509,s510),((s506,s508),((((s513,s504),s512),s511),s502))),s505),s507)),s501,s100);
   tree tree_18   [p = 0.024634] [P = 0.604527] = ((s501,((s509,s510),(((s506,s508),(((s511,s512),(s513,s504)),s502)),s505))),(s503,s507),s100);
   tree tree_19   [p = 0.024634] [P = 0.629161] = (((s503,s507),(((s509,s510),((s506,s508),((((s513,s504),s512),s511),s502))),s505)),s501,s100);
   tree tree_20   [p = 0.021971] [P = 0.651132] = ((s503,s507),(s501,((s509,s510),(((s506,s508),((((s513,s504),s511),s512),s502)),s505))),s100);
   tree tree_21   [p = 0.021971] [P = 0.673103] = ((s501,(((s509,s510),((s506,s508),(((s511,s512),(s513,s504)),s502))),s505)),(s503,s507),s100);
   tree tree_22   [p = 0.018642] [P = 0.691744] = ((s501,(((s509,s510),((s506,s508),(s502,(((s513,s504),s511),s512)))),s505)),(s503,s507),s100);
   tree tree_23   [p = 0.017976] [P = 0.709720] = (((s509,(s510,(((s506,s508),(((s513,s504),(s511,s512)),s502)),s505))),s507),(s503,s501),s100);
   tree tree_24   [p = 0.016644] [P = 0.726365] = ((s503,s501),((((s509,s510),s505),((s506,s508),(s502,(((s513,s504),s511),s512)))),s507),s100);
   tree tree_25   [p = 0.016644] [P = 0.743009] = ((s503,((((s509,s510),((s506,s508),((((s513,s504),s511),s512),s502))),s505),s507)),s501,s100);
   tree tree_26   [p = 0.015979] [P = 0.758988] = (((((s509,s510),s505),((s506,s508),(s502,(s511,((s513,s504),s512))))),s507),(s503,s501),s100);
   tree tree_27   [p = 0.015313] [P = 0.774301] = (((((s509,s510),((s506,s508),(s502,(s511,((s513,s504),s512))))),s505),s507),(s503,s501),s100);
   tree tree_28   [p = 0.015313] [P = 0.789614] = (((((s509,s510),((s506,s508),(s502,(((s513,s504),s511),s512)))),s505),s507),(s503,s501),s100);
   tree tree_29   [p = 0.013981] [P = 0.803595] = ((((s509,s510),(((s506,s508),((((s513,s504),s512),s511),s502)),s505)),s507),(s503,s501),s100);
end;

```

#### Consensus tree

This is the content of the file `rpoB.nexus.con`. The difference between the two trees is the information given as branch labels:

* First tree: labels are bipartition frequencies (= posterior probability of clade, if tree sample is from Bayesian MCMC analysis)
* Second tree: labels indicate the standard error of the mean of the corresponding branch length

```
#NEXUS

begin trees;
   [In this tree branch labels indicate the posterior
    probability of the bipartition corresponding to the branch.]
   tree prob = (s513:0.00544371,s504:0.00186322,((((s506:0.00858872,s508:0.0104387)0.976:0.00535366,((s501:0.0104984,(s503:0.0164008,s507:0.00855192)0.634:0.006356,s100:0.0225601)1.000:0.0132633,(s509:0.00208824,s510:0.00236803)0.947:0.00377207,s505:0.00548352)0.974:0.00371618)0.957:0.00351458,s502:0.0021187)1.000:0.00963544,s511:0.00900914,s512:0.00389581)1.000:0.00754406);

   [In this tree branch labels indicate the standard error
    of the mean for the branch length.]
   tree sem = (s513:0.00544371,s504:0.00186322,((((s506:0.00858872,s508:0.0104387)0.000082:0.00535366,((s501:0.0104984,(s503:0.0164008,s507:0.00855192)0.000123:0.006356,s100:0.0225601)0.000123:0.0132633,(s509:0.00208824,s510:0.00236803)0.000071:0.00377207,s505:0.00548352)0.000069:0.00371618)0.000072:0.00351458,s502:0.0021187)0.000120:0.00963544,s511:0.00900914,s512:0.00389581)0.000095:0.00754406);
end;
```

### Example 2: One large file from MCMC analysis, with verbose screen output

The command below causes sumt to do the following:

* Summarise the tree samples in the file `mhc.nexus.t`
* Discard 50% of tree samples as burn-in
* Report 75% credible set of topologies (i.e., all the most frequently seen topologies to a cumulated probability of 75%)
* Report mean, variance, and standard error of the mean of branch lengths for all bipartitions seen in more than 10% of input trees
* Overwrite any existing output files with no warning
* Print more verbose output to screen, including running count of distinct bipartitions and topologies seen in input trees

```
sumt -b 0.5 -t 0.75 -f 0.1 -n -v mhc.nexus.t
```

#### Screen output

This is printed to screen during run. Numbers in parentheses (after lines of dots) give a running count of how many distinct bipartitions and topologies we have seen so far (so after the first 5000 trees have been analyzed, we have seen 511 different tree topologies, and these contain a total of 118 different bipartitions). At the end of the run the actual and theoretical maximum for number of bipartitions in the consensus tree is reported. The total number of observed topologies and bipartitions (and the theoretical maximal possible number of bipartitions) is also reported. In this case we have seen a total of 141 bipartitions in the 912 distint tree topologies. If the 912 topologies had been completely different (in the sense of not sharing any bipartitions), then the number of distinct bipartitions would have been 39,216 (so 141 is a small fraction of that, indicating that some clades are observed in a large fraction of the tree samples).

```
   Counting trees in file 'mhc.nexus.t':         30,001

   Analyzing file: mhc.nexus.t (Weight: 1.000)
   Discarded 15,000 of 30,001 trees (burnin fraction=0.50)
   Processing trees ('.' signifies 100 trees):

   ..................................................     5000   (# bip:    118    # topo:    511)
   ..................................................    10000   (# bip:    130    # topo:    720)
   ..................................................    15000   (# bip:    141    # topo:    912)


   Bipartition list written to mhc.nexus.parts
   Consensus tree written to mhc.nexus.con
   Tree probabilities written to mhc.nexus.trprobs

   Done. 15,001 trees analyzed.
   Time spent: 0:00:13 (h:m:s)


   Different topologies seen:      912
   Different bipartitions seen:    141 (theoretical maximum: 39,216)
   Internal bipartitions in consensus tree:  18 (theoretical maximum: 20)
   Memory used: 1,030.95 MB.
```