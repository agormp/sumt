# sumt

![](https://img.shields.io/badge/version-3.0.0-blue)
[![PyPI downloads](https://static.pepy.tech/personalized-badge/sumt?period=total&units=none&left_color=black&right_color=blue&left_text=downloads&service=github)](https://pepy.tech/project/sumt)

The command-line program `sumt` computes consensus trees and other tree-summary statistics for sets of phylogenetic trees. The input trees can be in one or more input files, and will typically be from a Bayesian MCMC analysis (BEAST or MrBayes for instance) or from a bootstrap procedure. 

`sumt` can compute three different kinds of main tree summaries:

* Majority rule consensus tree
* Majority rule consensus tree, with all compatible bipartitions added
* Maximum bipartition credibility tree

Branch labels on these trees indicate clade support. 

`sumt` also produces a summary of observed bipartitions along with branch length means and variances and optionally, a list of tree topologies and how frequently they were observed. The name is taken from the `sumt` command in [MrBayes](https://nbisweden.github.io/MrBayes/index.html), whose functionality it also is meant to resemble. Clade support values and topology frequencies can be interpreted as posterior probabilities if the input trees are from a Bayesian MCMC analysis.


## Availability

The `sumt` source code is available on GitHub: https://github.com/agormp/sumt. The executable can be installed from PyPI: https://pypi.org/project/sumt/

## Installation

```
python3 -m pip install sumt
```

Upgrading to latest version:

```
python3 -m pip install --upgrade sumt
```

## Dependencies

`sumt` relies on the [phylotreelib library](https://github.com/agormp/phylotreelib) for phylogeny-related matters, and on [psutil](https://pypi.org/project/psutil/) for (optionally) monitoring memory usage. These are automatically included when using pip to install.

## Overview

* Input:
	* One or more files containing phylogenetic trees (all trees must have same leaf names), in NEXUS or Newick format.
	* Typically trees are from a Bayesian MCMC analysis, but could also be from a bootstrap procedure
	* Can read sample files from BEAST and MrBayes
* Output:
	* File containing summary tree with clade support values (= frequency of bipartition in input trees). Can be one of these:
		* Majority rule consensus tree
		* Majority rule consensus tree, with all compatible bipartitions added
		* Maximum bipartition credibility tree (similar to maximum clade credibility tree)
	* The file also contains a second consensus tree where branch labels indicate bipartition IDs, which can be used for interpreting bipartition file below.
	* File containing list of bipartitions present in input trees, along with mean and variance of corresponding branch lengths. This list includes both bipartitions that correspond to branches in the summary tree, and bipartitions not included in the summary
	* (Optionally) File containing list of observed tree topologies with posterior and cumulated probabilities
	* (Optionally) Progress indication is written to screen
* Optimized for speed and memory usage:
	* 100,000 trees with 41 leaves processed in 41 s, using max 50 MB memory on 2021 MacBook (4,440 distinct bipartitions seen)
	* Same file processed in 45 s, using max 4.0 GB memory when also keeping track of topologies (74,283 distinct topologies seen)
* Option to discard fraction of trees as burn-in (for Bayesian analyses)
* Option to compute average standard deviation of split frequencies when multiple input files are given. This can be used as a measure of convergence of Bayesian analyses, assuming that different files represent independent MCMC samples.
* Option to include all compatible bipartitions in consensus tree (in addition to those that are present in more than 50% of input trees).
* Option to root consensus tree using either outgroup, midpoint, or [minimum variance](https://pubmed.ncbi.nlm.nih.gov/28800608/) rooting.
* Option to assign specific weights to different input files.
* Option to automatically assign weights so all files have equal impact regardless of number of trees in them.
* Option to set basename of output files (default: basename will be stem of input file name)
* Option to get verbose progress indication:
	* Running counts of number of distinct bipartitions and topologies seen, along with theoretical maximum possible
	* Count of number of bipartitions in consensus tree, along with theoretical maximum
	* Report of maximum memory usage during processing

## Usage

```
usage: sumt    [-h] [-w WEIGHT INFILE] [--autow] [-i FORMAT] [--con | --all | --mbc]
               [-b NUM] [-t NUM] [-s] [-f NUM] [-n] [-v] [-q] [--basename NAME]
               [--rootmid | --rootminvar | -r TAXON [TAXON ...] | --rootfile FILE]
               [INFILE ...]

Computes summary tree and statistics from set of phylogenetic trees

options:
  -h, --help            show this help message and exit

Input tree files:
  INFILE                input FILE(s) containing phylogenetic trees (can list several files)
  -w WEIGHT INFILE      input FILEs with specified weights (repeat -w option for each input
                        file)
  --autow               automatically assign file weights based on tree counts, so all files
                        have equal impact (default is for all trees, not files, to be equally
                        important)
  -i FORMAT             format of input files: nexus, newick [default: nexus]

Type of summary tree:
  --con                 majority rule consensus tree [default]
  --all                 majority rule consensus tree with all compatible bipartitions added
  --mbc                 Maximum Bipartition Credibility (MBC) tree. MBC is similar to MCC
                        (Maximum Clade Credibility) tree but counting bipartitions instead
                        of clades, i.e. ignoring rooting. Additionally, branch lengths are
                        estimated from branch lengths of bipartitions and not from node
                        depths (i.e., again ignoring rooting)

Bayesian phylogeny options:
  -b NUM                burnin: fraction of trees to discard [0 - 1; default: 0.25]
  -t NUM                compute tree probabilities, report NUM percent credible interval [0
                        - 1]
  -s                    compute average standard deviation of split frequencies (ASDSF)
  -f NUM                Minimum frequency for including bipartitions in report and in
                        computation of ASDSF [default: 0.1]

Output to terminal and files:
  -n                    no warning when overwriting files
  -v                    verbose: more information, longer error messages
  -q                    quiet: don't print progress indication to terminal window. NOTE:
                        also turns on the -n option
  --basename NAME       base name of output files (default: derived from input file)

Rooting of summary tree:
  --rootmid             perform midpoint rooting of tree
  --rootminvar          perform minimum variance rooting of tree
  -r TAXON [TAXON ...]  root consensus tree on specified outgroup taxon/taxa
  --rootfile FILE       root consensus tree on outgroup taxa listed in file (one name per
                        line)
```

## Usage examples

### Example 1: 
#### Majority rule consensus tree, bipartition summary, topology summary, 
computation of average standard deviation of split frequencies, midpoint rooting

The command below causes `sumt` to do the following:

* Summarise the tree samples in the files `hiv.nexus.run1.t` and `hiv.nexus.run2.t`
* `--con`: Compute majority rule consensus tree (this is default and could have been omitted)
* `-b 0.25`: Discard 25% of tree samples as burn-in
* `-t 0.99`: Keep track of topology probabilities, report 99% credible set
* `-s`: Compute average standard deviation of split frequencies as a measure of MCMC convergence
* `-f 0.1`: Report mean, variance, and standard error of the mean of branch lengths for all bipartitions seen in more than 10% of input trees
* `--rootmid`: Perform midpoint rooting

```
sumt --con -b 0.25 -t 0.99 -f 0.1 --rootmid primates.nexus.run1.t primates.nexus.run2.t 
```

#### Screen output

This is printed to screen during run:

```
   Counting trees in file primates.nexus.run1.t          2,001
   Counting trees in file primates.nexus.run2.t          2,001

   Analyzing file: primates.nexus.run1.t (Weight: 0.500)
   Discarded 500 of 2,001 trees (burnin fraction=0.25)
   Processing trees ('.' signifies 100 trees):

   ...............


   Analyzing file: primates.nexus.run2.t (Weight: 0.500)
   Discarded 500 of 2,001 trees (burnin fraction=0.25)
   Processing trees ('.' signifies 100 trees):

   ...............

   Consensus tree written to primates.con
   Bipartition list written to primates.parts
   Tree probabilities written to primates.trprobs

   Done. 3,002 trees analyzed.
   Time spent: 0:00:00 (h:m:s)
```

#### Bipartition overview

Below are the contents of the file `primates.parts`, which lists information about bipartitions observed in input trees. A bipartition - or split - is a division of leaf-names into two groups: those on one side of an internal branch, and those on the other side. Two trees with different topologies can have the same bipartition, and that internal branch is then said to have been observed in both trees. `sumt` keeps track of bipartition frequencies, as well as means and variances of the lengths of the corresponding branches. For branches included on the summary tree, the clade support value is the same as the value in the column PROB (i.e., the frequency of the bipartion among the input trees).

Bipartitions are indicated using the "asterisk and dots" notation also used by e.g. MrBayes: Columns correspond to taxa (with column 1 = taxon 1). All taxa with "\*" (or ".") are in the same half of the bipartition. The ID column contains either a leaf name or a numerical branchID that corresponds to branch labels given in the second tree in the consensus tree file (in this case there are 5 leaves, and therefore only two internal branches). Branch IDs are numbered consecutively in order of bipartition frequency (so the branch with ID=7 has the 7th highest observed bipartition frequency).

```
List of bipartitions:

PART = Description of partition in .* format
PROB = Posterior probability of the partition
BLEN = Mean branch length
VAR  = Branch length variance
SEM  = Standard error of the mean for branch length
ID   = Leaf name or internal branch label, for those bipartitions that are included in consensus tree

PART    PROB      BLEN       VAR          SEM          ID
*....   1.000000  0.07585    (0.0003568)  (0.0003448)  Orangutan
.*...   1.000000  0.4324     (0.02003  )  (0.002583 )  Gibbon
..*..   1.000000  0.07269    (0.0005394)  (0.0004239)  Gorilla
...*.   1.000000  0.05678    (0.0002647)  (0.000297 )  Human
....*   1.000000  0.2834     (0.009653 )  (0.001793 )  Chimpanzee
.*..*   1.000000  0.1449     (0.003743 )  (0.001117 )  6
*..*.   0.965356  0.03702    (0.0004153)  (0.0003786)  7
```

#### Tree probabilities

This is the content of the file `primates.trprobs`. In this case there were only 5 leafs corresponding to a total of 15 possible trees, of which 3 were seen in the MCMC samples. Note: For data sets with more than about 15-20 taxa, each sampled tree will typically be unique and all topologies therefore have the same probability, meaning the credible set is not very useful. (Bipartitions on those trees will, however, not be unique, and clade probabilities carry useful information).

```
#NEXUS

[This file contains the 99% most probable trees found during the
MCMC search, sorted by posterior probability (the 99% HPD interval).
Lower case 'p' indicates the posterior probability of a tree.
Upper case 'P' indicates the cumulative posterior probability.]

begin trees;
    translate
        1     Chimpanzee,
        2     Gibbon,
        3     Gorilla,
        4     Human,
        5     Orangutan
    ;
    tree tree_1 [p = 0.965356] [P = 0.965356] = (((5,2),3),4,1);
    tree tree_2 [p = 0.018321] [P = 0.983678] = ((5,2),(3,1),4);
    tree tree_3 [p = 0.016322] [P = 1.000000] = (((5,2),1),4,3);
end;
```

#### Consensus tree

This is the content of the file `primates.nexus.con`. The difference between the two trees is the information given as branch labels:

* First tree: labels are bipartition frequencies (= posterior probability of bipartition, if tree samples are from Bayesian MCMC analysis)
* Second tree: labels are the branchIDs also indicated in the bipartition summary in the file `primates.parts`. This should make it simpler to understand what branch the bipartition corresponds to (open the tree file in a treeviewer such as FigTree and view the branch labels).

```
#NEXUS

begin trees;
   [In this tree branch labels indicate the posterior probability of the bipartition corresponding to the branch.]
   tree prob = (((Human:0.056782,Chimpanzee:0.0758462)0.965:0.0370172,Gorilla:0.0726937)1.000:0.144909,Orangutan:0.283409,Gibbon:0.432437);

   [In this tree branch labels indicate the bipartition ID listed in the file primates.parts.
    These branch labels can be used for interpreting the table of branch lenght info in that file]
   tree partID = (((Human:0.056782,Chimpanzee:0.0758462)7:0.0370172,Gorilla:0.0726937)6:0.144909,Orangutan:0.283409,Gibbon:0.432437);
end;
```

### Example 2: 
#### Maximum bipartition credibility tree, topology summary, verbose output, setting basename of output files

The command below causes sumt to do the following:

* Summarise the tree samples in the file `gp120.nexus.trees`
* `--mbc`: Compute maximum bipartition credibility tree (instead of majority rule consensus)
* `-b 0.1`: Discard 10% of tree samples as burn-in
* `-t 0.75`: Report 75% credible set of topologies (i.e., all the most frequently seen topologies to a cumulated probability of 75%)
* `-n`: Overwrite any existing output files with no warning
* `-v`: Print more verbose output to screen, including running count of distinct bipartitions and topologies seen in input trees
* `--basename /Users/bob/hiv`: produce output files with the indicated stem (/Users/bob/hiv.parts, /Users/bob/hiv.trprobs, /Users/bob/hiv.mbc)

```
sumt --mbc -b 0.1 -t 0.75 -nv --basename /Users/bob/hiv gp120.nexus.trees
```

#### Screen output

This is printed to screen during run. Numbers in parentheses (after lines of dots) give a running count of how many distinct bipartitions and topologies we have seen so far (so after the first 5000 trees have been analyzed, we have seen 4732 different tree topologies, and these contain a total of 842 different bipartitions). 

At the end of the run the actual and theoretical maximum for number of bipartitions in the consensus tree is reported. The total number of observed topologies and bipartitions (and the theoretical maximal possible number of bipartitions) is also reported. In this case we have seen 34,127 distint tree topologies among the 36,001 trees analyzed, meaning that about 95% of the sampled trees are unique. In these topologies we have seen 1,064 distinct bipartitions. If the 34,127 topologies had been completely different (in the sense of not sharing any bipartitions), then the number of distinct bipartitions would have been 1,262,699 (so 1,064 is a small fraction of that, indicating that some bipartitions are observed in a large fraction of the tree samples).

The ouput also indicates that the summary tree is fully resolved (has no polytomies). This will always be the case for a maximum bipartition credibility tree (because individual tree samples are typically fully resolved). It is further stated that the MBC tree has not been explicitly rooted (none of `sumt`'s rooting options were used).

Finally the Highest Log Bipartition Credibility is output (this is the sum of the logs of the bipartition frequencies, for those bipartitions that are present in the MBC tree).

```
   Counting trees in file gp120.nexus.trees         40,001

   Analyzing file: gp120.nexus.trees (Weight: 1.000)
   Discarded 4,000 of 40,001 trees (burnin fraction=0.10)
   Processing trees ('.' signifies 100 trees):

   ..................................................     5000   (# bip:    842    # topo:   4732)
   ..................................................    10000   (# bip:    920    # topo:   9483)
   ..................................................    15000   (# bip:    965    # topo:  14253)
   ..................................................    20000   (# bip:    989    # topo:  18985)
   ..................................................    25000   (# bip:   1022    # topo:  23720)
   ..................................................    30000   (# bip:   1037    # topo:  28433)
   ..................................................    35000   (# bip:   1063    # topo:  33175)
   ..........

   Maximum bipartition credibility tree written to /Users/bob/hiv.mbc
   Bipartition list written to /Users/bob/hiv.parts
   Tree probabilities written to /Users/bob/hiv.trprobs

   Done. 36,001 trees analyzed.
   Time spent: 0:00:20 (h:m:s)

   Max memory used: 1.77 GB.

   Different topologies seen:     34,127
   Different bipartitions seen:    1,064 (theoretical maximum: 1,262,699)
   Bipartitions in MBC tree:          37 (theoretical maximum: 37)
                                         (tree is fully resolved - no polytomies)

   MBC tree has not been explicitly rooted
   (Tree has been rooted at random internal node; root is at trifurcation)

   Highest Log Bipartition Credibility:  -28.96
```

### Example 3: 
#### Consensus tree with all compatible bipartitions, outgroup rooting

The command below causes sumt to do the following:

* Summarise the tree samples in the file `mhc_align.nexus.run1.t`
* `--all`: Compute majority rule consensus tree with all compatible bipartitions added (bipartitions with frequency < 50% are checked for compatitibiliy with tree iteratively, and added if possible, until tree is fully resolved or all bipartitions have been checked)
* `-b 0.1`: Discard 10% of tree samples as burn-in
* `-t 0.95`: Report 95% credible set of topologies (i.e., all the most frequently seen topologies to a cumulated probability of 95%)
* `-n`: Overwrite any existing output files with no warning
* `-v`: Print more verbose output to screen, including running count of distinct bipartitions and topologies seen in input trees
* `-r macaque olive_baboon yellow_baboon`: root consensus tree using outgroup consisting of the taxa "macaque", "olive_baboon", and "yellow_baboon".

```
sumt --all -b 0.1 -t 0.95 -nv -r macaque olive_baboon yellow_baboon mhc_align.nexus.run1.t
```

#### Screen output

```
   Counting trees in file mhc_align.nexus.run1.t         30,001

   Analyzing file: mhc_align.nexus.run1.t (Weight: 1.000)
   Discarded 3,000 of 30,001 trees (burnin fraction=0.10)
   Processing trees ('.' signifies 100 trees):

   ..................................................     5000   (# bip:     14    # topo:     12)
   ..................................................    10000   (# bip:     14    # topo:     13)
   ..................................................    15000   (# bip:     15    # topo:     14)
   ..................................................    20000   (# bip:     15    # topo:     16)
   ..................................................    25000   (# bip:     15    # topo:     16)
   ....................

   Consensus tree written to mhc.con
   Bipartition list written to mhc.parts
   Tree probabilities written to mhc.trprobs

   Done. 27,001 trees analyzed.
   Time spent: 0:00:02 (h:m:s)

   Max memory used: 32.92 MB.

   Different topologies seen:         16
   Different bipartitions seen:       15 (theoretical maximum: 96)
   Bipartitions in Consensus tree:     7 (theoretical maximum: 7)
                                         (tree is fully resolved - no polytomies)

   Consensus tree has been explicitly rooted
   (Root is at bifurcation)

   Log Bipartition Credibility:  -0.1825
```
