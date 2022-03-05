# sumt: command line program for computing consensus tree and other phylogenetic tree summaries

![](https://img.shields.io/badge/version-2.2.1-blue)

The command-line program `sumt` takes as input one or more files containing samples of phylogenetic trees (e.g., from a Bayesian MCMC analysis or a bootstrap procedure), and produces, as output, files containing (1) the majority-rule consensus tree with clade support values, (2) a summary of observed bipartitions along with branch length means and variances, and optionally (3) a list of tree topologies and how frequently they were observed. The name is taken from the `sumt` command in [MrBayes](https://nbisweden.github.io/MrBayes/index.html), whose functionality it also is meant to resemble. Clade support values and topology frequencies can be interpreted as posterior probabilities if the input trees are from a Bayesian MCMC analysis.


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
	* The file also contains a second consensus tree where branch labels = bipartition IDs, which can be used for interpreting bipartition file below.
	* File containing list of bipartitions present in input trees, along with mean and variance of corresponding branch lengths.
	* (Optionally) File containing list of observed tree topologies with posterior and cumulated probabilities
	* (Optionally) Progress indication is written to screen
* Reasonably fast and light on memory usage: 100,000 trees with 41 leaves processed in 1:18 (m:s), using max 4.9 GB memory on 2019 MacBook (99,034 topologies, 5,218 bipartitions)
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

* Summarise the tree samples in the files `primates.nexus.run1.t` and `primates.nexus.run2.t`
* `-b 0.25`: Discard 25% of tree samples as burn-in
* `-t 1.0`: Report 100% credible set of topologies (i.e., all topologies that were seen)
* `-s`: Compute average standard deviation of split frequencies as a measure of MCMC convergence
* `-f 0.1`: Report mean, variance, and standard error of the mean of branch lengths for all bipartitions seen in more than 10% of input trees

```
sumt -b 0.25 -t 1.0 -s -f 0.1 primates.nexus.run1.t primates.nexus.run2.t
```

#### Screen output

This is printed to screen during run:

```
   Counting trees in file 'primates.nexus.run1.t':          2,001
   Counting trees in file 'primates.nexus.run2.t':          2,001

   Analyzing file: primates.nexus.run1.t (Weight: 0.500)
   Discarded 500 of 2,001 trees (burnin fraction=0.25)
   Processing trees ('.' signifies 100 trees):

   ...............


   Analyzing file: primates.nexus.run2.t (Weight: 0.500)
   Discarded 500 of 2,001 trees (burnin fraction=0.25)
   Processing trees ('.' signifies 100 trees):

   ...............

   Average standard deviation of split frequencies: 0.000707

   Bipartition list written to primates.nexus.parts
   Consensus tree written to primates.nexus.con
   Tree probabilities written to primates.nexus.trprobs

   Done. 3,002 trees analyzed.
   Time spent: 0:00:00 (h:m:s)
```

#### Bipartition overview

This is the contents of the file `primates.nexus.parts`. Bipartitions are indicated using the "asterisk and dots" notation also used by e.g. MrBayes: Columns correspond to taxa (with column 1 = taxon 1). All taxa with "\*" (or ".") are in the same half of the bipartition. The ID column contains either a leaf name or a numerical branchID that corresponds to branch labels given in the second tree in the consensus tree file (in this case there are 5 leaves, and therefore only two internal branches).

```
List of bipartitions:

PART = Description of partition in .* format
PROB = Posterior probability of the partition
BLEN = Mean branch length
VAR  = Branch length variance
SEM  = Standard error of the mean for branch length
ID   = Leaf name or internal branch label, for those bipartitions that are included in consensus tree

PART    PROB      BLEN      VAR         SEM         ID
*....   1.000000  0.049286  (0.000068)  (0.000150)  Chimpanzee
.*...   1.000000  0.124501  (0.000181)  (0.000246)  Gibbon
..*..   1.000000  0.061740  (0.000085)  (0.000169)  Gorilla
...*.   1.000000  0.033342  (0.000055)  (0.000136)  Human
....*   1.000000  0.091739  (0.000134)  (0.000212)  Orangutan
.*..*   1.000000  0.055000  (0.000089)  (0.000173)  1
*.*..   0.943000  0.017117  (0.000026)  (0.000095)  2
```

#### Tree probabilities

This is the content of the file `primates.nexus.trprobs`. In this case there were only 5 leafs corresponding to a total of 15 possible trees (of which 3 were seen in the MCMC samples). Note: For data sets with more than about 15-20 taxa, each sampled tree will typically be unique and all topologies therefore have the same probability, meaning the credible set is not very useful. (Bipartitions on those trees will, however, not be unique, and clade probabilities carry useful information).

```
#NEXUS

[This file contains all trees that were found during the MCMC
search, sorted by posterior probability.
Lower case 'p' indicates the posterior probability of a tree.
Upper case 'P' indicates the cumulative posterior probability.]

begin trees;
   tree tree_1    [p = 0.943038] [P = 0.943038] = ((Orangutan,Gibbon),(Gorilla,Chimpanzee),Human);
   tree tree_2    [p = 0.056296] [P = 0.999334] = (((Orangutan,Gibbon),Gorilla),Human,Chimpanzee);
   tree tree_3    [p = 0.000666] [P = 1.000000] = (((Orangutan,Gibbon),Chimpanzee),Human,Gorilla);
end;
```

#### Consensus tree

This is the content of the file `primates.nexus.con`. The difference between the two trees is the information given as branch labels:

* First tree: labels are bipartition frequencies (= posterior probability of clade, if tree sample is from Bayesian MCMC analysis)
* Second tree: labels are the branchIDs also indicated in the bipartition summary in the file `primates.nexus.parts`. This should make it simpler to understand what branch the bipartition corresponds to (open the tree file in a treeviewer such as FigTree and view the branch labels).

```
#NEXUS

begin trees;
   [In this tree branch labels indicate the posterior probability of the bipartition corresponding to the branch.]
   tree prob = (Orangutan:0.0917385,Gibbon:0.124501,(Human:0.0333417,(Gorilla:0.06174,Chimpanzee:0.0492864)0.943:0.0171167)1.0:0.0549997);

   [In this tree branch labels indicate the bipartition ID listed in the file /Users/gorm/Documents/3_resources/example_data/primates.nexus.parts.
    These branch labels can be used for interpreting the table of branch lenght info in that file]
   tree partID = (Orangutan:0.0917385,Gibbon:0.124501,(Human:0.0333417,(Gorilla:0.06174,Chimpanzee:0.0492864)2:0.0171167)1:0.0549997);
end;
```

### Example 2: One large file from MCMC analysis, with verbose screen output

The command below causes sumt to do the following:

* Summarise the tree samples in the file `mhc.nexus.t`
* `-b 0.5`: Discard 50% of tree samples as burn-in
* `-t 0.75`: Report 75% credible set of topologies (i.e., all the most frequently seen topologies to a cumulated probability of 75%)
* `-f 0.1`: Report mean, variance, and standard error of the mean of branch lengths for all bipartitions seen in more than 10% of input trees
* `-n`: Overwrite any existing output files with no warning
* `-v`: Print more verbose output to screen, including running count of distinct bipartitions and topologies seen in input trees

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