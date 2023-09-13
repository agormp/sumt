# sumt

![](https://img.shields.io/badge/version-3.6.0-blue)
[![PyPI downloads](https://static.pepy.tech/personalized-badge/sumt?period=total&units=none&left_color=black&right_color=blue&left_text=PyPI%20downloads&service=github)](https://pepy.tech/project/sumt)

The command-line program `sumt` computes consensus trees and other tree-summary statistics for sets of phylogenetic trees. The input trees can be in one or more input files, and will typically be from a Bayesian MCMC analysis (BEAST or MrBayes for instance) or from a bootstrap procedure. 

`sumt` can compute four different kinds of main tree summaries:

* Majority rule consensus tree
* Majority rule consensus tree, with all compatible bipartitions added
* Maximum clade credibility tree
* Maximum bipartition credibility tree

Branch labels on these trees indicate clade support. 

`sumt` also produces a summary of observed bipartitions along with branch length means and variances and, optionally, a list of tree topologies and how frequently they were observed. The name is taken from the `sumt` command in [MrBayes](https://nbisweden.github.io/MrBayes/index.html), whose functionality it was originally meant to resemble. Clade support values and topology frequencies can be interpreted as posterior probabilities if the input trees are from a Bayesian MCMC analysis.


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
	* File containing summary tree with clade support values (= frequency of bipartition in input trees). 
	    * The summary tree can be one of these:
		    * Majority rule consensus tree
		    * Majority rule consensus tree, with all compatible bipartitions added
			* Maximum clade credibility (MCC) tree 
		    * Maximum bipartition credibility (MBC) tree - this is similar to MCC, but ignoring location of root.
	    * The file also contains a second consensus tree where branch labels indicate bipartition IDs, which can be used for interpreting bipartition file below.
	* File containing list of bipartitions (in "*." format) present in input trees, along with mean and variance of corresponding branch lengths. This list includes both bipartitions that correspond to branches in the summary tree, and bipartitions not included in the summary tree (e.g., low frequency bipartitions).
	* During run: progress bar showing percentage of file analyzed
	* (Optionally) File containing list of observed tree topologies with posterior and cumulated probabilities
* Optimized for speed and memory usage:
	* Consensus tree computed from 100,000 trees with 41 leaves in 28 s, using max 49 MB memory on 2021 MacBook (4,440 distinct bipartitions seen)
	* Same file processed in 30 s, using max 2.83 GB memory when also keeping track of topologies (74,283 distinct topologies seen)
* Option to discard fraction of trees as burn-in (for Bayesian analyses)
* Option to compute average standard deviation of split frequencies when multiple input files are given. This can be used as a measure of convergence of Bayesian analyses, assuming that different files represent independent MCMC chains.
* Option to include all compatible bipartitions in consensus tree (in addition to those that are present in more than 50% of input trees).
* Options to root consensus tree using either outgroup, midpoint, [minimum variance](https://pubmed.ncbi.nlm.nih.gov/28800608/) rooting, or based on where root is most frequently placed in input tree sample
* Option to set node depths to mean of those observed in input trees (only useful when input trees are based on clock model)
* Option to assign specific weights to different input files.
* Option to automatically assign weights so all files have equal impact regardless of number of trees in them.
* Option to set basename of output files (default: basename will be stem of input file name)
* Program prints information about run to stdout:
	* Number of leaves on tree
	* Number of different toplogies seen in input trees (if also using option -t)
	* Number of bipartitions or clades seen in input trees, along with theoretical maximum
	* Number of bipartitions in summary tree, along with theoretical maximum
	* Indication of whether summary tree is resolved
	* Indication of whether summary tree has been explicitly rooted and how
	* Log credibility (sum of logs of bipartition or clade frequencies in summary tree)
	* Indication of how branch lengths have been set
	* Root credibility (when input trees are clock trees)
	* Report of maximum memory usage during processing

## Usage

```
usage: sumt [-h] [--version] [--con | --all | --mcc | --mbc]
            [--rootmid | --rootminvar | -r TAXON [TAXON ...] | --rootfile FILE | --rootmaxfreq]
            [-b NUM] [-t NUM] [-s] [-f NUM] [-n] [-v] [-q] [--basename NAME] [--meandepth]
            [--autow] [--informat FORMAT] [-i FILE | -w WEIGHT FILE]

Computes summary tree and statistics from set of phylogenetic trees

options:
  -h, --help            show this help message and exit
  --version             show the program's version number and exit

Type of summary tree (pick one option):
  --con                 majority rule consensus tree
  --all                 majority rule consensus tree with all compatible bipartitions added
  --mcc                 Maximum Clade Credibility (MCC) tree. The MCC tree is determined by
                        inspecting tree samples and selecting the tree that has the highest
                        product of clade frequencies (= highest sum of log of clade frequencies).
                        The MCC tree is therefore a tree that has been observed in the pool of
                        tree samples, differing from the consensus tree which typically does not
                        match any individual sample. NOTE 1: only meaningful if input trees are
                        estimated using clock model. NOTE 2: by default, the MCC tree uses the
                        rooting of the specific tree sample. This will often (but not always)
                        correspond to the bipartition where the root is most commonly found in the
                        input trees.
  --mbc                 Maximum Bipartition Credibility (MBC) tree. The MBC tree is similar to the
                        MCC tree but counting bipartitions instead of clades, i.e. ignoring
                        rooting (two input trees can have the same set of bipartitions, but be
                        rooted in different locations).

Rooting of summary tree:
  --rootmid             perform midpoint rooting of summary tree
  --rootminvar          perform minimum variance rooting of summary tree
  -r TAXON [TAXON ...]  root summary tree on specified outgroup taxon/taxa
  --rootfile FILE       root summary tree on outgroup taxa listed in file (one name per line)
  --rootmaxfreq         root summary tree on bipartition where root is located most frequently in
                        input trees. NOTE: only meaningful if input trees are estimated using
                        clock model

Bayesian phylogeny options:
  -b NUM                burnin: fraction of trees to discard [0 - 1; default: 0.0]
  -t NUM                compute tree probabilities, report NUM percent credible interval [0 - 1]
  -s                    compute average standard deviation of split frequencies (ASDSF)
  -f NUM                Minimum frequency for including bipartitions in report and in computation
                        of ASDSF [default: 0.1]

Output to terminal and files:
  -n                    no warning when overwriting files
  -v                    verbose: show full traceback in the event of failed python execution
  -q                    quiet: don't print progress indication to terminal window. NOTE: also
                        turns on the -n option
  --basename NAME       base name of output files (default: derived from input file)

Estimation of node depths for clock trees:
  --meandepth           set node depth for each clade to mean node depth observed for that clade
                        among input trees (and branch lengths are then based on these depths).
                        NOTE 1: only meaningful if input trees are estimated using clock model.
                        NOTE 2: will only work if all clades in tree have been observed at least
                        once among input trees - the option will therefore fail for some rootings.

Other options:
  --autow               automatically assign file weights based on tree counts, so all files have
                        equal impact (default is for all trees, not files, to be equally
                        important)
  --informat FORMAT     format of input files: nexus, newick [default: nexus]

Input tree files:
  -i FILE               input FILE(s) containing phylogenetic trees (repeat -i FILE option for
                        each input file)
  -w WEIGHT FILE        input FILEs with specified weights (repeat -w WEIGHT FILE option for each
                        input file)
```

## Usage examples

### Example 1: 
#### Majority rule consensus tree, bipartition summary, topology summary, computation of average standard deviation of split frequencies, midpoint rooting

The command below causes `sumt` to do the following:

* `--con`: Compute majority rule consensus tree (this is default and could have been omitted)
* `-b 0.25`: Discard 25% of tree samples as burn-in
* `-t 0.99`: Keep track of topology probabilities, report 99% credible set
* `-s`: Compute average standard deviation of split frequencies as a measure of MCMC convergence (asdsf)
* `-f 0.1`: Include bipartitions seen in more than 10% of input trees for computations of (1) asdsf and of (2) branch lengt mean, variance, and standard error of the mean
* `--rootmid`: Perform midpoint rooting
* `-i primates.nexus.run1.t -i primates.nexus.run2.t `: Summarise the tree samples in the files `primates.nexus.run1.t` and `primates.nexus.run2.t`

```
sumt --con -b 0.25 -t 0.99 -f 0.1 --rootmid -i primates.nexus.run1.t -i primates.nexus.run2.t 
```

#### Screen output

This is printed to screen during run:

```
   Counting trees in file primates.nexus.run1.t                             2,501
   Counting trees in file primates.nexus.run2.t                             2,501

   Analyzing file: primates.nexus.run1.t (Weight: 0.500)
   Discarded 625 of 2,501 trees (burnin fraction=0.25)

   Processing trees:
   
   0      10      20      30      40      50      60      70      80      90     100
   v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v
   *********************************************************************************


   Analyzing file: primates.nexus.run2.t (Weight: 0.500)
   Discarded 625 of 2,501 trees (burnin fraction=0.25)

   Processing trees:
   
   0      10      20      30      40      50      60      70      80      90     100
   v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v
   *********************************************************************************

   Consensus tree written to primates.con
   Bipartition list written to primates.parts
   Tree probabilities written to primates.trprobs

   Done. 3,752 trees analyzed.
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
*....   1.000000  0.07601    (0.0003773)  (0.0003171)  Chimpanzee
.*...   1.000000  0.3582     (0.01896  )  (0.002248 )  Gibbon
..*..   1.000000  0.07316    (0.0005376)  (0.0003785)  Gorilla
...*.   1.000000  0.05669    (0.0002765)  (0.0002715)  Human
....*   1.000000  0.2833     (0.009456 )  (0.001588 )  Orangutan
.*..*   1.000000  0.1454     (0.003596 )  (0.000979 )  6
*..*.   0.969350  0.03699    (0.0003814)  (0.0003238)  7
```

#### Tree probabilities

This is the content of the file `primates.trprobs`. In this case there were only 5 leaves corresponding to a total of 15 possible trees, of which 3 were seen in the MCMC samples. Note: For data sets with more than about 15-20 taxa, each sampled tree will typically be unique and all topologies therefore have the same probability, meaning the credible set is not very useful. (Bipartitions on those trees will, however, not be unique, and clade probabilities carry useful information).

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
    tree tree_1 [p = 0.969350] [P = 0.969350] = (((5,2),3),4,1);
    tree tree_2 [p = 0.016258] [P = 0.985608] = (((5,2),1),4,3);
    tree tree_3 [p = 0.014392] [P = 1.000000] = ((5,2),(1,3),4);
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
   tree prob = ((((Human:0.0566895,Chimpanzee:0.0760097)0.969:0.0369884,Gorilla:0.0731637)1.000:0.145374,Orangutan:0.283342)1.000:0.0748381,Gibbon:0.35818);

   [In this tree branch labels indicate the bipartition ID listed in the file primates.parts.
    These branch labels can be used for interpreting the table of branch lenght info in that file]
   tree partID = ((((Human:0.0566895,Chimpanzee:0.0760097)7:0.0369884,Gorilla:0.0731637)6:0.145374,Orangutan:0.283342)Gibbon:0.0748381,Gibbon:0.35818);
end;
```

### Example 2: 
#### Maximum bipartition credibility tree, topology summary, verbose output, setting basename of output files

The command below causes sumt to do the following:

* `--mbc`: Compute maximum bipartition credibility tree (instead of majority rule consensus)
* `-b 0.1`: Discard 10% of tree samples as burn-in
* `-t 0.75`: Report 75% credible set of topologies (i.e., all the most frequently seen topologies to a cumulated probability of 75%)
* `-n`: Overwrite any existing output files with no warning
* `-v`: Print more verbose output to screen, including running count of distinct bipartitions and topologies seen in input trees
* `--basename /Users/bob/hiv`: produce output files with the indicated stem (/Users/bob/hiv.parts, /Users/bob/hiv.trprobs, /Users/bob/hiv.mbc)
* `-i gp120.nexus.trees`: Summarise the tree samples in the file `gp120.nexus.trees`

```
sumt --mbc -b 0.1 -t 0.75 -nv --basename /Users/bob/hiv -i gp120.nexus.trees
```

#### Screen output

This is printed to screen during run. The progress bar gradually fills up over the 20 seconds it takes to process the 36,001 post-burnin trees.

At the end of the run the actual and theoretical maximum for number of bipartitions in the MBC tree is reported. The total number of observed topologies and bipartitions (and the theoretical maximal possible number of bipartitions) is also reported. In this case we have seen 34,127 distint tree topologies among the 36,001 trees analyzed, meaning that about 95% of the sampled trees are unique. In these topologies we have seen 1,064 distinct bipartitions. If the 34,127 topologies had been completely different (in the sense of not sharing any bipartitions), then the number of distinct bipartitions would have been 1,262,699 (so 1,064 is a small fraction of that, indicating that some bipartitions are observed in a large fraction of the tree samples).

The ouput also indicates that the summary tree is fully resolved (has no polytomies). This will always be the case for a maximum bipartition credibility tree (because individual tree samples are typically fully resolved). It is further stated that the MBC tree has not been explicitly rooted (none of `sumt`'s rooting options were used).

Finally the Highest Log Bipartition Credibility is output (this is the sum of the logs of the bipartition frequencies, for those bipartitions that are present in the MBC tree).

```
   Counting trees in file gp120.DNA.align.nexus.run1.t                     40,001

   Analyzing file: gp120.DNA.align.nexus.run1.t (Weight: 1.000)
   Discarded 4,000 of 40,001 trees (burnin fraction=0.10)

   Processing trees:
   
   0      10      20      30      40      50      60      70      80      90     100
   v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v
   *********************************************************************************

   Maximum bipartition credibility tree written to bobhiv.mbc
   Bipartition list written to bobhiv.parts
   Tree probabilities written to bobhiv.trprobs

   Number of leaves on tree:          40
   Different topologies seen:     34,127
   Different bipartitions seen:    1,064 (theoretical maximum: 1,262,699)
   Bipartitions in MBC tree:          37 (theoretical maximum: 37)
                                         (tree is fully resolved - no polytomies)

   MBC tree has not been explicitly rooted
   (Tree has been rooted at random internal node; root is at trifurcation)

   Highest Log Bipartition Credibility:  -28.96

   Done. 36,001 trees analyzed.
   Time spent: 0:00:20 (h:m:s)
   Max memory used: 1.50 GB.
```

### Example 3: 
#### Consensus tree with all compatible bipartitions, outgroup rooting

The command below causes sumt to do the following:

* `--all`: Compute majority rule consensus tree with all compatible bipartitions added (bipartitions with frequency < 50% are checked for compatitibiliy with tree iteratively in order of decreasing frequencies, and added if possible. Iteration stops when the consensus tree is fully resolved or all bipartitions have been checked)
* `-b 0.1`: Discard 10% of tree samples as burn-in
* `-t 0.95`: Report 95% credible set of topologies (i.e., all the most frequently seen topologies to a cumulated probability of 95%)
* `-n`: Overwrite any existing output files with no warning
* `-v`: Print more verbose output to screen, including running count of distinct bipartitions and topologies seen in input trees
* `-r macaque olive_baboon yellow_baboon`: root consensus tree using outgroup consisting of the taxa "macaque", "olive_baboon", and "yellow_baboon".
* `-i mhc_align.nexus.run1.t`: Summarise the tree samples in the file `mhc_align.nexus.run1.t`

```
sumt --all -b 0.1 -t 0.95 -nv -r macaque olive_baboon yellow_baboon -i mhc_align.nexus.run1.t
```

#### Screen output

```

   Counting trees in file mhc_align.nexus.run1.t                           30,001

   Analyzing file: mhc_align.nexus.run1.t (Weight: 1.000)
   Discarded 3,000 of 30,001 trees (burnin fraction=0.10)

   Processing trees:
   
   0      10      20      30      40      50      60      70      80      90     100
   v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v
   *********************************************************************************

   Consensus tree written to mhc_align.con
   Bipartition list written to mhc_align.parts
   Tree probabilities written to mhc_align.trprobs

   Number of leaves on tree:           9
   Different topologies seen:         16
   Different bipartitions seen:       15 (theoretical maximum: 96)
   Bipartitions in Consensus tree:     6 (theoretical maximum: 6)
                                         (tree is fully resolved - no polytomies)

   Consensus tree has been explicitly rooted
   (Root is at bifurcation)

   Log Bipartition Credibility:  -0.1825

   Done. 27,001 trees analyzed.
   Time spent: 0:00:02 (h:m:s)
   Max memory used: 30.67 MB.
```
