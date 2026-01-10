# sumt

![](https://img.shields.io/badge/version-3.12.2-blue)
[![PyPI downloads](https://static.pepy.tech/personalized-badge/sumt?period=total&units=none&left_color=black&right_color=blue&left_text=PyPI%20downloads&service=github)](https://pepy.tech/project/sumt)
[![DOI](https://zenodo.org/badge/461522623.svg)](https://zenodo.org/doi/10.5281/zenodo.10148693)

The command-line program `sumt` computes consensus trees and other tree-summary statistics for sets of phylogenetic trees. The input trees can be in one or more input files, and will typically be from a Bayesian MCMC analysis (BEAST or MrBayes for instance) or from a bootstrap procedure. 

`sumt` can compute six different kinds of main tree summaries:

* Maximum clade credibility (MCC) tree
* HIPSTR (highest independent posterior subtree reconstruction)
* mrHIPSTR (Majority rule HIPSTR tree)
* Majority rule consensus tree
* Majority rule consensus tree, with all compatible bipartitions added
* Maximum bipartition credibility (MBC) tree

Branch lengths in summary tree can be set based on:

* Mean node depths across input trees (for clock trees)
* "Common Ancestor" node depths (for clock trees)
* Mean length of bipartitions (splits) across input trees across input trees

Rooting can be done using:

* Outgroup
* Midpoint
* Minimum variance
* Root of maximum credibility tree (for MCC)
* Root for most credible resolution of root bipartition (for HIPSTR and mrHIPSTR)

The name is taken from the `sumt` command in [MrBayes](https://nbisweden.github.io/MrBayes/index.html), whose functionality it was originally meant to resemble. Clade support values and topology frequencies can be interpreted as posterior probabilities if the input trees are from a Bayesian MCMC analysis.


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

## Citation

To cite sumt: use the link in the right sidebar under About --> Cite this repository.

## Dependencies

`sumt` relies on the [phylotreelib library](https://github.com/agormp/phylotreelib) for phylogeny-related matters, and on [psutil](https://pypi.org/project/psutil/) for (optionally) monitoring memory usage. These are automatically included when using pip to install.

## Overview

* Input:
	* One or more files containing phylogenetic trees (all trees must have same leaf names), in NEXUS or Newick format.
	* Typically trees are from a Bayesian MCMC analysis, but could also be from a bootstrap procedure
	* Can read sample files from BEAST and MrBayes
* Output:
	* File containing summary tree with clade support values (= frequency of bipartition or clade in input trees). 
	    * The summary tree can be one of these:
		    * Majority rule consensus tree
		    * Majority rule consensus tree, with all compatible bipartitions added
			* Maximum clade credibility (MCC) tree 
		    * Maximum bipartition credibility (MBC) tree - this is similar to MCC, but ignoring location of root.
			* HIPSTR tree ([highest independent posterior subtree reconstruction](https://academic.oup.com/bioinformatics/article/41/10/btaf488/8250098))
			* mrHIPSTR (majority rule HIPSTR)
		* By default, the summary tree contains [NEXUS metacomments](https://beast.community/nexus_metacomments) summarizing e.g. clade or bipartition frequencies, and mean and standard deviation for branch lengths or depths.
	* During run: progress bar showing percentage of file analyzed
	* (Optionally) File containing list of observed tree topologies with posterior and cumulated probabilities
* Optimized for speed and memory usage
* Option to discard fraction of trees as burn-in (for Bayesian analyses)
* Option to compute average standard deviation of split frequencies when multiple input files are given. This can be used as a measure of convergence of Bayesian analyses, assuming that different files represent independent MCMC chains.
* Option to include all compatible bipartitions in consensus tree (in addition to those that are present in more than 50% of input trees).
* Options to root consensus tree using either outgroup, midpoint, [minimum variance](https://pubmed.ncbi.nlm.nih.gov/28800608/) rooting, or based on root of maximum credibility tree (for MCC), or root for most credible resolution of root bipartition (for HIPSTR and mrHIPSTR)
* Option to set node depths to mean of those observed in input trees (only useful when input trees are based on clock model)
* Option to set node depths using "common ancestor depths" in input trees (same as "treeannotator -heights ca" in the BEAST package; only useful when input trees are based on clock model)
* When setting mean or CA node depths: depths of leaves will be the mean of those observed for leaves (only relevant when estimating some tip dates)
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
	* Report of time and maximum memory usage during processing

## Usage

```
usage: sumt [-h] [--version] [--informat FORMAT] [--outformat FORMAT]
            [-i FILE | -w WEIGHT FILE] [--autow] [--basename NAME] [-n] [-v] [-q]
            (--con | --all | --mcc | --mbc | --hip | --mrhip)
            (--noblen | --biplen | --meandepth | --cadepth)
            [--rootmid | --rootminvar | --rootog TAX [TAX ...] | --rootogfile FILE]
            [--rootcred] [-b NUM [NUM ...]] [-t NUM] [-s] [-f NUM]

Computes summary tree and statistics from set of phylogenetic trees

options:
  -h, --help            show this help message and exit
  --version             show the program's version number and exit

INPUT AND OUTPUT:
  --informat FORMAT     format of input tree files: newick, nexus [default: nexus]
  --outformat FORMAT    format of output tree file: newick, nexus, mcnexus [default:
                        mcnexus]. mcnexus: nexus with metacomments as used by BEAST (e.g.,
                        [&height=100.0])
  -i FILE               input FILE(s) containing phylogenetic trees (repeat -i FILE option
                        for each input file)
  -w WEIGHT FILE        input FILEs with specified weights (repeat -w WEIGHT FILE option for
                        each input file)
  --autow               automatically assign file weights based on tree counts, so all files
                        have equal impact (default is for all trees, not files, to be equally
                        important)
  --basename NAME       base name of output files (default: derived from input file)
  -n                    no warning when overwriting files
  -v                    verbose: show full traceback in the event of failed python execution
  -q                    quiet: don't print progress indication to terminal window. NOTE: also
                        turns on the -n option

TYPE OF SUMMARY TREE (pick one option):
  --con                 majority rule consensus tree
  --all                 majority rule consensus tree with all compatible bipartitions added
  --mcc                 Maximum Clade Credibility (MCC) tree. The MCC tree is determined by
                        inspecting tree samples and selecting the tree that has the highest
                        product of clade frequencies (= highest sum of log of clade
                        frequencies). The MCC tree is therefore a tree that has been observed
                        in the pool of tree samples, differing from the consensus tree which
                        typically does not match any individual sample. NOTE 1: only
                        meaningful if input trees are estimated using clock model or
                        otherwise rooted. NOTE 2: by default, the MCC tree uses the rooting
                        of the specific tree sample. This will often (but not always)
                        correspond to the bipartition where the root is most commonly found
                        in the input trees.
  --mbc                 Maximum Bipartition Credibility (MBC) tree. The MBC tree is similar
                        to the MCC tree but counting bipartitions instead of clades, i.e.
                        ignoring rooting (two input trees can have the same set of
                        bipartitions, but be rooted in different locations).
  --hip                 HIPSTR summary tree (Highest Independent Posterior SubTree). (see
                        Baele et al., Bioinformatics, 2025, 41(10), btaf488)The tree is built
                        by choosing, at each internal node, the child clade pair with the
                        highest combined posterior support, producing a fully resolved
                        summary tree not necessarily observed among the input trees.
  --mrhip               MrHIPSTR summary tree (majority rule HIPSTR tree). Like HIPSTR, but
                        only including clades with >= 50% support

ESTIMATION OF BRANCH LENGTHS (pick one option):
  --noblen              Do not set branch lengths (only the topology and branch- or clade-
                        support of the summary tree are estimated).
  --biplen              Set branch lengths in summary tree based on average for corresponding
                        leaf bipartitions:each branch in tree corresponds to a bipartition of
                        the leaves into two groups. Branch lenghts are set to the mean of the
                        length of thecorresponding bipartition across all input trees.
  --meandepth           set node depth for each clade to mean node depth observed for that
                        clade among input trees (and branch lengths are then based on these
                        depths). Warning: option is intended for input trees estimated using
                        a clock model. It requires that all clades in the summary tree have
                        been observed in the input trees, and may fail for some
                        rootings.NOTE: mean is computed across trees where the specific,
                        monophyletic clade is present, and may therefore be based on very few
                        (down to one) values. NOTE 2: may result in negative branch lengths.
  --cadepth             'Common Ancestor depth'. Same as option '--height ca' in
                        treeannotator. Uses all trees in input set when determining node-
                        depths. For a given clade: (1) Find the most recent common ancestor
                        of the leaves in that clade in each of the input trees. (2) Compute
                        node-depth of clade as the mean of the depths of these MRCAs. This is
                        different from --meandepth where only trees with that exact clade are
                        included when computing the mean. Warning: option is intended for
                        input trees estimated using a clock model. It requires that all
                        clades in the summary tree have been observed in the input trees, and
                        may fail for some rootings.

ROOTING OF SUMMARY TREE:
  --rootmid             perform midpoint rooting of summary tree
  --rootminvar          perform minimum variance rooting of summary tree
  --rootog TAX [TAX ...]
                        root summary tree on outgroup; specify outgroup taxon/taxa on
                        command-line
  --rootogfile FILE     root summary tree on outgroup; specify outgroup taxon/taxa in file
                        (one name per line)
  --rootcred            compute root credibility for all possible rooting locations and add a
                        'rootcred' attribute to branches in the summary tree. If an outgroup
                        is specified: track which branch (bipartition) the outgroup attaches
                        to in each input tree and report the frequency of rooting there. If
                        no outgroup is specified: assume input trees are rooted and track
                        root frequencies on branches. The cumulated root credibility may be
                        less than 100% if some root locations are not present in the summary
                        tree.

BAYESIAN PHYLOGENY OPTIONS:
  -b NUM [NUM ...]      burnin: fraction of trees to discard [0 - 1; default: [0]]. Either
                        one value (used on all input files), or one value per input file.
  -t NUM                compute tree probabilities, report NUM percent credible interval [0 -
                        1]
  -s                    compute average standard deviation of split frequencies (ASDSF)
  -f NUM                Minimum frequency for including bipartitions in computation of ASDSF
                        [default: 0.1]
```

## Usage examples

### Example 1: 
#### Majority rule consensus tree, bipartition summary, topology summary, computation of average standard deviation of split frequencies, midpoint rooting

The command below causes `sumt` to do the following:

* `--con`: Compute majority rule consensus tree
* `-b 0.25`: Discard 25% of tree samples as burn-in
* `-t 0.99`: Keep track of topology probabilities, report 99% credible set
* `-s`: Compute average standard deviation of split frequencies as a measure of MCMC convergence (asdsf)
* `-f 0.1`: Include bipartitions seen in more than 10% of input trees for computations of (1) asdsf and of (2) branch lengt mean, variance, and standard error of the mean
* `--rootmid`: Perform midpoint rooting
* `--biplen`: Set branch lengths to mean of those observed for the corresponding bipartitions in input trees
* `-i primates.run1.t -i primates.run2.t `: Summarise the tree samples in the files `primates.run1.t` and `primates.run2.t`

```
sumt --con -b 0.25 -t 0.99 -f 0.1 --rootmid --biplen -i primates.run1.t -i primates.run2.t
```

#### Screen output

This is printed to screen during run:

```
   Counting trees in file primates.run1.t                                      2,001
   Counting trees in file primates.run2.t                                      2,001
   
   Analyzing file: primates.run1.t (Weight: 1.000)
   Discarded 500 of 2,001 trees (burnin fraction=0.25)
   
   Processing trees:
   0      10      20      30      40      50      60      70      80      90     100
   v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v
   *********************************************************************************   
   
   Analyzing file: primates.run2.t (Weight: 1.000)
   Discarded 500 of 2,001 trees (burnin fraction=0.25)
   
   Processing trees:
   0      10      20      30      40      50      60      70      80      90     100
   v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v
   *********************************************************************************   
   
   Computing summary tree...done
   Tree probabilities written to primates.trprobs
   
   Number of leaves on input trees:       5
   Different topologies seen:             3
   Different bipartitions seen:           4 (theoretical maximum: 6)
   Bipartitions in Consensus tree:        2 (theoretical maximum: 2)
                                            (tree is fully resolved - no polytomies)
   
   Branch lengths set based on mean branch lengths for corresponding bipartitions
   
   Consensus tree has been midpoint-rooted
   
   Log bipartition credibility:  -0.028723
   
   Done. 3,002 trees analyzed.
   Time spent: 0:00:00 (h:m:s)
   Max memory used: 40.02 MB.
   
   Consensus tree written to primates.con
```

#### Tree probabilities

This is the content of the file `primates.trprobs`. In this case there were only 5 leaves corresponding to a total of 15 possible trees, of which 3 were seen in the MCMC samples. Note: For data sets with more than about 15-20 taxa, each sampled tree will typically be unique and all topologies therefore have the same probability, meaning the credible set is not very useful. (Bipartitions or clades on those trees will, however, not be unique, and clade probabilities carry useful information).

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
    tree tree_1 [p = 0.971686] [P = 0.971686] = ((3,(5,2)),1,4);
    tree tree_2 [p = 0.014990] [P = 0.986676] = ((5,2),(3,1),4);
    tree tree_3 [p = 0.013324] [P = 1.000000] = (((5,2),1),3,4);
end;
```

#### Consensus tree

This is the content of the file `primates.con`. NEXUS metacomments [&...] are used to summarise information across input trees. Depending on options there may be NEXUS metacomments for nodes (placed before the colons) and/or branches (placed after colons). In this particular case the summary tree contains the following attributes for each branch:

* bipartition_cred: credibility of this branch (bipartition) = frequency of bipartition in post-burnin trees
* length: mean branch length, taken across post-burnin trees where this bipartition was observed
* length_var and length_sd: variance and standard deviation of branch length across the same input trees

NEXUS metacomments can be visualised by some tree-viewers (e.g., figtree)

```
#NEXUS

begin trees;
	tree nexus_tree = ((((Human:0.05688946[&bipartition_cred=1,length=0.05688946,length_sd=0.01657206,length_var=0.000274633],Chimpanzee:0.07560542[&bipartition_cred=1,length=0.07560542,length_sd=0.01924216,length_var=0.0003702608])0.9716855:0.03674312[&bipartition_cred=0.9716855,length=0.03674312,length_sd=0.0200375,length_var=0.0004015013],Gorilla:0.07280193[&bipartition_cred=1,length=0.07280193,length_sd=0.0234786,length_var=0.0005512446])1:0.1464279[&bipartition_cred=1,length=0.1464279,length_sd=0.06030987,length_var=0.003637281],Orangutan:0.2792791[&bipartition_cred=1,length=0.2792791,length_sd=0.09485017,length_var=0.008996555])1:0.3550858[&bipartition_cred=1,length=0.3550858,length_sd=0.1376393,length_var=0.01894458],Gibbon:0.3550858[&bipartition_cred=1,length=0.3550858,length_sd=0.1376393,length_var=0.01894458]);
end;
```

### Example 2: 
#### Maximum bipartition credibility tree, topology summary,  setting basename of output files

The command below causes sumt to do the following:

* `--mbc`: Compute maximum bipartition credibility tree. The MBC tree is determined by inspecting tree samples and selecting the tree that has the highest product of bipartition frequencies (= highest sum of log of bipartition frequencies). This is similar to the MCC (maximum clade credibility) but counting bipartitions instead of monophyletic clades (i.e., ignoring rooting).
* `-b 0.1`: Discard 10% of tree samples as burn-in
* `-t 0.75`: Report 75% credible set of topologies (i.e., all the most frequently seen topologies to a cumulated probability of 75%)
* `-n`: Overwrite any existing output files with no warning
* `--basename /Users/bob/hiv`: produce output files with the indicated stem (/Users/bob/hiv.parts, /Users/bob/hiv.trprobs, /Users/bob/hiv.mbc)
* `--biplen`: Set branch lengths to mean of those observed for the corresponding bipartitions in input trees
* `-i gp120.trees`: Summarise the tree samples in the file `gp120.trees`

```
sumt --mbc -b 0.1 -t 0.75 -n --basename /Users/bob/hiv --biplen -i gp120.trees
```


### Example 3: 
#### Consensus tree with all compatible bipartitions, outgroup rooting

The command below causes sumt to do the following:

* `--all`: Compute majority rule consensus tree with all compatible bipartitions added (bipartitions with frequency < 50% are checked for compatitibiliy with tree iteratively in order of decreasing frequencies, and added if possible. Iteration stops when the consensus tree is fully resolved or all bipartitions have been checked)
* `-b 0.1`: Discard 10% of tree samples as burn-in
* `-t 0.95`: Report 95% credible set of topologies (i.e., all the most frequently seen topologies to a cumulated probability of 95%)
* `-n`: Overwrite any existing output files with no warning
* `--rootog macaque olive_baboon yellow_baboon`: root consensus tree using outgroup consisting of the taxa "macaque", "olive_baboon", and "yellow_baboon".
* `--biplen`: Set branch lengths to mean of those observed for the corresponding bipartitions in input trees
* `-i mhc_align.run1.t`: Summarise the tree samples in the file `mhc_align.run1.t`

```
sumt --all -b 0.1 -t 0.95 -n --rootog macaque olive_baboon yellow_baboon --biplen -i mhc_align.run1.t
```

### Example 4: 
#### Maximum clade credibility tree with "common ancestor" node depths, separate burnins

The command below causes sumt to do the following:

* `--mcc`: Compute maximum clade credibility tree. Note: input trees need to be based on a clock model for this to be meaningful. 
* `-b 0.25 0.4`: Discard 25% of tree samples in first file, and 40% of trees in second file as burn-in
* `--cadepth`: set node depth for each clade to mean node depth observed for MRCA of that clade among input trees (this is the same as `treeannotator -heights ca` in the BEAST package). Note: only meaningful when input trees are clock trees.
* `-n`: Overwrite any existing output files with no warning
* `-s`: Compute average standard deviation of clade frequencies as a measure of MCMC convergence
* `--basename beast_summary`: produce output files with the indicated stem 
* `-i beastrun1.trees -i beastrun2.trees`: Summarise the tree samples in the files `beastrun1.trees` and `beastrun2.trees`

```
sumt --mcc -b 0.25 0.4 --cadepth -ns --basename beast_summary -i beastrun1.trees -i beastrun2.trees 
```

#### Screen output

This is printed to screen during run:

```

   Counting trees in file beastrun1.trees                                      2,001
   Counting trees in file beastrun2.trees                                      2,001
   
   Analyzing file: beastrun1.trees (Weight: 1.000)
   Discarded 500 of 2,001 trees (burnin fraction=0.25)
   
   Processing trees:
   0      10      20      30      40      50      60      70      80      90     100
   v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v
   *********************************************************************************   
   
   Analyzing file: beastrun2.trees (Weight: 1.000)
   Discarded 800 of 2,001 trees (burnin fraction=0.40)
   
   Processing trees:
   0      10      20      30      40      50      60      70      80      90     100
   v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v
   *********************************************************************************   
   
   Computing summary tree...done
   
   Number of leaves on input trees:     508
   Different clades seen:           176,642 (theoretical maximum: 1,369,914)
   Bipartitions in MCC tree:            505 (theoretical maximum: 505)
                                            (tree is fully resolved - no polytomies)
   
   Branch lengths set based on common ancestor depths in input trees
   
   MCC tree rooted at original root of tree sample having highest clade credibility
   
   Highest log clade credibility:  -1245.29
   Average standard deviation of split frequencies: 0.023403
   
   Done. 2,702 trees analyzed.
   Time spent: 0:00:18 (h:m:s)
   Max memory used: 746.58 MB.
   
   Maximum clade credibility tree written to beast_summary.mcc
```

### Example 5: 
#### Maximum clade credibility tree with mean node depths

The command below causes sumt to do the following:

* `--mcc`: Compute maximum clade credibility tree. Note: input trees need to be based on a clock model for this to be meaningful. 
* `-b 0.25`: Discard 25% of tree samples as burn-in
* `--meandepth`: set node depth for each clade to mean node depth observed for that specific, monophyletic clade among input trees (this is the same as `treeannotator -heights mean` in the BEAST2 package). Note: only meaningful when input trees are clock trees.
* `-n`: Overwrite any existing output files with no warning
* `-s`: Compute average standard deviation of clade frequencies as a measure of MCMC convergence
* `--basename beast_summary`: produce output files with the indicated stem 
* `-i beastrun1.trees -i beastrun2.trees`: Summarise the tree samples in the files `beastrun1.trees` and `beastrun2.trees`

```
sumt --mcc -b 0.25 --meandepth -ns --basename beast_summary -i beastrun1.trees -i beastrun2.trees 
```

### Example 7: 
#### mrHIPSTR  tree with common ancestor node depths and root credibility

The command below causes sumt to do the following:

* `--mrhip`: Compute mrHIPSTR tree
* `-b 0.25`: Discard 25% of tree samples as burn-in
* `--cadepth`: set node depth for each clade to mean node depth observed for MRCA of that clade among input trees (this is the same as `treeannotator -heights ca` in the BEAST package). Note: only meaningful when input trees are clock trees. Node depths of leaves will be set to the mean leaf depth observed (only relevant when leaf dates are being estimated)
* `-n`: Overwrite any existing output files with no warning
* `-s`: Compute average standard deviation of clade frequencies as a measure of MCMC convergence
* `--basename beast_summary`: produce output files with the indicated stem 
* `--rootcred`: compute credibility of root location on summary tree; this is the fraction of post-burnin trees in which the used root position was observed. Also reported: the cumulated root credibility on the summary tree (may be less than 100% if some root splits are not in final summary).
* `-i beastrun1.trees -i beastrun2.trees`: Summarise the tree samples in the files `beastrun1.trees` and `beastrun2.trees`
```
sumt --mrhip -b 0.25 --cadepth -ns --basename beast_summary --rootcred -i beastrun1.trees -i beastrun2.trees 
```

#### Screen output

This is printed to screen during run:

```
   Counting trees in file beastrun1.trees                                      2,001
   Counting trees in file beastrun2.trees                                      2,001
   
   Analyzing file: beastrun1.trees (Weight: 1.000)
   Discarded 500 of 2,001 trees (burnin fraction=0.25)
   
   Processing trees:
   0      10      20      30      40      50      60      70      80      90     100
   v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v
   *********************************************************************************   
   
   Analyzing file: beastrun2.trees (Weight: 1.000)
   Discarded 500 of 2,001 trees (burnin fraction=0.25)
   
   Processing trees:
   0      10      20      30      40      50      60      70      80      90     100
   v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v
   *********************************************************************************   
   
   Computing summary tree...done
   
   Number of leaves on input trees:     508
   Different clades seen:           193,339 (theoretical maximum: 1,522,014)
   Bipartitions in mrHIPSTR tree:       505 (theoretical maximum: 505)
                                            (tree is fully resolved - no polytomies)
   
   Branch lengths set based on common ancestor depths in input trees
   
   mrHIPSTR tree rooted at most frequently observed root bipartition
   Root credibility (frequency of root bipartition in input trees):       86.3%
   Cumulated root credibility (sum of rootcred for all branches in tree): 100.0%
   
   Log clade credibility:  -1014.47
   Average standard deviation of split frequencies: 0.022985
   
   Done. 3,002 trees analyzed.
   Time spent: 0:00:22 (h:m:s)
   Max memory used: 2.30 GB.
   
   Majority rule HIPSTR (mrHIPSTR) tree written to beast_summary.mrhip
```

### Example 8: 
#### HIPSTR  tree with mean  node depths and root credibility

The command below causes sumt to do the following:

* `--hip`: Compute HIPSTR tree
* `-b 0.25`: Discard 25% of tree samples as burn-in
* `--meandepth`: set node depth for each clade to mean node depth observed for that clade among input trees (this is the same as `treeannotator -heights mean` in the BEAST package). Note: only meaningful when input trees are clock trees. Node depths of leaves will be set to the mean leaf depth observed (only relevant when leaf dates are being estimated)
* `-n`: Overwrite any existing output files with no warning
* `-s`: Compute average standard deviation of clade frequencies as a measure of MCMC convergence
* `--basename beast_summary`: produce output files with the indicated stem 
* `--rootcred`: compute credibility of root location on summary tree; this is the fraction of post-burnin trees in which the used root position was observed. Also reported: the cumulated root credibility on the summary tree (may be less than 100% if some root splits are not in final summary).
* `-i beastrun1.trees -i beastrun2.trees`: Summarise the tree samples in the files `beastrun1.trees` and `beastrun2.trees`
```
sumt --hip -b 0.25 --meandepth -ns --basename beast_summary --rootcred -i beastrun1.trees -i beastrun2.trees 
```

#### Screen output

This is printed to screen during run:

```

   Counting trees in file beastrun1.trees                                      2,001
   Counting trees in file beastrun2.trees                                      2,001
   
   Analyzing file: beastrun1.trees (Weight: 1.000)
   Discarded 500 of 2,001 trees (burnin fraction=0.25)
   
   Processing trees:
   0      10      20      30      40      50      60      70      80      90     100
   v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v
   *********************************************************************************   
   
   Analyzing file: beastrun2.trees (Weight: 1.000)
   Discarded 500 of 2,001 trees (burnin fraction=0.25)
   
   Processing trees:
   0      10      20      30      40      50      60      70      80      90     100
   v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v
   *********************************************************************************   
   
   Computing summary tree...done
   
   Number of leaves on input trees:     508
   Different clades seen:           193,339 (theoretical maximum: 1,522,014)
   Bipartitions in HIPSTR tree:         505 (theoretical maximum: 505)
                                            (tree is fully resolved - no polytomies)
   
   Branch lengths set based on mean node depths in input trees
   
   HIPSTR tree rooted at most frequently observed root bipartition
   Root credibility (frequency of root bipartition in input trees):       86.3%
   Cumulated root credibility (sum of rootcred for all branches in tree): 100.0%
   
   Log clade credibility:  -1013.9
   Average standard deviation of split frequencies: 0.022985
   
   Done. 3,002 trees analyzed.
   Time spent: 0:00:14 (h:m:s)
   Max memory used: 2.34 GB.
   
   HIPSTR tree written to beast_summary.hip
```