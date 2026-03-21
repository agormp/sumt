# sumt

[![PyPI version](https://img.shields.io/pypi/v/sumt)](https://pypi.org/project/sumt/)
[![PyPI downloads](https://static.pepy.tech/personalized-badge/sumt?period=total&units=none&left_color=black&right_color=blue&left_text=PyPI%20downloads&service=github)](https://pepy.tech/project/sumt)

`sumt` computes summary trees and associated statistics from one or more files of phylogenetic trees.

Input trees are typically posterior samples from Bayesian MCMC analyses (BEAST / MrBayes) or bootstrap replicates. `sumt` is also useful for any collection of plausible trees, e.g. equally parsimonious trees, near-optimal ML trees, or trees from heuristic searches. In these cases, support values should be interpreted as frequencies within the provided tree set (what structure is “in common”), rather than posterior probabilities.

Supported summary trees:

- Majority-rule consensus (`--con`)
- Majority-rule consensus + all compatible bipartitions (`--all`)
- Maximum clade credibility (`--mcc`)
- Maximum bipartition credibility (`--mbc`)
- HIPSTR (`--hip`) and majority-rule HIPSTR (`--mrhip`)
- Works with any set of trees (posterior samples, bootstrap replicates, equally parsimonious trees, near-optimal trees). Support values are computed as frequencies across the input trees.

Branch-length / node-depth options:

- `--noblen` (topology + support only)
- `--biplen` (mean bipartition lengths)
- `--meandepth` (mean clade depths, then derive branch lengths)
- `--cadepth` (TreeAnnotator-style “--height ca”: mean MRCA depths, then derive branch lengths)

Rooting options:

- `--rootmid`, `--rootminvar`
- `--rootog TAX[,TAX,...]` or `--rootogfile FILE`
- `--rootcred` to compute root credibility on the output tree

---

## Version 4.0.0

Version 4 is a **major** release (breaking CLI changes, plus new capabilities).

- Now with multi-processing and computation of credible intervals
- For a high-level overview, see **[What changed (4.0.0)](#what-changed-400)**.
- For migration help, see **[Upgrading from 3.x → 4.x](#upgrading-from-3x--4x)**.

---

## Installation

```bash
python3 -m pip install sumt
```

Upgrade:

```bash
python3 -m pip install --upgrade sumt
```

---

## Downloading example data used in this README

This README uses small example files that live in the repository under `examples/data/`.

### Option A: curl

macOS and most Linux distributions include `curl`. Windows 10/11 typically includes `curl` as well.

```bash
# BEAST example (NEXUS tree file, 1000 trees, 12 tips)
curl -L -o primate-mtDNA.trees \
  https://raw.githubusercontent.com/agormp/sumt/main/examples/data/primate-mtDNA.trees

# MrBayes example (two runs; .t tree files)
curl -L -o mrbayes.1.t \
  https://raw.githubusercontent.com/agormp/sumt/main/examples/data/mrbayes.1.t
curl -L -o mrbayes.2.t \
  https://raw.githubusercontent.com/agormp/sumt/main/examples/data/mrbayes.2.t
```

### Option B: Python

```bash
python -c "import urllib.request as u; u.urlretrieve('https://raw.githubusercontent.com/agormp/sumt/main/examples/data/primate-mtDNA.trees','primate-mtDNA.trees')"
python -c "import urllib.request as u; u.urlretrieve('https://raw.githubusercontent.com/agormp/sumt/main/examples/data/mrbayes.1.t','mrbayes.1.t')"
python -c "import urllib.request as u; u.urlretrieve('https://raw.githubusercontent.com/agormp/sumt/main/examples/data/mrbayes.2.t','mrbayes.2.t')"
```

---

## Quick start

For additional short copy/paste patterns (including at least one example for every option), jump to **[Examples covering all main options](#examples-covering-all-main-options)**.

For full option documentation, run `sumt -h` (or `sumt --help`).

### 1) Consensus tree from a BEAST posterior tree file

```bash
# Majority-rule consensus tree, branch lengths by mean bipartition length, midpoint rooted
sumt --con --biplen --rootmid -b 0.1 primate-mtDNA.trees
```

This writes (by default) a summary tree file named `primate-mtDNA.con` (suffix = summary-tree type).
Use `--outformat newick` if you prefer Newick output.

### 2) Credible intervals + parallel processing

```bash
# 95% central CI for branch lengths; automatic CPU selection (0 = automatic)
sumt --con --biplen --ci 0.95 --cpus 0 primate-mtDNA.trees
```

- `--ci 0.95` computes a central 95% credible interval for each estimated branch length (or node depth)
- Optional: --cik K adjusts quantile-approximation precision (details below)
- `--cpus 0` chooses a default number of worker processes (use `--cpus 1` to force single-process, or specify an exact number of processes)

### 3) Multi-file workflow (MrBayes example: two independent runs)

This demonstrates:

- multiple input files
- burn-in via -b (one value for all files; optional comma-separated per-file values)
- optional ASDSF (average standard deviation of split frequencies) computation (`-s`) as a convergence diagnostic
- optional parallel processing (`--cpus`)

```bash
# Typical: same burn-in for both independent runs
sumt --con --biplen --rootmid -b 0.25 -s --cpus 0 \
  mrbayes.1.t mrbayes.2.t

# If you really need different burn-in per file, use comma-separated values:
sumt --con --biplen --rootmid -b 0.25,0.4 -s --cpus 0 \
  mrbayes.1.t mrbayes.2.t

# Add 80% and 95% credible intervals on branch lengths
sumt --con --biplen --rootmid -b 0.25,0.4 --ci 0.8,0.95 \
  mrbayes.1.t mrbayes.2.t
```

Output files are written with suffixes matching the summary-tree type (`.con`, `.mcc`, `.mbc`, `.hip`, `.mrhip`).
Use `--basename NAME` to control the output prefix; otherwise it uses the stem of the first input file.

---

## Terminology note: “node depth” vs “node height” (tree orientation)

In `sumt` (and the underlying `phylotreelib`), trees are typically treated as **rooted at the bottom with tips at the top**.
Accordingly, the code and output use:

- **node depth** = distance from the **tips (leaves)** back to a node (i.e., “time before the most recent leaf”)

In many programming libraries, the same quantity is called **node height** (because the root is drawn at the top).
When you see “height” in other tools (or in BEAST/TreeAnnotator options), it usually corresponds to what `sumt` calls “depth”. I am in the process of changing terminology to match that standard.

---

## What the different summary trees mean

All summary trees represent a single “best” topology derived from a set of input trees, with support values (and optionally branch-length / node-depth summaries). If your input trees are not posterior samples (e.g., equally parsimonious trees), the reported supports are empirical frequencies in your set (how often a split/clade occurs), not Bayesian posterior probabilities.

### `--con` Majority-rule consensus tree

Includes all bipartitions (splits) observed in **≥ 50%** of the post-burnin trees.
Support is typically reported as bipartition frequency (posterior probability under a Bayesian interpretation).

### `--all` Majority-rule consensus + all compatible bipartitions

Starts with the majority-rule consensus tree and then considers additional bipartitions in descending frequency order.
A bipartition is added if it is **compatible** with the current partially resolved tree. This continues until the tree is
fully resolved or no more compatible bipartitions remain.

### `--mcc` Maximum clade credibility tree

Selects an **observed input tree** (not a newly constructed consensus topology) that maximizes the product of clade frequencies
(equivalently: maximizes the sum of log clade frequencies). This is mainly meaningful for **rooted** (often clock-like) trees.

### `--mbc` Maximum bipartition credibility tree

Like MCC, but uses bipartitions instead of clades, and therefore **ignores rooting**. Two trees can share the same bipartitions
but differ in root position; MBC treats them as equivalent with respect to credibility.

### `--hip` HIPSTR and `--mrhip` majority-rule HIPSTR

[HIPSTR (Highest Independent Posterior SubTree reconstruction)](https://academic.oup.com/bioinformatics/article/41/10/btaf488/8250098) builds a fully resolved summary tree by choosing, at each internal node,
the child-clade pair with the highest combined posterior support. A HIPSTR tree is typically **not** an observed input tree.

- `--hip`: includes clades even if < 50% (yields a fully resolved tree under the HIPSTR heuristic)
- `--mrhip`: includes only clades with ≥ 50% support (majority rule)

---

## What the branch-length / node-depth modes mean

You always pick exactly one of: `--noblen`, `--biplen`, `--meandepth`, `--cadepth`.

### `--noblen` (no branch lengths)

Computes **topology + support only**. Branch lengths in the output are set to 0 (or omitted if you choose Newick without lengths).

Use this when:

- you only care about the consensus topology/support, or
- your input trees do not have meaningful lengths (e.g., pure topology samples).

### `--biplen` (mean bipartition lengths)

For each branch in the summary tree, `sumt` identifies the corresponding **leaf bipartition** (split) and sets the branch length to the
**mean length** of that bipartition across the post-burnin input trees (where that bipartition occurs).

This works for unrooted summaries (e.g., `--con`, `--all`, `--mbc`) and does not assume clock-like trees.

With `--ci`, credible intervals (and median) are computed for **branch lengths**.

### `--meandepth` (mean clade depths, then derive branch lengths)

Intended for **rooted clock-like trees** (e.g., time trees).

For each clade in the *summary* tree, `sumt` sets the node depth to the **mean node depth observed for that exact monophyletic clade**,
computed only across those input trees where the clade occurs as a monophyletic group.

Then branch lengths are derived from depths (`blen = depth(parent) - depth(child)`).

Notes:

- This can be based on very few observations for rare clades.
- It may produce **negative branch lengths** in some cases (a known issue with mean-depth approaches).

With `--ci`, credible intervals (and median) are computed for **node depths**.

### `--cadepth` (common-ancestor depths; TreeAnnotator-style “heights ca”)

Also intended for **rooted clock-like trees**.

For each clade in the summary tree, `sumt` computes, in **every** post-burnin input tree, the depth of the MRCA of that clade’s tip set,
and then takes the mean of those MRCA depths across all trees.

This corresponds to TreeAnnotator’s “heights ca” approach.

With `--ci`, credible intervals (and median) are computed for **node depths**.

---

## What the rooting options mean

You can choose at most one of: `--rootmid`, `--rootminvar`, `--rootog`, `--rootogfile`.

If you do not specify a rooting option:

- For `--mcc`, the chosen sample tree’s root is retained.
- For `--con`, `--all`, and `--mbc`, the output should be treated as *unrooted* unless you root it explicitly (e.g. with `--rootmid`, `--rootminvar`, or an outgroup).
- For `--hip` and `--mrhip`, the output is **rooted by construction** (HIPSTR is defined on clades / child-clade pairs, i.e. rooted structure). In practice this means the output root reflects the rooted structure present in the input trees; if your input trees are not meaningfully rooted, HIPSTR-style summaries are usually not appropriate.

### `--rootmid` (midpoint rooting)

Places the root at the midpoint of the tree’s diameter (the longest tip-to-tip path).
Often useful as a quick heuristic when no outgroup is available.

### `--rootminvar` (minimum-variance rooting)

Chooses a root location that minimizes the variance in root-to-tip distances (aiming for the most “clock-like” rooting).
This follows [Mai, Sayyari & Mirarab (2017), *PLOS ONE* 12(8)](https://doi.org/10.1371/journal.pone.0182238).

### `--rootog TAX[,TAX,...]` and `--rootogfile FILE` (outgroup rooting)

Roots the summary tree on an outgroup.

- `--rootog` takes a comma-separated list of taxa on the command line.
- `--rootogfile` takes a file with one taxon name per line.

If you provide multiple outgroup taxa, `sumt` attempts to place the root on the branch separating those taxa from the remaining tips.

---

## Root credibility (`--rootcred`)

`--rootcred` annotates *branches* in the output tree with how often the root was observed to fall on that branch among the input trees.

There are two cases:

1) **If you use an outgroup** (`--rootog` / `--rootogfile`):

   - For each input tree, `sumt` identifies the branch where the outgroup attaches.
   - Root credibility for a branch is the fraction of input trees where the outgroup attached there.
   - This works even if the input trees are unrooted, because the outgroup attachment defines a rooting.

2) **If you do not use an outgroup**:

   - `sumt` assumes the input trees are already rooted (typical for BEAST clock trees).
   - It tracks the observed root location (root bipartition) directly across trees.

Why the “cumulated root credibility” may be < 100%:

- Root credibility is only reported on branches that exist in the final summary topology. If some root locations occur on branches (bipartitions) that are *not present* in the summary tree, their mass is not represented.

Example:

```bash
sumt --con --biplen --rootmid --rootcred primate-mtDNA.trees
```

---

## Credible intervals (`--ci`)

When you request `--ci`, `sumt` estimates **central credible intervals** for either:

- branch lengths (when using `--biplen`), or
- node depths (when using `--meandepth` or `--cadepth`)

Implementation note (approximate quantiles)

- Credible intervals are based on **estimated** quantiles (not exact order statistics).
- Internally, `sumt` uses a mergeable **log-bucket histogram** (see `QuantileAccumulator` in `phylotreelib`) to approximate quantiles in one pass.
- The approximation precision is controlled by `--cik K`, which sets the histogram resolution (2^K mantissa sub-bins per exponent bucket):
  - higher K ⇒ finer resolution but more memory/CPU
  - worst-case relative bucket midpoint error bound is about `2^-(K+1)` (e.g. K=7 ≈ 0.39%, K=8 ≈ 0.20%, K=9 ≈ 0.10%)
- Default is `--cik 7`.

Examples:

```bash
# One CI (default precision, --cik 7):
sumt --con --biplen --ci 0.95 primate-mtDNA.trees

# Several CIs:
sumt --con --biplen --ci 0.5,0.8,0.95 primate-mtDNA.trees

# Increase precision of the quantile approximation:
sumt --con --biplen --ci 0.95 --cik 9 primate-mtDNA.trees
```

---

## Examples covering all main options

Below are short “pattern” examples. They are deliberately small and repetitive so you can copy/paste a working starting point.

### Input/output controls

```bash
# Input format (autodetection usually works; use only if needed)
sumt --informat nexus --con --biplen primate-mtDNA.trees
sumt --informat newick --con --biplen mytrees.newick

# Output format
sumt --outformat newick --con --biplen primate-mtDNA.trees
sumt --outformat nexus  --con --biplen primate-mtDNA.trees

# Basename controls output stem
sumt --basename primates_summary --con --biplen primate-mtDNA.trees

# Suppress metacomments in NEXUS output
sumt --con --biplen --nometa primate-mtDNA.trees

# Overwrite without prompting; quiet mode implies -n
sumt -n --con --biplen primate-mtDNA.trees
sumt -q --con --biplen primate-mtDNA.trees

# Verbose tracebacks (useful for debugging)
sumt -v --con --biplen primate-mtDNA.trees
```

### Summary-tree types

```bash
sumt --con   --biplen primate-mtDNA.trees
sumt --all   --biplen primate-mtDNA.trees
sumt --mcc   --meandepth primate-mtDNA.trees
sumt --mbc   --biplen primate-mtDNA.trees
sumt --hip   --biplen primate-mtDNA.trees
sumt --mrhip --biplen primate-mtDNA.trees
```

### Branch-length / node-depth settings

```bash
# Topology + support only
sumt --con --noblen primate-mtDNA.trees

# Mean bipartition lengths
sumt --con --biplen primate-mtDNA.trees

# Mean clade depths (clock-like rooted trees)
sumt --mcc --meandepth primate-mtDNA.trees

# Common-ancestor depths (TreeAnnotator-style --height ca)
sumt --mcc --cadepth primate-mtDNA.trees
```

### Rooting and root credibility

```bash
# Midpoint / min-variance
sumt --con --biplen --rootmid    primate-mtDNA.trees
sumt --con --biplen --rootminvar primate-mtDNA.trees

# Outgroup on command line (comma-separated list; not biologically meaningful...)
sumt --con --biplen --rootog Macaca_fuscata,M._mulatta,M._fascicularis,M._sylvanus primate-mtDNA.trees

# Outgroup from file (one taxon per line)
printf "Macaca_fuscata\nM._mulatta\nM._fascicularis\nM._sylvanus\n" > outgroup.txt
sumt --con --biplen --rootogfile outgroup.txt primate-mtDNA.trees

# Root credibility on the output tree
sumt --con --biplen --rootmid --rootcred primate-mtDNA.trees
```

### Bayesian / diagnostic options

```bash
# Burn-in (one value for all files)
sumt --con --biplen -b 0.25 primate-mtDNA.trees

# Burn-in (one value per file, comma-separated; could be different per file)
sumt --con --biplen -b 0.25,0.25 mrbayes.1.t mrbayes.2.t

# Tree probabilities and credible set of topologies
sumt --con --biplen -t 0.95 primate-mtDNA.trees

# ASDSF across files (+ minimum frequency threshold)
sumt --con --biplen -s -f 0.1 mrbayes.1.t mrbayes.2.t

# Credible intervals with higher quantile-precision
sumt --con --biplen --ci 0.8,0.95 --cik 9 primate-mtDNA.trees
```

Notes on `-t PROB` (`.trprobs` output):

- When you run with `-t` (e.g. `-t 0.95`), `sumt` writes a file named `<basename>.trprobs`.
- The file is **NEXUS** format (`#NEXUS … begin trees; … end;`) containing the most probable topologies up to the
  requested cumulative probability (a “credible set”).
- Each tree entry is annotated with `p` (posterior probability of that topology) and `P` (cumulative probability so far).
- This is most useful for **small** numbers of taxa. With ~15–20+ taxa, almost every sampled tree is typically unique,
  so the “credible set” becomes long and less informative even though split/clade supports remain useful.


### Performance tuning

`sumt` can process large tree files faster by splitting work across multiple processes (`--cpus`).
Internally, each worker processes trees in **chunks** (`--chunksize`, default: 250 trees per chunk).
A larger chunk size reduces scheduling/serialization overhead, but can increase peak memory usage and sometimes
reduce load balancing across CPUs (especially if trees vary in processing cost).


```bash
# Parallel processing
sumt --con --biplen --cpus 8 primate-mtDNA.trees

# Force single-process
sumt --con --biplen --cpus 1 primate-mtDNA.trees

# Chunk size: larger reduces overhead, may increase memory usage
sumt --con --biplen --cpus 8 --chunksize 500 primate-mtDNA.trees
```

---

## Version 4.0.0 notes

### What changed (4.0.0)

This is a **major** release because it includes intentional breaking CLI changes:

- **Removed `-i` / `-w` / `--autow`**: input files are now positional, and file weights are no longer supported (each tree counts equally).
- **Added credible intervals** via `--ci`.
- **Added multiprocessing** via `--cpus` and `--chunksize`.
- **Simplified output formats**: `--outformat` is now `newick` or `nexus`; use `--nometa` to suppress metacomments.

### Upgrading from 3.x → 4.x

#### Input files

Old (3.x):

```bash
sumt ... -i file1.t -i file2.t
# or
sumt ... -w 0.5 file1.t -w 1.0 file2.t --autow
```

New (4.x):

```bash
sumt ... file1.t file2.t
```

Notes:

- Weighting was removed
- Filenames are now a positional argument (no longer requiring -i FILENAME)

#### Burn-in syntax

Old (3.x) accepted space-separated lists:

```bash
sumt ... -b 0.25 0.4 -i file1.t -i file2.t
```

New (4.x) uses one value or comma-separated values:

```bash
sumt ... -b 0.25,0.4 file1.t file2.t
```

#### Credible intervals (new)

```bash
sumt ... --ci 0.8,0.95
```

#### Multiprocessing (new)

```bash
sumt ... --cpus 8 --chunksize 500
```

---

## Citation

If you use `sumt` in academic work, the simplest option is to cite the GitHub repository (GitHub “Cite this repository” in the sidebar).

---

## Notes

- Some combinations (especially `--meandepth` / `--cadepth`) assume clock-like, rooted trees.
- Large tree files can be processed efficiently, but for best performance you may want to tune `--chunksize` and `--cpus`.
