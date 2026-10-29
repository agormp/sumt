# sumt

[![PyPI version](https://img.shields.io/pypi/v/sumt)](https://pypi.org/project/sumt/)
[![PyPI downloads](https://static.pepy.tech/personalized-badge/sumt?period=total&units=none&left_color=black&right_color=blue&left_text=PyPI%20downloads&service=github)](https://pepy.tech/project/sumt)

`sumt` computes summary trees and associated statistics from one or more files of phylogenetic trees.

Input trees are typically posterior samples from Bayesian MCMC analyses such as BEAST or MrBayes. They can also be other collections of plausible trees, such as bootstrap replicates.

## Installation

This README describes sumt 4.2.1, which requires Python 3.10 or newer and phylotreelib 3.x.

```bash
python3 -m pip install sumt
```

To upgrade an existing installation:

```bash
python3 -m pip install --upgrade sumt
```

## Quick start

Download the small BEAST example used below:

```bash
curl -L -o primate-mtDNA.trees \
  https://raw.githubusercontent.com/agormp/sumt/main/examples/data/primate-mtDNA.trees
```

### Consensus tree

Build a majority-rule consensus tree with mean bipartition lengths, discard 10% burn-in and midpoint-root the result:

```bash
sumt --con --biplen --rootmid -b 0.1 primate-mtDNA.trees
```

This writes `primate-mtDNA.con` in Nexus format. Use `--outformat newick` for Newick output or `--basename NAME` to choose another output prefix.

### MCC tree with common-ancestor heights and credible intervals

For rooted clock trees, build an MCC tree with common-ancestor heights and 80% and 95% credible intervals, using two worker processes:

```bash
sumt --mcc --caheight --ci 0.8,0.95 -b 0.1 --cpus 2 primate-mtDNA.trees
```

This writes `primate-mtDNA.mcc`. Credible intervals are central, equal-tailed intervals based on approximate quantiles. `--cik` controls their numerical resolution. See the [interval guide](https://github.com/agormp/sumt/blob/main/docs/recipes.md#credible-intervals---ci) for details.

### Combine independent runs

Input filenames are positional arguments. This example combines two MrBayes runs, discards 25% burn-in from each and reports the average standard deviation of split frequencies:

```bash
sumt --con --biplen --rootmid -b 0.25 -s mrbayes.1.t mrbayes.2.t
```

See the [recipes](https://github.com/agormp/sumt/blob/main/docs/recipes.md#downloading-example-data) for those example files and for per-file burn-in values.

## Choosing a summary

Choose one summary type and one branch-length/node-height method:

| Summary option | Result |
|---|---|
| `--con` | Majority-rule consensus |
| `--all` | Consensus with compatible bipartitions added |
| `--mcc` | Observed rooted topology with the highest clade credibility |
| `--mbc` | Observed unrooted topology with the highest bipartition credibility |
| `--hip` | Fully resolved rooted HIPSTR summary |
| `--mrhip` | HIPSTR favoring clades with more than 50% support, while allowing lower-support clades |

| Length/height option | Estimate |
|---|---|
| `--noblen` | Topology and support only |
| `--biplen` | Mean length of each matching bipartition |
| `--cladeheight` | Mean height where each exact clade occurs |
| `--caheight` | Mean MRCA height of each summary clade's tips across all input trees |

The height methods are intended for meaningfully rooted clock trees. The [option guide](https://github.com/agormp/sumt/blob/main/docs/recipes.md#what-the-branch-length--node-height-modes-mean) explains their assumptions and differences.

## Rooting and output

MCC retains the selected topology's root. HIPSTR and mrHIPSTR are rooted by construction. Use these methods with meaningfully rooted input trees.

Consensus and MBC output is unrooted unless a rooting option is requested. Choose midpoint rooting (`--rootmid`), minimum-variance rooting (`--rootminvar`) or an outgroup (`--rootog TAX[,TAX,...]`).

Nexus output declares `[&R]` for rooted results and `[&U]` for unrooted results. Summary statistics use seven significant digits. Input annotations are not summarized or copied to the output tree.

- `--ci 0.8,0.95` adds requested credible intervals.
- `--nometa` omits statistical annotations but retains the rootedness declaration.
- `--nolabel` omits internal support labels.
- `--noblen` omits branch lengths and their annotations.
- `--cpus 1` runs in one process. Set a larger value for parallel processing, or 0 for automatic selection.

Use `sumt --help` for all options. See the [recipes and option guide](https://github.com/agormp/sumt/blob/main/docs/recipes.md) for root credibility, topology-probability files, performance tuning and migration from earlier versions.

## Citation

If you use `sumt` in academic work, cite the GitHub repository using GitHub's “Cite this repository” link.
