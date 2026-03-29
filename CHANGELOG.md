## Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

### Unreleased (on GitHub, not yet in PyPi version)

#### Added
- Added option --usemedian: specifies the use of median as pointestimate for branch length / node-height options (default is to use mean).

---

### [4.1.1] - 2026 March 29

#### Changed
- Changed terminology from "depth" to "height" in help text and option names, to describe difference between date of most recent leaf to given other node (trees now perceived as being rooted at top, with leaves below, as in computer science).

### Deprecated
- Deprecated --meandepth option. Use --cladeheight instead; --meandepth will continue working for now, with deprecation warning. (First: terminology changed from depth to height. Secondly: user can in an upcoming version choose whether to use mean or median as point estimate when estimating branch lengths using any of the options in the BRANCHLENGTH group).


---

### [4.1.0] - 2026 March 23

#### Added
- Added `--cik K` to control the precision of the **approximate quantile** estimation used by `--ci` (credible intervals and median).
  Higher `K` gives finer resolution (more memory/CPU). See the README section **Credible intervals (`--ci`)** for details.

#### Changed
- Help text and README updated to explain that `--ci` endpoints are **estimated** (via a mergeable log-bucket histogram) and how `--cik` affects precision.

#### Fixed
- Ensured the `--cik` setting is consistently applied when computing common-ancestor depths (`--cadepth`) with multiprocessing, so partial results can be merged safely.

---

### [4.0.0] - 2026 March 20

#### Added
- Multiprocessing for tree processing via `--cpus` and `--chunksize`.
- Credible intervals (`--ci`) and median estimates for branch lengths / node depths.
- Optional peak memory monitoring and reporting (process tree RSS + system available memory at peak).
- Downloadable small example datasets referenced from the README (`examples/data/...`).
- Clearer, more structured CLI help text (`sumt -h`) and expanded README explanations and copy/paste recipes.

#### Changed
- **Breaking CLI change:** input tree files are now **positional arguments** (no longer `-i FILE`).
- Burn-in syntax changed: `-b` now accepts a single value or **comma-separated** per-file values (instead of space-separated lists).
- Output format options simplified: `--outformat` now uses `newick` or `nexus` (and Nexus metacomments are controlled by `--nometa`).
- Output/annotation model simplified: the former `mcnexus` mode is replaced by standard `nexus` output with optional metacomments (`--nometa` to suppress).

#### Removed
- File weighting support (`-w`, `--autow`); all trees now contribute equally.
- `-i` (repeated input flag); input files are now positional arguments.
- `mcnexus` output mode.

#### Fixed
- Various bugfixes and speedups
