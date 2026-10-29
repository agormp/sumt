## Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

### [4.2.1] - 2026 October 29

#### Fixed
- Require phylotreelib 3.0.1 for the Tree iteration and cache fixes.

---

### [4.2.0] - 2026 October 28

#### Added
- Nexus summary output declares `[&R]` for MCC/HIPSTR/mrHIPSTR and explicitly rooted consensus/MBC trees, or `[&U]` otherwise. Declarations remain with `--nometa`.

#### Fixed
- Reject the unfinished `--usemedian` option instead of reporting median estimates while computing means.
- Source distributions include the public documentation. README documentation links also work on PyPI.
- Corrected MCC, HIPSTR and mrHIPSTR descriptions in help text, status messages and the README.
- Memory monitoring stops promptly instead of delaying command exit until the next sample.
- `.trprobs` entries declare `[&R]` for clade topology probabilities, including when both topology types are tracked, or `[&U]` for bipartition-only probabilities. This is independent of main-summary rooting.
- Counting, burn-in and parallel processing handle comments, Nexus declarations and translation tables consistently.

#### Changed
- Reduced redundant rootedness test combinations while retaining workflow and output checks.
- Summary output retains requested interval annotations and seven-significant-digit precision with the new phylotreelib API. `--noblen` omits branch annotations, while support labels remain separately controlled.
- Requires Python 3.10 or newer and phylotreelib 3.x (`>=3.0.0,<4`).
- Shortened the README and moved detailed options, examples and older migration notes to `docs/recipes.md`.

---

### [4.1.5] - 2026 September 22

#### Fixed
- Added the option `--cadepth` as a hidden, deprecated compatibility alias for `--caheight`; using it now emits a deprecation warning.

---

### [4.1.4] - 2026 September 22

#### Fixed
- Updated the phylotreelib requirement to >= version 2.3.1. This fixes handling of situations where some nodes (typically leaves) have effectively constant node-height samples. These will now use the exact observed height for the median and credible-interval endpoints, rather than the midpoint of a log bucket.

---

### [4.1.3] - 2026 March 30

#### Fixed
- Fixed argparse bug introduced in 4.1.1 that caused sumt to crash.

---

### [4.1.2] - 2026 March 30

#### Changed
- Updated internal calls to use the renamed `CAHeightEstimator` class from phylotreelib 2.1.2 (no change to CLI behaviour; requires phylotreelib ≥ 2.1.2).

---

### [4.1.1] - 2026 March 29

#### Changed
- Changed terminology from "depth" to "height" in help text and option names, to describe difference between date of most recent leaf to given other node (trees now perceived as being rooted at top, with leaves below, as in computer science). The canonical options became `--cladeheight` (formerly `--meandepth`) and `--caheight` (formerly `--cadepth`).

#### Deprecated
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
