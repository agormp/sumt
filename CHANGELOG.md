## Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

### [4.1.0] - 2026 March 23

#### Added
- Added `--cik K` to control the precision of the **approximate quantile** estimation used by `--ci` (credible intervals and median).
  Higher `K` gives finer resolution (more memory/CPU). See the README section **Credible intervals (`--ci`)** for details.

#### Changed
- Help text and README updated to explain that `--ci` endpoints are **estimated** (via a mergeable log-bucket histogram) and how `--cik` affects precision.

#### Fixed
- Ensured the `--cik` setting is consistently applied when computing common-ancestor depths (`--cadepth`) with multiprocessing, so partial results can be merged safely.

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
