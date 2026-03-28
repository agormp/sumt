import phylotreelib as pt
import argparse, os, sys, time, math, copy, psutil, statistics, configparser, threading
from importlib.metadata import version, PackageNotFoundError
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from dataclasses import dataclass
from itertools import (takewhile,repeat,islice)
from operator import itemgetter
from pathlib import Path
import gc
gc.disable()        # Faster. Assume no cyclic references will ever be created

# CA-depth worker globals (set by _ca_worker_init)
_CA_PLAN = None
_CA_TRACKCI = False
_CA_QUANTILE_K = 7

####################################################################################

def main(commandlist=None):
    start=time.time()
    mem_mon = PeakMemoryMonitor(interval_s=1)
    mem_mon.start()
    args = parse_commandline(commandlist)

    try:
        output = OutputManager(args)
        setup_output_directory(args.outbase)
        n_trees_analyzed, count_burnin_filename_list = count_trees(args, output)

        if args.cpus == 1:
            treesummarylist = process_trees(count_burnin_filename_list, args, output)
            worker_pids = None
        else:
            treesummarylist, worker_pids = process_trees_concurrent(
                                count_burnin_filename_list, args, output, n_trees_analyzed)

        ave_std = compute_converge_stats(treesummarylist, args) if args.std else None
        treesummary = merge_treesummaries(treesummarylist)
        sumtree = compute_sumtree(treesummary, args, count_burnin_filename_list, output, n_trees_analyzed, worker_pids=None)

        if args.treeprobs:
            trproblist = compute_trprobs(treesummary, args)
            trprobs_status_message = print_trprobs(treesummary, trproblist, args, output)
            output.info(trprobs_status_message)

        print_result_summary(sumtree, treesummary, start, n_trees_analyzed,
                             ave_std, args, output, mem_mon, worker_pids)
        print_sumtree(sumtree, treesummary, args, output)

    except Exception as error:
        handle_error(error, args.verbose)

####################################################################################

@dataclass
class MemPeak:
    peak_total_rss_bytes = 0
    peak_available_bytes_at_peak = 0
    sample_count = 0

####################################################################################

class PeakMemoryMonitor:
    """
    Cross-platform (Linux/macOS/Windows) peak memory monitor.

    Measures:
      - peak total RSS across this process + all descendants
      - peak system memory used percent (overall machine pressure)

    Notes:
      - RSS sums may double-count shared pages (libraries, CoW)
      - Sampling interval trades overhead vs accuracy.
    """

    def __init__(self, interval_s: float = 0.2):
        self.interval_s = float(interval_s)
        self._stop = threading.Event()
        self._thread = None
        self.peak = MemPeak()

    def start(self):
        if self._thread is not None:
            return
        self._thread = threading.Thread(target=self._run, name="PeakMemoryMonitor", daemon=True)
        self._thread.start()

    def stop(self):
        if self._thread is None:
            return self.peak
        self._stop.set()
        self._thread.join()
        self._thread = None
        return self.peak

    def _run(self):
        root_pid = os.getpid()
        while not self._stop.is_set():
            total_rss = 0

            # System-level pressure (what you ultimately care about for “crash risk”)
            try:
                vm = psutil.virtual_memory()
                avail = int(vm.available)
            except Exception:
                avail = 0

            # Process tree RSS (parent + all descendants)
            try:
                root = psutil.Process(root_pid)
                procs = [root] + root.children(recursive=True)
            except Exception:
                procs = []

            for p in procs:
                try:
                    total_rss += p.memory_info().rss
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass

            if total_rss > self.peak.peak_total_rss_bytes:
                self.peak.peak_total_rss_bytes = total_rss
                self.peak.peak_available_bytes_at_peak = avail

            self.peak.sample_count += 1
            time.sleep(self.interval_s)

####################################################################################

def _get_pid(_=None):
    return os.getpid()

####################################################################################

def available_logical_cpus():
    n = psutil.cpu_count(logical=True)
    return n if n else 1

####################################################################################

def cpu_process_summary_line(processes_used):
    total = available_logical_cpus()
    pct = (processes_used / total * 100.0) if total else 0.0
    return f"Number of processors used: {processes_used} / {total} ({pct:.0f}%)"

####################################################################################

def parse_comma_separated_floats(parser, value, optname):
    vals = []
    for part in value.split(","):
        part = part.strip()
        if not part:
            continue
        try:
            vals.append(float(part))
        except ValueError:
            parser.error(f"{optname}: not a float: {part}")
    if not vals:
        parser.error(f"{optname}: no values provided")
    return vals

####################################################################################

def parse_comma_separated_strings(parser, value, optname):
    vals = [part.strip() for part in value.split(",") if part.strip()]
    if not vals:
        parser.error(f"{optname}: no values provided")
    return vals

####################################################################################

def parse_commandline(commandlist):
    # Python note: "commandlist" is to enable unit testing of argparse code
    # Will be "None" when run in script mode, and argparse will then automatically
    # take values from sys.argv[1:]
    # https://jugmac00.github.io/blog/testing-argparse-applications-the-better-way/

    parser = build_parser()
    args = parser.parse_args(commandlist)

    # Check for deprecated command --meandepth (now --cladedepth)
    if commandlist is None:
        commandlist = sys.argv
    if "--meandepth" in commandlist:
        print("WARNING: --meandepth is deprecated; use --cladedepth instead.", file=sys.stderr)

    # If output basename is not set: use stem of first input filename minus all suffixes
    if not args.outbase:
        infilepath = args.infilelist[0]
        args.outbase = Path(infilepath.stem.split(".")[0])

    # Parse burnin fractions from comma-separated string
    args.burninfrac = parse_comma_separated_floats(parser, args.burninfrac, "-b")

    if len(args.burninfrac) == 1:
        burnin_value = args.burninfrac[0]
        args.burninfrac = [burnin_value] * len(args.infilelist)
    elif len(args.burninfrac) != len(args.infilelist):
        parser.error("either provide one burnin value, or one value per input file")

    if any(x < 0 or x > 1 for x in args.burninfrac):
        parser.error("option -b: NUM must be between 0.0 and 1.0")

    # Parse CI probabilities from comma-separated string
    if args.ci is None:
        args.ci_probs = None
    else:
        args.ci_probs = parse_comma_separated_floats(parser, args.ci, "--ci")
        for p in args.ci_probs:
            if not (0.0 < p < 1.0):
                parser.error(f"--ci: prob must be in (0,1): {p}")
        args.ci_probs = sorted(set(args.ci_probs))   # deduplicate + sort for stable output

    if args.treeprobs and (args.treeprobs > 1 or args.treeprobs < 0):
        parser.error(f"option -t: NUM must be between 0.0 and 1.0 (provided value: {args.treeprobs})")

    nfiles = len(args.infilelist)
    if args.std and nfiles == 1:
        parser.error("cannot compute standard deviation (option -s) from one tree file")

    if args.quiet:
        args.nowarn = True

    if args.ogfile:
        args.outgroup = read_outgroup(args.ogfile)
    elif args.outgroup:
        args.outgroup = parse_comma_separated_strings(parser, args.outgroup, "--rootog")

    args.actively_rooted = (args.outgroup or args.rootmid or args.rootminvar)

    if (args.treetype in ["mcc", "hip", "mrhip"]) and args.actively_rooted:
        parser.error(
            f"Rooting method should not be specified for {args.treetype} trees "
            "(input trees are already assumed to be rooted)"
        )

    # Bipartitions need to be tracked in these situations
    args.trackbips = (
        args.treetype in ["con", "all", "mbc"]
        or args.biplen
    )

    # Clades need to be tracked in these situations
    args.trackclades = (
        args.treetype in ["mcc", "hip", "mrhip"]
        or args.cladedepth
        or args.cadepth
    )

    # Root needs to be tracked if rootcred==True, OR if using biplen with MCC
    args.trackroot = (
        args.rootcred
        or (args.biplen and (args.treetype == "mcc"))
    )

    # Branch lengths need to be tracked if biplen==True
    args.trackblen = args.biplen

    # Root-branch length fractions are only needed for MCC + biplen
    args.trackrootblen = (
        args.treetype == "mcc"
        and args.biplen
    )

    # Node depths need to be tracked if cladedepth==True
    args.trackdepth = args.cladedepth

    # Subclade-pairs need to be tracked if using HIPSTR or MrHIPSTR
    args.track_subcladepairs = args.treetype in ("hip", "mrhip")

    # Topologies need to be tracked if computing MCC or MBC trees, or tree probabilities
    args.tracktopo = (
        args.treetype in ("mcc", "mbc")
        or args.treeprobs
    )

    # Quantiles need to be tracked if one or more credible intervals are requested
    args.trackci = bool(args.ci_probs)

    # Check that quantile accumulator bin resolution is in reasonable range
    if not (1 <= args.ci_k <= 20):
        parser.error("--cik must be in range [1,20] (there will be 2^k bins in quantile estimator!)")

    if args.cpus < 0:
        parser.error("--cpus must be >= 0")
    if args.chunksize <= 0:
        parser.error("--chunksize must be >= 1")

    if args.cpus == 0:
        ncpu = psutil.cpu_count(logical=True) or 1
        args.cpus = ncpu

    return args

####################################################################################

def build_parser():

    try:
        sumt_version = version("sumt")
    except PackageNotFoundError:
        sumt_version = "Unknown"

    parser = argparse.ArgumentParser(description = "Computes summary tree and statistics from set of phylogenetic trees")

    parser.add_argument( "--version", action="version", version=f"%(prog)s {sumt_version}",
        help="show program version and exit")

    ####################################################################################

    inout_grp = parser.add_argument_group("INPUT / OUTPUT")

    inout_grp.add_argument("--informat", dest="informat", action="store", metavar="FORMAT",
        choices=["newick", "nexus"], default="nexus",
        help="input tree format: %(choices)s [default: %(default)s]")

    inout_grp.add_argument("--outformat", dest="outformat", action="store", metavar="FORMAT",
        choices=["newick", "nexus"], default="nexus",
        help="output tree format: %(choices)s [default: %(default)s]")

    inout_grp.add_argument("--nometa", action="store_true",
        help="omit Nexus metacomments with node and branch annotations in output")

    inout_grp.add_argument("--basename", action="store", type=Path, dest="outbase", metavar="NAME",
        help="base name for output files [default: derived from first input file]")

    inout_grp.add_argument("-n", action="store_true", dest="nowarn",
        help="no warning: overwrite existing output files without prompting")

    inout_grp.add_argument("-v", action="store_true", dest="verbose",
        help="verbose: show full Python traceback on error")

    inout_grp.add_argument("-q", action="store_true", dest="quiet",
        help="quiet: suppress progress output; also implies -n")

    inout_grp.add_argument("infilelist", nargs="+", metavar="FILE", type=Path,
        help="input phylogenetic tree file(s)")

    ####################################################################################

    sumtype_grp = parser.add_argument_group("SUMMARY TREE TYPE (choose one)")
    sumtype_excl = sumtype_grp.add_mutually_exclusive_group(required=True)

    sumtype_excl.add_argument("--con", dest="treetype", action="store_const", const="con",
        help="majority rule consensus tree")

    sumtype_excl.add_argument("--all", dest="treetype", action="store_const", const="all",
        help="majority rule consensus tree with all compatible bipartitions")

    sumtype_excl.add_argument("--mcc", dest="treetype", action="store_const", const="mcc",
        help= "maximum clade credibility (MCC) tree. "
              "The MCC tree is determined by inspecting tree samples and selecting the "
              "tree that has the highest product of clade frequencies (= highest sum of "
              "log of clade frequencies). The MCC tree is therefore a tree that has been "
              "observed in the pool of tree samples, differing from the consensus tree "
              "which typically does not match any individual sample. "
              "Meaningful mainly for rooted input trees, for "
              "example from clock-based analyses. By default, the rooting of the chosen "
              "sample tree is retained.")

    sumtype_excl.add_argument("--mbc", dest="treetype", action="store_const", const="mbc",
        help="maximum bipartition credibility (MBC) tree. "
             "Similar to the MCC tree, but uses bipartitions rather than clades and "
             "therefore ignores rooting. "
             "(two input trees can have the same set of bipartitions, but be rooted "
             "in different locations).")

    sumtype_excl.add_argument("--hip", dest="treetype", action="store_const", const="hip",
        help="HIPSTR summary tree (Highest Independent Posterior SubTree; "
             "Baele et al., Bioinformatics, 2025, 41(10)). "
             "Builds a fully resolved summary tree by choosing, at each internal node, "
             "the child-clade pair with the highest combined posterior support. "
             "Like a consensus tree, a HIPSTR tree has not necessarily been observed "
             "among the input trees.")

    sumtype_excl.add_argument("--mrhip", dest="treetype", action="store_const", const="mrhip",
        help="MrHIPSTR (majority-rule HIPSTR) summary tree. Like --hip, but includes only clades with "
                 "at least 50%% support.")

    ####################################################################################

    blen_grp = parser.add_argument_group(title= "BRANCH-LENGTH ESTIMATION (choose one)")
    blen_excl = blen_grp.add_mutually_exclusive_group(required=True)
    blen_excl.add_argument("--noblen", action="store_true",
        help="do not estimate branch lengths; compute topology and branch- or clade-support only")

    blen_excl.add_argument("--biplen", action="store_true",
        help="set branch lengths to the mean length of the corresponding leaf "
             "bipartition across input trees: "
             "each branch in tree corresponds to a bipartition of the leaves "
             "into two groups. Branch lenghts are set to the mean of the length of the "
             "corresponding bipartition across all input trees.")

    blen_excl.add_argument("--cladedepth", action="store_true",
        help="set node depths based on the depth of each observed clade, then derive "
                 "branch lengths from those depths. Intended for rooted, clock-like trees. "
                 "Requires all clades in the summary tree to have been observed in the input "
                 "trees and may fail for some rootings. "
                 "Mean is computed across trees where the specific, monophyletic clade "
                 "is present, and may therefore be based on very few (down to one) values. "
                 "May produce negative branch lengths.")

    # This is deprecated now
    blen_excl.add_argument("--meandepth", dest="cladedepth", action="store_true", help=argparse.SUPPRESS)

    blen_excl.add_argument("--cadepth", action="store_true",
        help="'common ancestor depth'; equivalent to TreeAnnotator --height ca."
             "Set node depths by mean MRCA depth across all trees, then derive branch "
             "lengths from those depths. Intended for rooted, clock-like trees.")

    ####################################################################################

    root_grp = parser.add_argument_group("ROOTING")

    root_excl = root_grp.add_mutually_exclusive_group()

    root_excl.add_argument("--rootmid", action="store_true",
        help="midpoint-root the summary tree")

    root_excl.add_argument("--rootminvar", action="store_true",
        help="root the summary tree by minimum-variance rooting")

    root_excl.add_argument("--rootog", dest="outgroup", metavar="TAX[,TAX,...]", type=str, default=None,
        help="root the summary tree on the specified outgroup taxon or taxa")

    root_excl.add_argument('--rootogfile', dest="ogfile", action="store", metavar="FILE", default=None,
        help="root the summary tree on outgroup taxa listed in FILE, one per line")

    root_grp.add_argument("--rootcred", action="store_true",
        help="compute root credibility for branches in the summary tree. "
             "With an outgroup, track the branch to which the outgroup attaches in each "
             "input tree. Otherwise assume rooted input trees and track observed root "
             "locations directly.")

    ####################################################################################

    bayes_grp = parser.add_argument_group("BAYESIAN OPTIONS")

    bayes_grp.add_argument("-b", dest="burninfrac", metavar="FRAC[,FRAC,...]", type=str, default="0",
        help="burn-in fraction(s) to discard [0,1]. Supply one value for all files, "
              "or one comma-separated value per input file [default: %(default)s]")

    bayes_grp.add_argument("--ci", metavar="PROB[,PROB,...]", type=str,
        help="compute one or more central credible intervals for branch lengths or node depths; "
             "for example, --ci 0.8,0.95. Credible interval endpoints are approximate quantiles "
             "(mergeable log-bucket histogram); use --cik to control precision.")

    bayes_grp.add_argument("--cik", dest="ci_k", type=int, default=7, metavar="K",
        help="set precision for --ci quantile estimation using a mergeable log-bucket histogram. "
             "Uses 2^K sub-bins per exponent bucket: higher K gives finer resolution but uses more "
             "memory/CPU. Worst-case relative bucket error is about 2^-(K+1): "
             "K=7 ~0.39%%, K=8 ~0.20%%, K=9 ~0.10%%. [default: %(default)s]")

    bayes_grp.add_argument("-t", type=float, dest="treeprobs", metavar="PROB",
        help="compute tree probabilities and report the PROB credible set [0,1]")

    bayes_grp.add_argument("-s", action="store_true", dest="std",
        help="compute average standard deviation of split frequencies (ASDSF) across "
             "input tree files")

    bayes_grp.add_argument("-f", type=float, dest="minfreq", metavar="NUM", default=0.1,
        help="minimum split frequency included in ASDSF computation [default: %(default)s]")

    ####################################################################################

    perf_grp = parser.add_argument_group("PERFORMANCE")

    perf_grp.add_argument("--cpus", type=int, default=0, metavar="N",
        help="number of CPUs to use for parallel processing. "
            "[default: 0 = automatic; 1 = run without parallel processing]")

    perf_grp.add_argument("--chunksize", type=int, default=250, metavar="N",
        help="number of trees per work chunk [default: %(default)s]. "
             "Larger values reduce overhead but use more memory and may reduce load balancing.")

    ####################################################################################

    other_grp = parser.add_argument_group("OTHER OPTIONS")

    other_grp.add_argument("--nolabel", action="store_true", dest="nolabel",
        help=argparse.SUPPRESS # do not print branch labels (=clade probabilities) on summary tree
        )

    ####################################################################################

    return parser

####################################################################################

def read_outgroup(ogfile):
    infile = open(ogfile, "r")
    outgroup = []
    for line in infile:
        leaf = line.strip()
        outgroup.append(leaf)
    infile.close()
    return outgroup

####################################################################################

class OutputManager:
    def __init__(self, args):
        self.quiet = args.quiet
        self.verbose = args.verbose

    def info(self, message="", padding=3, end="\n"):
        if not self.quiet:
            print(f"{padding * ' '}{message}", end=end)

    def force(self, message="", padding=3, end="\n"):
        if not self.quiet:
            print(f"{padding * ' '}{message}", end=end)
            sys.stdout.flush()

    def warning(self, message):
        if not self.quiet and self.verbose:
            print(f"{self.padding * ' '}[WARNING] {message}")

####################################################################################

def setup_output_directory(outbase):
    outbase.parent.mkdir(parents=True, exist_ok=True)

####################################################################################

def count_trees(args, output):

    count_burnin_filename_list = []
    n_trees_analyzed = 0
    for i,filename in enumerate(args.infilelist):
        treelist = []
        n_tot = fast_treecount(filename, args)
        burnin = int(args.burninfrac[i] * n_tot)
        count_burnin_filename_list.append((n_tot, burnin, filename))
        n_trees_analyzed += (n_tot - burnin)

    return (n_trees_analyzed, count_burnin_filename_list)

####################################################################################

def fast_treecount(filename, args):
    """Heuristic: count patterns ([;=\n] etc) to infer number of trees"""

    # Empirically: if ); is in file, then this == number of trees
    n_terminators = count_bytestring(filename, b");")
    if n_terminators > 0:
        return n_terminators

    # Count semicolon: if == 1 (possibly after nexus correction): 1 tree
    n_semicolons = count_bytestring(filename, b";")
    if n_semicolons == 1:
        return 1
    if args.informat == "nexus":
        n_other_semicolon_patterns  = count_bytestring(filename, b"egin taxa;")
        n_other_semicolon_patterns += count_bytestring(filename, b"egin trees;")
        n_other_semicolon_patterns += count_bytestring(filename, b"nd;")
        n_other_semicolon_patterns += count_bytestring(filename, b"imensions ntax")
        n_other_semicolon_patterns += count_bytestring(filename, b"axlabels")
        n_other_semicolon_patterns += count_bytestring(filename, b"ranslate")
        n_semicolons -= n_other_semicolon_patterns
    if n_semicolons == 1:
        return 1

    # If we got this far, and filetype is newick then bail out and use precise counting
    if args.informat == "newick":
        return count_trees_by_parsing(filename, args)

    # Final attempt to infer ntrees for nexus files:
    # count "= (",  and "tree "
    # Add the values that are not 0 to list, choose minimum as count
    # Should be robust to most variations, but should check at end of sumt...
    n_eqparen = count_bytestring(filename, b"= (")
    n_treestr = count_bytestring(filename, b"tree ")
    countlist = [n_semicolons, n_eqparen, n_treestr]
    notzero = [val for val in countlist if val>0]
    return min(notzero)

####################################################################################

def count_bytestring(filename, bytestring):
    """Fast counting of specific pattern. Bytestring argument must be given
    with b modifier (e.g., b');')"""

    # Modified from: https://stackoverflow.com/a/27517681/7836730
    with open(filename, 'rb') as f:
        bufsize = 1024*1024
        bufgen = takewhile(lambda x: x, (f.raw.read(bufsize) for _ in repeat(None)))

        prev_buf = b""
        count = 0

        for buf in bufgen:
            count += buf.count(bytestring)

            # For multi-byte patterns, consider overlaps between buffers
            if len(bytestring) > 1 and len(prev_buf) > 0:
                merged = prev_buf[-len(bytestring)+1:] + buf[:len(bytestring)-1]
                count += merged.count(bytestring)

            prev_buf = buf

    return count

####################################################################################

def count_trees_by_parsing(filename, args):
    # Open treefile. Discard (i.e., silently pass by) the requested number of trees
    if args.informat == "nexus":
        treefile = pt.Nexustreefile(filename)
    else:
        treefile = pt.Newicktreefile(filename)
    treecount = 0
    for tree in treefile:
        treecount += 1
    return treecount

####################################################################################

def process_trees(count_burnin_filename_list, args, output):

    treesummarylist = []
    for i, (count, burnin, filename) in enumerate(count_burnin_filename_list):
        output.force()
        output.force(f"Analyzing file: {filename}")

        # Open treefile. Discard (i.e., silently pass by) the requested number of trees
        if args.informat == "nexus":
            treefile = pt.Nexustreefile(filename)
        else:
            treefile = pt.Newicktreefile(filename)
        for j in range(burnin):
            treefile.readtree(returntree=False)
        output.info(f"Discarded {burnin:,} of {count:,} trees (burnin fraction={args.burninfrac[i]:.2f})")

        # Instantiate Treesummary.
        treesummary = pt.TreeSummary(
            trackbips=args.trackbips,
            trackclades=args.trackclades,
            trackroot=args.trackroot,
            trackblen=args.trackblen,
            trackrootblen=args.trackrootblen,
            trackdepth=args.trackdepth,
            tracktopo=args.tracktopo,
            track_subcladepairs=args.track_subcladepairs,
            store_trees=args.treeprobs,
            trackci=args.trackci,
            ci_probs=args.ci_probs,
            quantile_k=args.ci_k
        )

        # Initialize the progress bar
        progress = ProgressBar(ntot=count-burnin, output=output, quiet=args.quiet, text="Processing trees:")

        # Read post-burnin trees from file, add to treesummary, print progress bar
        for j in range(burnin, count):
            tree = treefile.readtree(returntree=True)
            treesummary.add_tree(tree)
            progress.update()
            del tree

        # Ensure the progress bar completes at the end
        progress.complete()

        treesummarylist.append(treesummary)
        output.info()

    return treesummarylist

####################################################################################

def process_trees_concurrent(count_burnin_filename_list, args, output, n_trees_analyzed):
    """Process trees using multiple processors."""

    output.info()
    filetable = file_overview_table(count_burnin_filename_list, output)
    output.info(filetable, padding=0)
    output.info()

    progress = ProgressBar(ntot=n_trees_analyzed, output=output, quiet=args.quiet,
                           text="Processing trees:")

    ncpus = args.cpus
    max_pending = 2 * ncpus
    max_chunk_size = args.chunksize

    # Global summaries (in parent only)
    treesummarylist = []
    for _ in range(len(count_burnin_filename_list)):
        ts_global = pt.TreeSummary(
            trackbips=args.trackbips,
            trackclades=args.trackclades,
            trackroot=args.trackroot,
            trackblen=args.trackblen,
            trackrootblen=args.trackrootblen,
            trackdepth=args.trackdepth,
            tracktopo=args.tracktopo,
            track_subcladepairs=args.track_subcladepairs,
            store_trees=args.treeprobs,
            trackci=args.trackci,
            ci_probs=args.ci_probs,
            quantile_k=args.ci_k
        )
        treesummarylist.append(ts_global)

    chunk_iterator = chunked_tree_strings_from_files(count_burnin_filename_list, max_chunk_size)
    worker_pids = set()

    with ProcessPoolExecutor(max_workers=ncpus) as ex:
        # Prime the pipeline
        pending = {
            ex.submit(worker_process_chunk, chunk, file_idx, parser_obj, args)
            for chunk, file_idx, parser_obj in islice(chunk_iterator, max_pending)
        }

        # Harvest + refill
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for fut in done:
                tree_summary, file_idx, pid = fut.result()
                worker_pids.add(pid)
                treesummarylist[file_idx].update(tree_summary)
                progress.update(len(tree_summary))

            # top up to max_pending (if input remains)
            n_to_submit = max_pending - len(pending)
            for chunk, file_idx, parser_obj in islice(chunk_iterator, n_to_submit):
                pending.add(ex.submit(worker_process_chunk, chunk, file_idx, parser_obj, args))

    output.info()
    return treesummarylist, worker_pids

##########################################################################################

def file_overview_table(count_burnin_filename_list, output, show_path=False, padding=3):
    """Prints table with columns for filename, treecount, burnin-fraction, burnin, kept"""

    filename_header = "Filename"
    treecount_header = "Treecount"
    burnfrac_header  = "BurninFrac"
    burnin_header    = "Burnin"
    kept_header  = "Kept"

    rows = []
    for count, burnin, filename in count_burnin_filename_list:
        fn = str(filename) if show_path else Path(filename).name
        kept = count - burnin
        burnfrac = (burnin / count) if count else 0.0  # avoid ZeroDivisionError
        burnfrac_str = f"{burnfrac:.2f}"
        rows.append((fn, count, burnfrac_str, burnin, kept))

    # Widths (at least as wide as the headers)
    fn_w = max([len(filename_header)] + [len(r[0]) for r in rows])
    tc_w = max([len(treecount_header)] + [len(str(r[1])) for r in rows])
    bf_w = max([len(burnfrac_header)]  + [len(r[2]) for r in rows])
    bi_w = max([len(burnin_header)]    + [len(str(r[3])) for r in rows])
    an_w = max([len(kept_header)]  + [len(str(r[4])) for r in rows])

    header = (
        f"{padding * ' '}{filename_header:<{fn_w}}  "
        f"{treecount_header:>{tc_w}}  "
        f"{burnfrac_header:>{bf_w}}  "
        f"{burnin_header:>{bi_w}}  "
        f"{kept_header:>{an_w}}"
    )
    divider = (
        f"{padding * ' '}{'-' * fn_w}  "
        f"{'-' * tc_w}  "
        f"{'-' * bf_w}  "
        f"{'-' * bi_w}  "
        f"{'-' * an_w}"
    )

    lines = [header, divider]
    for fn, count, burnfrac_str, burnin, kept in rows:
        lines.append(
            f"{padding * ' '}{fn:<{fn_w}}  "
            f"{count:>{tc_w}d}  "
            f"{burnfrac_str:>{bf_w}}  "
            f"{burnin:>{bi_w}d}  "
            f"{kept:>{an_w}d}"
        )

    return "\n".join(lines)

##########################################################################################

def chunked_tree_strings_from_files(count_burnin_filename_list, max_chunk_size):
    """yields (chunk_of_tree_strings, file_idx, parser_obj)"""

    for file_idx, (count, burnin, filename) in enumerate(count_burnin_filename_list):
        with pt.Treefile(filename, strip_comments=False) as tf:
            parser_obj = tf.parser_obj

            for _ in range(burnin):
                tf.readtree(returntree=False)

            chunk = []
            for _ in range(burnin, count):
                treestr = tf.readtree(returntree=False)
                chunk.append(treestr)
                if len(chunk) == max_chunk_size:
                    yield (chunk, file_idx, parser_obj)
                    chunk = []
            if chunk:
                yield (chunk, file_idx, parser_obj)

##########################################################################################

def worker_process_chunk(chunk, file_idx, parser_obj, args):
    """Worker: parse chunk and return (TreeSummary, file_idx, pid)."""

    ts = pt.TreeSummary(
        trackbips=args.trackbips,
        trackclades=args.trackclades,
        trackroot=args.trackroot,
        trackblen=args.trackblen,
        trackrootblen=args.trackrootblen,
        trackdepth=args.trackdepth,
        tracktopo=args.tracktopo,
        track_subcladepairs=args.track_subcladepairs,
        store_trees=args.treeprobs,
        trackci=args.trackci,
        ci_probs=args.ci_probs,
        quantile_k=args.ci_k
    )

    for treestr in chunk:
        treestr = pt.remove_comments(treestr)
        tree = pt.Tree._from_string_private(parser_obj, treestr)
        ts.add_tree(tree)

    return (ts, file_idx, os.getpid())

##########################################################################################

def _ca_worker_init(plan, trackci, quantile_k):
    global _CA_PLAN, _CA_TRACKCI, _CA_QUANTILE_K
    _CA_PLAN = plan
    _CA_TRACKCI = bool(trackci)
    _CA_QUANTILE_K = int(quantile_k)

####################################################################################

def worker_process_ca_chunk(chunk, parser_obj):
    """
    Worker: parse chunk -> update local CADepthEstimator -> return it.
    """
    est = pt.CADepthEstimator(_CA_PLAN, trackci=_CA_TRACKCI, quantile_k=_CA_QUANTILE_K)
    for treestr in chunk:
        treestr = pt.remove_comments(treestr)
        tree = pt.Tree._from_string_private(parser_obj, treestr)
        est.add_tree(tree)
    return (est, len(chunk), os.getpid())

####################################################################################

def set_ca_depths_concurrent(sumtree, count_burnin_filename_list, args, output, n_trees_analyzed):
    """
    Second pass over the tree samples, parallelized.
    Returns sumtree with CA depths written into nodes.
    """
    # Build plan once in parent
    plan = pt.CADepthEstimator.build_plan(
        sumtree,
        trackci=args.trackci and bool(args.ci_probs),
        ci_probs=args.ci_probs,
    )

    global_est = pt.CADepthEstimator(
        plan,
        trackci=args.trackci and bool(args.ci_probs),
        quantile_k=args.ci_k,
    )

    progress = ProgressBar(ntot=n_trees_analyzed, output=output, quiet=args.quiet,
                           text="Finding common ancestor node-depths:")

    ncpus = args.cpus
    max_pending = 2 * ncpus
    max_chunk_size = args.chunksize

    chunk_iterator = chunked_tree_strings_from_files(count_burnin_filename_list, max_chunk_size)
    worker_pids = set()

    with ProcessPoolExecutor(
        max_workers=ncpus,
        initializer=_ca_worker_init,
        initargs=(plan, args.trackci and bool(args.ci_probs), args.ci_k),
    ) as ex:

        # Prime
        pending = {
            ex.submit(worker_process_ca_chunk, chunk, parser_obj)
            for chunk, file_idx, parser_obj in islice(chunk_iterator, max_pending)
        }

        # Harvest + refill
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for fut in done:
                est_part, ntrees, pid = fut.result()
                worker_pids.add(pid)
                global_est.merge(est_part)
                progress.update(ntrees)

            # top up to max_pending (if input remains)
            n_to_submit = max_pending - len(pending)
            for chunk, file_idx, parser_obj in islice(chunk_iterator, n_to_submit):
                pending.add(ex.submit(worker_process_ca_chunk, chunk, parser_obj))

    output.info()
    global_est.write_into(sumtree)
    return sumtree, worker_pids

####################################################################################

class ProgressBar:
    def __init__(self, ntot, output, quiet=False, text="Processing items:"):
        self.output = output
        self.quiet = quiet
        self.progscale = "0      10      20      30      40      50      60      70      80      90     100"
        self.progticks = "v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v"
        self.ndots = len(self.progticks)
        self.ntot = ntot
        self.finished = 0
        self.items_per_dot = self.ntot / self.ndots
        self.n_dotsprinted = 0
        self.text = text

        # Print progress bar header
        if not self.quiet:
            self.output.force()
            self.output.force(self.text)
            self.output.force(f"{self.progscale}")
            self.output.force(f"{self.progticks}")
            self.output.force("", end="")

    def update(self, nfinished=1):
        """ Update the progress bar based on the number of processed trees. """
        if not self.quiet:
            self.finished += nfinished
            n_dots_expected = math.floor(self.finished / self.items_per_dot)
            if self.n_dotsprinted < n_dots_expected:
                n_missing = n_dots_expected - self.n_dotsprinted
                self.output.force("*" * n_missing, padding=0, end="")
                self.n_dotsprinted += n_missing

    def complete(self):
        """ Ensure all dots are printed at the end if they haven't been already. """
        if not self.quiet:
            if self.n_dotsprinted < self.ndots:
                n_missing = self.ndots - self.n_dotsprinted
                self.output.force("*" * n_missing)

####################################################################################

def compute_converge_stats(treesummarylist, args):
    """Compute average bipartition/clade frequency standard deviation between treesummaries"""

    # NOTES ON COMPUTATION:
    #   (1) all bipartitions/clades are included in computation, regardless of whether they are present
    #       in only one of the treesummaries (their freq is set to zero in those treesummaries
    #       where they are not present)
    #   (2) external branches (leading to leafs) are not included in computation since they
    #       will allways be present in all trees (=> p=1, std=0)
    #   (3) N-1 is used in the denominator when computing std. This has to do with sample vs.
    #       population estimates of std, and is in accordance with what MrBayes does. One could
    #       argue that N should be used instead of N-1 (to get the Max Likelihood estimate of std).

    if args.treetype in ("mcc", "hip", "mrhip"):
        ave_std = compute_converge_clades(treesummarylist, args)
    else:
        ave_std = compute_converge_biparts(treesummarylist, args)
    return ave_std

##########################################################################################

def compute_converge_clades(treesummarylist, args):
    """Compute average clade frequency standard deviation between treesummaries"""
    sum_std = 0
    N = float(len(treesummarylist))

    # Find combined set of clades (excluding leaves)
    # Only clades that have freq >= minfreq are kept
    cladeset = set()
    for treesummary in treesummarylist:
        for clade,node in treesummary.cladesummary.items():
            leafset = clade.get_clade()
            if len(leafset)>1 and node.posterior >= args.minfreq:
                cladeset.add(clade)

    # For each clade: compute std of freq of this clade across all treesummaries
    for clade in cladeset:
        freqlist = []
        for treesummary in treesummarylist:
            # If current clade not in current treesummary: set freq=0.0
            if clade in treesummary.cladesummary:
                freqlist.append(treesummary.cladesummary[clade].posterior)
            else:
                freqlist.append(0.0)
        sum_std += statistics.stdev(freqlist)

    ave_std = sum_std / len(cladeset)
    return ave_std

##########################################################################################

def compute_converge_biparts(treesummarylist, args):
    """Compute average bipartition frequency standard deviation between treesummaries"""
    sum_std = 0
    N = float(len(treesummarylist))

    # Find combined set of bipartitions (excluding external branches)
    # Only biparts that have freq >= minfreq are kept
    bipset = set()
    for treesummary in treesummarylist:
        for bipart,branch in treesummary.bipartsummary.items():
            (bip1, bip2) = bipart
            if len(bip1)>1 and len(bip2)>1 and branch.posterior >= args.minfreq:
                bipset.add(bipart)

    # For each internal bipart: compute std of freq of this bipart across all treesummaries
    for bipart in bipset:
        freqlist = []
        for treesummary in treesummarylist:
            # If current bipartition not in current treesummary: set freq=0.0
            if bipart in treesummary.bipartsummary:
                freqlist.append(treesummary.bipartsummary[bipart].posterior)
            else:
                freqlist.append(0.0)
        sum_std += statistics.stdev(freqlist)

    ave_std = sum_std / len(bipset)
    return ave_std

##########################################################################################

def  merge_treesummaries(treesummarylist):
    """Combine multiple TreeSummary objects into one"""

    treesummary = treesummarylist[0]
    for treesummary2 in treesummarylist[1:]:
        treesummary.update(treesummary2)
        del treesummary2
    return treesummary

##########################################################################################

def compute_sumtree(treesummary, args, count_burnin_filename_list, output, n_trees_analyzed, worker_pids=None):
    """Controls computation of summary tree, setting of node depths and branch lengths,
       and annotation of summary tree with relevant attributes.
       """

    # Rooting: only do something if user asked for it; otherwise keep "input"
    if args.outgroup:
        rooting = "og"
        og = args.outgroup
    elif args.rootmid:
        rooting = "mid"
        og = None
    elif args.rootminvar:
        rooting = "minvar"
        og = None
    else:
        rooting = None
        og = None

    output.info()
    output.force("Computing summary tree")

    if args.cadepth and args.cpus > 1:
        # 1) Build topology + root, but do NOT compute depths/blen
        sumtree = pt.build_sumtree(treesummary, treetype=args.treetype, rooting=rooting,
                                   blen="none", og=og, count_burnin_filename_list=None)

        # 2) Parallel CA depths (second pass)
        sumtree, ca_pids = set_ca_depths_concurrent(
            sumtree, count_burnin_filename_list, args, output, n_trees_analyzed)
        if worker_pids is not None:
            worker_pids |= ca_pids

        # 3) Branch lengths from depths
        sumtree.set_blens_from_depths()

        # 4) Annotate (support/rootcred etc)
        tpp = pt.TreePostProcessor(treesummary)
        sumtree = tpp.annotate_sumtree(sumtree)

    else:
        # Serial (non-parallel) processing
        blen = "cadepth" if args.cadepth else ("cladedepth" if args.cladedepth else ("biplen" if args.biplen else "none"))
        sumtree = pt.build_sumtree(treesummary, treetype=args.treetype, rooting=rooting,
                                   blen=blen, og=og, count_burnin_filename_list=count_burnin_filename_list)

    return sumtree

##########################################################################################

def print_sumtree(sumtree, treesummary, args, output):

    confilename = args.outbase.parent / (args.outbase.name + f".{args.treetype}")

    pt.configure_sumtree_printing(
        sumtree,
        treetype=args.treetype,
        blen=("cadepth" if args.cadepth else ("cladedepth" if args.cladedepth else ("biplen" if args.biplen else "none"))),
        trackci=args.trackci,
        ci_labels=(treesummary.ci_labels if args.trackci else None),
        precision=7,
        print_meta=not args.nometa,
        printlabels=not args.nolabel,
        printdist=not args.noblen,
    )

    with open_file_with_warning(confilename, args.nowarn, output) as confile:
        tree_str = sumtree.nexus() if args.outformat == "nexus" else sumtree.newick()
        confile.write(tree_str)
        confile.write("\n")

    output.info("")
    if args.treetype == "mbc":
        output.info(f"Maximum bipartition credibility tree written to {confilename}")
    elif args.treetype == "mcc":
        output.info(f"Maximum clade credibility tree written to {confilename}")
    elif args.treetype == "hip":
        output.info(f"HIPSTR tree written to {confilename}")
    elif args.treetype == "mrhip":
        output.info(f"Majority rule HIPSTR (mrHIPSTR) tree written to {confilename}")
    else:
        output.info(f"Consensus tree written to {confilename}")

##########################################################################################

def open_file_with_warning(filename, nowarn, output):
    if nowarn:
        return open(filename, "w")
    elif filename.is_file():
        overwrite = input(f"\n   File {filename} already exists.\n   Overwrite (y/n): ")
        if overwrite == "y":
            output.info(f"Overwriting file {filename}")
            output.info()
            return open(filename, "w")  # Overwrite
        else:
            output.info(f"Appending to file {filename}")
            output.info()
            return open(filename, "a")  # Append
    else:
        return open(filename, "w")

##########################################################################################

def compute_trprobs(treesummary, args):
    """Returns sorted list of [freq, tree] lists (highest freq first)"""

    # Python note: root trees in trprobs?
    if treesummary.trackclades:
        toposummary = treesummary.cladetoposummary
    elif treesummary.trackbips:
        toposummary = treesummary.biptoposummary
    trproblist = []
    for topology, topostruct in toposummary.items():
        trproblist.append((topostruct.posterior, topostruct.tree))

    # Sort report according to frequency (higher values first) and return
    trproblist = sorted(trproblist, key=itemgetter(0), reverse=True)

    return trproblist

##########################################################################################

def print_trprobs(treesummary, trproblist, args, output):
    topofilename = args.outbase.parent / (args.outbase.name + ".trprobs")
    with open_file_with_warning(topofilename, args.nowarn, output) as topofile:
        topofile.write("#NEXUS\n\n")
        if args.treeprobs < 1:
            topofile.write(f"[This file contains the {round(args.treeprobs*100)}% most probable trees found during the\n")
            topofile.write(f"MCMC search, sorted by posterior probability (the {round(args.treeprobs*100)}% HPD interval).\n")
        else:
            topofile.write("[This file contains all trees that were found during the MCMC\n")
            topofile.write("search, sorted by posterior probability. \n")
        topofile.write("Lower case 'p' indicates the posterior probability of a tree.\n")
        topofile.write("Upper case 'P' indicates the cumulative posterior probability.]\n\n")
        topofile.write("begin trees;\n")
        topofile.write(treesummary.translateblock)
        n = 1
        cum = 0.0
        for (freq, tree) in trproblist:
            cum += freq
            treestring = tree.newick(printdist=False, printlabels=False, transdict=treesummary.transdict)
            topofile.write(f"    tree tree_{n} [p = {freq:.6f}] [P = {cum:.6f}] = {treestring}\n")
            n += 1
            if cum > args.treeprobs:
                break
        topofile.write("end;\n")

    return f"Tree probabilities written to {topofilename}"

##########################################################################################

def fmt_mem(n_bytes):
    if n_bytes > 1E9:
        return f"{n_bytes / 1024**3:.1f} GB"
    else:
        return f"{n_bytes / 1024**2:.1f} MB"

def memory_summary_line(peak):
    total = psutil.virtual_memory().total
    pct = (peak.peak_total_rss_bytes / total * 100.0) if total else 0.0
    return (
        f"Peak memory used: "
        f"{fmt_mem(peak.peak_total_rss_bytes)} / {fmt_mem(total)} ({pct:.1f}%), "
        f"free at peak: {fmt_mem(peak.peak_available_bytes_at_peak)}"
    )

##########################################################################################

def print_result_summary(sumtree, treesummary, start, n_trees_analyzed,
                         ave_std, args, output, mem_mon, worker_pids=None):

    sumvar_tuple = compute_summary_variables(sumtree, treesummary, start, args)
    (n_leaves, grouptype, space, n_uniq_groupings, theo_max_groups, theo_max_biparts, n_topo_seen,
     treetype, n_internal_biparts, rootdegree, h, m, s) = sumvar_tuple

    # Information about bipartitions, clades and topologies
    output.info()
    output.info(f"Number of leaves on input trees: {n_leaves:>7,d}")
    if args.treeprobs:
        output.info(f"Different topologies seen: {n_topo_seen:>13,d}")
        output.info(f"Different {grouptype}s seen:{space}{n_uniq_groupings:>11,d} (theoretical maximum: {theo_max_groups * n_topo_seen:,d})")
    else:
        output.info(f"Different {grouptype}s seen:{space}{n_uniq_groupings:>11,d} (theoretical maximum: {theo_max_groups * n_trees_analyzed:,d})")
    tmpstr = f"Bipartitions in {treetype} tree:"
    output.info(f"{tmpstr:<34}{n_internal_biparts:>6,d} (theoretical maximum: {theo_max_biparts:,d})")

    if n_internal_biparts < theo_max_biparts:
        output.info("(tree contains polytomies)", padding=44)
    else:
        output.info("(tree is fully resolved - no polytomies)", padding=44)

    # Information about branch lengths
    if args.cladedepth:
        output.info()
        output.info(f"Branch lengths set based on mean node depths in input trees")
    if args.cadepth:
        output.info()
        output.info(f"Branch lengths set based on common ancestor depths in input trees")
    elif args.biplen:
        output.info()
        output.info(f"Branch lengths set based on mean branch lengths for corresponding bipartitions")
    elif args.noblen:
        output.info()
        output.info(f"Branch lengths have not been tracked")

    # Information about rooting
    output.info()
    if not args.actively_rooted:
        if args.treetype == "mcc":
            output.info(f"{treetype} tree rooted at original root of tree sample having highest clade credibility")
        elif args.treetype in ("hip", "mrhip"):
            output.info(f"{treetype} tree rooted at most frequently observed root bipartition")
        else:
            output.info(f"{treetype} tree has not been explicitly rooted")
            output.info(f"Tree has been rooted at random internal node; root is at {rootdegree}")
    else:
        if args.outgroup:
            output.info(f"{treetype} tree has been rooted based on outgroup")
        elif args.rootmid:
            output.info(f"{treetype} tree has been midpoint-rooted")
        elif args.rootminvar:
            output.info(f"{treetype} tree has been rooted using minimum variance-rooting")
        else:
            raise TreeError("rooting error") # Python note: Should never go here - remove when code tested

    if args.rootcred:
        if sumtree.rootcred:
            output.info(f"Root credibility (frequency of root bipartition in input trees):       {sumtree.rootcred * 100:.1f}%")
        output.info(f"Cumulated root credibility (sum of rootcred for all branches in tree): {sumtree.cumulated_rootcred * 100:.1f}%")


    # Information about log credibility
    if args.treetype == "mbc" or (args.treetype == "mcc" and not args.actively_rooted):
        output.info()
        output.info(f"Highest log {grouptype} credibility:  {sumtree.logcred:.6g}")
    else:
        output.info()
        output.info(f"Log {grouptype} credibility:  {sumtree.logcred:.6g}")

    if args.std:
        output.info(f"Average standard deviation of split frequencies: {ave_std:.6f}")

    output.info()
    output.info(f"Done. {n_trees_analyzed:,d} trees analyzed.\n   Time spent: {h:d}:{m:02d}:{s:02d} (h:m:s)")

    # Processor count reporting
    if worker_pids:
        procs_used = len(worker_pids)
        output.info(cpu_process_summary_line(procs_used))
    else:
        # Serial path: 1 process (this one)
        output.info(cpu_process_summary_line(1))

    peakmem = mem_mon.stop()
    output.info(memory_summary_line(peakmem))

##########################################################################################

def compute_summary_variables(sumtree, treesummary, start, args):

    if args.treeprobs:
        if args.trackclades:
            n_topo_seen = len(treesummary.cladetoposummary)
        elif args.trackbips:
            n_topo_seen = len(treesummary.biptoposummary)
    else:
        n_topo_seen = None

    nrootkids = len(sumtree.children(sumtree.root))
    if nrootkids == 2:
        rootdegree = "bifurcation"
    elif nrootkids == 3:
        rootdegree = "trifurcation"
    else:
        rootdegree = "multifurcation"

    n_leaves = len(treesummary.leaves)

    if args.treetype in ("mcc", "hip", "mrhip"):
        n_uniq_groupings = len(treesummary.cladesummary) - n_leaves
        theo_max_groups = n_leaves - 1
    else:
        n_uniq_groupings = len(treesummary.bipartsummary) - n_leaves
        theo_max_groups = n_leaves - 3

    theo_max_biparts = n_leaves - 3
    n_internal_biparts = sumtree.n_bipartitions()

    if args.treetype == "mcc":
        treetype, grouptype, space = "MCC", "clade", " " * 7
    elif args.treetype == "mbc":
        treetype, grouptype, space = "MBC", "bipartition", " " * 1
    elif args.treetype == "hip":
        treetype, grouptype, space = "HIPSTR", "clade", " " * 7
    elif args.treetype == "mrhip":
        treetype, grouptype, space = "mrHIPSTR", "clade", " " * 7
    else:
        treetype, grouptype, space = "Consensus", "bipartition", " " * 1

    time_spent = time.time() - start
    h = int(time_spent/3600)
    m = int((time_spent % 3600)/60)
    s = int(time_spent % 60)

    return   (n_leaves, grouptype, space, n_uniq_groupings, theo_max_groups, theo_max_biparts, n_topo_seen,
              treetype, n_internal_biparts, rootdegree, h, m, s)

##########################################################################################

def handle_error(error, verbose):
    print("\n\nExecution failed:")
    if verbose:
        import traceback
        traceback.print_exc(file=sys.stdout)
    else:
        print(error)
    sys.exit()

##########################################################################################

if __name__ == "__main__":
    main()
