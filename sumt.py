# By Anders Gorm Pedersen, agpe@dtu.dk

import phylotreelib as pt
import argparse, os, sys, time, math, copy, psutil, statistics, configparser
from itertools import (takewhile,repeat)
from operator import itemgetter
from pathlib import Path
import gc
gc.disable()        # Faster. Assume no cyclic references will ever be created

####################################################################################

def main(commandlist=None):
    start=time.time()
    pid = psutil.Process(os.getpid())
    args = parse_commandline(commandlist)
    output = OutputManager(args)

    try:
        setup_output_directory(args.outbase)        
        wt_file_list = parse_infilelist(args)
        n_trees_analyzed, wt_count_burnin_filename_list = count_trees(wt_file_list, args, output)
        memory1 = track_memory_usage(pid)
        
        treesummarylist = process_trees(wt_count_burnin_filename_list, args, output)
        if args.std:
            ave_std = compute_converge_stats(treesummarylist, args)
        else:
            ave_std = None
        treesummary = merge_treesummaries(treesummarylist)
        
        sumtree, logcred = compute_sumtree(treesummary, args, wt_count_burnin_filename_list)
        sumtree = root_sumtree(sumtree, args)
        sumtree = annotate_sumtree_root(sumtree, args)
        sumtree = set_sumtree_blen(sumtree, args)
        sumtree_status_message = print_sumtree(sumtree, args)
        output.info(sumtree_status_message)
        
        if args.treeprobs:
            trproblist = compute_trprobs(treesummary, args):
            trprobs_status_message = print_trprobs(trproblist, args)
            output.info(trprobs_status_message)
        
        sumdict = compute_summary_variables(sumtree, treesummary, start, memory1, ave_std, args)
        print_result_summary(sumdict, sumtree, args, output)

    except Exception as error:
        handle_error(error, args.verbose)

####################################################################################
####################################################################################

class ProgressBar:
    def __init__(self, total_trees, burnin, quiet=False):
        self.quiet = quiet
        self.total_trees = total_trees
        self.burnin = burnin
        self.processed_trees = 0
        self.progscale = "0      10      20      30      40      50      60      70      80      90     100"
        self.progticks = "v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v"
        self.ndots = len(self.progticks)
        self.n_tot = total_trees - burnin
        self.trees_per_dot = self.n_tot / self.ndots
        self.n_dotsprinted = 0

        # Print progress bar header
        if not self.quiet:
            print("\n\n   Processing trees:")
            print(f"   {self.progscale}")
            print(f"   {self.progticks}")
            print("   ", end="")
            sys.stdout.flush()

    def update(self):
        """ Update the progress bar based on the number of processed trees. """
        if not self.quiet:
            self.processed_trees += 1
            n_dots_expected = math.floor(self.processed_trees / self.trees_per_dot)
            if self.n_dotsprinted < n_dots_expected:
                n_missing = n_dots_expected - self.n_dotsprinted
                sys.stdout.write("*" * n_missing)
                sys.stdout.flush()
                self.n_dotsprinted += n_missing

    def complete(self):
        """ Ensure all dots are printed at the end if they haven't been already. """
        if not self.quiet:
            if self.n_dotsprinted < self.ndots:
                n_missing = self.ndots - self.n_dotsprinted
                sys.stdout.write("*" * n_missing)
                sys.stdout.flush()

####################################################################################
####################################################################################

class OutputManager:
    def __init__(self, args):
        self.quiet = args.quiet
        self.verbose = args.verbose

    def info(self, message="", margin=3, end="\n"):
        if not self.quiet:
            print(f"{margin * ' '}{message}", end=end)

    def force(self, message="", margin=3, end="\n"):
        if not self.quiet:
            print(f"{margin * ' '}{message}", end=end)
            sys.stdout.flush()

    def warning(self, message):
        if not self.quiet and self.verbose:
            print(f"{self.margin * ' '}[WARNING] {message}")
            
####################################################################################
####################################################################################

def setup_output_directory(outbase):
    outbase.parent.mkdir(parents=True, exist_ok=True)
    
####################################################################################

def track_memory_usage(pid):
    return pid.memory_full_info().rss
  
####################################################################################      

def parse_commandline(commandlist):
    # Python note: "commandlist" is to enable unit testing of argparse code
    # Will be "None" when run in script mode, and argparse will then automatically take values from sys.argv[1:]
    # https://jugmac00.github.io/blog/testing-argparse-applications-the-better-way/

    parser = build_parser()
    args = parser.parse_args(commandlist)

    if args.version:
        config = configparser.ConfigParser()
        config.read('setup.cfg')
        try:
            print(config['metadata']['version'])
            exit()
        except KeyError:
            print("Unknown")
            exit()

    if not (args.con or args.all or args.mcc or args.mbc):
        parser.error("\n    One of --con, --all, --mcc, or --mbc must be specified\n    (to select type of summary tree)")

    if not (args.noblen or args.biplen or args.meandepth or args.cadepth):
        parser.error("\n    One of --noblen, --biplen, --meandepth, or --cadepth must be specified\n    (to select how to estimate branch lengths)")

    if not args.infilelist and not args.fileweights:
        parser.error("Please list one or more tree files.")

    # If output basename is not set: use stem of infilenames minus all suffixes
    if not args.outbase:
        if args.infilelist:
            infilepath = args.infilelist[0]
        else:
            wt, infilepath = args.fileweights[0]
        args.outbase = Path(infilepath.stem.split('.')[0])

    if len(args.burninfrac) == 1:
        burnin_value = args.burninfrac[0]
        args.burninfrac = [burnin_value] * len(args.infilelist)
    elif len(args.burninfrac) != len(args.infilelist):
        parser.error("either provide one burnin value, or one value per input file")

    if any(x < 0 or x > 1 for x in args.burninfrac):
        parser.error("option -b: NUM must be between 0.0 and 1.0")

    if args.treeprobs and (args.treeprobs > 1 or args.treeprobs < 0):
        parser.error(f"option -t: NUM must be between 0.0 and 1.0 (provided value: -t {args.treeprobs})")

    if args.infilelist:
        nfiles = len(args.infilelist)
    else:
        nfiles = len(args.fileweights)
    if args.std and nfiles==1:
        parser.error("cannot compute standard deviation (option -s) from one tree file")

    if args.quiet:
        args.nowarn = True
        
    if args.rootog and not (args.outgroup or args.ogfile):
        parser.error("Option --rootog requires specifying outgroup using either --og or --ogfile")

    if args.mcc and (args.rootog or args.rootmid or args.rootminvar):
        parser.error("MCC tree is not compatible with any of these rooting methods: --rootmid, --rootminvar, --rootog")

    # Bipartitions need to be tracked in these situations
    if args.con or args.all or args.mbc or args.biplen:
        args.trackbips = True
    else:
        args.trackbips = False

    # Clades need to be tracked in these situations:
    if args.mcc or args.meandepth or args.cadepth:
        args.trackclades = True
    else:
        args.trackclades = False

    # Root needs to be tracked in these situations:
    if args.rootcred:   
        args.trackroot = True
    else:
        args.trackroot = False

    # Branch lengths need to be tracked if biplen==True
    if args.biplen:
        args.trackblen = True
    else:
        args.trackblen = False

    # Node depths need to be tracked if meandepth is set
    if args.meandepth:
        args.trackdepth = True
    else:
        args.trackdepth = False

    if args.ogfile:
        args.outgroup = read_outgroup(args.ogfile)

    if (args.outgroup or args.rootmid or args.rootminvar):
        args.actively_rooted = True
    else:
        args.actively_rooted = False

    return args

####################################################################################

def build_parser():

    parser = argparse.ArgumentParser(description = "Computes summary tree and statistics from set of phylogenetic trees")

    parser.add_argument('--version', action='store_true', dest="version",
                        help="show the program's version number and exit")

    ####################################################################################

    inout_grp = parser.add_argument_group("Input and output")
    infile_excl = inout_grp.add_mutually_exclusive_group()

    inout_grp.add_argument("--informat", action="store", metavar="FORMAT",
                      choices=["nexus", "newick"], default="nexus",
                      help="format of input tree files: %(choices)s [default: %(default)s]")

    inout_grp.add_argument("--outformat", action="store", metavar="FORMAT",
                      choices=["nexus", "newick"], default="nexus",
                      help="format of output tree file: %(choices)s [default: %(default)s]")

    infile_excl.add_argument("-i", action="append", dest='infilelist', metavar='FILE', type=Path,
                        help="input FILE(s) containing phylogenetic trees (repeat -i FILE option for each input file)")

    infile_excl.add_argument("-w", action="append", dest="fileweights",
                        nargs=2, metavar=("WEIGHT", "FILE"),
                        help="input FILEs with specified weights (repeat -w WEIGHT FILE option for each input file)")

    inout_grp.add_argument("--autow", action="store_true", dest="autoweight",
                     help="automatically assign file weights based on tree counts, so all files have equal impact "
                         + "(default is for all trees, not files, to be equally important)")

    inout_grp.add_argument("--basename", action="store", type=Path, dest="outbase", metavar="NAME",
                      help="base name of output files (default: derived from input file)")

    inout_grp.add_argument("-n", action="store_true", dest="nowarn",
                      help="no warning when overwriting files")

    inout_grp.add_argument("-v", action="store_true", dest="verbose",
                      help="verbose: show full traceback in the event of failed python execution")

    inout_grp.add_argument("-q", action="store_true", dest="quiet",
                      help="quiet: don't print progress indication to terminal window. NOTE: also turns on the -n option")


    ####################################################################################

    sumtype_grp = parser.add_argument_group("Type of summary tree (pick one option)")
    sumtype_excl = sumtype_grp.add_mutually_exclusive_group()

    sumtype_excl.add_argument("--con", action="store_true",
                              help="majority rule consensus tree")

    sumtype_excl.add_argument("--all", action="store_true",
                              help="majority rule consensus tree with all compatible bipartitions added")

    sumtype_excl.add_argument("--mcc", action="store_true",
                              help="Maximum Clade Credibility (MCC) tree. "
                              + "The MCC tree is determined by inspecting tree samples and selecting the "
                              + "tree that has the highest product of clade frequencies (= highest sum of "
                              + "log of clade frequencies). The MCC tree is therefore a tree that has been "
                              + "observed in the pool of tree samples, differing from the consensus tree "
                              + "which typically does not match any individual sample. "
                              + "NOTE 1: only meaningful if input trees are estimated using clock model. "
                              + "NOTE 2: by default, the MCC tree uses the rooting of the specific tree sample. "
                              + "This will often (but not always) correspond to the "
                              + "bipartition where the root is most commonly found in the input trees.")

    sumtype_excl.add_argument("--mbc", action="store_true",
                              help="Maximum Bipartition Credibility (MBC) tree. "
                              + "The MBC tree is similar to the MCC tree "
                              + "but counting bipartitions instead of clades, i.e. ignoring rooting "
                              + "(two input trees can have the same set of bipartitions, but be rooted "
                              + "in different locations).")

    ####################################################################################

    blen_grp = parser.add_argument_group(title= "Estimation of branch lengths (pick one option)")
    blen_excl = blen_grp.add_mutually_exclusive_group()
    blen_excl.add_argument("--noblen", action="store_true",
                      help="Do not set branch lengths (only the topology and branch- or clade-"
                          + "support of the summary tree are estimated). ")

    blen_excl.add_argument("--biplen", action="store_true",
                      help="Set branch lengths in summary tree based on average for corresponding "
                          + "leaf bipartitions:"
                          + "each branch in tree corresponds to a bipartition of the leaves "
                          + "into two groups. Branch lenghts are set to the mean of the length of the"
                          + "corresponding bipartition across all input trees.")

    blen_excl.add_argument("--meandepth", action="store_true",
                      help="set node depth for each clade to mean node depth observed for that "
                           + "clade among input trees "
                           + "(and branch lengths are then based on these depths). "
                           + "Warning: option is intended "
                            + "for input trees estimated using a clock model. "
                            + "It requires that all clades in the summary tree have "
                            + "been observed in the input trees, and may fail "
                            + "for some rootings."
                           + "NOTE: mean is computed across trees where the specific, monophyletic clade "
                           + "is present, and may therefore be based on very few (down to one) values. "
                           + "NOTE 2: may result in negative branch lengths. ")

    blen_excl.add_argument("--cadepth", action="store_true",
                      help="'Common Ancestor depth'. Same as option '--height ca' in treeannotator. "
                           + "Uses all trees in input set when determining node-depths. "
                           + "For a given clade: (1) Find the most recent "
                           + "common ancestor of the leaves in that clade in each of the input trees. "
                           + "(2) Compute node-depth of clade as the mean of the depths of these MRCAs. "
                           + "This is different from --meandepth where only "
                           + "trees with that exact clade are included when computing the mean. "
                           + "Warning: option is intended "
                            + "for input trees estimated using a clock model. "
                            + "It requires that all clades in the summary tree have "
                            + "been observed in the input trees, and may fail "
                            + "for some rootings.")

    ####################################################################################

    root_grp = parser.add_argument_group("Rooting of summary tree")
    
    root_excl = root_grp.add_mutually_exclusive_group()
    
    root_excl.add_argument("--rootmid", action="store_true",
                      help="perform midpoint rooting of summary tree")
                      
    root_excl.add_argument("--rootminvar", action="store_true",
                      help="perform minimum variance rooting of summary tree")

    root_excl.add_argument("--rootog", action="store_true",
                      help="root summary tree on outgroup (requires specifying outgroup using --og or --ogfile)")
                      
    og_group = root_grp.add_mutually_exclusive_group()
    
    og_group.add_argument('--og', dest="outgroup", metavar="NAME", nargs="+", default=None,
                      help='specify outgroup taxon/taxa on command-line')
                      
    og_group.add_argument('--ogfile', action="store", metavar="FILE", default=None,
                      help="specify outgroup taxon/taxa in file (one name per line)")
                      
    root_grp.add_argument("--rootcred", action="store_true",
                          help=("compute root credibility for all possible rooting locations and add a 'rootcred' "
                                "attribute to branches in the summary tree. If an outgroup is specified: track which "
                                "branch (bipartition) the outgroup attaches to in each input tree and report the "
                                "frequency of rooting there. If no outgroup is specified: assume input trees are "
                                "rooted and track root frequencies on branches. The cumulated root credibility may be "
                                "less than 100%% if some root locations are not present in the summary tree."))
                         
    ####################################################################################

    bayes_grp = parser.add_argument_group("Bayesian phylogeny options")

    bayes_grp.add_argument("-b", dest="burninfrac", metavar="NUM", type=float, default=[0], nargs='+',
                           help="burnin: fraction of trees to discard [0 - 1; default: %(default)s]. "
                           + "Either one value (used on all input files), or one value per input file.")

    bayes_grp.add_argument("-t", type=float, dest="treeprobs", metavar="NUM",
                      help="compute tree probabilities, report NUM percent credible interval [0 - 1]")

    bayes_grp.add_argument("-s", action="store_true", dest="std",
                      help="compute average standard deviation of split frequencies (ASDSF)")

    bayes_grp.add_argument("-f", type=float, dest="minfreq", metavar="NUM", default=0.1,
                      help="Minimum frequency for including bipartitions in report and in computation of ASDSF [default: %(default)s]")

    ####################################################################################

    outformat_grp = parser.add_argument_group("Output to terminal and files")


    ####################################################################################

    other_grp = parser.add_argument_group("Other options")

    other_grp.add_argument("--nolabel", action="store_true", dest="nolabel",
                     help=argparse.SUPPRESS # do not print branch labels (=clade probabilities) on summary tree
                     )

    ####################################################################################

    return parser

####################################################################################

def parse_infilelist(args):

    # If only unweighted filenames are given:
    # Reformat list of filenames into (weight, filename) tuple format expected by program
    # Set all weights to 1
    if args.infilelist:
        wt_file_list = [(1, filename) for filename in args.infilelist]

    # If only weighted filenames are listed:
    # Reformat list of tuples such that weight is in float (not string).
    # Normalize weights so their average is one.
    else:
        wt_file_list = []

        # Attempt to convert weight string to float. Print sensible error message if this fails
        for (wt_string, filename) in args.fileweights:
            try:
                wt = float(wt_string)
            except ValueError:
                msg = f'Invalid file weight: "{wt_string}" - value has to be a real number.'
                raise Exception(msg)
            wt_file_list.append((wt, filename))

        # Normalize weights, build final weight/file list:
        wtsum = 0.0
        n_files = len(wt_file_list)
        for (wt, filename) in wt_file_list:
            wtsum += wt
        wt_avg = wtsum / n_files
        wt_file_list = [(wt / wt_avg, filename) for (wt, filename) in wt_file_list]

    return wt_file_list

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
    # count "= (", "=  and "tree "
    # Add the values that are not 0 to list, choose minimum as count
    # Should be robust to most variations, but should check at end of sumt...
    n_eqparen = count_bytestring(filename, b"= (")
    n_treestr = count_bytestring(filename, b"tree ")
    countlist = [n_semicolons, n_eqparen, n_treestr]
    notzero = [val for val in countlist if val>0]
    return min(notzero)

####################################################################################

def count_trees(wt_file_list, args, output):

    count_list = []
    burnin_list = []
    n_postburnin = []
    for i,(wt, filename) in enumerate(wt_file_list):
        treelist = []
        output.force(f"Counting trees in file {str(filename):<40}", end="")
        n_tot = fast_treecount(filename, args)
        output.force(f"{n_tot:>15,d}")
        burnin = int(args.burninfrac[i] * n_tot)
        count_list.append(n_tot)
        burnin_list.append(burnin)
        n_postburnin.append(n_tot - burnin)

    # If automatic weighting requested: Compute new weights
    if args.autoweight:
        countsum = sum(n_postburnin)
        countavg = countsum / len(n_postburnin)
        new_wt_list = [countavg / count for count in n_postburnin]  # Normalized so avg=1

    # Construct final combined wt + count + burnin + filename list
    wt_count_burnin_filename_list = []
    for i in range(len(wt_file_list)):
        filename = wt_file_list[i][1]
        count = count_list[i]
        burnin = burnin_list[i]
        if args.autoweight:
            wt = new_wt_list[i]
        else:
            wt = wt_file_list[i][0]

        wt_count_burnin_filename_list.append((wt, count, burnin, filename))

    n_trees_analyzed = sum(count_list) - sum (burnin_list)
    return (n_trees_analyzed, wt_count_burnin_filename_list)

####################################################################################

def process_trees(wt_count_burnin_filename_list, args, output):

    treesummarylist = []
    interner = pt.Interner()

    for i, (weight, count, burnin, filename) in enumerate(wt_count_burnin_filename_list):
        output.force()
        output.force("Analyzing file: {} (Weight: {:5.3f})".format(filename, weight))

        # Open treefile. Discard (i.e., silently pass by) the requested number of trees
        if args.informat == "nexus":
            treefile = pt.Nexustreefile(filename, interner=interner)
        else:
            treefile = pt.Newicktreefile(filename, interner=interner)
        for j in range(burnin):
            treefile.readtree(returntree=False)
        output.info()
        output.info(f"Discarded {burnin:,} of {count:,} trees (burnin fraction={args.burninfrac[i]:.2f})")

        # Instantiate Treesummary.
        trackbips = args.trackbips
        trackclades = args.trackclades
        trackroot = args.trackroot
        trackblen = args.trackblen
        trackdepth = args.trackdepth
        if args.mcc or args.mbc or args.treeprobs:
            treesummary = pt.BigTreeSummary(store_trees=args.treeprobs,
                                            trackbips=trackbips, trackclades=trackclades, trackroot=trackroot,
                                            trackblen=trackblen, trackdepth=trackdepth)
        else:
            treesummary = pt.TreeSummary(trackbips=trackbips, trackclades=trackclades, trackroot=trackroot,
                                         trackblen=trackblen, trackdepth=trackdepth)

        # Initialize the progress bar
        progress = ProgressBar(total_trees=count, burnin=burnin, args.quiet)

        # Read post-burnin trees from file, add to treesummary, print progress bar
        for j in range(burnin, count):
            tree = treefile.readtree(returntree=True)
            treesummary.add_tree(tree, weight)
            progress.update()
            del tree

        # Ensure the progress bar completes at the end
        progress.complete()

        treesummarylist.append(treesummary)
        output.info()

    return treesummarylist

##########################################################################################

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

    if args.mcc:
        ave_std = compute_converge_clades(treesummarylist, args)
    else:
        ave_std = compute_converge_biparts(treesummarylist, args)
    return ave_std

##########################################################################################

def compute_converge_biparts(treesummarylist, args):
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

def compute_converge_clades(treesummarylist, args):
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

def  merge_treesummaries(treesummarylist):

    treesummary = treesummarylist[0]
    for treesummary2 in treesummarylist[1:]:
        treesummary.update(treesummary2)
        del treesummary2
    return treesummary

##########################################################################################

def compute_and_print_contree(treesummary, args, wt_count_burnin_filename_list):

    if args.mcc:
        sys.stdout.write("\n   Finding Maximum Clade Credibility tree...")
        sys.stdout.flush()
        contree, logcred = treesummary.max_clade_cred_tree()
    elif args.mbc:
        sys.stdout.write("\n   Finding Maximum Bipartition Credibility tree...")
        sys.stdout.flush()
        contree, logcred = treesummary.max_bipart_cred_tree()
    else:
        sys.stdout.write("\n   Computing consensus tree...")
        sys.stdout.flush()
        contree = treesummary.contree(allcompat=args.all)
        logcred = treesummary.log_bipart_credibility(contree.topology())
    sys.stdout.write("done.\n")
    sys.stdout.flush()

    if args.rootog:
        contree.rootout(args.outgroup)
    elif args.rootmid:
        contree.rootmid()
    elif args.rootminvar:
        contree.rootminvar()

    if args.rootcred:
        if args.actively_rooted or args.mcc:
            contree.rootcred = treesummary.compute_rootcred(contree)
        contree = treesummary.set_rootcredibility(contree)  # Add branch attribute with root credibilities

    if args.meandepth:
        contree = treesummary.set_mean_node_depths(contree)
    elif args.cadepth:
        sys.stdout.write("   Computing common ancestor depths...")
        sys.stdout.flush()
        contree = treesummary.set_ca_node_depths(contree, wt_count_burnin_filename_list)
        sys.stdout.write("done.\n")
    elif args.biplen and args.mcc:
        contree = treesummary.set_mean_biplen(contree)

    # If root is bifurcation and one child is leaf: Remove label from other child-branch
    # (In this case both branches from root are the same bipartition, so not useful to label internal branch part)
    # Python note: not sure this generalizes to all tree types. think...
    rootkids = contree.children(contree.root)
    if (len(rootkids) == 2) and (rootkids & contree.leaves):
        rootkids = list(rootkids)
        if rootkids[0] in contree.leaves:
            node2 = rootkids[1]
        else:
            node2 = rootkids[0]
        contree.set_branch_attribute(contree.root, node2, "label", "")

    if args.trackblen or args.trackdepth or args.cadepth:
        printdist = True
    else:
        printdist = False
        
def print_sumtree(sumtree, args):

    printdist = args.trackblen or args.trackdepth or args.cadepth
    printlabels = not args.nolabel

    metacomment_fields = ["posterior", "length"]
    if args.trackroot:
        metacomment_fields.append("rootcred")

    newick_prob_tree = contree.newick(labelfield="label", printdist=printdist, printlabels=printlabels)
    
    if args.mbc:
        confilename = args.outbase.parent / (args.outbase.name + ".mbc")
    elif args.mcc:
        confilename = args.outbase.parent / (args.outbase.name + ".mcc")
    elif args.all:
        confilename = args.outbase.parent / (args.outbase.name + ".all")
    else:
        confilename = args.outbase.parent / (args.outbase.name + ".con")
    if args.nowarn:
        confile = open(confilename, "w")
    elif confilename.is_file():
        overwrite = input(f"\n   File {confilename} already exists.\n   Overwrite (y/n): ")
        if overwrite== "y":
            confile = open(confilename, "w")            # Overwrite
            print("   Overwriting file {}\n".format(confilename))
        else:
            confile.write("#NEXUS\n\n")
            confile.write("begin trees;\n")
            confile.write("   [In this tree branch labels indicate the posterior probability of the bipartition corresponding to the branch.]\n")
            confile.write("   tree prob = ")
            confile.write(newick_prob_tree)
            confile.write("\nend;\n")
        
    if args.mbc:
        return f"Maximum bipartition credibility tree written to {confilename}"
    elif args.mcc:
        return f"Maximum clade credibility tree written to {confilename}"
    else:
        return f"Consensus tree written to {confilename}"
        
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

def print_trprobs(trproblist, args):
    topofilename = args.outbase.parent / (args.outbase.name + ".trprobs")
    if args.nowarn:
        topofile = open(topofilename, "w")
    elif topofilename.is_file():
        overwrite = input(f"\n   File {topofilename} already exists.\n   Overwrite (y/n): ")
        if overwrite== "y":
            topofile = open(topofilename, "w")            # Overwrite
            print(f"   Overwriting file {topofilename}\n")
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

def compute_summary_variables(sumtree, treesummary, start, memory1, ave_std, args):
    sumdict = {}
    
    memory2 = track_memory_usage(pid)
    sumdict["memorymax"] = max(memory1, memory2) 
    
    sumdict["ave_std"] = ave_std

    if args.treeprobs:
        if args.trackclades:
            sumdict["n_topo_seen"] = len(treesummary.cladetoposummary)
        elif args.trackbips:
            sumdict["n_topo_seen"] = len(treesummary.biptoposummary)

    nrootkids = len(sumtree.children(sumtree.root))
    if nrootkids == 2:
        sumdict['rootdegree'] = "bifurcation"
    elif nrootkids == 3:
        sumdict['rootdegree'] = "trifurcation"
    else:
        sumdict['rootdegree'] = "multifurcation"

    n_leaves = len(treesummary.leaves)
    sumdict['n_leaves'] = n_leaves

    if args.mcc:
        sumdict['n_uniq_groupings'] = len(treesummary.cladesummary) - n_leaves
        sumdict['theo_max_groups'] = n_leaves - 1
    else:
        sumdict['n_uniq_groupings'] = len(treesummary.bipartsummary) - n_leaves
        sumdict['theo_max_groups'] = n_leaves - 3

    sumdict['n_internal_biparts'] = sumtree.n_bipartitions()

    if args.mcc:
        sumdict['treetype'], sumdict['branchtype'], sumdict['space'] = "MCC", "clade", " " * 7
    elif args.mbc:
        sumdict['treetype'], sumdict['branchtype'], sumdict['space'] = "MBC", "bipartition", " " * 1
    else:
        sumdict['treetype'], sumdict['branchtype'], sumdict['space'] = "Consensus", "bipartition", " " * 1

    time_spent = time.time() - start
    sumdict["h"] = int(time_spent/3600)
    sumdict["m"] = int((time_spent % 3600)/60)
    sumdict["s"] = int(time_spent % 60)

    return sumdict

##########################################################################################

def print_result_summary(sumdict, sumtree, args):
    
    sd = sumdict
    
    # Information about bipartitions, clades and topologies
    output.info()
    output.info(f"Number of leaves on input trees: {sd['n_leaves']:>7,d}")
    if args.treeprobs:
        output.info(f"Different topologies seen: {sd['n_topo_seen']:>13,d}")
        output.info(f"Different {sd['branchtype']}s seen:{sd['space']}{sd['n_uniq_groupings']:>11,d} (theoretical maximum: {sd['theo_max_groups'] * sd['n_topo_seen']:,d})")
    else:
        output.info(f"Different {sd['branchtype']}s seen:{sd['space']}{sd['n_uniq_grouping']s:>11,d} (theoretical maximum: {sd['theo_max_groups'] * sd['n_trees_analyzed']:,d})")
    output.info(f"{'Bipartitions in ' + {sd['treetype']} + ' tree:':<34}{sd['n_internal_biparts']:>6,d} (theoretical maximum: {sd['theo_max_groups']:,d})")

    if sd['n_internal_biparts'] < sd['theo_max_bips']:
        output.info("(tree contains polytomies)", margin=44)
    else:
        output.info("(tree is fully resolved - no polytomies)", margin=44)

    # Information about rooting
    if not (args.actively_rooted or args.mcc):
        output.info()
        output.info(f"{sd['treetype']} tree has not been explicitly rooted")
        output.info(f"Tree has been rooted at random internal node; root is at {sd['rootdegree']}")
    else:
        if args.rootog:
            output.info()
            output.info(f"{sd['treetype']} tree has been rooted based on outgroup")
        elif args.rootmid:
            output.info()
            output.info(f"{sd['treetype']} tree has been midpoint-rooted")
        elif args.rootminvar:
            output.info()
            output.info(f"{sd['treetype']} tree has been rooted using minimum variance-rooting")
        elif args.mcc:
            output.info()
            output.info(f"MCC tree rooted at original root of tree sample having highest clade credibility")

    if args.rootcred:
        if args.actively_rooted or args.mcc:
            output.info(f"Root credibility (frequency of root bipartition in input trees):       {sumtree.rootcred * 100:.1f}%")
        output.info(f"Cumulated root credibility (sum of rootcred for all branches in tree): {sumtree.cumulated_rootcred * 100:.1f}%")
        

    # Information about branch lengths
    if args.meandepth:
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

    # Information about log credibility
    if args.mbc or (args.mcc and not args.actively_rooted):
        output.info()
        output.info(f"Highest log {sd['branchtype']} credibility:  {sd['logcred']:.6g}")
    else:
        output.info()
        output.info(f"Log {sd['branchtype']} credibility:  {sd['logcred']:.6g}")

    if args.std:
        output.info(f"Average standard deviation of split frequencies: {sd['ave_std']:.6f}")

    output.info()
    output.info(f"Done. {sd['n_trees_analyzed']:,d} trees analyzed.\n   Time spent: {sd['h']:d}:{sd['m']:02d}:{sd['s']:02d} (h:m:s)")

    if memorymax > 1E9:
        output.info("Max memory used: {:,.2f} GB.".format( sd['memorymax']  / (1024**3) ))
    else:
        output.info("Max memory used: {:,.2f} MB.".format( sd['memorymax']  / (1024**2) ))

##########################################################################################
##########################################################################################

if __name__ == "__main__":
    main()
    # import cProfile
    # cProfile.run('main()', 'tmp/profile.pstats')
