####################################################################################
####################################################################################

import phylotreelib as treelib
import argparse, os, sys, time, math, copy, psutil
from itertools import (takewhile,repeat)
from operator import itemgetter
from pathlib import Path
import gc

gc.disable()        # Faster. Assume no cyclic references will ever be created

def main(commandlist=None):
    # Python note: "commandlist" is to enable unit testing of argparse code
    # https://jugmac00.github.io/blog/testing-argparse-applications-the-better-way/
    start=time.time()
    args = parse_commandline(commandlist)
    try:
        args.outbase.parent.mkdir(parents=True, exist_ok=True) # Create intermediate dirs
        wt_file_list = parse_infilelist(args)
        if args.quiet:
            sys.stdout = open(os.devnull, 'w')
        n_trees_analyzed, wt_count_burnin_filename_list = count_trees(wt_file_list, args)
        if args.verbose:
            pid = psutil.Process(os.getpid())
            memory1 = pid.memory_full_info().rss
        treesummarylist = process_trees(wt_count_burnin_filename_list, args)
        if args.std:
            compute_and_print_converge_stats(treesummarylist, args.minfreq)
        treesummary = merge_treesummaries(treesummarylist)
        n_leafs = len(treesummary.leaves)
        total_unique_internal_biparts = len(treesummary.bipartsummary) - n_leafs
        treesummary.add_branchid()
        contree, logbipcred = compute_and_print_contree(treesummary, args)
        compute_and_print_biparts(treesummary, args)
        theo_maxbip_internal_unrooted = n_leafs - 3
        n_internal_biparts = contree.n_bipartitions()
        if args.treeprobs:
            compute_and_print_trprobs(treesummary, args)
            n_topo_seen = len(treesummary.toposummary)
        stop=time.time()

        if args.verbose:
            memory2 = pid.memory_full_info().rss
            memorymax = max(memory1, memory2)

        time_spent=stop-start
        h = int(time_spent/3600)
        m = int((time_spent % 3600)/60)
        s = int(time_spent % 60)
        print("\n   Done. {:,d} trees analyzed.\n   Time spent: {:d}:{:02d}:{:02d} (h:m:s)\n".format(n_trees_analyzed, h, m, s))
        if args.verbose:
            if args.mbc:
                treetype = "MBC"
            else:
                treetype = "Consensus"

            if memorymax > 1E9:
                print("   Max memory used: {:,.2f} GB.".format( memorymax  / (1024**3) ))
            else:
                print("   Max memory used: {:,.2f} MB.".format( memorymax  / (1024**2) ))

            if args.treeprobs:
                print("\n   Different topologies seen: {:10,d}".format(n_topo_seen))
            print("   Different bipartitions seen: {:8,d} (theoretical maximum: {:,d})".format(
                                                    total_unique_internal_biparts, theo_maxbip_internal_unrooted * n_topo_seen))
            if args.rooted:
                theo_maxbip_internal = theo_maxbip_internal_unrooted + 1
            else:
                theo_maxbip_internal = theo_maxbip_internal_unrooted
            print("   {:<34}".format(f"Bipartitions in {treetype} tree:"), end="")
            print(f"{n_internal_biparts:3,d} (theoretical maximum: {theo_maxbip_internal:,d})")

            if n_internal_biparts < theo_maxbip_internal:
                print("                                         (tree contains polytomies)")
            else:
                print("                                         (tree is fully resolved - no polytomies)")

            if args.rooted:
                print(f"\n   {treetype} tree has been explicitly rooted")
                print(f"   (Root is at bifurcation)")
            else:
                print(f"\n   {treetype} tree has not been explicitly rooted")
                print(f"   (Tree has been rooted at random internal node; root is at trifurcation)")

            if args.mbc:
                print(f"\n   Highest Log Bipartition Credibility:  {logbipcred:.4g}")
            else:
                print(f"\n   Log Bipartition Credibility:  {logbipcred:.4g}")


    except treelib.TreeError as error:
        print("Execution failed:\n")
        if args.verbose:
            import traceback
            traceback.print_exc(file=sys.stdout)
        else:
            print("Error: ", error)

        sys.exit()

####################################################################################
####################################################################################

def parse_commandline(commandlist):
    # Python note: "commandlist" is to enable unit testing of argparse code
    # Will be "None" when run in script mode, and argparse will then automatically take values from sys.argv[1:]

    parser = build_parser()
    args = parser.parse_args(commandlist)

    if args.infilelist and args.fileweights:
        parser.error("When using option -w all input files need to have a weight specified")

    if not args.infilelist and not args.fileweights:
        parser.error("Please list one or more tree files.")

    if not any([args.con, args.all, args.mbc]):
        parser.error("Please select type of summary tree using one of these options: --con, --all, --mbc")

    # If output basename is not set: determine semi-intelligently from infilenames:
    if not args.outbase:
        if args.infilelist:
            infilename = args.infilelist[0]
        else:
            wt, infilename = args.fileweights[0]
        if infilename.endswith(".run1.t"):
            args.outbase = outname[:-7]
        elif infilename.endswith(".t"):
            args.outbase = outname[:-2]
        else:
            args.outbase = infilename.stem

    if args.burninfrac > 1 or args.burninfrac < 0:
        parser.error("option -b: NUM must be between 0.0 and 1.0")

    if args.treeprobs and (args.treeprobs > 1 or args.treeprobs < 0):
        parser.error("option -t: NUM must be between 0.0 and 1.0")

    if args.infilelist:
        nfiles = len(args.infilelist)
    else:
        nfiles = len(args.fileweights)
    if args.std and nfiles==1:
        parser.error("cannot compute standard deviation from one tree file")

    if args.quiet:
        args.nowarn = True

    if args.rootfile:
        args.outgroup = read_outgroup(args.rootfile)

    return args

####################################################################################
####################################################################################

def build_parser():

    parser = argparse.ArgumentParser(description = "Computes summary tree and statistics from set of phylogenetic trees")

    infilegroup = parser.add_argument_group("Input tree files")

    infilegroup.add_argument('infilelist', nargs='*', metavar='INFILE', type=Path,
                        help="input FILE(s) containing phylogenetic trees (can list several files)")

    infilegroup.add_argument("-w", action="append", dest="fileweights",
                        nargs=2, metavar=("WEIGHT", "INFILE"),
                        help="input FILEs with specified weights (repeat -w option for each input file)")

    infilegroup.add_argument("--autow", action="store_true", dest="autoweight",
                     help="automatically assign file weights based on tree counts, so all files have equal impact"
                         + "(default is for all trees, not files, to be equally important)")

    infilegroup.add_argument("-i", action="store", dest="informat", metavar="FORMAT",
                      choices=["nexus", "newick"], default="nexus",
                      help="format of input files: %(choices)s [default: %(default)s]")

    ####################################################################################

    sumtypegroup = parser.add_argument_group("Type of summary tree")

    sumtypecommands = sumtypegroup.add_mutually_exclusive_group()

    sumtypecommands.add_argument("--con", action="store_true", default=True,
                              help="majority rule consensus tree")

    sumtypecommands.add_argument("--all", action="store_true",
                              help="majority rule consensus tree with all compatible bipartitions added")

    sumtypecommands.add_argument("--mbc", action="store_true",
                              help="Maximum Bipartition Credibility (MBC) tree. "
                              + "MBC is similar to MCC (Maximum Clade Credibility) tree "
                              + "but counting bipartitions instead of clades, i.e. ignoring rooting. "
                              + "Additionally, branch lengths are estimated from branch lengths of bipartitions "
                              + "and not from node depths (i.e., again ignoring rooting)")

    ####################################################################################

    bayesgroup = parser.add_argument_group("Bayesian phylogeny options")

    bayesgroup.add_argument("-b", type=float, dest="burninfrac", metavar="NUM", default=0.25,
                      help="burnin: fraction of trees to discard [0 - 1; default: %(default)s]")

    bayesgroup.add_argument("-t", type=float, dest="treeprobs", metavar="NUM",
                      help="compute tree probabilities, report NUM percent credible interval [0 - 1]")

    bayesgroup.add_argument("-s", action="store_true", dest="std",
                      help="compute average standard deviation of split frequencies (ASDSF)")

    bayesgroup.add_argument("-f", type=float, dest="minfreq", metavar="NUM", default=0.1,
                      help="Minimum frequency for including bipartitions in report and in computation of ASDSF [default: %(default)s]")

    ####################################################################################

    outformatgroup = parser.add_argument_group("Output to terminal and files")

    outformatgroup.add_argument("-n", action="store_true", dest="nowarn",
                      help="no warning when overwriting files")

    outformatgroup.add_argument("-v", action="store_true", dest="verbose",
                      help="verbose: more information, longer error messages")

    outformatgroup.add_argument("-q", action="store_true", dest="quiet",
                      help="quiet: don't print progress indication to terminal window. NOTE: also turns on the -n option")

    outformatgroup.add_argument("--basename", action="store", type=Path, dest="outbase", metavar="NAME",
                      help="base name of output files (default: derived from input file)")

    ####################################################################################

    rootgroup = parser.add_argument_group("Rooting of summary tree")
    rootcommands = rootgroup.add_mutually_exclusive_group()
    rootcommands.add_argument("--rootmid", action="store_true", dest="midpoint",
                      help="perform midpoint rooting of tree")
    rootcommands.add_argument("--rootminvar", action="store_true", dest="minvar",
                      help="perform minimum variance rooting of tree")

    rootcommands.add_argument("-r", dest="outgroup", metavar="TAXON", nargs="+", default=None,
                      help="root consensus tree on specified outgroup taxon/taxa")

    rootcommands.add_argument("--rootfile", action="store", metavar="FILE", default=None,
                      help="root consensus tree on outgroup taxa listed in file (one name per line)")

    return parser

####################################################################################
####################################################################################

def parse_infilelist(args):

    # If only unweighted filenames are given:
    # Reformat list of filenames into (weight, filename) tuple format expected by program
    # Set all weights to 1/n_files
    if args.infilelist:
        nfiles = len(args.infilelist)
        wt_file_list = [(1/nfiles, filename) for filename in args.infilelist]

    # If only weighted filenames are listed:
    # reformat list of tuples such that weight is in float (not string). Normalize weights so they sum to one.
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
        for (wt, filename) in wt_file_list:
            wtsum += wt
        wt_file_list = [(wt/wtsum, filename) for (wt, filename) in wt_file_list]

    return wt_file_list

####################################################################################
####################################################################################

def read_outgroup(rootfile):
    infile = open(rootfile, "r")
    outgroup = []
    for line in infile:
        leaf = line.strip()
        outgroup.append(leaf)
    infile.close()
    return outgroup

####################################################################################
####################################################################################

def count_trees(wt_file_list, args):

    # fast counting of specific pattern.
    # from: https://stackoverflow.com/a/27517681/7836730
    # Assumes all treestrings (and nothing else) ends in ");"
    def count_treestring_terminators(filename):
        f = open(filename, 'rb')
        bufsize = 1024*1024
        bufgen = takewhile(lambda x: x, (f.raw.read(bufsize) for _ in repeat(None)))
        return sum( buf.count(b');') for buf in bufgen if buf )

    def count_trees_by_parsing(filename, args):
        # Open treefile. Discard (i.e., silently pass by) the requested number of trees
        if args.informat == "nexus":
            treefile = treelib.Nexustreefile(filename)
        else:
            treefile = treelib.Newicktreefile(filename)
        treecount = 0
        for tree in treefile:
            treecount += 1
        return treecount

    count_list = []
    burnin_list = []
    sys.stdout.write("\n")
    for (wt, filename) in wt_file_list:
        treelist = []
        sys.stdout.write(f"   Counting trees in file {str(filename):<40}")
        sys.stdout.flush()
        n_tot = count_trees_by_parsing(filename, args)
        sys.stdout.write(f"{n_tot:>15,d}\n")
        sys.stdout.flush()
        burnin = int(args.burninfrac * n_tot)
        count_list.append(n_tot)
        burnin_list.append(burnin)

    # If automatic weighting requested: Compute new weights
    if args.autoweight:
        max_count = max(count_list)
        relwt_list = [max_count / count for count in count_list]
        relwtsum = sum(relwt_list)
        new_wt_list = [relwt / relwtsum for relwt in relwt_list]    # Normalized so sum=1

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
####################################################################################

def process_trees(wt_count_burnin_filename_list, args):

    treesummarylist = []
    for i, (weight, count, burnin, filename) in enumerate(wt_count_burnin_filename_list):
        sys.stdout.write("\n   Analyzing file: {} (Weight: {:5.3f})".format(filename, weight))
        sys.stdout.flush()

        # Open treefile. Discard (i.e., silently pass by) the requested number of trees
        if args.informat == "nexus":
            treefile = treelib.Nexustreefile(filename)
        else:
            treefile = treelib.Newicktreefile(filename)
        for j in range(burnin):
            treefile.readtree()
        sys.stdout.write(f"\n   Discarded {burnin:,} of {count:,} trees (burnin fraction={args.burninfrac:.2f})")

        # Instantiate Treesummary.
        # Re-use interner from first Treesummary to avoid duplication
        if i>0:
            interner = treesummarylist[0].interner
        else:
            interner = None
        if args.treeprobs or args.mbc:
            treesummary = treelib.BigTreeSummary(interner=interner,
                                                 store_trees=True)
        else:
            treesummary = treelib.TreeSummary(interner=interner)

        # Read remaining trees from file, add to treesummary
        sys.stdout.write("\n   Processing trees ('.' signifies 100 trees):\n")
        sys.stdout.flush()
        sys.stdout.write("\n   ")

        n_trees = 0
        for tree in treefile:
            n_trees += 1
            treesummary.add_tree(tree, weight)

            # Progress indicator
            if n_trees % 5000 == 0:
                sys.stdout.write(".  %7d" % n_trees)
                if args.verbose:
                    n_leaves = len(tree.leaves)
                    if args.treeprobs:
                        sys.stdout.write("   (# bip: %6d    # topo: %6d)\n   " % (len(treesummary.bipartsummary)-n_leaves, len(treesummary.toposummary)))
                    else:
                        sys.stdout.write("   (# bip: %6d)\n   " % (len(treesummary.bipartsummary)-n_leaves))
                else:
                    sys.stdout.write("\n   ")
                sys.stdout.flush()
            elif n_trees % 100 == 0:
                sys.stdout.write(".")
                sys.stdout.flush()

        treesummarylist.append(treesummary)
        print("\n")

    return treesummarylist

##########################################################################################
##########################################################################################

def compute_and_print_converge_stats(treesummarylist, minfreq):
    """Compute average bipartition frequency standard deviation between treesummaries"""

    # NOTES ON COMPUTATION:
    #   (1) all bipartitions are included in computation, regardless of whether they are present
    #       in only one of the treesummaries (their freq is set to zero in those treesummaries
    #       where they are not present)
    #   (2) external branches (leading to leafs) are not included in computation since they
    #       will allways be present in all trees (=> p=1, std=0)
    #   (3) N-1 is used in the denominator when computing std. This has to do with sample vs.
    #       population estimates of std, and is in accordance with what MrBayes does. One could
    #       argue that N should be used instead of N-1 (to get the Max Likelihood estimate of std).

    sum_std = 0
    N = float(len(treesummarylist))

    # Find combined set of bipartitions (excluding external branches)
    bipset = set()
    for treesummary in treesummarylist:
        for bipart in treesummary.bipartsummary:
            (bip1, bip2) = bipart
            if (len(bip1)==1 or len(bip2)==1):
                pass
            else:
                bipset.add(bipart)

    # Discard rare bipartitions: Only biparts that have freq >= minfreq are kept
    bipset_keep = set()
    for bipart in bipset:
        for treesummary in treesummarylist:
            try:
                if treesummary.bipartsummary[bipart].freq >= minfreq:
                    bipset_keep.add(bipart)
                    break
            except:
                pass
    del(bipset)

    # For each internal bipart: compute std of freq of this bipart across all treesummaries
    for bipart in bipset_keep:

        freqsum = 0
        sumsq = 0
        for treesummary in treesummarylist:
            # If current bipartition in treesummary, add freq etc. Otherwise add zero (=do nothing)
            if bipart in treesummary.bipartsummary:
                freqsum += treesummary.bipartsummary[bipart].freq
                sumsq += (treesummary.bipartsummary[bipart].freq)**2

        meansq = (freqsum/N)**2
        std = math.sqrt((sumsq-N*meansq)/(N-1))
        sum_std += std

    ave_std = sum_std / len(bipset_keep)

    print(("   Average standard deviation of split frequencies: {:.6f}\n".format(ave_std)))


##########################################################################################
##########################################################################################

def  merge_treesummaries(treesummarylist):

    treesummary = treesummarylist[0]
    for treesummary2 in treesummarylist[1:]:
        treesummary.update(treesummary2)
        del treesummary2
    return treesummary

##########################################################################################
##########################################################################################


def compute_and_print_biparts(treesummary, args):

    # Compute and retrieve results
    (leaflist, bipreslist) = bipart_report(treesummary, args)

    # Before printing results: check whether files already exist
    partsfilename = args.outbase.parent / (args.outbase.name + ".parts")
    if args.nowarn:
        partsfile = open(partsfilename, "w")
    elif partsfilename.is_file():
        overwrite = input(f"   File {partsfilename} already exists.\n   Overwrite (y/n): ")
        if overwrite== "y":
            partsfile = open(partsfilename, "w")            # Overwrite
            print(f"   Overwriting file {partsfilename}\n")
        else:
            partsfile = open(partsfilename, "a")            # Append
            print(f"   Appending to file {partsfilename}\n")
    else:
        partsfile = open(partsfilename, "w")

    # Print bipartitions
    partsfile.write("List of bipartitions:\n\n"
                    "PART = Description of partition in .* format\n"
                    "PROB = Posterior probability of the partition\n"
                    "BLEN = Mean branch length\n"
                    "VAR  = Branch length variance\n"
                    "SEM  = Standard error of the mean for branch length\n"
                    "ID   = Leaf name or internal branch label, for those bipartitions that are included in consensus tree\n\n")
    stringwidth = len(leaflist)
    partsfile.write("PART" + (stringwidth-1)*" " + "PROB      " + "BLEN       " + "VAR          " + "SEM          " + "ID\n")

    for (_, _, bipstring, freq, mean, var, sem, branchID) in bipreslist:
        if var == "NA":
            partsfile.write(f"{bipstring}   {freq:<8.6f}  {mean:<9.4g}  ({var:<9.4g})  ({sem:<9.4g})  {branchID}\n")
        else:
            partsfile.write(f"{bipstring}   {freq:<8.6f}  {mean:<9.4g}  ({var:<9.4g})  ({sem:<9.4g})  {branchID}\n")
    partsfile.close()
    print(f"   Bipartition list written to {partsfilename}")

##########################################################################################
##########################################################################################

def bipart_report(treesummary, args):
    """Return processed, almost directly printable, summary of all observed bipartitions"""

    leaflist = sorted(treesummary.leaves)
    position_dict = {}
    for position, leaf in enumerate(leaflist):
        position_dict[leaf] = position
    bipreport = []
    for _,bipart in treesummary.sorted_biplist:
        branch = treesummary.bipartsummary[bipart]
        freq = branch.freq
        if freq > args.minfreq:
            bipstring = bipart_to_string(bipart, position_dict, leaflist)
            bipsize = bipstring.count("*")              # Size of smaller set
            bipreport.append([1-freq, bipsize, bipstring,
                              freq, branch.length, branch.var, branch.sem, branch.branchID])

        bipreport = sorted(bipreport, key=itemgetter(0,1,2))
    # Return tuple of (leaflist, bipreport)
    return (leaflist, bipreport)

##########################################################################################

def bipart_to_string(bipartition, position_dict, leaflist):
    """Takes bipartition (set of two leaf sets) and returns string representation"""

    bipart1, bipart2 = bipartition

    # Bipartstring will be built as a list of chars first. Initialize with all "."
    stringwidth = len(leaflist)
    bipart_list = stringwidth * ["."]

    # Smaller set is represented by "*" characters. Larger set by "." characters (already set)
    if len(bipart1) < len(bipart2):
        smallset = bipart1
    else:
        smallset = bipart2

    for leaf in smallset:
        pos = position_dict[leaf]
        bipart_list[pos] = "*"

    return "".join(bipart_list)        # Concatenate into one string

##########################################################################################
##########################################################################################

def compute_and_print_contree(treesummary, args):

    if args.mbc:
        contree, logbipcred = treesummary.max_clade_cred_tree()
    else:
        contree = treesummary.contree(allcompat=args.all)
        logbipcred = treesummary.log_clade_credibility(contree.topology())

    if args.outgroup:
        contree.rootout(args.outgroup)
    elif args.midpoint:
        contree.rootmid()
    elif args.minvar:
        contree.rootminvar()
    newick_prob = contree.newick(labelfield="label")
    newick_branchID = contree.newick(labelfield="branchID")

    if args.mbc:
        confilename = args.outbase.parent / (args.outbase.name + ".mbc")
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
            confile = open(confilename, "a")            # Append
            print("   Appending to file {}\n".format(confilename))
    else:
        confile = open(confilename, "w")

    confile.write("#NEXUS\n")
    confile.write("\n")
    confile.write("begin trees;\n")
    confile.write("   [In this tree branch labels indicate the posterior probability of the bipartition corresponding to the branch.]\n")
    confile.write("   tree prob = ")
    confile.write(newick_prob)
    confile.write("\n\n   [In this tree branch labels indicate the bipartition ID listed in the file {}.\n".format(args.outbase.name + ".parts"))
    confile.write("    These branch labels can be used for interpreting the table of branch lenght info in that file]\n")
    confile.write("   tree partID = ")
    confile.write(newick_branchID)
    confile.write("\nend;\n")
    confile.close()

    if args.mbc:
        print(f"   Maximum bipartition credibility tree written to {confilename}")
    else:
        print("   Consensus tree written to {}".format(confilename))

    return contree, logbipcred

##########################################################################################
##########################################################################################

def compute_and_print_trprobs(treesummary, args):
    topolist = topo_report(treesummary)

    # Before printing results: check whether file already exist
    topofilename = args.outbase.parent / (args.outbase.name + ".trprobs")
    if args.nowarn:
        topofile = open(topofilename, "w")
    elif topofilename.is_file():
        overwrite = input(f"\n   File {topofilename} already exists.\n   Overwrite (y/n): ")
        if overwrite== "y":
            topofile = open(topofilename, "w")            # Overwrite
            print(f"   Overwriting file {topofilename}\n")
        else:
            topofile = open(topofilename, "a")            # Append
            print(f"   Appending to file {topofilename}\n")
    else:
        topofile = open(topofilename, "w")

    topofile.write("#NEXUS\n")
    topofile.write("\n")
    if args.treeprobs < 1:
        topofile.write(f"[This file contains the {round(args.treeprobs*100)}% most probable trees found during the\n")
        topofile.write(f"MCMC search, sorted by posterior probability (the {round(args.treeprobs*100)}% HPD interval).\n")
    else:
        topofile.write("[This file contains all trees that were found during the MCMC\n")
        topofile.write("search, sorted by posterior probability. \n")
    topofile.write("Lower case 'p' indicates the posterior probability of a tree.\n")
    topofile.write("Upper case 'P' indicates the cumulative posterior probability.]\n")
    topofile.write("\n")
    topofile.write("begin trees;\n")
    topofile.write(treesummary.translateblock)
    n=1
    cum = 0.0
    for (freq, tree) in topolist:
        cum += freq
        treestring = tree.newick(printdist=False, printlabels=False, transdict=treesummary.transdict)
        topofile.write(f"    tree tree_{n} [p = {freq:.6f}] [P = {cum:.6f}] = {treestring}\n")
        n += 1
        if cum > args.treeprobs:
            break

    topofile.write("end;\n")
    topofile.close()
    print(f"   Tree probabilities written to {topofilename}")

##########################################################################################
##########################################################################################

def topo_report(treesummary):
    """Returns list of [freq, treestring] lists"""

    # Python note: root trees in trprobs?
    toporeport = []
    for topology, topostruct in treesummary.toposummary.items():
        toporeport.append((topostruct.freq, topostruct.tree))

    # Sort report according to frequency (higher values first) and return
    toporeport = sorted(toporeport, key=itemgetter(0), reverse=True)

    return toporeport

##########################################################################################
##########################################################################################

if __name__ == "__main__":
    main()

