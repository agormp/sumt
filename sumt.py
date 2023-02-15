####################################################################################
####################################################################################

import phylotreelib as treelib
import os, sys, time, math, copy, psutil
from optparse import OptionParser
from itertools import (takewhile,repeat)
from operator import itemgetter
import gc

gc.disable()        # Faster. Assume no cyclic references will ever be created

def main():

    parser = build_parser()
    (options, wt_file_list) = parse_commandline(parser)
    start=time.time()

    try:
        # Silence output to stdout if requested
        # NOTE: automatically sets the "nowarn" option since interactivity is not possible
        if options.quiet:
            sys.stdout = open(os.devnull, 'w')
            options.nowarn = True

        if options.rootfile:
            outgroup = read_outgroup(options.rootfile)
        else:
            outgroup = options.outgroup

        (n_trees_analyzed, wt_count_burnin_filename_list) = count_trees(wt_file_list, options)

        if options.verbose:
            pid = psutil.Process(os.getpid())
            memory1 = pid.memory_full_info().rss

        treesummarylist = process_trees(wt_count_burnin_filename_list, options, outgroup)

        # Get basename of output files.
        # Use name supplied by user if present. Make sure all intermediate directories exist
        if options.outbase:
            filename = os.path.basename(options.outbase)     # NB: os.path.basename extracts filename from path
            dirname = os.path.dirname(options.outbase)
            if dirname and not os.path.isdir(dirname):
                os.makedirs(dirname)                        # construct any missing dirs from path
            outname = os.path.join(dirname, filename)

        # If basename is not set: determine semi-intelligently from infilenames:
        elif len(wt_file_list) > 1:
            for (wt, outname) in wt_file_list:
                if outname.endswith(".run1.t"):
                    outname = outname[:-7]
                    break
                elif outname.endswith(".t"):
                    outname = outname[:-2]
        else:
            outname = wt_file_list[0][1]
            if outname.endswith(".t"):
                outname = outname[:-2]

        if options.std:
            compute_and_print_converge_stats(treesummarylist, options.minfreq)

        # Collect all trees in single treesummary so final summary tree and bipart stats can be computed
        treesummary =  treesummarylist[0]
        for treesummary2 in treesummarylist[1:]:
            treesummary.update(treesummary2)
            del treesummary2
        n_leafs = len(treesummary.leaves)
        total_unique_internal_biparts = len(treesummary.bipartsummary) - n_leafs

        treesummary.add_branchid()

        contree = compute_and_print_contree(treesummarylist[0], options.allcomp, outgroup, outname,
                                                      options.midpoint, options.minvar, options.nowarn, options.outformat,
                                                      options.mbc)

        compute_and_print_biparts(treesummary, outname, options.nowarn, options.minfreq, options.mbc,
                                                      options.allcomp)

        theo_maxbip_internal = n_leafs - 2            # Maximum theoretical number of internal biparts = n-2 (rooted)
        n_internal_biparts = contree.n_bipartitions()

        if options.treeprobs:
            compute_and_print_trprobs(treesummarylist[0], options.treeprobs, outname, options.nowarn)
            n_topo_seen = len(treesummarylist[0].toposummary)
        stop=time.time()

        if options.verbose:
            memory2 = pid.memory_full_info().rss
            memorymax = max(memory1, memory2)

        time_spent=stop-start
        h = int(time_spent/3600)
        m = int((time_spent % 3600)/60)
        s = int(time_spent % 60)
        print("\n   Done. {:,d} trees analyzed.\n   Time spent: {:d}:{:02d}:{:02d} (h:m:s)\n".format(n_trees_analyzed, h, m, s))
        if options.verbose:
            if options.treeprobs:
                print("   Different topologies seen: {:8,d}".format(n_topo_seen))
                print("   Different bipartitions seen: {:6,d} (theoretical maximum: {:,d})".format(
                                                        total_unique_internal_biparts, theo_maxbip_internal * n_trees_analyzed))
            else:
                print("   Different bipartitions seen: {:6,d} (theoretical maximum: {:,d})".format(
                                                        total_unique_internal_biparts, theo_maxbip_internal * n_trees_analyzed))
            print("   Bipartitions in consensus tree: {:3,d} (theoretical maximum: {:,d})".format(
                                                        n_internal_biparts, theo_maxbip_internal))
            if memorymax > 1E9:
                print("   Max memory used: {:,.2f} GB.".format( memorymax  / (1024**3) ))
            else:
                print("   Max memory used: {:,.2f} MB.".format( memorymax  / (1024**2) ))
    except treelib.TreeError as exc:
        if options.verbose:
            import traceback
            traceback.print_exc(file=sys.stdout)
        else:
            print("Error: ", exc.errormessage)

        sys.exit()

####################################################################################
####################################################################################

def build_parser():
    use = """usage: %prog [options] FILE [FILE ...]\n       %prog [options] -w WEIGHT FILE -w WEIGHT FILE ..."""
    vers = "%prog 2.1.4"
    parser = OptionParser(usage=use, version=vers)
    parser.set_defaults(burninfrac=0.25, minfreq=0.1, allcomp=False, autoweight=False, outgroup=None,
                        rootfile=None, midpoint=False, informat="NEXUS", outformat="NEXUS", mbc=False,
                        nowarn=False, std=False, treeprobs=None, verbose=False, fileweights=None)

    parser.add_option("-I", type="choice", dest="informat",
                      choices=["NEXUS", "nexus", "NEWICK", "newick"], metavar="FORM",
                      help="format of input: nexus or newick [default: nexus]")

    parser.add_option("-O", type="choice", dest="outformat",
                      choices=["NEXUS", "nexus", "NEWICK", "newick"], metavar="FORM",
                      help="format of output: nexus or newick [default: nexus]")

    parser.add_option("--mbc", action="store_true", dest="mbc",
                      help="summarise trees using Maximum Bipartition Credibility (MBC) tree instead of majority rule consensus tree. "
                           + "MBC is similar to MCC (Maximum Clade Credibility) tree "
                           + "but counting bipartitions instead of clades, i.e. ignoring rooting. "
                           + "Additionally, branch lengths are estimated from branch lengths of bipartitions "
                           + "and not from node depths (i.e., again ignoring rooting)")

    parser.add_option("-q", action="store_true", dest="quiet",
                      help="quiet: don't print progress indication to terminal window. NOTE: also turns on the -n option")

    parser.add_option("-v", action="store_true", dest="verbose",
                      help="verbose: more information, longer error messages")

    parser.add_option("-n", action="store_true", dest="nowarn",
                      help="no warning when overwriting files")

    parser.add_option("--basename", action="store", type="string", dest="outbase", metavar="NAME",
                      help="base name of output files (default: derived from input file)")

    parser.add_option("-b", type="float", dest="burninfrac", metavar="NUM",
                      help="burnin: fraction of trees to discard [0 - 1; default: 0.25]")

    parser.add_option("-t", type="float", dest="treeprobs", metavar="NUM",
                      help="compute tree probabilities, report NUM percent credible interval [0 - 1]")

    parser.add_option("-s", action="store_true", dest="std",
                      help="compute average standard deviation of split frequencies (ASDSF)")

    parser.add_option("-f", type="float", dest="minfreq", metavar="NUM",
                      help="Min. frequency for including bipartitions in report and in computation of ASDSF [default: 0.1]")

    parser.add_option("-a", action="store_true", dest="allcomp",
                      help="add all compatible bipartitions to consensus tree")

    parser.add_option("-w", action="append",
                      type="string", dest="fileweights", nargs=2, metavar="WEIGHT FILE -w WEIGHT FILE ...",
                      help="put different weights on different FILEs")

    parser.add_option("--autow", action="store_true", dest="autoweight",
                    help="automatically assign file weights based on tree counts, so all files have equal impact")

    parser.add_option("-m", "--rootmid", action="store_true", dest="midpoint",
                      help="perform midpoint rooting of tree")

    parser.add_option("--rootminvar", action="store_true", dest="minvar",
                      help="perform minimum variance rooting of tree")

    parser.add_option("-r", action="append",
                      type="string", dest="outgroup", metavar="TAX [-r TAX ...]",
                      help="root consensus tree on specified outgroup taxon/taxa")

    parser.add_option("--rootfile", action="store", type="string", dest="rootfile", metavar="FILE",
                      help="root consensus tree on outgroup taxa listed in file (one name per line)")


    return parser

####################################################################################
####################################################################################

def parse_commandline(parser):

    # Parse commandline, check for errors
    (options, args) = parser.parse_args()

    # If neither unweighted nor weighted filenames are listed: abort with error message
    if len(args) == 0 and not options.fileweights:
        parser.error("Please list one or more tree files.")

    # If both weighted and unweighted filenames are listed: abort with error message
    elif len(args) > 0 and options.fileweights:
        print(("args: {}    op.fileweights: {}".format(args, options.fileweights)))
        parser.error("Weights must be given for all listed tree files")

    # If only unweighted filenames are given:
    # Reformat list of filenames into (weight, filename) tuple format expected by program
    # Set all weights to 1/n_files
    elif len(args) > 0 and not options.fileweights:
        wt_file_list = []
        wt = 1.0 / len(args)
        for filename in args:
            wt_file_list.append((wt, filename))

    # If only weighted filenames are listed:
    # reformat list of tuples such that weight is in float (not string). Normalize weights so they sum to one.
    elif len(args) == 0 and options.fileweights:

        tmp_wt_file_list = []

        # Attempt to convert weight string to float. Print sensible error message if this fails
        for (wt_string, filename) in options.fileweights:
            try:
                wt = float(wt_string)
            except ValueError as exc:
                badfloat = exc[0].split()[-1]     # Last word in "exc" is the offending value
                print(('Invalid file weight: "{}" - value has to be a real number.'.format(badfloat)))
                sys.exit()

            tmp_wt_file_list.append((wt, filename))

        # Normalize weights, build final weight/file list:
        wtsum = 0.0
        for (wt, filename) in tmp_wt_file_list:
            wtsum += wt
        wt_file_list = []
        for (wt, filename) in tmp_wt_file_list:
            wt_file_list.append((wt/wtsum, filename))


    for (wt, filename) in wt_file_list:
        if not os.path.isfile(filename):
            parser.error("File %s not found." % filename)

    if options.burninfrac > 1 or options.burninfrac < 0:
        parser.error("option -b: NUM must be between 0.0 and 1.0")

    if options.treeprobs and (options.treeprobs > 1 or options.treeprobs < 0):
        parser.error("option -t: NUM must be between 0.0 and 1.0")

    if options.outgroup and options.rootfile:
        parser.error("options -r and --rootfile are incompatible")

    if options.std and len(wt_file_list)==1:
        parser.error("cannot compute standard deviation from one tree file")

    n_root_options = sum(1 for opt in [options.midpoint, options.minvar, options.outgroup, options.rootfile] if opt)
    if n_root_options > 1:
        parser.error("More than one option for rooting was specified. Please use only one.")

    if options.rootfile:
        if not os.path.isfile(options.rootfile):
            parser.error("File %s not found." % options.rootfile)

    # return options, weights, and filenames for external use:
    return (options, wt_file_list)

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

def count_trees(wt_file_list, options):

    # fast counting of specific pattern.
    # from: https://stackoverflow.com/a/27517681/7836730
    # Assumes all treestrings (and nothing else) ends in ");"
    def count_treestring_terminators(filename):
        f = open(filename, 'rb')
        bufsize = 1024*1024
        bufgen = takewhile(lambda x: x, (f.raw.read(bufsize) for _ in repeat(None)))
        return sum( buf.count(b');') for buf in bufgen if buf )

    def count_trees_by_parsing(filename, options):
        # Open treefile. Discard (i.e., silently pass by) the requested number of trees
        if options.informat.lower() == "nexus":
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
        sys.stdout.write("   Counting trees in file {:<40}".format("'" + filename + "'" ":"))
        sys.stdout.flush()
        n_tot = count_trees_by_parsing(filename, options)
        sys.stdout.write("{:>15,d}\n".format(n_tot))
        sys.stdout.flush()
        burnin = int(options.burninfrac * n_tot)
        count_list.append(n_tot)
        burnin_list.append(burnin)

    # If automatic weighting requested: Compute new weights
    if options.autoweight:
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
        if options.autoweight:
            wt = new_wt_list[i]
        else:
            wt = wt_file_list[i][0]

        wt_count_burnin_filename_list.append((wt, count, burnin, filename))

    n_trees_analyzed = sum(count_list) - sum (burnin_list)
    return (n_trees_analyzed, wt_count_burnin_filename_list)

####################################################################################
####################################################################################

def process_trees(wt_count_burnin_filename_list, options, outgroup):

    treesummarylist = []
    for i, (weight, count, burnin, filename) in enumerate(wt_count_burnin_filename_list):
        sys.stdout.write("\n   Analyzing file: {} (Weight: {:5.3f})".format(filename, weight))
        sys.stdout.flush()

        # Open treefile. Discard (i.e., silently pass by) the requested number of trees
        if options.informat.lower() == "nexus":
            treefile = treelib.Nexustreefile(filename)
        else:
            treefile = treelib.Newicktreefile(filename)
        for j in range(burnin):
            treefile.readtree()
        sys.stdout.write(f"\n   Discarded {burnin:,} of {count:,} trees (burnin fraction={options.burninfrac:.2f})")

        # Instantiate Treesummary.
        # Re-use interner from first Treesummary to avoid duplication
        if i>0:
            interner = treesummarylist[0].interner
        else:
            interner = None
        if options.treeprobs or options.mbc:
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
                if options.verbose:
                    n_leaves = len(tree.leaves)
                    if options.treeprobs:
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

def compute_and_print_biparts(treesummary, filename, nowarn, minf, mbc, allcomp):

    # Compute and retrieve results
    (leaflist, bipreslist) = bipart_report(treesummary, minf, mbc, allcomp)

    # Before printing results: check whether files already exist
    partsfilename = filename + ".parts"
    if nowarn:
        partsfile = open(partsfilename, "w")
    elif os.path.isfile(partsfilename):
        overwrite = input("   File %s already exists.\n   Overwrite (y/n): " % partsfilename)
        if overwrite== "y":
            partsfile = open(partsfilename, "w")            # Overwrite
            print("   Overwriting file {}\n".format(partsfilename))
        else:
            partsfile = open(partsfilename, "a")            # Append
            print("   Appending to file {}\n".format(partsfilename))
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
    partsfile.write("PART" + (stringwidth-1)*" " + "PROB      " + "BLEN      " + "VAR         " + "SEM         " + "ID\n")

    for (_, _, bipstring, freq, mean, var, sem, branchID) in bipreslist:
        if var == "NA":
            partsfile.write(f"{bipstring}   {freq:8.6f}  {mean:8.6f}  ({var:>8})  ({sem:>8})  {branchID}\n")
        else:
            partsfile.write(f"{bipstring}   {freq:8.6f}  {mean:8.6f}  ({var:8.6f})  ({sem:8.6f})  {branchID}\n")

    print("   Bipartition list written to {}".format(partsfilename))

##########################################################################################
##########################################################################################

def bipart_report(treesummary, minfreq, mbc, allcomp):
    """Return processed, almost directly printable, summary of all observed bipartitions"""

    leaflist = sorted(treesummary.leaves)
    position_dict = {}
    for position, leaf in enumerate(leaflist):
        position_dict[leaf] = position
    bipreport = []
    for _,bipart in treesummary.sorted_biplist:
        branch = treesummary.bipartsummary[bipart]
        freq = branch.freq
        if freq > minfreq:
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

def compute_and_print_contree(treesummary, allcomp, outgroup, filename,
                              midpoint, minvar, nowarn, outformat, mbc):

    # Construct consensus tree with bipart freq labels
    if mbc:
        contree, logbipcred = treesummary.max_clade_cred_tree()
    else:
        contree = treesummary.contree(allcompat=allcomp)

    # If outgroup is given: attempt to root tree on provided outgroup.
    # If this is impossible then print warning and save midpoint rooted contree instead
    if outgroup:
        try:
            contree.rootout(outgroup)
        except treelib.TreeError as exc:
            print("Warning: ", exc.errormessage)
            print("Midpoint rooting used instead")
            contree.rootmid()
    elif midpoint:
        contree.rootmid()
    elif minvar:
        contree.rootminvar()
    newick_prob = contree.newick(labelfield="label")
    newick_branchID = contree.newick(labelfield="branchID")

    # Before printing results: check whether files already exist
    if mbc:
        confilename = filename + ".mbc"
    else:
        confilename = filename + ".con"
    if nowarn:
        confile = open(confilename, "w")
    elif os.path.isfile(confilename):
        overwrite = input("\n   File %s already exists.\n   Overwrite (y/n): " % confilename)
        if overwrite== "y":
            confile = open(confilename, "w")            # Overwrite
            print("   Overwriting file {}\n".format(confilename))
        else:
            confile = open(confilename, "a")            # Append
            print("   Appending to file {}\n".format(confilename))
    else:
        confile = open(confilename, "w")

    if outformat=="NEWICK" or outformat == "newick":
        confile.write(newick)
        confile.write("\n")
    else:
        confile.write("#NEXUS\n")
        confile.write("\n")
        confile.write("begin trees;\n")
        confile.write("   [In this tree branch labels indicate the posterior probability of the bipartition corresponding to the branch.]\n")
        confile.write("   tree prob = ")
        confile.write(newick_prob)
        confile.write("\n\n   [In this tree branch labels indicate the bipartition ID listed in the file {}.\n".format(filename + ".parts"))
        confile.write("    These branch labels can be used for interpreting the table of branch lenght info in that file]\n")
        confile.write("   tree partID = ")
        confile.write(newick_branchID)
        confile.write("\nend;\n")

    if mbc:
        print(f"   Maximum bipartition credibility tree written to {confilename}")
        print(f"   Highest Log Bipartition Credibility:  {logbipcred:.4g}")
    else:
        print("   Consensus tree written to {}".format(confilename))

    return contree

##########################################################################################
##########################################################################################

def compute_and_print_trprobs(treesummary, hpd_frac, filename, nowarn):
    topolist = topo_report(treesummary)

    # Before printing results: check whether file already exist
    topofilename = filename + ".trprobs"
    if nowarn:
        topofile = open(topofilename, "w")
    elif os.path.isfile(topofilename):
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
    if hpd_frac < 1:
        topofile.write(f"[This file contains the {round(hpd_frac*100)}% most probable trees found during the\n")
        topofile.write(f"MCMC search, sorted by posterior probability (the {round(hpd_frac*100)}% HPD interval).\n")
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
        treestring = tree.newick(printdist=False, printlabels=False)
        topofile.write(f"   tree tree_{n} [p = {freq:.6f}] [P = {cum:.6f}] = {treestring}\n")
        n += 1
        if cum > hpd_frac:
            break

    topofile.write("end;\n")
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

