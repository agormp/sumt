####################################################################################
####################################################################################

import phylotreelib as treelib
import os, sys, time, math, subprocess, copy
from optparse import OptionParser

# NOTE: maybe include this as option for aggressive optimization at runtime
# import gc
# gc.disable()        # Assume no cyclic references will ever be created

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

        (n_trees_analyzed, wt_count_burnin_file_list) = count_trees(wt_file_list, options)
        treesummarylist = process_trees(wt_count_burnin_file_list, options, outgroup)

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
            for treesummary in treesummarylist:
                treesummary.compute_bipfreq()
            compute_and_print_converge_stats(treesummarylist, options.minfreq)

        if options.verbose:
            # Check memory usage (this is presumably the point where most memory is used)
            pid= str(os.getpid())
            ps=subprocess.getoutput("ps up" + pid)
            memory=ps.split()[16]

        # Collect all trees in single treesummary so final contree and bipart stats can be computed
        for treesummary2 in treesummarylist[1:]:
            treesummarylist[0].update(treesummary2)
            del treesummary2
        total_unique_biparts = len(treesummarylist[0].bipartsummary)

        compute_and_print_biparts(treesummarylist[0], outname, options.nowarn, options.minfreq)
        n_biparts_contree = compute_and_print_contree(treesummarylist[0], options.allcomp,
                                                                        outgroup, outname, options.midpoint, options.nowarn,
                                                                         options.outformat)
        n_leafs = len(treesummarylist[0].leaves)
        theo_maxbip = 2 * n_leafs - 3     # Theoretical maximum number of bipartitions in tree = 2n-3
        n_internal_biparts = n_biparts_contree - n_leafs        # Number of internal bipartitions (excluding leaves)
        theo_maxbip_internal = theo_maxbip - n_leafs            # Maximum theoretical number of internal biparts = n-3

        if options.treeprobs:
            compute_and_print_trprobs(treesummarylist[0], options.treeprobs, outname, options.nowarn)
            n_trees_seen = len(treesummarylist[0].toposummary)
        stop=time.time()
        time_spent=stop-start
        h = int(time_spent/3600)
        m = int((time_spent % 3600)/60)
        s = int(time_spent % 60)
        print("\n   Done. {:,d} trees analyzed.\n   Time spent: {:d}:{:02d}:{:02d} (h:m:s)".format(n_trees_analyzed, h, m, s))
        if options.verbose:
            print("\n")
            if options.treeprobs:
                print("   Different topologies seen: {:8,d}".format(n_trees_seen))
                print("   Different bipartitions seen: {:6,d} (theoretical maximum: {:,d})".format(total_unique_biparts, theo_maxbip * n_trees_seen))
            else:
                print("   Different bipartitions seen: {:6,d}".format(total_unique_biparts))
            print("   Internal bipartitions in consensus tree: {:3,d} (theoretical maximum: {:,d})".format(n_internal_biparts, theo_maxbip_internal))
            print("   Memory used: {:,.2f} MB.".format(float(memory)/1024))
            print("")

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
    parser.set_defaults(burninfrac=0.25, minfreq=0.1, allcomp=False, autoweight=False, zeroterms=False, outgroup=None,
                        rootfile=None, midpoint=False, informat="NEXUS", outformat="NEXUS",
                        nowarn=False, std=False, treeprobs=None, verbose=False, fileweights=None)

    parser.add_option("-I", type="choice", dest="informat",
                      choices=["NEXUS", "nexus", "NEWICK", "newick"], metavar="FORM",
                      help="format of input: nexus or newick [default: nexus]")

    parser.add_option("-O", type="choice", dest="outformat",
                      choices=["NEXUS", "nexus", "NEWICK", "newick"], metavar="FORM",
                      help="format of output: nexus or newick [default: nexus]")

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
                      help="compute tree probabilities, report NUM percent credible interval [0 - 1; default: 1.00]")

    parser.add_option("-s", action="store_true", dest="std",
                      help="compute average standard deviation of split frequencies (ASDSF)")

    parser.add_option("-f", type="float", dest="minfreq", metavar="NUM",
                      help="Min. frequency for including bipartitions in report and in computation of ASDSF [default: 0.1]")

    parser.add_option("-a", action="store_true", dest="allcomp",
                      help="add all compatible bipartitions to consensus tree")

    parser.add_option("-z", action="store_true", dest="zeroterms",
                      help="include zero length terms when computing branch length and average standard deviation of split frequencies")

    parser.add_option("-w", action="append",
                      type="string", dest="fileweights", nargs=2, metavar="WEIGHT FILE -w WEIGHT FILE ...",
                      help="put different weights on different FILEs")



    parser.add_option("--autow", action="store_true", dest="autoweight",
                    help="automatically assign file weights based on tree counts, so all files have equal impact")

    parser.add_option("-m", action="store_true", dest="midpoint",
                      help="perform midpoint rooting of tree")



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

    if options.treeprobs > 1 or options.treeprobs < 0:
        parser.error("option -t: NUM must be between 0.0 and 1.0")

    if options.outgroup and options.rootfile:
        parser.error("options -r and --rootfile are incompatible")

    if options.std and len(wt_file_list)==1:
        parser.error("cannot compute standard deviation from one tree file")

    if options.midpoint and (options.outgroup or options.rootfile):
        parser.error("cannot perform midpoint (-m) and outgroup (-r/--rootfile) rooting simultaneously")

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

    # Count trees in all files
    count_list = []
    sys.stdout.write("\n")
    for (wt, filename) in wt_file_list:
        sys.stdout.write("   Counting trees in file {:<40}".format("'" + filename + "'" ":"))
        sys.stdout.flush()
        if options.informat.lower() == "nexus":
            treefile = treelib.Nexustreefile(filename)
        else:
            treefile = treelib.Newicktreefile(filename)
        n_tot = 0
        for tree in treefile:
            n_tot += 1
        sys.stdout.write("{:>15,d}\n".format(n_tot))
        sys.stdout.flush()
        count_list.append(n_tot)

    # Compute burnin numbers
    burnin_list = [int(options.burninfrac * count) for count in count_list]

    # If automatic weighting requested: Compute new weights
    if options.autoweight:
        max_count = max(count_list)
        relwt_list = [max_count / count for count in count_list]
        relwtsum = sum(relwt_list)
        new_wt_list = [relwt / relwtsum for relwt in relwt_list]    # Normalized so sum=1

    # Construct final combined wt + count + burnin + filename list
    wt_count_burnin_file_list = []
    for i in range(len(wt_file_list)):
        filename = wt_file_list[i][1]
        count = count_list[i]
        burnin = burnin_list[i]
        if options.autoweight:
            wt = new_wt_list[i]
        else:
            wt = wt_file_list[i][0]
        wt_count_burnin_file_list.append((wt, count, burnin, filename))

    n_trees_analyzed = sum(count_list) - sum (burnin_list)
    return (n_trees_analyzed, wt_count_burnin_file_list)

####################################################################################
####################################################################################

def process_trees(wt_count_burnin_file_list, options, outgroup):

    treesummarylist = []
    for (weight, count, burnin, filename) in wt_count_burnin_file_list:
        sys.stdout.write("\n   Analyzing file: {} (Weight: {:5.3f})".format(filename, weight))
        sys.stdout.flush()
        # Open treefile. Discard (i.e., silently pass by) the requested number of trees
        if options.informat.lower() == "nexus":
            treefile = treelib.Nexustreefile(filename)
        else:
            treefile = treelib.Newicktreefile(filename)
        i = 0
        while i < burnin:
            i += 1
            next(treefile)
        sys.stdout.write("\n   Discarded {:,} of {:,} trees (burnin fraction={:.2f})".format(burnin, count, options.burninfrac))
        sys.stdout.write("\n   Processing trees ('.' signifies 100 trees):\n")
        sys.stdout.flush()

        if options.treeprobs:
            treesummary = treelib.BigTreeSummary()
        else:
            treesummary = treelib.TreeSummary()

        # Read remaining trees from file, add to treesummary
        n_trees = 0
        sys.stdout.write("\n   ")
        for tree in treefile:
            n_trees += 1
            treesummary.add_tree(tree, weight)

            # Progress indicator
            if n_trees % 5000 == 0:
                sys.stdout.write(".  %7d" % n_trees)
                if options.verbose:
                    if options.treeprobs:
                        sys.stdout.write("   (# bip: %6d    # topo: %6d)\n   " % (len(treesummary.bipartsummary), len(treesummary.toposummary)))
                    else:
                        sys.stdout.write("   (# bip: %6d)\n   " % len(treesummary.bipartsummary))
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

            # Exclude if bipart is external (i.e., when one part of it has only one member)
            (bip1, bip2) = bipart
            if (len(bip1)==1 or len(bip2)==1):
                pass
            else:
                bipset.add(bipart)

    # Discard rare bipartitions: Only biparts that have freq >= minfreq are kept
    bipsetcopy = copy.deepcopy(bipset)
    for bipart in bipsetcopy:
        discard = True
        for treesummary in treesummarylist:
            bipsum = treesummary.bipartsummary
            if bipart in bipsum and bipsum[bipart].freq >= minfreq:
                discard = False
        if discard:
            bipset.remove(bipart)
    del(bipsetcopy)

    # For each internal bipart: compute std of freq of this bipart across all treesummaries
    for bipart in bipset:

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

    ave_std = sum_std / len(bipset)

    print(("   Average standard deviation of split frequencies: {:.6f}\n".format(ave_std)))


##########################################################################################
##########################################################################################

def compute_and_print_biparts(treesummary, filename, nowarn, minf):

    # Compute bipart freq + branch length var and sem for combined tree summary
    treesummary.compute_bipfreq()
    treesummary.compute_blen_var_and_sem()

    # Compute and retrieve results
    (leaflist, biplist) = bipart_report(treesummary, minfreq=minf)

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

    # Print list of taxa
    partsfile.write("List of taxa in bipartitions:\n\n")
    i = 1
    for leaf in leaflist:
        partsfile.write("%2d -- %s\n" % (i, leaf))
        i += 1

    # Print bipartitions
    partsfile.write("\n\nList of bipartitions:\n\n"
                    "PART = Description of partition in .* format\n"
                    "PROB = Posterior probability of the partition\n"
                    "BLEN = Mean branch length\n"
                    "VAR  = Branch length variance\n"
                    "SEM  = Standard error of the mean for branch length\n\n")
    stringwidth = len(leaflist)
    partsfile.write("PART" + (stringwidth-1)*" " + "PROB      " + "BLEN      " + "VAR         " + "SEM\n")

    for bipart in biplist:
        partsfile.write("%s   %8.6f  %8.6f  (%8.6f)  (%8.6f)\n" % (bipart[1], bipart[0], bipart[2], bipart[3], bipart[4]))

    print("   Bipartition list written to {}".format(partsfilename))

##########################################################################################
##########################################################################################

def bipart_report(treesummary, minfreq=0.05):
    """Return processed, almost directly printable, summary of all observed bipartitions"""

    # Bipart report consists of a tuple containing:
    #       (0) a sorted list of leaves (for interpreting bipartstring)
    #       (1) a sorted list of lists. Each item list is: [bipartstring, freq, mean, var, sem]
    #           entire list is sorted on bipartition frequency

    # Must first figure out which leaves correspond to which positions in bipartstring
    # Note: leaves are ordered alphabetically, meaning first char in bipstring corresponds
    # to first leaf in alphabetic sort
    leaflist = sorted(treesummary.leaves)

    position_dict = {}
    for position, leaf in enumerate(leaflist):
        position_dict[leaf] = position

    # Loop over all bipartitions in raw_result, build formatted result list in process
    bipreport = []
    for bipart in treesummary.bipartsummary:
        bipstring = bipart_to_string(bipart, position_dict, leaflist)
        bipsize = bipstring.count("*")              # Size of smaller set

        # Only report bipartitions that occur more often than "minfreq":
        if treesummary.bipartsummary[bipart].freq>minfreq:
            freq = treesummary.bipartsummary[bipart].freq
            mean = treesummary.bipartsummary[bipart].mean
            var = treesummary.bipartsummary[bipart].var
            sem = treesummary.bipartsummary[bipart].sem
            bipreport.append([freq, bipstring, mean, var, sem])

    # Sort bipreport according to (1) frequency (higher values first), (2) size of
    # smaller bipartition (external branches before internal branches), and
    # (3) bipartstring (*.. before .*. before ..*)
    # First construct temporary list of (1-freq, bipsize, bipstring, originial list-item)
    # tuples. Sort this list of tuples and then re-extract the original list-item again
    # (Example of Decorate, Sort, Undecorate idiom)
    tmplist = sorted([(1-bip[0], bip[1].count("*"), bip[1], bip) for bip in bipreport])
    bipreport = [tup[-1] for tup in tmplist]        # Last element of tuple is orig list

    # Return tuple of (leaflist, bipreport)
    return (leaflist, bipreport)

##########################################################################################
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
                              midpoint, nowarn, outformat):

    # Construct consensus tree with bipart freq labels
    contree = treesummary.contree(allcompat=allcomp)
    n_biparts = len(contree.bipdict())

    # If outgroup is given: attempt to root tree on provided outgroup.
    # If this is impossible then print warning and save midpoint rooted contree instead
    if outgroup:
        try:
            contree.rootout(outgroup)
        except treelib.TreeError as exc:
            print("Warning: ", exc.errormessage)
            print("Midpoint rooting used instead")
            contree.rootmid()

    # Perform midpoint rooting if requested
    elif midpoint:
        contree.rootmid()
    newick = contree.newick()

    # Before printing results: check whether files already exist
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
        confile.write("   [In this tree branch labels indicate the posterior\n")
        confile.write("    probability of the bipartition corresponding to the branch.]\n")
        confile.write("   tree prob = ")
        confile.write(newick)
        confile.write("\nend;\n")

    print("   Consensus tree written to {}".format(confilename))

    return n_biparts

##########################################################################################
##########################################################################################

def compute_and_print_trprobs(treesummary, hpd_frac, filename, nowarn):
    topolist = topo_report(treesummary)

    # Before printing results: check whether file already exist
    topofilename = filename + ".trprobs"
    if nowarn:
        topofile = open(topofilename, "w")
    elif os.path.isfile(topofilename):
        overwrite = input("\n   File %s already exists.\n   Overwrite (y/n): " % topofilename)
        if overwrite== "y":
            topofile = open(topofilename, "w")            # Overwrite
            print("   Overwriting file {}\n".format(topofilename))
        else:
            topofile = open(topofilename, "a")            # Append
            print("   Appending to file {}\n".format(topofilename))
    else:
        topofile = open(topofilename, "w")

    topofile.write("#NEXUS\n")
    topofile.write("\n")
    if hpd_frac < 1:
        topofile.write("[This file contains the {}% most probable trees found during the\n".format(round(hpd_frac*100)))
        topofile.write("MCMC search, sorted by posterior probability (the {}% HPD interval).\n".format(round(hpd_frac*100)))
    else:
        topofile.write("[This file contains all trees that were found during the MCMC\n")
        topofile.write("search, sorted by posterior probability. \n")
    topofile.write("Lower case 'p' indicates the posterior probability of a tree.\n")
    topofile.write("Upper case 'P' indicates the cumulative posterior probability.]\n")
    topofile.write("\n")
    topofile.write("begin trees;\n")

    n=1
    cum = 0.0
    for (freq, treestring) in topolist:
        cum += freq
        topofile.write("   tree tree_%-4d [p = %.6f] [P = %.6f] = %s\n" % (n, freq, cum, treestring))
        n += 1
        if cum > hpd_frac:
            break

    topofile.write("end;\n")
    print("   Tree probabilities written to {}".format(topofilename))

##########################################################################################
##########################################################################################

def topo_report(treesummary):
    """Returns list of [freq, treestring] lists"""

    # Python note: root trees in trprobs?
    treesummary.compute_topofreq()
    toporeport = []
    for topology, topostruct in treesummary.toposummary.items():
        treestring = topostruct.treestring
        freq = topostruct.freq
        toporeport.append((freq, treestring))

    # Sort report according to frequency (higher values first) and return
    toporeport.sort(reverse=True)

    return toporeport

##########################################################################################
##########################################################################################

if __name__ == "__main__":
    main()
    # import cProfile
    # cProfile.run('main()', 'tmp/profile.pstats')


