#!/usr/bin/python

import sys, getopt

from random import random

verbose = False
verbose_debug = False
suppress_warnings = False
output = sys.stdout

import base_types
from base_types import *

def generateSequences(number=1, means = [0.1, 0.3, 0.5, 0.7, 0.9], meanBlockLength=1000, sdBlockLength=200, regLength=100000):
    from random import gauss as normal

    for i in xrange(1,number+1):
        totalLength = 0;
        sequences = []
        for mean in means:
            loop = 0;

            while loop < regLength:
                length = max(0, int(normal( meanBlockLength, sdBlockLength )))
                if random() <= mean:
                    sequences.append([totalLength+loop, totalLength+loop+length])
                loop += ( length + 2 )
            sequences[-1][-1] = min(sequences[-1][-1], totalLength+regLength-1)

            totalLength += regLength

        of = open('sim_output_%i.bed' % i, 'w')
        for entry in sequences:
            of.write( "name\t%i\t%i\n" % tuple(entry) )
        of.close()

    of = open('sim_lengths.txt', 'w')
    of.write("name\t0\t%i" % totalLength )
    of.close()

def change_point_mle(data, lb, ub, si, bl):
    """Process 4.7 from the paper. 

    Data should be a cumm sum object.
    """
    # if the interval isn't long enough, return 0
    if ( ub - lb ) <= 2*bl: return 0

    # n is the total considered interval length
    # we add 1 because the interval is inclusive
    n = float(ub - lb + 1)
    # we add 1 to j because the interval includes data[si]
    j = float(si - lb + 1)

    s1_mean = float(data.value(si) - data.value(lb))/j
    s2_mean = float(data.value(ub) - data.value(si))/float(ub-si)
    tot_mean = float(data.value(ub) - data.value(lb))/n

    p1 = (j/n)*(s1_mean - tot_mean)**2
    p2 = ((n-j)/n)*(s2_mean - tot_mean)**2
    
    if verbose_debug:
        print "%i\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f" % \
            ( lb, ub, si, s1_mean, s2_mean, tot_mean, p1+p2 )
    return p1+p2

def split_var_estimate( data, lb, ub, si, bl ):
    """Process 4.8 from the paper. 

    Data should be a cumm sum object.
    """
    # n is the total considered interval length
    # we add 1 because the interval is inclusive
    n = float(ub - lb + 1)
    # we add 1 to t because the interval includes data[si]
    t = float(si - lb + 1)

    # the two sum's coefficients
    p1_coef = t/(n**2)
    p2_coef = (n-t)/(n**2)
    
    # the two sections means
    s1_mean = float(data.value(si) - data.value(lb))/t
    s2_mean = float(data.value(ub) - data.value(si+1))/(ub-si+1)
    
    def wmn(i, offset):
        return float(data.value(lb+i+offset) - data.value(lb+i))/offset

    # part 1 of equation 4.8
    offset = int(t*bl/n)
    sti, ei = ( 0, int(t-t*bl/n) )
    p1 = p1_coef*sum( [(wmn(i, offset) - s1_mean)**2 for i in xrange(sti, ei)] ) 

    # part 2 of equation 4.8
    sti, ei = ( int(t+1), int(n-(n-t)*bl/n) )
    offset = int((n-t)*bl/n)
    p2 = p2_coef*sum( [(wmn(i, offset) - s2_mean)**2 for i in xrange(sti, ei)] ) 
    
    return p1+p2

def split_p_value( data, lb, ub, si, bl ):
    """Probability that the means... ( What is this )

    Return: z_score, p_value ( z_score compared to chisq(df=1) )
    """
    mean_estimate = change_point_mle( data, lb, ub, si, bl )
    var_estimate = split_var_estimate( data, lb, ub, si, bl )
    test_stat = mean_estimate/sqrt( var_estimate )
    p_value = r.pchisq(test_stat, df=1)[0]
    if p_value > 0.5: p_value = 1.0-p_value
    return test_stat, p_value

def segment(data, minLength, b=0, maxIterations=100, M_TOL = 0):
    """A method to segment an annotation track.

    """

    # an object to keep track of the scores we've already calculated
    cached_scores = {}

    def find_best_window_score(data, lb, ub):
        if cached_scores.has_key((lb, ub)):
            return cached_scores[(lb, ub)]
        else:
            # the best MLE score ( in splitWindow )
            currBestScore = 0
            # the best split index ( in splitWindow )
            currBestIndex = None

            # loop through each possible index and calculate the scores
            # note that we're calling the iter_split_points method:
            #   this is an optimization for binary features - it only looks
            #   at the potential split points 
            if verbose_debug:
                print "lb\tub\tsi\tl_mean\tu_mean\tt_mean\tscore"
            for split_index in data.iter_split_points(lb+minLength, ub-minLength):
                # calculate the score of using this change point
                score = change_point_mle(data, lb, ub, split_index, minLength)

                # eliminate the occasional rounding error
                if score < M_TOL: score = 0

                # if this is the current best score, store a record of it
                if score > currBestScore:
                    currBestScore = score
                    currBestIndex = split_index

            cached_scores[(lb, ub)] = ( currBestScore, currBestIndex )
            return ( currBestScore, currBestIndex )

    def find_split(data, splitIndexes):
        bestWindow = None
        bestIndex = None
        bestScore = 0

        for splitWindowIndex in xrange(len(splitIndexes) - 1):
            # define the window that we will be looking for a split in
            splitWindow = (splitIndexes[splitWindowIndex], splitIndexes[splitWindowIndex+1])

            # if the window isnt at least 2L, it will be impossible to find a split of length L
            if ( splitWindow[1] - splitWindow[0] ) < 2*minLength: continue

            # find the best score and it's index
            score, index = find_best_window_score(data, splitWindow[0], splitWindow[1])
            
            # if there is no possible split, continue
            if index == None:
                continue
            
            if b > 0:
                # calculate J
                V = split_var_estimate(data, splitWindow[0], splitWindow[1], index, minLength) 
                B = score
                reg_len = float(splitWindow[1] - splitWindow[0])
                lam = minLength*reg_len/( splitIndexes[-1] - splitIndexes[0] )

                J = int( (reg_len*B)/sqrt(V*lam) > b )
                if verbose:
                    print >> output, "b: ", sqrt(reg_len*B)/(V*lam), " > ", b, "\tM: ", score, "\tBest Index: ", index

            # if the best score for this window is the best so far for all the windows
            # then make a record of that
            if score > bestScore and ( b == 0 or J > 0 ) :
                bestIndex = index
                bestScore = score
                bestWindow = splitWindowIndex

        if bestIndex != None:
            best_V_score = split_var_estimate(data, splitIndexes[bestWindow], splitIndexes[bestWindow+1], bestIndex, minLength) 
            return {'Best Window': bestWindow, 'Best Index': bestIndex, 'M': bestScore, 'V': best_V_score}
        # otherwise, return None
        else:
            return None
   
    ### Actually do stuff ##############################################################################################
    splitIndexes = [0,data.length]

    # if maxIterations isn't set, set it to the max number of splits that are
    # possible given the minLength restriction
    if maxIterations == None: 
        maxIterations = int( ( 2.0*len(data) )/minLength )

    for loop in xrange(maxIterations):
        split_data = find_split( data, splitIndexes )

        if verbose: print >> output,  loop, splitIndexes, "\n", split_data
        if split_data == None: break

        splitIndexes.insert(split_data['Best Window']+1, split_data['Best Index'])

    return splitIndexes

def merge_boundaries(boundaries, regions_list, min_boundary_len):
    """This merges boundaries.

    """

    base_region = regions_list[0]

    # a method to flatten the boundaries list
    def flatten(lst):
        for elem in lst:
            if type(elem) in (tuple, list):
                for i in flatten(elem):
                    yield i
            else:
                yield elem

    # make the boundaries list unique
    boundaries = list(set(flatten(boundaries)))
    # sort the boundaries list
    boundaries.sort()
    
    # test for regions that are too small
    loop = 0
    while loop < len(boundaries)-1:
        if verbose:
            print "Curr Boundaries: ", boundaries

        rs = boundaries[loop]
        re = boundaries[loop+1] 

        if re - rs + 1 < min_boundary_len:
            # if rs == 0, then we MUST put this in with the next bucket
            if rs == 0:
                boundaries.remove(re)
            # if re == length, then we MUST put this in with the next bucket
            elif re == boundaries[-1]:
                boundaries.remove(rs)
            else:
                # try and put this bucket with the bucket with the closest mean
                # with respect to base_region
                prev_reg = base_region[boundaries[loop-1]:boundaries[loop]]
                mean1 = prev_reg.featuresLength()/float( boundaries[loop] - boundaries[loop-1] + 1 )

                next_reg = base_region[boundaries[loop]:boundaries[loop+1]]
                mean2 = next_reg.featuresLength()/float( boundaries[loop+1] - boundaries[loop] + 1 )

                current_mean = base_region[rs:re].featuresLength()/float( re - rs + 1 )
                
                if verbose:
                    print "Lower Mean: %.5f   Current Mean: %.5f    Upper Mean: %.5f" % (mean1, current_mean, mean2)

                if abs(current_mean - mean1) <  abs(current_mean - mean2):
                    boundaries.remove(rs)
                else:
                    boundaries.remove(re)
            # retart the loop
            loop = 0
        else:
            #increment the loop
            loop += 1

    return boundaries


def usage():
  print >> output, """
  Split a bed file into more regions and writes the result into a new bed file.

  INPUTS:
  -i --input:  The regions to split, in a .bed file.

   alternatively, to split two bed files by their combined split points

  -1 --input1: region1 to split, in a .bed file.
  -2 --input2: region2 to split, in a .bed file.

  -d --domain : the domain of the data files, the portion of the genome over 
     which these features (-1 and -2) are defined.  The support of the 
     statistics. Usually determined by array coverage or "alignibility".  If the
     features are defined everywhere (e.g. such as may be the case in C. elegans
     data), then this file contains one line for each chromosome, with: 
     "chr_name 0 chr_length" on each line.

  -m --min_length: For the "double bootstrap" in any organism this should be at 
     least 5,000,000.  For human or mouse, it should be at least 1,000,000 in 
     general, for Drosophila at least 500,000, and for C. elegans, at least 
     250,000.  Larger values will cause down stream analysis to be more 
     conservative, so initial runs should use quite large values.  If 
     significance is obtained with larger values, it will certainly hold for 
     smaller ones.
    
  -s --max_splits: The maximum number of times to split a region. Defaults to
     no maximum. Recommended: no maximum.

  -p --plot: Plot the CDF of the region data and overlay the calculated split
     points. This option requires that the pylab module is part of your python 
     distribution ( usually part of matplotlib ) and that you are running this
     locally.

  -g --generate_test_data: Generates a region of test data and writes it to
     sim_output.bed and sim_lengths.txt. The region will be approximately 500k
     BP's long and consist of intervals with normally distributed lengths of 
     mean 1k BP's and SD's of 200 BP's. The region will consist of 5 ( unnamed )
     subregions of length 100k and with means (0.1, 0.3, 0.5, 0.7, 0.9). For 
     finer control over the generation procedure, load this as a module and
     call the function generateSequences directly.

  -o --output_file_prefix: A file prefix to name the output files. For instance, 
     if the input file names are test.bed and test.txt, and the prefix is split_
     the files will be written to split_test.bed and split_lengths.bed. The
     prefix defaults to split_. THIS WILL SILENTLY OVERWRITE ANY EXISTING FILES 
     OF THE SAME NAME !

  -v --verbose

  -w --supresswarnings: Supress merge warnings.
  """

def main():
    try:
        long_args = ["input=", "input1=", "input2=", "domain=", "min_length=", "max_splits=", "output_file_prefix=", "verbose", "plot", "generate_test_data", "help", "supresswarnings"]
        opts, args = getopt.getopt(sys.argv[1:], "i:1:2:d:m:s:o:vpghw", long_args)
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    # set the default options
    max_splits = None
    output_fname_prefix = "split_"
    do_plot = False

    for o, a in opts:
        if o in ("-v", "--verbose"):
            global verbose
            verbose = True
        elif o in ("-w", "--supresswarnings"):
            global supress_warnings
            supress_warnings = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-i", "--input"):
            split_file = open(a)
        elif o in ("-1", "--input1"):
            split_file1 = open(a)
        elif o in ("-2", "--input2"):
            split_file2 = open(a)
        elif o in ("-d", "--domain"):
            lengths_file = open(a)
        elif o in ("-m", "--min_length"):
            min_length = int(a)
        elif o in ("-s", "--max_splits"):
            max_splits = int(a)
        elif o in ("-p", "--plot"):
            try: 
                global pylab
                import pylab
            except ImportError:
                raise ImportError, "Must have pylab installed to use plot option."
            else:
                do_plot = True
        elif o in ("-g", "--generate_test_data"):
            print "Saving the simulated region as 'sim_output.bed' and 'sim_lengths.txt'..."
            generateSequences(number=2)
            print"\nfinished.\n"
            sys.exit()
        elif o in ("-o", "--output_file_prefix"):
            output_fname_prefix = a
        else:
            assert False, "unhandled option %s" % o

    try: 
        assert vars().has_key('split_file') or ( vars().has_key('split_file1') and vars().has_key('split_file2') )
        assert vars().has_key('lengths_file')
        assert vars().has_key('min_length') 
    except Exception, inst:
        print inst
        usage()
        sys.exit()

    if vars().has_key('split_file'):
        split_files = ( split_file, )
    else:
        split_files = (split_file1, split_file2)

    # build the regions list
    regions_s_to_split = [parse_bed_file(fp, lengths_file) for fp in split_files ]
    if verbose: 
        print >> output, "Region files parsed."

    # close all of the open files
    for fp in split_files:
        fp.close()
    lengths_file.close()

    if verbose:
        import time
        startTime = time.time()

    baseRegions = []
    for i in xrange(len(regions_s_to_split)):
        test = base_types.regions()
        baseRegions.append(test)

    # for each named region in all of the regions
    for key in regions_s_to_split[0].keys():
        # store the region boundaries for each area
        region_boundaries = []

        # a list of region's for this key
        region_list = [ regions[key] for regions in regions_s_to_split ]
        # make sure all of the lengths are the same
        # BAD python2.3 compat change
        assert len(set([ region.length for region in region_list ])) == 1


        # for each regions in the regions
        # note that we require the regions 'line up'
        # the named regions should have identical names and lengths
        for data in region_list:
            cd = data.build_cdf()
            if verbose:
                print >> output, "Cumulative Data object built for %s" % key
            region_boundaries.append( segment( cd, min_length, b=0, maxIterations=max_splits) )
            if verbose:
                print >> output, "\n", region_boundaries

            if do_plot:
                cd.plot( region_boundaries[-1] )
                print "Press enter to continue..."
                raw_input()

        # merge the region boundaries
        if len(region_list) == 1:
            boundaries = region_boundaries[0]
        else:
            boundaries = merge_boundaries( region_boundaries, region_list, min_length )

        if do_plot and len(split_files) > 1:
            for data in region_list:
                data.build_cdf().plot( boundaries )
                print "Press enter to continue..."
                raw_input()

        if verbose: 
            print boundaries

        for data, baseRegion in zip(region_list, baseRegions):
            # note that we dont put the beginning and end into the split argument
            split_data = data.split(boundaries[1:-1])
            baseRegion.extend(split_data)

    for fp, baseRegion in zip(split_files, baseRegions):
        ### FIX ME - i shouldnt be writing and overwriting the lengths file
        # write the output to a bed file
        of = open(output_fname_prefix + fp.name, 'w')
        olf = open(output_fname_prefix + lengths_file.name, 'w')
        baseRegion.writeBedFile(of, olf)
        of.close()
        olf.close()

    if verbose: print >> output, "Execution Time: ", time.time()-startTime

if __name__ == "__main__":
    main()
