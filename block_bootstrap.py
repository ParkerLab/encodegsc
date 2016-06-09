#!/usr/bin/python

import getopt, sys

import base_types
from base_types import *

verbose = False
output = sys.stdout

import string
import itertools

from random import uniform, randint
from math import sqrt
from operator import itemgetter

rpy = None
try:
    import rpy
except ImportError:
    pass

class sample(tuple):
    """Hold sample data from which we will calculate the test stat.

    """
    def __init__(self, *args, **kwargs):
        # tuple.__init__(self, *args, **kwargs)
        tuple.__init__(self)

    def __mul__(self, scalar):
        """Multiply a sample by a scalar.

        eg 0.5*(1,1,1,1,2,1) = (0.5,0.5,0.5,0.5,1.0,0.5)
        """
        # BAD python 2.3 compat change
        return type(self)( [scalar*item for item in self] )

    def __add__(self, other):
        """Add two samples together pointwise.

        eg (0,0,0,0,0,0) + (1,1,1,1,2,1) = (1,1,1,1,2,1)

        """
        if len(self) != len(other):
            raise ValueError, "Can only add samples of the same length."
        
        return type(self)( [ i1 + i2 for i1, i2 in itertools.izip( self, other ) ] )

class sample_dict(dict):
    """Store a dictionary of samples.
    
    This is generally used to store all of the samples for a region, categorized
    by their named partitions. An example of this would be looking at an entire,
    genome, but having it partitioned by chromosome. The keys would be the chr
    names.
    """
    def __setitem__(self, key, item):
        #assert isinstance(item, sample)
        dict.__setitem__(self, key, item)   

    def sum(self):
        """Returns the sum of all the overlap_sample objects.

        In example, if the values in this container were
            (0,0,0,0,2,4)
        and (1,4,5,7,2,4)

        this would return (1,4,5,7,4,8).
        """

        if len(self) == 0: 
            raise ValueError, "Can't sum an empty sequence"

        # define a container to hold the values
        empty = type(self.values()[0])(0,0,0,0,0,0)

        # iterate through the items, and add them to the empty container
        for item in self.values():
            empty += item

        return empty

    def weightedSum(self, weights):
        """Returns the sum of all overlap_sample objects weighted by weights.

        In example, if the values in this container were
            'R1': (0,0,0,0,2,4)
        and 'R2': (1,4,5,7,2,4)

        and the weights were 
            'R1': 0.1   
        and 'R2': 0.9

        then this would return 

           (0, 0, 0, 0, .2, .4) + (0.9, 3.6, 4.5, 6.3, 1.8, 3.6)
        =  (0.9, 3.6, 4.5, 6.3, 2.0, 4.0)

        """

        if len(self) == 0: 
            raise ValueError, "Can't sum an empty sequence"

        # make sure that the keys are identical
        assert set(weights.keys()) == set(self.keys())

        # make sure that the weights sum to 1
        assert round( sum( weights.values() ), 15) == 1.0

        # define a container to hold the values
        obj = self.values()[0]
        obj_type = type(obj)
        if obj_type in (int,float):
            empty = obj_type(0)
        else:
            empty = obj_type(args=[0]*len(obj))

        # BAD python2.3 compatability change
        return sum( [ self[key]*weights[key] for key in self.keys() ], empty )

def double_overlap( coveredRegion, coveringRegion, 
                    outer_regionFraction, inner_regionFraction, 
                    callback, number=1 ):
    """Double sample from a pair of regions and calculate aggregate stats.

    This is a wrapper for custom GSC statistics. This code samples from the
    regions with the specified region fractions, and then calls callback with
    s11,s12, s21, s22 for the ouer and inner samples. The callback should return 
    an object inherited from sample. Then, double overlap sample will return
    a sample dict containing the samples hashed by region name.
    
    It's probably easiest to quickly browse the code and an example.
    """
    osl = outer_regionFraction
    isl = inner_regionFraction

    # for loop in the number of samples to take
    for loop in xrange(number):
        # first, build a dict of the outer samples
        # these are what we use to calculate the expectation
        #   of the statistic
        osample = sample_dict()
        isample = sample_dict()

        # for each named region in the coveringRegion...
        for key in coveringRegion.keys():
            # rename the regions to be more readable
            cA = coveringRegion[key]
            eA = coveredRegion[key]
            if cA.length != eA.length:
                raise ValueError, "Regions must be the same length in bp's."
            
            ###############################
            # Take the outer samples 
            ###############################
            # find the length of the sample in bp's
            sample_length = int(eA.length*osl)

            # take the slices
            start_bp = randint(0, cA.length - sample_length - 1)
            os1 = eA.get_subregion( start_bp, start_bp+sample_length, shift_to_zero=True )
            start_bp = randint(0, cA.length - sample_length - 1)
            os2 = cA.get_subregion( start_bp, start_bp+sample_length, shift_to_zero=True )

            # save the outer sample stat
            osample[key] = callback( os1, os2 )
            
            ###############################
            # Take the (inner) sub samples 
            ###############################
            # find the length of the sample in bp's
            sample_length = int(sample_length*isl)
            
            # take the slices
            start_bp = randint(0, os1.length - sample_length - 1)
            is1 = eA.get_subregion( start_bp, start_bp+sample_length, shift_to_zero=True )
            start_bp = randint(0, os1.length - sample_length - 1)            
            is2 = cA.get_subregion( start_bp, start_bp+sample_length, shift_to_zero=True )
            
            # save the inner sample stat
            isample[key] = callback( is1, is2 )
        yield osample, isample

def single_overlap( coveredRegion, coveringRegion, 
                    regionFraction, callback, number=1 ):
    """Block sample from a pair of regions and calculate aggregate stats.

    This is a wrapper for custom GSC statistics. This code samples from the
    regions with the specified region fractions, and then calls callback with
    s11,s12, s21, s22 for the otuer and inner samples. The callback should return 
    an object inherited from sample. Then, single overlap sample will return
    a sample dict containing the samples hashed by region name.
    
    It's probably easiest to quickly browse the code and an example - correlation
    is a good placde to start.
    """
    sl = regionFraction
    
    # for loop in the number of samples to take
    for loop in xrange(number):
        sample = sample_dict()

        # for each named region in the coveringRegion...
        for key in coveringRegion.keys():
            # rename the regions to be more readable
            cA = coveringRegion[key]
            eA = coveredRegion[key]
            if cA.length != eA.length:
                raise ValueError, "Regions must be the same length in bp's."
            
            ###############################
            # Take the outer samples 
            ###############################
            # find the length of the sample in bp's
            sample_length = int(eA.length*sl)

            # take the slices
            start_bp = randint(0, cA.length - sample_length - 1)
            s1 = eA.get_subregion( start_bp, start_bp+sample_length, shift_to_zero=True )
            start_bp = randint(0, cA.length - sample_length - 1)
            s2 = cA.get_subregion( start_bp, start_bp+sample_length, shift_to_zero=True )

            # save the outer sample stat
            sample[key] = callback( s1, s2 )
        
        yield sample



def conditional_bp_overlap_stat( \
    covering_region, covered_region, region_fraction, num_samples ):
    #% Lastly, scale the bootstrap distribution to get the null distribution of
    #% the test statistic.
    #
    #% 2*frac = 2L/n.  So sqrt(2*frac)*T_n has the correct SD.  We can compute
    #% this by pulling the constant out:
    #
    #% store the samples from n and N
    #%SD = sqrt(2*frac)*std(Tn);
    #

    Tns = []
    Tns_2 = []
    # calculate the overlap stat num_samples times
    samples = random_regions_bp_overlap(covered_region, covering_region, region_fraction, num_samples)
    # BAD python2.3 compat change
    lengths = dict([ (key, covered_region[key].length) for key in covered_region.keys() ])
    
    if verbose:
        print >> output, "#### Sample Distribution Info"
        print >> output, "Sample #".rjust(8), "Tn".rjust(15)
    
    loop = 0
    for sample in samples:
        stats = conditional_overlap_sample_stat(sample, lengths, region_fraction)
        Tns.append(stats['Tn'])
        if verbose: print >> output, str(loop).rjust(8), str("%.8f" % stats['Tn']).rjust(15)
        loop += 1
    
    SD = sqrt(2*region_fraction)*std(Tns)
    if verbose:
        print >> output
        print >> output, 'mean: ', mean(Tns)
        print >> output, 'SD:   ', SD, '\n'

    real_stats = conditional_overlap_stat(covered_region, covering_region) 

    theoretical_stat_mean = real_stats['theoretical']
    observed_stat_mean = real_stats['observed']
    test_stat = real_stats['test_stat']

    if verbose:
        print >> output, 'NULL stat mean: ', theoretical_stat_mean
        print >> output, 'observed stat mean: ', observed_stat_mean      
        print >> output, 'test_stat: ', test_stat

    # Finally, refer the "mean zero" test statistic "test_stat" to the
    # distribution Tn.  Compute the SD and whatnot.
    #
    #%z_score = test_stat/SD;
    #%p_value = 1 - normcdf(test_stat,0,SD);

    z_score = test_stat/SD
    if verbose: print >> output, 'z_score: ', z_score

    p_value = min(1 - sn_cdf(z_score), sn_cdf(z_score))
    print >> output, 'p_value: ', p_value, "\n"

    return z_score, p_value



def double_bootstrap_stat(      covering_regions, covered_regions, 
                                outer_regionFraction, inner_regionFraction, 
                                aggregate_sample_callback,
                                regions_agg_callback, 
                                nsamples=1, 
                                real_stat_callback=None    ):
    """Wrapper for generic double bootstrap stat

    Input Parameters:
        covering_regions - the regions that num regions is calulated from
        covered_regions - the other regions

        outer_region_fraction - the ratio of the outer sample
        inner_region_fraction - the ratio of the subsample to the outer sample

        aggregate_samples_callback: this is what is called on the block samples 
        that are taken. the quintesential example of this is region overlap. 

        scale_sample_callback: scales the stat that was aggregated from 
        aggregate_samples_callback.
        
        regions_agg_callback: How to aggregate the samples. The two implemented 
        are the conditional and marginal version. Marginal adds up all of the 
        segments - the conditional weights them by the segment weights.

        nsamples - the number of samples to take
        
        real_stat_callback - sometimes we need to calculate the 'true' stat 
        differently than for the samples ( for instance, maybe it needsw to 
        be scaled differently. In such cases, pass this in. If it is none, the
        real stat is calculated by calling the sample functions on the whole 
        set of regions.

    """
   
    # store the samples
    o_samples = []
    i_samples = []

    if verbose:
        # TODO make this prettier, maybe with a sample header method?
        print >> output, "#### Samples\n"
        print >> output, "Sample #".rjust(8), "Outer Region Stat".rjust(23), \
        "Inner Region Stat".rjust(23)
        
    # take the samples and keep track of the relevant data
    for loop, item in \
        enumerate(double_overlap( 
            covered_regions, covering_regions, 
            outer_regionFraction, inner_regionFraction, 
            aggregate_sample_callback, nsamples 
        )):

        # unpack the samples tuple
        outer_samples, inner_samples = item                     
        # append the aggregated stats, should be a floats
        o_samples.append(regions_agg_callback(outer_samples))
        i_samples.append(regions_agg_callback(inner_samples))
        
        #print outer_samples, inner_samples
        #print agg_callback(outer_samples), agg_callback(inner_samples)
        if verbose:
            print >> output, str(loop+1).rjust(8), \
               str(o_samples[-1]).rjust(23), str(i_samples[-1]).rjust(23)
      
    # the expectation of the mean under the NULL
    # our estimate is the empirical mean of the outer sample
    Nmean = mean(o_samples)
    obsMean = mean(i_samples)

    if rpy != None and verbose:
        rpy.r.png( "dist.png", width=600, height=600  )
        rpy.r.par( mfrow=(2,2) )
        rpy.r.hist( o_samples, main="Outer Sample", xlab="", ylab=""  )
        rpy.r.hist( i_samples, main="Inner Sample", xlab="", ylab=""  )
        rpy.r.qqnorm( o_samples, main="Outer Sample", xlab="", ylab=""  )
        rpy.r.qqline( o_samples  )
        rpy.r.qqnorm( i_samples, main="Inner Sample", xlab="", ylab=""  )
        rpy.r.qqline( i_samples  )
        rpy.r.dev_off()

    dist = [ (sample - Nmean) for sample in i_samples ]
    dist_sd = sqrt(outer_regionFraction*inner_regionFraction)*std(dist)

    if verbose: print >> output, "\nSample Dist Info:\n"
    if verbose: print >> output, "\tThe expected mean under the NULL: ", Nmean
    if verbose: print >> output, "\tThe sample mean: ", obsMean
    if verbose: print >> output, "\tThe normalized bootstrap SD: ", dist_sd, "\n"

    # the stat on each of the segments
    if None == real_stat_callback:
        seg_stats = sample_dict([ (key, aggregate_sample_callback( covered_regions[key], covering_regions[key] )) for key in covered_regions.keys() ])
        real_stat_expectation  = regions_agg_callback( seg_stats )
    else:
        real_stat_expectation = real_stat_callback( covered_regions, covered_regions  )

    obs_mean = real_stat_expectation - Nmean
    if verbose: print >> output, "\tThe real stat: ", real_stat_expectation  
    if verbose: print >> output, "\tThe scaled real stat: ", obs_mean, "\n"    

    z_score = obs_mean/dist_sd
    print "\tZ Score: ", z_score

    p_value = min( 1 - sn_cdf( z_score ), sn_cdf( z_score ) )
    print "\tp-value: ", p_value
    
    return z_score, p_value

def single_bootstrap_stat(      covering_regions, covered_regions, 
                                regionFraction, 
                                exp_under_null_callback,
                                aggregate_sample_callback,
                                regions_agg_callback, 
                                nsamples=1,
                                real_stat_callback=None):
    """Wrapper for generic double bootstrap stat

    Input Parameters:
        covering_regions - the regions that num regions is calulated from
        covered_regions - the other regions

        region_fraction - the fraction of the region to sample

        exp_under_null_callback - determine the expectation of the stat under 
        the null, as a function of the covering and covered regions. 

        aggregate_samples_callback: this is what is called on the block samples 
        that are taken. the quintesential example of this is region overlap. 

        scale_sample_callback: scales the stat that was aggregated from 
        aggregate_samples_callback.
        
        regions_agg_callback: How to aggregate the samples. The two implemented 
        are the conditional and marginal version. Marginal adds up all of the 
        segments - the conditional weights them by the segment weights.

        nsamples - the number of samples to take

    """
   
    # store the samples
    samples = []

    if verbose:
        # TODO make this prettier, maybe with a sample header method?
        print >> output, "#### Samples\n"
        print >> output, "Sample #".rjust(8), "Stat".rjust(23)
        
    # take the samples and keep track of the relevant data
    for loop, sample in \
        enumerate(single_overlap( 
            covered_regions, covering_regions, 
            regionFraction,
            aggregate_sample_callback, nsamples 
        )):

        # append the aggregated stats, should be a floats
        samples.append(regions_agg_callback(sample))
        
        #print outer_samples, inner_samples
        #print agg_callback(outer_samples), agg_callback(inner_samples)
        if verbose:
            print >> output, str(loop+1).rjust(8), \
               str(samples[-1]).rjust(23)
    
    # the expectation of the mean under the NULL
    # our estimate is the empirical mean of the outer sample
    Nmean = exp_under_null_callback( covering_regions, covering_regions )
    obsMean = mean(samples)

    if rpy != None and verbose:
        rpy.r.png( "dist.png", width=1200, height=600  )
        rpy.r.par( mfrow=(1,2) )
        rpy.r.hist( samples, main="Sample Hist", xlab="", ylab="", n=nsamples/10.0  )
        rpy.r.qqnorm( samples, main="qqplot vs norm quantiles", xlab="", ylab=""  )
        rpy.r.qqline( samples )
        rpy.r.dev_off()

    dist = [ (sample - Nmean) for sample in samples ]
    dist_sd = sqrt(regionFraction)*std(dist)

    if verbose: print >> output, "\nSample Dist Info:\n"
    if verbose: print >> output, "\tThe expected mean under the NULL: ", Nmean
    if verbose: print >> output, "\tThe sample mean: ", obsMean
    if verbose: print >> output, "\tThe normalized bootstrap SD: ", dist_sd, "\n"

    # the stat on each of the segments
    # the stat on each of the segments
    if None == real_stat_callback:
        seg_stats = sample_dict([ (key, aggregate_sample_callback( covered_regions[key], covering_regions[key] )) for key in covered_regions.keys() ])
        real_stat_expectation  = regions_agg_callback( seg_stats )
    else:
        real_stat_expectation = real_stat_callback( covered_regions, covering_regions  )

    obs_mean = real_stat_expectation - Nmean
    if verbose: print >> output, "\tThe real stat: ", real_stat_expectation  
    if verbose: print >> output, "\tThe scaled real stat: ", obs_mean, "\n"    

    z_score = obs_mean/dist_sd
    print "\tZ Score: ", z_score

    p_value = min( 1 - sn_cdf( z_score ), sn_cdf( z_score ) )
    print "\tp-value: ", p_value
    
    return z_score, p_value

def fraction_basepair_overlap( coveredRegionSample, coveringRegionSample ):
    """ Calculate region overlap statistics.
    
    This is intended as a callback for single_overlap.
    """
    # stores the feature length of the covered region ( self )
    covered_feature_len = coveredRegionSample.featuresLength()

    # stores the feature length of the covering region ( coveringRegion )
    covering_feature_len = coveringRegionSample.featuresLength()

    # stores the observed total feature overlap between the region's
    overlap_feature_len = coveredRegionSample.overlap(coveringRegionSample)

    return sample( ( overlap_feature_len, covered_feature_len, covering_feature_len  ) )

def conditional_bp_overlap_stat( \
    covering_regions, covered_regions, \
    region_fraction, num_samples ):
    
    def agg_callback(sample):
        total_length = sum( item.length for item in covering_regions.values()  )
        Fn = 0.0
        Jn_num = 0.0
        Jn_denom = 0.0
        for r_name, values in sample.iteritems():
            # unpack the sample values
            overlap, covered_feature_len, covering_feature_len = values
            regionLength = float(covering_regions[r_name].length)
            # lambda in the paper
            length_frac = regionLength/total_length
            # the length of the sampled region
            sampleRegionLength = region_fraction*regionLength
            
            I = covered_feature_len/sampleRegionLength
            J = covering_feature_len/sampleRegionLength
            IJ = overlap/sampleRegionLength

            Fn += length_frac*(IJ/max(I, 1.0/sampleRegionLength))
            Jn_denom += length_frac*I
            Jn_num += length_frac*I*J
            
        return Fn

    
    def analytical_expectation( covering_regions, covered_regions ):
        totalLength = sum( item.length for item in covering_regions.values()  )
        
        Obs_num = 0.0
        Obs_den = 0.0
        I_num = 0.0
        for key in covered_regions.keys():
            covered_feature_len = covered_regions[key].featuresLength()
            covering_feature_len = covering_regions[key].featuresLength()
            overlap_feature_len = covered_regions[key].overlap(covering_regions[key])

            # stores the length of this region for both features
            region_length = float(covered_regions[key].length)
            
            # this is lambda_i in the paper
            regionFraction  = region_length/totalLength
            
            Obs_num += regionFraction*overlap_feature_len/region_length
            Obs_den += regionFraction*covered_feature_len/region_length
            
            I_num += regionFraction*covering_feature_len*covered_feature_len/(region_length**2)

        J_n = I_num/Obs_den;
        O_n = Obs_num/Obs_den;
        
        return { 'expected': J_n, 'observed': O_n  }
    
    def expectation_under_null( covering_regions, covered_regions ):
        return analytical_expectation( covering_regions, covered_regions )[ 'expected'  ]

    def real_stat( covering_regions, covered_regions ):
        return analytical_expectation( covering_regions, covered_regions )[ 'observed'  ]
    
    z_score, p_value = \
        single_bootstrap_stat( covering_regions, covered_regions, \
                               region_fraction, expectation_under_null, \
                               fraction_basepair_overlap, agg_callback, \
                               num_samples, real_stat )
    
    return z_score, p_value

def marginal_bp_overlap_stat( \
    covering_regions, covered_regions, \
    region_fraction, num_samples ):
    
    weights = covering_regions.regionFraction()  
    
    def agg_callback(sample):
        assert set( sample.keys() ) == set( covering_regions.keys() )
        assert set( sample.keys() ) == set( covered_regions.keys() )
        tot_overlap = 0
        tot_fl = 0
        for r_name, values in sample.iteritems():
            overlap_feature_len, covered_feature_len, covering_feature_len = values
            tot_overlap += float(overlap_feature_len)
            tot_fl += covered_feature_len

        # If we want to do this conditional on features, uncomment this
        """
        if tot_fl == 0: 
            return 0
        else:
            return tot_overlap/tot_fl
        """        
        return tot_overlap/tot_fl

    def expectation_under_null( covering_regions, covered_regions ):
        covered_fl = 0
        covering_fl = 0
        region_lengths = 0

        for r_name in covering_regions.keys():
            covered_fl += covered_regions[r_name].featuresLength()
            covering_fl += covering_regions[r_name].featuresLength()
            region_lengths += covered_regions[r_name].length
        
        # Note that this is really (covering*covered/rl**2)/covered,
        # which reduces algebraically to the below expression
        return covering_fl/float(region_lengths)
    
    z_score, p_value = \
        single_bootstrap_stat( covering_regions, covered_regions, \
                               region_fraction, expectation_under_null, \
                               fraction_basepair_overlap, agg_callback, \
                               num_samples )
    
    return z_score, p_value

def conditional_resample_region_overlap_stat(
                covering_regions, covered_regions, 
                outer_regionFraction, inner_regionFraction, 
                nsamples, min_overlap_fraction=0.0
        ):
    """Calculate the probability under the NULL of independence that the region
    overlap is purely do to chance.

    I use the machinery in double bootstrap stat to implement this.
    """
    # the region weights - defined as homogenous region length/ total length
    weights = covering_regions.regionFraction()  
    
    def agg_callback(sample):
        assert set( sample.keys() ) == set( covering_regions.keys() )
        assert set( sample.keys() ) == set( covered_regions.keys() )
        return sample.weightedSum( weights )

    def region_overlap( coveredRegionSample, coveringRegionSample ):
        """ Calculate region overlap statistics.

        This is intended as a callback for double_overlap.
        """
        percent_overlap = float(coveredRegionSample.regionOverlap(coveringRegionSample, min_overlap_fraction)) / coveredRegionSample.numRegions
        return percent_overlap

    return double_bootstrap_stat(covering_regions, covered_regions, 
                                outer_regionFraction, inner_regionFraction, 
                                region_overlap, agg_callback, 
                                nsamples=nsamples)
    
def marginal_resample_region_overlap_stat(
                covering_regions, covered_regions, 
                outer_regionFraction, inner_regionFraction, 
                nsamples, min_overlap_fraction=0.0
        ):
    """Calculate the probability under the NULL of independence that the region
    overlap is purely do to chance.

    I use the machinery in double bootstrap stat to implement this.
    """
    def agg_callback(sample):
        assert set( sample.keys() ) == set( covering_regions.keys() )
        assert set( sample.keys() ) == set( covered_regions.keys() )

        overlaps = sum([ item[0] for item in sample.values() ])
        num_regions = sum([ item[1] for item in sample.values() ])
        
        return float(overlaps)/num_regions

    def region_overlap( coveredRegionSample, coveringRegionSample ):
        """ Calculate region overlap statistics.

        This is intended as a callback for double_overlap.
        """
        region_overlap = ( float(coveredRegionSample.regionOverlap(coveringRegionSample, min_overlap_fraction)), coveredRegionSample.numRegions )
        return region_overlap

    return double_bootstrap_stat(covering_regions, covered_regions, 
                                outer_regionFraction, inner_regionFraction, 
                                region_overlap, agg_callback, 
                                nsamples=nsamples)

def conditional_pearson_correlation(
                covering_regions, covered_regions, 
                regionFraction, nsamples
        ):
    """Calculate the probability under the NULL of independence that the pearson
    correlation is purely do to chance.

    I use the machinery in double bootstrap stat to implement this.
    """
    
    def agg_callback(sample):
        """Correlation is asmytotically normal, so we weight the correlation 
           by 1/sqrt(region length), since length is the number of
           points that contribute to the correlation.
        
        """
        stat = 0
        for key, val in sample.iteritems():
            length = covering_regions[key].length
            stat += val/sqrt( length  )
        return stat

    def pearson_correlation( coveredRegionSample, coveringRegionSample ):
        """ Calculate region overlap statistics.

        This is intended as a callback for double_overlap.
        """
        correlation = coveredRegionSample.regionCorrelation(coveringRegionSample)
        return correlation

    def mean_under_null( covering_regions, covered_regions ):
        return 0.0

    return single_bootstrap_stat(covering_regions, covered_regions, 
                                 regionFraction, 
                                 mean_under_null,
                                 pearson_correlation, 
                                 agg_callback, 
                                 nsamples=nsamples)

def conditional_mean_fold_enrichment(
                covering_regions, covered_regions, 
                outer_regionFraction, inner_regionFraction, 
                nsamples
        ):
    """Calculate the probability under the NULL of independence that the fold
    of region 1 WRT region2 is due to chance.

    I use the machinery in double bootstrap stat to implement this.
    """
    # the region weights - defined as homogenous region length/ total length
    weights = covering_regions.regionFraction()  
    
    def agg_callback(sample):
        assert set( sample.keys() ) == set( covering_regions.keys() )
        assert set( sample.keys() ) == set( covered_regions.keys() )
        return sample.weightedSum( weights )

    def fold( coveredRegionSample, coveringRegionSample ):
        # special case for empty intervals        
        if len(coveredRegionSample) == 0 or len(coveringRegionSample) == 0: return 0

        totalEnrichment = 0

        currentMatches = []
        
        thisIter = coveredRegionSample.iter_intervals_and_values()
        nextCoveredMatch = thisIter.next()
        
        for cg_iv, cg_val in coveringRegionSample.iter_intervals_and_values():
            # first, get rid of any current matches that dont match
            # because of the ordering, these should start not matching at the begining
            for cd_iv, cd_val in currentMatches:
                if cd_iv.end >= cg_iv.start: break
                else: del currentMatches[0] 

            # next, add any new items to currentMatches
            while nextCoveredMatch != None and nextCoveredMatch[0].start <= cg_iv.end:
                currentMatches.append(nextCoveredMatch)
                try: nextCoveredMatch = thisIter.next()
                except StopIteration: nextCoveredMatch = None

            # finally, calculate the fold enrivchment of every item in currentMatches
            for cd_iv, cd_val in currentMatches:
                totalEnrichment += cd_iv.overlap(cg_iv)*( float(cg_val)/cd_val )

        return totalEnrichment/coveredRegionSample.featuresLength()

    return double_bootstrap_stat(   covering_regions, covered_regions, 
                                    outer_regionFraction, inner_regionFraction, 
                                    fold, agg_callback, 
                                    nsamples=nsamples   )

def usage():
  print >> output, """
  This calculates a p-value under the null hypothesis that the observed
  instance-coverage of Feature_1 by Feature_2 is due to chance.

  INPUTS:
  -1 : --feature_1 : the 'covered' input data file
  -2 : --feature_2 : the 'covering' input data file
  Input files are accepted in either the wiggle or bed format. If the filename
  extension is .wig, the wiggle file is parsed. Otherwise, the file is assumed
  to be a bed file. For asymmetric stats, the covered file is what the numerator
  in the test statistic should be. ie, for fold enrichment, the stat is the 
  quotient covering_region/covered_region. Details are in the doc strings.

  -d : --domain : the domain of the data files, the portion of the genome over 
  which these features (-1 and -2) are defined.  The support of the statistics.  
  Usually determined by array coverage or "alignibility".  If the features are 
  defined everywhere (e.g. such as may be the case in C. elegans data), then 
  this file contains one line for each chromosome, with: "chr_name 0 chr_length"
  on each line.

  The file formats are described in the README under INPUT FILE FORMATS.

  -r --region_fraction : the fraction of each region (e.g. chromosome) to take 
  in each block-wise sample.  The product of this number, and the minimum 
  segment length (e.g. for no segmentation, using whole human chromosomes, the  
  minimum segment length would be about 50Mb due to chromosome Y), should be  
  larger than the mixing distance of features -1 and -2.  For human, if the  
  minimum segment length is around 5Mb (some segmentation has been done, or the  
  features of interest are not defined everywhere throughout the genome, e.g.  
  due to mappability issues), then this number should be at least 0.01 for most   
  features. This number can be fitted under a stability criterion.
  
  -s --subregion_fraction: only used for the double bootstrap tests, 
  this specifies the fraction of the subsample to take. For instance, a 
  specified region fraction (-r) of 0.20 and a subregion fraction (-s) of 0.20
  would yield a net sample length of 0.04 of each named region. 

  -n --num_samples : the number of bootstrap samples to take.  100 will be 
  sufficient to get an idea of significance, but at least 10K should be used 
  for publication.

  -t --test: Determine the test to run. Accepts one of the following types
     basepair_overlap_conditional ( bc ) ( default )
     basepair_overlap_marginal ( bm )
     region_overlap_conditional ( rc )
     region_overlap_marginal ( rm )
     pearson_correlation_conditional ( cc )
     fold_enrichment_conditional ( fc )

     The particulars of the tests are discussed in the code ( docstrings ). 
  
  -fg --filter_covering: Filter covering bed file by a group name.
  -fd --filter_covered : Filter covered bed file by a group name. In the bed4 
   format, the 4th column stores a label for the given region.  When one of the
   above two options are chosen, the file will filter out regions that dont 
   match 'name' as given. Defaults to no filtering.

  -B --force_binary: Treat the input data as binary regardless of whether or not
  there is a value. This is useful for bed files which have value's for the 
  purpose of plotting at the UCSC genome browser ( ie all 1000 ).
  
  --min_overlap_fraction: This is only used for the region overlap test. If the 
  covered region isn't overlapped by at least (X*100)% of the covering region,
  then we dont consider this an overlap. Defaults to 0.0. That is, if the covered 
  region is overlapped by even 1 basepair, then we say it is covered.
  
  -o --output_file: A file to output all non-error output into. Will append to 
  this file if it already exists. defaults to standard out.

  """

def main():
    try:
        long_args = [   "help", "feature-1=", "feature-2=", "domain=", 
                        "region-fraction=", "subregion-fraction=", 
                        "num-samples=", "test=", "output-file=", 
                        "filter-covering=", "filter-covered=", "force-binary", "verbose", "min-overlap-fraction="    ]
        
        opts, args = getopt.getopt(sys.argv[1:], "1:2:d:r:s:n:t:o:fg:fd:Bv", long_args)
    except getopt.GetoptError, err:
        usage()
        print "Error: \n", str(err)
        print
        sys.exit(2)

    ## make percent_basepair_overlap the default test type
    test_type = 'bc'
    filter_covering = None
    filter_covered = None
    force_binary = False
    min_overlap_fraction = 0.0

    for o, a in opts:
        if o == "-v":
            global verbose
            verbose = True
            base_types.verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-1", "--feature-1"):
            covered_file = open(a)
        elif o in ("-2", "--feature-2"):
            covering_file = open(a)
        elif o in ("-d", "--domain"):
            lengths_file = open(a)
        elif o in ("-r", "--region-fraction"):
            region_fraction = float(a)
        elif o in ("-s", "--subregion-fraction"):
            subregion_fraction = float(a)
        elif o in ("-n", "--num-samples"):
            num_samples = int(a)
        elif o in ("-t", "--test"):
            test_type = a
        elif o in ("-B", "--force-binary"):
            force_binary = True
        elif o in ("-o", "--output_file"):
            global output
            output = open(a, 'a')
        elif o in ("-fg", "--filter-covering"):
            filter_covering = a
        elif o in ("-fd", "--filter-covered"):
            filter_covered = a
        elif o in ("--min-overlap-fraction", ):
            min_overlap_fraction = float( a )
            if min_overlap_fraction < 0 or min_overlap_fraction > 1.0:
                raise ValueError, "min_overlap_fraction must be between 0 and 1.0, inclusive."
        else:
            assert False, "unhandled option %s" % o

    try: covered_file, covering_file, lengths_file, region_fraction, num_samples
    except Exception, inst:
        usage()
        sys.exit()

    if verbose:
        import time
        startTime = time.time()

    if covered_file.name.endswith(".wig"):
        coveredAnnotations = parse_wiggle_file(covered_file, lengths_file, force_binary)
    else:
        coveredAnnotations = parse_bed_file(covered_file, lengths_file, filter_covered, force_binary)

    if covering_file.name.endswith(".wig"):
        coveringAnnotations = parse_wiggle_file(covering_file, lengths_file, force_binary)
    else:
        coveringAnnotations = parse_bed_file(covering_file, lengths_file, filter_covering, force_binary)

    covered_file.close()
    covering_file.close()
    lengths_file.close()

    if verbose:
        print "Input Files Parse Time: ", time.time() - startTime, "\n"
        startTime = time.time()

    if not vars().has_key('region_fraction'):
        raise ValueError, 'The region_fraction setting is not set - it is mandatory.'

    #### Tests that work as for continuous data *or* for binary data
    if test_type in ( 'pearson_correlation_conditional', 'cc' ): 
        conditional_pearson_correlation( 
            coveringAnnotations, coveredAnnotations, 
            region_fraction, num_samples 
        )
    elif test_type in ( 'fold_enrichment_conditional', 'fc' ): 
        if not vars().has_key('subregion_fraction'):
            raise ValueError, 'The fold_enrichment_conditional test requires that the subregion_fraction option be set.'
        conditional_mean_fold_enrichment( 
            coveringAnnotations, coveredAnnotations, 
            region_fraction, subregion_fraction, num_samples 
        )
    else:
    #### Tests that ONLY work for binary data
        # test for the correct object
        if type(coveredAnnotations.values()[0]) != binary_region \
           or type(coveringAnnotations.values()[0]) != binary_region:
            raise TypeError, "The %s test only works for binary data ( You can force this with the -B option )" % test_type
        if test_type in ( 'basepair_overlap_marginal', 'bm' ) : 
            marginal_bp_overlap_stat( coveringAnnotations, coveredAnnotations, 
                                      region_fraction, num_samples )
        elif test_type in ( 'basepair_overlap_conditional', 'bc' ) : 
            conditional_bp_overlap_stat( coveringAnnotations, coveredAnnotations, 
                                         region_fraction, num_samples )
        elif test_type in ( 'region_overlap_marginal', 'rm' ): 
            if not vars().has_key('subregion_fraction'):
                raise ValueError, 'The region_overlap_marginal test requires that the subregion_fraction option be set.'
            marginal_resample_region_overlap_stat( 
                coveringAnnotations, coveredAnnotations, 
                region_fraction, subregion_fraction, 
                num_samples, min_overlap_fraction
            )
        elif test_type in ( 'region_overlap_conditional', 'rc' ): 
            if not vars().has_key('subregion_fraction'):
                raise ValueError, 'The region_overlap_conditional test requires that the subregion_fraction option be set.'
            conditional_resample_region_overlap_stat( 
                coveringAnnotations, coveredAnnotations, 
                region_fraction, subregion_fraction, 
                num_samples, min_overlap_fraction 
            )
        else:
            raise ValueError, "Unrecognized test '%s'" % test_type     

    if verbose: print >> output, "\nExecution Time: ", time.time()-startTime

    output.close()

if __name__ == "__main__":
    main()


