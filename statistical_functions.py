##################### Standard Statistical Methods #############################
### Functions we need if numpy isn't installed
###
################################################################################

try: from numpy.random import multivariate_normal
except ImportError:
    def multivariate_normal(meanMat, covMat, n):
        """ this constructs 2D multivariate dist, given a mean and cov list 
        """

        from math import sqrt
        from random import gauss as normal

        # make sure that the mean and cov entries make sense
        assert len(meanMat) == 2
        assert len(covMat) == 2
        assert len(covMat[0]) == 2
        assert len(covMat[1]) == 2

        # make sure that cov matrix is symmetric
        assert covMat[0][1] == covMat[1][0]

        # first, calculate the cholesky decomposition of the cov matrix
        # I'll use an explicit Cholesky-Crout algorithm since I know it's 2x2
        L11 = sqrt( covMat[0][0] )
        L12 = covMat[0][1]/L11
        # I add +0.000001 to covMat[1][1] - L12**2 because rounding errors sometimes
        # give me very small negative values, and then sqrt raises a domain error
        # I could use a try-catch, but I like to avoid the overhead since this
        # is a potentially hot spot in the algorithm - and it doesnt really matter
        L22 = sqrt( covMat[1][1] - L12**2 + 0.000001)

        mvn = [ [], [] ]
        # then I'll build a 2xn vector of MVN random samples
        for loop in xrange(n):
            z1 = normal(0,1)
            z2 = normal(0,1)

            mvn[0].append( meanMat[0] + L11*z1 )
            mvn[1].append( meanMat[1] + L12*z1 + L22*z2 )

        return zip(mvn[0], mvn[1])

try:
    from numpy import mean, cov, std
except ImportError: 
    try:
        from scipy.stats.stats import mean, cov, std
    except ImportError:
        from math import sqrt
        def mean(data): return float(sum(data))/len(data)
        
        def var(data):
            ev = mean(data)
            return sum( [ (entry - ev)**2 for entry in data ] )/( len(data) - 1 )
        
        def std(data): return sqrt( var(data) )
        
        def cov(data):
            """ Takes an iterable of 'vectors' ( by which I mean iterables of the same length )
            """
            # make sure they are all the same length
            for entry in data: assert len(entry) == len(data[0])
            
            # its expecting to get the entries in the opposite order that we want them
            data = zip(*data)
            
            # pairwise covariance
            def _cov(X,Y):
                assert len(X) == len(Y)
                mx = mean(X)
                my = mean(Y)
                return sum( [ float((x-mx)*(y-my)) for x,y in zip(X,Y) ] )/( len(X) - 1 )
            
            # first, make a zeros covariance matrix
            covMat = [ [0]*len(data) for loop in xrange(len(data)) ]
            # fill in the covaraince matrix entries
            for oloop in xrange(len(data)):
                for iloop in xrange(oloop, len(data)):
                    tmpVar = _cov(data[oloop], data[iloop])
                    covMat[iloop][oloop] = tmpVar
                    covMat[oloop][iloop] = tmpVar
            return covMat

try:
    from scipy.stats.distributions import norm as norm_dist
    sn_cdf = norm_dist(loc=0, scale=1).cdf
except ImportError:
    def sn_cdf(x):
        # brute force integration of the sn pdf
        # probably accumualtes errors and has all sorts of other
        # undesirable properties
        # however, it is still accurate to 1e-8 for the worst case
        #
        # this could be much faster/more accurate if I weighted
        # the points better

        # if we're more than 20 sd over the mean, return 0.0
        if x > 20: return 0.0

        # if we're more than 20 sd less than the mean, return 0
        if x < -20: return 1.0

        def norm_pdf(x, u, s):
            from math import pi, sqrt, exp

            constTerm = 1/(s*sqrt(2*pi))
            expTerm = exp(-0.5*(s**-2)*((x-u)**2))

            return constTerm*expTerm

        # if we are past 20, the added value to the cdf
        # is unsubstantial: so just integrate to 20
        a = max(-20, float(x))
        b = 20.0

        estimate = 0
        num = 10000
        for loop in xrange(num):
            ln = a + loop*((b-a)/num)
            hn = a + (loop+1)*((b-a)/num)
            local_est = (hn-ln)*norm_pdf((hn+ln)/2.0, 0, 1)
            estimate += local_est

        # this integrates above x - so reverse it
        estimate = 1-estimate

        ## some code to test my approximation
        #from scipy.stats.distributions import norm as norm_dist
        #sci_sn_cdf = norm_dist(loc=0, scale=1).cdf
        #print estimate, sci_sn_cdf(x), sci_sn_cdf(x) - estimate

        return estimate

##################### End Standard Statistical Methods #########################

