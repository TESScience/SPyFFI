'''Tools to stack images, using a variety of algorithms
    (independent from Strategies, to optimize for 2D images).'''

from imports import *

class Stacker(Talker):
    '''Stack a cube of images, using some filter.'''

    def __init__(self, **kwargs):
        # decide whether or not this Stacker is chatty
        Talker.__init__(self, **kwargs)

class TruncatedMean(Stacker):
    '''Binning with TruncatedMean = break into subsets, reject the highest and lowest points from each and take the mean of the rest, sum these truncated means.'''
    def __init__(self, n=10, m=None):
        Stacker.__init__(self)
        self.n = n
        if m is None:
            self.m = self.n-2
        else:
            self.m = None
            assert(((self.n - self.m) % 2) == 0)
        self.name = "Central {m} out of {n}".format(m = self.m, n=self.n)

    def filter(self):
        assert(self.m == self.n-2)
        shape = self.timeseries.flux.shape
        reshapen = self.timeseries.flux.reshape(shape[0], shape[1]/self.n, self.n)
        sum = np.sum(reshapen, 2)
        max = np.max(reshapen, 2)
        min = np.min(reshapen, 2)
        corrected = np.mean((sum - max - min)/self.m,1)
        return corrected #- np.mean(corrected)

class outlierwithdecay(Strategy):
    '''Binning Strategy -- break into subsets, estimate standard deviation with weighted mean of current subset and previous estimate, reject from each and take the mean of the rest, sum these truncated means.'''
    def __init__(self, n=10, threshold=10.0, memory=0.90, safetybuffer=2.0, diagnostics=False):
        '''Initialize an outlierwith decay Strategy.

        n = the number of subexposures in each "chunk"
        threshold = how many sigma about the noise are required for a point to be an outlier?
        memory = the weight given to previous estimates of standard deviations (best estimate = memory[previous best estimate] + (1 - memory)[most recent chunk])
        safetybuffer = by what factor should we inflate the initial standard deviation to prevent overfitting?'''

        # initialize the basic Strategy class
        Strategy.__init__(self, n)

        # store the parameters of the filter
        self.threshold = threshold
        self.memory = memory
        self.safetybuffer = 2.0

        # define a name for this filter
        self.name = r'Rejecting {threshold}$\sigma$ Outliers; $\sigma$ from chunks of {n};  memory of {memory}%'.format(threshold=self.threshold, memory=self.memory*100, n=self.n)

        # for testing, keep a diagnostics flag to say whether to display the mean + std. estimates
        self.diagnostics = diagnostics

    def directory(self):
        '''Define directory in which to store plots and files related to this Strategy.'''
        try:
          self.workingDirectory = self.timeseries.cube.directory
        except:
          self.workingDirectory = '/Users/zkbt/Cosmos/Data/TESS/CR/'

        name = 'RejectionWith{0:.0f}percentMemory'.format(self.memory*100)

        dir = self.workingDirectory + name.replace(' ', '') + '/'
        zachopy.utils.mkdir(dir)

        try:
          dir = dir + '{0:.0f}_{1:.0f}/'.format(self.timeseries.pixel[0], self.timeseries.pixel[1])
          zachopy.utils.mkdir(dir)
          print dir
        except:
          pass
        return dir

    def strategyprefix(self):
        '''Custom Strategy prefix for this long named Strategy.'''
        name = 'RejectionWith{0:.0f}percentMemory'.format(self.memory*100)
        f = self.directory() + 'crdemo_filter{0}_n{1}'.format(name.replace(' ',''), self.n)
        return f

    def filter(self):
        '''Loop through the unbinned light curve, applying the realtime outlier-rejection filter.'''

        # make sure that the chunk size divides evenly into an exposure
        assert(self.timeseries.nsubexposures % self.n == 0)
        assert(self.n > 3)

        # define the number of chunks (and other convenience constants)
        nchunks = self.timeseries.nsubexposures/self.n
        nsubexposures = self.timeseries.nsubexposures
        nexposures = self.timeseries.nexposures
        n = self.n

        # create an array to store the final binned timeseries
        finaltimeseries = np.zeros(nexposures)

        # create arrays to store the per-chunk estimates of the mean and the standard deviation (these are arrays simply for plotting purposes)
        running_mean, running_std = np.zeros((nexposures, nchunks)), np.zeros((nexposures, nchunks))

        # initialize these estimates with the first chunk (rejecting no outliers)
        running_mean[0,0] = np.mean(self.timeseries.flux[0,0:n])
        running_std[0,0] = np.sqrt(np.sum((self.timeseries.flux[0,0:n] - running_mean[0,0])**2)/(n-1.0))

        # inflate the initial standard deviation measurement to prevent overfitting
        running_std[0,0] *= self.safetybuffer

        # set the first binned point to this mean estimate (there's no other choice)
        finaltimeseries[0] = running_mean[0,0]

        # loop over the binned exposures, and chunks within exposures
        count = 1
        for iexposure in np.arange(nexposures):
          for ichunk in np.arange(nchunks):

            # skip the very first point, because it's already been defined
            if (ichunk == 0)&(iexposure == 0):
              continue

            # pull out the light curve for this little chunk
            flux = self.timeseries.flux[iexposure, n*ichunk:n*(ichunk+1)]
            # pull out the last mean and standard deviation estimates
            best_mean = running_mean.flatten()[count-1]
            best_std = running_std.flatten()[count-1]

            # determine which points are not outliers, by being less than [threshold] sigma over the mean
            notoutlier = flux < (best_mean + self.threshold*best_std)
            assert(notoutlier.any())
            if notoutlier.any() == False:
              notoulier = np.ones_like(flux)
            # determine the mean and standard deviation of the good points in this chunk
            this_mean = np.mean(flux[notoutlier])
            this_std = np.sqrt(np.sum((flux[notoutlier] - this_mean)**2)/(n - 1.0))

            # store this binned exposure in the final array
            finaltimeseries[iexposure] += this_mean

            # mix this chunk into the running estimates, for the next chunk to use
            running_mean[iexposure, ichunk] = self.memory*best_mean + (1.0 - self.memory)*this_mean
            running_std[iexposure, ichunk] = self.memory*best_std + (1.0 - self.memory)*this_std
            # or would it be better to be doing this in variance space?, such as...
            # running_std[iexposure, ichunk] = np.sqrt(self.memory*best_std**2 + (1.0 - self.memory)*this_std**2)

            # advance the linear counter
            count += 1

        if self.diagnostics:

          # create a plot for showing the mean and standard deviation estimates
          plt.figure('diagnostics for {0}'.format(self.__class__), figsize=(7,4), dpi=150)
          plt.cla()

          # plot the two timeseries
          pargs = dict(linewidth=3, alpha=0.5)
          plt.plot(running_mean.flatten(), linestyle='--', label='mean', **pargs)
          plt.plot(running_std.flatten(), label='standard deviation', **pargs)
          plt.legend()

        # return the binned timeseries
        return finaltimeseries/nchunks
