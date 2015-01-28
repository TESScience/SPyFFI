'''Define, test, and compare potential cosmic ray mitigation strategies.'''
import matplotlib as mpl
mpl.rcParams['axes.labelsize'] = 8#'small'
mpl.rcParams['xtick.labelsize'] = 8#'small'
mpl.rcParams['ytick.labelsize'] = 8#'small'
mpl.rcParams['font.size'] = 8#'small'

#mpl.rcParams['mathtext.fontset'] = 'custom'
# mpl.rcParams['legend.handlelength'] = 4
#mpl.rcParams['mathtext.sf'] = 'sans'
#mpl.rcParams['mathtext.default'] = 'sf'
from imports import *
import transit
import scipy.stats
plt.ion()
#plt.rc('text', usetex=False)
#plt.rc('font', family='serif')



def drawCosmics(n):



	def weighted_values(values, probabilities, size):
	    bins = np.add.accumulate(probabilities)
	    return values[np.digitize(np.random.random_sample(size), bins)]


	x, y = np.load('/Users/zkbt/Cosmos/Data/TESS/CR/cosmics_histogram.npy')
	y /= np.sum(y)

	return weighted_values(x,y,n)

class timeseries():
	'''Object to store a light curve, both unbinned and binned.'''
	def __init__(self, cube=None, pixel=(0,0), nexposures=1324, nsubexposures=900, amplitude=None, subexposurecadence=2.0, mag=10.0):
		'''Initialize timeseries object, either as a toy model or from a pixel drawn from a simulated image cube.

		(for toy model):
		t = timeseries(nexposures=[how many binned exposures?], nsubexposures=[how many subexposures within each exposure?], amplitude=[what's the amplitude of cosmic rays relative to noise?])

		(for cube simulation):
		t = timeseries(cube=[an (pix)x(pix)x(time) image cube object], pixel=[tuple indicated which pixel to use, like '(4,3)'])'''

		# the number of subexposures is defined first of all
		self.nsubexposures = nsubexposures

		# either make a toy light curve based on input parameters
		if cube is None:
		 	self.nexposures = nexposures
			if amplitude is None:
				self.mag = mag
				self.containedflux = 0.5
				self.photonsfromstar = 1.7e6*2*73.0*10**(-0.4*mag)*self.containedflux
				skynoise = 10.0
				readoutnoise = 10.0
				signal = self.nsubexposures*self.photonsfromstar
				noise = np.sqrt(self.nsubexposures*(readoutnoise**2 + skynoise**2 + self.photonsfromstar))
				snr = signal/noise
				self.exposurenoise = 1.0/snr
				self.subexposurenoise = self.exposurenoise*np.sqrt(self.nsubexposures)
				self.cosmicsamplitude = 1000.0/self.photonsfromstar/self.nsubexposures#amplitude*self.exposurenoise
				self.scale =  self.timeseries.nsubexposures

			else:
				self.exposurenoise = 1.0/self.nsubexposures
				self.subexposurenoise = self.exposurenoise*np.sqrt(self.nsubexposures)
				self.cosmicsamplitude = amplitude*self.exposurenoise
				self.scale=1

			self.cosmicsperexposure = 1.0 - 0.999**self.nsubexposures
			self.cosmicspersubexposure = self.cosmicsperexposure/self.nsubexposures
			self.subexposurecadence=subexposurecadence
			self.exposurecadence = self.nsubexposures*self.subexposurecadence
			self.createSimple()
		else:
			self.createFromCube(cube, pixel)

	@property
	def shape(self):
		'''The shape of arrays in this timeseries.'''
		return (self.nexposures, self.nsubexposures)

	def createFromCube(self, cube, pixel):
		'''Populate this timeseries, using a pixel from a cube.'''

		self.nexposures = (np.int(cube.n/self.nsubexposures))
		x = pixel[0]
		y = pixel[1]
		# KLUDGE (possibly?) -- make sure that reshape is working the way we think it is
		self.flux = cube.photons[x,y,:self.nexposures*self.nsubexposures].reshape(self.shape)
		self.noiseless = cube.noiseless[x,y,:self.nexposures*self.nsubexposures].reshape(self.shape)
		self.cosmics = cube.cosmics[x,y,:self.nexposures*self.nsubexposures].reshape(self.shape)
		self.subexposurenoise = cube.sigma()[x,y]
		self.exposurenoise = self.subexposurenoise/np.sqrt(self.nsubexposures)
		iscosmic = cube.cosmics.flatten() > 0
		self.cosmicsamplitude = np.mean(cube.cosmics.flatten()[iscosmic])/self.nsubexposures/self.exposurenoise
		self.cosmicspersubexposure = np.sum(iscosmic).astype(np.float)/len(cube.cosmics.flatten())
		self.cosmicsperexposure = self.cosmicspersubexposure*self.nsubexposures

		self.x = np.arange(0, self.nexposures, 1.0/self.nsubexposures).reshape(self.shape)
		self.toy = False
		self.cube = cube
		self.pixel = pixel


	def createSimple(self, cosmics=True, noise=True):
		'''Populate this timeseries with a simple toy model.'''
		self.flux = np.ones(self.shape)
		if noise:
			self.addNoise()
		if cosmics:
			self.addCosmics()
		self.x = np.arange(0, self.nexposures, 1.0/self.nsubexposures).reshape(self.shape)
		self.toy = True

	def createTransit(self, cosmics=True, mag=10.0, period = 0.9351425135):

		# add in a transit

		self.createSimple(noise=False, cosmics=False)
		p = transit.Planet(period=period, t0=period/2.0)
		s = transit.Star()
		i = transit.Instrument()
		self.tm = transit.TM(planet=p, star=s, instrument=i)
		self.tlc = transit.TLC(self.x*self.exposurecadence/60.0/60.0/24.0, self.flux, self.subexposurenoise*np.ones_like(self.flux))
		self.tlc.linkModel(self.tm)
		self.flux = self.tm.model()

		# add in noise
		self.addNoise()
		if cosmics:
			self.addRealCosmics()


	def addNoise(self):
		'''For toy model, include Gaussian noise.'''
		self.noiseless = self.flux + 0.0
		self.flux += np.random.normal(0,self.subexposurenoise,self.shape)

	def addCosmics(self):
		'''For toy model, include cosmic rays as random impulses with fixed amplitude.'''
		self.cosmics = self.cosmicsamplitude*self.nsubexposures*np.random.poisson(self.cosmicspersubexposure, self.shape)
		self.flux += self.cosmics

	def addRealCosmics(self):
		self.cosmics = drawCosmics(len(self.flux.flatten())).reshape(self.flux.shape)/self.photonsfromstar
		self.flux += self.cosmics


	def plot(self):
		'''Simple plot of this timeseries.'''
		plt.figure('unbinned')
		plt.cla()
		try:
			x = self.x % self.tm.planet.period.value*60*60*24.0
		except:
			x = self.x
		plt.plot(x.flatten(), self.flux.flatten(), linewidth=0, marker='.', alpha=0.3)
		plt.draw()

	def __str__(self):
		'''Summary of this object.'''
		return "{nexposures} exposures, {nsubexposures} subexposures, cosmic rays {amplitude}X noise".format(nexposures=self.nexposures, nsubexposures=self.nsubexposures, amplitude=self.cosmicsamplitude/self.exposurenoise)

class strategy(object):
	'''Define structures and methods needed by all cosmic ray rejection strategies.'''
	def __init__(self, n=None):
		'''Initialize strategy object.'''

		# setup basic filtering parameters:
		self.name = 'Nothing'
		if n is None:
			self.n = 5
		else:
			self.n = n

		self.prefix = ''

		# define a dictionary of plotting parameters
		self.plotting = {'nsigma':10, 'toplot':96, 'flux':'mediumvioletred', 'nocosmics':'orange', 'naive':'blue'}

	def directory(self):
		'''Define directory in which to store plots and files related to this strategy.'''
		try:
			self.workingDirectory = self.timeseries.cube.directory
		except:
			self.workingDirectory = '/Users/zkbt/Cosmos/Data/TESS/CR/'


		dir = self.workingDirectory + self.name.replace(' ', '') + '/'
		zachopy.utils.mkdir(dir)

		try:
			dir = dir + '{0:.0f}_{1:.0f}/'.format(self.timeseries.pixel[0], self.timeseries.pixel[1])
			zachopy.utils.mkdir(dir)
			print dir
		except:
			pass
		return dir

	def fileprefix(self):

		f = self.directory() + 'crdemo_filter{0}_n{1}_cr{2}_{3}exposures'.format(self.name.replace(' ',''), self.n, self.timeseries.cosmicsamplitude/self.timeseries.exposurenoise, self.timeseries.nexposures)
		return f.replace('.', 'p')

	def strategyprefix(self):
		f = self.directory() + 'crdemo_filter{0}_n{1}'.format(self.name.replace(' ',''), self.n, self.timeseries.nexposures)
		return f

	def calculate(self, timeseries):
		'''Populate unbinned and binned arrays, using the filter defined for this strategy.'''

		# setup the binned and unbinned arrays
		self.timeseries = timeseries
		self.prefix = ''

		# set up the unbinned dictionary
		self.unbinned = {'x':self.timeseries.x.flatten(), 'flux':self.timeseries.flux.flatten(), 'nocosmics':(self.timeseries.flux - self.timeseries.cosmics).flatten()}

		# bin, using the strategy
		self.binned = self.combine()

		# calculate standard deviations of the binned light curves, with and without the cosmic ray filter
		self.unmititigated = np.std(self.binned['naive'])
		self.achieved = np.std(self.binned['flux'])

		self.prefix += '   '


	def mean(self):
		'''Collapse over the subexposures, using a simple mean.'''
		return np.mean(self.timeseries.flux, 1)

	def noiselessmean(self):
		'''Collapse over the subexposures, using a simple mean.'''
		return np.mean(self.timeseries.noiseless, 1)

	def filter(self):
		raise Exception("Uh-oh! It seems no filtering function was defined!")

	def combine(self):
		'''Bin the unbinned timeseries, using the appropriate filter.'''

		# define x coordinate for the binned timeseries
		x = np.arange(0.5, self.timeseries.nexposures +0.5, 1)

		# naive timeseries is binned using a simple mean
		naive = self.mean()

		# noirseless binned timeseries
		noiseless = self.noiselessmean()

		# flux timeseries is binned using this strategy's filter
		flux = self.filter()

		# nocosmics timeseries is binned using a simple mean, but assuming the cosmics were never injected
		nocosmics = np.mean(self.timeseries.flux - self.timeseries.cosmics, 1)

		# return the binned timeseries
		self.binned = {'x':x, 'flux':flux, 'nocosmics':nocosmics, 'naive':naive, 'noiseless':noiseless}
		return self.binned

	def plottimeseries(self, x,y, **kwargs):
		'''Plot a timeseries.'''

		# plot only a subset (for quickness and clarity's sake)
		ok = x <= self.plotting['toplot']

		# plot in the timeseries panel
		ax = self.ax['timeseries']
		ax.plot(x[ok], y[ok], **kwargs)
		ax.set_xlim(0, self.plotting['toplot'])

	def plothistogram(self, y,binwidth=0.1, ax=None, fixed=False, **kwargs):
		'''Plot a histogram, rotated clockwise by 90 degrees, to represent a projection of a timeseries plot.'''

		# set the binwidth
		self.binwidth = binwidth

		# create a histogram of the lightcurve values
		yhist, edges = np.histogram(y, bins=np.arange(np.min(y)-self.binwidth, np.max(y)+self.binwidth, self.binwidth), density=True)


		# define the "x"-axis at the centers of the histogram bins
		xhist = (edges[1:] + edges[0:-1])/2.0

		# plot in the histogram panel
		if ax is None:
			ax = self.ax['histbinned']
		ax.plot(yhist, xhist, **kwargs)

		ax.set_xscale('log')
		ax.set_xlim(np.min(yhist)*0.5, np.max(yhist)*1.5)


	def plotimage(self, ax=plt.gca()):
		'''Display the image from which this light curve came, and indicate which pixel.'''

		# make sure the light curve was actually derived from an image cube
		assert(self.timeseries.toy == False)

		# pull out the median image for plotting
		image = self.timeseries.cube.median()

		ax.imshow(np.log(np.transpose(image)), interpolation='nearest', cmap='gray_r')
		ax.plot(self.timeseries.pixel[0], self.timeseries.pixel[1],marker='o',markerfacecolor='None',markersize=10, alpha=0.5, markeredgecolor='gray')
		ax.text(self.timeseries.pixel[0], self.timeseries.pixel[1], '     {0}'.format(self.timeseries.pixel), horizontalalignment='left', verticalalignment='center', color='gray', alpha=0.5, weight='bold', size=4)


	def test(self, t, remake=False, niterations=20):
		'''Test strategy over a range of input cosmic ray amplitudes (using toy model light curves).'''
		self.timeseries = t
		print "testing {0}, with filename of {1}".format(self.name, self.strategyprefix())
		filename = self.strategyprefix() + '_{0}iterations_{1:03.0f}subexposures.npy'.format(niterations, self.timeseries.nsubexposures)
		try:
			# can we load a preprocessed file?
			assert(remake == False)
			print 'trying to load file from ', filename
			self.amplitudes, self.noise = np.load(filename)
			print self.strategyprefix() + " loaded from file"
		except:

			# create a grid of relevant cosmic ray amplitudes
			self.amplitudes = np.arange(0, 5, 0.5)
			self.noise = np.zeros((len(self.amplitudes), niterations))

			# loop over iterations, to improve precision of the noise estimates
			for iteration in np.arange(niterations):
				print "  iteration #{0}/{1}".format(iteration, niterations)
				# loop over amplitudes
				for i in range(len(self.amplitudes)):


					# create a new timeseries
					self.timeseries = timeseries(nexposures=t.nexposures, nsubexposures=t.nsubexposures, amplitude=self.amplitudes[i])

					# bin it, using the strategy
					self.calculate(self.timeseries)

					# create a demo plot of this light curve
					if iteration == 0:
						self.plot()

					# print this timeseries's summary
					print self.prefix + "{0}".format(self.timeseries)

					# store the achieve noise
					self.noise[i, iteration] = self.achieved/self.timeseries.exposurenoise

			self.noise = np.mean(self.noise, 1)

			# save this calculation, so we can use again if need be
			np.save(filename, (self.amplitudes, self.noise))
			print '   saved results to ', filename

		# return (x, y) for plotting
		return self.amplitudes, self.noise

	def plot(self,unbinned=False,**kwargs):
		'''Create a pretty plot comparing the binned timeseries and their histograms.'''


		# setup the plots
		includeimage = self.timeseries.toy == False
		self.figure = plt.figure(0, figsize=(9 + 3*includeimage,2), dpi=150)
		self.gs = plt.matplotlib.gridspec.GridSpec(1,6+includeimage,wspace=0,hspace=0,left=0.08, right=0.98, top=0.85, bottom=0.2)
		self.figure = plt.gcf()
		self.figure.clf()
		self.figure.suptitle("Cosmic Ray Rejection with [{0}]".format(self.name))
 		self.ax = {}
		self.ax['timeseries'] =	plt.subplot(self.gs[0,:-(1 + includeimage)])
		self.ax['timeseries'].set_autoscale_on(False)
		self.ax['timeseries'].set_xlabel(r'Exposure #')
		self.ax['timeseries'].set_ylabel(r'Flux')
		self.ax['histbinned'] = plt.subplot(self.gs[0,-(1 + includeimage)], sharey=self.ax['timeseries'])
		plt.setp(self.ax['histbinned'].get_yticklabels(), visible=False)
		if includeimage:
			self.ax['image'] = plt.subplot(self.gs[0,-1])
			plt.setp(self.ax['image'].get_xticklabels(), visible=False)
			plt.setp(self.ax['image'].get_yticklabels(), visible=False)




		# set up the zero-level for the plot, the width of the plotting windows, and the position of the text labels
		ylevel = np.median(self.timeseries.scale*self.binned['nocosmics'])
		ywidth = self.timeseries.scale*self.plotting['nsigma']*0.8*self.timeseries.exposurenoise
		y =  ylevel - ywidth


		# fiddle with the ylimits of the plots, depending on whether looking at binned or unbinned timeseries
		if unbinned:
			self.ax['histbinned'].set_ylim(np.min(self.timeseries.scale*self.unbinned['flux']), np.max(self.timeseries.scale*self.unbinned['flux']))
			#self.ax['histbinned'].set_ylim(-self.plotting['nsigma']*self.timeseries.subexposurenoise, self.plotting['nsigma']*self.timeseries.subexposurenoise)
		else:
			self.ax['histbinned'].set_ylim(-self.timeseries.scale*self.plotting['nsigma']*self.timeseries.exposurenoise + ylevel, self.timeseries.scale*self.plotting['nsigma']*self.timeseries.exposurenoise + ylevel)

		# plot the unbinned timeseries, if desired
		if unbinned:
			self.plottimeseries(self.unbinned['x'], self.timeseries.scale*self.unbinned['flux'], linewidth=1, alpha=0.1, color='black')# markersize=2, markeredgewidth=0, markerfacecolor='black', marker='o', linewidth=0, alpha=0.1)
			#self.plottimeseries(self.unbinned['x'], self.unbinned['nocosmics'], linewidth=1, alpha=0.1, color='orange')# markersize=2, markeredgewidth=0, markerfacecolor='black', marker='o', linewidth=0, alpha=0.1)

		# plot the cosmic ray hits as vertical lines
		cosmichit = ((self.unbinned['flux'] -  self.unbinned['nocosmics']) > 0)*(self.unbinned['x']<=self.plotting['toplot'])
		for x in self.unbinned['x'][cosmichit]:
			self.ax['timeseries'].axvline(x, color='black', alpha=0.1, linewidth=0.5)


		# plot the binned timeseries
		self.plottimeseries(self.binned['x'], self.timeseries.scale*self.binned['nocosmics'], markersize=4, marker='o', alpha=0.5, markerfacecolor='orange', color=self.plotting['nocosmics'], markeredgewidth=0, linewidth=1)
		self.plottimeseries(self.binned['x'], self.timeseries.scale*self.binned['naive'], markersize=4, marker='o', alpha=0.5, markerfacecolor=self.plotting['naive'], color=self.plotting['naive'], markeredgewidth=0, linewidth=1)
		self.plottimeseries(self.binned['x'], self.timeseries.scale*self.binned['flux'], markersize=4, marker='o', alpha=0.5, markerfacecolor=self.plotting['flux'], color=self.plotting['flux'], markeredgewidth=0, linewidth=1)

		# plot the binned histogram
		self.plothistogram(self.timeseries.scale*self.binned['nocosmics'], binwidth=self.timeseries.scale*0.2*self.timeseries.exposurenoise, alpha=0.5, color=self.plotting['nocosmics'])
		self.plothistogram(self.timeseries.scale*self.binned['naive'], binwidth=self.timeseries.scale*0.2*self.timeseries.exposurenoise, alpha=0.5, color=self.plotting['naive'])
		self.plothistogram(self.timeseries.scale*self.binned['flux'], binwidth=self.timeseries.scale*0.2*self.timeseries.exposurenoise, alpha=0.5, color=self.plotting['flux'])

		# plot the image and position, if this isn't just a toy model
		if includeimage:
			self.plotimage(ax=self.ax['image'])

		# labels describing the simulated timeseries
		self.ax['timeseries'].text(self.plotting['toplot']*0.02, y, "{rate:0.2f} cosmic rays per exposure, at {amplitude:0.1f}X noise".format( amplitude=self.timeseries.cosmicsamplitude/self.timeseries.exposurenoise, rate=self.timeseries.cosmicsperexposure), fontsize=6, color=self.plotting['flux'], alpha=0.7, horizontalalignment='left')
		self.ax['timeseries'].text(self.plotting['toplot']*0.98, y, "{nsubexposures} subexposures".format( nsubexposures=self.timeseries.nsubexposures), fontsize=6, color=self.plotting['flux'], alpha=0.7, horizontalalignment='right')

		# labels showing the unmitigated, achieved, and ideal noise values
		left, right = np.log(self.ax['histbinned'].get_xlim())
		span = right - left
		self.ax['histbinned'].text(np.exp(span*0.2 + left), y, "{1:.2f}".format(self.timeseries.scale*self.unmititigated, self.unmititigated/self.timeseries.exposurenoise), fontsize=6, color=self.plotting['naive'], horizontalalignment='center', alpha=0.7)
		self.ax['histbinned'].text(np.exp(span*0.5 + left), y, "{1:.2f}".format(self.timeseries.scale*self.achieved, self.achieved/self.timeseries.exposurenoise), fontsize=6, color=self.plotting['flux'], horizontalalignment='center', alpha=0.7)
		self.ax['histbinned'].text(np.exp(span*0.78 + left), y, "{1:.2f}".format(self.timeseries.scale*self.timeseries.exposurenoise, 1.0), fontsize=6, color=self.plotting['nocosmics'], horizontalalignment='center', alpha=0.7)
		self.ax['histbinned'].text(np.exp(span*0.2 + left), ylevel - ywidth*1.1, 'unmitigated', fontsize=4, color=self.plotting['naive'], horizontalalignment='center', alpha=0.7)
		self.ax['histbinned'].text(np.exp(span*0.5 + left), ylevel - ywidth*1.1, 'achieved', fontsize=4, color=self.plotting['flux'], horizontalalignment='center', alpha=0.7)
		self.ax['histbinned'].text(np.exp(span*0.8 + left), ylevel - ywidth*1.1, 'perfect', fontsize=4, color=self.plotting['nocosmics'], horizontalalignment='center', alpha=0.7)

		# draw and save the plot
		plt.draw()
		plt.savefig(self.fileprefix() + '.pdf')
		print "saved plot to " +  self.fileprefix() + '.pdf'


class mean(strategy):
	'''Binning strategy -- the simplest, taking the straight sum of all the subexposures.'''
	def __init__(self):
		strategy.__init__(self)
		self.name = "Mean"

	def filter(self):
		return np.mean(self.timeseries.flux, 1)


class median(strategy):
	'''Binning strategy -- break into non-overlapping subsets, take median of each, sum those medians.'''
	def __init__(self, n=None):
		strategy.__init__(self, n)
		self.name = "Median of {0}".format(self.n)

	def filter(self):
		shape = self.timeseries.flux.shape
		reshapen = self.timeseries.flux.reshape(shape[0], shape[1]/self.n, self.n)
		partial = np.median(reshapen, 2)
		return np.mean(partial, 1)

class shiftedmedian(strategy):
	'''Binning strategy -- break into overlapping subsets, take median of each, sum those medians.'''
	def __init__(self, n=None):
		strategy.__init__(self, n)
		self.name = "Shifted Median of {0}".format(self.n)

	def filter(self):
		shape = self.timeseries.flux.shape
		sumofmedians = np.zeros((shape[0], shape[1]/self.n))
		for i in range(self.n):
			# note, this wraps around to the start at the end of the array; but should work for testing
			reshapen = np.roll(self.timeseries.flux, i, 1).reshape(shape[0], shape[1]/self.n, self.n)
			sumofmedians += np.median(reshapen, 2)
			#print reshapen[0,0,:]
			#print reshapen[0,1,:]
			#print
		return np.mean(sumofmedians/self.n, 1)


class sullivan(strategy):
	'''Binning strategy -- idea Peter and I played around with; didn't work particularly well.'''

	def __init__(self, n=None):

		n = 3
		strategy.__init__(self, n)
		self.name = "Sullivan of {0}".format(self.n)

	def filter(self):
		shape = self.timeseries.flux.shape
		reshapen = self.timeseries.flux.reshape(shape[0], shape[1]/self.n, self.n)
		median = np.median(reshapen, 2)
		d = np.sum(reshapen, 2)/3.
		a = (3*d - reshapen[:,:,0])/2.
		b = (3*d - reshapen[:,:,1])/2.
		c = (3*d - reshapen[:,:,2])/2.

		diffs = np.array([np.abs(a-median), np.abs(b - median), np.abs(c -median), np.abs(d - median)])
		mask = diffs == np.min(diffs, 0).reshape(1,diffs.shape[1], diffs.shape[2])*np.ones_like(diffs)

		print
		print np.mean(mask[0,:,:])
		print np.mean(mask[1,:,:])
		print np.mean(mask[2,:,:])
		print np.mean(mask[3,:,:])

		values = np.array([a,b,c,d])
		partial = np.sum(mask*values,0)/np.sum(mask, 0)
		return np.mean(partial, 1)

class lowest(strategy):
	'''Binning strategy -- break into subsets, reject the highest point from each and take the mean of the rest, sum these truncated means.'''
	def __init__(self, n=None):
		strategy.__init__(self, n)
		m = None
		if m is None:
			m = self.n - 1
		self.m = m
		self.name = "Lowest {m} out of {n}".format(m = self.m, n=self.n)

	def filter(self):
		assert(self.m == self.n-1)
		shape = self.timeseries.flux.shape
		reshapen = self.timeseries.flux.reshape(shape[0], shape[1]/self.n, self.n)
		sum = np.sum(reshapen, 2)
		max = np.max(reshapen, 2)
		corrected = np.mean((sum - max)/self.m,1)
		return corrected/ np.mean(corrected)	# this isn't correct -- it should be an empirical correction factor!

class central(strategy):
	'''Binning strategy -- break into subsets, reject the highest and lowest points from each and take the mean of the rest, sum these truncated means.'''
	def __init__(self, n=None):
		strategy.__init__(self, n)
		m = None
		if m is None:
			m = self.n - 2
		self.m = m
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

class outlierwithdecay(strategy):
	'''Binning strategy -- break into subsets, estimate standard deviation with weighted mean of current subset and previous estimate, reject from each and take the mean of the rest, sum these truncated means.'''
	def __init__(self, n=10, threshold=10.0, memory=0.90, safetybuffer=2.0, diagnostics=False):
		'''Initialize an outlierwith decay strategy.

		n = the number of subexposures in each "chunk"
		threshold = how many sigma about the noise are required for a point to be an outlier?
		memory = the weight given to previous estimates of standard deviations (best estimate = memory[previous best estimate] + (1 - memory)[most recent chunk])
		safetybuffer = by what factor should we inflate the initial standard deviation to prevent overfitting?'''

		# initialize the basic strategy class
		strategy.__init__(self, n)

		# store the parameters of the filter
		self.threshold = threshold
		self.memory = memory
		self.safetybuffer = 2.0

		# define a name for this filter
		self.name = r'Rejecting {threshold}$\sigma$ Outliers; $\sigma$ from chunks of {n};  memory of {memory}%'.format(threshold=self.threshold, memory=self.memory*100, n=self.n)

		# for testing, keep a diagnostics flag to say whether to display the mean + std. estimates
		self.diagnostics = diagnostics

	def directory(self):
		'''Define directory in which to store plots and files related to this strategy.'''
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
		'''Custom strategy prefix for this long named strategy.'''
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


class hybrid(strategy):
	def __init__(self, n=None, threshold=7):
		strategy.__init__(self, n)
		m = None
		if m is None:
			m = self.n - 1
		self.m = m
		self.threshold = threshold
		self.name = r"Rejecting {threshold}$\sigma$ Outliers; at most {nminusm} out of {n}".format(threshold=self.threshold, nminusm = self.n - self.m, n=self.n)

	def filter(self):
		assert(self.m == self.n-1)
		shape = self.timeseries.flux.shape
		reshapen = self.timeseries.flux.reshape(shape[0], shape[1]/self.n, self.n)
		max = np.max(reshapen, 2)
		isntmax = reshapen < max.reshape(shape[0], shape[1]/self.n, 1)*np.ones_like(reshapen)
		nisntmax = np.sum(isntmax,2)
		nismax = (self.n - nisntmax)
		mean = np.sum(reshapen*isntmax, 2)/nisntmax
		meanofsquares = np.sum(reshapen**2*isntmax, 2)/nisntmax
		variance = meanofsquares
		maxisbad = (max - mean)**2/variance > self.threshold
		return np.mean((mean*(nisntmax + nismax*maxisbad) + nismax*(maxisbad == False)*max)/self.n, 1)


def plottransit(s):
	plt.figure('transit',figsize=(9,6), dpi=150)


	plt.clf()
	gs = plt.matplotlib.gridspec.GridSpec(2,2,height_ratios=[1, 0.3], width_ratios=[1,0.2], wspace=0, hspace=0)
	axtransit = plt.subplot(gs[0,0])
	axresiduals = plt.subplot(gs[1,0], sharex=axtransit)
	axhistogram = plt.subplot(gs[1,1], sharey=axresiduals)
	plt.setp(axtransit.get_xticklabels(), visible=False)

	period= s.timeseries.tm.planet.period.value
	t0 = s.timeseries.tm.planet.t0.value
	#
	#t_unbinned = s.unbinned['x']/24.0/60.0/60.0*s.timeseries.exposurecadence
	#flux_unbinned = s.unbinned['flux']
	#ax.plot((t_unbinned % period) - t0, flux_unbinned, linewidth=0, marker='o', color='gray', alpha=0.25, markersize=1)

	t_binned = s.binned['x']/24.0/60.0/60.0*s.timeseries.exposurecadence
	flux_binned = s.binned['flux']
	phasedtime = (t_binned % period) - t0
	sort = np.argsort(phasedtime)
	alphas = dict(naive=0.5, nocosmics=0.5, flux=0.5)
	axhistogram.set_autoscale_on(True)
	noises = {}
	positions = dict(naive = 0.2, flux = 0.5, nocosmics=0.8)
	for tag in ['naive','nocosmics', 'flux']:
		kw = dict(color=s.plotting[tag], markerfacecolor=s.plotting[tag],  markersize=2, marker='o', alpha=alphas[tag]*0.1, markeredgewidth=0, linewidth=0)
		axtransit.plot(phasedtime[sort], s.binned[tag][sort], **kw)
		binwidth = s.timeseries.nsubexposures*2.0/24.0/60.0/60.0/2.0
		bx, by, be = zachopy.oned.binto(phasedtime[sort], s.binned[tag][sort], binwidth)
		axtransit.plot(bx, by, linewidth=2, color=s.plotting[tag], alpha=alphas[tag])
		residuals = s.binned[tag][sort]-s.binned['noiseless'][sort]
		axresiduals.plot(phasedtime[sort], residuals, **kw)
		bx, by, be = zachopy.oned.binto(phasedtime[sort], residuals, binwidth)
		axresiduals.plot(bx, by, linewidth=2, color=s.plotting[tag], alpha=alphas[tag])
		noises[tag] = np.std(residuals)


		s.plothistogram(residuals, ax=axhistogram, binwidth=0.2*s.timeseries.exposurenoise, alpha=0.5, linewidth=2, color=s.plotting[tag])
		left, right = np.log(axhistogram.get_xlim())
		span = right - left
		ylevel = 0

		ywidth = s.timeseries.exposurenoise*10
		axhistogram.text(np.exp(span*positions[tag] + left), ylevel + ywidth*1.2, "{0:.0f}e-\n({1:.2f})\n{2:.0f}".format(noises[tag]*s.timeseries.photonsfromstar*s.timeseries.nsubexposures, noises[tag]/s.timeseries.exposurenoise, noises[tag]*1e6), fontsize=6, color=s.plotting[tag], horizontalalignment='center', alpha=0.7)
		axhistogram.text(np.exp(span*positions[tag] + left), ylevel + ywidth, tag, fontsize=4, color=s.plotting[tag], horizontalalignment='center', alpha=0.7)


	axresiduals.set_xlabel('Time from Mid-transit (days)')
	axresiduals.set_ylabel('O-C')
	axtransit.set_ylabel('Relative Flux')
	axtransit.set_title('{0}\non a {1:.1f} magnitude star with {2:.0f}% contained flux'.format(s.name, s.timeseries.mag, s.timeseries.containedflux*100))


	axresiduals.set_xlim(-0.05*period, 0.05*period)
	plt.setp(axhistogram.get_yticklabels(), visible=False)

	plt.draw()

def compare(nsubexposures=900, niterations=1, remake=True):

	nexposures = 1324*900/nsubexposures
	t = timeseries(nexposures=nexposures, nsubexposures=nsubexposures, amplitude=1.0)

	subsets = [[mean()],
				[mean(),outlierwithdecay(n=10,threshold=10, memory=0.1),outlierwithdecay(n=10,threshold=10, memory=0.5),outlierwithdecay(n=10,threshold=10, memory=0.9)],
				[mean(), median(n=3), median(n=4),  median( n=5), median( n=6)],
				[mean(),shiftedmedian(n=3), shiftedmedian(n=4), shiftedmedian(n=5), shiftedmedian(n=6)],
				[mean(),lowest(2), lowest(3), lowest(4), lowest(5),  lowest(10), lowest(20), lowest( 60)],
				[mean(),central(3), central(5), central(10),central(20), central( 60)],
				[mean(),  median(n=3), shiftedmedian(n=3), lowest(10),  central(10), outlierwithdecay(n=10,threshold=10, memory=0.5)]]
	labels = ['mean', 'outlierrejection', 'median', 'shiftedmedian', 'lowest', 'central',  'everything']
	for i in range(len(subsets)):
		strategies = subsets[i]
		names = []
		results = {}
		print labels[i]
		while(len(strategies)>0):
			print "  {0:2} strategies to go!".format(len(strategies))
			s = strategies.pop()
			s.test(t, niterations=niterations, remake=remake)
			plt.draw()
			names.append(s.name)
			results[s.name] = {'amplitudes':s.amplitudes, 'noise':s.noise, 'strategy':s}
			del s

		figure = plt.figure(figsize=(7,4), dpi=150)
		gs = plt.matplotlib.gridspec.GridSpec(1,4, wspace=0,hspace=0)

		ax = plt.subplot(gs[0:])
		plt.cla()
		medianlinestyle, lowestlinestyle, centrallinestyle, shiftedlinestyle, outlierstyle = 0,0,0, 0,0
		linestyles = ['-', '--']
		count =0
		for k in names:
			if "Median" in k and "Shifted" not in k:
				color = 'yellow'
				linestyle = medianlinestyle
				medianlinestyle += 1
			if  "Shifted" in k:
				color = 'blue'
				linestyle = shiftedlinestyle
				shiftedlinestyle  += 1

			if "Lowest" in k:
				color = 'orange'
				linestyle = lowestlinestyle
				lowestlinestyle += 1
			if "Central" in k:
				color = 'green'
				linestyle = centrallinestyle
				centrallinestyle += 1
			if "Reject" in k:
				color = 'red'
				linestyle = outlierstyle
				outlierstyle += 1
			print k
			n = np.float(results[k]['strategy'].n)
			scale = 2
			if "Mean" in k:
				color = 'black'
				linestyle=0
				linewidth=5
				alpha=0.2
			else:
				linewidth = (len(names) - linestyle + 1)/scale
				#linewidth = np.log(n/t.nsubexposures) - np.log(1.0/t.nsubexposures)
				alpha=1-(linewidth)/(len(names)+1.0)#(1 - np.log(1.0/t.nsubexposures))

			ax.plot(results[k]['amplitudes'], results[k]['noise'], label=k.replace(' out of ', '/').replace('; ', '\n '), color=color, linewidth=linewidth, alpha=alpha, linestyle=linestyles[linestyle % len(linestyles)])
			ax.legend(bbox_to_anchor=(0.95, 0.95), loc=1, borderaxespad=0., fontsize=7)
			ax.set_xlabel(r'Cosmic Ray Amplitude (in $\sigma$)')
			ax.set_ylabel(r'Achieved Noise (in $\sigma$)')
			ax.set_ylim(0.95, 1.25)
			ax.set_title('Cosmic Ray Rejection when Stacking {0} Exposures'.format(t.nsubexposures))
			count +=1
		plt.tight_layout()
		plt.draw()
		plt.savefig( '/Users/zkbt/Cosmos/Data/TESS/CR/comparisons/cosmicrayrejectioncomparisons_{2:03.0f}iterations_{1:.0f}subexposures_{0}.pdf'.format(labels[i], t.nsubexposures, niterations))


def testTransit(s=mean(), mag=10.0, short=True):

	if short:
		nexposures=1324*15
		nsubexposures=60
		period = 0.9351425135
	else:
		nexposures = 1324
		nsubexposures = 900
		period = 3.631246
	t = timeseries(nexposures=nexposures, nsubexposures=nsubexposures, mag=mag)
	t.createTransit(period=period)
	s.calculate(t)
	plottransit(s)
