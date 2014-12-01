'''Define and test potential cosmic ray mitigation strategies.'''
import matplotlib as mpl
mpl.rcParams['axes.labelsize'] = 8#'small'
mpl.rcParams['xtick.labelsize'] = 8#'small'
mpl.rcParams['ytick.labelsize'] = 8#'small'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['legend.handlelength'] = 4
#mpl.rcParams['mathtext.sf'] = 'sans'
#mpl.rcParams['mathtext.default'] = 'sf'

import zachopy.utils
import glob
import numpy as np
import matplotlib.pyplot as plt
import zachopy.star

plt.ion()
plt.rc('text', usetex=False)
plt.rc('font', family='serif')

workingDirectory = '/Users/zkbt/Cosmos/Data/TESS/CR/'


class timeseries():

	def __init__(self, nexposures=1000, nsubexposures=900, amplitude=2.0):
		self.nsubexposures = nsubexposures
	 	self.nexposures = nexposures
		self.shape = (self.nexposures, self.nsubexposures)
		self.exposurenoise = 1.0
		self.subexposurenoise = self.exposurenoise*np.sqrt(self.nsubexposures)
		self.cosmicsamplitude = amplitude*self.exposurenoise
		self.cosmicsperexposure = 0.5
		self.cosmicspersubexposure = self.cosmicsperexposure/self.nsubexposures
		self.createSimple()

	def createFromImages(self):
		pass

	def createSimple(self):
		self.flux = np.zeros(self.shape)
		self.addNoise()
		self.addCosmics()
		self.x = np.arange(0, self.nexposures, 1.0/self.nsubexposures).reshape(self.shape)

	def addNoise(self):
		self.flux += np.random.normal(0,self.subexposurenoise,self.shape)

	def addCosmics(self):
		self.cosmics = self.cosmicsamplitude*self.nsubexposures*np.random.poisson(self.cosmicspersubexposure, self.shape)
		self.flux += self.cosmics

	def plot(self):
		plt.cla()
		plt.plot(self.flux.flatten())
		plt.draw()

	def __str__(self):
		return "{nexposures} exposures, {nsubexposures} subexposures, cosmic rays {amplitude}X noise".format(nexposures=self.nexposures, nsubexposures=self.nsubexposures, amplitude=self.cosmicsamplitude)

class strategy(object):
	def __init__(self, t, n=None):



		# setup the timeseries
		self.timeseries=t

		# setup basic filtering parameters:
		self.name = 'Nothing'
		if n is None:
			self.n = self.timeseries.nsubexposures
		else:
			self.n = n

		self.prefix = ''

	def directory(self):

		dir = workingDirectory + self.name.replace(' ', '') + '/'
		zachopy.utils.mkdir(dir)
		return dir

	def fileprefix(self):
		f = self.directory() + 'crdemo_filter{0}_n{1}_cr{2}_{3}exposures'.format(self.name.replace(' ',''), self.n, self.timeseries.cosmicsamplitude, self.timeseries.nexposures)
		return f

	def strategyprefix(self):
		f = self.directory() + 'crdemo_filter{0}_n{1}_{2}exposures'.format(self.name.replace(' ',''), self.n, self.timeseries.nexposures)
		return f

	def calculate(self):
		# setup the binned and unbinned arrays
		self.prefix = ''
		self.unbinned = {'x':self.timeseries.x.flatten(), 'flux':self.timeseries.flux.flatten(), 'nocosmics':(self.timeseries.flux - self.timeseries.cosmics).flatten()}
		self.binned = self.combine()
		self.unmititigated = np.std(self.binned['naive'])
		self.achieved = np.std(self.binned['flux'])
		print self.prefix + self.name
		self.prefix += '   '


	def mean(self):
		return np.mean(self.timeseries.flux, 1)

	def filter(self):
		return self.mean()

	def combine(self):
		x = np.arange(0.5, self.timeseries.nexposures +0.5, 1)
		naive = self.mean()
		flux = self.filter()
		nocosmics = np.mean(self.timeseries.flux - self.timeseries.cosmics, 1)
		self.binned = {'x':x, 'flux':flux, 'nocosmics':nocosmics, 'naive':naive}
		return self.binned

	def plottimeseries(self, x,y, **kwargs):
		ok = x <= self.plotting['toplot']
		self.ax['timeseries'].plot(x[ok], y[ok], **kwargs)
		self.ax['timeseries'].set_xlim(0, self.plotting['toplot'])

	def plothistogram(self, y,binwidth=0.1,ax=plt.gca(),**kwargs):
		self.binwidth = binwidth
		yhist, edges = np.histogram(y, bins=np.arange(np.min(y)-self.binwidth, np.max(y)+self.binwidth, self.binwidth), density=True)
		xhist = (edges[1:] + edges[0:-1])/2.0
		ax.plot(yhist, xhist, **kwargs)
		ax.set_xscale('log')
		ax.set_xlim(np.min(yhist)*0.5, 1)

	def test(self):
		try:
			self.amplitudes, self.noise = np.load(self.strategyprefix() + '.npy')
		except:
			self.amplitudes = np.arange(0, 5, 0.5)
			self.noise = np.zeros_like(self.amplitudes)
			for i in range(len(self.amplitudes)):
				self.timeseries = timeseries(nexposures=self.timeseries.nexposures, nsubexposures=self.timeseries.nsubexposures, amplitude=self.amplitudes[i])
				#self.calculate()
				self.plot()
				print self.prefix + "{0}".format(self.timeseries)
				self.noise[i] = self.achieved
			np.save(self.strategyprefix() + '.npy', (self.amplitudes, self.noise))
		return self.amplitudes, self.noise

	def plot(self,unbinned=False,**kwargs):
		self.calculate()


		# setup the plots
		self.figure = plt.figure(0, figsize=(9,2), dpi=150)
		self.gs = plt.matplotlib.gridspec.GridSpec(1,6,wspace=0,hspace=0)
		self.figure = plt.gcf()
		self.figure.clf()
 		self.ax = {}
		self.ax['timeseries'] =	plt.subplot(self.gs[0,:-1])
		self.ax['timeseries'].set_autoscale_on(False)
		self.ax['timeseries'].set_xlabel(r'Exposure #')
		self.ax['timeseries'].set_ylabel(r'Flux/$\sigma$')

		self.ax['histbinned'] = plt.subplot(self.gs[0,-1], sharey=self.ax['timeseries'])
		plt.setp(self.ax['histbinned'].get_yticklabels(), visible=False)
		self.plotting = {'nsigma':10, 'toplot':42, 'flux':'mediumvioletred', 'nocosmics':'orange', 'naive':'blue'}
		if unbinned:
			self.ax['histbinned'].set_ylim(np.min(self.unbinned['flux']), np.max(self.unbinned['flux']))
			#self.ax['histbinned'].set_ylim(-self.plotting['nsigma']*self.timeseries.subexposurenoise, self.plotting['nsigma']*self.timeseries.subexposurenoise)

		else:
			self.ax['histbinned'].set_ylim(-self.plotting['nsigma']*self.timeseries.exposurenoise, self.plotting['nsigma']*self.timeseries.exposurenoise)
		self.figure.suptitle("Cosmic Ray Rejection with [{0}]".format(self.name))


		# plot the unbinned timeseries
		if unbinned:
			self.plottimeseries(self.unbinned['x'], self.unbinned['flux'], linewidth=1, alpha=0.1, color='black')# markersize=2, markeredgewidth=0, markerfacecolor='black', marker='o', linewidth=0, alpha=0.1)
			#self.plottimeseries(self.unbinned['x'], self.unbinned['nocosmics'], linewidth=1, alpha=0.1, color='orange')# markersize=2, markeredgewidth=0, markerfacecolor='black', marker='o', linewidth=0, alpha=0.1)
		cosmichit = ((self.unbinned['flux'] -  self.unbinned['nocosmics']) > 0)*(self.unbinned['x']<=self.plotting['toplot'])
		for x in self.unbinned['x'][cosmichit]:
			self.ax['timeseries'].axvline(x, color='black', alpha=0.1, linewidth=0.5)


		# plot the binned timeseries
		self.plottimeseries(self.binned['x'], self.binned['nocosmics'], markersize=4, marker='o', alpha=0.5, markerfacecolor='orange', color=self.plotting['nocosmics'], markeredgewidth=0, linewidth=1)
		self.plottimeseries(self.binned['x'], self.binned['naive'], markersize=4, marker='o', alpha=0.5, markerfacecolor=self.plotting['naive'], color=self.plotting['naive'], markeredgewidth=0, linewidth=1)
		self.plottimeseries(self.binned['x'], self.binned['flux'], markersize=4, marker='o', alpha=0.5, markerfacecolor=self.plotting['flux'], color=self.plotting['flux'], markeredgewidth=0, linewidth=1)

		# plot the binned histogram
		self.plothistogram(self.binned['nocosmics'], ax=self.ax['histbinned'], binwidth=0.2, alpha=0.5, color=self.plotting['nocosmics'])
		self.plothistogram(self.binned['naive'], ax=self.ax['histbinned'], binwidth=0.2, alpha=0.5, color=self.plotting['naive'])
		self.plothistogram(self.binned['flux'], ax=self.ax['histbinned'], binwidth=0.2, alpha=0.5, color=self.plotting['flux'])

		y =  -self.plotting['nsigma']*0.8
		self.ax['timeseries'].text(self.plotting['toplot']*0.02, y, "{rate:0.2f} cosmic rays per exposure, at {amplitude:0.1f}X noise".format( amplitude=self.timeseries.cosmicsamplitude, rate=self.timeseries.cosmicsperexposure), fontsize=6, color=self.plotting['flux'], alpha=0.7, horizontalalignment='left')
		self.ax['timeseries'].text(self.plotting['toplot']*0.98, y, "{nsubexposures} subexposures".format( nsubexposures=self.timeseries.nsubexposures), fontsize=6, color=self.plotting['flux'], alpha=0.7, horizontalalignment='right')

		left, right = np.log(self.ax['histbinned'].get_xlim())
		span = right - left
		self.ax['histbinned'].text(np.exp(span*0.25 + left), y, "{0:.2f}".format(self.unmititigated), fontsize=6, color=self.plotting['naive'], horizontalalignment='center', alpha=0.7)
		self.ax['histbinned'].text(np.exp(span*0.5 + left), y, "{0:.2f}".format(self.achieved), fontsize=6, color=self.plotting['flux'], horizontalalignment='center', alpha=0.7)
		self.ax['histbinned'].text(np.exp(span*0.75 + left), y, "{0:.2f}".format(self.timeseries.exposurenoise), fontsize=6, color=self.plotting['nocosmics'], horizontalalignment='center', alpha=0.7)

		self.ax['histbinned'].text(np.exp(span*0.25 + left), y*1.1, 'unmitigated', fontsize=4, color=self.plotting['naive'], horizontalalignment='center', alpha=0.7)
		self.ax['histbinned'].text(np.exp(span*0.5 + left), y*1.1, 'achieved', fontsize=4, color=self.plotting['flux'], horizontalalignment='center', alpha=0.7)
		self.ax['histbinned'].text(np.exp(span*0.75 + left), y*1.1, 'perfect', fontsize=4, color=self.plotting['nocosmics'], horizontalalignment='center', alpha=0.7)
		plt.draw()
		plt.tight_layout()
		plt.savefig(self.fileprefix() + '.pdf')


class mean(strategy):
	def __init__(self, t):
		strategy.__init__(self, t)
		self.name = "Mean"

	def filter(self):
		return np.mean(self.timeseries.flux, 1)


class median(strategy):
	def __init__(self, t, n=None):
		strategy.__init__(self, t, n)
		self.name = "Median of {0}".format(self.n)

	def filter(self):
		shape = self.timeseries.flux.shape
		reshapen = self.timeseries.flux.reshape(shape[0], shape[1]/self.n, self.n)
		partial = np.median(reshapen, 2)
		return np.mean(partial, 1)

class shiftedmedian(strategy):
	def __init__(self, t, n=None):
		strategy.__init__(self, t, n)
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
	def __init__(self, t, n=None):
		n = 3
		strategy.__init__(self, t, n)
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
	def __init__(self, t, n=None):
		strategy.__init__(self, t, n)
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
		return corrected - np.mean(corrected)

class central(strategy):
	def __init__(self, t, n=None):
		strategy.__init__(self, t, n)
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



class hybrid(strategy):
	def __init__(self, t, n=None, threshold=7):
		strategy.__init__(self, t, n)
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




def compare(nexposures=1000, t=None):
	if t is None:
		t = timeseries(nexposures=nexposures)

	subsets = [[mean(t)],
				[mean(t), median(t,n=3), median(t,n=4),  median(t, n=5), median(t, n=6)],
				[mean(t),shiftedmedian(t,n=3), shiftedmedian(t,n=4), shiftedmedian(t,n=5), shiftedmedian(t,n=6)],
				[mean(t),lowest(t,2), lowest(t,3), lowest(t,4), lowest(t,5),  lowest(t,10), lowest(t,20), lowest(t, 60),  lowest(t,900)],
				[mean(t),central(t,3), central(t,5), central(t,10),central(t,20), central(t, 60), central(t, 900)],
				[mean(t),  median(t,n=3), median(t,n=4), shiftedmedian(t,n=3), shiftedmedian(t,n=4),  lowest(t,10), lowest(t,20),  central(t,10),central(t,20) ]]
	labels = ['mean', 'median', 'shiftedmedian', 'lowest', 'central', 'everything']
	for i in range(len(subsets)):
		strategies = subsets[i]
		names = []
		results = {}
		print labels[i]
		while(len(strategies)>0):
			print "  {0:2} strategies to go!".format(len(strategies))
			s = strategies.pop()
			s.test()
			names.append(s.name)
			results[s.name] = {'amplitudes':s.amplitudes, 'noise':s.noise, 'strategy':s}
			del s

		figure = plt.figure(figsize=(7,4), dpi=150)
		gs = plt.matplotlib.gridspec.GridSpec(1,4, wspace=0,hspace=0)

		ax = plt.subplot(gs[0:3])
		plt.cla()
		medianlinestyle, lowestlinestyle, centrallinestyle, shiftedlinestyle = 0,0,0, 0
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

			ax.plot(results[k]['amplitudes'], results[k]['noise'], label=k.replace(' out of ', '/'), color=color, linewidth=linewidth, alpha=alpha, linestyle=linestyles[linestyle % len(linestyles)])
			ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=7)
			ax.set_xlabel(r'Cosmic Ray Amplitude (in $\sigma$)')
			ax.set_ylabel(r'Achieved Noise (in $\sigma$)')
			ax.set_ylim(1.0, 1.25)
			count +=1
		plt.tight_layout()
		plt.draw()
		plt.savefig(workingDirectory + 'cosmicrayrejectioncomparisons_{0}.pdf'.format(labels[i]))
