'''Generate pixel-lightcurve cubes with dimensions (xpix)x(ypix)x(time).'''
import zachopy.utils
import tess
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import zachopy.star

class cube(object):
	'''Cube to handle simulated postage stamp pixel light curves;
			has dimensions of (xpixels, ypixels, time).'''

	def __init__(self, name=None, size=32, n=900, cadence=2, jitter=False):
		'''Initialize a cube object.'''

		# define basic parameters
		self.size = size
		self.xpixels = self.size
		self.ypixels = self.size
		self.n = n
		self.cadence = cadence
		self.shape = (self.xpixels, self.ypixels, self.n)
		self.jitter = jitter
		# set up a tess camera object
		self.C=tess.Camera(stamp=self.size, cadence=self.cadence)

		# point that camera somewhere
		if name is None:
			self.ra = 258.828916667 # np.random.random()*360
			self.dec =  4.96380555556 # np.random.random()*180 -90
			print "   pointing at GJ1214"
		else:
			s = zachopy.star.SingleStar(name)
			self.ra = s.icrs.ra.degree
			self.dec = s.icrs.dec.degree
			print "   pointing at ", name

		self.C.point(self.ra, self.dec)

		# create a blank image with the camera
		self.I = tess.Image(self.C)
		self.I.image = self.I.zeros()

		# create empty (xpixels, ypixels, n)
		self.photons, self.cosmics, self.noiseless = np.zeros(self.shape), np.zeros(self.shape), np.zeros(self.shape)

		# populate the cube with simulated pixel data
		#self.simulate()
		#self.plot()

	def master(self):
		try:
			self.meanimage
		except:
			self.meanimage = np.mean(self.noiseless, 2)
		return self.meanimage.reshape(self.xpixels, self.ypixels, 1)

	def simulate(self):
		'''Use tess simulator to populate the image cube.'''
		for i in range(self.n):
			self.photons[:,:,i], self.cosmics[:,:,i], self.noiseless[:,:,i] = self.I.expose(jitter=self.jitter, write=False)

	def plot(self, normalization='None'):

		self.master()
		if normalization.lower() == 'none':
			normalizationimage = 1.0
			ylabel = 'Photons'
		elif normalization.lower() == 'master':
			normalizationimage = self.master()
			ylabel='Relative Flux'

		# create a relative light curve (dF/F)
		photonsnormalized = self.photons/normalizationimage
		cosmicsnormalized = self.cosmics/normalizationimage
		noiselessnormalized = self.noiseless/normalizationimage

		# set up a logarithmic color scale (going between 0 and 1)
		def color(x):
			zero = np.min(np.log(self.master()))
			span = np.max(np.log(self.master())) - zero
			normalized = (np.log(x) -  zero)/span
			return cm.YlGn(normalized)

		# create a plot
		scale = 1.5
		plt.figure(figsize = (np.minimum(self.xpixels*scale,10),np.minimum(self.ypixels*scale, 10)))
		gs = gridspec.GridSpec(self.xpixels,self.ypixels, wspace=0, hspace=0)


		for i in range(self.ypixels):
			for j in range(self.xpixels):

				try:
					sharex = ax
					sharey = ax
				except:
					sharex = None
					sharey = None


				ax = plt.subplot(gs[i,j], sharex=sharex, sharey=sharey)
				ax.patch.set_facecolor(color(self.meanimage[i,j]))
				ax.plot(photonsnormalized[i,j,:], color='black')
				ax.plot(noiselessnormalized[i,j,:], color='blue', alpha=0.3)
				ax.plot(cosmicsnormalized[i,j,:], color='green', alpha=0.3)

				if i == (self.ypixels-1) and j == 0:
					plt.setp(ax.get_xticklabels(), rotation=90)
					ax.set_xlabel('Time')
					ax.set_ylabel(ylabel)
				else:
					plt.setp(ax.get_xticklabels(), visible=False)
					plt.setp(ax.get_yticklabels(), visible=False)


		ax.set_ylim(np.min(photonsnormalized), np.max(photonsnormalized))

		plt.draw()
