'''Generate TESS pixel lightcurve cubes with dimensions (xpix)x(ypix)x(time).'''
import tess
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import zachopy.star

class Cube(object):
	'''Cube to handle simulated postage stamp pixel light curves;
			has dimensions of (xpixels, ypixels, time).'''

	def __init__(self, subject='test', size=32, n=900, cadence=2, jitter=False):
		'''Initialize a cube object.

		keyword arguments:
		[subject='test'] the patch of the sky to paint, either a star name that can be queried by Vizier or "test" for a test pattern grid of stars
		[size=32] the size of the image, in pixels
		[n=900] the number of exposures to simulate
		[cadence=2] the cadence of the simulated images, in seconds
		[jitter=False] should jitter be included between images?'''

		# define a prefix for output
		self.prefix = '   [cube] '

		# decide what to point the camera at...
		self.subject = subject
		if 'test' in self.subject.lower():
			# ...either a test pattern of stars...
			print self.prefix + 'pointing at a test pattern'
			self.testpattern = True
		else:
			self.testpattern = False
			# ...or a field centered on a particular star.
			try:
				print self.prefix + "trying to point at ", subject
				s = zachopy.star.SingleStar(name)
				self.ra = s.icrs.ra.degree
				self.dec = s.icrs.dec.degree
			except:
				print self.prefix + " but that failed, so we'll go with a test pattern"
				self.testpattern = True

		# define basic geometry
		self.size = size
		self.xpixels = self.size
		self.ypixels = self.size
		self.jitter = jitter

		# number of images and cadence
		self.n = n
		self.cadence = cadence
		self.shape = (self.xpixels, self.ypixels, self.n)

		# create a TESS camera and point it at the subject
		self.C=tess.Camera(stamp=self.size, cadence=self.cadence, testpattern=self.testpattern)

		# create a blank image with the camera
		self.I = tess.Image(self.C)
		self.I.image = self.I.zeros()

		# create empty (xpixels, ypixels, n)
		self.photons, self.cosmics, self.noiseless = np.zeros(self.shape), np.zeros(self.shape), np.zeros(self.shape)

		# create a dictionary to store a bunch of summaries
		self.summaries = {}

		# populate the cube with simulated pixel data
		# self.simulate()
		# self.plot()

	@property
	def directory(self):
		'''Return path to a directory where this cube's data can be stored.'''
		dir = self.I.directory + 'cubes/'
		zachopy.utils.mkdir(dir)
		return dir

	@property
	def filename(self):
		'''Return a filename for saving/loading this cube.'''
		if self.jitter:
			phrase = 'with'
		else:
			phrase = 'without'
		return self.directory + 'cube_{n:.0f}exp_at{cadence:.0f}s_{phrase}jitter.npy'.format(n=self.n, cadence=self.cadence, phrase=phrase)

	def simulate(self):
		'''Use TESS simulator to paint stars (and noise and cosmic rays) into the image cube.'''
		for i in range(self.n):
			self.photons[:,:,i], self.cosmics[:,:,i], self.noiseless[:,:,i] = self.I.expose(jitter=self.jitter, write=False, smear=False, remake=i==0)

	def save(self):
		'''Save this cube a 3D numpy array (as opposed to a series of FITS images).'''
		print " saving cube to " + self.filename
		np.save(self.filename, (self.photons, self.cosmics, self.noiseless))

	def load(self):
		'''Load this cube from a 3D numpy array, assuming one exists with the appropriate size, for this field, at this cadence, with the right number of exposures.'''
		print " loading simulated cube from " + self.filename
		try:
			self.photons, self.cosmics, self.noiseless = np.load(self.filename)
			print self.prefix + 'Loaded cube from ' + self.filename
		except:
			print self.prefix + 'No saved cube was found; generating a new one.\n          (looked in {0})'.format( self.filename)
			self.simulate()
			self.save()

	def cubify(self, image):
		'''Slightly reshape an image, so it can be cast into operations on the whole cube.'''
		return image.reshape(self.xpixels, self.ypixels, 1)

	@property
	def median(self):
		'''The median image.'''
		key = 'median'
		try:
			self.summaries[key]
		except:
			self.summaries[key] = np.median(self.photons, 2)
		return self.summaries[key]

	@property
	def mad(self):
		'''The median of the absolute deviation image.'''
		key = 'mad'
		try:
			self.summaries[key]
		except:
			self.summaries[key] = np.median(np.abs(self.photons - self.cubify(self.median)), 2)
		return self.summaries[key]

	@property
	def master(self):
		'''The (completely noiseless and unsaturated) raw master frame.'''
		key = 'master'
		try:
			self.summaries[key]
		except:
			self.summaries[key] = np.mean(self.noiseless, 2)
		return self.summaries[key]

	def nsigma(self):
		sigma = 1.48*self.cubify(self.mad)
		return (self.photons - self.cubify(self.median))/sigma

	def write(self, normalization='none'):
		'''Save all the images to FITS, inside the cube directory.'''

		# make a directory for the normalization used
		dir = self.directory + normalization + '/'
		zachopy.utils.mkdir(dir)

		if normalization == 'none':
			flux = self.photons
		if normalization == 'nsigma':
			flux = self.nsigma()

		# loop through the images in the cube
		for i in range(self.n):
			# pick some kind of normalization for the image
			image = flux[:,:,i]
			self.I.writeToFITS(image, dir + normalization + '_{0:05.0f}.fits'.format(i))

	def plot(self, normalization='median'):

		# choose how to normalize the lightcurves for plotting
		if normalization.lower() == 'none':
			normalizationarray = 1.0
			ylabel = 'Photons'
		elif normalization.lower() == 'master':
			normalizationarray = self.master().reshape(self.xpixels, self.ypixels, 1)
			ylabel='Relative Flux'
		elif normalization.lower() == 'median':
			normalizationarray = self.median().reshape(self.xpixels, self.ypixels, 1)
			ylabel='Relative Flux'

		# create a relative light curve (dF/F, in most cases)
		photonsnormalized = self.photons/normalizationarray
		cosmicsnormalized = self.cosmics/normalizationarray
		noiselessnormalized = self.noiseless/normalizationarray

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

		# loop over pixels (in x and y directions)
		for i in range(self.ypixels):
			for j in range(self.xpixels):

				# set up axis sharing, so zooming on one plot zooms the other
				try:
					sharex, sharey = ax, ax
				except:
					sharex, sharey = None, None
				ax = plt.subplot(gs[i,j], sharex=sharex, sharey=sharey)

				# color the plot panel based on the pixel's intensity
				ax.patch.set_facecolor(color(self.median()[i,j]))


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
