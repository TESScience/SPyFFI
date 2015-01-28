'''Generate TESS pixel lightcurve cubes with dimensions (xpix)x(ypix)x(time).'''
from imports import *
import SPyFFI, cr

class Cube(Talker):
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

		# decide whether or not this Cube is chatty
		Talker.__init__(self, mute=False, pithy=False)

		# decide what to point the camera at...
		self.subject = subject
		if 'test' in self.subject.lower():
			# ...either a test pattern of stars...
			self.speak( 'pointing at a test pattern')
			self.testpattern = True
		else:
			self.testpattern = False
			# ...or a field centered on a particular star.
			try:
				self.speak( "trying to point at {0}".format(subject))
				s = zachopy.star.SingleStar(name)
				self.ra = s.icrs.ra.degree
				self.dec = s.icrs.dec.degree
			except:
				self.speak( " but that failed, so we'll go with a test pattern")
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
		self.camera = SPyFFI.Camera(subarray=self.size, cadence=self.cadence, testpattern=self.testpattern)
		self.camera.cartographer.pithy = True

		# create a blank image with the camera
		self.ccd = self.camera.ccds[0]
		self.ccd.image = self.ccd.zeros()

		# create empty (xpixels, ypixels, n)
		if self.cadence == 2:
			bits = np.float32
		else:
			bits = np.float64
		self.photons, self.cosmics, self.noiseless = np.zeros(self.shape).astype(bits), np.zeros(self.shape).astype(bits), np.zeros(self.shape).astype(bits)

		# create a dictionary to store a bunch of summaries
		self.summaries = {}

		# populate the cube with simulated pixel data
		# self.simulate()
		# self.plot()

	def bin(self, nsubexposures=60, strategy=cr.central(n=10), plot=False):
		'''Bin together 2-second exposures, using some cosmic strategy.'''

		# make sure that we're starting with a 2-second exposure
		assert(self.cadence == 2)

		# create an empty binned cube object that has the right size and shape
		binned = Cube(subject=self.subject, size=self.size, cadence=self.cadence*nsubexposures, n=(np.int(self.n/nsubexposures)) )

		# add an additional array to keep track of what the unmitigated lightcurves would look like
		binned.unmitigated = np.zeros_like(binned.photons)

		# loop over the x and y pixels
		self.speak('binning {0} cube by {1} subexposures into a {2} cube'.format(self.shape, nsubexposures, binned.shape))
		for x in np.arange(self.shape[0]):
			for y in np.arange(self.shape[1]):
				self.speak('   {x}, {y} out of ({size}, {size})'.format(x=x, y=y, size=self.size))
				timeseries = cr.timeseries(self, (x,y), nsubexposures=nsubexposures)
				strategy.calculate(timeseries)
				if plot:
					strategy.plot()
				binned.photons[x,y] = strategy.binned['flux']
				binned.cosmics[x,y] = strategy.binned['naive'] - strategy.binned['nocosmics']
				binned.unmitigated[x,y] = strategy.binned['naive']
		return binned



	def display(self, cube=None, name='cube', limit=50):
		'''Use ds9 to display the image cube.'''
		self.speak('displaying (up to {0:.0f} exposures) of the pixel cube with ds9'.format(limit))
		if cube is None:
			cube = self.photons
		try:
			self.ds9
		except:
			self.ds9 = zachopy.display.ds9(name)
		self.ds9.many(cube, limit=limit)

	@property
	def directory(self):
		'''Return path to a directory where this cube's data can be stored.'''
		dir = self.ccd.directory + 'cubes/'
		zachopy.utils.mkdir(dir)
		return dir

	@property
	def filename(self):
		'''Return a filename for saving/loading this cube.'''
		if self.jitter:
			phrase = 'with'
		else:
			phrase = 'without'
		return self.directory + 'cube_{n:.0f}exp_at{cadence:.0f}s_{phrase}jitter_{intrapixel}.npy'.format(n=self.n, cadence=self.cadence, phrase=phrase, intrapixel=self.camera.psf.intrapixel.name)

	def simulate(self):
		'''Use TESS simulator to paint stars (and noise and cosmic rays) into the image cube.'''
		self.ccd.display =0
		self.camera.cartographer.pithy = False
		self.speak('populating the {0} cube with {1:.0f}-second exposures'.format(self.shape, self.cadence))
		for i in range(self.n):
			self.speak('filling exposure #{0:.0f}/{1:.0f}'.format(i, self.n))
			self.photons[:,:,i], self.cosmics[:,:,i], self.noiseless[:,:,i] = self.ccd.expose(jitter=self.jitter, write=False, smear=False, remake=i==0, terse=True, cosmics='fancy')
		self.background = self.ccd.backgroundimage

	def save(self):
		'''Save this cube a 3D numpy array (as opposed to a series of FITS images).'''
		self.speak( "Saving cube to " + self.filename)
		np.save(self.filename, (self.photons, self.cosmics, self.noiseless))

	def load(self):
		'''Load this cube from a 3D numpy array, assuming one exists with the appropriate size, for this field, at this cadence, with the right number of exposures.'''
		self.speak("Trying to load simulated cube from " + self.filename)
		try:
			self.photons, self.cosmics, self.noiseless = np.load(self.filename)
			self.speak( 'Loaded cube from ' + self.filename)
		except:
			self.speak('No saved cube was found; generating a new one.\n          (looked in {0})'.format( self.filename))
			self.simulate()
			self.save()

	def cubify(self, image):
		'''Slightly reshape an image, so it can be cast into operations on the whole cube.'''
		return image.reshape(self.xpixels, self.ypixels, 1)

	def median(self, which='photons'):
		'''The median image.'''
		array = self.__dict__[which].astype(np.float64)
		key = 'median'
		try:
			self.summaries[key+which]
		except:
			self.summaries[key+which] = np.median(array, 2)
		return self.summaries[key+which]

	def mean(self, which='photons'):
		'''The median image.'''
		array = self.__dict__[which].astype(np.float64)
		key = 'mean'
		try:
			self.summaries[key+which]
		except:
			self.summaries[key+which] = np.mean(array, 2)
		return self.summaries[key+which]


	def mad(self, which='photons'):
		'''The median of the absolute deviation image.'''
		array = self.__dict__[which].astype(np.float64)
		key = 'mad'
		try:
			self.summaries[key+which]
		except:
			self.summaries[key+which] = np.median(np.abs(array - self.cubify(self.median(which))), 2)
		return self.summaries[key+which]

	def std(self, which='photons'):
		'''The standard deviation image.'''
		array = self.__dict__[which].astype(np.float64)
		key = 'std'
		try:
			self.summaries[key+which]
		except:
			self.summaries[key+which] = np.std(array, 2)
		return self.summaries[key+which]

	def sigma(self, which='photons', robust=True):
		if robust:
			return 1.4826*self.mad(which)
		else:
			return self.std(which)

	def master(self, which='photons'):
		'''The (calculated) master frame.'''
		return self.median(which)

	def nsigma(self, which='photons', robust=True):
		array = self.__dict__[which].astype(np.float64)
		return (array - self.cubify(self.median(which)))/self.cubify(self.sigma(which, robust=robust))

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
			self.ccd.writeToFITS(image, dir + normalization + '_{0:05.0f}.fits'.format(i))

	def plot(self, normalization='none'):

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
			return plt.matplotlib.cm.YlGn(normalized)

		# create a plot
		scale = 1.5
		plt.figure(figsize = (np.minimum(self.xpixels*scale,10),np.minimum(self.ypixels*scale, 10)))
		gs = plt.matplotlib.gridspec.GridSpec(self.xpixels,self.ypixels, wspace=0, hspace=0)

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
