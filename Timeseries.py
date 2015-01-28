from imports import *
import settings
import transit

def drawCosmics(n):
	'''Draw random pixel values from a histogram of cosmic ray values (most of them are 0).'''
	def weighted_values(values, probabilities, size):
		bins = np.add.accumulate(probabilities)
		return values[np.digitize(np.random.random_sample(size), bins)]

	x, y = np.load(settings.dirs['intermediates']+'cosmics_histogram.npy')
	y /= np.sum(y)
	return weighted_values(x,y,n)


class Timeseries(Talker):
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
				self.scale =  self.nsubexposures

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

	def createTransit(self, cosmics=True, mag=10.0, period = 0.9351425135, rs_over_a = 0.1, k=0.1):

		# add in a transit
		self.createSimple(noise=False, cosmics=False)
		p = transit.Planet(period=period, t0=period/2.0, rs_over_a=rs_over_a, k=k)
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
