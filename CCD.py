'''Tools to generate simulated TESS images.'''
# import necessary packages
import settings
from imports import *
import Cosmics

# setup basic output options for this Python session
np.set_printoptions(threshold = 1e6, linewidth = 300)
plt.ion()

# define mapping between CCD number and quadrant
quadrants = {1:(1,1), 2:(-1,1), 3:(-1,-1), 4:(1,-1), 0:None}

class CCD(Talker):
	def __init__(self, number=1, camera=None, subarray=None, label='', display=False):
		'''Turn on a TESS CCD, which can be used to make simulated images.

			camera=None is parent TESS camera for this CCD, required for *any* conversion from (RA, Dec) to (x,y)
			number=1 is the CCD number of the detector is it? (1-4, like mathematical quadrants)
			subarray=None is the size of a square subarray that can optionally be defined for speediness
			label='' is a special label you can add for experimenting
			'''
		# decide whether or not this CCD is chatty
		Talker.__init__(self, mute=False, pithy=False)

		# keep track of where this image is
		self.number = number
		self.camera = camera
		self.subarray = subarray

		# some labeling for the image
		self.label = label
		self.note = ''


		# define the size and location of the detector within the field
		if subarray is None:
			# image size is the whole detector
			self.npix = 2048
			# image center (in focalxy coordinates) is set by the image size and the camera's CCD gap
			self.center = self.camera.gapinpixels/2.0 + self.npix/2.0*np.array(quadrants[number])
		else:
			# image size is the specified subarray size
			self.npix = subarray
			# image center is the center of the field (for now)
			self.center = np.array([0.0,0.0])

		# some numbers for the image
		self.xsize = self.npix
		self.ysize = self.npix
		self.xmin = 0
		self.xmax = self.xmin + self.xsize
		self.ymin = 0
		self.ymax = self.ymin + self.ysize

		# pull out the camera's position string, for saving images
		self.pos_string = self.camera.pos_string()

		# a few misc diagnostics
		self.display = display
		self.plot = False
		self.starcounter = 0
		self.nstars =0

		# start populating the image header (seems like we can't do this until we're sure camera has pointed somewhere)
		self.populateHeader()

	def show(self):
		'''Display the current (possibly in-progress image.)'''
		if self.display:
			try:
				self.ds9
			except:
				self.ds9 = zachopy.display.ds9(self.name.replace(' ',''))
			self.ds9.one(self.image)

	@property
	def name(self):
		'''Simple name for this CCD, for saving and printing.'''

		# define a base name, whether this is a subarray or a full ccd
		if self.number == 0:
			str = 'sub{0:d}x{0:d}'.format(self.subarray)
		else:
			str = 'ccd{0:d}'.format(self.number)

		# add label, if one exists
		if self.label != '':
			str += '_' + self.label

		return str


	@property
	def directory(self):
		'''Directory for saving all images from this CCD.'''
		d = self.camera.directory + self.name + '/'
		zachopy.utils.mkdir(d)
		return d

		#@property
		#def wcs(self):
		#	'''The WCS for this CCD, automatically derived from the Camera's WCS.'''
		#
		#	# what is the imagexy pixel coordinate of the center of the field?
		#	self.wcs.wcs.crpix =[self.npix/2.0,self.npix/2.0]


	def photons(self, mag):
		'''Use magnitude to calculate photons per second that will be recorded on a CCD.'''

		# this doesn't work great for M dwarfs, need to include multiband information at some point
		return self.camera.effective_area*10**(-0.4*mag)*1.7e6

	def populateHeader(self, ccd=None):
		'''Populate the header structure with information about the CCD image, and the inputs.'''

		# create an empty header
		try:
			self.header
		except:
			self.header = astropy.io.fits.Header()



		# fill it with some CCD details
		self.header['CCD'] = ''
		self.header['IMAGNOTE'] = ('', 'Details of this individual image')
		self.header['EXPTIME'] = (self.camera.cadence, '[s] exposure time ')
		self.header['NREADS'] = (np.round(self.camera.cadence/self.camera.singleread).astype(np.int), '# of reads summed')
		self.header['SUBEXPTI'] = (self.camera.singleread, '[s] exposure in a single subexposure')
		self.header['SATURATE'] = (self.camera.saturation*self.camera.cadence/self.camera.singleread, '[e-] saturation level in this image')
		self.header['READNOIS'] = (self.camera.read_noise, '[e-] read noise (per individual read)')
		self.header['READTIME'] = (self.camera.readouttime, '[s] time to transer to frame store')
		self.header['CCDNUM'] = (self.number, 'CCD number (1,2,3,4 or 0=fake subarrayor)')
		self.header['CCDSIZE'] = (self.npix, '[pix] size of one CCD')

		# fill in the timestamp for this CCD image
		self.setTime()

		# leave space to talk about the inputs to this image
		self.header['INPUTS'] = ''
		self.header['INPUNOTE'] = ('', 'Ingredients for simulated images.')

	def setTime(self):
		'''Based on the camera counter, apply a timestamp to this image.'''

		# add time to the image (in a very simplistic way -- no accounting for the spacecraft orbit)
		self.header['COUNTER'] = self.camera.counter, '# of exposures since start, for this field'
		self.header['BJD'] = self.camera.counter*self.camera.cadence/24.0/60.0/60.0, '[day] mid-exposure time - 2457827.0 (e.g.)'

	def pixels(self):
		'''Give grids of x and y values (2D arrays) of the image pixels.'''

		# use meshgrid to create (x,y) arrays, using default 'xy' indexing (x increases with column, y increases with row)
		x,y = np.meshgrid(np.arange(self.xsize) + self.xmin, np.arange(self.ysize) + self.ymin)#, indexing='ij')
		pix = self.camera.cartographer.point(x,y,'ccdxy')
		return pix

	def zeros(self):
		'''Create an image of zeros, the same size as the CCD.'''
		return np.zeros((self.xsize, self.ysize))

	def ones(self):
		'''Create an image of ones, the same size as the CCD.'''
		return np.ones((self.xsize, self.ysize))

	def zodicalBackground(self,elon, elat):
		'''Calcaulte the zodiacal background at a given ecliptic (lat, long).'''

		# from Josh and Peter's memo
		Vmax = 23.345
		DeltaV = 1.148
		b = np.abs(elat)
		V = Vmax  - DeltaV*((b-90.0)/90.0)**2
		assert((b < 90).all())
		return 10**(-0.4*(V-22.8))*(2.56e-3)*self.camera.effective_area*self.camera.pixel_solid_area

	def unresolvedBackground(self, glon, glat, complicated=False):
		'''Calculate the unresolved stellar background at a given galactic (lat, long).'''

		# from Josh and Peter's memo
		flip = glon > 180.0
		if np.sum(flip):
			glon[flip] -= 360.0

		if complicated:
			L = np.abs(glon/180.0)
			B = np.abs(glat/90.0)
			a1 = 18.7
			a2 = 4.3
			a3 = 0.52
			a4 = 10.2
			a5 = 0.46
			a6 = -3.74
			I_surface_brightness = a1 + a2*(1.0-np.exp(-L/a3)) + a4*(1.0-np.exp(-B/a5)) + a6*np.sqrt(L*B)
		else:
			a0 = 18.9733
			a1 = 8.833
			a2 = 4.007
			a3 = 0.805
			I_surface_brightness = a0 + a1*(np.abs(glat)/40.0) + a2*(np.abs(glon)/180.0)**a3
		return 10**(-0.4*I_surface_brightness)*1.7e6*self.camera.effective_area*self.camera.pixel_solid_area


	def writeToFITS(self, image, path, split=False, savetype=np.float32):
		'''General FITS writer for this CCD.'''

		# print status
		self.speak('saving {0}x{0} image to {1} with type {2}'.format(image.shape[0],path,savetype.__name__))
		cr = self.camera.cartographer.point(0.0, 0.0, 'focalxy')
		x, y = cr.ccdxy.tuple
		# modify the camera's WCS, based on the CCD number
		self.header['CRPIX1'] = x
		self.header['CRPIX2'] = y

		# write the file to FITS
		#astropy.io.fits.PrimaryHDU(np.transpose(savetype(image)), header=self.header).writeto(filename, clobber=True)
		astropy.io.fits.PrimaryHDU(savetype(image), header=self.header).writeto(path, clobber=True)


	def loadFromFITS(self, path):
		'''General FITS loader for this CCD.'''

		# print status
		self.speak('trying to load image from {0}'.format(path))


		try:
			# can we get by without the transposes?
			#image = np.transpose(astropy.io.fits.open(path)[0].data)
			image = astropy.io.fits.open(path)[0].data
			self.speak("       ...success!")
			return image
		except:
			self.speak("       ...failed.")
			assert(False)


	def writeFinal(self, lean=True):
		'''Write the final image from this CCD.'''

		self.header.extend(self.camera.header)
		self.header.extend(self.camera.psf.header)

		# make filename for this image
		self.note = 'final_{0:06.0f}'.format(self.camera.counter)
		finalfilename = self.directory + self.note + '.fits'

		# write the image to FITS
		self.speak('saving final TESS image')
		self.writeToFITS(self.image, finalfilename, savetype=np.int32)

		# optionally, write some other outputs too!
		if lean == False:
			self.note = 'withoutbackground_{0:06.0f}'.format(self.camera.counter)
			self.writeToFITS(self.image - self.image_background, self.directory + self.note + '.fits')



	def projectCatalog(self, write=True):
		'''Create using the camera's star catalog, and project stars using this CCD.'''
		self.speak('projecting the starmap onto CCD')

		# make sure the camera has a catalog defined
		try:
			self.camera.catalog
		except:
			self.camera.populateCatalog()

		# pull out positions, magnitudes, and temperatures
		ras, decs, tmag, temperatures = self.camera.catalog.arrays()
		self.camera.cartographer.ccd = self

		# create coordinate object for the stars
		stars = self.camera.cartographer.point(ras, decs, 'celestial')
		x,y = stars.ccdxy.tuple

		# trim everything down to only those stars that could be relevant for this camera
		buffer = 10
		ok = (x > -buffer) & (x < self.xsize + buffer) & (y > -buffer) & (y < self.ysize + buffer)

		# assign this CCD's stars
		self.starx = x[ok]
		self.stary = y[ok]
		self.starmag = np.array(tmag)[ok]
		self.startemp = np.array(temperatures)[ok]

		# keep track of which CCD we projected onto
		self.starsareon = self.name

		# write the catalog to a text file
		if write:
			outfile = self.directory + 'catalog_{pos}_{name}.txt'.format(pos=self.pos_string, name=self.name)
			np.savetxt(outfile, np.c_[ras[ok], decs[ok], self.starx, self.stary, self.starmag], fmt=['%.6f', '%.6f', '%.3f', '%.3f', '%.3f'])
			self.speak("save projected star catalog {0}".format(outfile))

	def addStar(self, ccdx, ccdy, mag, temp, verbose=True, plot=False):
		'''Add one star to an image, given position, magnitude, and effective temperature.'''


		self.speak("adding stars at ({0}, {1}) with magnitude of {2}".format(ccdx, ccdy, mag))
		self.speak(" ({0}/{1})".format(self.starcounter,self.nstars))

		# (this is probably a slow way of doing things -- speed it up!)
		ccdxy = self.camera.cartographer.point(ccdx + self.camera.nudge['x']/self.camera.pixelscale,ccdy + self.camera.nudge['y']/self.camera.pixelscale,'ccdxy')

		focalx, focaly = ccdxy.focalxy.tuple

		normalized, xindex, yindex = self.camera.psf.pixelizedPSF(ccdxy,temp)
		binned = normalized*self.camera.cadence*self.photons(mag)
		#binned = unnormed*self.camera.cadence*self.photons(mag)/np.sum(unnormed)

		'''if plot:
			try:
				self.ax_prnu.figure.clf()
			except:
				fi, ax = plt.subplots(1,4,figsize=(20,4), sharex=True, sharey=True)
				self.ax_prnu = ax[0]
				self.ax_psf = ax[1]
				self.ax_highres = ax[2]
				self.ax_prf = ax[3]
			extent = [self.psf.xgrid.min(), self.psf.xgrid.max(), self.psf.ygrid.min(), self.psf.ygrid.max()]
			self.ax_prnu.imshow(intrapixel,extent=extent,cmap='gray_r',interpolation='nearest')
			self.ax_psf.imshow(subgrid_psf,extent=extent,cmap='gray_r',interpolation='nearest')
			self.ax_highres.imshow(subgrid_psf*intrapixel,extent=extent,cmap='gray_r',interpolation='nearest')
			self.ax_prf.imshow(binned,extent=extent,cmap='gray_r',interpolation='nearest')
			self.ax_prf.set_xlim(-self.psf.dx_pixels, self.psf.dy_pixels)
			self.ax_prf.set_ylim(-self.psf.dx_pixels, self.psf.dy_pixels)
			self.ax_prf.set_aspect(1)
			self.ax_prf.figure.savefig(settings.prefix + 'plots/prnu_demo.pdf')
			plt.draw()'''

		ok = (xindex >= self.xmin) * (xindex < self.xsize) * (yindex >= self.ymin) * (yindex < self.ysize)
		self.starimage[yindex[ok], xindex[ok]] += binned[ok]

	def addStars(self, remake=False, jitter=False):
		self.speak("Adding stars.")

		if jitter:
			remake=True

		# define a grid of magnitude thresholds, will save an image containing all stars brighter than each
		dthreshold = 1
		magnitude_thresholds = np.arange(6,18,dthreshold)

		# if the final star image already exists, just load it
		try:
			assert(remake == False)
			self.note = 'starsbrighterthan{0}'.format(np.max(magnitude_thresholds))
			starsfilename = self.directory + self.note + '.fits'
			try:
				self.starimage
			except:
				self.starimage = self.loadFromFITS(starsfilename)

		# otherwise loop through thresholds, adding stars at each
		except:
			self.starimage = self.zeros()
			try:
				self.starx
				assert(self.starsareon == self.name)
			except:
				self.projectCatalog()

			for threshold in magnitude_thresholds:

				# define a filename for this magnitude range
				self.note = 'starsbrighterthan{0:02d}'.format(threshold)
				starsfilename = self.directory + self.note + '.fits'

				# load the existing stellar image, if possible
				try:
					assert(remake==False)
					self.starimage = self.loadFromFITS(starsfilename)
				except:
					# if this is the smallest threshold, include all the stars brighter than it
					if threshold == np.min(magnitude_thresholds):
						min = -100
					else:
						min = threshold - dthreshold
					# pick the stars to add to the image on this pass through
					ok =    (self.starx + self.camera.psf.dx_pixels_axis[-1] >= self.xmin) * \
							(self.starx + self.camera.psf.dx_pixels_axis[0] <= self.xmax) * \
							(self.stary + self.camera.psf.dy_pixels_axis[-1] >= self.ymin) * \
							(self.stary + self.camera.psf.dy_pixels_axis[0] <= self.ymax) * \
							(self.starmag < threshold)*(self.starmag >= min)
					x = self.starx[ok]
					y = self.stary[ok]
					mag = self.starmag[ok]
					temp = self.startemp[ok]

					if np.sum(ok) > 0:
						self.nstars += np.sum(ok)
						for i in range(len(x)):
							self.addStar(x[i], y[i], mag[i], temp[i])
							self.starcounter += 1

						if jitter == False:
							self.writeToFITS(self.starimage, starsfilename)
					self.show()

		self.image += self.starimage
		self.show()
		self.header['ISTARS'] = ('True', 'stars from UCAC4')
		if jitter:
			self.header['IJITTER'] = ('True', 'spacecraft jitter, motion between images')
		return self.starimage

	def addGalaxies(self):
		pass
		# looks like I should use http://vizier.cfa.harvard.edu/viz-bin/Cat?VII/155 for a source catalog?


	def addCosmics(self, gradient=False, version='fancy', diffusion=False, write=False, rate=5.0):
		'''Add cosmic rays to image.'''

		# print update
		self.speak('Adding cosmic rays')

		# filenames, in case saving is required


		# use Al's code to generate cosmic ray image of the correct size
		image = Cosmics.cosmicImage(exptime=self.camera.cadence, size=self.npix,
									gradient=gradient, version=version, diffusion=diffusion, rate=rate)

		# (optionally), write cosmic ray image
		if write:
			self.note = 'cosmics_{0:.1f}persecond_{1}seconds_{2:06.0f}'.format(rate,self.camera.cadence, self.camera.counter).replace('.', 'p')
			cosmicsfilename = self.directory + self.note + '.fits'
			self.writeToFITS(image, cosmicsfilename)

		# add the cosmics into the running image
		self.image += image
		# note that cosmics were included in the image header
		self.header['ICOSMICS'] = ('True', 'cosmic rays injected')
		self.show()
		return image

	def bleedSaturated(self, plot=False):
		'''Bleed saturated pixels in the image.'''

		self.speak('Bleeding saturated pixels.')

		# keep track of the original image
		untouched = self.image + 0.0
		original = np.sum(self.image)

		# set the saturation limit based on the number of individual reads included
		saturation_limit =  self.camera.saturation*self.camera.cadence/self.camera.singleread
		stilloversaturated = True

		# keep looping until all saturation problems are gone
		count =0
		while stilloversaturated:
			count +=1
			oversaturated = self.image > saturation_limit
			saturated = self.image >= saturation_limit
			for x in range(self.image.shape[1]):
				regions, nregions = scipy.ndimage.measurements.label(oversaturated[x,:])
				if nregions > 0:
					for i in np.arange(nregions)+1:
						y = (regions == i).nonzero()[0]
						if oversaturated[x,y].any():
							fluxtodistribute = np.sum(self.image[x,y])
							npixels = (fluxtodistribute/saturation_limit)
							center = np.int(np.mean(y))
							grow = np.int(np.floor(npixels/2 - 0.5))
							indices = np.arange(np.maximum(center - grow, 0),np.minimum(center + grow, self.image.shape[0]-1)+1)
							assert(y[0] in indices)

							existingflux = np.sum(self.image[x,indices])
							self.image[x,indices] = saturation_limit
							leftoverflux = existingflux - indices.shape[0]*saturation_limit
							#print "       distributing flux of {0:10.5f}X saturated pixels over {1:4} + {2:.5f} pixels".format(fluxtodistribute/saturation_limit, indices.shape[0], leftoverflux/saturation_limit)
							#print x, y
							#print indices
							leftedge = center - grow -1
							rightedge = center + grow +1
							try:
								try:
									self.image[x,leftedge] += leftoverflux/2.0
								except:
									self.image[x,rightedge] += leftoverflux/2.0
								try:
									self.image[x,rightege] += leftoverflux/2.0
								except:
									self.image[x,leftedge] += leftoverflux/2.0
							except:
								print "    this star seems to saturate the entire detector!"
			self.speak("on pass #{0} through saturation filter,  \n        the max saturation fraction is {1} \n        and the flux change over entire image is {2} electrons".format(count, np.max(self.image)/saturation_limit, np.sum(self.image) - original))

			# KLUDGE to prevent endless loops
			stilloversaturated = (self.image > saturation_limit).any() and count < 10

		self.note = 'saturation_{0}K'.format(self.camera.saturation).replace('.', 'p')
		saturationfilename = self.directory + self.note + '.fits'
		if not os.path.exists(saturationfilename):
			self.writeToFITS(self.image - untouched,saturationfilename)

		# update image header
		self.header['ISATURAT'] = ('True', 'bleed trails for saturated pixels')

		self.show()

	def addBackgrounds(self):
		'''Add smooth backgrounds (zodiacal light and unresolved stars) to background.'''

		# set up filenames for saving background, if need be
		self.note = 'backgrounds'
		backgroundsfilename = self.directory + self.note + '.fits'
		self.speak("Adding backgrounds.")

		# if the background image already exists, just load it
		try:
			self.backgroundimage
		except:

			try:
				self.backgroundimage = self.loadFromFITS(backgroundsfilename)
			# otherwise compute the backgrounds from scratch, and save them for next time
			except:
				# define a blank background image
				self.backgroundimage = self.zeros()
				# define coordinates (equatorial, Galactic, ecliptic) at every pixel in the image
				ra, dec = self.pixels().celestial.tuple
				elon, elat = self.pixels().ecliptic.tuple
				glon, glat = self.pixels().galactic.tuple

				# add the zodiacal light, using the simple model from Josh and Peter on the TESS wiki
				print "   including smooth model for zodiacal light."
				elon, elat
				self.backgroundimage += self.zodicalBackground(elon, elat)*self.camera.cadence
				self.header['IZODIACA'] = ('True', 'zodiacal light, treated as smooth')

				# add unresolved background light, using the simple model from Josh and Peter on the TESS wiki
				print "   including smooth model for unresolved stars in the Galaxy."
				self.backgroundimage += self.unresolvedBackground(glon, glat)*self.camera.cadence
				self.header['IUNRESOL'] = ('True', 'unresolved stars, treated as smooth background')

				# write the image, so it can just be loaded easily next time
				self.writeToFITS(self.backgroundimage, backgroundsfilename)

		# add the background image to the total image
		self.image += self.backgroundimage
		self.show()

	def addPhotonNoise(self):
		'''Add photon noise into an image.'''

		self.speak("Adding photon noise [sqrt(photons from stars and various backgrounds)].")
		noise_variance = self.image
		ok = noise_variance > 0

		noise = np.zeros_like(self.image)
		noise[ok] = np.sqrt(noise_variance[ok])

		assert(np.isfinite(noise).all())
		self.image += noise*np.random.randn(self.xsize, self.ysize)

		self.note = 'photonnoise'
		noisefilename = self.directory + self.note + '.fits'
		if not os.path.exists(noisefilename):
			self.writeToFITS(noise,noisefilename)
		self.header['IPHOTNOI'] = ('True', 'photon noise')
		self.show()

	def addReadNoise(self):
		'''Add read noise to image.'''
		self.speak("Adding read noise.")
		self.speak("    = quadrature sum of [{0} seconds]/[{1} seconds] = {2} reads with {3} e- each.".format(self.camera.cadence,self.camera.singleread, self.camera.cadence/self.camera.singleread, self.camera.read_noise))

		# calculate the variance due to read noise
		noise_variance = self.camera.cadence/self.camera.singleread*self.camera.read_noise**2

		# add noise into image
		self.image += np.sqrt(noise_variance)*np.random.randn(self.xsize, self.ysize)

		# update image header
		self.header['IREADNOI'] = ('True', 'read noise')
		self.show()


	def addSmear(self):
		'''Smear the image along the readout direction.'''
		self.speak("Adding readout smear.")
		self.speak("    assuming {0} second readout times on {1} second exposures.".format(self.camera.readouttime,self.camera.singleread))


		untouched = self.image + 0.0
		mean =np.mean(self.image, 1).reshape(self.image.shape[0],1)*self.ones()
		self.image += mean*self.camera.readouttime/self.camera.singleread

		self.note = 'readoutsmear'
		smearfilename = self.directory + self.note + '.fits'
		if not os.path.exists(smearfilename):
			self.writeToFITS(self.image - untouched,smearfilename)

		# update header
		self.header['ISMEAR'] = ('True', 'smearing during transer to frame store')
		self.show()

	def expose(self, plot=False, jitter=False, write=False, split=False, remake=False, smear=True, terse=False, cosmics='fancy', diffusion=False):
		'''Expose an image on this CCD.'''

		self.plot = plot
		self.terse = terse

		# create a blank image
		self.image = self.zeros()

		# jitter the camera, or at least update the
		if jitter:
			self.camera.jitter.jitter(self.camera.counter, header=self.header)

		# add stars to the image
		self.addStars(jitter=jitter, remake=remake)

		# add galaxies to the image
		self.addGalaxies()


		# add background to the image
		self.addBackgrounds()

		if write==False:
			stars = self.image + 0.0

		# add the photon noise from stars, galaxies, and backgrounds
		self.addPhotonNoise()

		# add cosmic rays to the image (after noise, because the *sub-Poisson* noise is already modeled with the Fano factor)
		cosmics = self.addCosmics(write=write, version=cosmics, diffusion=diffusion)

		# add smear from the finite frame transfer time
		if smear:
			self.addSmear()

		# create saturation bleed trails
		self.bleedSaturated()

		# add read noise, constant across detector
		self.addReadNoise()


		# finally, update the header for the image
		#self.populateHeader()

		# write the image
		if write:
			self.writeFinal()

		self.report("created image #{counter:07d} of {pos_string} with {cadence:.0f}s cadence".format(counter=self.camera.counter, pos_string=self.pos_string, cadence=self.camera.cadence))


		##### THIS WON"T WORK ON MULTIPLE CCD CAMERAS
		self.camera.advanceCounter()


		if write==False:
			return self.image, cosmics, stars

	def zip(self):
		'''Compress all the files in this directory.'''
		os.system('gzip -v {0}*.fits'.format(self.directory))



def gauss(x, y, xcenter, ycenter):
	rsquared = (x - xcenter)**2 + (y - ycenter)**2
	sigma = 1.0
	return np.exp(-0.5*rsquared/sigma**2)

def lorentz(x, y, xcenter, ycenter):
	rsquared = (x - xcenter)**2 + (y - ycenter)**2
	sigma = 1.0
	return 1.0/(rsquared/sigma**2 + 1.0)


def smoothedge( x,y,xcenter, ycenter,edge):
	rsquared = (x - xcenter)**2 + (y - ycenter)**2
	c = 1.0
	a = -c/edge**2
	return (a*rsquared + c)*(rsquared <= edge**2)
