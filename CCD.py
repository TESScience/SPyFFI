'''Tools to generate simulated TESS images.'''
# import necessary packages
import settings
from imports import *
from PSF import PSF
from Camera import Camera
import cosmical_original._cosmical
import cosmical_realistic._cosmical

# setup basic output options for this Python session
np.set_printoptions(threshold = 1e6, linewidth = 300)
plt.ion()




class CCD:
	def __init__(self, number=1, camera=None, subarray=None ):
		'''Turn on a TESS CCD:

			camera 			: the parent TESS camera for this CCD, required for *any* conversion from (RA, Dec) to (x,y)
			number = 1 		: which CCD in the field is it? (1-4, like mathematical quadrants)'''

		# keep track of where this image is
		self.number = number
		self.camera = camera

		if subarray is None:
			self.npix = 2048

		else:
			self.npix =
        self.npix = 4096 + self.gapinpixels

!!!! NEED TO TRANSFER ALL IMAGE GENERATION STUFF FROM CAMERA CLASS INTO HERE
		(give Camera an expose() method, which would expose on all four CCD's)

		self.xsize = self.camera.npix
		self.ysize = self.camera.npix
		self.gapinpixels = camera.gapinpixels
		self.subarraysize = (self.xsize - self.gapinpixels)/2
		self.xmin = 0
		self.xmax = self.xmin + self.xsize
		self.ymin = 0
		self.ymax = self.ymin + self.ysize
		self.xgapstart = (self.xsize - self.camera.gapinpixels)/2
		self.xgapend = (self.xsize + self.camera.gapinpixels)/2
		self.ygapstart, self.ygapend = self.xgapstart, self.xgapend
		self.pos_string = self.camera.pos_string()

		if self.camera.stamp is not None:
			self.prefix = settings.prefix + 'outputs/TESS_{stamp}pix_{pos}_{cadence:d}s'.format(stamp=self.camera.stamp, pos=self.pos_string, cadence=np.int(self.camera.cadence), counter=self.camera.counter).replace(' ','')
		else:
			self.prefix = settings.prefix + 'outputs/TESS_{pos}_{cadence:d}s'.format(pos=self.pos_string, cadence=np.int(self.camera.cadence), counter=self.camera.counter).replace(' ','')

		zachopy.utils.mkdir(self.prefix)
		self.directory = self.prefix + '/'
		self.image_note = ''
		self.nsubarrays = 32
		self.plot = False
		self.starcounter = 0
		self.nstars =0
		self.populateHeader()

	def getStars(self)
        if self.camera.testpattern:
            self.catalog = catalogs.TestPattern(size=self.npix*self.pixelscale)
        else:
            self.catalog = catalogs.UCAC4(ra=self.ra, dec=self.dec, radius=self.fov/np.sqrt(2)*1.01)

	def populateHeader(self, ccd=None):
		try:
			self.header
		except:
			self.header = astropy.io.fits.Header()
			self.header.extend(self.camera.header)
			self.header.extend(self.camera.psf.header)



		self.header['IMAGE'] = ''
		self.header['IMAGNOTE'] = ('', 'Details of this individual image')


		self.header['EXPTIME'] = (self.camera.cadence, '[s] exposure time ')
		self.header['NREADS'] = (np.round(self.camera.cadence/self.camera.singleread).astype(np.int), '# of reads summed')
		self.header['SUBEXPTI'] = (self.camera.singleread, '[s] exposure in a single subexposure')
		self.header['SATURATE'] = (self.camera.saturation*self.camera.cadence/self.camera.singleread, '[e-] saturation level in this image')
		self.header['READNOIS'] = (self.camera.read_noise, '[e-] read noise (per individual read)')
		self.header['READTIME'] = (self.camera.readouttime, '[s] time to transer to frame store')
		self.header['CCD'] = ('1,2,3,4', 'CCD number, using mathematical quadrants')
		self.header['NCCDS'] = (4, 'the number of CCDs going into this image')
		self.header['CCDSIZE'] = (self.subarraysize, '[pix] size of one CCD')
		self.setTime()

		self.header['INPUTS'] = ''
		self.header['INPUNOTE'] = ('', 'Ingredients for simulated images.')

	def setTime(self):
		self.header['COUNTER'] = self.camera.counter, '# of exposures since start, for this field'
		self.header['BJD'] = self.camera.counter*self.camera.cadence/24.0/60.0/60.0, '[day] mid-exposure time - 2457827.0 (e.g.)'






	def imageXY(self):
		x,y = np.meshgrid(np.arange(self.xsize) + self.xmin, np.arange(self.ysize) + self.ymin, indexing='ij')
		return x,y

	def imageAD(self):
		x,y = self.imageXY()
		ra, dec = self.camera.wcs.wcs_pix2world(x,y,1)
		return ra,dec

	def zeros(self):
		return np.zeros((self.xsize, self.ysize))

	def ones(self):
		return np.ones((self.xsize, self.ysize))

	def zodicalBackground(self,elat, elon):
			# from Josh and Peter's memo
			Vmax = 23.345
			DeltaV = 1.148
			b = np.abs(elat)
			V = Vmax  - DeltaV*((b-90.0)/90.0)**2
			assert((b < 90).all())
			return 10**(-0.4*(V-22.8))*(2.56e-3)*self.camera.effective_area*self.camera.pixel_solid_area

	def unresolvedBackground(self, glat, glon, complicated=False):
		# from Josh and Peter's memo
		flip = glon > 180
		if np.sum(flip):
			glon[flip] -= 360

		if complicated:
			L = np.abs(glon/180.0)
			B = np.abs(glat/90.0)
			a1 = 18.7
			a2 = 4.3
			a3 = 0.52
			a4 = 10.2
			a5 = 0.46
			a6 = -3.74
			I_surface_brightness = a1 + a2*(1-np.exp(-L/a3)) + a4*(1-np.exp(-B/a5)) + a6*np.sqrt(L*B)
		else:
			a0 = 18.9733
			a1 = 8.833
			a2 = 4.007
			a3 = 0.805
			I_surface_brightness = a0 + a1*(np.abs(glat)/40.0) + a2*(np.abs(glon)/180.0)**a3
		return 10**(-0.4*I_surface_brightness)*1.7e6*self.camera.effective_area*self.camera.pixel_solid_area

	def test(self):
		print "Creating test exposures."
		self.prefix = settings.prefix + 'outputs/TESS_test_image'
		self.createTestPattern()
		self.addStars()

	def createTestPattern(self):
		print "Replacing star catalog with a test pattern of stars."
		n = 25
		self.camera.starx, self.camera.stary = np.meshgrid( np.linspace(self.xmin, self.xmax),  np.linspace(self.ymin, self.ymax))
		self.camera.starmag = np.zeros_like(self.camera.starx) + 10.0

	def c1(self, image):

		return image[self.xsize - self.subarraysize:, self.ysize - self.subarraysize:]

	def c2(self, image):
		return image[0:self.subarraysize, self.ysize - self.subarraysize:]

	def c3(self, image):
		return image[0:self.subarraysize,0:self.subarraysize]

	def c4(self, image):
		return image[self.xsize - self.subarraysize:,0:self.subarraysize]

	def writeToFITS(self, image, path, split=False, savetype=np.float32):
		print "    writing to " + path


		fov = self.gapMask()*image
		if split:
			functions = [self.c1, self.c2, self.c3, self.c4]
			for q in [1,2,3,4]:
				image = functions[q-1](fov)
				filename = path.replace('.fits', '_c{0}.fits'.format(q))

				if q == 1:
					self.header['CRPIX1'] = -self.gapinpixels/2
					self.header['CRPIX2'] = -self.gapinpixels/2
				if q == 2:
					self.header['CRPIX1'] = self.subarraysize + self.gapinpixels/2
					self.header['CRPIX2'] = -self.gapinpixels/2
				if q == 3:
					self.header['CRPIX1'] = self.subarraysize + self.gapinpixels/2
					self.header['CRPIX2'] = self.subarraysize + self.gapinpixels/2
				if q == 4:
					self.header['CRPIX1'] = -self.gapinpixels/2
					self.header['CRPIX2'] =  self.subarraysize + self.gapinpixels/2

				self.header['CCD'] = (q, 'CCD number, using mathematical quadrants')
				self.header['NCCDS'] = (1, 'the number of CCDs going into this image')
				#assert((savetype(image) > 0).all())
				astropy.io.fits.PrimaryHDU(np.transpose(savetype(image)), header=self.header).writeto(filename, clobber=True)
		else:
			image = fov
			self.header['CRPIX1'] =  self.subarraysize + self.gapinpixels/2
			self.header['CRPIX2'] =  self.subarraysize + self.gapinpixels/2
			self.header['CCD'] = ('1,2,3,4', 'CCD number, using mathematical quadrants')
			self.header['NCCDS'] = (4, 'the number of CCDs going into this image')
			#assert((savetype(image) > 0).all())
			astropy.io.fits.PrimaryHDU(np.transpose(savetype(image)), header=self.header).writeto(path, clobber=True)


	def loadFromFITS(self, path):
		print "    trying to load image from ", path
		try:
			image = np.transpose(astropy.io.fits.open(path)[0].data)
			print "       ...success!"
			return image
		except:
			print "       ...failed."
			assert(False)


	def writeFinal(self, split=False, lean=True):
		self.note = 'final_{0:06.0f}'.format(self.camera.counter)
		finalfilename = self.directory + self.note + '.fits'
		print "Writing final TESS images."

		self.writeToFITS(self.image, finalfilename, split=split, savetype=np.int32)

		if lean == False:
			self.note = 'withoutbackground_{0:06.0f}'.format(self.camera.counter)
			self.writeToFITS(self.image - self.image_background, self.directory + self.note + '.fits')

	def setupPlots(self):
		fi, ax = plt.subplots(1,3,figsize=(30,10),sharex=True, sharey=True)
		self.ax_imstars = ax[0]
		self.ax_imperfect = ax[1]
		self.ax_imelectronics = ax[2]
		self.ax_imstars.set_aspect(1)
		plt.ion()
		plt.show()



	def relativeposition(self, x, y):
		scale = self.camera.npix/2.0
		return (x - scale)/scale, (y-scale)/scale

	def xy_from_center(self, x, y):
		scale = self.camera.npix/2.0
		return (x - scale), (y-scale)

	def addStar(self, x, y, mag, temp, verbose=True, plot=False):

		#xsubpixel = self.camera.psf.xsubgrid + np.long(x)
		#ysubpixel = self.camera.psf.ysubgrid + np.long(y)

		xrel, yrel = self.relativeposition(x,y)
		if verbose and self.starcounter % 10 == 0:
			print "         adding star at ({0}, {1}) with magnitude of {2}".format(x, y, mag)
			print "            assuming ({0}, {1}) relitve position and temperature of {2}".format(xrel, yrel, temp)
			print "               ({0}/{1})".format(self.starcounter,self.nstars)
		#subgrid_psf = self.camera.psf.psf(xrel, yrel, 5000.0)
		#subgrid_psf = lorentz(xsubpixel, ysubpixel, x, y)*smoothedge(xsubpixel, ysubpixel, x, y,self.gridpixels)
		#intrapixel = prnu(xsubpixel, ysubpixel)
		#assert(subgrid_psf.shape == intrapixel.shape)
		#normalized_psf = subgrid_psf*self.camera.photons(mag)*self.camera.cadence/np.sum(subgrid_psf)


		#binned = zachopy.utils.rebin_total(normalized_psf*intrapixel,2*self.psf.gridpixels,2*self.psf.gridpixels)

		unnormed = self.camera.psf.binnedpsf(x- self.camera.npix/2.0,y-self.camera.npix/2.0,temp)
		binned = unnormed*self.camera.cadence*self.camera.photons(mag)/np.sum(unnormed)

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


		xindex =  -self.camera.psf.dx_pixels + np.long(x) - self.xmin
		yindex =  -self.camera.psf.dy_pixels + np.long(y) - self.ymin

		ok = (xindex >= self.xmin) * (xindex < self.xsize) * (yindex >= self.ymin) * (yindex < self.ysize)
		print " a major kludge in addStar()! "
		ok = ok[:-1,:-1]
		self.starimage[xindex[ok], yindex[ok]] += binned[ok]

	def addStars(self, remake=False, jitter=False):
		try:
			assert(self.terse)
		except:
			print "Adding stars."

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
			except:
				self.camera.project()

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
					ok =    (self.camera.starx + self.camera.psf.dx_pixels_axis[-1] >= self.xmin) * \
							(self.camera.starx + self.camera.psf.dx_pixels_axis[0] <= self.xmax) * \
							(self.camera.stary + self.camera.psf.dy_pixels_axis[-1] >= self.ymin) * \
							(self.camera.stary + self.camera.psf.dy_pixels_axis[0] <= self.ymax) * \
							(self.camera.starmag < threshold)*(self.camera.starmag >= min)
					x = self.camera.starx[ok]
					y = self.camera.stary[ok]
					mag = self.camera.starmag[ok]
					temp = self.camera.startemp[ok]

					if np.sum(ok) > 0:
						self.nstars += np.sum(ok)
						for i in range(len(x)):
							self.addStar(x[i], y[i], mag[i], temp[i])
							self.starcounter += 1

						if jitter == False:
							self.writeToFITS(self.starimage, starsfilename)
		self.image += self.starimage
		self.header['ISTARS'] = ('True', 'stars from UCAC4')
		if jitter:
			self.header['IJITTER'] = ('True', 'spacecraft jitter, motion between images')
		return self.starimage

	def addGalaxies(self):
		pass
		# looks like I should use http://vizier.cfa.harvard.edu/viz-bin/Cat?VII/155 for a source catalog?


	def addCosmicsAl(self, split=False, gradient=False, write=False, version='fancy', diffusion=False):

		try:
			assert(self.terse)
		except:
			print "Adding cosmic rays."
		rate = 5.0 # cts cm^-2 s^-1
		self.note = 'cosmics_{0:.1f}persecond_{1}seconds_{2:06.0f}'.format(rate,self.camera.cadence, self.camera.counter).replace('.', 'p')
		cosmicsfilename = self.directory + self.note + '.fits'


		if gradient:
			smallexptime = self.camera.cadence
			bigexptime = self.camera.cadence*2
		else:
			smallexptime = self.camera.cadence*1.5
			bigexptime = smallexptime


		if version == 'original':
			# you'll also need to change the import statement up at the start
			nexpected = (self.camera.physicalpixelsize*self.camera.npix)**2*(0.5*smallexptime + 0.5*bigexptime)*rate
			ndrawn = np.random.poisson(nexpected)
			image = np.transpose(cosmical_original._cosmical.cosmical(smallexptime, bigexptime, ndrawn, self.camera.npix, self.camera.npix))
			if diffusion:
				kernal = np.array([	[0.0034, 0.0516, 0.0034],
									[0.0516, 0.7798, 0.0516],
									[0.0034, 0.0516, 0.0034]])
				image = scipy.signal.convolve2d(image, kernal, mode='same')
				print "convolved with {0}".format(kernal)
		elif version == 'fancy':
			intdiffusion=np.int(diffusion)
			image = np.transpose(cosmical_realistic._cosmical.cosmical(rate, smallexptime, bigexptime, self.camera.npix, self.camera.npix, intdiffusion))
			print '=========================='
			print smallexptime, bigexptime

		if write:
			self.writeToFITS(image, cosmicsfilename, split=split)

		self.image += image
		self.header['ICOSMICS'] = ('True', 'cosmic rays injected')
		return image

	def addCosmics(self, split=False):
		try:
			assert(self.terse)
		except:
			print "Adding cosmic rays."


		rate = 5 # cts cm^-2 s^-1
		depth = self.camera.physicalpixeldepth

		# Roland says a typical 10 pixel cosmic ray deposits 15000 electrons per trace
		#electron_deposition = 15000.0/np.sqrt((10*self.camera.physicalpixelsize)**2 + self.camera.physicalpixeldepth**2)
		# assume 100 electrons deposited per micron of silicon traveled
		electron_deposition = 100.0*1e4
		detectorarea = (self.camera.npix*self.camera.physicalpixelsize)**2		# cm^2
		ncosmics = 2*np.int(np.round(self.camera.cadence*detectorarea*rate))
		print "    {0} cosmic rays in a {1} second exposure".format(ncosmics, self.camera.cadence)
		print "    assuming {0} e- deposited per micron (= {1} e- per 15 microns = {2} e- per {3} microns)".format(electron_deposition/1e4, electron_deposition*15/1e4,  electron_deposition*depth,  depth*1e4,)

		self.note = 'cosmics_{0:.1f}persecond_{1:.1f}epermicron_{2}seconds_{3:06.0f}'.format(rate,electron_deposition/1e4,self.camera.cadence, self.camera.counter).replace('.', 'p')
		cosmicsfilename = self.directory + self.note + '.fits'

		try:
			image_cosmics = self.loadFromFITS(cosmicsfilename)
		except:
			x0, y0, vx, vy, vz = np.random.random((5,ncosmics))
			x0 *= self.xsize
			y0 *= self.ysize
			vx -= 0.5
			vy -= 0.5
			vs = np.sqrt(vx**2 + vy**2 + vz**2)
			vr = np.sqrt(vx**2 + vy**2)
			image_cosmics = self.zeros()
			for i in range(ncosmics):

				z = np.linspace(0, depth/self.camera.physicalpixelsize, 1+ np.maximum(depth/self.camera.physicalpixelsize/(vz[i]/vr[i]),1))
				x = np.round(x0[i] + vx[i]/vz[i]*z).astype(np.int)
				y = np.round(y0[i] + vy[i]/vz[i]*z).astype(np.int)
				silicon = (vs[i]/vz[i]*(z[-1] - z[0]))*self.camera.physicalpixelsize	# microns
				pixels = np.maximum(np.int(np.round(vr[i]/vz[i]*(z[-1] - z[0]) )),1)
				#if silicon*1e4 > 100:
				#	print z
				assert(pixels > 0)
				brightness = electron_deposition*silicon/pixels
				ok = (x >= self.xmin)*(x < self.xmax)*(y >= self.ymin)*(y < self.ymax)
				image_cosmics[x[ok],y[ok]] += brightness
				if i % 10000 == 0:
					print "  {0:12}/{1}".format(i,ncosmics) + " -- cosmic travelled {0:7.1f} microns through {1:6d} pixels  ".format(silicon*1e4, pixels) + "; the path length was {0:6.1f} microns per pixel, depositing {1:8.1f} electrons per pixel".format(silicon/pixels*1e4, brightness)
				#print "       x = {0}\n       y = {1}\n       r = {2}".format(x,y,np.sqrt(x**2 + y**2))

			self.writeToFITS(image_cosmics, cosmicsfilename, split=split)
		self.image += image_cosmics
		self.header['ICOSMICS'] = ('True', 'cosmic rays injected')

	def bleedSaturated(self, plot=False):
		try:
			assert(self.terse)
		except:
			print "Bleeding saturated pixels."
		untouched = self.image + 0.0
		original = np.sum(self.image)
		# set the saturation limit based on the number of individual reads included
		saturation_limit =  self.camera.saturation*self.camera.cadence/self.camera.singleread
		stilloversaturated = True

		count =0
		# keep looping until all saturation problems are gone
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
			try:
				assert(self.terse)
			except:
				print "      on pass #{0} through saturation filter,  \n        the max saturation fraction is {1} \n        and the flux change over entire image is {2} electrons".format(count, np.max(self.image)/saturation_limit, np.sum(self.image) - original)
			# KLUDGE!
			stilloversaturated = (self.image > saturation_limit).any() and count < 10
		self.note = 'saturation_{0}K'.format(self.camera.saturation).replace('.', 'p')

		saturationfilename = self.directory + self.note + '.fits'
		if not os.path.exists(saturationfilename):
			self.writeToFITS(self.image - untouched,saturationfilename)
		self.header['ISATURAT'] = ('True', 'bleed trails for saturated pixels')

	def addBackgrounds(self):

		# set up filenames for saving background, if need be
		self.note = 'backgrounds'
		backgroundsfilename = self.directory + self.note + '.fits'
		try:
			assert(self.terse)
		except:
			print "Adding backgrounds."

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
				ra, dec = self.imageAD()
				gal = astropy.coordinates.ICRS(ra=ra, dec=dec, unit=(astropy.units.deg,astropy.units.deg)).galactic
				elon, elat = zachopy.borrowed.crossfield.euler(ra, dec, select=3)

				# add the zodiacal light, using the simple model from Josh and Peter on the TESS wiki
				print "   including smooth model for zodiacal light."
				self.backgroundimage += self.zodicalBackground(elat, elon)*self.camera.cadence
				self.header['IZODIACA'] = ('True', 'zodiacal light, treated as smooth')

				# add unresolved background light, using the simple model from Josh and Peter on the TESS wiki
				print "   including smooth model for unresolved stars in the Galaxy."
				self.backgroundimage += self.unresolvedBackground(gal.b.degree, gal.l.degree)*self.camera.cadence
				self.header['IUNRESOL'] = ('True', 'unresolved stars, treated as smooth background')

				# write the image, so it can just be loaded easily next time
				self.writeToFITS(self.backgroundimage, backgroundsfilename)

		# add the background image to the total image
		self.image += self.backgroundimage

	def addPhotonNoise(self):
		try:
			assert(self.terse)
		except:
			print "Adding photon noise."
			print "    sqrt(photons from stars and various backgrounds)"
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

	def addReadNoise(self):
		try:
			assert(self.terse)
		except:
			print "Adding read noise."
			print "    = quadrature sum of [{0} seconds]/[{1} seconds] = {2} reads with {3} e- each.".format(self.camera.cadence,self.camera.singleread, self.camera.cadence/self.camera.singleread, self.camera.read_noise)
		noise_variance = self.camera.cadence/self.camera.singleread*self.camera.read_noise**2
		self.image += np.sqrt(noise_variance)*np.random.randn(self.xsize, self.ysize)
		self.note = 'readnoise'

		noisefilename = self.directory + self.note + '.fits'
		if not os.path.exists(noisefilename):
			self.writeToFITS(np.sqrt(noise_variance*np.ones_like(self.image)),noisefilename)
		self.header['IREADNOI'] = ('True', 'read noise')

	def gapMask(self):
		mask = self.ones()
		mask[self.xgapstart:self.xgapend,:] = 0
		mask[:, self.ygapstart:self.ygapend] = 0
		return mask

	def addGaps(self):
		try:
			assert(self.terse)
		except:
			print "Removing the pixels in the CCD gaps."
		self.image[self.xgapstart:self.xgapend,:] = 0
		self.image[:, self.ygapstart:self.ygapend] = 0


	def addSmear(self):
		try:
			assert(self.terse)
		except:
			print "Adding readout smear."
			print "    assuming {0} second readout times on {1} second exposures.".format(self.camera.readouttime,self.camera.singleread)
		height = 2048
		untouched = self.image + 0.0
		chunk = self.image[:,:2048]
		mean =np.mean(chunk, 1).reshape(chunk.shape[0],1)*np.ones_like(chunk)
		chunk += mean*self.camera.readouttime/self.camera.singleread


		chunk = self.image[:,-2048:]
		mean =np.mean(chunk, 1).reshape(chunk.shape[0],1)*np.ones_like(chunk)
		chunk += mean*self.camera.readouttime/self.camera.singleread



		self.note = 'readoutsmear'
		smearfilename = self.directory + self.note + '.fits'
		if not os.path.exists(smearfilename):
			self.writeToFITS(self.image - untouched,smearfilename)

		self.header['ISMEAR'] = ('True', 'smearing during transer to frame store')


	def advanceCounter(self):
		self.camera.counter +=1
		'''try:
			counter = 0
			files = glob.glob(self.directory + 'final*.fits')
			for f in files:
				number = np.in(f.split('/')[-1].split('.')[0][-4:])
				counter = np.maximum(number, counter)
			self.camera.counter = counter
		except:
			self.camera.counter = 0'''

	def expose(self, plot=False, jitter=False, write=False, split=False, remake=False, smear=True, terse=False, cosmics='fancy', diffusion=False):
		self.plot = plot
		self.terse = terse

		# create a blank image
		self.image = self.zeros()

		# jitter the camera, or at least update the
		if jitter:
			self.camera.jitter(header=self.header)

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
		cosmics = self.addCosmicsAl(split=split, write=write, version=cosmics, diffusion=diffusion)

		# add smear from the finite frame transfer time
		if smear:
			self.addSmear()

		# create saturation bleed trails
		self.bleedSaturated()

		# add read noise, constant across detector
		self.addReadNoise()

		# zero out the gaps in the CCD
		self.addGaps()

		# finally, update the header for the image
		#self.populateHeader()

		# write the image
		if write:
			self.writeFinal(split=split)

		if terse:
			print "[tess] created image #{counter:07d} of {pos_string} with {cadence:.0f}s cadence".format(counter=self.camera.counter, pos_string=self.camera.pos_string(), cadence=self.camera.cadence)
		else:
			print "-------_-------_-------_-------_-------_-------_-------_-------"
			print "CREATING IMAGE #{0:07d}".format(self.camera.counter)
			print "----_----_----_----_----_----_----_----_----_----_----_----_---"

		self.advanceCounter()

		#
		if write==False:
			return self.image, cosmics, stars

	def zip(self):
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


#T = Camera()
#T.point(15,60)	# - Cassiopeia
#T.point(82, 1) # Orion
#T.jitter()
#T.point(285,35) #- Lyra
#T.point(280,30)	# - Sagittarius
#T.project()
#I = Image(T)
#I.populateHeader()
#I.expose()
#I.raster()
