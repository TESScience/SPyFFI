'''Tools to generate simulated TESS images.'''
# import necessary packages
import settings
from imports import *
import Cosmics, Catalogs

# setup basic output options for this Python session
np.set_printoptions(threshold = 1e6, linewidth = 300)

zipsuffix = ''

# define mapping between CCD number and quadrant
quadrants = {1:(1,1), 2:(-1,1), 3:(-1,-1), 4:(1,-1), 0:None}


def dndmag(m):
	# a rough fit to a random (small) field
	f = np.exp(-3.293 + 0.697*m)
	return f

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


		# start populating the image header (seems like we can't do this until we're sure camera has pointed somewhere)
		#self.populateHeader()

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

	def photons(self, mag):
		'''Use magnitude to calculate photons per second that will be recorded on a CCD.'''

		# this doesn't work great for M dwarfs, need to include multiband information at some point
		return self.camera.effective_area*10**(-0.4*mag)*1.7e6

	def populateHeader(self, ccd=None):
		'''Populate the header structure with information about the CCD image, and the inputs.'''

		# create an empty header
		try:
			self.camera.header
		except:
			self.camera.populateHeader()#astropy.io.fits.Header()
		self.header = self.camera.header


		# fill it with some CCD details
		self.header['CCD'] = ''
		self.header['CCDNOTE'] = ('', 'Details of this individual image')
		self.header['EXPTIME'] = (self.camera.cadence, '[s] exposure time ')
		self.header['NREADS'] = (np.round(self.camera.cadence/self.camera.singleread).astype(np.int), '# of reads summed')
		self.header['SUBEXPTI'] = (self.camera.singleread, '[s] exposure in a single subexposure')
		self.header['SATURATE'] = (self.camera.saturation*self.camera.cadence/self.camera.singleread, '[e-] saturation level in this image')
		self.header['READNOIS'] = (self.camera.read_noise, '[e-] read noise (per individual read)')
		self.header['READTIME'] = (self.camera.readouttime, '[s] time to transer to frame store')
		self.header['CCDNUM'] = (self.number, 'CCD number (1,2,3,4 or 0=fake subarray)')
		self.header['CCDSIZE'] = (self.npix, '[pix] size of one CCD')

		# fill in the timestamp for this CCD image
		self.setTime()

	def addInputLabels(self):
		try:
			self.header['INPUTS']
		except KeyError:
			# leave space to talk about the inputs to this image
			self.header['INPUTS'] = ''
			self.header['INPUNOTE'] = ('', 'Ingredients for simulated images.')

	def setTime(self):
		'''Based on the camera counter, apply a timestamp to this image.'''

		# add time to the image (in a very simplistic way -- no accounting for the spacecraft orbit)
		self.header['COUNTER'] = self.camera.counter, '# of exposures since start, for this field'
		self.bjd0 = 2457827.0
		self.bjd = self.bjd0 + self.camera.counter*self.camera.cadence/24.0/60.0/60.0
		self.bjdantisun = self.bjd0 + 13.7
		self.header['BJD0'] = self.bjd0, '[day] base time subtracted from all BJD'
		self.header['BJD'] = self.bjd - self.bjd0, '[day] mid-exposure time - BJD0'
		self.header['ANTISUN'] = self.bjdantisun - self.bjd0, '[day] time of antisun - BJD0'
		self.epoch = (self.bjd - 2451544.5)/365.25 + 2000.0
		self.header['EPOCH'] = self.epoch, '[years] epoch of mid-exposure time'

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
		'''Calcaulte the zodiacal background at a given celestial (lat, long).'''

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


	def writeToFITS(self, image, path, split=False, savetype=np.float32, zip=False):
		'''General FITS writer for this CCD.'''

		# print status
		self.speak('saving {0}x{0} image to {1} with type {2}'.format(image.shape[0],path,savetype.__name__))
		cr = self.camera.cartographer.point(0.0, 0.0, 'focalxy')
		x, y = cr.ccdxy.tuple

		#quad = quadrants[self.number]
		#if quad is not None:
		#	offset = self.npix
		#	x += quad[0]*

		# modify the camera's WCS, based on the CCD number
		self.header['CRPIX1'] = x + 0.5
		self.header['CRPIX2'] = y + 0.5

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

	@property
	def fileidentifier(self):
		return '{pos_string}_{name}_{counter:06.0f}'.format(pos_string=self.camera.pos_string(), name=self.name, counter=self.camera.counter)

	def writeFinal(self, lean=True):
		'''Write the final image from this CCD.'''

		#self.header.extend(self.camera.header)
		#self.header.extend(self.camera.psf.header)

		# make filename for this image
		self.note = 'simulated_'+self.fileidentifier
		finalfilename = self.directory + self.note + '.fits' + zipsuffix

		# write the image to FITS
		self.speak('saving final TESS image')
		self.writeToFITS(self.image, finalfilename, savetype=np.int32)

		# optionally, write some other outputs too!
		if lean == False:
			self.note = 'withoutbackground_{0:06.0f}'.format(self.camera.counter)
			self.writeToFITS(self.image - self.image_background, self.directory + self.note + '.fits')

	def createStamps(self, x, y, radius=5):

		# create an empty array
		mask = self.zeros()

		# create a (currently fixed size and shape) aperture
		diameter = radius*2 + 1
		xgrid,ygrid = np.meshgrid(np.arange(-radius,radius+1),np.arange(-radius,radius+1))
		aperture = np.zeros_like(xgrid)

		#make a circular aperture
		r = np.sqrt(xgrid**2 + ygrid**2)
		aperture[r <= radius] = 1


		# loop over stars
		for i in range(len(x)):
			thisx = np.minimum(np.maximum(np.round(x[i]).astype(np.int) + xgrid, 0), self.xsize-1)
			thisy = np.minimum(np.maximum(np.round(y[i]).astype(np.int) + ygrid, 0), self.ysize-1)
			mask[thisy,thisx] += aperture

		self.stampimage = mask
		self.note = 'stampdefinition'
		stampfilename = self.directory + self.note + '.fits'
		self.writeToFITS(self.stampimage,stampfilename)

	def stampify(self):
		self.image *= self.stampimage != 0

	def setupCatalog(self, write=True):
		'''setup the CCD's catalog, by creating it from the camera and then trimming'''

		self.speak("setting up this CCD's starmap")

		# make sure the camera has a catalog defined
		try:
			self.camera.catalog
			self.speak("the camera already had a catalog of {0} elements defined; using it!".format(len(self.camera.catalog.ra)))
		except AttributeError:
			self.speak("populating a new catalog for the camera")
			self.camera.populateCatalog()

		# we want to make a trimmed catalog for this CCD.
		# first figure out which ones are on the CCD

		# pull out positions, magnitudes, and temperatures at the time of the first exposure
		self.speak('taking an intial snapshot at {0} = {1}'.format(self.bjd, self.epoch))
		ras, decs, tmag, temperatures = self.camera.catalog.snapshot(self.bjd,
								exptime=self.camera.cadence/60.0/60.0/24.0)
		self.speak('  done!')
		assert(ras.shape == tmag.shape)

		# assign the cartrographer's CCD to this one
		self.camera.cartographer.ccd = self

		# create coordinate object for the stars
		stars = self.camera.cartographer.point(ras, decs, 'celestial')
		x,y = stars.ccdxy.tuple

		# trim everything down to only those stars that could be relevant for this ccd
		buffer = 10
		ok = (x > -buffer) & (x < self.xsize + buffer) & (y > -buffer) & (y < self.ysize + buffer)

		# trim to postage stamps, if desired
		if self.camera.stamps is not None:

			# select (randomly?) some target stars
			np.random.seed(self.number)


			# currently a kludge -- eventually replace with actual target selection
			'''# this is how the weights were estimated, from a small field
			hy, edges = np.histogram(tmag)
			hx = 0.5*(edges[1:] + edges[:-1])
			plt.figure()
			plt.plot(hy, hx)
			np.polyfit(hx, np.log(hy), 1)'''


			# weight stars inversely to their abundance, to give roughly uniform distribution of magnitudes



			# start with targets with centers on the chip
			onccd = (np.round(x) > 0) & \
					(np.round(x) < self.xsize) & \
					(np.round(y) > 0) & \
					(np.round(y) < self.ysize)
			weights = (1.0/dndmag( self.camera.catalog.tmag)*( self.camera.catalog.tmag >= 6)*( self.camera.catalog.tmag <= 16))[onccd]
			weights /= np.sum(weights)
		 	itargets = np.random.choice(onccd.nonzero()[0],
										size=np.minimum(self.camera.stamps, np.sum(weights != 0)),
										replace=False,
										p=weights)


			# pull out the target stars, and write them to disk
			self.targetstarcatalog = Catalogs.Trimmed(self.camera.catalog, itargets)
			targetsoutfile = self.directory + 'postagestamptargets_{pos}_{name}.txt'.format(pos=self.pos_string, name=self.name)
			self.targetstarcatalog.writeProjected(ccd=self, outfile=targetsoutfile)

			#for t in range(len(self.targetstarcatalog.ra)):
			#	i = (self.camera.catalog.ra == self.targetstarcatalog.ra[t])*(self.camera.catalog.dec == self.targetstarcatalog.dec[t])
			#	assert(np.array(self.camera.catalog.lightcurvecodes)[i] == np.array(self.targetstarcatalog.lightcurvecodes)[t])

			# make a stamp image centered on target stars
			self.createStamps(x[itargets], y[itargets])

			# test which stars fall in the mask
			ix, iy = np.round(x[onccd]).astype(np.int), np.round(y[onccd]).astype(np.int)

			# mask out those stars that fall outside stamps (or detector)
			ok *= onccd
			ok[onccd] *= self.stampimage[iy,ix]

		# keep track of whether this is a stamp catalog or not
		self.stamps = self.camera.stamps

		# create the CCD catalog
		self.catalog = Catalogs.Trimmed(self.camera.catalog, ok)

		# write the catalog out to a text file
		outfile = self.directory + 'catalog_{pos}_{name}.txt'.format(pos=self.pos_string, name=self.name)
		self.catalog.writeProjected(ccd=self, outfile=outfile)

	def projectCatalog(self, write=True):
		'''Create using the camera's star catalog, and project stars using this CCD.'''
		self.speak('projecting the starmap onto CCD')

		# make sure the camera has a catalog defined
		try:
			self.catalog
			assert(self.stamps == self.camera.stamps)
		except (AttributeError,AssertionError):
			self.setupCatalog()


		# pull out positions, magnitudes, and temperatures
		self.speak('taking a snapshot at {0} = {1}'.format(self.bjd, self.epoch))
		ras, decs, tmag, temperatures = self.catalog.snapshot(self.bjd, exptime=self.camera.cadence/60.0/60.0/24.0)
		self.speak('  done!')
		assert(ras.shape == tmag.shape)
		self.camera.cartographer.ccd = self

		# create coordinate object for the stars
		stars = self.camera.cartographer.point(ras, decs, 'celestial')
		x,y = stars.ccdxy.tuple

		# apply differential velocity aberration, based on the time offset from antisun
		dt = self.bjd - self.bjdantisun
		dx,dy = self.aberrations(stars, dt)
		fieldcenter = self.camera.cartographer.point(0,0,'focalxy')
		centerx,centery = self.aberrations(fieldcenter,dt)
		x+=dx-np.mean(dx)
		y+=dy-np.mean(dy)

		# trim everything down to only those stars that could be relevant for this camera
		buffer = 10
		ok = (x > -buffer) & (x < self.xsize + buffer) & (y > -buffer) & (y < self.ysize + buffer)

		# assign this CCD's stars
		self.starx = x[ok]
		self.stary = y[ok]
		self.starmag = np.array(tmag)[ok]
		self.startemp = np.array(temperatures)[ok]
		self.starlc = np.array(self.camera.catalog.lightcurvecodes)[ok]
		self.starbasemag = np.array(self.camera.catalog.tmag)[ok]
		# keep track of which CCD we projected onto
		self.starsareon = self.name

	def aberrations(self, stars, dt):

		# make sure the cartographer is defined
		try:
			assert(self.aberrator.ccd == self)
		except AttributeError:
			self.aberrator = Aberrator(self.camera.cartographer)
			self.header['ABERRATE'] = ''
			self.header['AB_NOTE'] = '', 'Velocity aberration ingredients.'
			self.header['AB_DEF'] = 'd?=BETA*cos(L-FCLON+DLON)*AB?FUNC(x,y)-FCD?', "[L=stars' ec. lon.]"

			for pix in ['x','y']:
				self.header['ABD{0}FUNC'.format(pix.upper())] = 'C + CX*x + CY*y + CXX*x**2 + CYY*y**2 + CXY*x*y'#'{0:+.3f}{1:+.3f}*x{2:+.3f}*y{3:+.3f}*x**2{4:+.3f}*y**2{5:+.3f}*x*y'.format(*self.aberrator.coefs[pix])
				self.header['ABD{0}_C'.format(pix.upper())] = self.aberrator.coefs[pix][0]
				self.header['ABD{0}_CX'.format(pix.upper())] =  self.aberrator.coefs[pix][1]
				self.header['ABD{0}_CY'.format(pix.upper())] =  self.aberrator.coefs[pix][2]
				self.header['ABD{0}_CXX'.format(pix.upper())] =  self.aberrator.coefs[pix][3]
				self.header['ABD{0}_CYY'.format(pix.upper())] =  self.aberrator.coefs[pix][4]
				self.header['ABD{0}_CXY'.format(pix.upper())] = self.aberrator.coefs[pix][5]

		# sign will definitely be wrong on this
		beta = 29.8*zachopy.units.km/zachopy.units.c # unitless (radians)
		if self.camera.warpspaceandtime:
			warp=self.camera.warpspaceandtime
			beta/=warp
			self.header['AB_WARP'] = 'speed of light is {0}X what it should be'.format(warp), '(for testing)'
		self.header['AB_BETA'] = beta, '[radians] v/c (from orbit tangential velocity)'

		# what is the celestial longitude of the field center?
		fieldcenter = self.camera.cartographer.point(0,0,type='focalxy')
		self.header['AB_FCLON'] = fieldcenter.ecliptic.elon, '[deg] ecliptic lon. of focal plane center'

		# how much is each star offset from antisun, at this time?
		dtheta = 360.0*dt/365.25
		theta = stars.ecliptic.elon - fieldcenter.ecliptic.elon + dtheta
		self.header['AB_DLON'] = dtheta, '[deg] motion of Earth (in ecliptic lon.)'

		# calculate the current positions and the nudge of celestial longitude
		x, y = stars.ccdxy.tuple	# pixels
		delon = beta*np.cos(theta*np.pi/180.0)*180/np.pi # degrees

		# calculate the aberration at the center of the camera FOV
		fcx, fcy = fieldcenter.ccdxy.tuple
		dfcelon = beta*np.cos(dtheta*np.pi/180.0)*180/np.pi
		fcdx, fcdy = self.aberrator.derivatives['x'](fcx, fcy)*dfcelon, self.aberrator.derivatives['y'](fcx,fcy)*dfcelon
		self.header['AB_FCDX'] = fcdx, '[pix] dx of FOV center'
		self.header['AB_FCDY'] = fcdy, '[pix] dy of FOV center'

		self.delon = delon # for testing
		self.fcdx, self.fcdy = fcdx, fcdy

		# calculate the absolute velocity aberration
		# the linear plane model is probably going to break down at the pole!
		vax, vay = self.aberrator.derivatives['x'](x,y)*delon, self.aberrator.derivatives['y'](x,y)*delon

		# subtract the field centers
		dx, dy = vax-fcdx, vay-fcdy
		return dx, dy





	def addStar(self, ccdx, ccdy, mag, temp, verbose=False, plot=False):
		'''Add one star to an image, given position, magnitude, and effective temperature.'''


		#self.speak("adding stars at ({0}, {1}) with magnitude of {2}".format(ccdx, ccdy, mag))
		#self.speak(" ({0}/{1})".format(self.starcounter,self.nstars))

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
		self.starcounter = 0
		self.nstars =0
		if jitter:
			remake=True
		self.camera.cartographer.pithy = True
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

			# propagate proper motions and project onto the detector
			self.projectCatalog()


			#for threshold in magnitude_thresholds:
			if True:
				self.speak('adding {0} stars to image'.format(len(self.starx)))
				threshold = np.max(magnitude_thresholds)
				# define a filename for this magnitude range
				self.note = 'starsbrighterthan{0:02d}'.format(threshold)
				starsfilename = self.directory + self.note + '.fits'

				# load the existing stellar image, if possible
				try:
					assert(remake==False)
					self.starimage = self.loadFromFITS(starsfilename)
				except:
					# if this is the smallest threshold, include all the stars brighter than it
					#if threshold == np.min(magnitude_thresholds):
					minimum = -100
					#else:
					#	minimum = threshold - dthreshold
					# pick the stars to add to the image on this pass through
					ok =    (self.starx + self.camera.psf.dx_pixels_axis[-1] >= self.xmin) * \
							(self.starx + self.camera.psf.dx_pixels_axis[0] <= self.xmax) * \
							(self.stary + self.camera.psf.dy_pixels_axis[-1] >= self.ymin) * \
							(self.stary + self.camera.psf.dy_pixels_axis[0] <= self.ymax) * \
							(self.starmag < threshold)*(self.starmag >= minimum)
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
		self.addInputLabels()
		if self.camera.testpattern:
			self.header['ISTARS'] = ('True', 'stars from a test pattern')
		else:
			self.header['ISTARS'] = ('True', 'stars from UCAC4')

		if jitter:
			self.header['IJITTER'] = ('True', 'spacecraft jitter, motion between images')
		else:
			self.header['IJITTER'] = ('False', 'no spacecraft jitter apply')

		if self.camera.aberrate:
			self.header['IVELABER'] = ('True'), 'differential velocity aberration applied'
		else:
			self.header['IVELABER'] = ('False'), 'no differential velocity aberration'
		return self.starimage

	def addGalaxies(self):
		pass
		# looks like I should use http://vizier.cfa.harvard.edu/viz-bin/Cat?VII/155 for a source catalog?


	def addCosmics(self, gradient=False, version='fancy', diffusion=False, write=False, rate=5.0, correctcosmics=True):
		'''Add cosmic rays to image.'''

		# print update
		self.speak('Adding cosmic rays')

		# filenames, in case saving is required


		# use Al's code to generate cosmic ray image of the correct size
		image = Cosmics.cosmicImage(exptime=self.camera.cadence, size=self.npix,
									gradient=gradient, version=version, diffusion=diffusion, rate=rate)

		# (optionally), write cosmic ray image
		if write:
			self.note = 'cosmics_'+self.fileidentifier
			cosmicsfilename = self.directory + self.note + '.fits' + zipsuffix
			self.writeToFITS(image, cosmicsfilename)

		# add the cosmics into the running image
		self.addInputLabels()
		if (correctcosmics == False) or self.camera.cadence <= 2:
			self.image += image
			self.header['ICOSMICS'] = ('True', 'cosmic rays injected')
		else:
			self.header['ICOSMICS'] = ('False', 'cosmic rays injected')
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

			# keep track of iterations, to prevent infinite loops!
			count +=1

			# identify saturated pixels
			oversaturated = self.image > saturation_limit
			saturated = self.image >= saturation_limit

			# loop over columns, treating each separately
			for x in range(self.image.shape[1]):

				# identify continuous saturated regions
				regions, nregions = scipy.ndimage.measurements.label(saturated[:,x])
				if nregions > 0:
					for i in np.arange(nregions)+1:
						y = (regions == i).nonzero()[0]
						if oversaturated[y,x].any():
							#self.speak('')
							#self.speak('in column {0}'.format(x))
							# in this saturated region, how much flux needs to be redistributed?
							fluxtodistribute = np.sum(self.image[y,x])
							#self.speak('must distribute {0} electrons'.format(fluxtodistribute))
							# how many pixels would this correspond to? (including the original pixels)
							npixels = (fluxtodistribute/saturation_limit)
							#self.speak('which we could do over {0} pixels'.format(npixels))

							# find the center of the saturated region
							center = np.mean(y)

							# find how far we away from center we can totally saturate pixel
							grow = (npixels-1.0)/2.0
							indices = np.arange(np.maximum(np.ceil(center - grow).astype(np.int), 0),np.minimum(np.floor(center + grow).astype(np.int), self.image.shape[0]-1)+1)
							#self.speak('with indices of {0}'.format(indices))

							assert(y[0] in indices)

							# record the flux we're starting with in this region
							existingflux = np.sum(self.image[indices,x])

							# saturate the pixels needed
							self.image[indices,x] = saturation_limit
							leftoverflux = existingflux - indices.shape[0]*saturation_limit
							#self.speak('leaving {0} behind'.format(leftoverflux))

							'''if leftoverflux > 0:
								leftedge = indices.min() -1
								rightedge = indices.max() +1
							else:'''
							leftedge = indices.min()-1
							rightedge = indices.max()+1
							try:
								try:
									self.image[leftedge,x] += leftoverflux/2.0
								except:
									self.image[rightedge,x] += leftoverflux/2.0
								try:
									self.image[rightedge,x] += leftoverflux/2.0
								except:
									self.image[leftedge,x] += leftoverflux/2.0
							except:
								self.speak("this star seems to saturate the entire detector!")
			self.speak("on pass #{0} through saturation filter: the max saturation fraction is {1}; and the flux change over entire image is {2} electrons".format(count, np.max(self.image)/saturation_limit, np.sum(self.image) - original))
			#self.display=True
			#self.show()
			#self.display=False
			#self.input('saturation?')

			# KLUDGE to prevent endless loops
			stilloversaturated = (self.image > saturation_limit).any() and count < 10

		self.note = 'saturation_{0}K'.format(self.camera.saturation).replace('.', 'p')
		saturationfilename = self.directory + self.note + '.fits'
		if not os.path.exists(saturationfilename):
			self.writeToFITS(self.image - untouched,saturationfilename)

		# update image header
		self.addInputLabels()
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
				# define coordinates (equatorial, Galactic, celestial) at every pixel in the image
				ra, dec = self.pixels().celestial.tuple
				elon, elat = self.pixels().celestial.tuple
				glon, glat = self.pixels().galactic.tuple

				# add the zodiacal light, using the simple model from Josh and Peter on the TESS wiki
				self.speak( "   including smooth model for zodiacal light.")
				elon, elat
				self.backgroundimage += self.zodicalBackground(elon, elat)*self.camera.cadence
				self.addInputLabels()
				self.header['IZODIACA'] = ('True', 'zodiacal light, treated as smooth')

				# add unresolved background light, using the simple model from Josh and Peter on the TESS wiki
				self.speak( "   including smooth model for unresolved stars in the Galaxy.")
				self.backgroundimage += self.unresolvedBackground(glon, glat)*self.camera.cadence
				self.header['IUNRESOL'] = ('True', 'unresolved stars, treated as smooth background')

				# write the image, so it can just be loaded easily next time
				self.writeToFITS(self.backgroundimage, backgroundsfilename)

		# add the background image to the total image
		self.image += self.backgroundimage
		self.show()

	def addPhotonNoise(self):
		'''Add photon noise into an image.'''

		self.noiseimage = self.zeros()
		self.speak("Adding photon noise [sqrt(photons from stars and various backgrounds)].")
		noise_variance = self.image
		ok = noise_variance > 0

		noise = np.zeros_like(self.image)
		noise[ok] = np.sqrt(noise_variance[ok])

		assert(np.isfinite(noise).all())
		self.image += noise*np.random.randn(self.xsize, self.ysize)
		self.noiseimage = noise

		self.note = 'photonnoise'
		noisefilename = self.directory + self.note + '.fits'
		if not os.path.exists(noisefilename):
			self.writeToFITS(noise,noisefilename)
		self.addInputLabels()
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
		try:
			self.noiseimage = np.sqrt(self.noiseimage**2 + noise_variance)
		except:
			self.noiseimage = np.sqrt(noise_variance)

		# update image header
		self.addInputLabels()
		self.header['IREADNOI'] = ('True', 'read noise')
		self.show()


	def addSmear(self):
		'''Smear the image along the readout direction.'''
		self.speak("Adding readout smear.")
		self.speak("    assuming {0} second readout times on {1} second exposures.".format(self.camera.readouttime,self.camera.singleread))


		untouched = self.image + 0.0
		mean =np.mean(self.image, 0).reshape(1,self.image.shape[0])*self.ones()
		self.image += mean*self.camera.readouttime/self.camera.singleread

		self.note = 'readoutsmear'
		smearfilename = self.directory + self.note + '.fits'
		if not os.path.exists(smearfilename):
			self.writeToFITS(self.image - untouched,smearfilename)

		# update header
		self.addInputLabels()
		self.header['ISMEAR'] = ('True', 'smearing during transer to frame store')
		self.show()


	def expose(self, plot=False, jitter=False, write=False, split=False, remake=False, smear=True, terse=False, cosmics='fancy', diffusion=False, correctcosmics=True, writenoiseless=True, skipcosmics=False, **kwargs):
		'''Expose an image on this CCD.'''
		self.plot = plot
		self.terse = terse

		# create a blank image
		self.image = self.zeros()

		self.populateHeader()
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

		if writenoiseless:
			# make filename for this image
			self.note = 'noiseless_'+self.fileidentifier
			noiselessfilename = self.directory + self.note + '.fits' + zipsuffix

			# write the image to FITS
			self.speak('saving noiseless TESS image')
			self.writeToFITS(self.image, noiselessfilename, savetype=np.int32)


		# add the photon noise from stars, galaxies, and backgrounds
		self.addPhotonNoise()

		if skipcosmics == False:
			# add cosmic rays to the image (after noise, because the *sub-Poisson* noise is already modeled with the Fano factor)
			cosmics = self.addCosmics(write=write, version=cosmics, diffusion=diffusion, correctcosmics=correctcosmics)

		# add smear from the finite frame transfer time
		if smear:
			self.addSmear()

		# create saturation bleed trails
		self.bleedSaturated()

		# add read noise, constant across detector
		self.addReadNoise()


		# finally, update the header for the image
		#self.populateHeader()

		if self.camera.stamps is not None:
			self.stampify()

		# write the image
		if write:
			self.writeFinal()

		self.report("created image #{counter:07d} of {pos_string} with {cadence:.0f}s cadence".format(counter=self.camera.counter, pos_string=self.pos_string, cadence=self.camera.cadence))

		# advance the camera's counter (and therefore timestep) if this is the last of the CCDs
		if self == self.camera.ccds[-1]:
			self.camera.advanceCounter()

		if write==False:
			return self.image, cosmics, stars


class Aberrator(Talker):
	'''object to keep track of how to apply velocity abberation; must be reset for each CCD'''
	def __init__(self, cartographer):
		Talker.__init__(self)

		# create a grid of stars spanning the CCD
		ccd = cartographer.ccd
		self.ccd = ccd
		ngrid = 20
		xgrid,ygrid = np.meshgrid(np.linspace(ccd.xmin, ccd.xmax, ngrid),
							 np.linspace(ccd.ymin, ccd.ymax, ngrid))
		x, y = xgrid.flatten(), ygrid.flatten()


		# estimate dx/delon and dy/delon


		self.strings, self.derivatives, self.coefs, self.inputs = {}, {}, {}, {}
		self.raw = {}
		template = '{0:+.10f}*ccdx{1:+.10f}*ccdy{2:+.10f}'
		A = np.vstack([np.ones(len(x)), x,y,x**2,y**2, x*y]).T

		delta=1.0/60.0/60.0 # step, in degrees, for calculating numerical derivative
		notnudged = cartographer.point(x,y,type='ccdxy')
		elon,elat = notnudged.celestial.tuple

		nudgedinlon = cartographer.point(elon+delta,elat,type='celestial')
		dxdelon = (nudgedinlon.ccdxy.x - x)/delta
		self.coefs['x'] = np.linalg.lstsq(A,dxdelon)[0]
		def model_dxdelon(x,y):
			return self.coefs['x'][0] +  self.coefs['x'][1]*x + self.coefs['x'][2]*y  + self.coefs['x'][3]*x**2 + self.coefs['x'][4]*y**2  + self.coefs['x'][5]*x*y # pixels/degree
		self.derivatives['x'] = model_dxdelon
		self.strings['x'] = template.format(*self.coefs['x'])
		self.raw['x'] = dxdelon
		self.inputs['x'] = x

		dydelon = (nudgedinlon.ccdxy.y - y)/delta
		self.coefs['y'] = np.linalg.lstsq(A,dydelon)[0]
		def model_dydelon(x,y):
			return self.coefs['y'][0] +  self.coefs['y'][1]*x + self.coefs['y'][2]*y  + self.coefs['y'][3]*x**2 + self.coefs['y'][4]*y**2  + self.coefs['y'][5]*x*y # pixels/degree

		self.derivatives['y'] = model_dydelon
		self.strings['y'] = template.format(*self.coefs['y'])
		self.raw['y'] = dydelon
		self.inputs['y'] = y

		#plt.ion()
		'''
		plt.figure()
		plt.scatter(x, 3600.0/dxdelon, c=y)
		plt.ylabel('arcsec of longitude/xpix')
		plt.figure()

		plt.scatter(x, 3600.0/dxdelat, c=y)
		plt.ylabel('arcsec of longitude/xpix')

		assert(False)
		'''
		plt.figure(figsize=(20,10))
		gs = plt.matplotlib.gridspec.GridSpec(2,3, hspace=0.3, bottom=.2)
		for i,k in enumerate(['x','y']):
			plt.subplot(gs[i,0])
			plt.scatter(self.inputs[k], self.raw[k], c=y, edgecolor='none')
			plt.ylabel('dccd{0}/delon (pixels/degrees)'.format(k))
			if i == 1:
				plt.xlabel('directly\ncalculated')
			plt.subplot(gs[i,1])
			plt.scatter(self.inputs[k], self.derivatives[k](x,y), c=y, edgecolor='none')
			if i == 1:
				plt.xlabel('model')

			plt.subplot(gs[i,2])
			plt.scatter(self.inputs[k], self.raw[k] - self.derivatives[k](x,y), c=y, edgecolor='none')
			if i == 1:
				plt.xlabel('residuals')

		plt.savefig(self.ccd.directory + 'aberrationgeometry.pdf')

	def plotPossibilities(self, n=100):
		x, y = np.random.uniform(0, self.ccd.xsize, n), np.random.uniform(0, self.ccd.ysize, n)
		stars = self.ccd.camera.cartographer.point(x,y,type='ccdxy')

		dx,dy,delon = [],[], []
		fcdx, fcdy = [], []
		bjds = np.linspace(0,365,1000)+self.ccd.bjd0
		for bjd in bjds:
			nudges = self.ccd.aberrations(stars,bjd)
			dx.append(nudges[0])
			dy.append(nudges[1])
			fcdx.append(self.ccd.fcdx)
			fcdy.append(self.ccd.fcdy)
			delon.append(self.ccd.delon)

		dx = np.array(dx)
		dy = np.array(dy)
		fcdx = np.array(fcdx)
		fcdy = np.array(fcdy)
		plt.figure(figsize=(8,8))
		gs = plt.matplotlib.gridspec.GridSpec(2,2, left=0.15, wspace=0.3)
		bjds -= min(bjds)
		# top row is uncorrected
		ax = plt.subplot(gs[0,0])
		plt.axvline(27.4, color='gray', alpha=1)
		ax.plot(bjds,dx+fcdx.reshape(len(bjds),1), alpha=0.3)
		plt.ylabel('velocity\naberration (pixels)')
		plt.title('x')
		ax = plt.subplot(gs[0,1],sharey=ax,sharex=ax)
		ax.plot(bjds,dy+fcdy.reshape(len(bjds),1), alpha=0.3)
		plt.axvline(27.4, color='gray', alpha=1)
		plt.xlim(-1+ min(bjds), max(bjds)+1)
		plt.title('y')

		ax = plt.subplot(gs[1,0])
		ax.plot(bjds,dx, alpha=0.3)
		plt.axvline(27.4, color='gray', alpha=1)
		plt.xlabel('Time (days)')
		plt.ylabel('differential velocity\naberration (pixels)')

		ax = plt.subplot(gs[1,1],sharey=ax,sharex=ax)
		ax.plot(bjds,dy, alpha=0.3)
		plt.axvline(27.4, color='gray', alpha=1)
		plt.xlim(-1 + min(bjds), max(bjds)+1)
		plt.xlabel('Time (days)')
		plt.savefig(self.ccd.directory+'aberrationoveroneyear.pdf')

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
