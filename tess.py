'''Tools to generate simulated TESS images.'''
# import necessary packages
import settings, catalogs
import numpy as np, matplotlib.pyplot as plt
import scipy.ndimage, scipy.signal, scipy.interpolate
import astropy.io, astropy.units, astropy.coordinates, astropy.wcs
import zachopy.utils, crossfield
import os, copy
import subprocess
from cosmical._cosmical import cosmical

# setup basic output options for this Python session
np.set_printoptions(threshold = 1e6, linewidth = 300)
plt.ion()





# define everything related to PSFs
class PSF:

	# initialize the PSF class
	def __init__(self, camera=None):

		# how much to bin the original (very high resolution) PSFs before integrating over the pixels
		self.initial_binning = 4


		self.subpixsize = 0.25/15.0*self.initial_binning
		self.pixsize = 1.0
		self.dtemperature = 2000
		self.noffset = 21
		self.drotation = 10
		self.nradii = 10
		self.nsubpix = 960/self.initial_binning

		self.dxax = np.arange(self.nsubpix)*self.subpixsize
		self.dxax -= np.mean(self.dxax)
		self.dyax = self.dxax
		self.dx, self.dy = np.meshgrid(self.dxax, self.dyax)
		self.gridpixels = self.dx.shape[0]*self.subpixsize/self.pixsize/2
		self.subgridresolution=self.subpixsize/self.pixsize
		self.subgrid = np.arange(-self.gridpixels,self.gridpixels , self.subgridresolution)
		self.xsubgrid, self.ysubgrid = np.meshgrid(self.subgrid, self.subgrid)
		self.grid = np.arange(-self.gridpixels,self.gridpixels,1,dtype=np.int)
		self.xgrid, self.ygrid = np.meshgrid(self.grid, self.grid)
		self.intrapixel = prnu(self.dx, self.dy)

		self.setCamera(camera)
		#self.populateBinned()
		self.populateHeader()

	def populateHeader(self):
		self.header = astropy.io.fits.Header()
		self.header[''] = ('', '')
		self.header['PRF'] = ''
		self.header['PRFNOTE'] = ('','Pixel Response Function library parameters')
		self.header['PLIBDTEM'] = (self.dtemperature, '[K] d(effective temperature)')
		self.header['PLIBNOFF'] = (self.noffset, '# of subpixel shifts, in both x and y')
		self.header['PLIBDROT'] = (self.drotation, '[deg] d(rotation around center)')
		self.header['PLIBNRAD'] = (self.nradii, '# of radial distances from field center')
		self.header['PNSUBPIX'] = (self.nsubpix, '[pix] # of initial subpixels')
		#self.header['PSUBSIZE'] = (self.subpixsize, '[pix] initial subpixel grid-size')
		#self.header['PPIXSIZE'] = (self.pixsize, '[pix] pixel size')


	# load the raw high resolution PSFs from Deb (no jitter yet)
	def populateMonochromaticLibrary(self, plot=True, inspect=False):
		print "Populating the raw, monochromatic, PSF library."
		try:
			(self.monochromaticlibrary, self.positions, self.wavelengths, self.dx, self.dy, self.dxax, self.dyax)  = np.load(settings.prefix + 'intermediates/psf_library.npy')
		except:
			# load the high resolution PSF's from a FITs file
			data = astropy.io.fits.open(settings.prefix + 'inputs/psfs_3p33.fits')[0].data
			# trim the high resolution PSF, so that an even number of whole pixels is included
			assert(data.shape[2] == data.shape[3])

			binned = zachopy.utils.rebin(data,data.shape[0], data.shape[1], data.shape[2]/self.initial_binning, data.shape[3]/self.initial_binning)
			ntrim = (binned.shape[2] - self.nsubpix)/2
			data = binned[:,:,ntrim:-ntrim,ntrim:-ntrim]
			npositions, nwavelengths, nxs, nys = data.shape
			self.monochromaticlibrary = {}

			positions = np.array([0,0.5,1.0,np.sqrt(2)])*2048
			wavelengths = 0.675 + np.arange(nwavelengths)*0.05
			self.positions = positions
			self.wavelengths = wavelengths
			assert(nxs == nys)

			extent=[self.dxax.min(), self.dxax.max(), self.dyax.min(), self.dyax.max()]
			if plot:
				fi, ax_psf_fits  = plt.subplots(data.shape[0],data.shape[1], figsize=(25.95, 9.8), sharex=True, sharey=True)

			for i in range(data.shape[0]):

				self.monochromaticlibrary[positions[i]] = {}
				for j in range(data.shape[1]):

					#interpolator = scipy.interpolate.interp2d(rho.flatten()[::10], phi.flatten()[::10], data[i,j,:,:].flatten()[::10], kind='linear', fill_value=0.0)
					#assert(False)
					self.monochromaticlibrary[positions[i]][wavelengths[j]] = data[i,j,:,:]
					if plot:
						if inspect:
							try:
								for a in ax_psf_zoom:
									a.cla()
							except:
								fig = plt.figure(figsize=(10,10))
								plt.subplots_adjust(hspace=0, wspace=0)
							ax_map = fig.add_subplot(2,2,3)
							ax_vert = fig.add_subplot(2,2,4, sharey=ax_map)
							ax_hori = fig.add_subplot(2,2,1, sharex=ax_map)
							ax_psf_zoom = [ax_map, ax_vert, ax_hori]
							ax_hori.plot(self.dxax, np.sum(self.monochromaticlibrary[positions[i]][wavelengths[j]], 0))
							ax_vert.plot(np.sum(self.monochromaticlibrary[positions[i]][wavelengths[j]], 1), self.dyax)
							ax_map.imshow(np.log(self.monochromaticlibrary[positions[i]][wavelengths[j]]),extent=extent, cmap='gray_r', interpolation='nearest', vmin=0  )
							plt.draw()
						ax_psf_fits[i,j].imshow(np.log(self.monochromaticlibrary[positions[i]][wavelengths[j]]), extent=extent, cmap='gray_r', interpolation='nearest', vmin=0 )
						if i == data.shape[0]-1:
							ax_psf_fits[i,j].set_xlabel("{0} micron".format(wavelengths[j]))

						if j == 0:
							ax_psf_fits[i,j].set_ylabel("{0:4.1f} pix.".format(positions[i]))
						plt.subplots_adjust(hspace=0, wspace=0)
			if plot:
				plt.draw()
				plt.savefig(settings.prefix + 'plots/psf_library.pdf')


			np.save(settings.prefix + 'intermediates/psf_library.npy',(self.monochromaticlibrary, self.positions, self.wavelengths, self.dx, self.dy, self.dxax, self.dyax))

	# integrate the monochromatic PSFs into a temperature library (no jitter yet)
	def populateTemperatureLibrary(self, plot=True):
		print "Populating the wavelength-integrated (unjittered) PSF library."
		temperature_filename = settings.prefix + 'intermediates/temperature_psfs.npy'.format(self.camera.cadence)
		try:
			self.temperaturelibrary, self.temperatures, self.positions = np.load(temperature_filename)
		except:
			try:
				self.monochromaticlibrary
			except:
				self.populateMonochromaticLibrary()
			weighting = astropy.io.fits.open(settings.prefix + 'inputs/ph_filt.fits')[0].data
			self.temperatures = np.arange(3000, 12001, self.dtemperature)
			for position in self.positions:
				try:
					self.temperaturelibrary
				except:
					self.temperaturelibrary = {}

				# loop over stellar temperatures
				for temperature in self.temperatures:
					try:
						self.temperaturelibrary[position]
					except:
						self.temperaturelibrary[position] = {}
					print "   modeling a PSF at {0}".format( position)
					print "      with a temperature of {0}K".format(temperature)
					psf = np.zeros_like(self.dx)
					# loop over wavelengths
					for w in np.arange(weighting.shape[1]):
						# loop over powers in the polynomial for weight vs. teff
						this_weight = 0.0
						for p in np.arange(weighting.shape[0]):
							this_weight += weighting[p,w]*(4000.0/temperature)**p
						psf += this_weight*self.monochromaticlibrary[position][self.wavelengths[w]]
						print "         {0} X [{1}]".format(this_weight, self.wavelengths[w])

					filename = "integrated_psf_for_{0}Kat{1}".format(temperature, position*100).replace(' ', '').replace('.','p') + '.pdf'
					self.temperaturelibrary[position][temperature] = psf/np.sum(psf)
					if plot:
						self.plot(psf, title="{0}K star at position (0,{1})".format(temperature, position),output=settings.prefix + 'plots/' + filename)
					print ''
			np.save(temperature_filename, (self.temperaturelibrary, self.temperatures, self.positions))

	# convolve the high resolution PSFs with the jitter expected for the camera's cadence
	def populateJitteredLibrary(self):
		print "Populating the high-resolution PSF library for an expected jitter for a {0:04.0f}s cadence.".format(self.camera.cadence)
		jittered_filename = settings.prefix + 'intermediates/jittered_psfs_{0:04.0f}s.npy'.format(self.camera.cadence)
		try:
			self.jitteredlibrary, self.temperatures, self.positions, self.jitteredlibrarytime = np.load(jittered_filename)
		except:
			try:
				self.temperaturelibrary
			except:
				self.populateTemperatureLibrary()

			self.jitteredlibrary = copy.copy(self.temperaturelibrary)
			for position in self.positions:
				for temperature in self.temperatures:
					self.jitteredlibrary[position][temperature] = scipy.signal.convolve2d(self.temperaturelibrary[position][temperature], self.camera.jittermap[0]/np.sum(self.camera.jittermap[0]), 'same', 'fill', 0)
			self.jitteredlibrarytime = self.camera.cadence
			np.save(jittered_filename, (self.jitteredlibrary, self.temperatures, self.positions, self.jitteredlibrarytime))
		assert(self.jitteredlibrarytime == self.camera.cadence)

	# populate a library of binned PRFS, using the jittered high-resolution library
	def populateBinned(self):
		binned_filename = settings.prefix + 'intermediates/binned_psf_library_{0:04.0f}s.npy'.format(self.camera.cadence)
		try:
			print "Trying to load ", binned_filename
			(self.binned, self.radii, self.temperatures, self.rotations, self.xoffsets, self.yoffsets) = np.load(binned_filename)
		except:
			print "Didn't find {0}, so creating a new binned library."
			self.populateMonochromaticLibrary()
			self.populateTemperatureLibrary()
			self.rotations = np.arange(0, 360, self.drotation)
			self.xoffsets = np.linspace(0,1,self.noffset)
			self.yoffsets = np.linspace(0,1,self.noffset)
			self.radii = np.linspace(0, np.max(self.positions),self.nradii)
			try:
				self.binned
			except:
				self.binned = {}
			for position in self.radii:
				try:
					self.binned[position]
				except:
					self.binned[position] = {}

				for temperature in self.temperatures:
					try:
						self.binned[position][temperature]
					except:
						self.binned[position][temperature] = {}

					for rotation in self.rotations:
						try:
							self.binned[position][temperature][rotation]
						except:
							self.binned[position][temperature][rotation] = {}

						for xoffset in self.xoffsets:
							try:
								self.binned[position][temperature][rotation][xoffset]
							except:
								self.binned[position][temperature][rotation][xoffset] = {}

							for yoffset in self.yoffsets:
								radius = position
								# make sure the sign is right!
								x = np.round(radius*np.cos(-rotation*np.pi/180)) + xoffset
								y = np.round(radius*np.sin(-rotation*np.pi/180)) + yoffset
								self.binned[position][temperature][rotation][xoffset][yoffset] = self.binnedhighrespsf( x, y, temperature)
								#self.plot(subgrid_psf*intrapixel,title='({0},{1},{2},{3},{4})'.format(position, temperature, rotation, xoffset, yoffset))
						print '({0},{1},{2},{3},{4})'.format(position, temperature, rotation, xoffset, yoffset)
				np.save(binned_filename, (self.binned, self.radii, self.temperatures, self.rotations, self.xoffsets, self.yoffsets))

	def binnedhighrespsf(self, x, y, temperature):
		subgrid_psf = self.highrespsf(x, y, temperature)
		assert(subgrid_psf.shape == self.intrapixel.shape)
		thisbinned = zachopy.utils.rebin_total(subgrid_psf*self.intrapixel,2*self.gridpixels,2*self.gridpixels)
		return thisbinned

	def binnedpsf(self, x, y, temperature, verbose=False, interpolation='linear'):

		try:
			self.binned
		except:
			self.populateBinned()
		r = np.sqrt(x**2 + y**2)
		position = zachopy.utils.find_nearest(self.radii, r, verbose=verbose)
		temperature = zachopy.utils.find_nearest(self.temperatures, temperature, verbose=verbose)
		rot = self.rotation(x,y)*180/np.pi
		if verbose:
			print ""
			print "{0},{1},{2},{3}".format(x, y, temperature, rot)
		while rot > 360 or rot < 0:
			if rot < 0:
				rot += 360
			if rot > 360:
				rot -= 360

		rotation = zachopy.utils.find_nearest(self.rotations, rot , verbose=verbose)
		if interpolation == 'nearest':
			xoffset = zachopy.utils.find_nearest(self.xoffsets, x % 1, verbose=verbose)
			yoffset = zachopy.utils.find_nearest(self.yoffsets, y % 1, verbose=verbose)
			return self.binned[position][temperature][rotation][xoffset][yoffset]
		if interpolation == 'linear':
			xfrac, yfrac = x % 1, y % 1
			xbelow, xabove = zachopy.utils.find_two_nearest(self.xoffsets, xfrac, verbose=verbose)
			xspan = xabove - xbelow
			xbelow_weight, xabove_weight = (xabove - xfrac)/xspan, (xfrac - xbelow)/xspan
			assert(xbelow_weight > 0)
			if verbose:
				print "  weights are ", xbelow_weight, xabove_weight
			ybelow, yabove = zachopy.utils.find_two_nearest(self.yoffsets, yfrac, verbose=verbose)
			yspan = yabove - ybelow
			ybelow_weight, yabove_weight =  (yabove - yfrac)/yspan, (yfrac - ybelow)/yspan
			if verbose:
				print "  weights are ", ybelow_weight, yabove_weight
			assert(xbelow_weight + xabove_weight ==1)
			assert(ybelow_weight + yabove_weight ==1)

			return 	xbelow_weight*ybelow_weight*self.binned[position][temperature][rotation][xbelow][ybelow] + \
					xabove_weight*ybelow_weight*self.binned[position][temperature][rotation][xabove][ybelow] + \
					xabove_weight*yabove_weight*self.binned[position][temperature][rotation][xabove][yabove] + \
					xbelow_weight*yabove_weight*self.binned[position][temperature][rotation][xbelow][yabove]


	def setCamera(self, camera):
		self.camera = camera

	def plotTemperatures(self):
		try:
			self.temperaturelibrary
		except:
			self.populateTemperatureLibrary()

		positions = self.temperaturelibrary.keys()
		temperatures = self.temperaturelibrary[positions[0]].keys()
		positions.sort()
		temperatures.sort()
		fi, ax = plt.subplots(len(positions), len(temperatures), figsize=(len(temperatures)*4, len(positions)*4), sharex=True, sharey=True)
		extent=[self.dxax.min(), self.dxax.max(), self.dyax.min(), self.dyax.max()]
		fi.suptitle("PSFs convolved with 2D jitter in a {0}s exposure.".format(self.camera.cadence))
		for p in range(len(positions)):
			for t in range(len(temperatures)):
				ax[p,t].imshow(np.log(psf), extent=extent, cmap='gray_r', interpolation='nearest', vmin=0 )

				if p == len(positions)-1:
					ax[p,t].set_xlabel("{0}K".format(temperatures[t]))

				if t == 0:
					ax[p,t].set_ylabel("{0:4.1f} pix.".format(positions[p]))
		plt.show()
		plt.savefig(settings.prefix + 'plots/psfs_with_temperature_{0:04.0f}.pdf'.format(self.camera.cadence))


	def plot(self, psf, title=None, output=None):
		try:
			for a in self.ax_psf_zoom:
				a.cla()
		except:
			fig = plt.figure(figsize=(10,10))
			plt.subplots_adjust(hspace=0, wspace=0)
			ax_map = fig.add_subplot(2,2,3)
			ax_vert = fig.add_subplot(2,2,4, sharey=ax_map)
			ax_hori = fig.add_subplot(2,2,1, sharex=ax_map)
			self.ax_psf_zoom = (ax_map, ax_vert, ax_hori)

		(ax_map, ax_vert, ax_hori) = self.ax_psf_zoom
		ax_hori.plot(self.dxax, np.sum(psf, 0)/np.sum(psf,0).max(), color='black', linewidth=3)
		ax_vert.plot(np.sum(psf, 1)/np.sum(psf,1).max(),self.dyax, color='black', linewidth=3)
		ax_vert.semilogx()
		ax_hori.semilogy()
		ax_hori.set_ylim(1e-6,1e0)
		ax_vert.set_xlim(1e-6,1e0)

		ax_vert.tick_params(labelleft=False)
		ax_hori.tick_params(labelbottom=False)
		if title is not None:
			ax_hori.set_title(title)
		ax_map.imshow(np.log(psf), cmap='gray_r', extent=[self.dxax.min(), self.dxax.max(), self.dyax.min(), self.dyax.max()],vmin=0.0,interpolation='nearest')
		plt.draw()
		if output is not None:
			ax_map.figure.savefig(output)
			print "   saved figure to: ", output

	def rotation(self, xpos, ypos):
		return -np.arctan2(ypos, xpos)

	def highrespsf(self, xpos, ypos, temperature):

		try:
			self.jitteredlibrary
		except:
			self.populateJitteredLibrary()

		possible_positions = self.jitteredlibrary.keys()
		r = np.sqrt(xpos**2 + ypos**2)
		indices = np.arange(len(self.positions))
		i_fraction = np.interp(r, self.positions, indices)
		below = self.jitteredlibrary[self.positions[np.floor(i_fraction)]]
		above = self.jitteredlibrary[self.positions[np.ceil(i_fraction)]]
		frac = i_fraction % 1.0
		#print "   {0}".format( frac)
		rounded_temperature = zachopy.utils.find_nearest(self.temperatures, temperature)
		unrotated_psf = below[rounded_temperature]*(1.0 - frac) + above[rounded_temperature]*frac
		rotated_psf = scipy.ndimage.interpolation.rotate(unrotated_psf, self.rotation(xpos, ypos)*180/np.pi, reshape=False)
		nudged_psf = np.zeros_like(rotated_psf)
		xnudge = np.int((xpos % 1)/self.subpixsize)
		ynudge = np.int((ypos % 1)/self.subpixsize)
		#print xnudge, ynudge, ' nudges'
		nudged_psf[xnudge:, ynudge:]= rotated_psf[:rotated_psf.shape[0]-xnudge, :rotated_psf.shape[1]-ynudge]
		#print [xnudge, ynudge]
		#print [rotated_psf.shape[0]-xnudge, rotated_psf.shape[1]-ynudge]
		return nudged_psf

	def xy2rhophi(self, x, y, rotation=0.0):
		rho = np.sqrt(x**2 + y**2)
		phi = np.arctan2(y, x) + rotation
		return rho, phi

	def rhophi2xy(self, rho, phi, rotation=0.0):
		x = np.cos(phi - rotation)*rho
		y = np.sin(phi - rotation)*rho
		return x,y

# use matrices to handle the rotations due to jitter
def rx(theta):
	return np.mat([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
def ry(theta):
	return np.mat([[np.cos(theta), 0, np.sin(theta)],[0, 1, 0], [-np.sin(theta), 0, np.cos(theta)]])
def rz(theta):
	return np.mat([[np.cos(theta), -np.sin(theta), 0],[np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
def rotate(ra, dec, x, y):
	cart = np.transpose(np.mat(astropy.coordinates.spherical_to_cartesian(1,dec*np.pi/180, ra*np.pi/180)))
	#assert(x==0.0)
	rotated = rx(x*np.pi/180)*ry(y*np.pi/180)*cart
	#print cart
	#print ' --> '
	#print rotated
	sphe = np.array(astropy.coordinates.cartesian_to_spherical(rotated[0,0], rotated[1,0], rotated[2,0]))*180/np.pi
	print "      {0:f}, {1:f} --> {2:f},{3:f}".format(ra, dec, sphe[2], sphe[1])
	return sphe[2], sphe[1]

# define a camera class
class Camera:
	def __init__(self, stamp=None, cadence=1800, ra=270,dec=66.56070833333332):
		print "Initializing a new TESS camera object."
		self.pixelscale = 24.0/4096*60*60							# arcsec!!
		self.offset_from_ecliptic = 6.0									# degrees
		self.number_of_segments = 13
		self.entrance_pupil_diameter = 10.5								# cm
		self.effective_area = 67.5										# cm^2
		self.physicalpixelsize = 15.0/1e4								# cm
		self.physicalpixeldepth = 100.0/1e4								# cm
		self.read_noise = 10.0											# electrons per read
		self.saturation = 150000.0										# electrons per read
		self.shortcadence = 2.0*60.0									# seconds
		self.longcadence = cadence										# seconds
		self.cadence = self.longcadence
		self.singleread = 2.0											# seconds
		self.readouttime = 0.005										# seconds

		if stamp is None:
			self.physicalgap = 0.2											# cm
			self.gapinpixels = np.round(self.physicalgap/self.physicalpixelsize).astype(np.int)		# pixels
			self.npix = 4096 + self.gapinpixels
		else:
			self.physicalgap = 0.0
			self.gapinpixels = 0
			self.npix = stamp

		self.fov = self.npix*self.pixelscale/60.0/60.0					# degrees
		self.nudge = {'x':0.0, 'y':0.0, 'z':0.0}						# nudge relative to nominal spacecraft pointing (arcsec)
		self.counter = 0
		self.pixel_solid_area = (self.pixelscale)**2 		# arcsec^2
		self.psf = PSF(camera=self)
		self.setCadence(self.longcadence)
		self.psf.setCamera(self)
		self.stamp = stamp
		self.point(ra, dec)

	def populateHeader(self):
		self.header = astropy.io.fits.Header()
		self.header['CAMERA'] = ''
		self.header['CAMNOTE'] = ('', 'properties of the Camera')
		self.header['FOV'] = (self.fov, '[deg] field of view')
		self.header['SCALE'] = (self.pixelscale, '["/pix] pixel scale')
		self.header['DIAMETER'] = (self.entrance_pupil_diameter, '[cm] entrace pupil diameter')
		self.header['EFFAREA'] = (self.effective_area, '[cm^2] effective area (at ref. wavelength)')
		self.header['PIXSIZE'] = (self.physicalpixelsize, '[cm] physical pixel size')
		self.header['PIXDEPTH'] = (self.physicalpixeldepth, '[cm] physical pixel depth')
		self.header['PHYSIGAP'] = (self.physicalgap, '[cm] gap between CCDs')
		self.header['PIXELGAP'] = (self.gapinpixels, '[pix] gap size in pixels (rough)')

		self.header['WCS'] = ''
		self.header['WCSNOTE'] = ('', 'World Cooridinate System for this image')
		if self.stamp is not None:
			self.header['STAMP'] = (self.stamp, 'THIS IMAGE IS JUST {0}x{1} POSTAGE STAMP!'.format(self.stamp, self.stamp))

		self.header.extend(self.wcs.to_header())

		self.header['MOTION'] = ''
		self.header['MOTNOTE'] = ('', 'properties of the image motion applied')
		self.header['JITTERX'] = (0, '["] jitter-induced nudge')
		self.header['JITTERY'] = (0, '["] jitter-induced nudge')

	def setCadence(self, cadence=1800.0):
		self.cadence = cadence
		print "Setting cadence to {0} seconds = {1} reads.".format(self.cadence, self.cadence/self.singleread)
		self.loadJitterBall()

	def point(self, ra=None, dec=None):
		print "Pointing TESS at the sky."
		# point TESS, to be used both for initial pointing, and for jitters
		if ra is not None and dec is not None:
			self.ra = ra
			self.dec = dec

		try:
			self.ra
			self.dec
		except:
			print "Please point your telescope somewhere. No RA or DEC defined."


		# the number of dimensions
		w = astropy.wcs.WCS(naxis=2)

		# the pixel coordinates of the reference position (taken to be center of field)
		w.wcs.crpix = [self.npix/2.0,self.npix/2.0]

		# the pixel scale, in degrees, taken to be FOV/number of pixels
		w.wcs.cdelt = [-self.fov/self.npix,self.fov/self.npix]

		# the celestial coordinates at the reference position (input by user)
		nudged_ra, nudged_dec = rotate(self.ra, self.dec,  self.nudge['x']/60.0/60.0, self.nudge['y']/60.0/60.0)
		print "!!!!!!!!!!!!!!!"
		w.wcs.crval = [nudged_ra, nudged_dec]

		# the rotation of the field (currently just a pure rotation, no shear)
		#rot = self.nudge['z']/60.0/60.0*np.pi/180.0
		#w.wcs.pc = [[np.cos(rot), -np.sin(rot)],[np.sin(rot), np.cos(rot)]]

		# the coordinate system type - what should I use?
		w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

		# set this to be the WCS
		self.wcs = w
		self.populateHeader()

	def loadJitterBall(self):
		print "Loading jitterball for {0} second cadences.".format(self.cadence)
		try:
			# check to see if the jitterball is already loaded
			self.jitterball
			# make sure the we're using the right jitterball for this cadence
			assert(self.jittercadence == self.cadence)
		except:
			# if the jitterball isn't already loaded, then do so
			jitterfile = settings.prefix + 'intermediates/jitter_{0:04.0f}.npy'.format(self.cadence)
			try:
				# if a processed jitterfile already exists, then just load it up:
				self.jitterball, self.jittermap = np.load(jitterfile)
				self.jittercadence = self.cadence
			except:
				print "   taking the raw jitter input file AttErrTimeArcsec_80k.dat"
				# otherwise, create a binned jitterfile
				# load simulated jitterball that Roland got from Orbital
				try:
					data = np.load(settings.prefix + 'intermediates/raw_jitter_file.npy')
				except:
					data = astropy.io.ascii.read(settings.prefix + "inputs/AttErrTimeArcsec_80k.dat", names=['t','x', 'y', 'z'])
					np.save(settings.prefix + 'intermediates/raw_jitter_file.npy', data)

				# subtract means
				data['x'] -= np.mean(data['x'])
				data['y'] -= np.mean(data['y'])
				data['z'] -= np.mean(data['z'])
				# scale jitterball to requirements (should be inflation by ~1.5)
				required_jitter_rms = 2.0/3
				original_fifthsec = np.sqrt(np.mean(data['x']**2 + data['y']**2))
				data['x'] *= required_jitter_rms/original_fifthsec
				data['y'] *= required_jitter_rms/original_fifthsec
				data['z'] *= required_jitter_rms/original_fifthsec

				# smooth them to the required cadence
				print "   smoothing the jitter to {0}s cadence.".format(self.cadence)
				spacing = data['t'][1] - data['t'][0]
				n = np.long(self.cadence/spacing)
				filter = np.ones(n)/n
				smoothed_t = np.convolve(data['t'], filter, mode='valid')
				smoothed_x = np.convolve(data['x'], filter, mode='valid')
				smoothed_y = np.convolve(data['y'], filter, mode='valid')
				smoothed_z = np.convolve(data['z'], filter, mode='valid')


				t = np.convolve(data['t'], filter, mode='valid')[::n]
				x = np.convolve(data['x'], filter, mode='valid')[::n]
				y = np.convolve(data['y'], filter, mode='valid')[::n]
				z = np.convolve(data['z'], filter, mode='valid')[::n]

				# plot each dimension separately
				print "   plotting the binned jitter timeseries to " + settings.prefix + 'plots/jitter_timeseries_{0:04.0f}.pdf'.format(self.cadence)
				fi, ax = plt.subplots(3,1, sharey=True, sharex=True)
				plt.ion()
				ax[0].plot(data['t'], data['x'], alpha=0.5, color='black')
				ax[0].plot(t, x, linewidth=2, alpha=0.5, marker='o', color='red')
				ax[1].plot(data['t'], data['y'], alpha=0.5, color='black')
				ax[1].plot(t, y, linewidth=2, alpha=0.5, marker='o', color='red')
				ax[2].plot(data['t'], data['z'], alpha=0.5, color='black')
				ax[2].plot(t, z, linewidth=2, alpha=0.5, marker='o', color='red')
				ax[0].set_xlim(0,self.cadence*10)
				ax[0].set_title('Expected TESS Pointing Jitter for {0}s Cadence'.format(self.cadence))
				ax[0].set_ylabel('x (")')
				ax[1].set_ylabel('y (")')
				ax[2].set_ylabel('z (")')
				ax[2].set_xlabel('Time (seconds)')
				plt.show()
				fi.savefig(settings.prefix + 'plots/jitter_timeseries_{0:04.0f}.pdf'.format(self.cadence))

				narcsec = 3
				bins = np.ceil(narcsec/self.pixelscale/self.psf.subpixsize).astype(np.int)*2 +1	# bins are in units of subpixels
				range = [[-(bins-1)/2*self.psf.subpixsize,(bins-1)/2*self.psf.subpixsize],[-(bins-1)/2*self.psf.subpixsize,(bins-1)/2*self.psf.subpixsize]] # range is in units of pixels

				x_interpolator = scipy.interpolate.interp1d(smoothed_t, smoothed_x,'nearest',fill_value=0,bounds_error=False)
				y_interpolator = scipy.interpolate.interp1d(smoothed_t, smoothed_y,'nearest',fill_value=0,bounds_error=False)
				jittermap_to_plot = np.histogram2d((data['x'] - x_interpolator(data['t']))/self.pixelscale, (data['y']  - y_interpolator(data['t']))/self.pixelscale, bins=50, range=range,normed=True)

				print "   plotting the binned jitter map to " + settings.prefix + 'plots/jitter_map_{0:04.0f}.pdf'.format(self.cadence)

				# plot map
				'''fi, ax = plt.subplots(1,1, sharey=True, sharex=True)
				plt.ion()
				ax.imshow(self.jittermap[0] , cmap='gray_r', extent=[-3,3,-3,3],interpolation='nearest')

				ax.set_xlim(-3,3)
				ax.set_ylim(-3,3)
				ax.set_title('TESS Pointing Jitter over {0}s'.format(self.cadence))
				ax.set_ylabel('x (")')
				ax.set_ylabel('y (")')
				plt.show()'''

				self.plothist2d(jittermap_to_plot,title='TESS Pointing Jitter over {0}s'.format(self.cadence),xtitle='Pixels', ytitle='Pixels')
				plt.savefig(settings.prefix + 'plots/jitter_map_{0:04.0f}.pdf'.format(self.cadence))




				print range
				print bins
				self.jittermap = np.histogram2d((data['x'] - x_interpolator(data['t']))/self.pixelscale, (data['y']  - y_interpolator(data['t']))/self.pixelscale, bins=bins, range=range,normed=True)
				self.plothist2d(self.jittermap ,title='TESS Pointing Jitter over {0}s'.format(self.cadence),xtitle='Pixels', ytitle='Pixels')
				plt.savefig(settings.prefix + 'plots/jitter_map_adopted_{0:04.0f}.pdf'.format(self.cadence))

				self.jitterball = (x,y,z)
				self.jittercadence = self.cadence
				np.save(jitterfile,( self.jitterball, self.jittermap))

	def plothist2d(self, hist, title=None, log=False, xtitle=None, ytitle=None):
		map = hist[0]
		x = hist[1][1:] + (hist[1][0] - hist[1][1])/2.0
		y = hist[2][1:]+ (hist[2][0] - hist[2][1])/2.0
		fig = plt.figure(figsize=(10,10))
		plt.clf()
		plt.subplots_adjust(hspace=0, wspace=0)
		ax_map = fig.add_subplot(2,2,3)
		ax_vert = fig.add_subplot(2,2,4, sharey=ax_map)
		ax_hori = fig.add_subplot(2,2,1, sharex=ax_map)

		ax_hori.plot(x, np.sum(map, 0)/np.sum(map), marker='o', color='black', linewidth=3)
		ax_vert.plot(np.sum(map, 1)/np.sum(map), y, marker='o', color='black', linewidth=3)
		if log:
			ax_vert.semilogx()
			ax_hori.semilogy()
		if log:
			bottom = np.min(map[map > 0])/np.maximum(np.sum(map,0).max(),np.sum(map,1).max())
		else:
			bottom = 0
		top = 1
		ax_hori.set_ylim(bottom,top)
		ax_vert.set_xlim(bottom,top)

		ax_vert.tick_params(labelleft=False)
		ax_hori.tick_params(labelbottom=False)
		if title is not None:
			ax_hori.set_title(title)
		if xtitle is not None:
			ax_map.set_xlabel(xtitle)
		if ytitle is not None:
			ax_map.set_ylabel(ytitle)
		if log:
			ax_map.imshow(np.log(map), cmap='gray_r', extent=[x.min(), x.max(), y.min(), y.max()],interpolation='nearest')
		else:
			ax_map.imshow(map, cmap='gray_r', extent=[x.min(), x.max(), y.min(), y.max()],interpolation='nearest')
		plt.draw()

	def jitter(self, header=None):
		try:
			self.jitterball
		except:
			self.loadJitterBall()
		print "Jittering TESS a little bit."
		self.nudge['x'] = self.counter*0.0#self.jitterball[0][self.counter]#*0.002
		self.nudge['y'] = -self.counter*3.0#self.jitterball[1][self.counter]#*0.002
		self.nudge['z'] = 0.0#self.jitterball[2][self.counter]#*0.002

		try:
			header['MOTION'] = ''
			header['MOTNOTE'] = ('', 'properties of the image motion applied')
			header['JITTERX'] = (self.nudge['x'], '["] jitter-induced nudge')
			header['JITTERY'] = (self.nudge['y'], '["] jitter-induced nudge')
		except:
			pass

		self.point()

	def pos_string(self):
		coords = astropy.coordinates.ICRS(ra=self.ra*astropy.units.degree, dec=self.dec*astropy.units.degree)
		return "{0:02}h{1:02}m{2:02}s{3:+03}d{4:02}m{5:02}s".format(np.int(coords.ra.hms[0]),np.int(coords.ra.hms[1]),np.int(coords.ra.hms[2].round()), np.int(coords.dec.dms[0]),np.int(np.abs(coords.dec.dms[1])),np.int(np.abs(coords.dec.dms[2].round())))

	def project(self, plot=False, write=False):
		print "Creating a starmap for this TESS field."
		ras, decs, rmag, jmag, imag, temperatures = catalogs.stars(ra=self.ra, dec=self.dec, radius=self.fov/np.sqrt(2)*1.01, catalog='UCAC4')
		deltamag = np.max(imag) - imag
		size = deltamag**2/16.0

		#try:
		#	self.ax_radec.cla()
		#except:
		#	fi, self.ax_radec = plt.subplots(1,1)
		#self.ax_radec.scatter(ras, decs, s=size, marker='.', color='black', alpha=0.3)

		x, y = self.wcs.wcs_world2pix(ras, decs, 1)
		starsfilename = settings.prefix + 'plots/' +  "starfield_{ra}_{dec}".format(ra=self.ra, dec=self.dec).replace(' ', '') + '.pdf'
		if not os.path.exists(starsfilename):
			if plot:
				try:
					self.ax_projected.cla()
				except:
					fi, self.ax_projected = plt.subplots(1,1)
				self.ax_projected.scatter(x, y, s=size, marker='o', color='black', alpha=0.3, edgecolors='none')
				self.ax_projected.plot([0, 4096, 4096, 0, 0], [0, 0, 4096, 4096, 0], color='red', alpha=0.5, linewidth=3)
				self.ax_projected.set_aspect(1)
				self.ax_projected.set_xlim(-128, 4096+128)
				self.ax_projected.set_ylim(-128, 4096+128)
				self.ax_projected.set_xlabel('Pixel Coordinates')
				self.ax_projected.set_ylabel('Pixel Coordinates')
				self.ax_projected.figure.savefig(starsfilename)

		if write:
			outfile = settings.prefix + 'intermediates/stars_{pos}_{fov:2.0f}deg.txt'.format(pos=self.pos_string(), fov=self.fov)
			np.savetxt(outfile, np.c_[ras, decs, x, y, imag], fmt=['%.6f', '%.6f', '%.3f', '%.3f', '%.3f'])
			print " saved stars to {0}".format(outfile)


		self.starx = x
		self.stary = y
		self.starmag = np.array(imag)
		self.startemp = np.array(temperatures)

	def photons(self, mag):
		# this doesn't work for M dwarfs, need to include multiband information
		return self.effective_area*10**(-0.4*mag)*1.7e6


class Image:
	def __init__(self, camera):
		self.camera = camera
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

		xsubpixel = self.camera.psf.xsubgrid + np.long(x)
		ysubpixel = self.camera.psf.ysubgrid + np.long(y)

		xrel, yrel = self.relativeposition(x,y)
		if verbose and self.starcounter % 10000 == 0:
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
			self.ax_prf.set_xlim(-self.psf.gridpixels, self.psf.gridpixels)
			self.ax_prf.set_ylim(-self.psf.gridpixels, self.psf.gridpixels)
			self.ax_prf.set_aspect(1)
			self.ax_prf.figure.savefig(settings.prefix + 'plots/prnu_demo.pdf')
			plt.draw()'''


		xindex =  -self.camera.psf.xgrid + np.long(x) - self.xmin
		yindex =  -self.camera.psf.ygrid + np.long(y) - self.ymin

		ok = (xindex >= self.xmin) * (xindex < self.xsize) * (yindex >= self.ymin) * (yindex < self.ysize)
		self.starimage[xindex[ok], yindex[ok]] += binned[ok]

	def addStars(self, remake=False, jitter=False):
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
					ok =    (self.camera.starx + self.camera.psf.gridpixels >= self.xmin) * \
							(self.camera.starx - self.camera.psf.gridpixels <= self.xmax) * \
							(self.camera.stary + self.camera.psf.gridpixels >= self.ymin) * \
							(self.camera.stary - self.camera.psf.gridpixels <= self.ymax) * \
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

	def addCosmicsAl(self, split=False, gradient=False, write=False):


		print "Adding cosmic rays."
		rate = 5 # cts cm^-2 s^-1
		self.note = 'cosmics_{0:.1f}persecond_{1}seconds_{2:06.0f}'.format(rate,self.camera.cadence, self.camera.counter).replace('.', 'p')
		cosmicsfilename = self.directory + self.note + '.fits'


		if gradient:
			smallexptime = self.camera.cadence
			bigexptime = self.camera.cadence*2
		else:
			smallexptime = self.camera.cadence*1.5
			bigexptime = smallexptime
		nexpected = (self.camera.physicalpixelsize*self.camera.npix)**2*(0.5*smallexptime + 0.5*bigexptime)*rate
		ndrawn = np.random.poisson(nexpected)
		image = np.transpose(cosmical(smallexptime, bigexptime, ndrawn, self.camera.npix, self.camera.npix))

		if write:
			self.writeToFITS(image, cosmicsfilename, split=split)

		self.image += image
		self.header['ICOSMICS'] = ('True', 'cosmic rays injected')
		return image

	def addCosmics(self, split=False):
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
		print "Adding backgrounds."

		# if the background image already exists, just load it
		try:
			self.backgroundimage = self.loadFromFITS(backgroundsfilename)
		# otherwise compute the backgrounds from scratch, and save them for next time
		except:
			# define a blank background image
			self.backgroundimage = self.zeros()
			# define coordinates (equatorial, Galactic, ecliptic) at every pixel in the image
			ra, dec = self.imageAD()
			gal = astropy.coordinates.ICRS(ra=ra, dec=dec, unit=(astropy.units.deg,astropy.units.deg)).galactic
			elon, elat = crossfield.euler(ra, dec, select=3)

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
		print "Removing the pixels in the CCD gaps."
		self.image[self.xgapstart:self.xgapend,:] = 0
		self.image[:, self.ygapstart:self.ygapend] = 0


	def addSmear(self):
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

	def expose(self, plot=False, jitter=False, write=False, split=False):
		self.plot = plot

		print "-------_-------_-------_-------_-------_-------_-------_-------"
		print "CREATING IMAGE #{0:06d}".format(self.camera.counter)
		print "----_----_----_----_----_----_----_----_----_----_----_----_---"

		# create a blank image
		self.image = self.zeros()

		# jitter the camera, or at least update the
		if jitter:
			self.camera.jitter(header=self.header)

		# add stars to the image
		self.addStars(jitter=jitter)

		# add galaxies to the image
		self.addGalaxies()


		# add background to the image
		self.addBackgrounds()

		if write==False:
			stars = self.image + 0.0

		# add the photon noise from stars, galaxies, and backgrounds
		self.addPhotonNoise()

		# add cosmic rays to the image (after noise, because the *sub-Poisson* noise is already modeled with the Fano factor)
		cosmics = self.addCosmicsAl(split=split, write=write)

		# add smear from the finite frame transfer time
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

def prnu(x,y,type='boxcar'):
	x_withinpixel = x % 1
	y_withinpixel = y % 1
	if type == 'boxcar':
		onedge = (np.abs(x_withinpixel) <= 0.1)+(np.abs(x_withinpixel) >= 0.9)+(np.abs(y_withinpixel) <= 0.1)+(np.abs(y_withinpixel) >= 0.9)
		response = np.ones_like(x)
		response[onedge] *= 0.95
	return response


plt.close('all')
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
