from imports import *
import settings, Catalogs
from PSF import PSF
from Cartographer import Cartographer
from CCD import CCD
from Jitter import Jitter

# define a camera class
class Camera(Talker):
    '''Keep track of one camera's entire field of view.'''
    def __init__(self, cadence=1800, ra=270,dec=66.56070833333332, testpattern=False, subarray=None, label='', number=1, magnitudes=[10]):
        '''Initialize camera, fill it with CCDs, and point it at the sky or at a testpattern.'''

        # decide whether or not this Camera is chatty
        Talker.__init__(self, mute=False, pithy=False)

        # keep track of which exposure is being simulated
        self.counter = 0

        # [self.speak will only report for this object if mute and pithy are turned off]
        self.speak("Turning on a new TESS camera object.")

        # keep track of any special features with this Camera
        self.subarray = subarray
        self.label = label
        self.testpattern = testpattern

        # decide whether to point the Camera either at real stars (from the sky) or at a test pattern (a grid of stars)
        if self.testpattern:
            # pretend a test pattern is pointed at (ra, dec) = (0.0, 0.0)
            self.ra = 0.0
            self.dec = 0.0
        else:
            # if real stars, use the input (ra, dec)
            self.ra = ra
            self.dec = dec

        # assign the cadence for this camera
        self.singleread = 2.0											# seconds
        self.readouttime = 0.005										# seconds
        self.setCadence(cadence)                                        # seconds

        # define scales for the Camera
        self.pixelscale = 21.1#24.0/4096*60*60							# arcsec!! (from Peter's paper)
        self.entrance_pupil_diameter = 10.5								# cm (from Peter's paper)
        self.effective_area = 69.1										# cm^2 (from Peter's paper)
        self.physicalpixelsize = 15.0/1e4								# cm
        self.physicalpixeldepth = 100.0/1e4								# cm
        self.read_noise = 10.0											# electrons per read
        self.saturation = 150000.0										# electrons per read

        # define the gaps between the CCDs (numbers still in flux, as of late 2014)
        self.physicalgap = 0.2											# cm
        self.gapinpixels = self.physicalgap/self.physicalpixelsize  	# pixels (not necessarily integer)

        # define the field of view of the Camera (one side of the CCD square)
        self.fov = (4096 + self.gapinpixels)*self.pixelscale/60.0/60.0					# degrees
        self.pixel_solid_area = (self.pixelscale)**2 		             # arcsec^2

        # turn on the necessary CCDs in this Camera
        if self.subarray is None:
            # if we're not dealing with a subarray, then turn on CCD's 1,2,3,4
            self.ccdnumbers = np.arange(4) + 1
        else:
            # if this is a subarray, then turn on one (imaginary) CCD and call it 0
            self.ccdnumbers = np.arange(1)
        self.ccds = [CCD(n,subarray=self.subarray,camera=self) for n in self.ccdnumbers]

        # assign a cartographer to this Camera and start it out on the first CCD
        self.cartographer = Cartographer(camera=self, ccd=self.ccds[0])



        # start the camera out unjittered from its nominal position
        self.nudge = {'x':0.0, 'y':0.0, 'z':0.0}						# nudge relative to nominal spacecraft pointing (arcsec)


        # point the Camera
        self.point(self.ra, self.dec)

        # load the PSF for this Camera
        self.psf = PSF(camera=self)

        # load the jitterball for this camera
        self.jitter = Jitter(camera=self)



    @property
    def fielddirectory(self):
        if self.label == '':
            d = settings.prefix + 'outputs/{pos}/'.format(pos=self.pos_string())
        else:
            d = settings.prefix + 'outputs/{pos}_{label}/'.format(pos=self.pos_string(), label=self.label)
        zachopy.utils.mkdir(d)
        return d

    @property
    def directory(self):
        d = self.fielddirectory + '{0:d}s/'.format(self.cadence)
        zachopy.utils.mkdir(d)
        return d

    def expose(self, **kwargs):
        '''Take an exposure on all the available CCD's.'''
        for c in self.ccds:
            c.expose(**kwargs)

    def populateHeader(self):
        '''Populate the header structure with information about the Camera, and its WCS.'''

        # create an empty header
        self.header = astropy.io.fits.Header()

        # fill it with some Camera details
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
        if self.subarray is not None:
            self.header['SUBARRAY'] = (self.subarray, 'THIS IMAGE IS JUST {0}x{1} POSTAGE STAMP!'.format(self.subarray, self.subarray))

        # fill it with the WCS information
        self.header['WCS'] = ''
        self.header['WCSNOTE'] = ('', 'World Cooridinate System for this image')
        self.header.extend(self.wcs.to_header())

        # note the jitter motions that have been applied to this image
        self.header['MOTION'] = ''
        self.header['MOTNOTE'] = ('', 'properties of the image motion applied')
        self.header['JITTERX'] = (0, '["] jitter-induced nudge')
        self.header['JITTERY'] = (0, '["] jitter-induced nudge')

    def setCadence(self, cadence=1800.0):
        '''Set the cadence of this Camera, in seconds.'''

        # set the cadence
        self.cadence = cadence
        self.speak("Setting cadence to {0} seconds = {1} reads.".format(self.cadence, self.cadence/self.singleread))

    def point(self, ra=None, dec=None):
        '''Point this Camera at the sky, by using the field-specified (ra,dec) and (if active) the jitter nudge for this exposure.'''

        # if (ra,dec) given, then point the whole Camera there
        if ra is not None and dec is not None:
            self.ra, self.dec = ra, dec

        # make sure an (ra,dec) are defined
        try:
            self.ra
            self.dec
            self.speak('Pointing the camera at (ra,dec) = {0:.6f},{1:.6f}'.format(self.ra, self.dec))
        except:
            self.report("Please point your telescope somewhere. No RA or DEC defined.")

        # create a blank WCS object, for converting between (ra,dec) and (x,y)
        self.wcs = astropy.wcs.WCS(naxis=2)

        # the focalxy coordinates of the reference position [taken to be center of field, which is (x,y) = (0.0, 0.0) in focalxy coordinates]
        self.wcs.wcs.crpix = [0.0, 0.0]

        # the pixel scale, in degrees
        self.wcs.wcs.cdelt = [-self.pixelscale/60.0/60.0,self.pixelscale/60.0/60.0]

        # the celestial coordinates at the reference position (input by user)
        nudged_ra, nudged_dec = zachopy.spherical.rotate(self.ra, self.dec,  self.nudge['x']/60.0/60.0, self.nudge['y']/60.0/60.0)
        self.wcs.wcs.crval = [nudged_ra, nudged_dec]

        # the rotation of the field (currently just a pure rotation, no shear)
        #rot = self.nudge['z']/60.0/60.0*np.pi/180.0
        #w.wcs.pc = [[np.cos(rot), -np.sin(rot)],[np.sin(rot), np.cos(rot)]]

        # the coordinate system type - what should I use?
        self.wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

        # set this to be the WCS
        #self.populateHeader()

    def pos_string(self):
        '''Return the position string for this field.'''
        if self.testpattern:
            try:
                return self.catalog.name
            except:
                return 'testpattern'
        else:
            coords = astropy.coordinates.ICRS(ra=self.ra*astropy.units.degree, dec=self.dec*astropy.units.degree)
            return "{0:02}h{1:02}m{2:02}s{3:+03}d{4:02}m{5:02}s".format(np.int(coords.ra.hms[0]),np.int(coords.ra.hms[1]),np.int(coords.ra.hms[2].round()), np.int(coords.dec.dms[0]),np.int(np.abs(coords.dec.dms[1])),np.int(np.abs(coords.dec.dms[2].round())))

    def advanceCounter(self):
        '''Take one step forward in time with this Camera.'''
        self.counter +=1

    def populateCatalog(self, **kwargs):
        '''Create a catalog of stars that are visible with this Camera.'''

        # figure out how wide we need to search for stars
        if self.subarray is None:
            size = self.fov
        else:
            size = self.subarray*self.pixelscale/3600.0

        # make a test pattern or a catalog of real stars
        if self.testpattern:
            # figure out how big a test pattern to create (in arcsec)

            self.catalog = Catalogs.TestPattern(size=size*3600.0, **kwargs)
        else:
            self.catalog = Catalogs.UCAC4(ra=self.ra, dec=self.dec, radius=size/np.sqrt(2)*1.01, **kwargs)


	'''def c1(self, image):

		return image[self.xsize - self.subarraysize:, self.ysize - self.subarraysize:]

	def c2(self, image):
		return image[0:self.subarraysize, self.ysize - self.subarraysize:]

	def c3(self, image):
		return image[0:self.subarraysize,0:self.subarraysize]

	def c4(self, image):
		return image[self.xsize - self.subarraysize:,0:self.subarraysize]'''
