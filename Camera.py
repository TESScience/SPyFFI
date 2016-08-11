import logging
import os

import numpy as np
import zachopy.utils
import astropy.io.fits
import zachopy.spherical

import settings
import Catalogs
from PSF import PSF
from Cartographer import Cartographer
from CCD import CCD
from Jitter import Jitter
from Focus import Focus
from settings import log_file_handler

logger = logging.getLogger(__name__)
logger.addHandler(log_file_handler)


# define a camera class
class Camera(object):
    """Keep track of one camera's entire field of view."""

    def __init__(self,
                 cadence=1800,  # what cadence for exposures?
                 ra=270, dec=66.56070833333332,  # field center,
                 testpattern=False,  # is this using a test pattern?
                 subarray=None,  # define a subarray at field center?
                 label='',  # a special label for this field?
                 cameranumber=1,  # which camera is this? (not used yet)
                 warpspaceandtime=False,  # slow the speed of light?
                 counterstep=1,  # jump by multiple exposures each time?
                 aberrate=True,  # apply velocity aberration?
                 positionangle=None,  # position angle of the field
                 stamps={2: None, 120: None, 1800: None},  # how many postage stamps?
                 variablefocus=False,
                 dirprefix='',
                 psfkw={},
                 jitterkw={},
                 focuskw={}
                 ):
        """Initialize camera, fill it with CCDs, and point it at the sky or at a testpattern."""

        self.psfkw, self.jitterkw, self.focuskw = psfkw, jitterkw, focuskw

        # in case you want everything stored in its own directory
        self.dirprefix = dirprefix
        zachopy.utils.mkdir(os.path.join(settings.outputs, self.dirprefix))

        # KLUDGE (or is it?)
        self.stamps = stamps

        # keep track of which exposure is being simulated
        self.counter = 0

        # [logger.info will only report for this object if mute and pithy are turned off]
        logger.info("turning on a new TESS camera object.")

        # keep track of any special features with this Camera
        self.subarray = subarray
        self.label = label

        # will this camera see a test pattern?
        self.testpattern = testpattern
        logger.info('the camera expects the stars to be drawn from {}'.format(
            {True: 'the sky', False: 'a test pattern'}[self.testpattern]))

        # warping time (for enhanced aberration)
        self.warpspaceandtime = warpspaceandtime
        if self.warpspaceandtime:
            logger.info('the camera will warp space and time by slowing '
                        'the speed of light to {}c'.format(self.warpspaceandtime))

        # should positions be aberrated?
        self.aberrate = aberrate

        # should the focus be allowed to vary?
        self.variablefocus = variablefocus

        # speeding up time
        self.counterstep = counterstep
        if self.counterstep > 1:
            logger.info('the camera will speed up time'
                        ' by a factor of {}'.format(self.counterstep))


        # if real stars, use the input (ra, dec)
        self.ra = ra
        self.dec = dec
        logger.info('the camera FOV is centered at (ra,dec) = {:.2f}, {:.2f} deg.'.format(
            self.ra, self.dec))


        # assign the cadence for this camera
        self.singleread = 2.0  # seconds
        self.readouttime = 0.05  # seconds

        # define scales for the Camera
        self.pixelscale = 21.1  # arcsec!! (from Peter's paper)
        self.entrance_pupil_diameter = 10.5  # cm (from Peter's paper)
        self.effective_area = 69.1  # cm^2 (from Peter's paper) should be 63.0
        self.physicalpixelsize = 15.0 / 1e4  # cm
        self.physicalpixeldepth = 100.0 / 1e4  # cm
        self.read_noise = 10.0  # electrons per read
        self.saturation = 150000.0  # electrons per read

        # define the gaps between the CCDs (numbers still in flux, as of late 2014)
        self.physicalgap = 0.2  # cm
        self.gapinpixels = self.physicalgap / self.physicalpixelsize  # pixels (not necessarily integer)

        # define the field of view of the Camera (one side of the CCD square)
        self.fov = (4096 + self.gapinpixels) * self.pixelscale / 60.0 / 60.0  # degrees
        self.pixel_solid_area = (self.pixelscale) ** 2  # arcsec^2

        # turn on the necessary CCDs in this Camera
        if self.subarray is None:
            # if we're not dealing with a subarray, then turn on CCD's 1,2,3,4
            self.ccdnumbers = np.arange(4) + 1
            logger.info("populating camera with 4 CCDs")
        else:
            # if this is a subarray, then turn on one (imaginary) CCD and call it 0
            self.ccdnumbers = np.arange(1)
            logger.info("populating camera with one, centered, CCD subarray")
        self.ccds = [CCD(n, subarray=self.subarray, camera=self) for n in self.ccdnumbers]

        # assign a cartographer to this Camera and start it out on the first CCD
        self.cartographer = Cartographer(camera=self, ccd=self.ccds[0])

        # start the camera out unjittered from its nominal position
        self.nudge = {'x': 0.0, 'y': 0.0, 'z': 0.0}  # nudge relative to nominal spacecraft pointing (arcsec)

        self.setCadence(cadence)  # seconds

        # point the Camera
        self.point(self.ra, self.dec)

    @property
    def fielddirectory(self):
        """define the field directory for this camera"""
        if self.label == '':
            d = os.path.join(settings.outputs,
                             '{dirprefix}{pos}'.format(pos=self.pos_string(), dirprefix=self.dirprefix))
        else:
            d = os.path.join(settings.outputs,
                             '{dirprefix}{pos}_{label}'.format(pos=self.pos_string(),
                                                               label=self.label,
                                                               dirprefix=self.dirprefix))
        zachopy.utils.mkdir(d)
        return d

    @property
    def directory(self):
        d = os.path.join(self.fielddirectory, '{0:d}s/'.format(self.cadence))
        zachopy.utils.mkdir(d)
        return d

    def expose(self, **kwargs):
        """Take an exposure on all the available CCD's."""
        for c in self.ccds:
            c.expose(**kwargs)

    def populateHeader(self):
        """Populate the header structure with information about the Camera, and its WCS."""

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
        self.header['FOCUS'] = (None, 'distance from optimal focus (microns)')
        if self.subarray is not None:
            self.header['SUBARRAY'] = (
                self.subarray, 'THIS IMAGE IS JUST {0}x{1} POSTAGE STAMP!'.format(self.subarray, self.subarray))

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
        """Set the cadence of this Camera, in seconds."""

        # set the cadence
        self.cadence = cadence
        logger.info(
            "setting cadence to {0} seconds = {1:.0f} reads.".format(self.cadence, self.cadence / self.singleread))

        #
        self.focus = Focus(camera=self, **self.focuskw)

        # load the jitterball for this camera
        self.jitter = Jitter(camera=self, nsubpixelsperpixel=self.psfkw['nsubpixelsperpixel'], **self.jitterkw)

        # load the PSF for this Camera
        self.psf = PSF(camera=self, **self.psfkw)

        # make sure the background image gets reset
        for c in self.ccds:
            try:
                del c.backgroundimage
            except AttributeError:
                pass

    def counterToBJD(self, counter):
        self.bjd0 = 2457827.0
        return self.bjd0 + counter * self.cadence / 24.0 / 60.0 / 60.0

    def point(self, ra=None, dec=None):
        """Point this Camera at the sky, by using the field-specified (ra,dec) and (if active) the jitter nudge for this exposure."""

        # if (ra,dec) given, then point the whole Camera there
        if ra is not None and dec is not None:
            self.ra, self.dec = ra, dec

        # make sure an (ra,dec) are defined
        try:
            self.ra
            self.dec
            logger.info('pointing the camera at (ra,dec) = {0:.6f},{1:.6f}'.format(self.ra, self.dec))
        except:
            self.report("Please point your telescope somewhere. No RA or DEC defined.")

        # create a blank WCS object, for converting between (ra,dec) and (x,y)
        self.wcs = astropy.wcs.WCS(naxis=2)

        # the focalxy coordinates of the reference position [taken to be center of field, which is (x,y) = (0.0, 0.0) in focalxy coordinates]
        self.wcs.wcs.crpix = [0.0, 0.0]

        # the pixel scale, in degrees
        self.wcs.wcs.cdelt = [-self.pixelscale / 60.0 / 60.0, self.pixelscale / 60.0 / 60.0]



        # the celestial coordinates at the reference position (input by user)
        nudged_ra, nudged_dec = zachopy.spherical.rotate(self.ra, self.dec, self.nudge['x'] / 60.0 / 60.0,
                                                         self.nudge['y'] / 60.0 / 60.0)
        self.wcs.wcs.crval = [nudged_ra, nudged_dec]

        # the rotation of the field (currently just a pure rotation, no shear)
        # rot = self.nudge['z']/60.0/60.0*np.pi/180.0
        # w.wcs.pc = [[np.cos(rot), -np.sin(rot)],[np.sin(rot), np.cos(rot)]]

        # the coordinate system type - what should I use?
        self.wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

        # set this to be the WCS
        # self.populateHeader()

    def pos_string(self):
        """Return the position string for this field."""
        if self.testpattern:
            try:
                return self.catalog.name
            except:
                return 'testpattern'
        else:
            coords = astropy.coordinates.ICRS(ra=self.ra * astropy.units.degree, dec=self.dec * astropy.units.degree)
            return "{0:02}h{1:02}m{2:02}s{3:+03}d{4:02}m{5:02}s".format(np.int(coords.ra.hms[0]),
                                                                        np.int(coords.ra.hms[1]),
                                                                        np.int(coords.ra.hms[2].round()),
                                                                        np.int(coords.dec.dms[0]),
                                                                        np.int(np.abs(coords.dec.dms[1])),
                                                                        np.int(np.abs(coords.dec.dms[2].round())))

    def advanceCounter(self):
        """Take one step forward in time with this Camera."""
        self.counter += self.counterstep

    def populateCatalog(self, **kwargs):
        """Create a catalog of stars that are visible with this Camera."""

        # figure out how wide we need to search for stars
        if self.subarray is None:
            size = self.fov
        else:
            size = self.subarray * self.pixelscale / 3600.0

        # make a test pattern or a catalog of real stars
        if self.testpattern:
            # figure out how big a test pattern to create (in arcsec)
            self.catalog = Catalogs.TestPattern(size=size * 3600.0, **kwargs)
        else:
            self.catalog = Catalogs.UCAC4(ra=self.ra, dec=self.dec, radius=size / np.sqrt(2) * 1.01, **kwargs)

    @property
    def effective_fov(self):
        """return the radius (degrees) needed to reach the detector corners"""

        # is it a subarray or no?
        if self.subarray is None:
            # set the radius to be the actual radial FOV of the camera
            radius = self.fov / np.sqrt(2)
        else:
            # set the radius to the corner of the subarray
            radius = self.subarray * self.pixelscale / 60.0 / 60.0

        return radius

        '''def c1(self, image):

            return image[self.xsize - self.subarraysize:, self.ysize - self.subarraysize:]

        def c2(self, image):
            return image[0:self.subarraysize, self.ysize - self.subarraysize:]

        def c3(self, image):
            return image[0:self.subarraysize,0:self.subarraysize]

        def c4(self, image):
            return image[self.xsize - self.subarraysize:,0:self.subarraysize]'''
