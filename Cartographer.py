import numpy as np
import astropy.coordinates
import zachopy.borrowed.crossfield as crossfield
import logging
from settings import log_file_handler

logger = logging.getLogger(__name__)
logger.addHandler(log_file_handler)


class Cartographer(object):
    """An object to handle all conversions between coordinate systems."""

    def __init__(self, camera=None, ccd=None):
        """Initialize Cartographer, telling it where it is and how to act."""

        # These will be necessary for more complicated conversions
        self.camera, self.ccd = None, None
        self.setCamera(camera)
        self.setCCD(ccd)

    def updatePossibilities(self):
        """Update which kinds of cartographic coordinates are valid.
         Should be run every time camera or ccd are updated."""

        # no matter what, Cartographer should be able to deal with these
        possibilities = ['focalxy', 'focalrtheta']

        # if possible, figure out how to convert between focal plane coordinates and the sky
        if self.camera is not None:
            possibilities.extend(['celestial', 'ecliptic', 'galactic'])

        # if possible, figure out how to convert betweeen focal plane and pixel coordinates
        if self.ccd is not None:
            possibilities.extend(['ccdxy'])

        # update the possibility lists
        self.possibleinputs = possibilities
        self.possibleoutputs = possibilities

    def setCCD(self, ccd):
        """Tell this Cartographer to think about a particular CCD object."""
        self.ccd = ccd
        self.updatePossibilities()

    def setCamera(self, camera):
        """Tell this Cartographer to think about a particular Camera object."""
        self.camera = camera
        self.updatePossibilities()

    def point(self, a, b=None, type='focalxy'):
        """Ask Cartographer to point somewhere, using any (valid) kind of coordinate.

            inputs are *either*:
                point(position) where coord is a bonafide position object
                or
                point(a, b, type) where (a,b) are the two coordinate values and type is the kind of position they should be interpreted as"""

        # allow us to input either a coordinate object, or a pair of coordinates
        # logger.info('Cartographer pointing a new coordinate:', 1)
        try:
            assert (a.__class__.__base__.__name__ == 'position')
            type = a.__class__.__name__
            temp = a
            a, b = temp.a, temp.b
            del temp
        except:
            pass
        # logger.info('using [{a}, {b}, {type}] as input'.format(a=a,b=b, type=type), 2)


        # make sure the input coordinates can be understood
        # assert(type in self.possibleinputs)
        # logger.info('input coordinate is parsable', 2)

        # remove any previous coordinate definitions
        # logger.info('clearing previous coordinates', 2)
        for coord in self.possibleinputs:
            try:
                del self.__dict__[coord]
            except:
                pass

        # assign the coordinates to Cartographer's memory
        self.input = type
        self.coordinate = eval('{type}(a,b,cartographer=self)'.format(a=a, b=b, type=type))
        self.__dict__[self.input] = self.coordinate
        logger.debug('pointing at {0}'.format(self.coordinate))
        return self.coordinate

    def quote(self, type='focalrtheta'):
        """Ask Cartographer to say where it's pointing, using any (valid) kind of coordinate."""

        # make sure the output is possible
        # assert(type in self.possibleinputs)

        # use the coordinates property definitions to return the desired output
        output = eval('self.coordinate.{type}'.format(type=type))
        logger.info('{0} is {1}'.format(self.coordinate, output))
        return output


class position(object):
    """General (a,b) coordinate object. Each coordinate can be either a scalar, or an N-dimensional array."""

    def __init__(self, a, b, cartographer=None):
        self.a = a
        self.b = b
        assert (np.size(a) == np.size(b))
        self.cartographer = cartographer
        self.ccd = cartographer.ccd
        self.camera = cartographer.camera
        self.aname = 'a'
        self.bname = 'b'
        self.name = 'position'

    def __str__(self):
        """How should these coordinates be represented as a string?"""
        if np.size(self.a) == 1:
            return '{name} ({aname}, {bname}) = ({a},{b})'.format(**self.__dict__)
        elif np.size(self.a) > 1:
            n = np.size(self.a)
            return '{name} ({aname}, {bname}) = arrays of {n} elements'.format(n=n, **self.__dict__)
        else:
            return 'nothing!'

    @property
    def tuple(self):
        return self.a, self.b

    @property
    def arrays(self):
        """For ease of plotting, return the coordinates at a tuple of arrays."""
        return self.a, self.b

    # each type of coordinate will overwrite at least two of these base definitions (in particular focalxy and celestial)
    @property
    def focalxy(self):
        return self.celestial.focalxy

    @property
    def celestial(self):
        return self.focalxy.celestial

    @property
    def focalrtheta(self):
        return self.focalxy.focalrtheta

    @property
    def ccdxy(self):
        return self.focalxy.ccdxy

    @property
    def galactic(self):
        return self.celestial.galactic

    @property
    def ecliptic(self):
        return self.celestial.ecliptic


class focalxy(position):
    """Cartesian coordinates in the focal plane. Measured in pixels from center of field."""

    def __init__(self, a, b, cartographer=None):
        position.__init__(self, a, b, cartographer)
        self.aname, self.bname, self.name = 'x', 'y', 'focalplane'
        self.x, self.y = self.a, self.b

    # define conversions to all cartesian coordinates and to celestial
    @property
    def focalxy(self):
        return self

    @property
    def focalrtheta(self):
        radius = np.sqrt(self.x ** 2 + self.y ** 2)
        theta = (np.arctan2(self.y, self.x) + 2 * np.pi) % (2 * np.pi)
        return focalrtheta(radius, theta, self.cartographer)

    @property
    def ccdxy(self):
        xcenter, ycenter = self.ccd.center
        x = self.x - xcenter + (self.ccd.xsize - 1) / 2.0
        y = self.y - ycenter + (self.ccd.ysize - 1) / 2.0
        return ccdxy(x, y, self.cartographer)

    @property
    def celestial(self):
        ra, dec = self.camera.wcs.wcs_pix2world(self.x, self.y, 1)
        return celestial(ra, dec, self.cartographer)


class focalrtheta(position):
    """Polar coordinates in the focal plane. Measured in pixels from center of field."""

    def __init__(self, a, b, cartographer=None):
        position.__init__(self, a, b, cartographer)
        self.aname, self.bname, self.name = 'r', 'theta', 'focalplane'
        self.r, self.theta = self.a, self.b

    # use properties to define conversions (at least focalxy and self)
    @property
    def focalxy(self):
        x = self.r * np.cos(self.theta)
        y = self.r * np.sin(self.theta)
        return focalxy(x, y, self.cartographer)

    @property
    def focalrtheta(self):
        return self


class ccdxy(position):
    """Cartesian coordinates on the CCD plane. Measured in pixels from the lower left corner of the CCD."""

    def __init__(self, a, b, cartographer=None):
        position.__init__(self, a, b, cartographer)
        self.aname, self.bname, self.name = 'x', 'y', 'ccd{0:.0f}'.format(self.ccd.number)
        self.x, self.y = self.a, self.b

    @property
    def integerpixels(self):
        return np.round(self.x).astype(np.int), np.round(self.y).astype(np.int)

    @property
    def fractionalpixels(self):
        return self.x - np.round(self.x), self.y - np.round(self.y)

    # use properties to define conversions (at least focalxy and self)
    @property
    def focalxy(self):
        xcenter, ycenter = self.ccd.center
        x = self.x + xcenter - (self.ccd.xsize - 1) / 2.0
        y = self.y + ycenter - (self.ccd.ysize - 1) / 2.0

        return focalxy(x, y, self.cartographer)

    @property
    def ccdxy(self):
        return self


class celestial(position):
    """R.A. and Dec. on the sky."""

    def __init__(self, a, b, cartographer=None):
        position.__init__(self, a, b, cartographer)
        self.aname, self.bname, self.name = 'ra', 'dec', 'celestial'
        self.ra, self.dec = self.a, self.b

    # define conversions to all sky coordinates and to focalxy
    @property
    def focalxy(self):
        x, y = self.camera.wcs.wcs_world2pix(self.ra, self.dec, 1)
        xcenter, ycenter = self.camera.wcs.wcs.crpix
        return focalxy(x - xcenter, y - ycenter, self.cartographer)

    @property
    def celestial(self):
        return self

    @property
    def galactic(self):
        gal = astropy.coordinates.ICRS(ra=self.ra, dec=self.dec, unit=(astropy.units.deg, astropy.units.deg)).galactic
        l, b = gal.l.degree, gal.b.degree
        return galactic(l, b, self.cartographer)

    @property
    def ecliptic(self):
        elon, elat = crossfield.euler(self.ra, self.dec, select=3)
        return ecliptic(elon, elat, self.cartographer)


class galactic(position):
    """Galactic longitude and latitude on the sky."""

    def __init__(self, a, b, cartographer=None):
        position.__init__(self, a, b, cartographer)
        self.aname, self.bname, self.name = 'glon', 'glat', 'galactic'
        self.glon, self.glat = self.a, self.b

    # use properties to define conversions (at least celestial and self)
    @property
    def celestial(self):
        cel = astropy.coordinates.Galactic(l=self.glon, b=self.glat, unit=(astropy.units.deg, astropy.units.deg)).icrs
        ra, dec = cel.ra.degree, cel.dec.degree
        return celestial(ra, dec, self.cartographer)

    @property
    def galactic(self):
        return self


class ecliptic(position):
    """Ecliptic longitude and latitude on the sky."""

    def __init__(self, a, b, cartographer=None):
        position.__init__(self, a, b, cartographer)
        self.aname, self.bname, self.name = 'elon', 'elat', 'ecliptic'
        self.elon, self.elat = self.a, self.b

    # use properties to define conversions (at least celestial and self)
    @property
    def celestial(self):
        ra, dec = crossfield.euler(self.elon, self.elat, select=4)
        return celestial(ra, dec, self.cartographer)

    @property
    def ecliptic(self):
        return self


class ds9(position):
    def __init__(self, a, b, cartographer=None):
        position.__init__(self, a, b, Cartographer)
        self.aname, self.bname, self.name = 'x', 'y', 'ds9'

    @property
    def imshow(self):
        return imshow(self.a - 1, self.b - 1)


class imshow(position):
    def __init__(self, a, b, cartographer=None):
        position.__init__(self, a, b, Cartographer)
        self.aname, self.bname = 'x', 'y'
        self.name = 'imshow'
