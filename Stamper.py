import os
import Catalogs
import numpy as np
import astropy
import logging
from settings import log_file_handler

logger = logging.getLogger(__name__)
logger.addHandler(log_file_handler)


def dndmag(m):
    # a rough fit to a random (small) field
    f = np.exp(-3.293 + 0.697 * m)
    return f


class Stamper(object):
    def __init__(self, specifier, ccd):
        self.specifier = specifier
        self.ccd = ccd
        self.camera = self.ccd.camera


        # create stamps randomly, just with a given number per CCD
        if specifier is None:
            self.fromNothing()
        elif type(specifier) == int:
            self.fromRandom(nstamps=specifier)
        elif type(specifier) == str:
            self.fromFile(filename=specifier)

    def xy(self):
        # assign the cartrographer's CCD to this one
        self.camera.cartographer.ccd = self.ccd

        # create coordinate object for the stars
        stars = self.camera.cartographer.point(self.ra, self.dec, 'celestial')


        # return the CCD xy coordinates
        return stars.ccdxy.tuple

    def write(self):
        # create an astropy table
        t = astropy.table.Table([self.ra, self.dec, self.radii], names=('ra', 'dec', 'radius'))

        # write it out
        f = os.path.join(
            self.ccd.directory, 'postagestamptargets_{pos}_{name}.txt'.format(pos=self.ccd.pos_string,
                                                                              name=self.ccd.name))

        t.write(f, format='ascii.fixed_width', delimiter='|', bookend=False)
        logger.info('wrote stamp definition catalog to {}'.format(f))

    def fromNothing(self):
        logger.info('no postage stamps defined')
        self.ccd.stampimage = self.ccd.ones()

    def fromRandom(self, nstamps, radius=10):
        logger.info('randomly populating {} with {} postage stamps'.format(self.ccd.name, nstamps))

        # select (randomly) some target stars, seeded by the CCD number
        np.random.seed(self.ccd.number)

        ras, decs, tmag, temperatures = self.camera.catalog.snapshot(self.camera.bjd,
                                                                     exptime=self.camera.cadence / 60.0 / 60.0 / 24.0)

        # assign the cartrographer's CCD to this one
        self.camera.cartographer.ccd = self.ccd

        # create coordinate object for the stars
        stars = self.camera.cartographer.point(ras, decs, 'celestial')
        x, y = stars.ccdxy.tuple

        # start with targets with centers on the chip
        onccd = (np.round(x) > 0) & \
                (np.round(x) < self.ccd.xsize) & \
                (np.round(y) > 0) & \
                (np.round(y) < self.ccd.ysize)

        # weight stars inversely to their abundance, to give roughly uniform distribution of magnitudes
        weights = \
            (1.0 / dndmag(self.camera.catalog.tmag) * (self.camera.catalog.tmag >= 6) * (
            self.camera.catalog.tmag <= 16))[
                onccd]
        weights /= np.sum(weights)
        itargets = np.random.choice(onccd.nonzero()[0],
                                    size=np.minimum(nstamps, np.sum(weights != 0)),
                                    replace=False,
                                    p=weights)

        # populate position arrays
        self.ra = self.camera.catalog.ra[itargets]
        self.dec = self.camera.catalog.dec[itargets]
        self.radii = np.ones_like(self.ra) * radius

        self.finishFromStars()

    def finishFromStars(self):

        # convert to ccd coordinates
        self.x, self.y = self.xy()

        # write out the stamp catalog
        self.write()

        # populate the stamp image
        self.populateStampImage()

    def fromFile(self, filename):
        logger.info('loading RA and Dec stamp centers from {}'.format(filename))
        self.table = astropy.io.ascii.read(filename, names=['ra', 'dec', 'radius'])

        self.ra = self.table['ra'].data
        self.dec = self.table['dec'].data
        self.radii = self.table['radius'].data

        self.finishFromStars()

    def trimCatalog(self, catalog):
        '''trim a catalog to contain only stars that fall on the CCD and in a stamp'''

        # assuming stars don't move in and out of postage stamps over time
        ras, decs, tmag, temperatures = self.camera.catalog.snapshot(self.camera.bjd,
                                                                     exptime=self.camera.cadence / 60.0 / 60.0 / 24.0)

        # assign the cartrographer's CCD to this one
        self.camera.cartographer.ccd = self.ccd

        # create coordinate object for the stars
        stars = self.camera.cartographer.point(ras, decs, 'celestial')
        x, y = stars.ccdxy.integerpixels

        onccd = (x >= 0) * (x < self.ccd.xsize) * (y >= 0) * (y < self.ccd.ysize)
        ok = np.ones_like(onccd)
        ok[onccd] *= self.ccd.stampimage[y[onccd], x[onccd]].astype(np.bool)

        return Catalogs.Trimmed(self.camera.catalog, ok)

    def populateStampImage(self):

        x = self.x
        y = self.y
        radii = self.radii

        # create an empty array
        mask = self.ccd.zeros()

        # create a (currently fixed size and shape) aperture
        diameter = max(radii) * 2 + 1
        xgrid, ygrid = np.meshgrid(np.arange(-max(radii), max(radii) + 1), np.arange(-max(radii), max(radii) + 1))

        # make a circular aperture
        r = np.sqrt(xgrid ** 2 + ygrid ** 2)

        # loop over stars
        for i in range(len(x)):
            if (x[i] >= 0) * (x[i] < self.ccd.xsize) * (y[i] >= 0) * (y[i] < self.ccd.ysize):
                thisx = np.minimum(np.maximum(np.round(x[i]).astype(np.int) + xgrid, 0), self.ccd.xsize - 1).astype(
                    np.int)
                thisy = np.minimum(np.maximum(np.round(y[i]).astype(np.int) + ygrid, 0), self.ccd.ysize - 1).astype(
                    np.int)
                aperture = r <= radii[i]
                mask[thisy, thisx] += aperture

        self.ccd.stampimage = mask
        self.ccd.note = 'stampdefinition'
        stampfilename = os.path.join(self.ccd.directory, self.ccd.note + '.fits')
        self.ccd.writeToFITS(self.ccd.stampimage, stampfilename)
