from imports import *
from Aperture import Aperture

class Photometer(Talker):
    '''Measure the brightness of stars in an image cube.'''

    def __init__(self, cube, **kwargs):
        # decide whether or not this CCD is chatty
        Talker.__init__(self, **kwargs)


        # make handy shortcuts for all the things a photometer might like
        self.cube = cube
        self.catalog = self.cube.catalog
        self.camera = self.cube.camera
        self.cartographer = self.camera.cartographer

    def drawApertures(self, plot=False, level=5):
        '''Take all stars in the catalog, and draw apertures around them.'''


        # figure out the CCD x and y positions of the stars, and their magnitudes
        self.x, self.y = self.cartographer.point(self.catalog.ra, self.catalog.dec, 'celestial').ccdxy.tuple
        self.mag = self.catalog.tmag

        # create a noiseless master image
        self.images = {}
        self.images['median'] = self.cube.median('photons')
        self.images['stars'] = self.cube.median('photons') - self.cube.background
        self.images['noise'] = np.maximum(self.cube.noise, np.sqrt(self.images['median']))

        self.cube.display(self.images['stars'])
        self.cube.ds9.one(self.images['noise'], clobber=False)
        self.cube.ds9.one(self.images['stars']/self.images['noise'], clobber=False)

        # define regions from the star image (these are starting points for the aperture finder)
        threshold = level*self.images['noise']
        self.images['labeled'], nlabels = scipy.ndimage.measurements.label(self.images['stars'] > threshold)
        self.cube.ds9.one(self.images['labeled'], clobber=False)
        self.cube.ds9.scale(scale='log', mode='minmax')
        self.cube.ds9.match()

        # create one aperture for
        self.apertures = []
        for i in range(len(self.x)):
            ap = Aperture(i)
            try:
                ap.create(self.images, self.x[i], self.y[i], self.mag[i], plot=plot)
                self.apertures.append(ap)
            except:
                pass

    def measure(self):
        for a in self.apertures:
            a.measure(self.cube)

        noises = {}
        for k in self.apertures[0].noises.keys():
            noises[k] = np.array([a.noises[k]/a.median for a in self.apertures])
        mag = np.array([a.mag for a in self.apertures])
        return mag, noises
