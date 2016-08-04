"""Keep track of the intrapixel sensitivity variations."""

import numpy as np
import scipy.signal
import scipy.interpolate
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import logging
from settings import log_file_handler

logger = logging.getLogger(__name__)
logger.addHandler(log_file_handler)


class Intrapixel(object):
    def __init__(self, nsubpixels=240, supersample=8, name=None):
        # decide whether or not this PSF is chatty
        self.name = name
        logger.info("creating an intrapixel sensitivity map called {0}".format(self.name))
        self.nsubpixels = nsubpixels
        self.supersample = supersample

    def prnu(self, x, y):
        raise RuntimeError("No implemented for base Intrapixel class")

    def createSmoothedArray(self):
        xaxis, yaxis = np.linspace(-0.5, 0.5, self.supersample * self.nsubpixels), np.linspace(-0.5, 0.5,
                                                                                               self.supersample * self.nsubpixels)
        x, y = np.meshgrid(xaxis, yaxis)

        raw = self.prnu(x, y)
        kernal = np.ones((self.supersample, self.supersample)).astype(np.float)
        kernal /= np.sum(kernal)
        smoothed = scipy.signal.convolve2d(raw, kernal, boundary='symm', mode='same')

        interpolator = scipy.interpolate.RectBivariateSpline(xaxis, yaxis, smoothed, kx=1, ky=1)

        plt.figure('intrapixel')
        gs = gridspec.GridSpec(1, 2)
        axRaw = plt.subplot(gs[0])
        axSmoothed = plt.subplot(gs[1], sharex=axRaw, sharey=axRaw)

        kw = dict(interpolation='nearest', cmap='gray', vmin=np.min(raw), vmax=1.0, extent=[-0.5, 0.5, -0.5, 0.5])
        axRaw.imshow(raw, **kw)
        axSmoothed.imshow(interpolator(x, y), **kw)
        plt.draw()
        a = raw_input('?\n')


class Perfect(Intrapixel):
    """Perfect pixel has uniform sensitivity everywhere."""

    def __init__(self):
        super(Perfect, self).__init__(nsubpixels=0, supersample=0, name='perfectpixels')

    def prnu(self, x, y):
        """Calculate the pixel response non-uniformity for a pixel."""

        # return response
        return np.ones_like(x)


class Boxcar(Intrapixel):
    """Boxcar pixels have fixed width borders that are less sensitive."""

    def __init__(self, edgewidth=0.10, edgesensitivity=0.95, nsubpixels=240, supersample=8):
        self.edgewidth = edgewidth
        self.edgesensitivity = edgesensitivity
        super(Boxcar, self).__init__(
            nsubpixels=nsubpixels,
            supersample=supersample,
            name='{0:02.0f}border{1:02.0f}sensitive'.format(self.edgewidth * 100,
                                                            self.edgesensitivity * 100).replace('.', 'p'))

    def prnu(self, x, y):
        """Calculate the pixel response non-uniformity for a pixel."""

        # start out with a uniform response
        response = np.ones_like(x)

        # calculate distance from pixel edges
        x_distancefrompixeledge = np.abs((x % 1) - 0.5)
        y_distancefrompixeledge = np.abs((y % 1) - 0.5)

        # figure out which of input coordinates are on a pixel's border
        on_border = (x_distancefrompixeledge <= self.edgewidth) | (y_distancefrompixeledge <= self.edgewidth)

        # make the border pixels less sensitive
        response[on_border] *= self.edgesensitivity

        # return response
        return response
