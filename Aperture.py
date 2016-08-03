import numpy as np
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import logging
from settings import log_file_handler

logger = logging.getLogger(__name__)
logger.addHandler(log_file_handler)


class Aperture(object):

    def __init__(self, name=None, **kwargs):
        self.name = name

    def create(self, images, x, y, mag, plot=False):
        """Using some input images and a catalog entry, define the pixel aperture for a star."""

        # keep track of the basics of this aperture
        self.x, self.y, self.mag = x, y, mag

        # figure out which label in the labeled image is relevant to this star
        label = images['labeled'][np.round(y), np.round(x)]
        if label == 0:
            self.row, self.col = np.array([np.round(y).astype(np.int)]), np.array([np.round(x).astype(np.int)])
        else:

            # identify and sort the pixels that might contribute
            ok = images['labeled'] == label
            okr, okc = ok.nonzero()
            sorted = np.argsort(images['stars'][ok]/images['noise'][ok])[::-1]

            signal = np.cumsum(images['stars'][ok][sorted])
            noise = np.sqrt(np.cumsum(images['noise'][ok][sorted]**2))
            snr = signal/noise
            mask = ok*0
            toinclude = sorted[:np.argmax(snr)+1]
            mask[okr[toinclude], okc[toinclude]] = 1

            self.row, self.col = okr[toinclude], okc[toinclude]

            if plot:
                fi = plt.figure('selecting an aperture', figsize=(10,3), dpi=100)
                fi.clf()
                gs = gridspec.GridSpec(3,2,width_ratios=[1,.1], wspace=0.1, hspace=0.01, top=0.9, bottom=0.2)
                ax_line = plt.subplot(gs[:,0])
                ax_image = plt.subplot(gs[0,1])
                ax_ok = plt.subplot(gs[1,1])
                ax_mask = plt.subplot(gs[2,1])

                shape = ok.shape
                row, col = ok.nonzero()
                left, right = np.maximum(np.min(col)-1, 0), np.minimum(np.max(col)+2, shape[1])
                bottom, top = np.maximum(np.min(row)-1, 0),  np.minimum(np.max(row)+2, shape[1])

                kwargs = dict( extent=[left, right, bottom, top],  interpolation='nearest')
                ax_image.imshow(np.log(images['median'][bottom:top, left:right]), cmap='gray_r', **kwargs)
                ax_ok.imshow(ok[bottom:top, left:right], cmap='Blues', **kwargs)
                ax_mask.imshow(mask[bottom:top, left:right], cmap='YlOrRd', **kwargs)

                for a in [ax_ok, ax_mask, ax_image]:
                    plt.setp(a.get_xticklabels(), visible=False)
                    plt.setp(a.get_yticklabels(), visible=False)

                ax_line.plot(signal.flatten(), signal.flatten()/noise.flatten(), marker='o', linewidth=1, color='black',alpha=0.5, markersize=10)
                ax_line.set_xlabel('Total Star Flux in Aperture')
                ax_line.set_ylabel('Total Signal-to-Noise Ratio')
                ax_line.set_title('{0:.1f} magnitude star'.format(mag))
                plt.draw()
                self.input('hmmm?')
        self.n = len(self.row)
        logger.info('created {0}'.format(self))

    def __str__(self):
        return 'aperture for a {mag:.2f} magnitude star at ({x:.1f},{y:.1f})'.format(**self.__dict__)

    def measure(self, cube, plot=True):
        logger.info('calculating lightcurve for {0}'.format(self))
        shape = (len(self.row), 1)
        self.mitigated = np.sum(cube.photons[self.row, self.col, :] - cube.background[self.row, self.col].reshape(shape), 0)
        self.unmitigated = np.sum(cube.unmitigated[self.row, self.col, :] - cube.background[self.row, self.col].reshape(shape), 0)
        self.cosmics = np.sum(cube.cosmics[self.row, self.col, :], 0)
        self.withoutcosmics = self.unmitigated - self.cosmics

        # kludge for outlier rejection statistics
        self.npixelsaffected = np.sum((cube.unmitigated[self.row, self.col, :] - cube.photons[self.row, self.col, :]) > 0.1*cube.noise[self.row, self.col].reshape(shape), 0)

        # calculate various noises
        self.noises = {}
        self.noises['expected'] = np.sqrt(np.sum(cube.noise[self.row, self.col]**2))
        self.noises['achieved'] = np.std(self.mitigated)
        self.noises['unmitigated'] = np.std(self.unmitigated)
        self.noises['withoutcosmics'] = np.std(self.withoutcosmics)
        self.noises['cosmics'] = np.std(self.cosmics)
        self.median = np.median(self.withoutcosmics)
