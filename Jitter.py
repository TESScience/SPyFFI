"""Keep track of spacecraft jitter."""
import settings
import numpy as np
import astropy.table
import os.path
import zachopy.utils
import matplotlib.pylab as plt
import scipy.interpolate
import logging
import matplotlib.gridspec as gridspec
from settings import log_file_handler

logger = logging.getLogger(__name__)
logger.addHandler(log_file_handler)


def makeCartoon(seed=1):
    """Create a cartoon jitter timeseries"""

    np.random.seed(seed)
    rmsat2s = 2.0 / 3.0
    rmsat120s = 0.21

    # share across two dimensions
    rmsat120s1d = rmsat120s / np.sqrt(2)
    rmsat2s1d = rmsat2s / np.sqrt(2)

    cadence = 2.0
    smoothscale = 10.0
    nsmooth = int(smoothscale / cadence)
    t = np.arange(0, 30 * 24 * 60 * 60, 2)
    n = len(t)
    d = {}
    d['t'] = t
    for k in ['x', 'y']:
        v = np.random.normal(0, 1, n)
        for i in range(2):
            v = np.convolve(v, np.ones(nsmooth), mode='same')
        d[k] = v / np.std(v) * rmsat2s1d

    table = astropy.table.Table(d, names=['t', 'x', 'y'])
    table.write(os.path.join(settings.inputs, 'cartoon.jitter'), format='ascii.fixed_width', bookend=False)


class Jitter(object):
    def __init__(self, camera=None, jitterrms=None, rawjitterbasename="AttErrTimeArcsec_80k.dat",
                 nsubpixelsperpixel=None,
                 amplifyinterexposurejitter=1.0):


        # set an extra directory specifically for this object
        self.directory = 'jitter/'

        # store the input camera
        self.camera = camera

        # set up the initial raw jitter file
        # (this one cam from Roland, some time ago)
        self.rawfile = os.path.join(settings.inputs, rawjitterbasename)

        # what do you want the RMS to be rescaled to?
        self.jitterrms = jitterrms

        # how much should the exposure to exposure jitter be amplified
        self.amplifyinterexposurejitter = amplifyinterexposurejitter

        # (for creating map that will be used to convolve with the psf)
        self.nsubpixelsperpixel = nsubpixelsperpixel

        # update the jitterball to one that has been binned to this cadence
        self.load()

    def load(self, remake=False):
        """make sure that a jitterball (timeseries of roll,pitch,yaw) has
      been loaded and binned to the appropriate exposure times"""

        try:
            # if the jitterball is already loaded into memory
            #  *and* of the correct cadence, we're all set!
            self.jitterball

            # make sure the we're using the right jitterball for this cadence
            assert (self.jittercadence == self.camera.cadence)

            # make sure we're not trying to remake the jitterball
            assert (remake == False)

        except (AttributeError, AssertionError):
            # load the processed jitterball
            self.loadProcessedJitterball()

    @property
    def basename(self):
        cadencestatement = '.cadence{:.0f}s'.format(self.camera.cadence)

        if self.jitterrms is not None:
            jitterstatement = '.rescaledto{:0.2f}arcsec'.format(self.jitterrms)
        else:
            jitterstatement = '.unscaled'

        return os.path.basename(self.rawfile) + cadencestatement + jitterstatement

    @property
    def processedfile(self):
        """determine what the processed filename should be (based on rawfile)"""

        # store the processed files in the intermediates directory
        directory = settings.intermediates + self.directory

        # make sure a jitter directory actually exists
        zachopy.utils.mkdir(directory)

        # define the filename
        return os.path.join(directory, self.basename + '.processed.npy')

    def loadProcessedJitterball(self):
        """load a pre-processed jitterball, from the intermediates directory"""

        # if not, populate the jitterball for this cadence
        logger.info(
            'populating the jitterball for {0:.0f} second cadence, '
            'based on the raw jitter file {1}.'.format(
                self.camera.cadence,
                self.basename))

        try:
            # if a processed file already exists, load it
            self.jitterball, self.jittermap = np.load(self.processedfile)
            self.jittercadence = self.camera.cadence
        except IOError:

            logger.info('no processed jitter file was found for {}'.format(
                self.basename))

            self.loadUnprocessedJitterball()

    def loadUnprocessedJitterball(self):
        """load from a raw jitter file, process, and save"""

        # otherwise, create a binned jitter structure
        logger.info('loading raw jitter from {}'.format(self.rawfile))

        # load the raw file
        if 'AttErrTimeArcsec' in self.rawfile:
            self.rawdata = astropy.io.ascii.read(self.rawfile,
                                                 names=['t', 'x', 'y', 'z'])
        else:
            self.rawdata = astropy.io.ascii.read(self.rawfile, names=['t', 'x', 'y'])

        # subtract means
        self.rawdata['x'] -= np.mean(self.rawdata['x'])
        self.rawdata['y'] -= np.mean(self.rawdata['y'])

        # scale jitterball to requirements (should be inflation by ~1.5)
        if self.jitterrms is not None:
            # STILL A KLUDGE! NEED ROLL, PITCH, YAW!
            original_rms = np.sqrt(np.mean(self.rawdata['x'] ** 2 +
                                           self.rawdata['y'] ** 2))
            self.rawdata['x'] *= self.jitterrms / original_rms
            self.rawdata['y'] *= self.jitterrms / original_rms

        # smooth them to the required cadence
        logger.info("smoothing the jitter to {0}s cadence".format(
            self.camera.cadence))

        #  figure out the time-spacing of the jitter timeseries
        spacings = self.rawdata['t'][1:] - self.rawdata['t'][:-1]
        spacing = np.median(spacings)

        # make sure that the jitter timeseries is evenly spaced
        aboutright = 0.01
        assert ((np.abs(spacings - spacing) < spacing * aboutright).all())

        # create a convolution filter, to smooth to camera's cadence
        n = np.long(self.camera.cadence / spacing)
        filter = np.ones(n) / n

        # construct smoothed timeseries, sampled at raw time resolution
        smoothed_t = np.convolve(self.rawdata['t'], filter, mode='valid')
        smoothed_x = np.convolve(self.rawdata['x'], filter, mode='valid')
        smoothed_y = np.convolve(self.rawdata['y'], filter, mode='valid')

        # sample smoothed timeseries at the camera's cadence
        t = smoothed_t[::n]
        x = smoothed_x[::n]
        y = smoothed_y[::n]

        # plot each dimension separately
        logger.info('saving binned jitter timeseries plot')


        # create the plot of the timeseries
        plotdirectory = os.path.join(settings.plots, self.directory)
        zachopy.utils.mkdir(plotdirectory)
        bkw = dict(alpha=0.5, color='black')
        rkw = dict(linewidth=2, alpha=0.5, marker='o', color='red')
        fi, ax = plt.subplots(2, 1, sharey=True, sharex=True)
        ax[0].plot(self.rawdata['t'], self.rawdata['x'], **bkw)
        ax[0].plot(t, x, **rkw)
        ax[1].plot(self.rawdata['t'], self.rawdata['y'], **bkw)
        ax[1].plot(t, y, **rkw)
        ax[0].set_xlim(0, self.camera.cadence * 10)
        ax[0].set_title('TESS Pointing Jitter for \n{}\nfor {}s Cadence'.format(self.basename, self.camera.cadence),
                        fontsize=6)
        ax[0].set_ylabel('x (")')
        ax[1].set_ylabel('y (")')
        ax[1].set_xlabel('Time (seconds)')
        fi.savefig(os.path.join(plotdirectory,
                                self.basename + '_timeseries.pdf'))

        # make interpolators to keep track of the running smooth means
        ikw = dict(kind='nearest', fill_value=0, bounds_error=False)
        xip = scipy.interpolate.interp1d(smoothed_t, smoothed_x, **ikw)
        yip = scipy.interpolate.interp1d(smoothed_t, smoothed_y, **ikw)

        # assign the jittermap here, in units of subpixels, to be used for convolution in the PSF code
        arcsectosubpixels = 1.0 / self.camera.pixelscale * self.nsubpixelsperpixel
        xoff = (self.rawdata['x'] - xip(self.rawdata['t'])) * arcsectosubpixels
        yoff = (self.rawdata['y'] - yip(self.rawdata['t'])) * arcsectosubpixels

        npixelsfromcenter = 1
        nbins = np.maximum(np.max(np.abs(np.sqrt(xoff ** 2 + yoff ** 2))),
                           1)  # npixelsfromcenter*self.nsubpixelsperpixel
        limits = [[-nbins, nbins], [-nbins, nbins]]

        # define the jittermap as a 2D histrogram for convolution within exps
        self.jittermap = np.histogram2d(xoff, yoff,
                                        bins=nbins,
                                        range=limits,
                                        normed=True)

        # define the binned jitterball, for nudges between exps
        self.jitterball = (x, y)

        # keep track of the jitter cadence associated with this
        self.jittercadence = self.camera.cadence

        logger.info('saving jittermap plots')

        # plot the adopted jitterball, as more useful binning
        plothist2d(self.jittermap, scale=1.0 / self.nsubpixelsperpixel,
                   title='TESS Pointing Jitter over {0}s'.format(
                       self.camera.cadence),
                   xtitle='Pixels', ytitle='Pixels',
                   filename=os.path.join(plotdirectory, self.basename + '_jittermap.pdf'))

        # save the necessary jitter files
        logger.info('saving the jitter files to {0}'.format(self.processedfile))
        np.save(self.processedfile, (self.jitterball, self.jittermap))

    @property
    def x(self):
        return self.amplifyinterexposurejitter * self.jitterball[0]

    @property
    def y(self):
        return self.amplifyinterexposurejitter * self.jitterball[1]

    def writeNudges(self, outfile='jitter.txt'):

        counters = np.arange(len(self.x))
        bjds = self.camera.counterToBJD(counters)
        time = bjds - np.min(bjds)
        plt.figure('jitter timeseries')
        gs = gridspec.GridSpec(2, 1, hspace=0.15)
        kw = dict(linewidth=2)
        ax = None

        for i, what in enumerate((self.x, self.y)):
            ax = plt.subplot(gs[i], sharex=ax, sharey=ax)
            ax.plot(time, what, **kw)
            ax.set_ylabel(['dRA (arcsec)', 'dDec (arcsec)'][i])
            if i == 0:
                ax.set_title('Jitter Timeseries from\n{}'.format(self.basename))

        plt.xlabel('Time from Observation Start (days)')
        plt.xlim(np.min(time), np.max(time))
        plt.draw()
        plt.savefig(outfile.replace('.txt', '.pdf'))

        data = [counters, bjds, self.x, self.y]
        names = ['imagenumber', 'bjd', 'arcsecnudge_ra', 'arcsecnudge_dec']

        t = astropy.table.Table(data=data, names=names)
        t.write(outfile.replace('.txt', '_amplifiedby{}.txt'.format(self.amplifyinterexposurejitter)),
                format='ascii.fixed_width', delimiter=' ')
        logger.info("save jitter nudge timeseries to {0}".format(outfile))

    def applyNudge(self,
                   counter=None,  # which row to use from jitterball?
                   dx=None, dy=None,  # custom nudges, in arcsec
                   header=None,  # the FITS header in which to record nudges
                   ):

        """jitter the camera by a little bit,
      by introducing nudges draw from a
      (cadence-appropriate) jitterball timeseries."""

        # make sure the jitterball has been populated
        self.load()
        n = len(self.x)

        # should we be applying a custom offset?
        usecustom = (counter is None)
        if usecustom:
            self.camera.nudge['x'] = dx
            self.camera.nudge['y'] = dy
        else:
            # if we're over the counter, loop back
            i = counter % n
            self.camera.nudge['x'] = self.x[i]
            self.camera.nudge['y'] = self.y[i]

        # if possible, write the details to the supplied FITS header
        try:
            header['MOTION'] = ''
            header['MOTNOTE'] = ('',
                                 'properties of the image motion applied')
            header['JITTERX'] = (self.camera.nudge['x'],
                                 '["] jitter-induced nudge')
            header['JITTERY'] = (self.camera.nudge['y'],
                                 '["] jitter-induced nudge')
            header['JITPFILE'] = (self.basename,
                                  'processed jitter filename')
            header['JITSCALE'] = (self.amplifyinterexposurejitter, 'jitter magnified by ? relative to file')
            header['JITCOUNT'] = (i, 'which row of jitter file was applied?')

            logger.info('updated header keywords')
        except TypeError:
            logger.info('no header was found to update')

        # move the camera, using the updated nudge values
        logger.info("nudged the camera to {x},{y}"
                   " away from nominal pointing.".format(**self.camera.nudge))


def plothist2d(hist, title=None, log=False, scale=1.0,
               xtitle=None, ytitle=None, filename=None):
    """Plot a 2D histogram."""
    map = hist[0]
    x = (hist[1][1:] + (hist[1][0] - hist[1][1]) / 2.0) * scale
    y = (hist[2][1:] + (hist[2][0] - hist[2][1]) / 2.0) * scale
    fig = plt.figure(figsize=(10, 10))
    plt.clf()
    plt.subplots_adjust(hspace=0, wspace=0)
    ax_map = fig.add_subplot(2, 2, 3)
    ax_vert = fig.add_subplot(2, 2, 4, sharey=ax_map)
    ax_hori = fig.add_subplot(2, 2, 1, sharex=ax_map)

    ax_hori.plot(x, np.sum(map, 0) / np.sum(map), marker='o', color='black', linewidth=3)
    ax_vert.plot(np.sum(map, 1) / np.sum(map), y, marker='o', color='black', linewidth=3)
    if log:
        ax_vert.semilogx()
        ax_hori.semilogy()
    if log:
        bottom = np.min(map[map > 0]) / np.maximum(np.sum(map, 0).max(), np.sum(map, 1).max())
    else:
        bottom = 0
    top = 1
    ax_hori.set_ylim(bottom, top)
    ax_vert.set_xlim(bottom, top)

    ax_vert.tick_params(labelleft=False)
    ax_hori.tick_params(labelbottom=False)
    if title is not None:
        ax_hori.set_title(title)
    if xtitle is not None:
        ax_map.set_xlabel(xtitle)
    if ytitle is not None:
        ax_map.set_ylabel(ytitle)

    try:
        xhalf, yhalf = (x[1] - x[0]) / 2.0, (y[1] - y[0]) / 2.0
    except IndexError:
        xhalf, yhalf = 0.5, 0.5

    kw = dict(cmap='gray_r',
              extent=[x.min() - xhalf, x.max() + xhalf,
                      y.min() - yhalf, y.max() + yhalf],
              interpolation='nearest')
    if log:
        y = np.log(map)
    else:
        y = map
    ax_map.imshow(y, **kw)
    if filename is not None:
        fig.savefig(filename)
