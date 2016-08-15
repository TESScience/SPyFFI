import multiprocessing
import glob
import os.path
import logging

from scipy.io import loadmat
import zachopy.utils
import numpy as np
import astropy.io.fits
import scipy.signal
import matplotlib.pylab as plt

import settings
import Intrapixel
from CCD import CCD
from Cartographer import Cartographer
from settings import log_file_handler

logger = logging.getLogger(__name__)
logger.addHandler(log_file_handler)


# define everything related to PSFs
class PSF(object):
    # initialize the PSF class
    def __init__(self, camera=None,
                 version=None,
                 debprefix='woods_prf_feb2016/RAYS_ptSrc_wSi_Oct27model_AP40.6_75C_F3p314adj',
                 focus_toinclude=(0, 10), stellartemp_toinclude=(4350,),
                 nsubpixelsperpixel=101, npixels=21,
                 npositions_toinclude=21, noffsets_toinclude=11):

        # link this PSF to a Camera (hopefully one with a Cartographer)
        self.setCamera(camera)

        # the debprefix sets the basic input of data to the PSFs
        self.debprefix = debprefix
        self.version = version
        logger.info("initializing PSF painter, based on {}".format(self.debprefix))

        # the basic geometry of the unbinned pixels
        self.nsubpixelsperpixel = nsubpixelsperpixel
        self.npixels = npixels
        self.noffsets = noffsets_toinclude
        self.npositions = npositions_toinclude

        # limit the library, to keep things manageable memory-wise
        self.focus_toinclude = focus_toinclude
        self.stellartemp_toinclude = stellartemp_toinclude

        # fill in an intrapixel sensitivity
        self.intrapixel = Intrapixel.Perfect()

        # self.display = ds9('PSF')
        # self.populateJitteredPSFLibrary()
        self.populateBinned()
        # self.populateHeader()

    def findAvailable(self):
        """find the available matlab structure files from Deb"""
        self.debfiles = []
        for focus in self.focus_toinclude:
            self.debfiles.extend(
                glob.glob(os.path.join(settings.inputs,
                                       self.debprefix + '_hx*_hy*_foc{:.0f}umPRFs.mat'.format(focus))))
        logger.info(
            'there are {} files from Deb, spanning focus of {}'.format(len(self.debfiles), self.focus_toinclude))
        assert (len(self.debfiles) > 0)

    @property
    def deblibrarydirectory(self):
        f = '{:.0f}'.format
        d = os.path.join(self.versiondirectory,
                         'focus{}_stellartemp{}/'.format(
                             'and'.join(map(f, self.focus_toinclude)),
                             'and'.join(map(f, self.stellartemp_toinclude))))
        zachopy.utils.mkdir(d)
        return d

    def populateUnjitteredPSFLibrary(self):
        """populate (from scratch), a library of Deb's PRFs"""

        filename = os.path.join(self.deblibrarydirectory, 'originaldeblibrary.npy')
        try:
            self.psflibrary = np.load(filename)[0]
            logger.info('loaded original Deb library from {0}'.format(filename))
        except IOError:

            # figure out the available files to load
            self.findAvailable()

            # create an empty dictionary
            self.psflibrary = {}
            # will be indexed as focus, stellar stellartemp, fieldx, fieldy

            # loop through the files, and populate the dictionary
            for i, debfile in enumerate(self.debfiles):
                logger.info('ingesting file {} of {}'.format(i, len(self.debfiles)))
                self.ingestDebFile(debfile)

            logger.info('saving original Deb library to {}'.format(filename))
            np.save(filename, (self.psflibrary,))

        # this library hasn't been jittered by anything
        self.jitteredby = None

        # summarize what exists in the library (and define the arrays)
        self.summarizeLibrary()

        # to deal with edge effects near x=0 and y=0
        # rows (first indices are y, columns are x!)
        for focus in self.unbinned_axes['focus']:
            for stellartemp in self.unbinned_axes['stellartemp']:
                # MAKE SURE I HAVE THE TRANSPOSES RIGHT!
                # add entries that are flipped through x=0
                fieldx = np.min(self.unbinned_axes['fieldx_mm'])
                self.psflibrary[focus][stellartemp][-fieldx] = {}
                for fieldy in self.unbinned_axes['fieldy_mm']:
                    self.psflibrary[focus][stellartemp][-fieldx][fieldy] = self.psflibrary[focus][stellartemp][fieldx][
                                                                               fieldy][::1, ::-1]

                # add entries that are flipped through y=0
                fieldy = np.min(self.unbinned_axes['fieldy_mm'])
                self.psflibrary[focus][stellartemp][-fieldx][-fieldy] = self.psflibrary[focus][stellartemp][fieldx][
                                                                            fieldy][::-1, ::-1]
                for fieldx in self.unbinned_axes['fieldx_mm']:
                    self.psflibrary[focus][stellartemp][fieldx][-fieldy] = self.psflibrary[focus][stellartemp][fieldx][
                                                                               fieldy][::-1, ::1]

        self.summarizeLibrary()

    def snaptogrid(self, xvalue, yvalue):
        """Deb's models are almost but not exactly on a perfect grid in x and y
            detector coordinates. So, we have to define a grid that they can
            snap to, for the sake of making indexing work.

            Be careful, this might be introducing some overall distortions
            to the field geometry! (Or does it?)"""

        # self.physicalpixelsize is in cm, positions are in mm
        # have a position tolerance of about 2 pixels
        tolerance = 2 * self.camera.physicalpixelsize * 10.0
        # tolerance = 1.0/30.0 (kludge)
        # tolerance = 3*21.0/60.0/60.0 (from degrees)

        # find the closest in the grid, if it's within tolerance, snap to it!
        #  otherwise, add new gridpoint, which subsequent points will snap to


        try:
            self.xgrid, self.ygrid
        except AttributeError:
            self.xgrid, self.ygrid = [xvalue], [yvalue]

        def addorsnap(grid, value):
            nearest = zachopy.utils.find_nearest(np.sort(grid), value)
            offset = np.abs(value - nearest)
            if offset < tolerance:
                logger.info('snapped {:.5f} to {:.5f}'.format(value, nearest))
                return nearest
            else:
                grid.append(value)
                logger.info('added {:.5f} to grid')
                for g in grid:
                    logger.info('  {0:.5f}'.format(g))
                return value

        logger.info('snapping the x-coordinate')
        xsnapped = addorsnap(self.xgrid, xvalue)
        logger.info('snapping the y-coordinate')
        ysnapped = addorsnap(self.ygrid, yvalue)

        return xsnapped, ysnapped

    def summarizeLibrary(self):
        """print out all the keys in the Deb library"""
        keys = ['focus', 'stellartemp', 'fieldx_mm', 'fieldy_mm']
        counts = {}
        for k in keys:
            counts[k] = 0

        # create empty dictionary to store the axes along with the raw PSFs are sampled
        self.unbinned_axes = {}

        # loop through the nest of the dictionary
        self.unbinned_axes['focus'] = np.sort(self.psflibrary.keys())
        for focus in self.unbinned_axes['focus']:
            counts['focus'] += 1
            logger.info('focus={:.0f}um'.format(focus))
            self.unbinned_axes['stellartemp'] = np.sort(self.psflibrary[focus].keys())
            for stellartemp in self.unbinned_axes['stellartemp']:
                counts['stellartemp'] += 1
                logger.info('  stellartemp={:.0f}K'.format(stellartemp))
                self.unbinned_axes['fieldx_mm'] = np.sort(self.psflibrary[focus][stellartemp].keys())
                for fieldx in self.unbinned_axes['fieldx_mm']:
                    counts['fieldx_mm'] += 1
                    logger.info('    fieldx={:.3f}mm'.format(fieldx))
                    self.unbinned_axes['fieldy_mm'] = np.sort(self.psflibrary[focus][stellartemp][fieldx].keys())
                    for fieldy in self.unbinned_axes['fieldy_mm']:
                        counts['fieldy_mm'] += 1
                        logger.info('      fieldy={:.3f}mm'.format(fieldy))

        # give a quick numerical summary
        logger.info('library contains:')
        for k in keys:
            logger.info('   {} {} entries'.format(counts[k], k))

        #
        logger.info('it has been jittered by {}'.format(self.jitteredby))

    @property
    def quickloadingdirectory(self):
        d = os.path.join(self.versiondirectory, 'scratch')
        zachopy.utils.mkdir(d)
        return d

    def ingestDebFile(self, debfile):
        """load a matlab structure containing PRFs at various subpixel offsets,
            at a particular position on the detector (and focus)"""

        logger.info('loading PSF from {}'.format(debfile))

        quicksummaryfilename = os.path.join(self.quickloadingdirectory,
                                            os.path.basename(debfile) + '.summary.npy')
        quickPSFfilename = os.path.join(self.quickloadingdirectory,
                                        os.path.basename(debfile) + '.PSF.npy')
        try:
            logger.info('trying to load quickfile from {}'.format(quickPSFfilename))
            stellartemps, focuss, fieldxs, fieldys = np.load(quicksummaryfilename)
            PSFs = np.load(quickPSFfilename)
            logger.info(' ...succeeded!')
        except IOError:
            logger.info(' ...failed!')
            logger.info('loading original from {}'.format(debfile))
            mat = loadmat(debfile)['PRF_stellar']

            # pull out the stuff we're interested in
            n = mat['stellar_temp'][0, :].size
            stellartemps = [mat['stellar_temp'][0, i][0][0] for i in range(n)]
            focuss = [mat['focus_mm'][0, i][0][0] for i in range(n)]
            fieldxs = [mat['field_position'][0, i][0][0] for i in range(n)]
            fieldys = [mat['field_position'][0, i][0][1] for i in range(n)]
            PSFs = np.array([mat['PSFimage'][0, i].astype(np.float32) for i in range(n)])
            assert (np.isfinite(PSFs).all())

            # save for faster future loading
            np.save(quicksummaryfilename, (stellartemps, focuss, fieldxs, fieldys))
            np.save(quickPSFfilename, PSFs)

        for i in range(len(stellartemps)):
            stellartemp = stellartemps[i]
            if stellartemp in self.stellartemp_toinclude:
                focus = focuss[i]
                fieldx, fieldy = self.snaptogrid(fieldxs[i], fieldys[i])
                # transpose the PSF, so rows become columsn
                PSF = PSFs[i].T

                # add it to the dictionary (with cascade of trys, to define nests)
                try:
                    self.psflibrary[focus]
                except KeyError:
                    self.psflibrary[focus] = {}
                finally:
                    try:
                        self.psflibrary[focus][stellartemp]
                    except KeyError:
                        self.psflibrary[focus][stellartemp] = {}
                    finally:
                        try:
                            self.psflibrary[focus][stellartemp][fieldx]
                        except KeyError:
                            self.psflibrary[focus][stellartemp][fieldx] = {}
                        finally:
                            try:
                                self.psflibrary[focus][stellartemp][fieldx][fieldy]
                            except KeyError:
                                self.psflibrary[focus][stellartemp][fieldx][fieldy] = {}
                            finally:
                                self.psflibrary[focus][stellartemp][fieldx][fieldy] = PSF

                logger.info('added PSF at:')
                logger.info('   focus = {:.0f}'.format(focus))
                logger.info('   stellartemp = {:.0f}'.format(stellartemp))
                logger.info('   fieldx = {:.3f}'.format(fieldx))
                logger.info('   fieldy = {:.3f}'.format(fieldy))

        # self.display.one(PSF.T)
        self.setupPixelArrays()

    def setupPixelArrays(self):
        """set up pixel coordinate arrays, for raw unbinned, for intermediate binned, and for pixel-integrated"""

        try:
            self.dx_subpixels
            return
        except AttributeError:
            logger.info('setting up the pixel arrays')


            # how many subpixels are available per pixel
            self.unbinned_nsubpixelsperpixel = self.nsubpixelsperpixel
            # the size of an individual subpixel, in units of pixels
            self.unbinned_subpixelsize = 1.0 / self.unbinned_nsubpixelsperpixel
            # how big is the stamp around the PSF, in pixels
            self.unbinned_npixels = self.npixels
            # the expected total size of the PSF in subpixels
            self.unbinned_nsubpixels = self.unbinned_nsubpixelsperpixel * self.unbinned_npixels
            # the actual size of the PSF image
            self.dx_subpixelsize, self.dy_subpixelsize = self.npixels * self.nsubpixelsperpixel, self.npixels * self.nsubpixelsperpixel
            # make sure the PSF image size is what we expect
            assert (self.dx_subpixelsize == self.unbinned_nsubpixels)

            # define arrays and grid of subpixel x-y coordinates (in pixels) for the PSF image
            self.dx_subpixels_axis = np.arange(self.unbinned_nsubpixels) * self.unbinned_subpixelsize
            self.dx_subpixels_axis -= np.mean(self.dx_subpixels_axis)
            self.dy_subpixels_axis = self.dx_subpixels_axis + 0.0
            self.dx_subpixels, self.dy_subpixels = np.meshgrid(self.dx_subpixels_axis, self.dy_subpixels_axis)

            # create a grid of (rounded) full pixels that would contain at least some subpixels, centered on 0.0, 0.0
            left, right = self.subpixel2pixel([np.min(self.dx_subpixels), np.max(self.dx_subpixels)])
            bottom, top = self.subpixel2pixel([np.min(self.dy_subpixels), np.max(self.dy_subpixels)])
            self.dx_pixels_axis = np.arange(left, right + 1)
            self.dy_pixels_axis = np.arange(bottom, top + 1)
            self.dx_pixels, self.dy_pixels = np.meshgrid(self.dx_pixels_axis, self.dy_pixels_axis)

            # define arrays for the edges of the pixels
            self.dx_pixels_edges = np.zeros(self.dx_pixels_axis.size + 1)
            self.dx_pixels_edges[0:-1] = self.dx_pixels_axis - 0.5
            self.dx_pixels_edges[-1] = self.dx_pixels_edges[-2] + 1.0
            self.dy_pixels_edges = np.zeros(self.dy_pixels_axis.size + 1)
            self.dy_pixels_edges[0:-1] = self.dy_pixels_axis - 0.5
            self.dy_pixels_edges[-1] = self.dy_pixels_edges[-2] + 1.0

            logger.info('created pixel coordinate arrays')

            # create a subarray CCD with these parameters, to aid pixelization calculations
            self.ccd = CCD(camera=self.camera, subarray=self.dx_pixels.shape[0], number=1, label='PSF')
            self.cartographer = Cartographer(camera=self.ccd.camera, ccd=self.ccd)
            self.cartographer.pithy = True

    @property
    def basedirectory(self):
        """the directory in which all processed PSF data will be stored"""
        d = os.path.join(settings.intermediates, 'psfs')
        zachopy.utils.mkdir(d)
        return d

    @property
    def versiondirectory(self):
        """the directory for this particular version of the PSFs"""
        d = os.path.join(self.basedirectory, self.version)
        zachopy.utils.mkdir(d)
        return d

    @property
    def plotdirectory(self):
        """the directory where plots should be stored"""
        d = os.path.join(self.versiondirectory, 'plots')
        zachopy.utils.mkdir(d)
        return d

    def subpixel2pixel(self, anysubpixel):
        """convert from fractional pixel to an integer pixel (pixel centers are at 0.0's; edges are at 0.5's)."""
        return np.round(anysubpixel).astype(np.int)

    def setCamera(self, camera):
        """Associated a camera structure with this PSF painter."""
        self.camera = camera

    def populateHeader(self):
        """Create a PSF header structure, to store the details of the PSF simulation."""
        self.header = astropy.io.fits.Header()
        self.header[''] = ('', '')
        self.header['PRF'] = ''
        self.header['PRFNOTE'] = ('', 'Pixel Response Function library parameters')
        # self.header['PLIBDTEM'] = (self.dstellartemp, '[K] d(effective stellartemp)')
        # self.header['PLIBNOFF'] = (self.noffset, '# of subpixel offsets, in both x and y')
        # self.header['PLIBDROT'] = (self.drotation, '[deg] d(rotation around center)')
        # self.header['PLIBNRAD'] = (self.nradii, '# of radial distances from field center')
        # self.header['PSUBSIZE'] = (unbinned_subpixelsize, '[pix] subpixel size in PSF integral')
        # self.header['PNSUBPIX'] = (unbinned_nsubpixels, '[pix] # of subpixels used for initial PSF integration')
        # self.header['PPIXSIZE'] = (self.pixsize, '[pix] pixel size')

    def populateJitteredPSFLibrary(self):
        """convolve Deb's PSFs with a jittermap, at the camera's cadence"""
        logger.info('populating the jittered PSF library')
        jitteredfilename = os.path.join(
            self.deblibrarydirectory, 'jitteredlibrary_{}.npy'.format(self.camera.jitter.basename))

        try:
            self.psflibrary
            assert (self.jitteredby == self.camera.jitter.basename)
            logger.info('already loaded into memory')
            return
        except (AttributeError, AssertionError):
            logger.info('trying to load it from disk')
        try:
            # try to load the dictionary from memory
            self.psflibrary, self.jitteredby = np.load(jitteredfilename)

            # make sure the proper jittering has been applied
            assert (self.jitteredby == self.camera.jitter.basename)

            logger.info('loaded from {}'.format(jitteredfilename))

            self.summarizeLibrary()
        except (IOError,):
            logger.info('could not load it, remaking it')
            # make sure the unjittered library has been populated
            self.populateUnjitteredPSFLibrary()

            # jitter that library
            logger.info(
                'jittering all the PSFs in the library, using {} at cadence {}s'.format(self.camera.jitter.basename,
                                                                                        self.camera.cadence))
            for focus in self.unbinned_axes['focus']:
                for stellartemp in self.unbinned_axes['stellartemp']:
                    for fieldx in self.unbinned_axes['fieldx_mm']:
                        for fieldy in self.unbinned_axes['fieldy_mm']:
                            unjittered = self.psflibrary[focus][stellartemp][fieldx][fieldy]
                            kernel = self.camera.jitter.jittermap[0] / np.sum(self.camera.jitter.jittermap[0])
                            jittered = scipy.signal.convolve2d(unjittered, kernel, 'same', 'fill', 0).astype(np.float32)
                            self.psflibrary[focus][stellartemp][fieldx][fieldy] = jittered
                            logger.info(
                                'jittered focus={focus:.1f}, stellartemp={stellartemp:.0f}, fieldx={fieldx:.2f}, fieldy={fieldy:.2f}'.format(
                                    **locals()))


            # make sure to repopulate the summaries
            self.summarizeLibrary()

            # this library hasn't been jittered by anything
            self.jitteredby = self.camera.jitter.basename

            # save the jittered library
            np.save(jitteredfilename, (self.psflibrary, self.jitteredby))
            logger.info('saved jittered PSF library to {0}'.format(jitteredfilename))

    @property
    def pixelstomm(self):
        return self.camera.physicalpixelsize * 10.0

    @property
    def pixelstodeg(self):
        return self.camera.pixelscale / 3600.0

    # @profile
    def highResolutionPSF(self, position, stellartemp=5000, focus=0, chatty=False):
        """Paint a high resolution PSF at a given (cartographer-style) position,
            at a particular focus and stellar stellartemp."""

        logger.info("generating high-resolution PSF for "
                    "focus={focus}um, "
                    "T={stellartemp}K, "
                    "position={position}".format(
            position=position,
            stellartemp=stellartemp,
            focus=focus))

        # make sure a jittered library has been loaded
        try:
            self.psflibrary
            # assert(self.jitteredby == self.camera.jitter.basename)
        except (AttributeError, AssertionError):
            self.populateJitteredPSFLibrary()

        # get a cartographic object, so we can translate among coordinates
        focalx, focaly = position.focalxy.tuple

        # which available stellartemp is closest to ours?
        rounded_focus = zachopy.utils.find_nearest(self.unbinned_axes['focus'], focus)
        if chatty:
            logger.info('rounding focus={focus}um to {rounded_focus}um'.format(**locals()))

        # which available stellartemp is closest to ours?
        rounded_stellartemp = zachopy.utils.find_nearest(self.unbinned_axes['stellartemp'], stellartemp)
        if chatty:
            logger.info('rounding stellartemp={stellartemp}K to {rounded_stellartemp}K'.format(**locals()))

        library = self.psflibrary[rounded_focus][rounded_stellartemp]

        # which four focal plane positions are closest to ours?
        x = np.abs(focalx * self.pixelstomm)
        if focalx >= 0:
            xdirection = 1
        else:
            xdirection = -1
        xbounds = zachopy.utils.find_two_nearest(self.unbinned_axes['fieldx_mm'], x)
        xweightbelow, xweightabove = zachopy.utils.interpolation_weights(xbounds, x)
        xbelow, xabove = xbounds

        y = np.abs(focaly * self.pixelstomm)
        if focaly >= 0:
            ydirection = 1
        else:
            ydirection = -1
        ybounds = zachopy.utils.find_two_nearest(self.unbinned_axes['fieldy_mm'], y)
        yweightbelow, yweightabove = zachopy.utils.interpolation_weights(ybounds, y)
        ybelow, yabove = ybounds

        psf = xweightbelow * yweightbelow * library[xbelow][ybelow] + \
              xweightbelow * yweightabove * library[xbelow][yabove] + \
              xweightabove * yweightbelow * library[xabove][ybelow] + \
              xweightabove * yweightabove * library[xabove][yabove]

        if chatty:
            logger.info("treating {position} as".format(**locals()))
            logger.info(' {xweightbelow}*{yweightbelow}*library[{xbelow}][{ybelow}]'.format(**locals()))
            logger.info(' {xweightbelow}*{yweightabove}*library[{xbelow}][{yabove}]'.format(**locals()))
            logger.info(' {yweightbelow}*{yweightbelow}*library[{ybelow}][{ybelow}]'.format(**locals()))
            logger.info(' {yweightbelow}*{yweightabove}*library[{ybelow}][{yabove}]'.format(**locals()))

        # MAKE SURE THE TRANSPOSE IS RIGHT!
        return psf[::ydirection, ::xdirection] / np.sum(psf)

        # return a PSF that has the right shape, but is still centered at (0.0, 0.0)
        # return psf

    # @profile
    def binHighResolutionPSF(self, position, stellartemp=5000, dx=0.0, dy=0.0, focus=0.0, plot=True, chatty=False,
                             figure=None):
        """Pixelize a high resolution PSF at a given coordinate (cannot be used to put stars directly on images)."""

        # center the CCD subarray at the *rounded* pixel of the quoted position
        if chatty:
            logger.info('pixelizing a PSF at {0}'.format(position))

        # make sure the pixel arrays are already set up
        self.setupPixelArrays()

        self.ccd.center = np.array(np.round(position.focalxy.tuple))
        if chatty:
            logger.info('moved PSF subarray to {0}'.format(self.cartographer.ccd.center))

        # create the high-resolution PSF appropriate for this position
        zeroCenteredSubgridPSF = self.highResolutionPSF(position, stellartemp=stellartemp, focus=focus, chatty=chatty)

        # calculate the x and y subpixel grids, for both the unshifted and shifted pixels
        unshiftedx, unshiftedy = self.dx_subpixels, self.dy_subpixels
        shiftedx, shiftedy = unshiftedx + dx, unshiftedy + dy

        # MAKE SURE THERE'S NOT AN ACCIDENTAL TRANSPOSE IN HERE!
        prnu = self.intrapixel.prnu

        zeroCenteredBinnedPSF, xedges, yedges = \
            np.histogram2d(unshiftedy.flatten(), unshiftedx.flatten(), \
                           bins=[self.dy_pixels_edges, self.dx_pixels_edges], \
                           weights=zeroCenteredSubgridPSF.flatten() * prnu(unshiftedy.flatten(), unshiftedx.flatten()))

        recenteredBinnedPSF, xedges, yedges = \
            np.histogram2d(shiftedy.flatten(), shiftedx.flatten(), \
                           bins=[self.dy_pixels_edges, self.dx_pixels_edges], \
                           weights=zeroCenteredSubgridPSF.flatten() * prnu(shiftedy.flatten(), shiftedx.flatten()))

        if plot:
            plt.ioff()
            try:
                self.pixelizingfigure
            except AttributeError:
                self.pixelizingfigure = plt.figure('Pixelizing the PSF', figsize=(7.5, 7.5), dpi=100)
            plt.clf()
            gs = plt.matplotlib.gridspec.GridSpec(2, 2, wspace=0.05, hspace=0.05, top=0.80)
            self.axHighResolution = plt.subplot(gs[0, 0])
            self.axHighResolutionNudged = plt.subplot(gs[0, 1], sharex=self.axHighResolution,
                                                      sharey=self.axHighResolution)
            self.axPixelized = plt.subplot(gs[1, 0], sharex=self.axHighResolution, sharey=self.axHighResolution)
            self.axPixelizedNudged = plt.subplot(gs[1, 1], sharex=self.axHighResolution, sharey=self.axHighResolution)

            axes = [self.axHighResolution, self.axHighResolutionNudged, self.axPixelized, self.axPixelizedNudged]

            for a in axes:
                a.imshow(prnu(unshiftedx, unshiftedy), extent=extent(unshiftedx, unshiftedy), origin='lower',
                         cmap='Oranges_r', alpha=0.25, interpolation='nearest')

            kw = dict(interpolation='nearest', cmap='gray_r', alpha=0.75, origin='lower')
            # vmax=np.maximum(np.max(zeroCenteredBinnedPSF), np.max(recenteredBinnedPSF))
            vmin, vmax = np.percentile(zeroCenteredSubgridPSF, [1, 99.9])
            self.axHighResolution.imshow((zeroCenteredSubgridPSF), extent=extent(unshiftedx, unshiftedy), vmin=vmin,
                                         vmax=vmax, **kw)
            self.axHighResolutionNudged.imshow((zeroCenteredSubgridPSF), extent=extent(shiftedx, shiftedy), vmin=vmin,
                                               vmax=vmax, **kw)
            vmax = np.maximum(np.max(zeroCenteredBinnedPSF), np.max(recenteredBinnedPSF))
            self.axPixelized.imshow(zeroCenteredBinnedPSF, extent=extent(self.dx_pixels_edges, self.dy_pixels_edges),
                                    vmax=vmax, **kw)
            self.axPixelizedNudged.imshow(recenteredBinnedPSF,
                                          extent=extent(self.dx_pixels_edges, self.dy_pixels_edges), vmax=vmax, **kw)
            pixeledgekw = dict(alpha=0.1, color='green')

            # draw the pixel boundries
            '''for a in axes:
                for x in self.dx_pixels_edges:
                    a.axvline(x, **pixeledgekw)
                for y in self.dy_pixels_edges:
                    a.axhline(y, **pixeledgekw)'''
            plotsize = 3.5
            self.axHighResolution.set_xlim(-plotsize, plotsize)
            self.axHighResolution.set_ylim(-plotsize, plotsize)
            self.axHighResolution.set_title('Unnudged', fontsize=12)
            self.axHighResolutionNudged.set_title('Nudged', fontsize=12)
            self.axHighResolution.set_ylabel('High Resolution', fontsize=12)
            self.axPixelized.set_ylabel('Pixelized', fontsize=12)
            self.axPixelized.set_xlabel('(in pixels)', fontsize=9)
            self.axPixelizedNudged.set_xlabel('(in pixels)', fontsize=9)

            plt.setp(self.axHighResolution.get_xticklabels(), visible=False)
            plt.setp(self.axHighResolutionNudged.get_xticklabels(), visible=False)
            plt.setp(self.axHighResolutionNudged.get_yticklabels(), visible=False)
            plt.setp(self.axPixelizedNudged.get_yticklabels(), visible=False)
            plt.suptitle(
                "Pixelizing the TESS Point Spread Function\n{stellartemp:.0f}K star, {focus:.2f}um focus,\njittered by {jitter},\nintrapixel of {intrapixel}\n({focalx:.0f},{focaly:.0f}) pixels from focal plane center\n({dx:.2f},{dy:.2f}) from pixel center".format(
                    focus=focus, focalx=self.ccd.center[0], focaly=self.ccd.center[1], dx=dx, dy=dy,
                    stellartemp=stellartemp, jitter=self.jitteredby, intrapixel=self.intrapixel.name))
            plt.draw()

            # self.display.one(zeroCenteredSubgridPSF, frame=0)
        centralx, centraly = position.ccdxy.integerpixels

        # return the pixelized, binned, PSF
        return recenteredBinnedPSF, centralx + self.dx_pixels, centraly + self.dy_pixels

    # populate a library of binned PRFS, using the jittered high-resolution library
    # @profile
    def parallelPopulateBinned(self, plot=False, chatty=True):
        """Populate a library of binned PRFs, using the jittered, wavelength-integrated, high-resolution library."""
        self.setupPixelArrays()

        # TODO: Make function
        binned_filename = \
            os.path.join(self.deblibrarydirectory,
                         'pixelizedlibrary'
                         '_{jitter}'
                         '_{intrapixel}'
                         '_{npositions:02.0f}positions'
                         '_{noffsets:02.0f}offsets.npy'.format(
                             jitter=self.camera.jitter.basename, intrapixel=self.intrapixel.name,
                             npositions=self.npositions, noffsets=self.noffsets))

        try:
            logger.info('trying to load PSFs from {0}'.format(binned_filename))
            self.binned, self.binned_axes = np.load(binned_filename)
            logger.info('...success!')
        except IOError:
            logger.info('creating a new library of binned PSFs')
            self.populateJitteredPSFLibrary()

            # set up the grid of values over which the
            self.binned_axes = {}
            self.binned_axes['focus'] = self.unbinned_axes['focus']
            self.binned_axes['stellartemp'] = self.unbinned_axes['stellartemp']
            self.binned_axes['fieldx_px'] = np.round(
                np.linspace(-np.max(self.unbinned_axes['fieldx_mm']), np.max(self.unbinned_axes['fieldx_mm']),
                            self.npositions) / self.pixelstomm).astype(np.int)
            self.binned_axes['fieldy_px'] = np.round(
                np.linspace(-np.max(self.unbinned_axes['fieldy_mm']), np.max(self.unbinned_axes['fieldy_mm']),
                            self.npositions) / self.pixelstomm).astype(np.int)
            self.binned_axes['xoffset'] = np.linspace(-0.5, 0.5, self.noffsets)  # 11)
            self.binned_axes['yoffset'] = np.linspace(-0.5, 0.5, self.noffsets)  # 11)

            pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
            self.numberofbinnedentries = 1
            for k, v in self.binned_axes.items():
                l = len(v)
                logger.info('including {} entries for {}'.format(l, k))
                self.numberofbinnedentries *= l
            self.countthroughbinnedentries = 0
            try:
                self.binned
            except AttributeError:
                self.binned = {}
                for focus in self.binned_axes['focus']:
                    logger.info('adding focus {0:.1f}um'.format(focus), 1)
                    try:
                        self.binned[focus]
                    except KeyError:
                        self.binned[focus] = {}

                        for stellartemp in self.binned_axes['stellartemp']:
                            logger.info('adding stellartemp of {0:.0f}K'.format(stellartemp), 2)

                            try:
                                self.binned[focus][stellartemp]
                            except KeyError:
                                self.binned[focus][stellartemp] = {}

                                for fieldx in self.binned_axes['fieldx_px']:
                                    logger.info('adding focal plane x of {0:.0f} pixels'.format(fieldx), 3)
                                    try:
                                        self.binned[focus][stellartemp][fieldx]
                                    except KeyError:
                                        self.binned[focus][stellartemp][fieldx] = {}

                                        for fieldy in self.binned_axes['fieldy_px']:
                                            logger.info('adding focal plane y of {0:.0f} pixels'.format(fieldy), 4)
                                            try:
                                                self.binned[focus][stellartemp][fieldx][fieldy]
                                            except KeyError:
                                                self.binned[focus][stellartemp][fieldx][fieldy] = {}

                                            position = self.cartographer.point(fieldx, fieldy, 'focalxy')

                                            for xoffset in self.binned_axes['xoffset']:
                                                logger.info('adding xoffset of {0:.2f} pixels'.format(xoffset), 6)
                                                try:
                                                    self.binned[focus][stellartemp][fieldx][fieldy][xoffset]
                                                except KeyError:
                                                    manager = multiprocessing.Manager()
                                                    self.binned[focus][stellartemp][fieldx][fieldy][
                                                        xoffset] = manager.dict()

                                                jobs = [multiprocessing.Process(target=self.addPixelized, args=(
                                                    self.binned[focus][stellartemp][fieldx][fieldy][xoffset], position,
                                                    focus, stellartemp, xoffset, yoffset)) for yoffset in
                                                        self.binned_axes['yoffset']]
                                                # for yoffset in self.binned_axes['yoffset']:
                                                #    pool.apply_async(self.addPixelized, (self.binned[focus][stellartemp][fieldx][fieldy][xoffset], position, focus, stellartemp, xoffset, yoffset))
                                                for j in jobs:
                                                    j.start()
                                                for j in jobs:
                                                    j.join()

                                                self.countthroughbinnedentries += len(self.binned_axes['yoffset'])
                                                logger.info('')
                                                logger.info('{}/{} PSFs binned'.format(self.countthroughbinnedentries,
                                                                                       self.numberofbinnedentries))
                                                logger.info('')

                np.save(binned_filename, (self.binned, self.binned_axes))
                logger.info('saved binned PSF library to {0}'.format(binned_filename))

    def addPixelized(self, d, position, focus, stellartemp, xoffset, yoffset, plot=False, chatty=True):
        logger.info('adding yoffset of {0:.2f} pixels'.format(yoffset), 6)
        d[yoffset] = \
            self.binHighResolutionPSF(position, stellartemp=stellartemp, dx=xoffset, dy=yoffset, focus=focus, plot=plot,
                                      chatty=chatty)[0]
        if plot:
            # TODO: Make function
            binned_filename = \
                os.path.join(self.deblibrarydirectory,
                             'pixelizedlibrary'
                             '_{jitter}'
                             '_{intrapixel}'
                             '_{npositions:02.0f}positions'
                             '_{noffsets:02.0f}offsets.npy'.format(
                                 jitter=self.camera.jitter.basename, intrapixel=self.intrapixel.name,
                                 npositions=self.npositions, noffsets=self.noffsets))
            plotdir = binned_filename.replace('pixelizedlibrary_', 'plotsfor_') + '/'
            zachopy.utils.mkdir(plotdir)
            plotfile = plotdir + 'f{focus:.0f}t{stellartemp:.0f}xf{fieldx:.0f}yf{fieldy:.0f}xo{xoffset:.2f}yo{yoffset:.2f}.pdf'.format(
                **locals())
            plt.savefig(plotfile)
            logger.info('saved plot to {}'.format(plotfile))
        logger.info(
            '(focus={focus:.2f}um, T={stellartemp:.0f}K, pos={position}, dx={xoffset:.2f} pixels, dy={yoffset:.2f} pixels)'.format(
                **locals()), 6)
        return None
        # self.input('thoughts?')

    def populateBinned(self, plot=False, chatty=True):
        """Populate a library of binned PRFs, using the jittered, wavelength-integrated, high-resolution library."""

        # KLUDGE! for testing!
        # self.parallelPopulateBinned()
        # return

        self.setupPixelArrays()

        binned_filename = \
            os.path.join(self.deblibrarydirectory,
                         'pixelizedlibrary_{jitter}_{intrapixel}_'
                         '{npositions:02.0f}positions_{noffsets:02.0f}offsets.npy'.format(
                             jitter=self.camera.jitter.basename, intrapixel=self.intrapixel.name,
                             npositions=self.npositions, noffsets=self.noffsets))
        try:
            logger.info('trying to load PSFs from {0}'.format(binned_filename))
            self.binned, self.binned_axes = np.load(binned_filename)
            logger.info('...success!')
        except IOError:
            logger.info('creating a new library of binned PSFs')
            self.populateJitteredPSFLibrary()

            # set up the grid of values over which the
            self.binned_axes = {}
            self.binned_axes['focus'] = self.unbinned_axes['focus']
            self.binned_axes['stellartemp'] = self.unbinned_axes['stellartemp']
            self.binned_axes['fieldx_px'] = np.round(
                np.linspace(-np.max(self.unbinned_axes['fieldx_mm']), np.max(self.unbinned_axes['fieldx_mm']),
                            self.npositions) / self.pixelstomm).astype(np.int)
            self.binned_axes['fieldy_px'] = np.round(
                np.linspace(-np.max(self.unbinned_axes['fieldy_mm']), np.max(self.unbinned_axes['fieldy_mm']),
                            self.npositions) / self.pixelstomm).astype(np.int)
            self.binned_axes['xoffset'] = np.linspace(-0.5, 0.5, self.noffsets)  # 11)
            self.binned_axes['yoffset'] = np.linspace(-0.5, 0.5, self.noffsets)  # 11)

            self.numberofbinnedentries = 1
            for k, v in self.binned_axes.items():
                l = len(v)
                logger.info('including {} entries for {}'.format(l, k))
                self.numberofbinnedentries *= l
            self.countthroughbinnedentries = 0
            try:
                self.binned
            except AttributeError:
                self.binned = {}
                for focus in self.binned_axes['focus']:
                    logger.info('adding focus {0:.1f}um'.format(focus), 1)
                    try:
                        self.binned[focus]
                    except KeyError:
                        self.binned[focus] = {}

                        for stellartemp in self.binned_axes['stellartemp']:
                            logger.info('adding stellartemp of {0:.0f}K'.format(stellartemp), 2)

                            try:
                                self.binned[focus][stellartemp]
                            except KeyError:
                                self.binned[focus][stellartemp] = {}

                                for fieldx in self.binned_axes['fieldx_px']:
                                    logger.info('adding focal plane x of {0:.0f} pixels'.format(fieldx), 3)
                                    try:
                                        self.binned[focus][stellartemp][fieldx]
                                    except KeyError:
                                        self.binned[focus][stellartemp][fieldx] = {}

                                        for fieldy in self.binned_axes['fieldy_px']:
                                            logger.info('adding focal plane y of {0:.0f} pixels'.format(fieldy), 4)
                                            try:
                                                self.binned[focus][stellartemp][fieldx][fieldy]
                                            except KeyError:
                                                self.binned[focus][stellartemp][fieldx][fieldy] = {}

                                            position = self.cartographer.point(fieldx, fieldy, 'focalxy')

                                            for xoffset in self.binned_axes['xoffset']:
                                                logger.info('adding xoffset of {0:.2f} pixels'.format(xoffset), 6)
                                                try:
                                                    self.binned[focus][stellartemp][fieldx][fieldy][xoffset]
                                                except KeyError:
                                                    self.binned[focus][stellartemp][fieldx][fieldy][xoffset] = {}

                                                for yoffset in self.binned_axes['yoffset']:
                                                    logger.info('adding yoffset of {0:.2f} pixels'.format(yoffset), 6)
                                                    self.binned[focus][stellartemp][fieldx][fieldy][xoffset][yoffset] = \
                                                        self.binHighResolutionPSF(position, stellartemp=stellartemp,
                                                                                  dx=xoffset, dy=yoffset, focus=focus,
                                                                                  plot=plot, chatty=chatty)[0]
                                                    if plot:
                                                        plotdir = binned_filename.replace('pixelizedlibrary_',
                                                                                          'plotsfor_') + '/'
                                                        zachopy.utils.mkdir(plotdir)
                                                        plotfile = plotdir + 'f{focus:.0f}t{stellartemp:.0f}xf{fieldx:.0f}yf{fieldy:.0f}xo{xoffset:.2f}yo{yoffset:.2f}.pdf'.format(
                                                            **locals())
                                                        plt.savefig(plotfile)
                                                        logger.info('saved plot to {}'.format(plotfile))
                                                    logger.info(
                                                        '(focus={focus:.2f}um, T={stellartemp:.0f}K, pos={position}, dx={xoffset:.2f} pixels, dy={yoffset:.2f} pixels)'.format(
                                                            **locals()), 6)
                                                    self.countthroughbinnedentries += 1
                                                    logger.info(
                                                        '{}/{} PSFs binned'.format(self.countthroughbinnedentries,
                                                                                   self.numberofbinnedentries))
                                                    # self.input('thoughts?')
                np.save(binned_filename, (self.binned, self.binned_axes))
                logger.info('saved binned PSF library to {0}'.format(binned_filename))

    def comparePSFs(self, position, stellartemp=4000, verbose=False, plot=True, justnew=False, center=None):
        """Compare the PSF pulled out of the library to a newly pixelized one."""

        newpsf, newx, newy = self.newlyPixelizedPSF(position, stellartemp)
        librarypsf, libraryx, libraryy = self.pixelizedPSF(position, stellartemp)

        if plot:
            # set up the plotting figure
            plt.figure('How Good is the PSF Library?', figsize=(7.5, 7.5), dpi=100)
            plt.clf()
            gs = plt.matplotlib.gridspec.GridSpec(2, 2, wspace=0.05, hspace=0.05, top=0.8, height_ratios=[1, 0.75])
            self.axLibrary = plt.subplot(gs[0, 0])
            self.axNew = plt.subplot(gs[0, 1], sharex=self.axLibrary, sharey=self.axLibrary)
            self.axLightcurve = plt.subplot(gs[1, :])
            axes = [self.axLibrary, self.axNew]

            # if center is None:
            #    center = position

            # plot the intrapixel sensitivity
            unshiftedx, unshiftedy = self.dx_subpixels, self.dy_subpixels
            integerx, integery = center.ccdxy.integerpixels
            for a in axes:
                a.imshow(self.intrapixel.prnu(unshiftedx, unshiftedy),
                         extent=extent(unshiftedx + integerx, unshiftedy + integery), cmap='Oranges_r', alpha=0.25,
                         interpolation='nearest')

            kw = dict(interpolation='nearest', cmap='gray_r', alpha=0.75, vmin=0, vmax=0.5)
            self.axLibrary.imshow(librarypsf, extent=extent(libraryx, libraryy), **kw)
            self.axNew.imshow(newpsf, extent=extent(newx, newy), **kw)


            # draw the pixel boundries
            '''pixeledgekw = dict(alpha=0.1, color='green')
            for a in axes:
                for x in self.dx_pixels_edges:
                    a.axvline(x, **pixeledgekw)
                for y in self.dy_pixels_edges:
                    a.axhline(y, **pixeledgekw)'''

            ccdx, ccdy = center.ccdxy.tuple
            plotsize = 3.5
            self.axLibrary.set_autoscale_on(False)
            self.axLibrary.set_xlim(-plotsize + ccdx, plotsize + ccdx)
            self.axLibrary.set_ylim(-plotsize + ccdy, plotsize + ccdy)
            self.axLibrary.set_title('Interpolated from Library', fontsize=12)
            self.axNew.set_title('Just Pixelized', fontsize=12)
            for a in axes:
                plt.setp(a.get_xticklabels(), visible=False)
                plt.setp(a.get_yticklabels(), visible=False)
            dx, dy = position.ccdxy.fractionalpixels
            plt.suptitle(
                "Comparing Recently Recalculated PSF to the Library\n{stellartemp:.0f}K star, jittered for {jitter:.0f}s, intrapixel of {intrapixel}\n({focalx:.0f},{focaly:.0f}) pixels from focal plane center\n({dx:.2f},{dy:.2f}) from pixel center".format(
                    focalx=self.ccd.center[0], focaly=self.ccd.center[1], dx=dx, dy=dy, stellartemp=stellartemp,
                    jitter=self.jitteredlibrarytime, intrapixel=self.intrapixel.name))
            # plt.draw()
        return np.sum(librarypsf), np.sum(newpsf)

    def newlyPixelizedPSF(self, position, stellartemp=4000, verbose=False):
        """Wrapper for binHighResolutionPSF, to directly compare with pixelizedPSF."""
        logger.info('re-pixelizing a PSF at {0}'.format(position))
        dx, dy = position.ccdxy.fractionalpixels
        logger.info('at subpixel offsets of {0}; CCD center at {1} with size of {2}'.format((dx, dy), self.ccd.center,
                                                                                            self.ccd.npix))
        return self.binHighResolutionPSF(position, stellartemp, dx=dx, dy=dy, plot=False)

    # @profile
    def pixelizedPSF(self, position, focus=0.0, stellartemp=4000, verbose=False):
        """Drop a pixelized PSF, drawn from the library, at a particular position."""

        # make sure the binned PSF library is already loaded
        try:
            self.binned
        except AttributeError:
            self.populateBinned()


        # need to determine [radius][stellartemp][theta][xoffset][yoffset] to pull out of library
        fieldx, fieldy = position.focalxy.tuple
        key_stellartemp = zachopy.utils.find_nearest(self.binned_axes['stellartemp'], stellartemp, verbose=verbose)
        key_fieldx = zachopy.utils.find_nearest(self.binned_axes['fieldx_px'], fieldx, verbose=verbose)
        key_fieldy = zachopy.utils.find_nearest(self.binned_axes['fieldx_px'], fieldy, verbose=verbose)

        ccdx, ccdy = position.ccdxy.tuple
        xoffset, yoffset = position.ccdxy.fractionalpixels
        centralx, centraly = position.ccdxy.integerpixels

        # interpolate
        xoffsetbounds = zachopy.utils.find_two_nearest(self.binned_axes['xoffset'], xoffset, verbose=verbose)
        xbelow, xabove = xoffsetbounds
        xbelow_weight, xabove_weight = zachopy.utils.interpolation_weights(xoffsetbounds, xoffset)
        yoffsetbounds = zachopy.utils.find_two_nearest(self.binned_axes['yoffset'], yoffset, verbose=verbose)
        ybelow, yabove = yoffsetbounds
        ybelow_weight, yabove_weight = zachopy.utils.interpolation_weights(yoffsetbounds, yoffset)

        try:
            atfocus = self.binned[focus]
            needtofocusinterpolate = False
        except KeyError:
            # logger.info('could not find focus entry exactly at {}'.format(focus))
            needtofocusinterpolate = True

        if needtofocusinterpolate:
            # logger.info('interpolating in focus')
            focusbounds = zachopy.utils.find_two_nearest(self.binned_axes['focus'], focus, verbose=verbose)
            focusbelow, focusabove = focusbounds
            focusbelow_weight, focusabove_weight = zachopy.utils.interpolation_weights(focusbounds, focus)
            prfbelow = self.binned[focusbelow][key_stellartemp][key_fieldx][key_fieldy]
            prfabove = self.binned[focusabove][key_stellartemp][key_fieldx][key_fieldy]

            ibelow = xbelow_weight * ybelow_weight * prfbelow[xbelow][ybelow] + \
                     xabove_weight * ybelow_weight * prfbelow[xabove][ybelow] + \
                     xabove_weight * yabove_weight * prfbelow[xabove][yabove] + \
                     xbelow_weight * yabove_weight * prfbelow[xbelow][yabove]

            iabove = xbelow_weight * ybelow_weight * prfabove[xbelow][ybelow] + \
                     xabove_weight * ybelow_weight * prfabove[xabove][ybelow] + \
                     xabove_weight * yabove_weight * prfabove[xabove][yabove] + \
                     xbelow_weight * yabove_weight * prfabove[xbelow][yabove]

            interpolated = ibelow * focusbelow_weight + iabove * focusabove_weight
        else:
            # logger.info('evaluating at focus gridpoint')
            prf = atfocus[key_stellartemp][key_fieldx][key_fieldy]
            interpolated = xbelow_weight * ybelow_weight * prf[xbelow][ybelow] + \
                           xabove_weight * ybelow_weight * prf[xabove][ybelow] + \
                           xabove_weight * yabove_weight * prf[xabove][yabove] + \
                           xbelow_weight * yabove_weight * prf[xbelow][yabove]

        assert (interpolated is not None)

        return interpolated, centralx + self.dx_pixels, centraly + self.dy_pixels

    def magnifiedPSF(self, position, focus=0.0, stellartemp=4000, verbose=False, binby=2):
        """create a magnified PSF to plop into an image, for Kalo's reference PSFs
            right now, only puts PSFs centered at the centers of pixels (no subpixel offsets)"""

        # pull out the (zero-centered) high resolution PSF for this position
        zeroCenteredSubgridPSF = self.highResolutionPSF(position, stellartemp=stellartemp, focus=focus, chatty=False)

        originalsize = self.nsubpixelsperpixel * self.npixels
        newsize = np.round(originalsize / binby).astype(np.int)
        start = ((originalsize - newsize) % binby) / 2
        trimmed = zeroCenteredSubgridPSF[start:start + newsize * binby, start:start + newsize * binby]
        newshape = (newsize, binby, newsize, binby)
        zeroCenteredBinnedPSF = trimmed.reshape(newshape).sum(3).sum(1)

        logger.info(zeroCenteredBinnedPSF.shape)
        logger.info(newshape)
        logger.info('effective magnification is {}'.format(np.float(self.nsubpixelsperpixel) / binby))
        dx, dy = np.meshgrid(np.arange(newsize), np.arange(newsize))
        dx -= dx.mean().astype(np.int)
        dy -= dy.mean().astype(np.int)

        ccdx, ccdy = position.ccdxy.integerpixels

        return zeroCenteredBinnedPSF, dx + ccdx, dy + ccdy

    def plotPSF(self, psf, title=None, output=None):
        """Plot a PSF."""
        try:
            for a in self.ax_psf_zoom:
                a.cla()
        except:
            fig = plt.figure(figsize=(10, 10))
            plt.subplots_adjust(hspace=0, wspace=0)
            ax_map = fig.add_subplot(2, 2, 3)
            ax_vert = fig.add_subplot(2, 2, 4, sharey=ax_map)
            ax_hori = fig.add_subplot(2, 2, 1, sharex=ax_map)
            self.ax_psf_zoom = (ax_map, ax_vert, ax_hori)

        (ax_map, ax_vert, ax_hori) = self.ax_psf_zoom
        ax_hori.plot(self.dx_subpixels_axis, np.sum(psf, 0) / np.sum(psf, 0).max(), color='black', linewidth=3)
        ax_vert.plot(np.sum(psf, 1) / np.sum(psf, 1).max(), self.dy_subpixels_axis, color='black', linewidth=3)
        ax_vert.semilogx()
        ax_hori.semilogy()
        ax_hori.set_ylim(1e-6, 1e0)
        ax_vert.set_xlim(1e-6, 1e0)

        ax_vert.tick_params(labelleft=False)
        ax_hori.tick_params(labelbottom=False)
        if title is not None:
            ax_hori.set_title(title)
        ax_map.imshow(np.sqrt(psf), cmap='gray_r',
                      extent=[self.dx_subpixels_axis.min(), self.dx_subpixels_axis.max(), self.dy_subpixels_axis.min(),
                              self.dy_subpixels_axis.max()], vmin=0.0, interpolation='nearest')
        plt.draw()
        if output is not None:
            ax_map.figure.savefig(output)
            logger.info("   saved PSF plot to {0}".format(output))


def extent(x, y):
    return [x.min(), x.max(), y.min(), y.max()]
