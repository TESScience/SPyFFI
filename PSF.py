from imports import *
import settings
import Intrapixel
from CCD import CCD
from Cartographer import Cartographer
# define everything related to PSFs
class PSF(Talker):

    # initialize the PSF class
    def __init__(self, camera=None):

        # decide whether or not this PSF is chatty
        Talker.__init__(self, mute=False, pithy=False)
        self.speak("Initializing the TESS point spread function painter.")

        # link this PSF to a Camera (hopefully one with a Cartographer)
        self.setCamera(camera)


        # parameters used for generating the binned PSF library
        self.version = 'original'
        self.dtemperature = 10000    # the step in temperature, stamper will choose nearest
        self.noffset = 21#5           # number of sub-pixel offsets, stamper will interpolate between offsets
        self.drotation = 10#90         # number of rotation angles for PSF's (0-360), stamper will ???
        self.nradii = 10           # number of radial focalplaneradius for PSF's, stamper will ???

        # create the necessary x,y pixel grids for painting PSFs (at both high and low resolution)
        self.setupPixelArrays()

        # fill in an intrapixel sensitivity
        self.intrapixel = Intrapixel.Perfect(nsubpixels=self.numberof_subpixelsforintegrating)

        # if possible, associate a Camera object with this PSF
        self.setCamera(camera)
        #self.populateBinned()
        self.populateHeader()


    @property
    def versiondirectory(self):
        d = settings.prefix + 'intermediates/' + self.version + '/'
        zachopy.utils.mkdir(d)
        return d


    @property
    def plotdirectory(self):
        d = settings.prefix + 'plots/' + self.version + '/'
        zachopy.utils.mkdir(d)
        return d


    def setupPixelArrays(self):
        '''Set up the pixel coordinate arrays, neede both for integrating the high resolution PSFs and for painting on binned PSFs.'''

        # try to load the pixel coordinate arrays from saved file (this done not for speedup, but to ensure that all subsequent work will be using the correct arrays for this version)
        arraysfilename = self.versiondirectory + 'pixelarrays.npy'
        try:
            # try loading these from a file
            self.initial_binning, self.subpixelsforintegrating_size, self.numberof_subpixelsforintegrating  = np.load(arraysfilename)
            self.speak('loaded PSF subpixel array definitions from {0}'.format(arraysfilename))
        except:
            # otherwise, create them from scratch
            if self.version == 'original':
                # details related to the orginal PSF files from Deb, by way of Peter.
                self.initial_binning = 4                            # how much to bin the original (very high resolution, from Deb) PSFs before integrating over the pixels
                self.subpixelsforintegrating_size = 0.25/15.0*self.initial_binning    # with this binning, how big is each subpixelsforintegratingel in the high resolution PSF?
                self.numberof_subpixelsforintegrating = 960/self.initial_binning             # how many subpixelsforintegratingels are there in the high resolution PSF image?
            # save these arrays for future use
            np.save(arraysfilename, (self.initial_binning, self.subpixelsforintegrating_size, self.numberof_subpixelsforintegrating))
            self.speak('save PSF subpixel array definitions to {0}'.format(arraysfilename))

        # define grids of subpixel offsets (positive and negative, centered on 0.0, 0.0), over which the PSF integration will be carried out
        # (how far are the *centers* of the subpixels from the star's focalplaneradius?)
        self.dx_subpixelsforintegrating_axis = np.arange(self.numberof_subpixelsforintegrating)*self.subpixelsforintegrating_size
        self.dx_subpixelsforintegrating_axis -= np.mean(self.dx_subpixelsforintegrating_axis)
        self.dy_subpixelsforintegrating_axis = self.dx_subpixelsforintegrating_axis
        self.dx_subpixelsforintegrating, self.dy_subpixelsforintegrating = np.meshgrid(self.dx_subpixelsforintegrating_axis, self.dy_subpixelsforintegrating_axis)

        # create a grid of full pixels that would contain at least some subpixels, centered on 0.0, 0.0
        left, right = self.subpixel2pixel((np.min(self.dx_subpixelsforintegrating), np.max(self.dx_subpixelsforintegrating)))
        bottom, top = self.subpixel2pixel((np.min(self.dy_subpixelsforintegrating), np.max(self.dy_subpixelsforintegrating)))
        self.dx_pixels_axis = np.arange(left, right+1)
        self.dy_pixels_axis = np.arange(bottom, top+1)
        self.dx_pixels, self.dy_pixels = np.meshgrid(self.dx_pixels_axis, self.dy_pixels_axis)

        # define arrays for the edges of the pixels
        self.dx_pixels_edges = np.zeros(self.dx_pixels_axis.size+1)
        self.dx_pixels_edges[0:-1] = self.dx_pixels_axis - 0.5
        self.dx_pixels_edges[-1] = self.dx_pixels_edges[-2] + 1.0
        self.dy_pixels_edges = np.zeros(self.dy_pixels_axis.size+1)
        self.dy_pixels_edges[0:-1] = self.dy_pixels_axis - 0.5
        self.dy_pixels_edges[-1] = self.dy_pixels_edges[-2] + 1.0


        # create a subarray CCD with these parameters, to aid pixelization calculations
        self.ccd = CCD(camera=self.camera, subarray=self.dx_pixels.shape[0], number=1, label='PSF')
        self.cartographer = Cartographer(camera=self.ccd.camera, ccd=self.ccd)
        self.cartographer.pithy = True

    def subpixel2pixel(self,anysubpixel):
        '''Convert from fractional pixel to an integer pixel (pixel centers are at 0.0's; edges are at 0.5's).'''
        return np.round(anysubpixel).astype(np.int)


    def setCamera(self, camera):
        '''Associated a camera structure with this PSF painter.'''
        self.camera = camera


    def populateHeader(self):
        '''Create a PSF header structure, to store the details of the PSF simulation.'''
        self.header = astropy.io.fits.Header()
        self.header[''] = ('', '')
        self.header['PRF'] = ''
        self.header['PRFNOTE'] = ('','Pixel Response Function library parameters')
        self.header['PLIBDTEM'] = (self.dtemperature, '[K] d(effective temperature)')
        self.header['PLIBNOFF'] = (self.noffset, '# of subpixel offsets, in both x and y')
        self.header['PLIBDROT'] = (self.drotation, '[deg] d(rotation around center)')
        self.header['PLIBNRAD'] = (self.nradii, '# of radial distances from field center')
        self.header['PSUBSIZE'] = (self.subpixelsforintegrating_size, '[pix] subpixel size in PSF integral')
        #self.header['PNSUBPIX'] = (self.numberof_subpixelsforintegrating, '[pix] # of subpixels used for initial PSF integration')
        #self.header['PPIXSIZE'] = (self.pixsize, '[pix] pixel size')


    def populateHighResolutionMonochromaticLibrary(self, plot=True, inspect=False):
        '''Load the high resolution PSFs from Deb, which are monochromatic and unjittered.'''

        self.speak("populating the raw, monochromatic, PSF library")
        monochromatic_filename = self.versiondirectory + 'monochromatic_psfs.npy'
        try:
            # maybe we can load from a preprocessed file?
            (self.monochromaticlibrary, self.focalplaneradii, self.wavelengths)  = np.load(monochromatic_filename)
            self.speak("loaded high-resolution monochromatic library from {0}".format(monochromatic_filename))
        except:
            # load the high resolution PSF's from a FITs file
            deb_filename = settings.prefix + 'inputs/psfs_3p33.fits'
            data = astropy.io.fits.open(deb_filename)[0].data
            self.speak("loaded high-resolution monochromatic library from {0}".format(deb_filename))

            # trim the high resolution PSF, so that an even number of whole pixels is included
            assert(data.shape[2] == data.shape[3])
            binned = zachopy.utils.rebin(data,data.shape[0], data.shape[1], data.shape[2]/self.initial_binning, data.shape[3]/self.initial_binning)
            ntrim = (binned.shape[2] - self.numberof_subpixelsforintegrating)/2
            data = binned[:,:,ntrim:-ntrim,ntrim:-ntrim]

            # define some useful counts
            nfocalplaneradii, nwavelengths, nxs, nys = data.shape
            assert(nxs == nys)

            # create an empty monochromatic PSF library
            self.monochromaticlibrary = {}

            # input the focalplaneradius and wavelength coordinates that Deb used
            focalplaneradii = np.array([0,0.5,1.0,np.sqrt(2)])*2048
            wavelengths = 0.675 + np.arange(nwavelengths)*0.05

            # keep track of these
            self.focalplaneradii = focalplaneradii
            self.wavelengths = wavelengths

            if plot:
                # create a grid of axes and set up the coordinate system
                fi, ax_psf_fits  = plt.subplots(data.shape[0],data.shape[1], figsize=(25.95, 9.8), sharex=True, sharey=True)
                extent=[self.dx_subpixelsforintegrating_axis.min(), self.dx_subpixelsforintegrating_axis.max(), self.dy_subpixelsforintegrating_axis.min(), self.dy_subpixelsforintegrating_axis.max()]

            # loop over focalplaneradii
            for i in range(nfocalplaneradii):
                # populate a dictionary of dictionaries
                self.monochromaticlibrary[focalplaneradii[i]] = {}

                # loop over wavelengths
                for j in range(nwavelengths):

                    # store monochromatic PSF at right place in library
                    self.monochromaticlibrary[focalplaneradii[i]][wavelengths[j]] = data[i,j,:,:]


                    if plot:
                        if inspect:
                            try:
                                # maybe we can clear and use existing axes
                                for a in ax_psf_zoom:
                                    a.cla()
                            except:
                                # create figure
                                fig = plt.figure(figsize=(10,10))
                                plt.subplots_adjust(hspace=0, wspace=0)
                                ax_map = fig.add_subplot(2,2,3)
                                ax_vert = fig.add_subplot(2,2,4, sharey=ax_map)
                                ax_hori = fig.add_subplot(2,2,1, sharex=ax_map)
                                ax_psf_zoom = [ax_map, ax_vert, ax_hori]

                            # plot a big zoom of each monochromatic PSF, with histograms
                            ax_hori.plot(self.dx_subpixelsforintegrating_axis, np.sum(self.monochromaticlibrary[focalplaneradii[i]][wavelengths[j]], 0))
                            ax_vert.plot(np.sum(self.monochromaticlibrary[focalplaneradii[i]][wavelengths[j]], 1), self.dy_subpixelsforintegrating_axis)
                            ax_map.imshow(np.sqrt(self.monochromaticlibrary[focalplaneradii[i]][wavelengths[j]]),extent=extent, cmap='gray_r', interpolation='nearest', vmin=0  )
                            plt.draw()

                        # populate a plot in the library grid
                        ax_psf_fits[i,j].imshow(np.sqrt(self.monochromaticlibrary[focalplaneradii[i]][wavelengths[j]]), extent=extent, cmap='gray_r', interpolation='nearest', vmin=0 )
                        if i == data.shape[0]-1:
                            ax_psf_fits[i,j].set_xlabel("{0} micron".format(wavelengths[j]))

                        if j == 0:
                            ax_psf_fits[i,j].set_ylabel("{0:4.1f} pix.".format(focalplaneradii[i]))
                        plt.subplots_adjust(hspace=0, wspace=0)
            if plot:
                plt.draw()
                plotfilename = self.plotdirectory + 'monochromatic_psf_library.pdf'
                plt.savefig(plotfilename)
                self.speak('saved plot of high-resolution, monochromatic libarary to {0}'.format(plotfilename))

            # save the library, so it'd be easier to load next time
            np.save(monochromatic_filename,(self.monochromaticlibrary, self.focalplaneradii, self.wavelengths))
            self.speak('saved high-resolution monochromatic to {0}'.format(monochromatic_filename))


    def populateHighResolutionWavelengthIntegratedLibrary(self, plot=True):
        '''Integrate the monochromatic PSFs into a (still unjittered) wavelength-integrated PSF library.'''

        self.speak("populating the wavelength-integrated wavelength-integrated, high-resolution PSF library")
        temperature_filename = self.versiondirectory + 'wavelengthintegrated_psfs.npy'
        try:
            # maybe we can load a preprocessed file?
            self.temperaturelibrary, self.temperatures, self.focalplaneradii = np.load(temperature_filename)
            self.speak('loaded wavelength-integrated, high-resolution PSF library from {0}'.format(temperature_filename))
        except:
            self.speak('creating a wavelength-integrated, high-resolution PSF library')
            # otherwise, we have to recalculate from scratch
            try:
                # make sure the monochromatic library has been populated
                self.monochromaticlibrary
            except:
                self.populateHighResolutionMonochromaticLibrary()

            # load weight coefficients calculated by Peter
            weighting = astropy.io.fits.open(settings.prefix + 'inputs/ph_filt.fits')[0].data

            # create a grid of temperatures at which PSFs will be calculated
            self.temperatures = np.arange(4000, 12001, self.dtemperature)

            # define a handy function for plotting a PSF
            if plot:
                def plotPSF(psf, title=None, output=None):
                    '''Plot a PSF.'''
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
                    ax_hori.plot(self.dx_subpixelsforintegrating_axis, np.sum(psf, 0)/np.sum(psf,0).max(), color='black', linewidth=3)
                    ax_vert.plot(np.sum(psf, 1)/np.sum(psf,1).max(),self.dy_subpixelsforintegrating_axis, color='black', linewidth=3)
                    ax_vert.semilogx()
                    ax_hori.semilogy()
                    ax_hori.set_ylim(1e-6,1e0)
                    ax_vert.set_xlim(1e-6,1e0)

                    ax_vert.tick_params(labelleft=False)
                    ax_hori.tick_params(labelbottom=False)
                    if title is not None:
                        ax_hori.set_title(title)
                    ax_map.imshow(np.sqrt(psf), cmap='gray_r', extent=[self.dx_subpixelsforintegrating_axis.min(), self.dx_subpixelsforintegrating_axis.max(), self.dy_subpixelsforintegrating_axis.min(), self.dy_subpixelsforintegrating_axis.max()],vmin=0.0,interpolation='nearest')
                    plt.draw()
                    if output is not None:
                        ax_map.figure.savefig(output)
                        self.speak("   saved PSF plot to {0}".format(output))


            # loop over focal plane radii
            self.speak('calculating wavelength-integrated PSFs at focalplaneradii of {0}'.format(self.focalplaneradii))
            for focalplaneradius in self.focalplaneradii:
                try:
                    self.temperaturelibrary
                except:
                    self.temperaturelibrary = {}

                # loop over stellar temperatures
                self.speak('calculating wavelength-integrated PSFs at temperatures of {0}'.format(self.temperatures))
                for temperature in self.temperatures:
                    try:
                        self.temperaturelibrary[focalplaneradius]
                    except:
                        self.temperaturelibrary[focalplaneradius] = {}

                    # create an empty psf image
                    psf = np.zeros_like(self.dx_subpixelsforintegrating)

                    # loop over wavelengths, adding in monochromatic psfs according to their weights
                    for w in np.arange(weighting.shape[1]):
                        # loop over powers in the polynomial for weight vs. teff
                        this_weight = 0.0

                        # loop over the terms of the temperature polynomial (this is stupidly slow, but it only has to happen once)
                        for p in np.arange(weighting.shape[0]):
                            this_weight += weighting[p,w]*(4000.0/temperature)**p
                            psf += this_weight*self.monochromaticlibrary[focalplaneradius][self.wavelengths[w]]
                            #self.speak("         {0} X [{1}]".format(this_weight, self.wavelengths[w]))

                    # after having calculated PSF image, store it in the library
                    self.temperaturelibrary[focalplaneradius][temperature] = psf/np.sum(psf)

                    # save a plot of the PSF
                    if plot:
                        filename = self.plotdirectory + "wavelengthintegrated_{0:05.0f}Kat{1:.0f}".format(temperature, focalplaneradius) + '.pdf'
                        plotPSF(psf, title="{0}K star at focalplaneradius (0,{1})".format(temperature, focalplaneradius),output= filename)

            np.save(temperature_filename, (self.temperaturelibrary, self.temperatures, self.focalplaneradii))
            self.speak('saved wavelength-integrated, high-resolution PSF library to {0}'.format(temperature_filename))
            self.plotTemperatureDependence()


    def plotTemperatureDependence(self):
        '''Plot the unjittered PSF library, after binning over wavelengths for stars of different temperatures.'''
        try:
            self.temperaturelibrary
        except:
            self.populateHighResolutionWavelengthIntegratedLibrary()

        # which radii and temperatures do we loop over?
        focalplaneradii = self.temperaturelibrary.keys()
        temperatures = self.temperaturelibrary[focalplaneradii[0]].keys()
        focalplaneradii.sort()
        temperatures.sort()

        # create a grid of plotting axes
        fi, ax = plt.subplots(len(focalplaneradii), len(temperatures), figsize=(len(temperatures)*4, len(focalplaneradii)*4), sharex=True, sharey=True)
        extent=[self.dx_subpixelsforintegrating_axis.min(), self.dx_subpixelsforintegrating_axis.max(), self.dy_subpixelsforintegrating_axis.min(), self.dy_subpixelsforintegrating_axis.max()]
        fi.suptitle("PSFs convolved with 2D jitter in a {0}s exposure.".format(self.camera.cadence))

        for p in range(len(focalplaneradii)):
            for t in range(len(temperatures)):

                psf = self.temperaturelibrary[focalplaneradii[p]][temperatures[t]]

                # plot this entry in the library
                ax[p,t].imshow(np.sqrt(psf), extent=extent, cmap='gray_r', interpolation='nearest', vmin=0 )

                # for the last row, include temperature xlabels
                if p == len(focalplaneradii)-1:
                    ax[p,t].set_xlabel("{0}K".format(temperatures[t]))

                # for the first column, include position labels
                if t == 0:
                        ax[p,t].set_ylabel("{0:4.1f} pix.".format(focalplaneradii[p]))

        # plot and save the figure
        plt.subplots_adjust(hspace=0, wspace=0)
        plt.savefig(self.plotdirectory + 'wavelengthintegratedpsfs{0:04.0f}_temperaturedependence.pdf'.format(self.camera.cadence))


    def populateHighResolutionJitteredLibrary(self):
        '''Convolve the high-resolution PSF with the jitter expected for the camera's cadence.'''
        self.speak("populating the high-resolution PSF library for an expected jitter for a {0:04.0f}s cadence.".format(self.camera.cadence))
        jittered_filename = self.versiondirectory + 'jittered{0:04.0f}s_psfs.npy'.format(self.camera.cadence)
        try:
            # maybe the jittered library can be loaded straightaway?
            self.jitteredlibrary, self.temperatures, self.focalplaneradii, self.jitteredlibrarytime = np.load(jittered_filename)
            self.speak('loaded jittered, high-resolution PSF library from {0}'.format(jittered_filename))

        except:
            try:
                self.temperaturelibrary
            except:
                self.populateHighResolutionWavelengthIntegratedLibrary()

            self.jitteredlibrary = copy.copy(self.temperaturelibrary)
            self.speak('calculating jittered PSFs at focalplaneradii of {0}'.format(self.focalplaneradii))
            for focalplaneradius in self.focalplaneradii:
                self.speak('calculating jittered PSFs at temperatures of {0}'.format(self.temperatures))
                for temperature in self.temperatures:
                    self.jitteredlibrary[focalplaneradius][temperature] = scipy.signal.convolve2d(self.temperaturelibrary[focalplaneradius][temperature], self.camera.jittermap[0]/np.sum(self.camera.jittermap[0]), 'same', 'fill', 0)
                    self.jitteredlibrarytime = self.camera.cadence
                    np.save(jittered_filename, (self.jitteredlibrary, self.temperatures, self.focalplaneradii, self.jitteredlibrarytime))
                    self.speak('saved jittered, high-resolution PSF library to {0}'.format(jittered_filename))
                    assert(self.jitteredlibrarytime == self.camera.cadence)


    def highResolutionPSF(self, position, temperature=5000, chatty=False):
        '''Paint a high resolution PSF at a given (cartographer-style) position, for a star of a given temperature.'''

        if chatty:
            self.speak("generating high-resolution PSF for {position} and T={temperature}K".format(position=position, temperature=temperature))
        # make sure a jittered library has been loaded
        try:
            self.jitteredlibrary
        except:
            self.populateHighResolutionJitteredLibrary()

        # get a cartographic object, so we can translate among coordinates
        focalx, focaly = position.focalxy.tuple
        focalplaneradius, focalplanetheta = position.focalrtheta.tuple

        # which available temperature is closest to ours?
        rounded_temperature = zachopy.utils.find_nearest(self.temperatures, temperature)
        if chatty:
            self.speak('  assuming {temperature} can be approximated as {rounded_temperature}'.format(**locals()))

        # which two available library focal plane radii are closest to ours?
        radiusbounds = zachopy.utils.find_two_nearest(self.focalplaneradii, focalplaneradius)
        radius_below, radius_above = radiusbounds
        fraction_below, fraction_above = zachopy.utils.interpolation_weights(radiusbounds, focalplaneradius)

        below = self.jitteredlibrary[radius_below]
        above = self.jitteredlibrary[radius_above]
        if chatty:
            self.speak("  treating focal plane radius of {focalplaneradius} as {fraction_below}x{radius_below} + {fraction_above}x{radius_above}".format(**locals()))
        unrotated_psf = below[rounded_temperature]*fraction_below + above[rounded_temperature]*fraction_above

        # rotate that PSF to the appropriate angle (IS THE ROTATION IN THE RIGHT DIRECTION???)
        rotated_psf = scipy.ndimage.interpolation.rotate(unrotated_psf, -focalplanetheta*180/np.pi, reshape=False)

        # return a PSF that has the right shape, but is still centered at (0.0, 0.0)
        return rotated_psf

    def binHighResolutionPSF(self, position, temperature=5000, dx=0.0, dy=0.0, plot=True, chatty=False):
        '''Pixelize a high resolution PSF at a given coordinate (cannot be used to put stars directly on images).'''

        # center the CCD subarray at the *rounded* pixel of the quoted position
        if chatty:
            self.speak('pixelizing a PSF at {0}'.format(position))
        self.ccd.center = np.array(np.round(position.focalxy.tuple))
        if chatty:
            self.speak('moved PSF subarray to {0}'.format(self.cartographer.ccd.center))

        # create the high-resolution PSF appropriate for this position
        zeroCenteredSubgridPSF = self.highResolutionPSF(position, temperature)

        # calculate the x and y subpixel grids, for both the unshifted and shifted pixels
        unshiftedx, unshiftedy = self.dx_subpixelsforintegrating, self.dy_subpixelsforintegrating
        shiftedx, shiftedy = unshiftedx + dx, unshiftedy + dy


        # MAKE SURE THERE'S NOT AN ACCIDENTAL TRANSPOSE IN HERE!
        prnu = self.intrapixel.prnu

        zeroCenteredBinnedPSF, xedges, yedges = \
            np.histogram2d(unshiftedy.flatten(), unshiftedx.flatten(), \
                            bins=[self.dy_pixels_edges, self.dx_pixels_edges], \
                            weights=zeroCenteredSubgridPSF.flatten()*prnu(unshiftedy.flatten(), unshiftedx.flatten()))
        recenteredBinnedPSF, xedges, yedges = \
            np.histogram2d(shiftedy.flatten(), shiftedx.flatten(), \
                            bins=[self.dy_pixels_edges, self.dx_pixels_edges], \
                            weights=zeroCenteredSubgridPSF.flatten()*prnu(shiftedy.flatten(), shiftedx.flatten()))

        if plot:

            plt.figure('Pixelizing the PSF', figsize=(7.5,7.5), dpi=100)
            plt.clf()
            gs = plt.matplotlib.gridspec.GridSpec(2,2,wspace=0.05,hspace=0.05, top=0.83)
            self.axHighResolution = plt.subplot(gs[0,0])
            self.axHighResolutionNudged = plt.subplot(gs[0,1], sharex=self.axHighResolution, sharey=self.axHighResolution)
            self.axPixelized = plt.subplot(gs[1,0], sharex=self.axHighResolution, sharey=self.axHighResolution)
            self.axPixelizedNudged = plt.subplot(gs[1,1], sharex=self.axHighResolution, sharey=self.axHighResolution)

            axes = [self.axHighResolution, self.axHighResolutionNudged, self.axPixelized, self.axPixelizedNudged]


            for a in axes:
                a.imshow(prnu(unshiftedx, unshiftedy), extent=extent(unshiftedx, unshiftedy), cmap='Oranges_r', alpha=0.25, interpolation='nearest')

            kw = dict( interpolation='nearest', cmap='gray_r', alpha=0.75)
            vmax=np.maximum(np.max(zeroCenteredBinnedPSF), np.max(recenteredBinnedPSF))
            self.axHighResolution.imshow(zeroCenteredSubgridPSF, extent=extent(unshiftedx, unshiftedy), **kw)
            self.axHighResolutionNudged.imshow(zeroCenteredSubgridPSF,  extent=extent(shiftedx, shiftedy),  **kw)
            vmax=np.maximum(np.max(zeroCenteredBinnedPSF), np.max(recenteredBinnedPSF))
            self.axPixelized.imshow(zeroCenteredBinnedPSF, extent=extent(self.dx_pixels_edges, self.dy_pixels_edges), vmax=vmax, **kw)
            self.axPixelizedNudged.imshow(recenteredBinnedPSF, extent=extent(self.dx_pixels_edges, self.dy_pixels_edges), vmax=vmax,**kw)
            pixeledgekw = dict(alpha=0.1, color='green')

            # draw the pixel boundries
            '''for a in axes:
                for x in self.dx_pixels_edges:
                    a.axvline(x, **pixeledgekw)
                for y in self.dy_pixels_edges:
                    a.axhline(y, **pixeledgekw)'''
            plotsize = 3.5
            self.axHighResolution.set_xlim(-plotsize,plotsize)
            self.axHighResolution.set_ylim(-plotsize,plotsize)
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
            plt.suptitle("Pixelizing the TESS Point Spread Function\n{temperature:.0f}K star, jittered for {jitter:.0f}s, intrapixel of {intrapixel}\n({focalx:.0f},{focaly:.0f}) pixels from focal plane center\n({dx:.2f},{dy:.2f}) from pixel center".format(focalx=self.ccd.center[0], focaly=self.ccd.center[1], dx=dx, dy=dy, temperature=temperature, jitter=self.jitteredlibrarytime, intrapixel=self.intrapixel.name))
            plt.draw()

        centralx, centraly = position.ccdxy.integerpixels

        # return the pixelized,
        return recenteredBinnedPSF, centralx + self.dx_pixels, centraly + self.dy_pixels


    # populate a library of binned PRFS, using the jittered high-resolution library
    def populateBinned(self, plot=False):
        '''Populate a library of binned PRFs, using the jittered, wavelength-integrated, high-resolution library.'''
        binned_filename = self.versiondirectory + 'pixelizedlibrary_{cadence:04.0f}s_{nradii:.0f}radii_{noffset:.0f}offsets_{drotation:.0f}degrees_{dtemperature:.0f}K_{intrapixel}.npy'.format(nradii=self.nradii, noffset=self.noffset, drotation=self.drotation, dtemperature=self.dtemperature, cadence=self.camera.cadence, intrapixel=self.intrapixel.name)
        try:
            (self.binned, self.focalr, self.temperatures, self.focaltheta, self.xoffsets, self.yoffsets) = np.load(binned_filename)
            self.speak('loaded binned PSFs from {0}'.format(binned_filename))
        except:
            self.speak('creating a new library of binned PSFs')
            self.populateHighResolutionJitteredLibrary()

            # set up the grid of values over which the
            self.focaltheta = np.arange(0, 360, self.drotation)
            self.xoffsets = np.linspace(-0.5,0.5,self.noffset)
            self.yoffsets = np.linspace(-0.5,0.5,self.noffset)
            self.focalr = np.linspace(0, np.max(self.focalplaneradii),self.nradii)
            try:
                self.binned
            except:
                self.binned = {}
                for radius in self.focalr:
                    self.speak('adding focal plane radius of {0:.1f} pixels'.format(radius), 1)
                    try:
                        self.binned[radius]
                    except:
                        self.binned[radius] = {}

                        for temperature in self.temperatures:
                            self.speak('adding temperature of {0:.0f}K'.format(temperature), 2)

                            try:
                                self.binned[radius][temperature]
                            except:
                                self.binned[radius][temperature] = {}

                                for theta in self.focaltheta:
                                    self.speak('adding focal plane theta of {0:.0f} degrees'.format(theta), 3)
                                    try:
                                        self.binned[radius][temperature][theta]
                                    except:
                                        self.binned[radius][temperature][theta] = {}
                                        position = self.cartographer.point(radius, theta/180.0*np.pi, 'focalrtheta')
                                        for xoffset in self.xoffsets:
                                            self.speak('adding xoffset of {0:.2f} pixels'.format(xoffset), 4)
                                            try:
                                                self.binned[radius][temperature][theta][xoffset]
                                            except:
                                                self.binned[radius][temperature][theta][xoffset] = {}

                                                for yoffset in self.yoffsets:
                                                    self.speak('adding yoffset of {0:.2f} pixels'.format(yoffset), 5)
                                                    self.binned[radius][temperature][theta][xoffset][yoffset] = self.binHighResolutionPSF(position, temperature=temperature, dx=xoffset, dy=yoffset, plot=plot)[0]
                                                    self.speak('(r={radius:.0f}, T={temperature:.0f}K, theta={theta:.0f} degrees , dx={xoffset:.2f} pixels, dy={yoffset:.2f} pixels)'.format(radius=radius, temperature=temperature, theta=theta, xoffset=xoffset, yoffset=yoffset), 6)
                np.save(binned_filename, (self.binned, self.focalr, self.temperatures, self.focaltheta, self.xoffsets, self.yoffsets))
                self.speak('saved binned PSF library to {0}'.format(binned_filename))



    def comparePSFs(self, position, temperature=4000, verbose=False, plot=True, justnew=False, center=None):
        '''Compare the PSF pulled out of the library to a newly pixelized one.'''

        newpsf, newx, newy = self.newlyPixelizedPSF(position, temperature)
        librarypsf, libraryx, libraryy = self.pixelizedPSF(position, temperature)

        if plot:
            # set up the plotting figure
            plt.figure('How Good is the PSF Library?', figsize=(7.5,7.5), dpi=100)
            plt.clf()
            gs = plt.matplotlib.gridspec.GridSpec(2,2,wspace=0.05,hspace=0.05, top=0.8, height_ratios=[1, 0.75])
            self.axLibrary = plt.subplot(gs[0,0])
            self.axNew = plt.subplot(gs[0,1], sharex=self.axLibrary, sharey=self.axLibrary)
            self.axLightcurve = plt.subplot(gs[1,:])
            axes = [self.axLibrary, self.axNew]

            #if center is None:
            #    center = position

            # plot the intrapixel sensitivity
            unshiftedx, unshiftedy = self.dx_subpixelsforintegrating, self.dy_subpixelsforintegrating
            integerx, integery = center.ccdxy.integerpixels
            for a in axes:
                a.imshow(self.intrapixel.prnu(unshiftedx, unshiftedy), extent=extent(unshiftedx+integerx, unshiftedy+integery), cmap='Oranges_r', alpha=0.25, interpolation='nearest')

            kw = dict( interpolation='nearest', cmap='gray_r', alpha=0.75, vmin=0, vmax=0.5)
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
            self.axLibrary.set_xlim(-plotsize+ccdx,plotsize+ccdx)
            self.axLibrary.set_ylim(-plotsize+ccdy,plotsize+ccdy)
            self.axLibrary.set_title('Interpolated from Library', fontsize=12)
            self.axNew.set_title('Just Pixelized', fontsize=12)
            for a in axes:
                plt.setp(a.get_xticklabels(), visible=False)
                plt.setp(a.get_yticklabels(), visible=False)
            dx, dy = position.ccdxy.fractionalpixels
            plt.suptitle("Comparing Recently Recalculated PSF to the Library\n{temperature:.0f}K star, jittered for {jitter:.0f}s, intrapixel of {intrapixel}\n({focalx:.0f},{focaly:.0f}) pixels from focal plane center\n({dx:.2f},{dy:.2f}) from pixel center".format(focalx=self.ccd.center[0], focaly=self.ccd.center[1], dx=dx, dy=dy, temperature=temperature, jitter=self.jitteredlibrarytime, intrapixel=self.intrapixel.name))
            #plt.draw()
        return np.sum(librarypsf), np.sum(newpsf)

    def newlyPixelizedPSF(self, position, temperature=4000, verbose=False):
        '''Wrapper for binHighResolutionPSF, to directly compare with pixelizedPSF.'''
        self.speak('re-pixelizing a PSF at {0}'.format(position))
        dx, dy = position.ccdxy.fractionalpixels
        self.speak('at subpixel offsets of {0}; CCD center at {1} with size of {2}'.format((dx,dy), self.ccd.center, self.ccd.npix))
        return self.binHighResolutionPSF(position, temperature, dx=dx, dy=dy, plot=False)

    def pixelizedPSF(self, position, temperature=4000, verbose=False):
        '''Drop a pixelized PSF, drawn from the library, at a particular position.'''

        # make sure the binned PSF library is already loaded
        try:
            self.binned
        except:
            self.populateBinned()

        # need to determine [radius][temperature][theta][xoffset][yoffset] to pull out of library
        focalr, focaltheta = position.focalrtheta.tuple
        key_radius = zachopy.utils.find_nearest(self.focalr, focalr, verbose=verbose)
        key_theta = zachopy.utils.find_nearest(self.focaltheta, focaltheta*180/np.pi, verbose=verbose)
        key_temperature = zachopy.utils.find_nearest(self.temperatures, temperature, verbose=verbose)

        ccdx, ccdy = position.ccdxy.tuple
        xoffset, yoffset = position.ccdxy.fractionalpixels
        centralx, centraly = position.ccdxy.integerpixels


        xoffsetbounds = zachopy.utils.find_two_nearest(self.xoffsets, xoffset, verbose=verbose)
        xbelow, xabove = xoffsetbounds
        xbelow_weight, xabove_weight = zachopy.utils.interpolation_weights(xoffsetbounds, xoffset)
        yoffsetbounds = zachopy.utils.find_two_nearest(self.yoffsets, yoffset, verbose=verbose)
        ybelow, yabove = yoffsetbounds
        ybelow_weight, yabove_weight = zachopy.utils.interpolation_weights(yoffsetbounds, yoffset)

        interpolatedPSF = \
            xbelow_weight*ybelow_weight*self.binned[key_radius][key_temperature][key_theta][xbelow][ybelow] + \
            xabove_weight*ybelow_weight*self.binned[key_radius][key_temperature][key_theta][xabove][ybelow] + \
            xabove_weight*yabove_weight*self.binned[key_radius][key_temperature][key_theta][xabove][yabove] + \
            xbelow_weight*yabove_weight*self.binned[key_radius][key_temperature][key_theta][xbelow][yabove]

        assert( interpolatedPSF is not None)
        return interpolatedPSF, centralx + self.dx_pixels, centraly + self.dy_pixels

def extent(x,y):
    return [x.min(), x.max(), y.min(), y.max()]
