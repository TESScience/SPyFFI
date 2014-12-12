from imports import *
import pixels

# define everything related to PSFs
class PSF:

    # initialize the PSF class
    def __init__(self, camera=None):
        print "Initializing the TESS point spread function painter."

        # details related to the orginal PSF files from Deb, by way of Peter.
        self.initial_binning = 4                            # how much to bin the original (very high resolution, from Deb) PSFs before integrating over the pixels
        self.subpixelsforintegrating_size = 0.25/15.0*self.initial_binning    # with this binning, how big is each subpixelsforintegratingel in the high resolution PSF?
        self.numberof_subpixelsforintegrating = 960/self.initial_binning             # how many subpixelsforintegratingels are there in the high resolution PSF image?

        # parameters used for generating the binned PSF library (the "stamper" will use this library to choose how to paint each star on the image)
        self.dtemperature = 2000    # the step in temperature, stamper will choose nearest
        self.noffset = 21           # number of sub-pixel offsets, stamper will interpolate between offsets
        self.drotation = 10         # number of rotation angles for PSF's (0-360), stamper will ???
        self.nradii = 10            # number of radial focalplanetradius for PSF's, stamper will ???

        # define grids of subpixel offsets (positive and negative, centered on 0.0, 0.0), over which the PSF integration will be carried out
        # (how far are the *centers* of the subpixels from the star's focalplanetradius?)
        self.dx_subpixelsforintegrating_axis = np.arange(self.numberof_subpixelsforintegrating)*self.subpixelsforintegrating_size
        self.dx_subpixelsforintegrating_axis -= np.mean(self.dx_subpixelsforintegrating_axis)
        self.dy_subpixelsforintegrating_axis = self.dx_subpixelsforintegrating_axis
        self.dx_subpixelsforintegrating, self.dy_subpixelsforintegrating = np.meshgrid(self.dx_subpixelsforintegrating_axis, self.dy_subpixelsforintegrating_axis)

        if False:   # (do I actually use any of these?)
            # create a grid of full pixels that would contain at least some subpixels, centered on 0.0, 0.0
            left, right = self.subpixel2pixel((np.min(self.dx_subpixelsforintegrating), np.max(self.dx_subpixelsforintegrating)))
            bottom, top = self.subpixel2pixel((np.min(self.dy_subpixelsforintegrating), np.max(self.dy_subpixelsforintegrating)))
            self.dx_pixels_axis = np.arange(left, right+1)
            self.dx_pixels_axis = np.arange(right, top+1)
            self.dx_pixels, self.dy_pixels = np.meshgrid(self.dx_pixels_axis, self.dx_pixels_axis)

        # fill in an intrapixel sensitivity
        # self.intrapixel = pixels.prnu(self.dx_subpixelsforintegrating, self.dx_subpixelsforintegrating)

        # if possible, associate a Camera object with this PSF
        self.setCamera(camera)
        #self.populateBinned()
        self.populateHeader()

    def subpixel2pixel(self,anysubpixel):
        '''Convert from fractional pixel to an integer pixel (pixel centers are at 0.0's; edges are at 0.5's).'''
        return np.round(anysubpixel).astype(np.int)

    def populateHeader(self):
        '''Create a PSF header structure, to store the details of the PSF simulation.'''
        self.header = astropy.io.fits.Header()
        self.header[''] = ('', '')
        self.header['PRF'] = ''
        self.header['PRFNOTE'] = ('','Pixel Response Function library parameters')
        self.header['PLIBDTEM'] = (self.dtemperature, '[K] d(effective temperature)')
        self.header['PLIBNOFF'] = (self.noffset, '# of subpixelsforintegratingel shifts, in both x and y')
        self.header['PLIBDROT'] = (self.drotation, '[deg] d(rotation around center)')
        self.header['PLIBNRAD'] = (self.nradii, '# of radial distances from field center')
        self.header['PSUBSIZE'] = (self.subpixelsforintegrating_size, '[pix] subpixel grid size used in initial PSF integration')
        #self.header['PNSUBPIX'] = (self.numberof_subpixelsforintegrating, '[pix] # of subpixels used for initial PSF integration')
        #self.header['PPIXSIZE'] = (self.pixsize, '[pix] pixel size')

    def populateHighResolutionMonochromaticLibrary(self, plot=True, inspect=False):
        '''Load the high resolution PSFs from Deb, which are monochromatic and unjittered.'''
        print " -populating the raw, monochromatic, PSF library"
        monochromatic_filename = settings.prefix + 'intermediates/monochromatic_psf_library.npy'
        try:
            # maybe we can load from a preprocessed file?
            (self.monochromaticlibrary, self.focalplaneradii, self.wavelengths, self.dx, self.dy, self.dx_subpixelsforintegrating_axis, self.dy_subpixelsforintegrating_axis)  = np.load(monochromatic_filename)
        except:
            # load the high resolution PSF's from a FITs file
            data = astropy.io.fits.open(settings.prefix + 'inputs/psfs_3p33.fits')[0].data

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

            # input the focalplanetradius and wavelength coordinates that Deb used
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
                            ax_map.imshow(np.log(self.monochromaticlibrary[focalplaneradii[i]][wavelengths[j]]),extent=extent, cmap='gray_r', interpolation='nearest', vmin=0  )
                            plt.draw()

                        # populate a plot in the library grid
                        ax_psf_fits[i,j].imshow(np.log(self.monochromaticlibrary[focalplaneradii[i]][wavelengths[j]]), extent=extent, cmap='gray_r', interpolation='nearest', vmin=0 )
                        if i == data.shape[0]-1:
                            ax_psf_fits[i,j].set_xlabel("{0} micron".format(wavelengths[j]))

                        if j == 0:
                            ax_psf_fits[i,j].set_ylabel("{0:4.1f} pix.".format(focalplaneradii[i]))
                        plt.subplots_adjust(hspace=0, wspace=0)
            if plot:
                plt.draw()
                plt.savefig(settings.prefix + 'plots/monochromatic_psf_library.pdf')

            # save the library, so it'd be easier to load next time
            np.save(monochromatic_filename,(self.monochromaticlibrary, self.focalplaneradii, self.wavelengths, self.dx, self.dy, self.dx_subpixelsforintegrating_axis, self.dy_subpixelsforintegrating_axis))

    def populateHighResolutionTemperatureLibrary(self, plot=True):
        '''Integrate the monochromatic PSFs into a (still unjittered) wavelength-integrated PSF library.'''

        print " -populating the wavelength-integrated (unjittered) PSF library."
        temperature_filename = settings.prefix + 'intermediates/temperature_psfs.npy'.format(self.camera.cadence)
        try:
            # maybe we can load a preprocessed file?
            self.temperaturelibrary, self.temperatures, self.focalplaneradii = np.load(temperature_filename)
        except:
            # otherwise, we have to recalculate from scratch

            try:
                # make sure the monochromatic library has been populated
                self.monochromaticlibrary
            except:
                self.populateHighResolutionMonochromaticLibrary()

            # load weight coefficients calculated by Peter
            weighting = astropy.io.fits.open(settings.prefix + 'inputs/ph_filt.fits')[0].data

            # create a grid of temperatures at which PSFs will be calculated
            self.temperatures = np.arange(3000, 12001, self.dtemperature)

            # loop over focal plane radii
            for focalplanetradius in self.focalplaneradii:
                try:
                    self.temperaturelibrary
                except:
                    self.temperaturelibrary = {}

                    # loop over stellar temperatures
                    for temperature in self.temperatures:
                        try:
                            self.temperaturelibrary[focalplanetradius]
                        except:
                            self.temperaturelibrary[focalplanetradius] = {}

                        # create an empty psf image
                        psf = np.zeros_like(self.dx)

                        # loop over wavelengths, adding in monochromatic psfs according to their weights
                        for w in np.arange(weighting.shape[1]):
                            # loop over powers in the polynomial for weight vs. teff
                            this_weight = 0.0
                            for p in np.arange(weighting.shape[0]):
                                this_weight += weighting[p,w]*(4000.0/temperature)**p
                                psf += this_weight*self.monochromaticlibrary[focalplanetradius][self.wavelengths[w]]
                                print "         {0} X [{1}]".format(this_weight, self.wavelengths[w])

                                filename = "integrated_psf_for_{0}Kat{1}".format(temperature, focalplanetradius*100).replace(' ', '').replace('.','p') + '.pdf'
                                self.temperaturelibrary[focalplanetradius][temperature] = psf/np.sum(psf)
                                if plot:
                                    self.plot(psf, title="{0}K star at focalplanetradius (0,{1})".format(temperature, focalplanetradius),output=settings.prefix + 'plots/' + filename)
                                    print ''
                                    np.save(temperature_filename, (self.temperaturelibrary, self.temperatures, self.focalplaneradii))

                                    # convolve the high resolution PSFs with the jitter expected for the camera's cadence
                                    def populateJitteredLibrary(self):
                                        print "Populating the high-resolution PSF library for an expected jitter for a {0:04.0f}s cadence.".format(self.camera.cadence)
                                        jittered_filename = settings.prefix + 'intermediates/jittered_psfs_{0:04.0f}s.npy'.format(self.camera.cadence)
                                        try:
                                            self.jitteredlibrary, self.temperatures, self.focalplaneradii, self.jitteredlibrarytime = np.load(jittered_filename)
                                        except:
                                            try:
                                                self.temperaturelibrary
                                            except:
                                                self.populateHighResolutionTemperatureLibrary()

                                                self.jitteredlibrary = copy.copy(self.temperaturelibrary)
                                                for focalplanetradius in self.focalplaneradii:
                                                    for temperature in self.temperatures:
                                                        self.jitteredlibrary[focalplanetradius][temperature] = scipy.signal.convolve2d(self.temperaturelibrary[focalplanetradius][temperature], self.camera.jittermap[0]/np.sum(self.camera.jittermap[0]), 'same', 'fill', 0)
                                                        self.jitteredlibrarytime = self.camera.cadence
                                                        np.save(jittered_filename, (self.jitteredlibrary, self.temperatures, self.focalplaneradii, self.jitteredlibrarytime))
                                                        assert(self.jitteredlibrarytime == self.camera.cadence)

    # populate a library of binned PRFS, using the jittered high-resolution library
    def populateBinned(self):
        binned_filename = settings.prefix + 'intermediates/binned_psf_library_{0:04.0f}s.npy'.format(self.camera.cadence)
        try:
            print "Trying to load ", binned_filename
            (self.binned, self.radii, self.temperatures, self.rotations, self.xoffsets, self.yoffsets) = np.load(binned_filename)
        except:
            print "Didn't find {0}, so creating a new binned library."
            self.populateHighResolutionMonochromaticLibrary()
            self.populateHighResolutionTemperatureLibrary()
            self.rotations = np.arange(0, 360, self.drotation)
            self.xoffsets = np.linspace(0,1,self.noffset)
            self.yoffsets = np.linspace(0,1,self.noffset)
            self.radii = np.linspace(0, np.max(self.focalplaneradii),self.nradii)
            try:
                self.binned
            except:
                self.binned = {}
                for focalplanetradius in self.radii:
                    try:
                        self.binned[focalplanetradius]
                    except:
                        self.binned[focalplanetradius] = {}

                        for temperature in self.temperatures:
                            try:
                                self.binned[focalplanetradius][temperature]
                            except:
                                self.binned[focalplanetradius][temperature] = {}

                                for rotation in self.rotations:
                                    try:
                                        self.binned[focalplanetradius][temperature][rotation]
                                    except:
                                        self.binned[focalplanetradius][temperature][rotation] = {}

                                        for xoffset in self.xoffsets:
                                            try:
                                                self.binned[focalplanetradius][temperature][rotation][xoffset]
                                            except:
                                                self.binned[focalplanetradius][temperature][rotation][xoffset] = {}

                                                for yoffset in self.yoffsets:
                                                    radius = focalplanetradius
                                                    # make sure the sign is right!
                                                    x = np.round(radius*np.cos(-rotation*np.pi/180)) + xoffset
                                                    y = np.round(radius*np.sin(-rotation*np.pi/180)) + yoffset
                                                    self.binned[focalplanetradius][temperature][rotation][xoffset][yoffset] = self.binnedhighrespsf( x, y, temperature)
                                                    #self.plot(subgrid_psf*intrapixel,title='({0},{1},{2},{3},{4})'.format(focalplanetradius, temperature, rotation, xoffset, yoffset))
                                                    print '({0},{1},{2},{3},{4})'.format(focalplanetradius, temperature, rotation, xoffset, yoffset)
                                                    np.save(binned_filename, (self.binned, self.radii, self.temperatures, self.rotations, self.xoffsets, self.yoffsets))

    def binnedhighrespsf(self, x, y, temperature):
        subgrid_psf = self.highrespsf(x, y, temperature)
        assert(subgrid_psf.shape == self.intrapixel.shape)
        thisbinned = zachopy.utils.rebin_total(subgrid_psf*self.intrapixel,2*self.gridpixels,2*self.gridpixels)
        return thisbinned

    def binnedpsf(self, x, y, temperature, verbose=False, interpolation='linear'):

        try:
            self.binned
        except:
            self.populateBinned()
            r = np.sqrt(x**2 + y**2)
            focalplanetradius = zachopy.utils.find_nearest(self.radii, r, verbose=verbose)
            temperature = zachopy.utils.find_nearest(self.temperatures, temperature, verbose=verbose)
            rot = self.rotation(x,y)*180/np.pi
            if verbose:
                print ""
                print "{0},{1},{2},{3}".format(x, y, temperature, rot)
                while rot > 360 or rot < 0:
                    if rot < 0:
                        rot += 360
                        if rot > 360:
                            rot -= 360

                            rotation = zachopy.utils.find_nearest(self.rotations, rot , verbose=verbose)
                            if interpolation == 'nearest':
                                xoffset = zachopy.utils.find_nearest(self.xoffsets, x % 1, verbose=verbose)
                                yoffset = zachopy.utils.find_nearest(self.yoffsets, y % 1, verbose=verbose)
                                return self.binned[focalplanetradius][temperature][rotation][xoffset][yoffset]
                                if interpolation == 'linear':
                                    xfrac, yfrac = x % 1, y % 1
                                    xbelow, xabove = zachopy.utils.find_two_nearest(self.xoffsets, xfrac, verbose=verbose)
                                    xspan = xabove - xbelow
                                    xbelow_weight, xabove_weight = (xabove - xfrac)/xspan, (xfrac - xbelow)/xspan
                                    assert(xbelow_weight > 0)
                                    if verbose:
                                        print "  weights are ", xbelow_weight, xabove_weight
                                        ybelow, yabove = zachopy.utils.find_two_nearest(self.yoffsets, yfrac, verbose=verbose)
                                        yspan = yabove - ybelow
                                        ybelow_weight, yabove_weight =  (yabove - yfrac)/yspan, (yfrac - ybelow)/yspan
                                        if verbose:
                                            print "  weights are ", ybelow_weight, yabove_weight
                                            assert(xbelow_weight + xabove_weight ==1)
                                            assert(ybelow_weight + yabove_weight ==1)

                                            return 	xbelow_weight*ybelow_weight*self.binned[focalplanetradius][temperature][rotation][xbelow][ybelow] + \
                                            xabove_weight*ybelow_weight*self.binned[focalplanetradius][temperature][rotation][xabove][ybelow] + \
                                            xabove_weight*yabove_weight*self.binned[focalplanetradius][temperature][rotation][xabove][yabove] + \
                                            xbelow_weight*yabove_weight*self.binned[focalplanetradius][temperature][rotation][xbelow][yabove]


    def setCamera(self, camera):
        '''Associated a camera structure with this PSF painter.'''
        self.camera = camera

    def plotTemperatures(self):
        '''Plot the unjittered PSF library, after binning over wavelengths for stars of different temperatures.'''
        try:
            self.temperaturelibrary
        except:
            self.populateHighResolutionTemperatureLibrary()

        # which radii and temperatures do we loop over?
        focalplaneradii = self.temperaturelibrary.keys()
        temperatures = self.temperaturelibrary[focalplaneradii[0]].keys()
        focalplaneradii.sort()
        temperatures.sort()

        # create a grid of plotting axes
        fi, ax = plt.subplots(len(focalplaneradii), len(temperatures), figsize=(len(temperatures)*4, len(focalplaneradii)*4), sharex=True, sharey=True)
        extent=[self.dx_subpixelsforintegrating_axis.min(), self.dx_subpixelsforintegrating_axis.max(), self.dy_subpixelsforintegrating_axis.min(), self.dy_subpixelsforintegrating_axis.max()]
        fi.suptitle("PSFs convolved with 2D jitter in a {0}s efocalxure.".format(self.camera.cadence))

        for p in range(len(focalplaneradii)):
            for t in range(len(temperatures)):

                # plot this entry in the library
                ax[p,t].imshow(np.log(psf), extent=extent, cmap='gray_r', interpolation='nearest', vmin=0 )

                # for the last row, include temperature xlabels
                if p == len(focalplaneradii)-1:
                    ax[p,t].set_xlabel("{0}K".format(temperatures[t]))

                # for the first column, include position labels
                if t == 0:
                        ax[p,t].set_ylabel("{0:4.1f} pix.".format(focalplaneradii[p]))

        # plot and save the figure
        plt.draw()
        plt.savefig(settings.prefix + 'plots/psfs_with_temperature_{0:04.0f}.pdf'.format(self.camera.cadence))


    def plot(self, psf, title=None, output=None):
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
        ax_map.imshow(np.log(psf), cmap='gray_r', extent=[self.dx_subpixelsforintegrating_axis.min(), self.dx_subpixelsforintegrating_axis.max(), self.dy_subpixelsforintegrating_axis.min(), self.dy_subpixelsforintegrating_axis.max()],vmin=0.0,interpolation='nearest')
        plt.draw()
        if output is not None:
            ax_map.figure.savefig(output)
            print "   saved figure to: ", output

    def rotation(self, focalx, focaly):
        '''Convert an (x,y) position (in focal plane coordinates) into an angle relative to the field center.'''
        return -np.arctan2(focaly, focalx)

    def highrespsf(self, focalx, focaly, temperature, verbose=False):
        '''Paint a high resolution PSF at a given x and y position (in focal plane coordinates), for a star of a given temperature.'''

        if verbose:
            print "  -generating a high resolution PSF for focal plane positions ({x}, {y}) and T={temperature}K".format(x=focalx, y=focaly, temperature=temperature)

        try:
            # make sure a jittered library has been loaded
            self.jitteredlibrary
        except:
            self.populateJitteredLibrary()

        # which available temperature is closest to ours?
        rounded_temperature = zachopy.utils.find_nearest(self.temperatures, temperature)
        if verbose:
            "     treating stellar temperature of {temperature} as {rounded_temperature}".format(**locals())

        # which two available library focal plane radii are closest to ours?
        focalplaneradius = np.sqrt(focalx**2 + focaly**2)
        bounds = zachopy.utils.find_two_nearest(self.focalplaneradii, focalplaneradius)
        span = bounds[1] - bounds[0]
        fraction_above = (focalplaneradius - bounds[0])/span
        fraction_below = 1.0 - fraction_above
        radius_above = self.focalplaneradii[bounds[1]]
        radius_below = self.focalplaneradii[bounds[0]]
        assert((fraction_below >= 0)&(fraction_above >=0))
        if verbose:
            "     treating focal plane radius of {focalplaneradius} as {fraction_below}x{radius_below} + {fraction_above}x{radius_above}".format(**locals())
        below = self.jitteredlibrary[radius_below]
        above = self.jitteredlibrary[radius_above]

        # intepolate to get a PSF at the right radius
        unrotated_psf = below[rounded_temperature]*fraction_below + above[rounded_temperature]*frac

        # rotate that PSF to the appropriate angle
        rotated_psf = scipy.ndimage.interpolation.rotate(unrotated_psf, self.rotation(focalx, focaly)*180/np.pi, reshape=False)

        # nudge the PSF from its (0.0, 0.0) position relative
        nudged_psf = np.zeros_like(rotated_psf)
        xnudge = np.int((focalx % 1)/self.subpixelsforintegrating_size)
        ynudge = np.int((focaly % 1)/self.subpixelsforintegrating_size)
        #print xnudge, ynudge, ' nudges'
        nudged_psf[xnudge:, ynudge:]= rotated_psf[:rotated_psf.shape[0]-xnudge, :rotated_psf.shape[1]-ynudge]
        #print [xnudge, ynudge]
        #print [rotated_psf.shape[0]-xnudge, rotated_psf.shape[1]-ynudge]
        return nudged_psf

    def xy2radiusrotation(self, focalx, focaly, rotationoffset=0.0):
        '''Convert from (x,y) focal plane positions to (radius,rotation) focal plane positions.'''
        radius = np.sqrt(x**2 + y**2)
        rotation = np.arctan2(focaly, focalx) + rotationoffset
        return radius, rotation

    def radiusrotation2xy(self, radius, rotation, rotationoffset=0.0):
        '''Convert from (radius,rotation) focal plane positions to (x,y) focal plane positions.'''

        focalx = np.cos(rotation - rotationoffset)*radius
        focalx = np.sin(rotation - rotationoffset)*radius
        return focalx, focaly
