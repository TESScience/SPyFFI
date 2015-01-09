from imports import *
import settings, catalogs
from PSF import PSF
from Cartographer import Cartographer
from CCD import CCD

# define a camera class
class Camera(Talker):
    '''Keep track of one camera's entire field of view.'''
    def __init__(self, cadence=1800, ra=270,dec=66.56070833333332, testpattern=False, subarray=None, label='', number=1):
        '''Initialize camera, fill it with CCDs, and point it at the sky or at a testpattern.'''

        # decide whether or not this Camera is chatty
        Talker.__init__(self, mute=False, pithy=False)

        # keep track of which exposure is being simulated
        self.counter = 0

        # [self.speak will only report for this object if mute and pithy are turned off]
        self.speak("Turning on a new TESS camera object.")

        # keep track of any special features with this Camera
        self.subarray = subarray
        self.label = label
        self.testpattern = testpattern

        # decide whether to point the Camera either at real stars (from the sky) or at a test pattern (a grid of stars)
        if self.testpattern:
            # pretend a test pattern is pointed at (ra, dec) = (0.0, 0.0)
            self.ra = 0.0
            self.dec = 0.0
        else:
            # if real stars, use the input (ra, dec)
            self.ra = ra
            self.dec = dec

        # assign the cadence for this camera
        self.singleread = 2.0											# seconds
        self.readouttime = 0.005										# seconds
        self.setCadence(cadence)                                        # seconds

        # define scales for the Camera
        self.pixelscale = 21.1#24.0/4096*60*60							# arcsec!! (from Peter's paper)
        self.entrance_pupil_diameter = 10.5								# cm (from Peter's paper)
        self.effective_area = 69.1										# cm^2 (from Peter's paper)
        self.physicalpixelsize = 15.0/1e4								# cm
        self.physicalpixeldepth = 100.0/1e4								# cm
        self.read_noise = 10.0											# electrons per read
        self.saturation = 150000.0										# electrons per read

        # define the gaps between the CCDs (numbers still in flux, as of late 2014)
        self.physicalgap = 0.2											# cm
        self.gapinpixels = self.physicalgap/self.physicalpixelsize  	# pixels (not necessarily integer)

        # define the field of view of the Camera (one side of the CCD square)
        self.fov = (4096 + self.gapinpixels)*self.pixelscale/60.0/60.0					# degrees
        self.pixel_solid_area = (self.pixelscale)**2 		             # arcsec^2

        # turn on the necessary CCDs in this Camera
        if self.subarray is None:
            # if we're not dealing with a subarray, then turn on CCD's 1,2,3,4
            self.ccdnumbers = np.arange(4) + 1
        else:
            # if this is a subarray, then turn on one (imaginary) CCD and call it 0
            self.ccdnumbers = np.arange(1)
        self.ccds = [CCD(n,subarray=self.subarray,camera=self) for n in self.ccdnumbers]

        # assign a cartographer to this Camera and start it out on the first CCD
        self.cartographer = Cartographer(camera=self, ccd=self.ccds[0])



        # start the camera out unjittered from its nominal position
        self.nudge = {'x':0.0, 'y':0.0, 'z':0.0}						# nudge relative to nominal spacecraft pointing (arcsec)

        # point the Camera
        self.point(self.ra, self.dec)


        # load the PSF for this Camera
        self.psf = PSF(camera=self)

    @property
    def fielddirectory(self):
        if self.label == '':
            d = settings.prefix + 'outputs/{pos}/'.format(pos=self.pos_string())
        else:
            d = settings.prefix + 'outputs/{pos}_{label}/'.format(pos=self.pos_string(), label=self.label)
        zachopy.utils.mkdir(d)
        return d

    @property
    def directory(self):
        d = self.fielddirectory + '{0:d}s/'.format(self.cadence)
        zachopy.utils.mkdir(d)
        return d

    def expose(self):
        '''Take an exposure on all the available CCD's.'''
        for c in self.ccds:
            c.expose()

    def populateHeader(self):
        '''Populate the header structure with information about the Camera, and its WCS.'''

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
        if self.subarray is not None:
            self.header['SUBARRAY'] = (self.subarray, 'THIS IMAGE IS JUST {0}x{1} POSTAGE STAMP!'.format(self.subarray, self.subarray))

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
        '''Set the cadence of this Camera, in seconds.'''

        # set the cadence
        self.cadence = cadence
        self.speak("Setting cadence to {0} seconds = {1} reads.".format(self.cadence, self.cadence/self.singleread))

        # update the jitterball to one that has been binned to this cadence
        self.loadJitterBall()

    def point(self, ra=None, dec=None):
        '''Point this Camera at the sky, by using the field-specified (ra,dec) and (if active) the jitter nudge for this exposure.'''

        # if (ra,dec) given, then point the whole Camera there
        if ra is not None and dec is not None:
            self.ra, self.dec = ra, dec

        # make sure an (ra,dec) are defined
        try:
            self.ra
            self.dec
            self.speak('Pointing the camera at (ra,dec) = {0:.6f},{1:.6f}'.format(self.ra, self.dec))
        except:
            self.report("Please point your telescope somewhere. No RA or DEC defined.")

        # create a blank WCS object, for converting between (ra,dec) and (x,y)
        self.wcs = astropy.wcs.WCS(naxis=2)

        # the focalxy coordinates of the reference position [taken to be center of field, which is (x,y) = (0.0, 0.0) in focalxy coordinates]
        self.wcs.wcs.crpix = [0.0, 0.0]

        # the pixel scale, in degrees
        self.wcs.wcs.cdelt = [-self.pixelscale/60.0/60.0,self.pixelscale/60.0/60.0]

        # the celestial coordinates at the reference position (input by user)
        nudged_ra, nudged_dec = zachopy.spherical.rotate(self.ra, self.dec,  self.nudge['x']/60.0/60.0, self.nudge['y']/60.0/60.0)
        self.wcs.wcs.crval = [nudged_ra, nudged_dec]

        # the rotation of the field (currently just a pure rotation, no shear)
        #rot = self.nudge['z']/60.0/60.0*np.pi/180.0
        #w.wcs.pc = [[np.cos(rot), -np.sin(rot)],[np.sin(rot), np.cos(rot)]]

        # the coordinate system type - what should I use?
        self.wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

        # set this to be the WCS
        self.populateHeader()

    def loadJitterBall(self):
        '''Load the jitterball for this camera, binned to the appropriate cadence.'''
        self.speak('Populating the jitterball for {0:.0f} second cadence.'.format(self.cadence))


        try:
            # if the jitterball is already loaded and of the correct cadence, we're all set!
            self.jitterball

            # make sure the we're using the right jitterball for this cadence
            assert(self.jittercadence == self.cadence)
        except:
            # if the jitterball isn't already loaded (or has wrong cadence), then load/create a new one!
            jitterfile = settings.prefix + 'intermediates/jitter_{0:04.0f}.npy'.format(self.cadence)
            try:
                # if a processed jitterfile already exists, then just load it up:
                self.jitterball, self.jittermap = np.load(jitterfile)
                self.jittercadence = self.cadence
            except:
                # otherwise, create a binned jitterfile

                # load simulated jitterball that Roland got from Orbital
                self.speak("Using the raw jitter input file AttErrTimeArcsec_80k.dat (which for convenience may be saved in raw_jitter_file.npy)")
                try:
                  data = np.load(settings.prefix + 'intermediates/raw_jitter_file.npy')
                except:
                  data = astropy.io.ascii.read(settings.prefix + "inputs/AttErrTimeArcsec_80k.dat", names=['t','x', 'y', 'z'])
                  np.save(settings.prefix + 'intermediates/raw_jitter_file.npy', data)

                # subtract means
                data['x'] -= np.mean(data['x'])
                data['y'] -= np.mean(data['y'])
                data['z'] -= np.mean(data['z'])

                # scale jitterball to requirements (should be inflation by ~1.5)
                required_jitter_rms = 2.0/3
                original_fifthsec = np.sqrt(np.mean(data['x']**2 + data['y']**2))
                data['x'] *= required_jitter_rms/original_fifthsec
                data['y'] *= required_jitter_rms/original_fifthsec
                data['z'] *= required_jitter_rms/original_fifthsec

                # smooth them to the required cadence
                self.speak("Smoothing the jitter to {0}s cadence.".format(self.cadence))

                # define convolution filter to smooth timeseries over as many samples as necessary to match the cadence
                spacing = data['t'][1] - data['t'][0]
                n = np.long(self.cadence/spacing)
                filter = np.ones(n)/n

                # construct smoothed timeseries, sampled at native (high) time resolution
                smoothed_t = np.convolve(data['t'], filter, mode='valid')
                smoothed_x = np.convolve(data['x'], filter, mode='valid')
                smoothed_y = np.convolve(data['y'], filter, mode='valid')
                smoothed_z = np.convolve(data['z'], filter, mode='valid')

                # sample smoothed timeseries at the camera's cadence
                t = smoothed_t[::n]
                x = smoothed_x[::n]
                y = smoothed_y[::n]
                z = smoothed_z[::n]

                # plot each dimension separately
                plotfilename =  settings.prefix + 'plots/jitter_timeseries_{0:04.0f}.pdf'.format(self.cadence)
                self.speak("Saving plot of the binned jitter timeseries to " + plotfilename)

                # create the plot
                fi, ax = plt.subplots(3,1, sharey=True, sharex=True)
                ax[0].plot(data['t'], data['x'], alpha=0.5, color='black')
                ax[0].plot(t, x, linewidth=2, alpha=0.5, marker='o', color='red')
                ax[1].plot(data['t'], data['y'], alpha=0.5, color='black')
                ax[1].plot(t, y, linewidth=2, alpha=0.5, marker='o', color='red')
                ax[2].plot(data['t'], data['z'], alpha=0.5, color='black')
                ax[2].plot(t, z, linewidth=2, alpha=0.5, marker='o', color='red')
                ax[0].set_xlim(0,self.cadence*10)
                ax[0].set_title('Expected TESS Pointing Jitter for {0}s Cadence'.format(self.cadence))
                ax[0].set_ylabel('x (")')
                ax[1].set_ylabel('y (")')
                ax[2].set_ylabel('z (")')
                ax[2].set_xlabel('Time (seconds)')
                plt.show()
                fi.savefig(plotfilename)

                # create a 2D jittermap
                narcsec = 3
                bins = np.ceil(narcsec/self.pixelscale/self.psf.subpixsize).astype(np.int)*2 +1	# bins are in units of subpixels
                range = [[-(bins-1)/2*self.psf.subpixsize,(bins-1)/2*self.psf.subpixsize],[-(bins-1)/2*self.psf.subpixsize,(bins-1)/2*self.psf.subpixsize]] # range is in units of pixels

                # make interpolators to keep track of the running smooth means at every moment
                x_interpolator = scipy.interpolate.interp1d(smoothed_t, smoothed_x,'nearest',fill_value=0,bounds_error=False)
                y_interpolator = scipy.interpolate.interp1d(smoothed_t, smoothed_y,'nearest',fill_value=0,bounds_error=False)


                # assign the jittermap here, which will be used for convolution in the PSF code
                self.jittermap = np.histogram2d((data['x'] - x_interpolator(data['t']))/self.pixelscale, (data['y']  - y_interpolator(data['t']))/self.pixelscale, bins=bins, range=range,normed=True)
                self.jitterball = (x,y,z)
                self.jittercadence = self.cadence

                self.speak('Saving plots of the intra-exposure jitter map (which which PSFs should be convolved).')
                # plot an easier to view histogram of the jitterball
                jittermap_to_plot = np.histogram2d((data['x'] - x_interpolator(data['t']))/self.pixelscale, (data['y']  - y_interpolator(data['t']))/self.pixelscale, bins=50, range=range,normed=True)
                self.plothist2d(jittermap_to_plot,title='TESS Pointing Jitter over {0}s'.format(self.cadence),xtitle='Pixels', ytitle='Pixels')
                plt.savefig(settings.prefix + 'plots/jitter_map_{0:04.0f}.pdf'.format(self.cadence))

                # plot the adopted jitterball, as more useful binning
                self.plothist2d(self.jittermap ,title='TESS Pointing Jitter over {0}s'.format(self.cadence),xtitle='Pixels', ytitle='Pixels')
                plt.savefig(settings.prefix + 'plots/jitter_map_adopted_{0:04.0f}.pdf'.format(self.cadence))

                # save the necessary jitter files so we don't have to go through this again
                self.speak('Saving the jitter files for this cadence to {0}'.format(jitterfile))
                np.save(jitterfile,( self.jitterball, self.jittermap))

    def plothist2d(self, hist, title=None, log=False, xtitle=None, ytitle=None):
        '''Plot a 2D histogram. (Used for loadJitterBall plots -- should probably move all this to new Jitterball object.)'''
        map = hist[0]
        x = hist[1][1:] + (hist[1][0] - hist[1][1])/2.0
        y = hist[2][1:]+ (hist[2][0] - hist[2][1])/2.0
        fig = plt.figure(figsize=(10,10))
        plt.clf()
        plt.subplots_adjust(hspace=0, wspace=0)
        ax_map = fig.add_subplot(2,2,3)
        ax_vert = fig.add_subplot(2,2,4, sharey=ax_map)
        ax_hori = fig.add_subplot(2,2,1, sharex=ax_map)

        ax_hori.plot(x, np.sum(map, 0)/np.sum(map), marker='o', color='black', linewidth=3)
        ax_vert.plot(np.sum(map, 1)/np.sum(map), y, marker='o', color='black', linewidth=3)
        if log:
          ax_vert.semilogx()
          ax_hori.semilogy()
        if log:
          bottom = np.min(map[map > 0])/np.maximum(np.sum(map,0).max(),np.sum(map,1).max())
        else:
          bottom = 0
        top = 1
        ax_hori.set_ylim(bottom,top)
        ax_vert.set_xlim(bottom,top)

        ax_vert.tick_params(labelleft=False)
        ax_hori.tick_params(labelbottom=False)
        if title is not None:
          ax_hori.set_title(title)
        if xtitle is not None:
          ax_map.set_xlabel(xtitle)
        if ytitle is not None:
          ax_map.set_ylabel(ytitle)
        if log:
          ax_map.imshow(np.log(map), cmap='gray_r', extent=[x.min(), x.max(), y.min(), y.max()],interpolation='nearest')
        else:
          ax_map.imshow(map, cmap='gray_r', extent=[x.min(), x.max(), y.min(), y.max()],interpolation='nearest')
        plt.draw()

    def jitter(self, dx=None, dy=None, dz=None, header=None ):
        '''Jitter the cameras by a little bit, by introducing nudges draw from a (cadence-appropriate) jitterball timeseries.'''

        # make sure the jitterball has been populated
        self.loadJitterBall()

        # assign the nudges in two translations and one rotation
        if dx is None:
            self.nudge['x'] = self.jitterball[0][self.counter]#*0.002#self.counter*0.0#
        else:
            self.nudge['x'] = dx

        if dy is None:
            self.nudge['y'] = self.jitterball[1][self.counter]#*0.002#-self.counter*3.0#
        else:
            self.nudge['y'] = dy

        if dz is None:
            self.nudge['z'] = self.jitterball[2][self.counter]#*0.002#-self.counter*3.0#
        else:
            self.nudge['z'] = dz

        # if possible, write the details to the supplied FITS header
        try:
            header['MOTION'] = ''
            header['MOTNOTE'] = ('', 'properties of the image motion applied')
            header['JITTERX'] = (self.nudge['x'], '["] jitter-induced nudge')
            header['JITTERY'] = (self.nudge['y'], '["] jitter-induced nudge')
        except:
            pass

        # move the camera, using the updated nudge values
        self.speak("Jittering the camera to {x},{y},{z} away from nominal pointing.".format(**self.nudge))
        self.point()

    def pos_string(self):
        '''Return the position string for this field.'''
        if self.testpattern:
            return "testpattern"
        else:
            coords = astropy.coordinates.ICRS(ra=self.ra*astropy.units.degree, dec=self.dec*astropy.units.degree)
            return "{0:02}h{1:02}m{2:02}s{3:+03}d{4:02}m{5:02}s".format(np.int(coords.ra.hms[0]),np.int(coords.ra.hms[1]),np.int(coords.ra.hms[2].round()), np.int(coords.dec.dms[0]),np.int(np.abs(coords.dec.dms[1])),np.int(np.abs(coords.dec.dms[2].round())))

    def advanceCounter(self):
        '''Take one step forward in time with this Camera.'''
        self.camera.counter +=1

    def populateCatalog(self):
        '''Create a catalog of stars that are visible with this Camera.'''

        # figure out how wide we need to search for stars
        if self.subarray is None:
            size = self.fov
        else:
            size = self.subarray*self.pixelscale/3600.0

        # make a test pattern or a catalog of real stars
        if self.testpattern:
            # figure out how big a test pattern to create (in arcsec)

            self.catalog = catalogs.TestPattern(size=size*3600.0)
        else:
            self.catalog = catalogs.UCAC4(ra=self.ra, dec=self.dec, radius=size/np.sqrt(2)*1.01)



	def c1(self, image):

		return image[self.xsize - self.subarraysize:, self.ysize - self.subarraysize:]

	def c2(self, image):
		return image[0:self.subarraysize, self.ysize - self.subarraysize:]

	def c3(self, image):
		return image[0:self.subarraysize,0:self.subarraysize]

	def c4(self, image):
		return image[self.xsize - self.subarraysize:,0:self.subarraysize]
