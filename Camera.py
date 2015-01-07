from imports import *
import settings, catalogs
from PSF import PSF
from Cartographer import Cartographer
from CCD import CCD

# define a camera class
class Camera(Talker):
    '''Keep track of one camera's entire field of view.'''
    def __init__(self, cadence=1800, ra=270,dec=66.56070833333332, testpattern=False, subarray=None):
        '''Initialize camera, fill it with CCDs, and point it at the sky or at a testpattern.'''

        # decide whether or not this Camera is chatty
        Talker.__init__(self, mute=False, pithy=False)

        # [self.speak will only report for this object if mute and pithy are turned off]
        self.speak("Turning on a new TESS camera object.")

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
        if subarray is None:
            # if we're not dealing with a subarray, then turn on CCD's 1,2,3,4
            self.ccdnumbers = np.arange(4) + 1
        else:
            # if this is a subarray, then turn on one (imaginary) CCD and call it 0
            self.ccdnumbers = np.arange(1)
        self.ccds = [CCD(n,subarray=subarray,camera=self) for n in self.ccdnumbers]

        # point the Camera either at real stars (from the sky) or at a test pattern (a grid of stars)
        self.testpattern = testpattern
        if self.testpattern:
            # pretend a test pattern is pointed at (ra, dec) = (0.0, 0.0)
            self.ra = 0.0
            self.dec = 0.0
        else:
            # if real stars, use the input (ra, dec)
            self.ra = ra
            self.dec = dec
        self.point(self.ra, self.dec)

        # assign a cartographer to this Camera and start it out on the first CCD
        self.cartographer = Cartographer(camera=self, ccd=self.ccds[0])

        # start the camera out unjittered from its nominal position
        self.nudge = {'x':0.0, 'y':0.0, 'z':0.0}						# nudge relative to nominal spacecraft pointing (arcsec)

        # keep track of which exposure is being simulated
        self.counter = 0

        # load the PSF for this Camera
        self.psf = PSF(camera=self)


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

        # the pixel scale, in degrees, taken to be FOV/number of pixels
        w.wcs.cdelt = [-self.fov/self.npix,self.fov/self.npix]

        # the celestial coordinates at the reference position (input by user)
        nudged_ra, nudged_dec = zachopy.spherical.rotate(self.ra, self.dec,  self.nudge['x']/60.0/60.0, self.nudge['y']/60.0/60.0)
        print "!!!!!!!!!!!!!!!"
        w.wcs.crval = [nudged_ra, nudged_dec]

        # the rotation of the field (currently just a pure rotation, no shear)
        #rot = self.nudge['z']/60.0/60.0*np.pi/180.0
        #w.wcs.pc = [[np.cos(rot), -np.sin(rot)],[np.sin(rot), np.cos(rot)]]

        # the coordinate system type - what should I use?
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

        # set this to be the WCS
        self.wcs = w
        self.populateHeader()

  def loadJitterBall(self):
    print "Loading jitterball for {0} second cadences.".format(self.cadence)
    try:
      # check to see if the jitterball is already loaded
      self.jitterball
      # make sure the we're using the right jitterball for this cadence
      assert(self.jittercadence == self.cadence)
    except:
      # if the jitterball isn't already loaded, then do so
      jitterfile = settings.prefix + 'intermediates/jitter_{0:04.0f}.npy'.format(self.cadence)
      try:
        # if a processed jitterfile already exists, then just load it up:
        self.jitterball, self.jittermap = np.load(jitterfile)
        self.jittercadence = self.cadence
      except:
        print "   taking the raw jitter input file AttErrTimeArcsec_80k.dat"
        # otherwise, create a binned jitterfile
        # load simulated jitterball that Roland got from Orbital
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
        print "   smoothing the jitter to {0}s cadence.".format(self.cadence)
        spacing = data['t'][1] - data['t'][0]
        n = np.long(self.cadence/spacing)
        filter = np.ones(n)/n
        smoothed_t = np.convolve(data['t'], filter, mode='valid')
        smoothed_x = np.convolve(data['x'], filter, mode='valid')
        smoothed_y = np.convolve(data['y'], filter, mode='valid')
        smoothed_z = np.convolve(data['z'], filter, mode='valid')


        t = np.convolve(data['t'], filter, mode='valid')[::n]
        x = np.convolve(data['x'], filter, mode='valid')[::n]
        y = np.convolve(data['y'], filter, mode='valid')[::n]
        z = np.convolve(data['z'], filter, mode='valid')[::n]

        # plot each dimension separately
        print "   plotting the binned jitter timeseries to " + settings.prefix + 'plots/jitter_timeseries_{0:04.0f}.pdf'.format(self.cadence)
        fi, ax = plt.subplots(3,1, sharey=True, sharex=True)
        plt.ion()
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
        fi.savefig(settings.prefix + 'plots/jitter_timeseries_{0:04.0f}.pdf'.format(self.cadence))

        narcsec = 3
        bins = np.ceil(narcsec/self.pixelscale/self.psf.subpixsize).astype(np.int)*2 +1	# bins are in units of subpixels
        range = [[-(bins-1)/2*self.psf.subpixsize,(bins-1)/2*self.psf.subpixsize],[-(bins-1)/2*self.psf.subpixsize,(bins-1)/2*self.psf.subpixsize]] # range is in units of pixels

        x_interpolator = scipy.interpolate.interp1d(smoothed_t, smoothed_x,'nearest',fill_value=0,bounds_error=False)
        y_interpolator = scipy.interpolate.interp1d(smoothed_t, smoothed_y,'nearest',fill_value=0,bounds_error=False)
        jittermap_to_plot = np.histogram2d((data['x'] - x_interpolator(data['t']))/self.pixelscale, (data['y']  - y_interpolator(data['t']))/self.pixelscale, bins=50, range=range,normed=True)

        print "   plotting the binned jitter map to " + settings.prefix + 'plots/jitter_map_{0:04.0f}.pdf'.format(self.cadence)

        # plot map
        '''fi, ax = plt.subplots(1,1, sharey=True, sharex=True)
        plt.ion()
        ax.imshow(self.jittermap[0] , cmap='gray_r', extent=[-3,3,-3,3],interpolation='nearest')

        ax.set_xlim(-3,3)
        ax.set_ylim(-3,3)
        ax.set_title('TESS Pointing Jitter over {0}s'.format(self.cadence))
        ax.set_ylabel('x (")')
        ax.set_ylabel('y (")')
        plt.show()'''

        self.plothist2d(jittermap_to_plot,title='TESS Pointing Jitter over {0}s'.format(self.cadence),xtitle='Pixels', ytitle='Pixels')
        plt.savefig(settings.prefix + 'plots/jitter_map_{0:04.0f}.pdf'.format(self.cadence))




        print range
        print bins
        self.jittermap = np.histogram2d((data['x'] - x_interpolator(data['t']))/self.pixelscale, (data['y']  - y_interpolator(data['t']))/self.pixelscale, bins=bins, range=range,normed=True)
        self.plothist2d(self.jittermap ,title='TESS Pointing Jitter over {0}s'.format(self.cadence),xtitle='Pixels', ytitle='Pixels')
        plt.savefig(settings.prefix + 'plots/jitter_map_adopted_{0:04.0f}.pdf'.format(self.cadence))

        self.jitterball = (x,y,z)
        self.jittercadence = self.cadence
        np.save(jitterfile,( self.jitterball, self.jittermap))

  def plothist2d(self, hist, title=None, log=False, xtitle=None, ytitle=None):
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

  def jitter(self, header=None):
    try:
      self.jitterball
    except:
      self.loadJitterBall()
    print "Jittering TESS a little bit."
    self.nudge['x'] = self.counter*0.0#self.jitterball[0][self.counter]#*0.002
    self.nudge['y'] = -self.counter*3.0#self.jitterball[1][self.counter]#*0.002
    self.nudge['z'] = 0.0#self.jitterball[2][self.counter]#*0.002

    try:
      header['MOTION'] = ''
      header['MOTNOTE'] = ('', 'properties of the image motion applied')
      header['JITTERX'] = (self.nudge['x'], '["] jitter-induced nudge')
      header['JITTERY'] = (self.nudge['y'], '["] jitter-induced nudge')
    except:
      pass

    self.point()

  def pos_string(self):
    if self.testpattern:
        return "testpattern"
    else:
        coords = astropy.coordinates.ICRS(ra=self.ra*astropy.units.degree, dec=self.dec*astropy.units.degree)
        return "{0:02}h{1:02}m{2:02}s{3:+03}d{4:02}m{5:02}s".format(np.int(coords.ra.hms[0]),np.int(coords.ra.hms[1]),np.int(coords.ra.hms[2].round()), np.int(coords.dec.dms[0]),np.int(np.abs(coords.dec.dms[1])),np.int(np.abs(coords.dec.dms[2].round())))

  def project(self, plot=False, write=False):
    print "Creating a starmap for this TESS field."
    #ras, decs, rmag, jmag, imag, temperatures = catalogs.stars(ra=self.ra, dec=self.dec, radius=self.fov/np.sqrt(2)*1.01, catalog='UCAC4')
    ras, decs, imag, temperatures = self.catalog.arrays()
    deltamag = np.max(imag) - imag
    size = deltamag**2/16.0

    #try:
    #	self.ax_radec.cla()
    #except:
    #	fi, self.ax_radec = plt.subplots(1,1)
    #self.ax_radec.scatter(ras, decs, s=size, marker='.', color='black', alpha=0.3)


    #focalxy = self.cartographer.focalxy(ras, decs, 'radec')
    #xcenter, ycenter = self.wcs.crpix
    #x, y = focalxy.x + xcenter, focalxy.y + ycenter
    x, y = self.wcs.wcs_world2pix(ras, decs, 1)
    starsfilename = settings.prefix + 'plots/' +  "starfield_{ra}_{dec}".format(ra=self.ra, dec=self.dec).replace(' ', '') + '.pdf'
    if not os.path.exists(starsfilename):
      if plot:
        try:
          self.ax_projected.cla()
        except:
          fi, self.ax_projected = plt.subplots(1,1)
        self.ax_projected.scatter(x, y, s=size, marker='o', color='black', alpha=0.3, edgecolors='none')
        self.ax_projected.plot([0, 4096, 4096, 0, 0], [0, 0, 4096, 4096, 0], color='red', alpha=0.5, linewidth=3)
        self.ax_projected.set_aspect(1)
        self.ax_projected.set_xlim(-128, 4096+128)
        self.ax_projected.set_ylim(-128, 4096+128)
        self.ax_projected.set_xlabel('Pixel Coordinates')
        self.ax_projected.set_ylabel('Pixel Coordinates')
        self.ax_projected.figure.savefig(starsfilename)

    if write:
      outfile = settings.prefix + 'intermediates/stars_{pos}_{fov:2.0f}deg.txt'.format(pos=self.pos_string(), fov=self.fov)
      np.savetxt(outfile, np.c_[ras, decs, x, y, imag], fmt=['%.6f', '%.6f', '%.3f', '%.3f', '%.3f'])
      print " saved stars to {0}".format(outfile)


    self.starx = x
    self.stary = y
    self.starmag = np.array(imag)
    self.startemp = np.array(temperatures)

  def photons(self, mag):
    # this doesn't work for M dwarfs, need to include multiband information
    return self.effective_area*10**(-0.4*mag)*1.7e6
