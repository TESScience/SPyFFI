# define everything related to PSFs
class PSF:

  # initialize the PSF class
  def __init__(self, camera=None):

    '''Details related to the orginal PSF files from Deb, by way of Peter.'''
    # how much to bin the original (very high resolution, from Deb) PSFs before integrating over the pixels
    self.initial_binning = 4
    # with this binning, how big is each subpixel in the high resolution PSF?
    self.subpixsize = 0.25/15.0*self.initial_binning
    # how many subpixels are there in the high resolution PSF image?
    self.nsubpix = 960/self.initial_binning

    # how big is a pixel, in pixel units (perhaps silly, but might come in handy)
    self.pixsize = 1.0

    # parameters used for generating the binned PSF library
    #  (the "stamper" will use this library to choose how to paint each star on the image)
    self.dtemperature = 2000    # the step in temperature, stamper will choose nearest
    self.noffset = 21           # number of sub-pixel offsets, stamper will interpolate between offsets
    self.drotation = 10         # number of rotation angles for PSF's (0-360), stamper will ???
    self.nradii = 10            # number of radial position for PSF's, stamper will ???

    self.dxax = np.arange(self.nsubpix)*self.subpixsize
    self.dxax -= np.mean(self.dxax)
    self.dyax = self.dxax
    self.dx, self.dy = np.meshgrid(self.dxax, self.dyax)
    self.gridpixels = self.dx.shape[0]*self.subpixsize/self.pixsize/2
    self.subgridresolution=self.subpixsize/self.pixsize
    self.subgrid = np.arange(-self.gridpixels,self.gridpixels , self.subgridresolution)
    self.xsubgrid, self.ysubgrid = np.meshgrid(self.subgrid, self.subgrid)
    self.grid = np.arange(-self.gridpixels,self.gridpixels,1,dtype=np.int)
    self.xgrid, self.ygrid = np.meshgrid(self.grid, self.grid)
    self.intrapixel = prnu(self.dx, self.dy)

    self.setCamera(camera)
    #self.populateBinned()
    self.populateHeader()

  def populateHeader(self):
    self.header = astropy.io.fits.Header()
    self.header[''] = ('', '')
    self.header['PRF'] = ''
    self.header['PRFNOTE'] = ('','Pixel Response Function library parameters')
    self.header['PLIBDTEM'] = (self.dtemperature, '[K] d(effective temperature)')
    self.header['PLIBNOFF'] = (self.noffset, '# of subpixel shifts, in both x and y')
    self.header['PLIBDROT'] = (self.drotation, '[deg] d(rotation around center)')
    self.header['PLIBNRAD'] = (self.nradii, '# of radial distances from field center')
    self.header['PNSUBPIX'] = (self.nsubpix, '[pix] # of initial subpixels')
    #self.header['PSUBSIZE'] = (self.subpixsize, '[pix] initial subpixel grid-size')
    #self.header['PPIXSIZE'] = (self.pixsize, '[pix] pixel size')


  # load the raw high resolution PSFs from Deb (no jitter yet)
  def populateMonochromaticLibrary(self, plot=True, inspect=False):
    print "Populating the raw, monochromatic, PSF library."
    try:
      (self.monochromaticlibrary, self.positions, self.wavelengths, self.dx, self.dy, self.dxax, self.dyax)  = np.load(settings.prefix + 'intermediates/psf_library.npy')
    except:
      # load the high resolution PSF's from a FITs file
      data = astropy.io.fits.open(settings.prefix + 'inputs/psfs_3p33.fits')[0].data
      # trim the high resolution PSF, so that an even number of whole pixels is included
      assert(data.shape[2] == data.shape[3])

      binned = zachopy.utils.rebin(data,data.shape[0], data.shape[1], data.shape[2]/self.initial_binning, data.shape[3]/self.initial_binning)
      ntrim = (binned.shape[2] - self.nsubpix)/2
      data = binned[:,:,ntrim:-ntrim,ntrim:-ntrim]
      npositions, nwavelengths, nxs, nys = data.shape
      self.monochromaticlibrary = {}

      positions = np.array([0,0.5,1.0,np.sqrt(2)])*2048
      wavelengths = 0.675 + np.arange(nwavelengths)*0.05
      self.positions = positions
      self.wavelengths = wavelengths
      assert(nxs == nys)

      extent=[self.dxax.min(), self.dxax.max(), self.dyax.min(), self.dyax.max()]
      if plot:
        fi, ax_psf_fits  = plt.subplots(data.shape[0],data.shape[1], figsize=(25.95, 9.8), sharex=True, sharey=True)

      for i in range(data.shape[0]):

        self.monochromaticlibrary[positions[i]] = {}
        for j in range(data.shape[1]):

          #interpolator = scipy.interpolate.interp2d(rho.flatten()[::10], phi.flatten()[::10], data[i,j,:,:].flatten()[::10], kind='linear', fill_value=0.0)
          #assert(False)
          self.monochromaticlibrary[positions[i]][wavelengths[j]] = data[i,j,:,:]
          if plot:
            if inspect:
              try:
                for a in ax_psf_zoom:
                  a.cla()
              except:
                fig = plt.figure(figsize=(10,10))
                plt.subplots_adjust(hspace=0, wspace=0)
              ax_map = fig.add_subplot(2,2,3)
              ax_vert = fig.add_subplot(2,2,4, sharey=ax_map)
              ax_hori = fig.add_subplot(2,2,1, sharex=ax_map)
              ax_psf_zoom = [ax_map, ax_vert, ax_hori]
              ax_hori.plot(self.dxax, np.sum(self.monochromaticlibrary[positions[i]][wavelengths[j]], 0))
              ax_vert.plot(np.sum(self.monochromaticlibrary[positions[i]][wavelengths[j]], 1), self.dyax)
              ax_map.imshow(np.log(self.monochromaticlibrary[positions[i]][wavelengths[j]]),extent=extent, cmap='gray_r', interpolation='nearest', vmin=0  )
              plt.draw()
            ax_psf_fits[i,j].imshow(np.log(self.monochromaticlibrary[positions[i]][wavelengths[j]]), extent=extent, cmap='gray_r', interpolation='nearest', vmin=0 )
            if i == data.shape[0]-1:
              ax_psf_fits[i,j].set_xlabel("{0} micron".format(wavelengths[j]))

            if j == 0:
              ax_psf_fits[i,j].set_ylabel("{0:4.1f} pix.".format(positions[i]))
            plt.subplots_adjust(hspace=0, wspace=0)
      if plot:
        plt.draw()
        plt.savefig(settings.prefix + 'plots/psf_library.pdf')


      np.save(settings.prefix + 'intermediates/psf_library.npy',(self.monochromaticlibrary, self.positions, self.wavelengths, self.dx, self.dy, self.dxax, self.dyax))

  # integrate the monochromatic PSFs into a temperature library (no jitter yet)
  def populateTemperatureLibrary(self, plot=True):
    print "Populating the wavelength-integrated (unjittered) PSF library."
    temperature_filename = settings.prefix + 'intermediates/temperature_psfs.npy'.format(self.camera.cadence)
    try:
      self.temperaturelibrary, self.temperatures, self.positions = np.load(temperature_filename)
    except:
      try:
        self.monochromaticlibrary
      except:
        self.populateMonochromaticLibrary()
      weighting = astropy.io.fits.open(settings.prefix + 'inputs/ph_filt.fits')[0].data
      self.temperatures = np.arange(3000, 12001, self.dtemperature)
      for position in self.positions:
        try:
          self.temperaturelibrary
        except:
          self.temperaturelibrary = {}

        # loop over stellar temperatures
        for temperature in self.temperatures:
          try:
            self.temperaturelibrary[position]
          except:
            self.temperaturelibrary[position] = {}
          print "   modeling a PSF at {0}".format( position)
          print "      with a temperature of {0}K".format(temperature)
          psf = np.zeros_like(self.dx)
          # loop over wavelengths
          for w in np.arange(weighting.shape[1]):
            # loop over powers in the polynomial for weight vs. teff
            this_weight = 0.0
            for p in np.arange(weighting.shape[0]):
              this_weight += weighting[p,w]*(4000.0/temperature)**p
            psf += this_weight*self.monochromaticlibrary[position][self.wavelengths[w]]
            print "         {0} X [{1}]".format(this_weight, self.wavelengths[w])

          filename = "integrated_psf_for_{0}Kat{1}".format(temperature, position*100).replace(' ', '').replace('.','p') + '.pdf'
          self.temperaturelibrary[position][temperature] = psf/np.sum(psf)
          if plot:
            self.plot(psf, title="{0}K star at position (0,{1})".format(temperature, position),output=settings.prefix + 'plots/' + filename)
          print ''
      np.save(temperature_filename, (self.temperaturelibrary, self.temperatures, self.positions))

  # convolve the high resolution PSFs with the jitter expected for the camera's cadence
  def populateJitteredLibrary(self):
    print "Populating the high-resolution PSF library for an expected jitter for a {0:04.0f}s cadence.".format(self.camera.cadence)
    jittered_filename = settings.prefix + 'intermediates/jittered_psfs_{0:04.0f}s.npy'.format(self.camera.cadence)
    try:
      self.jitteredlibrary, self.temperatures, self.positions, self.jitteredlibrarytime = np.load(jittered_filename)
    except:
      try:
        self.temperaturelibrary
      except:
        self.populateTemperatureLibrary()

      self.jitteredlibrary = copy.copy(self.temperaturelibrary)
      for position in self.positions:
        for temperature in self.temperatures:
          self.jitteredlibrary[position][temperature] = scipy.signal.convolve2d(self.temperaturelibrary[position][temperature], self.camera.jittermap[0]/np.sum(self.camera.jittermap[0]), 'same', 'fill', 0)
      self.jitteredlibrarytime = self.camera.cadence
      np.save(jittered_filename, (self.jitteredlibrary, self.temperatures, self.positions, self.jitteredlibrarytime))
    assert(self.jitteredlibrarytime == self.camera.cadence)

  # populate a library of binned PRFS, using the jittered high-resolution library
  def populateBinned(self):
    binned_filename = settings.prefix + 'intermediates/binned_psf_library_{0:04.0f}s.npy'.format(self.camera.cadence)
    try:
      print "Trying to load ", binned_filename
      (self.binned, self.radii, self.temperatures, self.rotations, self.xoffsets, self.yoffsets) = np.load(binned_filename)
    except:
      print "Didn't find {0}, so creating a new binned library."
      self.populateMonochromaticLibrary()
      self.populateTemperatureLibrary()
      self.rotations = np.arange(0, 360, self.drotation)
      self.xoffsets = np.linspace(0,1,self.noffset)
      self.yoffsets = np.linspace(0,1,self.noffset)
      self.radii = np.linspace(0, np.max(self.positions),self.nradii)
      try:
        self.binned
      except:
        self.binned = {}
      for position in self.radii:
        try:
          self.binned[position]
        except:
          self.binned[position] = {}

        for temperature in self.temperatures:
          try:
            self.binned[position][temperature]
          except:
            self.binned[position][temperature] = {}

          for rotation in self.rotations:
            try:
              self.binned[position][temperature][rotation]
            except:
              self.binned[position][temperature][rotation] = {}

            for xoffset in self.xoffsets:
              try:
                self.binned[position][temperature][rotation][xoffset]
              except:
                self.binned[position][temperature][rotation][xoffset] = {}

              for yoffset in self.yoffsets:
                radius = position
                # make sure the sign is right!
                x = np.round(radius*np.cos(-rotation*np.pi/180)) + xoffset
                y = np.round(radius*np.sin(-rotation*np.pi/180)) + yoffset
                self.binned[position][temperature][rotation][xoffset][yoffset] = self.binnedhighrespsf( x, y, temperature)
                #self.plot(subgrid_psf*intrapixel,title='({0},{1},{2},{3},{4})'.format(position, temperature, rotation, xoffset, yoffset))
            print '({0},{1},{2},{3},{4})'.format(position, temperature, rotation, xoffset, yoffset)
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
    position = zachopy.utils.find_nearest(self.radii, r, verbose=verbose)
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
      return self.binned[position][temperature][rotation][xoffset][yoffset]
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

      return 	xbelow_weight*ybelow_weight*self.binned[position][temperature][rotation][xbelow][ybelow] + \
          xabove_weight*ybelow_weight*self.binned[position][temperature][rotation][xabove][ybelow] + \
          xabove_weight*yabove_weight*self.binned[position][temperature][rotation][xabove][yabove] + \
          xbelow_weight*yabove_weight*self.binned[position][temperature][rotation][xbelow][yabove]


  def setCamera(self, camera):
    self.camera = camera

  def plotTemperatures(self):
    try:
      self.temperaturelibrary
    except:
      self.populateTemperatureLibrary()

    positions = self.temperaturelibrary.keys()
    temperatures = self.temperaturelibrary[positions[0]].keys()
    positions.sort()
    temperatures.sort()
    fi, ax = plt.subplots(len(positions), len(temperatures), figsize=(len(temperatures)*4, len(positions)*4), sharex=True, sharey=True)
    extent=[self.dxax.min(), self.dxax.max(), self.dyax.min(), self.dyax.max()]
    fi.suptitle("PSFs convolved with 2D jitter in a {0}s exposure.".format(self.camera.cadence))
    for p in range(len(positions)):
      for t in range(len(temperatures)):
        ax[p,t].imshow(np.log(psf), extent=extent, cmap='gray_r', interpolation='nearest', vmin=0 )

        if p == len(positions)-1:
          ax[p,t].set_xlabel("{0}K".format(temperatures[t]))

        if t == 0:
          ax[p,t].set_ylabel("{0:4.1f} pix.".format(positions[p]))
    plt.show()
    plt.savefig(settings.prefix + 'plots/psfs_with_temperature_{0:04.0f}.pdf'.format(self.camera.cadence))


  def plot(self, psf, title=None, output=None):
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
    ax_hori.plot(self.dxax, np.sum(psf, 0)/np.sum(psf,0).max(), color='black', linewidth=3)
    ax_vert.plot(np.sum(psf, 1)/np.sum(psf,1).max(),self.dyax, color='black', linewidth=3)
    ax_vert.semilogx()
    ax_hori.semilogy()
    ax_hori.set_ylim(1e-6,1e0)
    ax_vert.set_xlim(1e-6,1e0)

    ax_vert.tick_params(labelleft=False)
    ax_hori.tick_params(labelbottom=False)
    if title is not None:
      ax_hori.set_title(title)
    ax_map.imshow(np.log(psf), cmap='gray_r', extent=[self.dxax.min(), self.dxax.max(), self.dyax.min(), self.dyax.max()],vmin=0.0,interpolation='nearest')
    plt.draw()
    if output is not None:
      ax_map.figure.savefig(output)
      print "   saved figure to: ", output

  def rotation(self, xpos, ypos):
    return -np.arctan2(ypos, xpos)

  def highrespsf(self, xpos, ypos, temperature):

    try:
      self.jitteredlibrary
    except:
      self.populateJitteredLibrary()

    possible_positions = self.jitteredlibrary.keys()
    r = np.sqrt(xpos**2 + ypos**2)
    indices = np.arange(len(self.positions))
    i_fraction = np.interp(r, self.positions, indices)
    below = self.jitteredlibrary[self.positions[np.floor(i_fraction)]]
    above = self.jitteredlibrary[self.positions[np.ceil(i_fraction)]]
    frac = i_fraction % 1.0
    #print "   {0}".format( frac)
    rounded_temperature = zachopy.utils.find_nearest(self.temperatures, temperature)
    unrotated_psf = below[rounded_temperature]*(1.0 - frac) + above[rounded_temperature]*frac
    rotated_psf = scipy.ndimage.interpolation.rotate(unrotated_psf, self.rotation(xpos, ypos)*180/np.pi, reshape=False)
    nudged_psf = np.zeros_like(rotated_psf)
    xnudge = np.int((xpos % 1)/self.subpixsize)
    ynudge = np.int((ypos % 1)/self.subpixsize)
    #print xnudge, ynudge, ' nudges'
    nudged_psf[xnudge:, ynudge:]= rotated_psf[:rotated_psf.shape[0]-xnudge, :rotated_psf.shape[1]-ynudge]
    #print [xnudge, ynudge]
    #print [rotated_psf.shape[0]-xnudge, rotated_psf.shape[1]-ynudge]
    return nudged_psf

  def xy2rhophi(self, x, y, rotation=0.0):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x) + rotation
    return rho, phi

  def rhophi2xy(self, rho, phi, rotation=0.0):
    x = np.cos(phi - rotation)*rho
    y = np.sin(phi - rotation)*rho
    return x,y
