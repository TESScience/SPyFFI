'''Keep track of spacecraft jitter.'''
from imports import *
import settings

class Jitter(Talker):
  def __init__(self, camera=None, jitterrms=None, rawjitterbasename="AttErrTimeArcsec_80k.dat", nsubpixelsperpixel=None,
      amplifyinterexposurejitter=1.0):
      # decide whether or not this CCD is chatty
      Talker.__init__(self, mute=False, pithy=False)

      # set an extra directory specifically for this object
      self.directory = 'jitter/'

      # store the input camera
      self.camera = camera

      # set up the initial raw jitter file
      # (this one cam from Roland, some time ago)
      self.rawfile = settings.prefix + "inputs/" + rawjitterbasename

      # what do you want the RMS to be rescaled to?
      self.jitterrms = jitterrms

      # how much should the exposure to exposure jitter be amplified
      self.amplifyinterexposurejitter = amplifyinterexposurejitter

      # (for creating map that will be used to convolve with the psf)
      self.nsubpixelsperpixel = nsubpixelsperpixel

      # update the jitterball to one that has been binned to this cadence
      self.load()


  def load(self, remake=False):
    '''make sure that a jitterball (timeseries of roll,pitch,yaw) has
      been loaded and binned to the appropriate exposure times'''

    try:
      # if the jitterball is already loaded into memory
      #  *and* of the correct cadence, we're all set!
      self.jitterball

      # make sure the we're using the right jitterball for this cadence
      assert(self.jittercadence == self.camera.cadence)

      # make sure we're not trying to remake the jitterball
      assert(remake == False)

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
    '''determine what the processed filename should be (based on rawfile)'''

    # store the processed files in the intermediates directory
    directory = settings.intermediates + self.directory


    # make sure a jitter directory actually exists
    zachopy.utils.mkdir(directory)

    # define the filename
    return directory + self.basename + '.processed.npy'

  def loadProcessedJitterball(self):
    '''load a pre-processed jitterball, from the intermediates directory'''

    # if not, populate the jitterball for this cadence
    self.speak(
      'populating the jitterball for {0:.0f} second cadence, '
      'based on the raw jitter file {1}.'.format(
                          self.camera.cadence,
                          self.basename))

    try:
      # if a processed file already exists, load it
      self.jitterball, self.jittermap = np.load(self.processedfile)
      self.jittercadence = self.camera.cadence
    except IOError:


      self.speak('no processed jitter file was found for {}'.format(
                          self.basename))

      self.loadUnprocessedJitterball()

  def loadUnprocessedJitterball(self):
    '''load from a raw jitter file, process, and save'''

    # otherwise, create a binned jitter structure
    self.speak('loading raw jitter from {}'.format(self.rawfile))

    # load the raw file
    self.rawdata = astropy.io.ascii.read(self.rawfile,
                    names=['t','x', 'y', 'z'])

    # subtract means
    self.rawdata['x'] -= np.mean(self.rawdata['x'])
    self.rawdata['y'] -= np.mean(self.rawdata['y'])
    self.rawdata['z'] -= np.mean(self.rawdata['z'])

    # scale jitterball to requirements (should be inflation by ~1.5)
    if self.jitterrms is not None:
      # STILL A KLUDGE! NEED ROLL, PITCH, YAW!
      original_rms = np.sqrt(np.mean( self.rawdata['x']**2 +
                      self.rawdata['y']**2))
      self.rawdata['x'] *= self.jitterrms/original_rms
      self.rawdata['y'] *= self.jitterrms/original_rms
      self.rawdata['z'] *= self.jitterrms/original_rms

    # smooth them to the required cadence
    self.speak("smoothing the jitter to {0}s cadence".format(
                          self.camera.cadence))

    #  figure out the time-spacing of the jitter timeseries
    spacings = self.rawdata['t'][1:] - self.rawdata['t'][:-1]
    spacing = np.median(spacings)

    # make sure that the jitter timeseries is evenly spaced
    aboutright = 0.01
    assert((np.abs(spacings - spacing) < spacing*aboutright).all())

    # create a convolution filter, to smooth to camera's cadence
    n = np.long(self.camera.cadence/spacing)
    filter = np.ones(n)/n

    # construct smoothed timeseries, sampled at raw time resolution
    smoothed_t = np.convolve(self.rawdata['t'], filter, mode='valid')
    smoothed_x = np.convolve(self.rawdata['x'], filter, mode='valid')
    smoothed_y = np.convolve(self.rawdata['y'], filter, mode='valid')
    smoothed_z = np.convolve(self.rawdata['z'], filter, mode='valid')

    # sample smoothed timeseries at the camera's cadence
    t = smoothed_t[::n]
    x = smoothed_x[::n]
    y = smoothed_y[::n]
    z = smoothed_z[::n]

    # plot each dimension separately
    self.speak('saving binned jitter timeseries plot')


    # create the plot of the timeseries
    plotdirectory = settings.plots + self.directory
    zachopy.utils.mkdir(plotdirectory)
    bkw = dict(alpha=0.5, color='black')
    rkw = dict(linewidth=2, alpha=0.5, marker='o', color='red')
    fi, ax = plt.subplots(3,1, sharey=True, sharex=True)
    ax[0].plot(self.rawdata['t'], self.rawdata['x'], **bkw)
    ax[0].plot(t, x, **rkw)
    ax[1].plot(self.rawdata['t'], self.rawdata['y'], **bkw)
    ax[1].plot(t, y, **rkw)
    ax[2].plot(self.rawdata['t'], self.rawdata['z'], **bkw)
    ax[2].plot(t, z, **rkw)
    ax[0].set_xlim(0,self.camera.cadence*10)
    ax[0].set_title('TESS Pointing Jitter for \n{}\nfor {}s Cadence'.format(self.basename, self.camera.cadence), fontsize=6)
    ax[0].set_ylabel('x (")')
    ax[1].set_ylabel('y (")')
    ax[2].set_ylabel('z (")')
    ax[2].set_xlabel('Time (seconds)')
    fi.savefig(plotdirectory + self.basename + '_timeseries.pdf')

    # make interpolators to keep track of the running smooth means
    ikw = dict(kind='nearest',fill_value=0,bounds_error=False)
    xip = scipy.interpolate.interp1d(smoothed_t,smoothed_x,**ikw)
    yip = scipy.interpolate.interp1d(smoothed_t,smoothed_y,**ikw)

    # assign the jittermap here, to be used for convolution in the PSF code
    xoff = (self.rawdata['x']-xip(self.rawdata['t']))/self.camera.pixelscale
    yoff = (self.rawdata['y']-yip(self.rawdata['t']))/self.camera.pixelscale

    # create a 2D jittermap (MAKE SURE THIS IS OKAY!)
    narcsec = np.maximum(np.max(np.abs(xoff)), np.max(np.abs(yoff)))#3

    # calculate the number of bins (in units of subpixels)
    bins = np.ceil(narcsec/self.camera.pixelscale*self.nsubpixelsperpixel).astype(np.int)*2 +1

    # calculate the range of the map (in units of pixels)
    range = [   [-(bins-1)/2/self.nsubpixelsperpixel,(bins-1)/2/self.nsubpixelsperpixel],
                [-(bins-1)/2/self.nsubpixelsperpixel,(bins-1)/2/self.nsubpixelsperpixel]]

    # define the jittermap as a 2D histrogram for convolution within exps
    self.jittermap = np.histogram2d(xoff, yoff,
                      bins=bins,
                      range=range,
                      normed=True)

    # define the binned jitterball, for nudges between exps
    self.jitterball = (x,y,z)

    # keep track of the jitter cadence associated with this
    self.jittercadence = self.camera.cadence

    self.speak('saving jittermap plots')
    # make an easier-to-view histogram of the jitterball for plotting
    jittermap_to_plot = np.histogram2d(xoff, yoff,
                      bins=50,
                      range=range,
                      normed=True)

    # plot this 2D histogram
    plothist2d(jittermap_to_plot,
            title='TESS Pointing Jitter over {0}s'.format(
              self.camera.cadence),
            xtitle='Pixels', ytitle='Pixels',
          filename=plotdirectory+self.basename+'_jittermaphires.pdf')

    # plot the adopted jitterball, as more useful binning
    plothist2d(self.jittermap,
            title='TESS Pointing Jitter over {0}s'.format(
              self.camera.cadence),
            xtitle='Pixels', ytitle='Pixels',
          filename=plotdirectory+self.basename+'_jittermapadopt.pdf')

    # save the necessary jitter files
    self.speak('saving the jitter files to {0}'.format(self.processedfile))
    np.save(self.processedfile, (self.jitterball, self.jittermap))

  @property
  def x(self):
    return self.amplifyinterexposurejitter*self.jitterball[0]

  @property
  def y(self):
    return self.amplifyinterexposurejitter*self.jitterball[1]

  @property
  def z(self):
    return self.amplifyinterexposurejitter*self.jitterball[2]

  def writeNudges(self, outfile='jitter.txt'):

    counters = np.arange(len(self.x))
    bjds = self.camera.counterToBJD(counters)
    time = bjds - np.min(bjds)
    plt.figure('jitter timeseries')
    gs = plt.matplotlib.gridspec.GridSpec(2,1,hspace=0.15)
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


    data= [counters, bjds, self.x, self.y]
    names=['imagenumber', 'bjd', 'arcsecnudge_ra', 'arcsecnudge_dec']

    t = astropy.table.Table(data=data, names=names)
    t.write(outfile.replace('.txt', '_amplifiedby{}.txt'.format(self.amplifyinterexposurejitter)), format='ascii.fixed_width', delimiter=' ')
    self.speak("save jitter nudge timeseries to {0}".format(outfile))

  def applyNudge(self,
                    counter=None, # which row to use from jitterball?
                    dx=None, dy=None, dz=None, # custom nudges, in arcsec
                    header=None, # the FITS header in which to record nudges
                    ):

    '''jitter the camera by a little bit,
      by introducing nudges draw from a
      (cadence-appropriate) jitterball timeseries.'''

    # make sure the jitterball has been populated
    self.load()
    n = len(self.x)

    # should we be applying a custom offset?
    usecustom = (counter is None)
    if usecustom:
      self.camera.nudge['x'] = dx
      self.camera.nudge['y'] = dy
      self.camera.nudge['z'] = dz
    else:
      # if we're over the counter, loop back
      i = counter % n
      self.camera.nudge['x'] = self.x[i]
      self.camera.nudge['y'] = self.y[i]
      self.camera.nudge['z'] = self.z[i]

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

      self.speak('updated header keywords')
    except TypeError:
      self.speak('no header was found to update')

    # move the camera, using the updated nudge values
    self.speak("nudged the camera to {x},{y}"
          " away from nominal pointing.".format(**self.camera.nudge))


def plothist2d(hist,  title=None, log=False,
            xtitle=None, ytitle=None, filename=None):

  '''Plot a 2D histogram.'''
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

  xhalf, yhalf = (x[1]-x[0])/2.0, (y[1] - y[0])/2.0
  kw =  dict(cmap='gray_r',
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
