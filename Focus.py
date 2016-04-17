'''Keep track of camera focus (to make cartoon drifts).'''
from imports import *
import settings

class Focus(Talker):
  def __init__(self, camera=None, span=[0.0,10.0]):
      Talker.__init__(self, mute=False, pithy=False)

      # store the input camera
      self.camera = camera

      # what do you want the absolute focus range to be?
      self.span = span


      self.orbit = 13.7
      self.sigma = 0.1
      self.jump = 0.5
      self.best, self.worst = self.span
      self.scale = self.worst - self.best

  def model(self, counter):
      bjd = self.camera.counterToBJD(counter)
      phase = (((bjd - self.camera.bjd0)/self.orbit + 0.5) % 1) - 0.5
      best, worst = self.span

      return np.exp(-0.5*phase**2/self.sigma**2)*self.scale + best

  def writeModel(self, outfile):
      time = np.linspace(0, 2*self.orbit, 1000)
      counter = 24*60*60/self.camera.cadence*time
      plt.figure('focus timeseries')
      assert(time.shape == counter.shape)
      plt.plot(time, self.model(counter), linewidth=2)
      plt.xlabel('Time from Observation Start (days)')
      plt.ylabel('Focus (um)')
      plt.xlim(np.min(time), np.max(time))
      plt.ylim(*self.span)
      plt.draw()
      plt.savefig(outfile.replace('.txt', '.pdf'))

      with open(outfile, 'w') as f:
          f.write('x = (((bjd - {})/{} + 0.5) % 1) - 0.5\n'.format(self.camera.bjd0, self.orbit))
          f.write('focus = {best} + {scale}*exp(-0.5*(x/{sigma})**2)'.format(**self.__dict__))


  def applyNudge(self,
                    counter=None, # which row to use from jitterball?
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
      self.camera.nudge['x'] = scale*self.x[i]
      self.camera.nudge['y'] = scale*self.y[i]
      self.camera.nudge['z'] = scale*self.z[i]

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
      header['JITSCALE'] = (scale, 'jitter magnified by ? relative to file')
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
