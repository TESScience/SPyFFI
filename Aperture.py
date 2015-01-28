from imports import *


class Aperture(Talker):
    def __init__(self):
        # decide whether or not this CCD is chatty
        Talker.__init__(self, mute=False, pithy=False)

    def define(self, image):
        '''Using a master image and an input catalog, define the pixel aperture for a star.'''

def test():
  pix = 10
  nstar = np.sort(np.random.rand(pix)*10)[::-1]
  nsky = 10.0*np.ones_like(nstar)

  signal = np.cumsum(nstar)
  noise = np.sqrt(np.cumsum((nstar + nsky)))

  plt.plot(nstar, signal/noise)
  return signal, noise
