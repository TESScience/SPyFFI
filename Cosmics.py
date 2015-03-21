from imports import *
import cosmical_realistic._cosmical

def test(n):
    for i in range(n):
        cosmicImage(size=100)

def cosmicImage(exptime=1800.0, size=2048, rate=5.0, gradient=False, version='fancy', diffusion=False, buffer=100):
    '''Generate a cosmic ray image, using Al Levine's fast C code.'''

    # if a gradient is set, allow the exposure times to be different
    if gradient:
      smallexptime = exptime
      bigexptime = exptime*2
    else:
      smallexptime = exptime*1.5
      bigexptime = exptime*1.5

    # add margin around image, because Al's code doesn't start cosmic ray images off the screen
    bufferedsize = size + 2*buffer

    # call different codes, depending
    if version == 'fancy':

        # set the diffusion flag to be 0 if False, 1 if True (as integer type)
        intdiffusion=0#np.int(diffusion)

        # call the fancy cosmic ray code
        image = cosmical_realistic._cosmical.cosmical(rate, smallexptime, bigexptime, bufferedsize, bufferedsize, intdiffusion)

    # if we need to diffuse the image, use Al's kernal (from the fancy code)
    if diffusion:
        kernal = np.array([	[0.0034, 0.0516, 0.0034],
                            [0.0516, 0.7798, 0.0516],
                            [0.0034, 0.0516, 0.0034]])
        image = scipy.signal.convolve2d(image, kernal, mode='same')

    if buffer > 0:
        goodimage = image[buffer:-buffer,buffer:-buffer]
    else:
        goodimage = image[:,:]
    del(image)
    # return the image
    return goodimage

import zachopy.display
def display(**kwargs):
    i = cosmicImage(**kwargs)
    d = zachopy.display.ds9('cosmics')
    d.one(i)
    plt.ion()
    plt.figure('sum over columns')
    plt.plot(np.sum(i,1))
    plt.draw()
    return i
