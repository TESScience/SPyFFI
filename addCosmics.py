from imports import *
import cosmical_original._cosmical
import cosmical_realistic._cosmical

def cosmicImage(exptime=1800.0, size=2048, gradient=False, version='fancy', diffusion=False):
    '''Generate a cosmic ray image, using Al Levine's fast C code.'''
    # the cosmic ray rate
    rate = 5.0

    # if a gradient is set, allow the exposure times to be different
    if gradient:
      smallexptime = exptime
      bigexptime = exptime*2
    else:
      smallexptime = exptime*1.5
      bigexptime = exptime*1.5

    # call different codes, depending
    if version == 'original':

        # draw an Poisson number, for how many cosmic rays should be included in this image
        physicalpixelsize = 15.0/1e4    # cm
        nexpected = (physicalpixelsize*size)**2*(0.5*smallexptime + 0.5*bigexptime)*rate
        ndrawn = np.random.poisson(nexpected)

        # call the original cosmic ray code
        image = np.transpose(cosmical_original._cosmical.cosmical(smallexptime, bigexptime, ndrawn, size, size))

        # if we need to diffuse the image, use Al's kernal (from the fancy code)
        if diffusion:
            kernal = np.array([	[0.0034, 0.0516, 0.0034],
                                [0.0516, 0.7798, 0.0516],
                                [0.0034, 0.0516, 0.0034]])
            image = scipy.signal.convolve2d(image, kernal, mode='same')
    elif version == 'fancy':

        # set the diffusion flag to be 0 if False, 1 if True (as integer type)
        intdiffusion=np.int(diffusion)

        # call the fancy cosmic ray code
        image = np.transpose(cosmical_realistic._cosmical.cosmical(rate, smallexptime, bigexptime, size, size, intdiffusion))
    print smallexptime, bigexptime
    # return the image
    return image

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
