# noinspection PyUnresolvedReferences
from cosmical_realistic import cosmical
import scipy.signal
import numpy as np


def cosmicImage(exptime=1800.0, size=2048, rate=5.0, gradient=False, diffusion=False, buffer_size=100):
    """Generate a cosmic ray image, using Al Levine's fast C code."""

    # if a gradient is set, allow the exposure times to be different
    if gradient:
        smallexptime = exptime
        bigexptime = exptime * 2
    else:
        smallexptime = exptime * 1.5
        bigexptime = exptime * 1.5

    # add margin around image, because Al's code doesn't start cosmic ray images off the screen
    bufferedsize = size + 2 * buffer_size

    # set the diffusion flag to be 0 if False, 1 if True (as integer type)
    intdiffusion = 0  # np.int(diffusion)

    # call the fancy cosmic ray code
    image = cosmical(rate, smallexptime, bigexptime, bufferedsize, bufferedsize, intdiffusion)

    # if we need to diffuse the image, use Al's kernal (from the fancy code)
    if diffusion:
        kernal = np.array([[0.0034, 0.0516, 0.0034],
                           [0.0516, 0.7798, 0.0516],
                           [0.0034, 0.0516, 0.0034]])
        image = scipy.signal.convolve2d(image, kernal, mode='same')

    good_image = image[buffer_size:-buffer_size, buffer_size:-buffer_size] if buffer_size > 0 else image[:, :]

    # return the image
    return good_image