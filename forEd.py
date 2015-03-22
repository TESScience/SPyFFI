import sys
import numpy as np
from SPyFFI import Observation

# figure out the TESS magnitude for the test pattern stars
try:
    magnitude = np.float(sys.argv[1])
except:
    magnitude =10

# create an observation cube, initializing the magnitude range, the subarray size, the number of exposures, and whether or not to nudge each star by a random fraction of a pixel (necessary to prevent aliasing)
o = Observation.TestPattern(magnitudes=[magnitude], subarray=4296, nexposures=25, random=True)

# create the images and store them to output directory
o.expose(jitter=False)
