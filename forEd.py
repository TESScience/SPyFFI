import sys
import numpy as np
from SPyFFI import Observation

try:
    magnitude = np.float(sys.argv[1])
except:
    magnitude =10
o = Observation.TestPattern(magnitudes=[magnitude], subarray=4296, nexposures=25, random=True)
o.expose(jitter=True)
