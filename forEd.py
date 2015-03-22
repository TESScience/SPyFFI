import sys
import numpy as np
from SPyFFI import Observation

magnitude =10.0; np.float(sys.argv[1])
o = Observation.TestPattern(magnitudes=[magnitude], subarray=100, n=10)
o.expose()

