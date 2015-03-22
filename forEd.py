#!/corscorpii/d1/zkbt/Ureka/variants/common/bin/python
#SBATCH --share
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition RegNodes,HyperNodes
#SBATCH --output=/corscorpii/d1/zkbt/TESS/output.log
#SBATCH	--time=48:00:00
#SBATCH --job-name=SPyFFI

import sys
import numpy as np
from SPyFFI import Observation

magnitude = np.float(sys.argv[1])
o = Observation.TestPattern(magnitudes=[magnitude], subarray=4096, n=10, random=True)
o.expose(jitter=True)
