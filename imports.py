'''Things that will probably need to be imported into most of the TESS code chunks.'''

# some basics
import numpy as np, matplotlib.pyplot as plt
import os, copy, subprocess, glob

# some scipy tools for interpolation and image filtering
import scipy.ndimage, scipy.signal, scipy.interpolate

# lots from astropy
import astropy.io, astropy.units, astropy.coordinates, astropy.wcs, astropy.table

# general tools from my library
import zachopy.utils, zachopy.borrowed.crossfield, zachopy.spherical, zachopy.oned, zachopy.star

# a parent class that allows text reporting to be muted/unmuted
from zachopy.Talker import Talker

#from memory_profiler import profile

# what is the path to the main code directory?
codedir = os.path.abspath(os.path.dirname(__file__))

# what is the base data output directory?
