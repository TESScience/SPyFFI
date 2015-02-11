from SPyFFI.imports import *
import SPyFFI.settings as settings
from SPyFFI.Cube import Cube
from SPyFFI.Photometer import Photometer
from SPyFFI.Noise import noise
from Strategies import *
import textwrap
import os

#cube = Cube(cadence=120, size=32, n=100, stacker='Central 8 out of 10')
cadence = 120
cube = Cube(cadence=cadence, size=256, n=25, stacker='Rejecting Outliers')
cube.camera.populateCatalog(random=True, magnitudes=[6,16])
cube.load()
phot = Photometer(cube)
phot.drawApertures()
for a in phot.apertures:
    a.measure(cube)

mag = np.array([a.mag for a in phot.apertures])
n = np.array([a.n for a in phot.apertures])
frac = np.array([np.mean(a.npixelsaffected.astype(np.float)/a.n) for a in phot.apertures])

b = 1.0
bx, bf, be = zachopy.oned.binto(mag, frac, b)
bx, bn, be = zachopy.oned.binto(mag, n, b)

for i in range(len(bx)):
    print '{0:2.0f} to {1:2.0f}: {2:4.1f} pixels in aperature, {3:.3f} affected'.format(bx[i] - b/2.0, bx[i]+b/2.0, bn[i], bf[i])

plt.figure('fraction of pixels affected')
plt.cla()
plt.plot(mag, frac, linewidth=0, marker='o', alpha=0.5, color='gray')
plt.plot(bx, bf, linewidth=5, alpha=0.75, color='red')
plt.title('{0}s exposures'.format(cadence))
plt.xlabel('TESS magnitude')
plt.ylabel('Fraction of Pixels Affect by Outlier Rejector')
plt.savefig(settings.dirs['plots'] + 'fractionaffected{0}.pdf'.format(cadence))
