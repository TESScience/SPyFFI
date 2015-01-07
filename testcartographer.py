from imports import *

import Camera
reload(Camera)
c = Camera.Camera(stamp=100)
c.cartographer.point(c.catalog.ra, c.catalog.dec, 'celestial')
plt.ion()
plt.cla()
xy = c.cartographer.quote('focalxy')
plt.plot(*xy.arrays, marker='o', linewidth=0)
