from imports import *
import Image

d = zachopy.display.ds9('Image')
C = Image.Camera(testpattern=True, cadence=1800.0, stamp=100)
I = Image.Image(C)
I.image = I.zeros()
for n in range(100):
    I.addCosmicsAl(version='fancy', diffusion=False, gradient=False)
d.one(I.image)
