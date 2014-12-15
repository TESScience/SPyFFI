'''A few simple tests of the TESS image simulation tools.'''
import tess
from imports import *




def sub(ra=270, dec=66.56070833333332, cadence=2, n=900, stamp=100, jitter=False):
	C = tess.Camera(stamp=stamp)
	C.setCadence(cadence)
	C.point(ra, dec)
	I = tess.Image(C)
	for i in range(n):
		I.expose(write=True, split=False, jitter=jitter)


def create(ra=270, dec=66.56070833333332, cadence=2, n=900):
	C = tess.Camera()
	C.setCadence(cadence)
	C.point(ra, dec)
	I = tess.Image(C)
	for i in range(n):
		I.expose(write=True)

def subtractPairs(ra=82, dec=1,cadence=2):
	C = tess.Camera()
	C.setCadence(cadence)
	C.point(ra, dec)
	I = tess.Image(C)
	g = glob.glob(I.directory + 'final*')
	for i in range(len(g)-1):
		this = I.loadFromFITS(g[i])
		that = I.loadFromFITS(g[i+1])
		filename = I.directory + g[i+1].split('/')[-1].split('.')[0] + '-' + g[i].split('/')[-1].split('.')[0] + '.fits'
		I.writeToFITS(that - this, filename)


def threefields():
	# north celestial pole
	create(ra=270, dec=66.56070833333332, cadence=1800, n=1)
	# north galactic pole
	create(ra=192.25, dec=27.4, cadence=1800, n=1)
	# galatic center
	create(ra=266.4166666666667, dec=-29.00777777777778, cadence=1800, n=1)

def forPT():
	pointings = [(270,66.56070833333332), (192.25,27.4), (266.4166666666667,-29.00777777777778)]
	C = tess.Camera()
	for pos in pointings:
		ra, dec = pos
		C.setCadence(1800)
		C.point(ra, dec)
		C.project(write=True)
