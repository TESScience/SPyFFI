'''Wrappers to create a stack of observations,
of any size
of any cadence
of any length
of either the sky or a test pattern.'''


import Camera, Catalogs
from imports import *

class Observation(Talker):
	def __init__(self, nexposures=1, **kwargs):

		Talker.__init__(self)
		self.nexposures = nexposures
		self.createCamera(**kwargs)
		self.createCatalog(**kwargs)


class TestPattern(Observation):
	def __init__(self, **kwargs):
		Observation.__init__(self, **kwargs)
		self.camera.testpattern=True

	def createCamera(self, cadence=2, ra=0.0, dec=0.0, subarray=100, **kwargs):
		self.camera = Camera.Camera(cadence=cadence, ra=ra, dec=dec, subarray=subarray)

	def createCatalog(self, spacing=200.0, magnitudes=[10], ra=0.0, dec=0.0, random=True, nudge=21.1, pm=0.0, **kwargs):

		# determine the size of the catalog from the camera's subarray size (in pixels)
		assert(self.camera.subarray is not None)
		size = self.camera.pixelscale*self.camera.subarray # in arcsec
		self.camera.catalog = Catalogs.makeCatalog(name='testpattern', size=size, spacing=spacing, magnitudes=magnitudes, random=random, nudge=nudge, pm=pm)

	def expose(self, remake=False, write=True, display=False, jitter=True):
		self.ccd = self.camera.ccds[0]
		self.ccd.display = display
		for i in range(self.nexposures):
			self.ccd.expose(write=write, remake=remake, jitter=jitter)
