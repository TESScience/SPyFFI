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

	def expose(self, remake=False, write=True, display=False, jitter=True):
		self.ccd = self.camera.ccds[0]
		self.ccd.display = display
		for i in range(self.nexposures):
			self.ccd.expose(write=write, remake=remake, jitter=jitter)

	def createCamera(self, cadence=2, ra=0.0, dec=0.0, subarray=100, **kwargs):
		self.camera = Camera.Camera(cadence=cadence, ra=ra, dec=dec, subarray=subarray)

	def create(self, **kwargs):
		todo = {2:3, 120:3, 1800:int(27.4*48)}
		for k in todo.keys():
			self.camera.setCadence(k)
			self.nexposures = todo[k]
			self.expose(**kwargs)

class Sky(Observation):
	def __init__(self, subarray=None, **kwargs):
		Observation.__init__(self, subarray=subarray, **kwargs)
		self.camera.testpattern=False

	def createCatalog(self, radius=0.2, **kwargs):
		# determine the size of the catalog from the camera's subarray size (in pixels)
		ra, dec = self.camera.ra, self.camera.dec
		radius = self.camera.subarray*self.camera.pixelscale/60.0/60.0
		self.camera.catalog = Catalogs.makeCatalog(name='UCAC4', ra=ra, dec=dec, radius=radius)

class SkySubarray(Sky):
	def __init__(self, subarray=200, **kwargs):
		kwargs['subarray'] = subarray
		Sky.__init__(self, **kwargs)


class SkyFFI(Sky):
	def __init__(self, **kwargs):
		kwargs['subarray'] = None
		Sky.__init__(self, **kwargs)

	def createCamera(self, subarray=None, cadence=2, ra=0.0, dec=0.0, **kwargs):
		self.camera = Camera.Camera(cadence=cadence, ra=ra, dec=dec, subarray=subarray)

	def createCatalog(self, fast=False, **kwargs):
		# determine the size of the catalog from the camera's subarray size (in pixels)
		ra, dec = self.camera.ra, self.camera.dec
		radius = self.camera.fov/np.sqrt(2)*1.01
		if fast:
			radius *= 0.1
		self.camera.catalog = Catalogs.makeCatalog(name='UCAC4', ra=ra, dec=dec, radius=radius)

	def expose(self, remake=False, write=True, display=False, jitter=True):
		self.ccd = self.camera.ccds[0]
		self.ccd.display = display
		for i in range(self.nexposures):
			self.camera.expose(write=write, remake=remake, jitter=jitter)


class TestPattern(Observation):
	def __init__(self, **kwargs):
		Observation.__init__(self, **kwargs)
		self.camera.testpattern=True

	def createCatalog(self, spacing=200.0, magnitudes=[10], ra=0.0, dec=0.0, random=True, nudge=21.1, pm=0.0, **kwargs):

		# determine the size of the catalog from the camera's subarray size (in pixels)
		assert(self.camera.subarray is not None)
		size = self.camera.pixelscale*self.camera.subarray # in arcsec
		self.camera.catalog = Catalogs.makeCatalog(name='testpattern', size=size, spacing=spacing, magnitudes=magnitudes, random=random, nudge=nudge, pm=pm)
