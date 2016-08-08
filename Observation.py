"""Wrappers to create a stack of observations,
of any size
of any cadence
of any length
of either the sky or a test pattern."""
import os

import Camera
import Catalogs
import numpy as np
from debug import DebugDict
from defaults import inputs as default
import logging
from settings import log_file_handler

logger = logging.getLogger(__name__)
logger.addHandler(log_file_handler)


class Observation(object):
    """an observation object handles a simulated group of observations,
        using the same general pointing, the same stars, the same lightcurves"""

    def __init__(self, inputs=default, debug=False):
        """initialize a basic observation object"""
        from documentation.input import documented_keywords

        # store the dictionary of dictionaries of inputs
        self.inputs = inputs if not debug else DebugDict(inputs, documented_keywords)

        # print the inputs for this observation
        logger.info("  creating a new observation, with the following inputs:")
        for k in inputs.keys():
            logger.info('   {:>20s}'.format('inputs[{}] = '.format(k)))
            for l in inputs[k].keys():
                logger.info('     {:>30s}:{}'.format(l, inputs[k][l]))

        # setup the basics of the observation
        self.cadencestodo = self.inputs['observation']['cadencestodo']
        self.collate = self.inputs['observation']['collate']
        self.testpattern = self.inputs['catalog']['name'].lower() == 'testpattern'

        # create the camera
        self.createCamera()

        # create the catalog
        self.createCatalog()

    def expose(self):
        """execute one exposure of this observation, looping through all CCDs"""

        # create one exposure, by looping over the CCDs
        #   (the last CCD will update the counter)
        for i, c in enumerate(self.camera.ccds):
            c.expose(**self.inputs['expose'])

    def create(self):
        """make *all* the exposures for this observation,
           looping through all CCD's, either collating or not"""

        # loop over the cadences that need to be done
        for k in self.cadencestodo.keys():
            # set the cadence to this one
            self.camera.setCadence(k)
            np.save(os.path.join(self.camera.directory, 'observationdictionary.npy'), self.inputs)
            # reset the counter
            self.camera.counter = 0

            if self.collate:
                # if collating, loop through exposure numbers, exposing all CCDs
                for i in range(self.cadencestodo[k]):
                    for c in self.camera.ccds:
                        c.expose(**self.inputs['expose'])

            else:
                # if not collating, loop through CCDs, exposing all for each
                for c in self.camera.ccds:
                    # reset the counter to 0
                    self.camera.counter = 0
                    for i in range(self.cadencestodo[k]):
                        # expose this CCD
                        c.expose(advancecounter=False, **self.inputs['expose'])
                        # advance the counter by hand
                        c.camera.advanceCounter()

    def createCamera(self):

        logger.info('setting up the camera for this Observation.')
        kw = self.inputs['camera']
        self.camera = Camera.Camera(testpattern=self.testpattern, **kw)
        # try:
        #    self.camera.label = kw['label']
        # except KeyError:
        #    pass

    def createCatalogFromStars(self):
        # determine the catalog purview from the camera object
        ra, dec = self.camera.ra, self.camera.dec
        radius = self.camera.effective_fov * 1.01
        kw = self.inputs['catalog']['skykw']

        logger.info('creating catalog from "real" stars')
        self.camera.catalog = Catalogs.UCAC4(
            ra=ra, dec=dec, radius=radius,
            lckw=self.inputs['catalog']['lckw'],
            **kw)

    def createCatalogWithTestPattern(self):
        # determine the catalog purview from the camera object
        ra, dec = self.camera.ra, self.camera.dec
        size = 2 * self.camera.effective_fov * 1.01 * 3600.0

        logger.info('creating catalog representing a test pattern of stars')
        kw = self.inputs['catalog']['testpatternkw']
        self.camera.catalog = Catalogs.TestPattern(
            ra=ra, dec=dec, size=size,
            lckw=self.inputs['catalog']['lckw'],
            **kw)

    def createCatalog(self):
        try:
            logger.info('setting up catalog based on '
                        'camera centered at {x.ra:.2f}, {x.dec:.2f}'.format(x=self.camera))
        except AttributeError:
            raise RuntimeError('createCamera must be run before createCatalog')

        if self.testpattern:
            self.createCatalogWithTestPattern()
        else:
            self.createCatalogFromStars()
