'''Wrappers to create a stack of observations,
of any size
of any cadence
of any length
of either the sky or a test pattern.'''


import Camera, Catalogs
from imports import *
from defaults import inputs as default

class Observation(Talker):
    '''an observation object handles a simulated group of observations,
        using the same general pointing, the same stars, the same lightcurves'''

    def __init__(self, inputs=default):
        '''initialize a basic observation object'''

        # store the dictionary of dictionaries of inputs
        self.inputs = inputs

        # initialize the talker
        Talker.__init__(self)

        # print the inputs for this observation
        self.speak("  creating a new observation, with the following inputs:")
        for k in inputs.keys():
            self.speak('   {:>20s}'.format('inputs[{}] = '.format(k)))
            for l in inputs[k].keys():
                self.speak('     {:>30s}:{}'.format(l, inputs[k][l]))

        # setup the basics of the observation
        self.setupObservation()

        # create the camera
        self.createCamera()

        # create the catalog
        self.createCatalog()

    def setupObservation(self):
        '''set up the basics of this observation set'''
        kw = self.inputs['observation']
        self.cadencestodo = kw['cadencestodo']
        self.collate = kw['collate']

        name = self.inputs['catalog']['name'].lower()
        self.testpattern = name == 'testpattern'

    def expose(self):
        '''execute one exposure of this observation, looping through all CCDs'''

        # create one exposure, by looping over the CCDs
        #   (the last CCD will update the counter)
        for i, c in enumerate(self.camera.ccds):
            c.expose(**self.inputs['expose'])

    def create(self, **kwargs):
        '''make *all* the exposures for this observation,
           looping through all CCD's, either collating or not'''

        # loop over the cadences that need to be done
        for k in self.cadencestodo.keys():
            # set the cadence to this one
            self.camera.setCadence(k)
            np.save(self.camera.directory + 'observationdictionary.npy', self.inputs)
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
                    self.counter = 0
                    for i in range(self.cadencestodo[k]):
                        # expose this CCD
                        c.expose(advancecounter=False, **self.inputs['expose'])
                        # advance the counter by hand
                        c.camera.advanceCounter()

    def createCamera(self):

        self.speak('setting up the camera for this Observation.')
        kw = self.inputs['camera']
        self.camera = Camera.Camera(testpattern=self.testpattern, **kw)
        #try:
        #    self.camera.label = kw['label']
        #except KeyError:
        #    pass

    def createCatalogFromStars(self):
        # determine the catalog purview from the camera object
        ra, dec = self.camera.ra, self.camera.dec
        radius = self.camera.effective_fov*1.01
        kw = self.inputs['catalog']['skykw']

        self.speak('creating catalog from "real" stars')
        self.camera.catalog = Catalogs.UCAC4(
                                ra=ra, dec=dec, radius=radius,
                                lckw=self.inputs['catalog']['lckw'],
                                **kw)

    def createCatalogWithTestPattern(self):
        # determine the catalog purview from the camera object
        ra, dec = self.camera.ra, self.camera.dec
        size = 2*self.camera.effective_fov*1.01*3600.0

        self.speak('creating catalog representing a test pattern of stars')
        kw = self.inputs['catalog']['testpatternkw']
        self.camera.catalog = Catalogs.TestPattern(
                                ra=ra, dec=dec, size=size,
                                lckw=self.inputs['catalog']['lckw'],
                                **kw)

    def createCatalog(self):
        try:
            self.speak('setting up catalog based on '
                'camera centered at {x.ra:.2f}, {x.dec:.2f}'.format(x=self.camera))
        except AttributeError:
            raise RuntimeError('createCamera must be run before createCatalog')


        if self.testpattern:
            self.createCatalogWithTestPattern()
        else:
            self.createCatalogFromStars()
