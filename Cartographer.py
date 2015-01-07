from imports import *

class Cartographer(Talker):
    '''An object to handle all conversions between coordinate systems.'''
    def __init__(self, camera=None, ccd=None, verbose=True):
        '''Initialize Cartographer, telling it where it is and how to act.'''

        # decide whether or not this Cartographer is chatty
        Talker.__init__(self, mute=False, pithy=False)

        # these will be necessary for more complicated conversions
        self.camera, self.ccd = None, None
        self.setCamera(camera)
        self.setCCD(ccd)


    def updatePossibilities(self):
        '''Update which kinds of cartographic coordinates are valid. Should be run every time camera or ccd are updated.'''

        # no matter what, Cartographer should be able to deal with these
        possibilities = ['focalxy','focalrtheta']

        # if possible, figure out how to convert between focal plane coordinates and the sky
        if self.camera is not None:
            possibilities.extend(['celestial','ecliptic','galactic'])

        # if possible, figure out how to convert betweeen focal plane and pixel coordinates
        if self.ccd is not None:
            possibilities.extend(['ccdxy','intrapixelxy'])

        # update the possibility lists
        self.possibleinputs = possibilities
        self.possibleoutputs = possibilities

    def setCCD(self, ccd):
        '''Tell this Cartographer to think about a particular CCD object.'''
        self.ccd = ccd
        self.updatePossibilities()

    def setCamera(self, camera):
        '''Tell this Cartographer to think about a particular Camera object.'''
        self.camera = camera
        self.updatePossibilities()



    def point(self, a, b=None, type='focalxy'):
        '''Ask Cartographer to point somewhere, using any (valid) kind of coordinate.

            inputs are *either*:
                point(position) where coord is a bonafide position object
                or
                point(a, b, type) where (a,b) are the two coordinate values and type is the kind of position they should be interpreted as'''

        # allow us to input either a coordinate object, or a pair of coordinates
        self.speak('Cartographer pointing a new coordinate:', 1)
        try:
            assert(a.__class__.__base__.__name__ == 'position')
            type = a.__class__.__name__
            temp = a
            a, b = temp.a, temp.b
            del temp
        except:
            pass
        #self.speak('using [{a}, {b}, {type}] as input'.format(a=a,b=b, type=type), 2)


        # make sure the input coordinates can be understood
        assert(type in self.possibleinputs)
        self.speak('input coordinate is parsable', 2)

        # remove any previous coordinate definitions
        self.speak('clearing previous coordinates', 2)
        for coord in self.possibleinputs:
            try:
                del self.__dict__[coord]
            except:
                pass

        # assign the coordinates to Cartographer's memory
        self.input = type
        self.coordinate = eval('{type}(a,b,Cartographer=self)'.format(a=a,b=b, type=type))
        self.__dict__[self.input] = self.coordinate
        self.speak('now set to {0}\n'.format(self.coordinate), 2)

    def quote(self, type='focalrtheta'):
        '''Ask Cartographer to say where it's pointing, using any (valid) kind of coordinate.'''

        # make sure the output is possible
        assert(type in self.possibleinputs)

        # use the coordinates property definitions to return the desired output
        output = eval('self.coordinate.{type}'.format(type=type))
        self.speak('Output coordinate is {0}'.format(output), 2)
        return output

    def focalxy(self, *args):
        '''Wrapper to quickly return the focalxy coordinate.'''
        self.point(*args)
        return self.quote('focalxy')

    def focalxy(self, *args):
        '''Wrapper to quickly return the focalxy coordinate.'''
        self.point(*args)
        return self.quote('focalxy')


class position(object):
    '''General (a,b) coordinate object. Each coordinate can be either a scalar, or an N-dimensional array.'''
    def __init__(self, a, b, Cartographer=None):
        self.a = a
        self.b = b
        assert(np.size(a) == np.size(b))
        self.Cartographer = Cartographer
        self.aname = 'a'
        self.bname = 'b'
        self.name = 'position'

    def __str__(self):
        '''How should these coordinates be represented as a string?'''
        if np.size(self.a) == 1:
            return '{name} ({aname}, {bname}) = ({a},{b})'.format(**self.__dict__)
        elif np.size(self.a) > 1:
            n = np.size(self.a)
            return '{name} ({aname}, {bname}) = arrays of {n} elements'.format(n=n, **self.__dict__)
        else:
            return 'nothing!'

    @property
    def arrays(self):
        '''For ease of plotting, return the coordinates at a tuple of arrays.'''
        return self.a, self.b

class focalxy(position):
    '''Cartesian coordinates in the focal plane. Measured from center of field.'''
    def __init__(self,a,b,Cartographer=None):
        position.__init__(self,a,b,Cartographer)
        self.aname, self.bname, self.name = 'x','y','focalplane'
        self.x, self.y = self.a, self.b

    # use properties to define the most efficient conversion pathways to get from this coordinate to others
    @property
    def focalxy(self):
        return self

    @property
    def focalrtheta(self):
        radius = np.sqrt(self.x**2 + self.y**2)
        theta = np.arctan2(self.y, self.x)
        return focalrtheta(radius, theta)



class focalrtheta(position):
    '''Polar coordinates in the focal plane. Measured from center of field.'''
    def __init__(self,a,b,Cartographer=None):
        position.__init__(self,a,b,Cartographer)
        self.aname, self.bname, self.name = 'r','theta','focalplane'
        self.r, self.theta = self.a, self.b

    # use properties to define the most efficient conversion pathways to get from this coordinate to others
    @property
    def focalxy(self):
        x = self.r*np.cos(self.theta)
        y = self.r*np.sin(self.theta)
        return focalxy(x, y)

    @property
    def focalrtheta(self):
        return self



class celestial(position):
    '''R.A. and Dec. on the sky.'''
    def __init__(self,a,b,Cartographer=None):
        position.__init__(self,a,b,Cartographer)
        self.aname, self.bname, self.name = 'ra', 'dec','sky'
        self.ra, self.dec = self.a, self.b

    # use properties to define the most efficient conversion pathways to get from this coordinate to others
    @property
    def focalxy(self):
        x, y = self.Cartographer.camera.wcs.wcs_world2pix(self.ra, self.dec, 1)
        xcenter, ycenter = self.Cartographer.camera.wcs.wcs.crpix
        return focalxy(x - xcenter, y-ycenter)

    @property
    def focalrtheta(self):
        return self.focalxy.focalrtheta

    @property
    def celestial(self):
        return self












class ds9(position):
    def __init__(self, a, b,Cartographer=None):
        position.__init__(self,a,b,Cartographer)
        self.aname, self.bname, self.name = 'x', 'y', 'ds9'

    @property
    def imshow(self):
        return imshow(self.a-1, self.b-1)

class imshow(position):
    def __init__(self, a, b,Cartographer=None):
        position.__init__(self,a,b,Cartographer)
        self.aname, self.bname = 'x', 'y'
        self.name = 'imshow'
