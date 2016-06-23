from imports import *
import Camera, Cube, Photometer, settings

class Track(object):
    '''Quick way to make a path for a pixel to follow. Must be created with a Cartographer.'''
    def __init__(self, cartographer=None, n=100):
        assert(cartographer is not None)
        self.n = n
        self.cartographer = cartographer
        self.indices = np.arange(n)
        self.steps = np.linspace(0, 1, n)

        self.make()
        self.findCenter()

    def xy(self):
        return self.focalxy.tuple

    def findCenter(self):
        x, y = self.positions.focalxy.tuple
        self.center = self.cartographer.point(np.mean(x), np.mean(y), 'focalxy')


class Spiral(Track):
    def __init__(self, x=0, y=0, **kwargs):
        Track.__init__(self, **kwargs)

        #self.center(x,y)

    def make(self, rmax=2.5, thetamax=7*np.pi):
        r, theta = self.steps*rmax, self.steps*thetamax
        self.positions = self.cartographer.point(r, theta, 'focalrtheta')

class Line(Track):
    def __init__(self, x=0, y=0, **kwargs):
        Track.__init__(self, **kwargs)
        self.make()
        #self.center(x,y)

    def make(self, xmax=2.5, ymax=1.0):
        x, y = xmax*np.linspace(-1, 1, self.n), ymax*np.linspace(-1, 1, self.n)
        self.positions = self.cartographer.point(x, y, 'focalxy')


class Tester(Talker):
    '''Tester object, to test various SPyFFI components
        (and also leave variables and attributes availbe to play with interactively.)'''
    def __init__(self, subarray=200, testpattern=True, cadence=1800):
        Talker.__init__(self)
        self.camera = Camera.Camera(subarray=subarray, testpattern=testpattern, cadence=cadence)
        self.buster = self.camera.cartographer
        self.jitter = self.camera.jitter

    def photometer(self, remake=False, size=100, cadence=2, n=100, jitter=False):
        self.cube = Cube.Cube(n=n,size=size, cadence=cadence, jitter=jitter)
        self.cube.camera.populateCatalog(random=True, magnitudes=[6,16])
        self.cube.load(remake=remake)
        self.phot = Photometer.Photometer(self.cube)

    def cartographer(self):
        '''Make sure the Cartographer can convert coordinates back and forth.'''
        B = self.buster
        for i in B.possibleinputs:
            B.point(0,0,i)
            for o in B.possibleoutputs:
                B.quote(o)

    def jitter(self):
        self.jitter = self.camera.jitter

    def image(self, remake=True, jitter=False):
        '''Make sure the CCDs can exposure new images.'''
        self.ccd = self.camera.ccds[0]
        self.ccd.expose(write=True, remake=remake, jitter=jitter)

    def shiftPSF(self):
        '''Make sure we can shift PSFs around.'''
        C = self.camera
        P = C.psf

        directory = settings.prefix + 'plots/psfs/'
        zachopy.utils.mkdir(directory)
        resolution = 60

        # nudgeds
        this = directory + 'shiftByNudges/'
        print "testing " + this
        zachopy.utils.mkdir(this)
        counter = 0
        for theta in np.linspace(0,7*np.pi,resolution):
            position = P.cartographer.point(0, 0,'focalxy')
            radius = theta/7/np.pi
            P.binHighResolutionPSF(position,dx=radius*np.cos(theta),dy=radius*np.sin(theta),temperature=5000)
            filename = this + '{0:04.0f}.pdf'.format(counter)
            plt.savefig(filename)
            print " saved plot to " + filename
            counter += 1


        # translation
        this = directory + 'shiftByTranslation/'
        print "testing " + this
        zachopy.utils.mkdir(this)
        counter = 0
        for x in np.linspace(-2500,2500,resolution):
            position = P.cartographer.point(x,0.0,'focalxy')
            P.binHighResolutionPSF(position,dx=0,dy=0,temperature=5000)
            filename = this + '{0:04.0f}.pdf'.format(counter)
            plt.savefig(filename)
            print " saved plot to " + filename
            counter += 1

        # rotation
        this = directory + 'shiftByRotation/'
        print "testing " + this
        zachopy.utils.mkdir(this)
        counter = 0
        for theta in np.linspace(0,2*np.pi,resolution):
            position = P.cartographer.point(1500, theta,'focalrtheta')
            P.binHighResolutionPSF(position,dx=0,dy=0,temperature=5000)
            filename = this + '{0:04.0f}.pdf'.format(counter)
            plt.savefig(filename)
            print " saved plot to " + filename
            counter += 1


    def lightcurveFromPSF(self):
        '''Make sure we get a reasonable light curve when we shift a star around on the detector (without PSF library).'''
        C = self.camera
        P = C.psf
        resolution=300
        flux = np.zeros(resolution)
        plt.figure('lightcurve test')
        x = np.linspace(0.0, 2.0,resolution)
        for i in np.arange(resolution) :
            position = P.cartographer.point(1 - (1-x[i])**2,np.pi*2*x[i],'focalrtheta')
            pixelized = P.newlyPixelizedPSF(position, verbose=True)
            flux[i] = np.sum(pixelized)

        plt.plot(x, flux, marker='o', alpha=0.5)

    def libraryPSF(self, n=100):
        '''Make sure our PSF library calculations give nearly the same answer as redoing pixelization integrals from scratch.'''
        directory = settings.prefix + 'plots/psfs/'

        C = self.camera
        P = C.psf

        tracks = [Line, Spiral]
        for t in tracks:
            track = t(cartographer=P.cartographer, n=n)
            this = directory + 'library{0}Offsets/'.format(track.__class__.__name__)
            zachopy.utils.mkdir(this)
            r, theta = track.positions.focalrtheta.tuple

            library, new = np.zeros(track.n).astype(np.float), np.zeros(track.n).astype(np.float)

            for i in np.arange(track.n):
                self.speak('testing the library PSFs at position {0:.0f}/{1:.0f} along the track'.format(i, track.n))
                position = P.cartographer.point(r[i],theta[i], 'focalrtheta')
                library[i], new[i] = P.comparePSFs(position, center=track.center)
                trackx, tracky = track.positions.ccdxy.tuple
                P.axLightcurve.cla()
                P.axLightcurve.plot(track.steps[:i+1], library[:i+1], linewidth=2, alpha=0.5, color='orange', label='library')
                P.axLightcurve.plot(track.steps[:i+1], new[:i+1], linewidth=2, alpha=0.5, color='blue', label='recalculated')
                P.axLightcurve.set_xlim(0,1)
                P.axLightcurve.set_ylim(0.971, 1.001)
                for ax in [P.axNew, P.axLibrary]:
                    ax.plot(trackx, tracky, alpha=0.25, linewidth=2, color='red')
                    ax.plot(trackx[i], tracky[i], linewidth=0, marker='+', color='red', alpha=0.5, markersize=20, markeredgewidth=4)
                P.axLightcurve.plot(track.steps[i], library[i], linewidth=0, marker='+', markeredgecolor='orange', alpha=0.5, markersize=20, markeredgewidth=4)
                P.axLightcurve.plot(track.steps[i], new[i], linewidth=0, marker='+', color='blue', alpha=0.5, markersize=20, markeredgewidth=4)
                plt.legend(frameon=False)
                plt.setp(P.axLightcurve.get_xticklabels(), visible=False, fontsize=8)

                plt.draw()

                filename = this + '{0:04.0f}.pdf'.format(i)
                plt.savefig(filename)
