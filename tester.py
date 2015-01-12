from imports import *
import Camera, settings


def cartographer():
    C = Camera.Camera(testpattern=True, subarray=100)
    I = C.ccds[0]
    B = C.cartographer
    for i in B.possibleinputs:
        B.point(0,0,i)
        for o in B.possibleoutputs:
            B.quote(o)

def image():
    C = Camera.Camera(testpattern=True, subarray=100)
    I = C.ccds[0]
    I.expose(write=True)

def psfshifts():
    C = Camera.Camera()
    P = C.psf

    directory = settings.prefix + 'plots/psfs/'
    zachopy.utils.mkdir(directory)
    resolution = 60

    # nudgeds
    this = directory + 'nudges/'
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
    this = directory + 'translation/'
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
    this = directory + 'rotation/'
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


def psflightcurve():

    C = Camera.Camera()
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

def psflibrary():
    directory = settings.prefix + 'plots/psfs/'
    this = directory + 'spiral/'
    zachopy.utils.mkdir(this)
    C = Camera.Camera()
    P = C.psf
    resolution = 200
    x = np.linspace(0, 1, resolution)
    ints = np.arange(resolution)
    r, theta = x*2.5, np.pi*5*x
    track = P.cartographer.point(r,theta,'focalrtheta')
    library, new = np.zeros_like(x), np.zeros_like(x)
    trackx, tracky = track.focalxy.tuple
    trackx -= track.ccd.center[0]
    tracky -= track.ccd.center[1]
    for i in np.arange(resolution):

        position = P.cartographer.point(r[i],theta[i], 'focalrtheta')
        library[i], new[i] = P.comparePSFs(position)

        P.axLightcurve.cla()
        P.axLightcurve.plot(ints[:i+1], library[:i+1], linewidth=2, alpha=0.5, color='orange', label='library')
        P.axLightcurve.plot(ints[:i+1], new[:i+1], linewidth=2, alpha=0.5, color='blue', label='recalculated')
        P.axLightcurve.set_xlim(0, resolution)
        for ax in [P.axNew, P.axLibrary]:
            ax.plot(trackx, tracky, alpha=0.25, linewidth=2, color='red')
            ax.plot(trackx[i], tracky[i], linewidth=0, marker='+', color='red', alpha=0.5, markersize=20, markeredgewidth=4)
        P.axLightcurve.plot(ints[i], library[i], linewidth=0, marker='+', markeredgecolor='orange', alpha=0.5, markersize=20, markeredgewidth=4)
        P.axLightcurve.plot(ints[i], new[i], linewidth=0, marker='+', color='blue', alpha=0.5, markersize=20, markeredgewidth=4)
        plt.legend()
        plt.draw()

        filename = this + '{0:04.0f}.pdf'.format(i)
        plt.savefig(filename)
