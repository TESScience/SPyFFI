import Cosmics
import zachopy.display
import matplotlib.pyplot as plt
import numpy as np

def demoimages():
    '''Create example images at three cadences.'''

    d = zachopy.display.ds9('Image')

    for cadence in [2.0, 120.0, 1800.0]:
        image = Cosmics.cosmicImage(exptime=cadence, size=512, rate=5.0, gradient=False, version='fancy', diffusion=False)
        d.one(image, clobber=cadence==2.0)
    print "fiddle with the images until you like them, and then save as example_cosmics.png"

# estimate the PDF of cosmic ray flux per pixel
def histogram(remake=False, diffusion=False, n=10, size=4096, xpeak = 700, version='fancy'):

    plt.ion()
    gs = plt.matplotlib.gridspec.GridSpec(2,1,wspace=0, hspace=0.05)
    ax = []
    ax.append(plt.subplot(gs[0]))
    ax.append(plt.subplot(gs[1], sharex=ax[0]))
    plt.setp(ax[0].get_xticklabels(), visible=False)

    ax[0].set_title('Effect of Cosmic Rays on 2s Subexposures')
    ax[0].set_yscale('log')
    ax[0].set_ylabel('probability (density) of\ncosmic ray/pixel/subexposure')
    ax[1].set_xlim(0,2500)
    ax[0].set_ylim(3e-9,1e-3)
    ax[1].set_ylabel('cumulative probability')
    ax[1].set_xlabel('# of electrons from cosmic ray')


    ax[1].set_yscale('log')
    ax[1].set_ylim(1e-6, 1e-2)

    def pdf(y):
        return y.astype(np.float)/np.sum(y)

    def ccdf(y):
        return 1.0 - np.cumsum(y.astype(np.float))/np.sum(y)

    filename='cosmics_histogram.npy'
    try:
        assert(remake == False)
        x, y = np.load(filename)
    except:

        sum = np.zeros(2500)



        for count in range(n):
            image = Cosmics.cosmicImage(exptime=2.0, size=size, rate=5.0, gradient=False, version=version, diffusion=diffusion)
            y, edges = np.histogram(image, bins=np.arange(2501))
            x = 0.5*(edges[:-1] + edges[1:])
            kw = dict(color='orange', alpha=0.1)
            ax[0].plot(x, pdf(y), **kw)
            ax[1].plot(x, ccdf(y), **kw)
            #
            sum += y
            print count
            plt.draw()

        y = sum/(count + 1.0)
        np.save(filename, (x,y))

    kw = dict(linewidth=4, color='black', alpha=0.75)
    ax[0].plot(x, pdf(y), **kw)
    ax[1].plot(x, ccdf(y), **kw)
    for a in ax:
        kw = dict(color='green', alpha=0.5)
        a.axvline(xpeak, linewidth=4, linestyle='--', **kw)
    ax[1].axhline(ccdf(y)[xpeak], linewidth=4,linestyle='--', **kw)
    ax[1].text(1300, 4e-5, '{0:.4f} of pixels are hit by\nat least {1} electrons \nin a 2s subexposure'.format(ccdf(y)[xpeak],xpeak), weight='bold', va='top', ha='center', **kw)
    plt.draw()
    plt.savefig('histogram_of_cosmics.pdf')

def test():
    image = Cosmics.cosmicImage(exptime=2.0, size=16, rate=5.0, gradient=False, version='fancy', diffusion=False)
    return image
