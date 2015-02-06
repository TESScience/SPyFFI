from SPyFFI.imports import *
import SPyFFI.settings as settings
from SPyFFI.Cube import Cube
from SPyFFI.Photometer import Photometer
from SPyFFI.Noise import noise
from Strategies import *

class Noises(Talker):

    def __init__(self, cadence=120, **kwargs):

        # decide whether or not this CCD is chatty
        Talker.__init__(self, **kwargs)

        self.cadence = cadence

    def create(self, n=100, size=100, jitter=False, remake=False):
        self.cube = Cube(cadence=self.cadence, n=n, size=size, jitter=jitter)
        filename = self.cube.directory + 'photometryon{0}s{1}px{2}im.npy'.format(self.cadence, size, n)
        try:
            self.unmitigated, self.mitigated = np.load(filename)
            self.speak('loaded from {0}'.format(filename))
        except:
            self.cube.camera.populateCatalog(random=True, magnitudes=[6,16])
            self.cube.load(remake=remake)
            self.phot = Photometer(self.cube)
            self.phot.drawApertures()
            self.unmitigated = self.phot.measure()
            self.phot.cube.oneWeirdTrickToTrimCosmics(threshold=4)
            self.mitigated = self.phot.measure()
            np.save(filename, (self.unmitigated, self.mitigated))
            self.speak('saved to {0}'.format(filename))

    def plot(self):

        # create figure to compare achieved noises
        fi = plt.figure('{0} seconds'.format(self.cadence), figsize=(10,8), dpi=100)
        keys = ['No Mitigation', 'Ground Mitigation']
        gs = plt.matplotlib.gridspec.GridSpec(2,len(keys),height_ratios=[1,0.3],hspace=0.05,wspace=0.05)
        fi.suptitle('Cosmic Ray Impact for {0}s Exposures'.format(self.cadence), weight='extra bold')
        ax_noise, ax_comp = [], []
        for j in range(len(keys)):
            if keys[j] == 'No Mitigation':
                mag, exp, ach = self.unmitigated
            elif keys[j] == 'Ground Mitigation':
                mag, exp, ach = self.mitigated

            # sort stars by magnitude
            i = np.argsort(mag)

            # define axes for this subset
            try:
                noisey = ax_noise[0]
                compy = ax_comp[0]
                sharex = ax_noise[0]
            except:
                noisey = None
                compy = None
                sharex = None

            ax_noise.append(plt.subplot(gs[0,j], sharey=noisey, sharex=sharex))
            ax_comp.append(plt.subplot(gs[1,j], sharex=ax_noise[0], sharey=compy))


            scale= 1e6
            pkw = dict(alpha=0.25,marker='o',linewidth=0, markeredgewidth=0, markersize=5)
            ax_noise[j].plot(mag[i], scale*exp[i], color='black',**pkw)
            ax_noise[j].plot(mag[i], scale*ach[i], color='Sienna', **pkw)
            bin=0.2
            bach = zachopy.oned.binto(mag[i], scale*ach[i], bin)
            bexp = zachopy.oned.binto(mag[i], scale*exp[i], bin)
            ax_noise[j].plot(bach[0], bach[1], color='SaddleBrown', alpha=1, linewidth=3)

            x = np.linspace(6,16,100)
            #ax_noise[j].plot(x, scale*noise(imag=x, exptime=self.cadence, ra=self.cube.camera.ra, dec=self.cube.camera.dec, verbose=False), linewidth=3, linestyle='--', color='gray', alpha=0.5)
            ax_noise[j].set_yscale('log')
            plt.setp(ax_noise[j].get_xticklabels(), visible=False)
            if j > 0:
                plt.setp(ax_noise[j].get_yticklabels(), visible=False)
                plt.setp(ax_comp[j].get_yticklabels(), visible=False)
            else:
                ax_noise[j].set_ylabel('Noise in a {0}s Binned Exposure (ppm)'.format(self.cadence))
                ax_comp[j].set_ylabel('Ratio of\nNoises')#np.min(ach[i]/exp[i]), np.max(ach[i]/exp[i]))

            ax_comp[j].plot(mag[i], ach[i]/exp[i], color='Sienna', **pkw)
            ax_comp[j].plot(bexp[0], bach[1]/bexp[1], color='SaddleBrown', alpha=1, linewidth=3)
            #ax_comp[j].plot(mag[i], noise(imag=mag[i], exptime=self.cadence,  ra=self.cube.camera.ra, dec=self.cube.camera.dec, verbose=False)/exp[i],  linewidth=3, linestyle='--', color='gray', alpha=0.5)
            ax_comp[j].axhline(1.0, linewidth=3, color='gray')

            ax_comp[j].set_xlabel('TESS magnitude')
            ax_noise[j].set_title(keys[j])
            ax_comp[j].set_ylim(np.min(ach[i]/exp[i]), np.max(ach[i]/exp[i]))
            ax_comp[j].set_xlim(6, 16)#np.min(ach[i]/exp[i]), np.max(ach[i]/exp[i]))

        filename = settings.dirs['plots'] + 'cosmicsin{0}.pdf'.format(self.cadence)
        fi.savefig(filename)
        self.speak('saved plot to {0}'.format(filename))
