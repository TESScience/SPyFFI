from SPyFFI.imports import *
import SPyFFI.settings as settings
from SPyFFI.Cube import Cube
from SPyFFI.Photometer import Photometer
from SPyFFI.Noise import noise
from Strategies import *
import textwrap
import os

def final():
    directory = '/Users/zkbt/Cosmos/Data/TESS/FFIs/outputs/testpattern/120s/sub400x400/cubes/'
    #spacefilename = '/Users/zkbt/Cosmos/Data/TESS/FFIs/outputs/testpattern/120s/sub256x256/cubes/spacemitigated_photometry_400exp_at120s_withjitter_perfectpixels_Rejecting10sigmaoutlierswithchunksof10andmemory0.5.npy'
    spacefilename = directory + 'spacemitigated_photometry_1000exp_at120s_withjitter_perfectpixels_Central8outof10.npy'
    groundfilename = directory +  'groundmitigated_photometry_1000exp_at120s_withoutjitter_perfectpixels_Sum.npy'
    space = np.load(spacefilename)
    ground = np.load(groundfilename)[1]




    pkw = dict(alpha=0.25,marker='o',linewidth=0, markeredgewidth=0, markersize=5)

    scale= 1e6
    bin=0.5
    ax_noise, ax_comp = [], []

    keys = ['Ground Mitigation', 'Space Mitigation']
    colors = ['DarkOrange', 'SlateBlue']
    fi = plt.figure('120 seconds', figsize=(10,16), dpi=75)
    gs = plt.matplotlib.gridspec.GridSpec(2,1,height_ratios=[0.75,1],hspace=0.05,wspace=0.05)
    ax_noise = plt.subplot(gs[0])
    ax_comp = plt.subplot(gs[1])
    ax_noise.set_ylabel('Noise in a 120s Binned Exposure (ppm)')
    ax_comp.set_ylabel('Ratio of Noises')
    ax_comp.set_xlabel('TESS Magnitude')
    ax_noise.set_yscale('log')
    plt.setp(ax_noise.get_xticklabels(), visible=False)
    syslimit = 60.0*np.sqrt(30)
    for j in range(len(keys)):
        if keys[j] == 'Ground Mitigation':
            mag, noises = ground
        elif keys[j] == 'Space Mitigation':
            mag, noises = space
        exp, ach = noises['expected'], noises['achieved']
        # sort stars by magnitude
        i = np.argsort(mag)


        bach = zachopy.oned.binto(mag[i], scale*ach[i], bin)
        bexp = zachopy.oned.binto(mag[i], scale*exp[i], bin)
        x= bexp[0]
        # plot the raw noises
        ax_noise.plot(bach[0], bach[1], color=colors[j], alpha=1, linewidth=3, label=keys[j])
        ax_comp.plot(mag[i], ach[i]/exp[i], color=colors[j], **pkw)
        ax_comp.plot(bexp[0], bach[1]/bexp[1], color=colors[j], alpha=1, linewidth=3)
    ax_noise.plot(bach[0], bexp[1], color='gray', alpha=1, linewidth=3, label='(if no cosmics existed)')


    # plot the comparison
    ax_comp.axhline(1.0, linewidth=3, color='gray')

    ax_comp.set_ylim(0.95, 1.35)#np.min(ach[i]/exp[i]), np.max(ach[i]/exp[i]))
    ax_comp.set_xlim(6, 16)#np.min(ach[i]/exp[i]), np.max(ach[i]/exp[i]))
    inputs = dict(imag=x, exptime=120, ra=0.0, dec=0.0, verbose=False)
    ax_noise.plot(x, np.sqrt(1 + (syslimit/bexp[1])**2)*bexp[1], linewidth=3, linestyle='--', color='gray', alpha=0.5, label='("incorrigible" noise)')
    ax_comp.plot(x, np.sqrt(1 + (syslimit/bexp[1])**2), linewidth=3, linestyle='--', color='gray', alpha=0.5)
    ax_noise.legend(loc='upper left')
    plt.savefig('/Users/zkbt/Dropbox/TESS/memos/cosmics/bestgroundvsspace.pdf')
    print noise(**inputs)/bexp[1]

class Noises(Talker):

    def __init__(self, cadence=120, **kwargs):

        # decide whether or not this CCD is chatty
        Talker.__init__(self, **kwargs)

        self.cadence = cadence

    def create(self, n=100, size=100, remake=False, strategy='Central 8 out of 10'):
        self.ground = Cube(cadence=self.cadence, n=n, size=size, jitter=False)
        # calculate ground-mitigation
        groundfilename = self.ground.filename.replace('cube_', 'groundmitigated_photometry_')
        threshold = 4
        self.descriptions = {}
        line = 40
        self.descriptions['No Mitigation'] = textwrap.fill('doing nothing to correct cosmic rays', line)

        self.descriptions['Ground Mitigation'] = textwrap.fill('rejecting {0}sigma pixel outliers on the ground, assuming stellar variability and jitter are perfectly known'.format(threshold), line)
        try:
            assert(remake==False)
            self.unmitigated, self.groundmitigated = np.load(groundfilename)
            self.speak('loaded from {0}'.format(groundfilename))
        except:
            self.ground.camera.populateCatalog(random=True, magnitudes=[6,16])
            self.ground.load(remake=remake)
            self.phot = Photometer(self.ground)
            self.phot.drawApertures()
            self.unmitigated = self.phot.measure()
            self.phot.cube.oneWeirdTrickToTrimCosmics(threshold=threshold)
            self.groundmitigated = self.phot.measure()
            np.save(groundfilename, (self.unmitigated, self.groundmitigated))
            self.speak('saved to {0}'.format(groundfilename))

        self.descriptions['Space Mitigation'] = textwrap.fill('mitigating cosmic rays onboard, using a [{0}] strategy, including jitter (and uniform intrapixel sensitivity) but not stellar variability'.format(strategy), line)
        self.space = Cube(cadence=self.cadence, n=n, size=size, jitter=True, stacker=strategy)
        spacefilename = self.space.filename.replace('cube_', 'spacemitigated_photometry_')
        try:
            assert(remake==False)
            self.spacemitigated = np.load(spacefilename)
            self.speak('loaded from {0}'.format(spacefilename))
        except:
            self.space.camera.catalog = self.ground.camera.catalog
            self.space.load(remake=remake)
            self.phot = Photometer(self.space)
            self.phot.drawApertures()
            self.spacemitigated = self.phot.measure()
            np.save(spacefilename, (self.spacemitigated))
            self.speak('saved to {0}'.format(spacefilename))

    def plot(self):

        # create figure to compare achieved noises
        fi = plt.figure('{0} seconds'.format(self.cadence), figsize=(15,10), dpi=75)
        fi.clf()
        keys = ['No Mitigation', 'Ground Mitigation', 'Space Mitigation']
        gs = plt.matplotlib.gridspec.GridSpec(2,len(keys),height_ratios=[1,0.3],hspace=0.05,wspace=0.05)
        fi.suptitle('Cosmic Ray Impact for {0}s Exposures'.format(self.cadence), weight='extra bold', size=20)
        ax_noise, ax_comp = [], []
        for j in range(len(keys)):
            if keys[j] == 'No Mitigation':
                mag, noises = self.unmitigated
            elif keys[j] == 'Ground Mitigation':
                mag, noises = self.groundmitigated
            elif keys[j] == 'Space Mitigation':
                mag, noises = self.spacemitigated
            exp, ach = noises['expected'], noises['achieved']
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
            bin=0.5
            bach = zachopy.oned.binto(mag[i], scale*ach[i], bin)
            bexp = zachopy.oned.binto(mag[i], scale*exp[i], bin)
            ax_noise[j].plot(bach[0], bach[1], color='SaddleBrown', alpha=1, linewidth=3)

            x = np.linspace(6,16,100)
            #ax_noise[j].plot(x, scale*noise(imag=x, exptime=self.cadence, ra=self.ground.camera.ra, dec=self.ground.camera.dec, verbose=False), linewidth=3, linestyle='--', color='gray', alpha=0.5)
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
            #ax_comp[j].plot(mag[i], noise(imag=mag[i], exptime=self.cadence,  ra=self.ground.camera.ra, dec=self.ground.camera.dec, verbose=False)/exp[i],  linewidth=3, linestyle='--', color='gray', alpha=0.5)
            ax_comp[j].axhline(1.0, linewidth=3, color='gray')

            ax_comp[j].set_xlabel('TESS magnitude')
            ax_noise[j].set_title(keys[j])
            ax_comp[j].set_ylim(0.8, 1.5)#np.min(ach[i]/exp[i]), np.max(ach[i]/exp[i]))
            ax_comp[j].set_xlim(6, 16)#np.min(ach[i]/exp[i]), np.max(ach[i]/exp[i]))
            ax_noise[j].text(11, np.max(ax_noise[j].get_ylim())*0.5, self.descriptions[keys[j]], va='top', ha='center', alpha=0.6)

        #filename = settings.dirs['plots'] + 'cosmicsin{0}.pdf'.format(self.cadence)
        filename = self.space.filename.replace('cube_','rejectioncomparison_').replace('.npy', '.pdf')
        fi.savefig(filename)
        self.speak('saved plot to {0}'.format(filename))
        os.system('cp {0} ~/Dropbox/TESS/memos/cosmics/.'.format(filename))
