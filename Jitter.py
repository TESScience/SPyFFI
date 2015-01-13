'''Keep track of spacecraft jitter.'''
from imports import *
import settings

class Jitter(Talker):

    def __init__(self, camera=None):

            # decide whether or not this CCD is chatty
            Talker.__init__(self, mute=False, pithy=False)

            # store the input camera
            self.camera = camera

            # update the jitterball to one that has been binned to this cadence
            self.load()


    def load(self, remake=False):
        '''Load the jitterball for this camera, binned to the appropriate cadence.'''
        self.speak('Populating the jitterball for {0:.0f} second cadence.'.format(self.camera.cadence))

        try:
            # if the jitterball is already loaded and of the correct cadence, we're all set!
            self.jitterball

            # make sure the we're using the right jitterball for this cadence
            assert(self.jittercadence == self.camera.cadence)
        except:
            # if the jitterball isn't already loaded (or has wrong cadence), then load/create a new one!
            jitterfile = settings.prefix + 'intermediates/jitter_{0:04.0f}.npy'.format(self.camera.cadence)
            try:
                assert(remake == False)
                # if a processed jitterfile already exists, then just load it up:
                self.jitterball, self.jittermap = np.load(jitterfile)
                self.jittercadence = self.camera.cadence
            except:
                # otherwise, create a binned jitterfile

                # load simulated jitterball that Roland got from Orbital
                self.speak("Using the raw jitter input file AttErrTimeArcsec_80k.dat (which for convenience may be saved in raw_jitter_file.npy)")
                try:
                  data = np.load(settings.prefix + 'intermediates/raw_jitter_file.npy')
                except:
                  data = astropy.io.ascii.read(settings.prefix + "inputs/AttErrTimeArcsec_80k.dat", names=['t','x', 'y', 'z'])
                  np.save(settings.prefix + 'intermediates/raw_jitter_file.npy', data)

                # subtract means
                data['x'] -= np.mean(data['x'])
                data['y'] -= np.mean(data['y'])
                data['z'] -= np.mean(data['z'])

                # scale jitterball to requirements (should be inflation by ~1.5)
                required_jitter_rms = 2.0/3
                original_fifthsec = np.sqrt(np.mean(data['x']**2 + data['y']**2))
                data['x'] *= required_jitter_rms/original_fifthsec
                data['y'] *= required_jitter_rms/original_fifthsec
                data['z'] *= required_jitter_rms/original_fifthsec

                # smooth them to the required cadence
                self.speak("Smoothing the jitter to {0}s cadence.".format(self.camera.cadence))

                # define convolution filter to smooth timeseries over as many samples as necessary to match the cadence
                spacing = data['t'][1] - data['t'][0]
                n = np.long(self.camera.cadence/spacing)
                filter = np.ones(n)/n

                # construct smoothed timeseries, sampled at native (high) time resolution
                smoothed_t = np.convolve(data['t'], filter, mode='valid')
                smoothed_x = np.convolve(data['x'], filter, mode='valid')
                smoothed_y = np.convolve(data['y'], filter, mode='valid')
                smoothed_z = np.convolve(data['z'], filter, mode='valid')

                # sample smoothed timeseries at the camera's cadence
                t = smoothed_t[::n]
                x = smoothed_x[::n]
                y = smoothed_y[::n]
                z = smoothed_z[::n]

                # plot each dimension separately
                plotfilename =  settings.prefix + 'plots/jitter_timeseries_{0:04.0f}.pdf'.format(self.camera.cadence)
                self.speak("Saving plot of the binned jitter timeseries to " + plotfilename)

                # create the plot
                fi, ax = plt.subplots(3,1, sharey=True, sharex=True)
                ax[0].plot(data['t'], data['x'], alpha=0.5, color='black')
                ax[0].plot(t, x, linewidth=2, alpha=0.5, marker='o', color='red')
                ax[1].plot(data['t'], data['y'], alpha=0.5, color='black')
                ax[1].plot(t, y, linewidth=2, alpha=0.5, marker='o', color='red')
                ax[2].plot(data['t'], data['z'], alpha=0.5, color='black')
                ax[2].plot(t, z, linewidth=2, alpha=0.5, marker='o', color='red')
                ax[0].set_xlim(0,self.camera.cadence*10)
                ax[0].set_title('Expected TESS Pointing Jitter for {0}s Cadence'.format(self.camera.cadence))
                ax[0].set_ylabel('x (")')
                ax[1].set_ylabel('y (")')
                ax[2].set_ylabel('z (")')
                ax[2].set_xlabel('Time (seconds)')
                plt.show()
                fi.savefig(plotfilename)

                # create a 2D jittermap
                narcsec = 3
                bins = np.ceil(narcsec/self.camera.pixelscale/self.camera.psf.subpixsize).astype(np.int)*2 +1	# bins are in units of subpixels
                range = [[-(bins-1)/2*self.psf.subpixsize,(bins-1)/2*self.camera.psf.subpixsize],[-(bins-1)/2*self.camera.psf.subpixsize,(bins-1)/2*self.camera.psf.subpixsize]] # range is in units of pixels

                # make interpolators to keep track of the running smooth means at every moment
                x_interpolator = scipy.interpolate.interp1d(smoothed_t, smoothed_x,'nearest',fill_value=0,bounds_error=False)
                y_interpolator = scipy.interpolate.interp1d(smoothed_t, smoothed_y,'nearest',fill_value=0,bounds_error=False)


                # assign the jittermap here, which will be used for convolution in the PSF code
                self.jittermap = np.histogram2d((data['x'] - x_interpolator(data['t']))/self.camera.pixelscale, (data['y']  - y_interpolator(data['t']))/self.camera.pixelscale, bins=bins, range=range,normed=True)
                self.jitterball = (x,y,z)
                self.jittercadence = self.camera.cadence

                self.speak('Saving plots of the intra-exposure jitter map (which which PSFs should be convolved).')
                # plot an easier to view histogram of the jitterball
                jittermap_to_plot = np.histogram2d((data['x'] - x_interpolator(data['t']))/self.camera.pixelscale, (data['y']  - y_interpolator(data['t']))/self.camera.pixelscale, bins=50, range=range,normed=True)
                self.plothist2d(jittermap_to_plot,title='TESS Pointing Jitter over {0}s'.format(self.camera.cadence),xtitle='Pixels', ytitle='Pixels')
                plt.savefig(settings.prefix + 'plots/jitter_map_{0:04.0f}.pdf'.format(self.camera.cadence))

                # plot the adopted jitterball, as more useful binning
                self.plothist2d(self.jittermap ,title='TESS Pointing Jitter over {0}s'.format(self.camera.cadence),xtitle='Pixels', ytitle='Pixels')
                plt.savefig(settings.prefix + 'plots/jitter_map_adopted_{0:04.0f}.pdf'.format(self.camera.cadence))

                # save the necessary jitter files so we don't have to go through this again
                self.speak('Saving the jitter files for this cadence to {0}'.format(jitterfile))
                np.save(jitterfile,( self.jitterball, self.jittermap))

    def plothist2d(self, hist, title=None, log=False, xtitle=None, ytitle=None):
        '''Plot a 2D histogram. (Used for load plots -- should probably move all this to new Jitterball object.)'''
        map = hist[0]
        x = hist[1][1:] + (hist[1][0] - hist[1][1])/2.0
        y = hist[2][1:]+ (hist[2][0] - hist[2][1])/2.0
        fig = plt.figure(figsize=(10,10))
        plt.clf()
        plt.subplots_adjust(hspace=0, wspace=0)
        ax_map = fig.add_subplot(2,2,3)
        ax_vert = fig.add_subplot(2,2,4, sharey=ax_map)
        ax_hori = fig.add_subplot(2,2,1, sharex=ax_map)

        ax_hori.plot(x, np.sum(map, 0)/np.sum(map), marker='o', color='black', linewidth=3)
        ax_vert.plot(np.sum(map, 1)/np.sum(map), y, marker='o', color='black', linewidth=3)
        if log:
          ax_vert.semilogx()
          ax_hori.semilogy()
        if log:
          bottom = np.min(map[map > 0])/np.maximum(np.sum(map,0).max(),np.sum(map,1).max())
        else:
          bottom = 0
        top = 1
        ax_hori.set_ylim(bottom,top)
        ax_vert.set_xlim(bottom,top)

        ax_vert.tick_params(labelleft=False)
        ax_hori.tick_params(labelbottom=False)
        if title is not None:
          ax_hori.set_title(title)
        if xtitle is not None:
          ax_map.set_xlabel(xtitle)
        if ytitle is not None:
          ax_map.set_ylabel(ytitle)
        if log:
          ax_map.imshow(np.log(map), cmap='gray_r', extent=[x.min(), x.max(), y.min(), y.max()],interpolation='nearest')
        else:
          ax_map.imshow(map, cmap='gray_r', extent=[x.min(), x.max(), y.min(), y.max()],interpolation='nearest')
        plt.draw()

    @property
    def x(self):
        return self.jitterball[0]

    @property
    def y(self):
        return self.jitterball[1]

    @property
    def z(self):
        return self.jitterball[2]

    def jitter(self, counter, dx=None, dy=None, dz=None, header=None, scale=1.0 ):
        '''Jitter the cameras by a little bit, by introducing nudges draw from a (cadence-appropriate) jitterball timeseries.'''

        # make sure the jitterball has been populated
        self.load()

        # assign the nudges in two translations and one rotation
        if dx is None:
            self.camera.nudge['x'] = scale*self.jitterball[0][counter]#*0.002#self.camera.counter*0.0#
        else:
            self.camera.nudge['x'] = dx

        if dy is None:
            self.camera.nudge['y'] = scale*self.jitterball[1][counter]#*0.002#-self.camera.counter*3.0#
        else:
            self.camera.nudge['y'] = dy

        if dz is None:
            self.camera.nudge['z'] = scale*self.jitterball[2][counter]#*0.002#-self.camera.counter*3.0#
        else:
            self.camera.nudge['z'] = dz

        # if possible, write the details to the supplied FITS header
        try:
            header['MOTION'] = ''
            header['MOTNOTE'] = ('', 'properties of the image motion applied')
            header['JITTERX'] = (self.camera.nudge['x'], '["] jitter-induced nudge')
            header['JITTERY'] = (self.camera.nudge['y'], '["] jitter-induced nudge')
        except:
            pass

        # move the camera, using the updated nudge values
        self.speak("Jittering the camera to {x},{y},{z} away from nominal pointing.".format(**self.camera.nudge))
        #self.camera.point()
