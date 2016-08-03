"""Calculate the expected photometric precision for TESS.

(translated from Josh Winn's IDL TESS signal-to-noise calculator
 on the TESS wiki and updated to include calculations
 published with Peter Sullivan's simulation paper)."""
import pkgutil

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.ascii
import scipy.interpolate
import logging
from settings import log_file_handler

logger = logging.getLogger(__name__)
logger.addHandler(log_file_handler)

import Cartographer


# create a cartographer for managing conversions between ecliptic and galactic coordinates
carto = Cartographer.Cartographer()

# create an interpolator to estimate the best number of pixels in a photometric aperture
optimalpixelsdata = astropy.io.ascii.read(pkgutil.get_data(__name__, 'relations/optimalnumberofpixels.txt'))
optimalpixelsinterpolator = scipy.interpolate.interp1d(optimalpixelsdata['tmag'], optimalpixelsdata['npix'],
                                                       kind='linear', bounds_error=True)


def optimal_npix(imag):
    """Return the average best number of pixels to use in a photometric aperture, given a Cousins I magnitude."""
    return optimalpixelsinterpolator(imag)


# create an interpolator to estimate the average enclosed flux fraction for a aperture size
encloseddata = astropy.io.ascii.read(pkgutil.get_data(__name__, 'relations/enclosedfractionofflux.txt'))
enclosedinterpolator = scipy.interpolate.interp1d(encloseddata['npix'], encloseddata['fraction'])


def enclosed_fraction(npix):
    """Return the average fraction of enclosed energy in a given number of pixels."""
    return enclosedinterpolator(npix)


# create an interpolator to estimate the TESS magnitude for a star,
# given is Cousin I magnitude and effective temperature
fluxesdata = astropy.io.ascii.read(pkgutil.get_data(__name__, 'relations/fluxes.txt'))
fluxesinterpolator = scipy.interpolate.interp1d(fluxesdata['teff'], fluxesdata['icminust'] / 1000.0)


def IctoT(Icmag, teff=5000.0):
    """Return the TESS magnitude of a star, given its Cousin I magnitude and effective temperature."""
    tmag = Icmag - fluxesinterpolator(teff)
    return tmag


def noise(imag=10.0, exptime=1800.0, teff=5000.0,
          elon=0.0, elat=30.0, glon=None, glat=None, ra=None, dec=None,
          subexptime=2.0, npix_aper=4, frac_aper=0.76, e_pix_ro=10.0,
          effective_area=73.0, pix_scale=21.1, sys_limit=60.0,
          verbose=False):
    """Calculate noise, given input Ic magnitude, returing the fractional rms (= 1/snr)

        Mandatory inputs

           imag, $                           apparent mag in Cousins I band

        Optional inputs

           exptime                           total exposure time in seconds
           teff                              effective temperature in Kelvins
           elon, elat                        ecliptic coordinates in degrees
           subexptime                        subexposure time (n_exp = exptime/subexptime)
           npix_aper                         number of pixels in photometric aperture
           frac_aper                         fraction of flux enclosed in photometric aperture
           e_pix_ro                          rms in no. photons/pixel from readout noise
           effective_area                    geometric collecting area
           pix_scale                         arcsec per pixel
           sys_limit                         minimum uncertainty in 1 hr of data, in ppm
           verbose                           request verbose output
    """

    # convert from Cousins I to TESS magnitude
    tmag = IctoT(imag, teff=teff)

    # pick the optimal number pixels in the photometric aperture
    npix_aper = optimal_npix(imag)

    # determine the fraction of the stellar flux included
    frac_aper = enclosed_fraction(npix_aper)

    # solid area of a pixel
    omega_pix = pix_scale ** 2.

    # how many subexposures composed this one exposure?
    n_exposures = exptime / subexptime

    # the TESS zeropint
    tmag0 = 1.514e6

    # photoelectrons from the star
    e_star = 10.0 ** (-0.4 * tmag) * tmag0 * effective_area * exptime * frac_aper

    if ra is not None and dec is not None:
        elon, elat = carto.point(ra, dec, 'celestial').ecliptic.tuple

    logger.debug('imag = {}'.format(imag))
    logger.debug('tmag = {}'.format(tmag))
    logger.debug('tmag0 = {}'.format(tmag0))
    logger.debug('exptime = {}'.format(exptime))
    logger.debug('teff = {}'.format(teff))
    logger.debug('elon = {}'.format(elon))
    logger.debug('elat = {}'.format(elat))
    logger.debug('npix_aper = {}'.format(npix_aper))
    logger.debug('frac_aper = {}'.format(frac_aper))
    logger.debug('subexptime = {}'.format(subexptime))
    logger.debug('n_exposures = {}'.format(n_exposures))
    logger.debug('e_pix_ro = {}'.format(e_pix_ro))
    logger.debug('effective_area = {}'.format(effective_area))
    logger.debug('pix_scale = {}'.format(pix_scale))
    logger.debug('omega_pix = {}'.format(omega_pix))
    logger.debug('sys_limit = {}'.format(sys_limit))
    logger.debug('e_star = {}'.format(e_star))

    # photoelectrons/pixel from zodiacal light
    dlat = (np.abs(elat) - 90.) / 90.
    vmag_zodi = 23.345 - 1.148 * dlat ** 2.
    e_pix_zodi = 10.0 ** (-0.4 * (vmag_zodi - 22.8)) * 2.39e-3 * effective_area * omega_pix * exptime

    logger.debug('vmag_zodi = {}'.format(vmag_zodi))
    logger.debug('e_pix_zodi = {}'.format(e_pix_zodi))

    # photoelectrons/pixel from background stars
    try:
        coord = carto.point(elon, elat, 'ecliptic')
        glon, glat = coord.galactic.tuple
    except:
        glon, glat = 96.36079818, -30.18846954
    glon = np.array([glon])
    glat = np.array([glat])

    logger.debug('glon = {GLON}, glat = {GLAT}'.format(GLON=glon, GLAT=glat))

    dlat = np.abs(glat) / 40.0
    dlon = glon
    q = (dlon > 180.)
    dlon[q] = 360. - dlon[q]
    dlon = np.abs(dlon) / 180.0
    p = [18.9733, 8.833, 4.007, 0.805]
    imag_bgstars = p[0] + p[1] * dlat + p[2] * dlon ** (p[3])
    e_pix_bgstars = 10.0 ** (-0.4 * imag_bgstars) * 1.7e6 * effective_area * omega_pix * exptime


    logger.debug('imag_bgstars = {}'.format(imag_bgstars))
    logger.debug('e_pix_bgstars = {}'.format(e_pix_bgstars))

    noise_star = np.sqrt(e_star) / e_star
    noise_sky = np.sqrt(npix_aper * (e_pix_zodi + e_pix_bgstars)) / e_star
    noise_ro = np.sqrt(npix_aper * n_exposures) * e_pix_ro / e_star
    noise_sys = 0.0 * noise_star + sys_limit / 1e6 / np.sqrt(exptime / 3600.)
    noise = np.sqrt(noise_star ** 2. + noise_sky ** 2. + noise_ro ** 2. + noise_sys ** 2.)

    '''if verbose:
        logger.debug('noise_star [ppm] = ', noise_star*1e6
        logger.debug('noise_sky  [ppm] = ', noise_sky*1e6
        logger.debug('noise_ro   [ppm] = ', noise_ro*1e6
        logger.debug('noise_sys  [ppm] = ', noise_sys*1e6
        logger.debug('noise      [ppm] = ', noise*1e6
        '''
    return noise


def demo(span=27.4, period=12.345678, mean=17, amplitude=1.0):
    """Demonstration of the TESS noise calculator, on a faint and highly-variable star."""

    # create times at a half hour spacing
    t = np.arange(0, span, 0.5 / 24.0)
    n = len(t)

    # create a (perfectly smooth) noiseless model
    noiselessmodel = mean + amplitude * np.sin(2 * np.pi * t / period)

    # calculate the per-point photometric uncertainty for each point of that model
    # noinspection PyTypeChecker
    perpointuncertainty = noise(imag=noiselessmodel)

    # create one random realization of the noise
    noiserealization = np.random.normal(0, 1, n) * perpointuncertainty

    # simulate measurements
    simulated = noiselessmodel + noiserealization

    # create plot showing the demonstration
    plt.ion()
    plt.figure('demonstration', figsize=(10, 3), dpi=200)
    plt.cla()
    # plot the simulated measurements
    plt.errorbar(t, simulated, perpointuncertainty, marker='o', elinewidth=2, linewidth=0, color='black', alpha=0.5)
    # plot the noiseless model
    plt.plot(t, noiselessmodel, color='green', linewidth=2, alpha=0.5)
    # clean up the look of the plot
    # noinspection PyTypeChecker
    plt.ylim(mean + amplitude + np.max(perpointuncertainty) * 5, mean - amplitude - np.min(perpointuncertainty) * 5)
    plt.xlim(np.min(t), np.max(t))
    plt.xlabel('Time (in days)')
    plt.ylabel('Flux (magnitudes)')
    plt.tight_layout()

    return t, simulated, perpointuncertainty
