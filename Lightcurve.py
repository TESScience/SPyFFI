"""Simulate some toy model light curves."""
import numpy as np
import matplotlib.pyplot as plt
import zachopy.units as u
import logging
from settings import log_file_handler

logger = logging.getLogger(__name__)
logger.addHandler(log_file_handler)


rotation_table = None


def get_rotationtable():
    """Lazily load the McQuillan rotation periods table"""
    import astropy.io.ascii as ascii
    import pkgutil
    global rotation_table
    if rotation_table is None:
        logger.info("Reading McQuillan rotation table")
        rotation_table = ascii.read(pkgutil.get_data(__name__, 'data/rotation_McQuillan.txt'))
    return rotation_table


transit_table = None


def get_transit_table():
    """Lazily load the Kepler TCE table"""
    # from Seader et al.
    # "The search includes a total of $198,646$ targets,
    # of which $112,001$ were observed in every quarter
    # and $86,645$ were observed in a subset of the 17 quarters."
    import astropy.io.ascii as ascii
    import pkgutil
    global transit_table
    if transit_table is None:
        logger.info("Reading Kepler TCE table")
        transit_table = ascii.read(pkgutil.get_data(__name__, 'data/keplerTCE_DR24.txt'))
    return transit_table


def draw_rotation(prng=np.random):
    """Return a random sine curve, drawn from McQuillan rotation periods"""
    row = prng.choice(get_rotationtable())
    # noinspection PyTypeChecker,PyUnresolvedReferences
    return Sinusoid(P=row['PRot'],
                    E=prng.uniform(0, row['PRot']),
                    A=row['Rper'] / 1.0e6)


# noinspection PyUnresolvedReferences
def draw_transit(prng=np.random):
    """Draw a random transit from the Kepler TCE list"""
    row = prng.choice(get_transit_table())
    T14 = row['tce_duration'] / 24.0
    return Trapezoid(P=row['tce_period'],
                     E=row['tce_time0bk'] + 2545833.0,
                     T14=T14,
                     T23=max(T14 - 2 * row['tce_ingress'] / 24.0, 0),
                     D=row['tce_depth'] / 1.0e6)


def parseCode(code):
    """extract name and traits from a lightcurve code"""

    traits = {}
    name, traitstrings = code.split('|')
    for t in traitstrings.split(','):
        k, v = t.split('=')
        traits[k] = np.float(v)
    return name, traits


def generate(code):
    """generate a lightcurve object from a lightcurve code string"""
    name, traits = parseCode(code)
    return globals()[name](**traits)


# TODO: Support Custom
# TODO: This should not have a **kw splat
def random(options=('trapezoid', 'sin'),
           fractionwithextremelc=0.01, fractionwithrotation=None,
           fractionwithtrapezoid=None, fractionwithcustom=0.0, **kw):
    """
    random() returns random Lightcurve.

    random() makes use of these keyword arguments:

        options=['trapezoid', 'sin'] (a list of the kinds of a variability to choose from)
        fractionwithextremelc=0.01 (what fraction should be drawn from a separate "extreme" population of light curves?)
        fractionwithtrapezoid=None (what fraction should get trapezoid? if None, defaults to 18%, from Kepler)
        fractionwithrotation=None (what fraction should get rotation? if None, defaults to 26%, from Kepler)

        when injecting, a star gets (in this order):
            a chance to be extreme
            a chance to be a trapezoid
            a chance to be a rotator
    """

    # first of all, give preference to try to be extreme
    if np.random.uniform(0, 1) < fractionwithextremelc:
        return cartoonrandom(options=options, extreme=True)
    else:
        # by default, set the fraction that get rotation to that from Kepler
        if fractionwithrotation is None:
            fractionwithrotation = 34030.0 / 133030.0
        # by default, set the fraction that get transits to that from Kepler
        if fractionwithtrapezoid is None:
            fractionwithtrapezoid = 20152 / 112001.0

        # if 'custom' in options:
        #     if np.random.uniform(0,1) < fractionwithcustom:
        #         return drawCustom()

        # give trapezoids preference over sin curves
        if 'trapezoid' in options:
            if np.random.uniform(0, 1) < fractionwithtrapezoid:
                return draw_transit()

        # then, try to include a sine curve
        if 'sin' in options:
            if np.random.uniform(0, 1) < fractionwithrotation:
                return draw_rotation()


    # if nothing else, make the light curve a constant
    return constant()


def cartoonrandom(options=('trapezoid', 'sin'), extreme=False, prng=np.random):
    """Generate a random lightcurve, from a cartoonish population"""

    name = prng.choice(options)
    if extreme or name == 'sin':
        if extreme:
            p = [0.1, 30.0]
            a = [0.1, 1]
        else:
            p = [0.1, 30.0]
            a = [0.0001, 0.02]

        P = 10 ** prng.uniform(*np.log10(p))
        E = prng.uniform(0, P)
        A = 10 ** prng.uniform(*np.log10(a))

        # noinspection PyTypeChecker
        return Sinusoid(P=P, E=E, A=A)

    elif name == 'trapezoid':
        if extreme:
            p = [0.1, 30.0]
            d = [0.1, 1]
        else:
            p = [0.1, 30.0]
            d = [0.0001, 0.01]
        P = 10 ** prng.uniform(*np.log10(p))
        E = prng.uniform(0, P)

        mass = prng.uniform(0.1, 1.5)
        radius = mass
        stellar_density = 3 * mass * u.Msun / (4 * np.pi * (radius * u.Rsun) ** 3)
        rsovera = (3 * np.pi / u.G / (P * u.day) ** 2 / stellar_density) ** (1.0 / 3.0)
        T14 = rsovera * P / np.pi
        # noinspection PyTypeChecker
        T23 = prng.uniform(0, T14)
        D = 10 ** prng.uniform(*np.log10(d))

        # noinspection PyTypeChecker
        return Trapezoid(P=P, E=E, D=D, T23=T23, T14=T14)


class LightCurve(object):
    """The LightCurve class defines the basics of a light curve object, which can
        injected into TESS simulations. It handles basic functionality, like
        (importantly), integrating a light curve over a finite exposure time."""

    def __init__(self):
        super(LightCurve, self).__init__()

    def demo(self, tmin=0, tmax=27.4, cadence=30.0 / 60.0 / 24.0, offset=0, raw=False, ax=None):
        """make a plot of a light curve, to show what it looks like"""
        t = np.arange(tmin, tmax, cadence)
        if ax is None:
            plt.figure('demo', figsize=(8, 3))
        else:
            plt.sca(ax)
        y = self.model(t)
        if raw:
            plt.plot(t, y + offset, alpha=0.25, linewidth=1, color='royalblue')
        plt.plot(t, self.integrated(t) + offset, alpha=0.5, linewidth=1, color='darkorange')
        plt.xlim(tmin, tmax)
        # plt.ylim(np.max(y)+0.01, np.min(y)-0.01)
        plt.xlabel('Time (days)')
        plt.ylabel('Flux (mag.)')

    @property
    def code(self):
        """returns a string describing the LightCurve lightcurve"""
        t = ''
        for k, v in self.traits.iteritems():
            t += '{k}={v},'.format(k=k, v=v)

        return '{0}({1})'.format(self.__class__.__name__, t[:-1])

    def integrated(self, t, exptime=30.0 / 60.0 / 24.0, resolution=100):
        """integrate the flux over a finite exposure time.
            for proper averageing, model() must be defined in magnitudes"""

        # deal with the edge case of only a single time point being passed
        try:
            t.shape
        except AttributeError:
            t = np.array([t])

        # don't waste time on this if the light curve is a constant
        if self.__class__.__name__ == 'constant':
            return self.model(np.array(t))

        # create a high-resolution subsampled timeseries
        nudges = np.linspace(-exptime / 2.0, exptime / 2.0, resolution)
        subsampled = t.reshape(1, t.shape[0]) + nudges.reshape(nudges.shape[0], 1)

        # make sure the average is photon-weighted (as opposed to magnitude weighted)
        flux = 10 ** (-0.4 * self.model(subsampled))
        mag = -2.5 * np.log10(flux.mean(0))
        assert (mag.shape == t.shape)
        return mag

    def __repr__(self):
        """string representation of the light curve"""
        return '<{0}>'.format(self.code)


class constant(LightCurve):
    """a "constant" LightCurve light curve is just that, a constant brightness"""

    def __init__(self, **kw):
        self.traits = {}
        LightCurve.__init__(self)

    def model(self, t):
        """model is returned in magnitudes, relative to a baseline level"""
        return np.zeros_like(t)


class Sinusoid(LightCurve):
    """a "sin" LightCurve light curve is just one sine (period, phase, amplitude)"""

    def __init__(self, P=3.1415926, E=0.0, A=0.1, **kw):
        self.traits = dict(P=P, E=E, A=A)
        LightCurve.__init__(self)

    def model(self, t):
        """model is returned in magnitudes, relative to a baseline level"""
        return self.traits['A'] * np.sin(2 * np.pi * (t - self.traits['E']) / self.traits['P'])


class Trapezoid(LightCurve):
    """A "Trapezoid" LightCurve light curve is a generic simplified eclipse model:
                P = period
                E = epoch of one eclipse center
                D = depth
                T14 = duration between contacts 1-4
                T23 = duration between contacts 2-3"""

    def __init__(self, P=3.1415926, E=0.0, D=0.01, T23=0.1, T14=0.1):
        self.traits = dict(P=P, E=E, D=D, T23=T23, T14=T14)
        LightCurve.__init__(self)

    def closesttransit(self, t):
        return np.round((t - self.traits['E']) / self.traits['P']) * self.traits['P'] + self.traits['E']

    def timefrommidtransit(self, t):
        return t - self.closesttransit(t)

    def model(self, t):
        """model is returned in magnitudes, relative to a baseline level"""
        flux = np.zeros_like(t)
        dt = np.abs(self.timefrommidtransit(t))
        start, finish = self.traits['T23'] / 2.0, self.traits['T14'] / 2.0
        depth = self.traits['D']

        i = dt <= start
        flux[i] = depth
        i = (dt <= finish) * (dt > start)
        flux[i] = (finish - dt[i]) / (finish - start) * depth

        return flux

# class custom(LightCurve):
#     def __init__(self, filepath):
#         ''' load a light curve from a file path'''
#         # if you want faster, load the whole light curve into memory here
#
#         # elements of traits will appear in catalog file
#         self.traits = dict(filepath=filepath)
#
#         #
#         self.discretelightcurve = "load"(filename)
#         self.interpolator = "create an interpolator"
#
#     def model(self, t):
#         ''' interpolate a light curve to an infitesimal timepoint '''
#         # if you want less memory intensive, you'll have to keep reloading here
#
#         return self.interpolator(t)
#         # return MAGNITUDES
#
# def drawCustom():
#     '''this function returns a random "custom" light curve object'''
#     good luck!
