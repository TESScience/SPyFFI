'''Simulate some toy model light curves.'''

from imports import *
import zachopy.units as u

# kludge to figure out path to data resources
import os
datadir = os.path.abspath(os.path.join(codedir, 'data'))

rotationtable = None
transittable = None

# read the McQuillan table (once)
def drawRotation():
    global rotationtable
    if rotationtable is None:
        print "[lightcurve.py] reading McQuillan rotation table"
        rotationtable = astropy.io.ascii.read(datadir+'/rotation_McQuillan.txt')
    row = rotationtable[np.random.randint(0, len(rotationtable))]
    P=row['PRot']
    E=np.random.uniform(0,P)
    A=row['Rper']/1.0e6
    return sin(**locals())


# read the Kepler TCE table (once)
def drawTransit():
    #from Seader et al. "The search includes a total of $198,646$ targets, of which $112,001$ were observed in every quarter and $86,645$ were observed in a subset of the 17 quarters.""
    global transittable
    if transittable is None:
        print "[lightcurve.py] reading Kepler TCE table"
        transittable = astropy.io.ascii.read(datadir+'/keplerTCE_DR24.txt')
        transittable = transittable[transittable['tce_depth']>1.0]

    row = transittable[np.random.randint(0, len(transittable))]
    P=row['tce_period']
    E=row['tce_time0bk']+2545833.0
    T14=row['tce_duration']/24.0
    T23=T14 - 2*row['tce_ingress']/24.0
    T23=np.maximum(T23, 0)
    D=row['tce_depth']/1.0e6
    return trapezoid(**locals())

def parseCode(code):
    '''extract name and traits from a lightcurve code'''

    traits = {}
    name, traitstrings = code.split('|')
    for t in traitstrings.split(','):
        k, v = t.split('=')
        traits[k] = np.float(v)
    return name, traits

def generate(code):
    name, traits = parseCode(code)
    return globals()[name](**traits)

def random(options=['trapezoid', 'sin'], extreme=False, **kw):
    '''
    random() returns random Lightcurve.

    random() makes use of these keyword arguments:
        options=['trapezoid', 'sin'] (a list of the kinds of a variability to choose from)
        extreme=False (should we allow extreme variability [good for movies] or no?)
    '''

    if extreme:
        return cartoonrandom(options=options, extreme=extreme)

    fractionrotators = len(rotationtable)/133030.0
    fractiontransiting = len(transittable)/112001.0

    if 'trapezoid' in options:
        if np.random.uniform(0,1) < fractiontransiting:
            return drawTransit()

    if 'sin' in options:
        if np.random.uniform(0,1) < fractionrotators:
            return drawRotation()

    return constant()

def cartoonrandom(options=['trapezoid', 'sin'], extreme=False):
    if extreme:
        opens = ['sin']
    name = np.random.choice(options)
    if name == 'sin':
        if extreme:
            p = [0.001, 0.1]
            a = [0.1, 1]
        else:
            p = [0.1, 30.0]
            a = [0.0001, 0.02]

        P=10**np.random.uniform(*np.log10(p))
        E=np.random.uniform(0,P)
        A=10**np.random.uniform(*np.log10(a))
        return sin(**locals())

    if name == 'trapezoid':
        if extreme:
            p = [0.001, 0.3]
            d = [0.1, 1]
        else:
            p = [0.1, 30.0]
            d = [0.0001, 0.01]
        P=10**np.random.uniform(*np.log10(p))
        E=np.random.uniform(0,P)

        mass = np.random.uniform(0.1, 1.5)
        radius = mass
        stellar_density = 3*mass*u.Msun/(4*np.pi*(radius*u.Rsun)**3)
        rsovera = (3*np.pi/u.G/(P*u.day)**2/stellar_density)**(1.0/3.0)
        T14 = rsovera*P/np.pi
        T23=np.random.uniform(0, T14)
        D=10**np.random.uniform(*np.log10(d))

        return trapezoid(**locals())

def test(n=100,step=0.01, **kw):
    plt.cla()
    f = plt.figure(figsize=(5,10))
    ax = plt.subplot()
    o =0
    for i in range(n):
        lc = random(**kw)
        lc.demo(offset=o, ax=ax)
        o += step

    ax.set_ylim(step*n, 0)
    plt.draw()

class Cartoon(Talker):
    def __init__(self):
        Talker.__init__(self)

    def demo(self, tmin=0, tmax=27.4, cadence=30.0/60.0/24.0, offset=0, raw=False, ax=None):
        t = np.arange(tmin, tmax, cadence)
        if ax is None:
            plt.figure('demo', figsize=(8,3))
        else:
            plt.sca(ax)
        y = self.model(t)
        if raw:
            plt.plot(t, y+offset, alpha=0.25, linewidth=1, color='royalblue')
        plt.plot(t, self.integrated(t)+offset, alpha=0.5, linewidth=1, color='darkorange')
        plt.xlim(tmin,tmax)
        #plt.ylim(np.max(y)+0.01, np.min(y)-0.01)
        plt.xlabel('Time (days)')
        plt.ylabel('Flux (mag.)')

    @property
    def code(self):
        '''returns a string describing the cartoon lightcurve'''
        t = ''
        for k, v in self.traits.iteritems():
            t += '{k}={v},'.format(k=k,v=v)

        return '{0}({1})'.format(self.__class__.__name__, t[:-1])

    def integrated(self, t, exptime=30.0/60.0/24.0, resolution=100):

        # don't waste time on this if the light curve is a constant
        if self.__class__.__name__ == 'constant':
            return self.model(t)

        # deal with the edge case of only a single time point being passed
        try:
            t.shape
        except AttributeError:
            t = np.array([t])

        # create a high-resolution subsampled timeseries
        nudges = np.linspace(-exptime/2.0, exptime/2.0, resolution)
        subsampled = t.reshape(1, t.shape[0]) + nudges.reshape(nudges.shape[0], 1)

        # make sure the average is photon-weighted (as opposed to magnitude weighted)
        flux = 10**(-0.4*self.model(subsampled))
        mag = -2.5*np.log10(flux.mean(0))
        assert(mag.shape == t.shape)
        return mag

    def __repr__(self):
        return '<{0}>'.format(self.code)


class constant(Cartoon):
    def __init__(self, **kw):
        self.traits = {}
        Cartoon.__init__(self)

    def model(self, t):
        return np.zeros_like(t)

class sin(Cartoon):
    def __init__(self, P=3.1415926, E=0.0, A=0.1, **kw):
        self.traits = dict(P=P, E=E, A=A)
        Cartoon.__init__(self)

    def model(self,t):
        return self.traits['A']*np.sin(2*np.pi*(t - self.traits['E'])/self.traits['P'])

class trapezoid(Cartoon):
    def __init__(self, P=3.1415926, E=0.0, D=0.01, T23=0.1, T14=0.1, **kw):
        self.traits = dict(P=P, E=E, D=D, T23=T23, T14=T14)
        Cartoon.__init__(self)

    def closesttransit(self, t):
        return np.round((t-self.traits['E'])/self.traits['P'])*self.traits['P'] + self.traits['E']

    def timefrommidtransit(self, t):
        return t-self.closesttransit(t)

    def model(self,t):
        flux = np.zeros_like(t)
        dt = np.abs(self.timefrommidtransit(t))
        start, finish = self.traits['T23']/2.0, self.traits['T14']/2.0
        depth = self.traits['D']

        i = dt <= start
        flux[i] = depth
        i = (dt<=finish)*(dt>start)
        flux[i] = (finish-dt[i])/(finish-start)*depth

        return flux
