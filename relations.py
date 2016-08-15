'''Handy relations, mostly interpolations out of paper tables.'''
import os
import settings
import astropy.io.ascii
import scipy.interpolate
import numpy as np
import numpy.polynomial.polynomial as polynomial


# convert colors, using the Sloan stellar locus
def davenport(color, input='r-J', output='r-i'):
    data = astropy.io.ascii.read(os.path.join(settings.inputs, 'davenport_table1.txt'),
                                 names=['g-i', '#', 'u-g', 'eu-g', 'g-r', 'eg-r', 'r-i', 'er-i', 'i-z', 'ei-z', 'z-J',
                                        'ez-J', 'J-H', 'eJ-H', 'H-K', 'eH-K'])
    if input == 'r-J':
        x = data['r-i'] + data['i-z'] + data['z-J']
    if output == 'r-i':
        y = data['r-i']
    # plt.figure()
    # plt.scatter(x, y)
    interpolator = scipy.interpolate.interp1d(x, y, 'linear', fill_value=0, bounds_error=False)
    new = interpolator(color)
    new[color > np.max(x)] = y[np.argmax(x)]
    new[color < np.min(x)] = y[np.argmin(x)]
    return new


# spit out the Pickles temperature associated with an input color
def pickles(color, input='R-J'):
    # assuming
    data = astropy.io.ascii.read(os.path.join(settings.inputs, 'pickles_table2.txt'))
    ok = data['[Fe/H]'] == 0
    if input == 'R-J':
        x = data['V-J'][ok] - data['V-R'][ok]
    y = data['logT'][ok]
    # plt.figure()
    # plt.scatter(x, 10**y)

    coeff, stats = polynomial.polyfit(x, y, 4, full=True)
    interpolator = polynomial.Polynomial(coeff)

    # interpolator = scipy.interpolate.interp1d(x,y,'linear', bounds_error=False)
    new = interpolator(color)
    new[color > np.max(x)] = y[np.argmax(x)]
    new[color < np.min(x)] = y[np.argmin(x)]
    return 10 ** new
