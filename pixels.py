'''Keep track of what pixels do.'''
from imports import *

def prnu(x,y,type='boxcar'):
    x_distancefrompixeledge = (x - 0.5) % 1.0
    y_distancefrompixeledge = (y - 0.5) % 1.0
    if type == 'boxcar':
        onedge = (x_distancefrompixeledge <= 0.1)|(x_distancefrompixeledge >=0.9)|(y_distancefrompixeledge <= 0.1)|(y_distancefrompixeledge >=0.9)
        response = np.ones_like(x)
        response[onedge] *= 0.95
    return response
