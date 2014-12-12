'''Keep track of what pixels do.'''

def prnu(x,y,type='boxcar'):
  x_withinpixel = x % 1
  y_withinpixel = y % 1
  if type == 'boxcar':
    onedge = (np.abs(x_withinpixel) <= 0.1)+(np.abs(x_withinpixel) >= 0.9)+(np.abs(y_withinpixel) <= 0.1)+(np.abs(y_withinpixel) >= 0.9)
    response = np.ones_like(x)
    response[onedge] *= 0.95
  return response
