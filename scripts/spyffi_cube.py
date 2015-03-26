
import Cube

# plot timeseries of all the pixels
c = Cube.Cube(n=60, size=6, jitter=True, cadence=2, magnitudes=[10], random=True, stacker='Sum')
c.simulate()
c.plot('median')


# display a cube of images in ds9
d = Cube.Cube(n=180, size=32, jitter=True, cadence=2, magnitudes=[10], random=True, stacker='Sum')
d.simulate()
d.display()


# display a cube of images in ds9
e = Cube.Cube(n=180, size=256, jitter=True, cadence=2, magnitudes=[6, 16], random=True, stacker='Sum')
e.simulate(); e.display(limit=180)
