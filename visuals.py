import Cube

c = Cube.Cube(subject='HD 36486', cadence=120, n=60, size=170, jitter=True)
c.load()
c.display()
'''p = Photometer.Photometer(c)
p.drawApertures(level=5)
m = [a.mag for a in p.apertures]
s = np.random.choice(np.argsort(m)[50], size=13)
apertures = np.array(p.apertures)[s]
i = np.zeros_like(c.master())
for a in apertures:
    i[a.row, a.col] += 1'''

import zachopy.display
ds9 = zachopy.display.ds9('apertures')
ds9.one(c.master(), clobber=True)
