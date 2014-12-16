'''Keep track of catalogs of objects, usually stars.'''
from imports import *
import settings, relations

class Star(object):
	'''a Star object, containing at least RA + Dec + magnitude'''
	def __init__(self, ra=0.0, dec=0.0, tmag=10.0, **kwargs):

		self.coord = astropy.coordinates.ICRS(ra=ra, dec=dec, unit=(astropy.units.deg,astropy.units.deg))
		self.ra = ra
		self.dec = dec
		self.tmag = tmag
		for k in kwargs.keys():
			self.__dict__[k] = kwargs[k]

class Catalog(object):
	'''an object to keep track of lots of stars'''
	def __init__(self):
		pass

	def arrays(self):
		return self.ra, self.dec, self.tmag, self.temperature

	def plot(self):
		plt.ion()
		plt.figure('star chart')
		self.ax = plt.subplot()
		self.ax.cla()
		deltamag = 20.0 - self.tmag
		size = deltamag**2*5
		self.ax.scatter(self.ra, self.dec, s=size, marker='o', color='grey', alpha=0.3, edgecolors='black')
		for i in range(len(self.ra)):
			self.ax.text(self.ra[i], self.dec[i], '{0:.2f}'.format(self.tmag[i]),horizontalalignment='center', verticalalignment='center', alpha=0.5, size=8, color='green',weight='bold')
		self.ax.set_aspect(1)
		self.ax.set_xlabel('Right Ascension')
		self.ax.set_ylabel('Declination')
		plt.show()

class TestPattern(Catalog):
	'''a test pattern catalog, creating a grid of stars to fill an image'''
	def __init__(self, size=3000.0, spacing=150.0, magnitudes=[7,12], ra=0.0, dec=0.0):
		'''create a size x size square (in arcsecs) test pattern of stars,
		with spacing (in arcsecs) between each element and
		magnitudes spanning the range of magnitudes'''

		# how many stars do we need?
		pixels = np.maximum(np.int(size/spacing), 1)
		n = pixels**2

		# construct a linear grid of magnitudes
		self.tmag = np.linspace(np.min(magnitudes), np.max(magnitudes),n)[::-1]
		ras, decs = np.meshgrid(np.arange(pixels)*spacing, np.arange(pixels)*spacing)
		self.dec = ((decs - np.mean(decs))/3600.0 + dec).flatten()
		self.ra = (ras - np.mean(ras)).flatten()/np.cos(self.dec*np.pi/180.0)/3600.0 + ra
		self.temperature = 5800.0 + np.zeros_like(self.ra)

class UCAC4(Catalog):
	def __init__(self, ra=0.0, dec=90.0, radius=0.2, write=True):
		self.load(ra=ra, dec=dec, radius=radius, write=write)

	def load(self, ra=0.0, dec=90.0, radius=0.2, write=True):
		catalog = 'UCAC4'
		ratag = '_RAJ2000'
		dectag = '_DEJ2000'
		if catalog=='UCAC4':
			vcat = 'I/322A/out'
			rmagtag ='f.mag'
			jmagtag = 'Jmag'
			vmagtag = 'Vmag'
			columns = ['_RAJ2000','_DECJ2000','f.mag','Jmag', 'Vmag']
		#if catalog=='Tycho-2':
		#	vcat = 'I/259/tyc2'
		#	rmagtag = 'VTmag'
		#	columns = ['_RAJ2000','_DECJ2000','VTmag']
		v = astroquery.vizier.Vizier(catalog=vcat,columns=columns)
		v.ROW_LIMIT = -1
		bcatalog = 'SIMBAD'
		starsfilename = settings.prefix + 'intermediates/' +  "{catalog}_{ra}_{dec}_{radius}".format(catalog=catalog, ra=ra, dec=dec, radius=radius).replace(' ', '') + '.npy'
		brightstarsfilename = settings.prefix + 'intermediates/' +  "{catalog}_{ra}_{dec}_{radius}".format(catalog=bcatalog, ra=ra, dec=dec, radius=radius).replace(' ', '') + '.npy'
		try:
			t = np.load(starsfilename)
			print "   loading a catalog of stars from ", starsfilename
		except:
			print "   querying {catalog} for ra = {ra}, dec = {dec}, radius = {radius}".format(catalog=catalog, ra=ra, dec=dec, radius=radius)
			t = v.query_region(astropy.coordinates.ICRS(ra=ra, dec=dec, unit=(astropy.units.deg,astropy.units.deg)), radius='{:f}d'.format(radius))[0]
			np.save(starsfilename, t)

		'''try:
			bt = np.load(brightstarsfilename)
			print "      and loading extra bright stars from ", brightstarsfilename

		except:
			bv =Vizier(catalog='I/131A/sao', columns= ['_RAJ2000','_DECJ2000','Vmag'],column_filters={"Vmag":"<3"})
			bv.ROW_LIMIT = -1
			print "   querying bright stars from {catalog} for ra = {ra}, dec = {dec}, radius = {radius}".format(catalog=bcatalog, ra=ra, dec=dec, radius=radius)
			bt = bv.query_region(astropy.coordinates.ICRS(ra=ra, dec=dec, unit=(astropy.units.deg,astropy.units.deg)), radius='{:f}d'.format(radius))[0]
			np.save(brightstarsfilename, bt)'''


		ras = np.array(t[:][ratag])
		decs = np.array(t[:][dectag])
		rmag = np.array(t[:][rmagtag])
		jmag = np.array(t[:][jmagtag])
		vmag = np.array(t[:][vmagtag])

		rbad = (np.isfinite(rmag) == False)*(np.isfinite(vmag))
		rmag[rbad] = vmag[rbad]
		rbad = (np.isfinite(rmag) == False)*(np.isfinite(jmag))
		rmag[rbad] = jmag[rbad]

		jbad = (np.isfinite(jmag) == False)*(np.isfinite(vmag))
		jmag[jbad] = vmag[jbad]
		jbad = (np.isfinite(jmag) == False)*(np.isfinite(rmag))
		jmag[jbad] = rmag[jbad]

		vbad = (np.isfinite(vmag) == False)*(np.isfinite(rmag))
		vmag[vbad] = rmag[vbad]
		vbad = (np.isfinite(vmag) == False)*(np.isfinite(jmag))
		vmag[vbad] = jmag[vbad]



		temperatures = relations.pickles(rmag-jmag)
		imag = rmag - relations.davenport(rmag-jmag)

		'''bad = np.isfinite(imag) == False
		bad = imag < 5
		plt.cla()
		plt.scatter(ras[bad], decs[bad], alpha=0.5, color='blue')
		plt.scatter(bt[:][ratag], bt[:][dectag], s=100/bt[:][vmagtag]**2,color='red',alpha=0.5)


		assert(False)'''
		#assert(np.sum(np.isfinite(imag)==False) == 0)
		ok = np.isfinite(imag)
		print "      found {0} stars with {1} < V < {2}".format(np.sum(ok), np.min(rmag[ok]), np.max(rmag[ok]))

		self.ra = ras[ok]
		self.dec = decs[ok]
		self.tmag = imag[ok]
		self.temperature = temperatures[ok]
		#return ras[ok], decs[ok], rmag[ok], jmag[ok], imag[ok], temperatures[ok]
