'''Define, test, and compare potential cosmic ray mitigation strategies.'''
#import matplotlib as mpl
#mpl.rcParams['axes.labelsize'] = 8#'small'
#mpl.rcParams['xtick.labelsize'] = 8#'small'
#mpl.rcParams['ytick.labelsize'] = 8#'small'
#mpl.rcParams['font.size'] = 8#'small'


from ..imports import *
from ..Timeseries import Timeseries
from Strategies import *
plt.ion()
#plt.rc('text', usetex=False)
#plt.rc('font', family='serif')

def compare(nsubexposures=900, niterations=1, remake=False):
	'''Compare performance of cosmic ray rejection strategies on toy model (individual pixel) light curves.'''
	nexposures = 1324*900/nsubexposures
	t = Timeseries(nexposures=nexposures, nsubexposures=nsubexposures, amplitude=1.0)

	subsets = [[mean()],
				[mean(),outlierwithdecay(n=10,threshold=10, memory=0.1),outlierwithdecay(n=10,threshold=10, memory=0.5),outlierwithdecay(n=10,threshold=10, memory=0.9)],
				[mean(), median(n=3), median(n=4),  median( n=5), median( n=6)],
				[mean(),shiftedmedian(n=3), shiftedmedian(n=4), shiftedmedian(n=5), shiftedmedian(n=6)],
				[mean(),lowest(2), lowest(3), lowest(4), lowest(5),  lowest(10), lowest(20), lowest( 60)],
				[mean(),central(3), central(5), central(10),central(20), central( 60)],
				[mean(),  median(n=3), shiftedmedian(n=3), lowest(10),  central(10), outlierwithdecay(n=10,threshold=10, memory=0.5)]]
	labels = ['mean', 'outlierrejection', 'median', 'shiftedmedian', 'lowest', 'central',  'everything']
	for i in range(len(subsets)):
		strategies = subsets[i]
		names = []
		results = {}
		print labels[i]
		while(len(strategies)>0):
			print "  {0:2} strategies to go!".format(len(strategies))
			s = strategies.pop()
			s.test(t, niterations=niterations, remake=remake)
			plt.draw()
			names.append(s.name)
			results[s.name] = {'amplitudes':s.amplitudes, 'noise':s.noise, 'strategy':s}
			del s

		figure = plt.figure(figsize=(7,4), dpi=150)
		gs = plt.matplotlib.gridspec.GridSpec(1,4, wspace=0,hspace=0)

		ax = plt.subplot(gs[0:])
		plt.cla()
		medianlinestyle, lowestlinestyle, centrallinestyle, shiftedlinestyle, outlierstyle = 0,0,0, 0,0
		linestyles = ['-', '--']
		count =0
		for k in names:
			if "Median" in k and "Shifted" not in k:
				color = 'yellow'
				linestyle = medianlinestyle
				medianlinestyle += 1
			if  "Shifted" in k:
				color = 'blue'
				linestyle = shiftedlinestyle
				shiftedlinestyle  += 1

			if "Lowest" in k:
				color = 'orange'
				linestyle = lowestlinestyle
				lowestlinestyle += 1
			if "Central" in k:
				color = 'green'
				linestyle = centrallinestyle
				centrallinestyle += 1
			if "Reject" in k:
				color = 'red'
				linestyle = outlierstyle
				outlierstyle += 1
			print k
			n = np.float(results[k]['strategy'].n)
			scale = 2
			if "Mean" in k:
				color = 'black'
				linestyle=0
				linewidth=5
				alpha=0.2
			else:
				linewidth = (len(names) - linestyle + 1)/scale
				#linewidth = np.log(n/t.nsubexposures) - np.log(1.0/t.nsubexposures)
				alpha=1-(linewidth)/(len(names)+1.0)#(1 - np.log(1.0/t.nsubexposures))

			ax.plot(results[k]['amplitudes'], results[k]['noise'], label=k.replace(' out of ', '/').replace('; ', '\n '), color=color, linewidth=linewidth, alpha=alpha, linestyle=linestyles[linestyle % len(linestyles)])
			ax.legend(bbox_to_anchor=(0.95, 0.95), loc=1, borderaxespad=0., fontsize=7)
			ax.set_xlabel(r'Cosmic Ray Amplitude (in $\sigma$)')
			ax.set_ylabel(r'Achieved Noise (in $\sigma$)')
			ax.set_ylim(0.95, 1.25)
			ax.set_title('Cosmic Ray Rejection when Stacking {0} Exposures'.format(t.nsubexposures))
			count +=1
		plt.tight_layout()
		plt.draw()
		plt.savefig( '/Users/zkbt/Cosmos/Data/TESS/CR/comparisons/cosmicrayrejectioncomparisons_{2:03.0f}iterations_{1:.0f}subexposures_{0}.pdf'.format(labels[i], t.nsubexposures, niterations))
