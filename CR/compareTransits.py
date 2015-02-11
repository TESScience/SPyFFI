import transit
from SPyFFI.imports import *
from SPyFFI.Timeseries import Timeseries
from Strategies import *


def compare():
    totest = [ lowest(10),  central(10), outlierwithdecay(n=10,threshold=10, memory=0.5)]
    for s in totest:
        for n in [60,900]:
            for l in ['Hot Jupiter around Sun', 'Habitable Planet around WD']:
                for m in [15]:
                    testTransit(s, mag=m, nsubexposures=n, label=l)


def plottransit(s, label=''):
    plt.figure('transit',figsize=(12, 8), dpi=150)


    plt.clf()
    gs = plt.matplotlib.gridspec.GridSpec(2,2,height_ratios=[1, 0.3], width_ratios=[1,0.2], wspace=0, hspace=0, top=0.85)
    axtransit = plt.subplot(gs[0,0])
    axresiduals = plt.subplot(gs[1,0], sharex=axtransit)
    axhistogram = plt.subplot(gs[1,1], sharey=axresiduals)
    plt.setp(axtransit.get_xticklabels(), visible=False)

    period= s.timeseries.tm.planet.period.value
    t0 = s.timeseries.tm.planet.t0.value
    #
    #t_unbinned = s.unbinned['x']/24.0/60.0/60.0*s.timeseries.exposurecadence
    #flux_unbinned = s.unbinned['flux']
    #ax.plot((t_unbinned % period) - t0, flux_unbinned, linewidth=0, marker='o', color='gray', alpha=0.25, markersize=1)

    t_binned = s.binned['x']/24.0/60.0/60.0*s.timeseries.exposurecadence
    flux_binned = s.binned['flux']
    phasedtime = (t_binned % period) - t0
    sort = np.argsort(phasedtime)
    alphas = dict(naive=0.5, nocosmics=0.5, flux=0.5)
    axhistogram.set_autoscale_on(True)
    noises = {}
    positions = dict(naive = 0.2, flux = 0.5, nocosmics=0.8)
    for tag in ['naive','nocosmics', 'flux']:
        kw = dict(color=s.plotting[tag], markerfacecolor=s.plotting[tag],  markersize=2, marker='o', alpha=alphas[tag]*0.25, markeredgewidth=0, linewidth=0)
        axtransit.plot(phasedtime[sort], s.binned[tag][sort], **kw)
        binwidth = s.timeseries.nsubexposures*2.0/24.0/60.0/60.0/2.0
        bx, by, be = zachopy.oned.binto(phasedtime[sort], s.binned[tag][sort], binwidth)
        axtransit.plot(bx, by, linewidth=2, color=s.plotting[tag], alpha=alphas[tag])
        residuals = s.binned[tag][sort]-s.binned['noiseless'][sort]
        axresiduals.plot(phasedtime[sort], residuals, **kw)
        bx, by, be = zachopy.oned.binto(phasedtime[sort], residuals, binwidth)
        axresiduals.plot(bx, by, linewidth=2, color=s.plotting[tag], alpha=alphas[tag])
        noises[tag] = np.std(residuals)


        s.plothistogram(residuals, ax=axhistogram, binwidth=0.2*s.timeseries.exposurenoise, alpha=0.5, linewidth=2, color=s.plotting[tag])
        left, right = np.log(axhistogram.get_xlim())
        span = right - left
        ylevel = 0

        ywidth = s.timeseries.exposurenoise*7
        axhistogram.text(np.exp(span*positions[tag] + left), ylevel + ywidth*1.2, "{0:.0f}e-\n({1:.2f})\n{2:.0f}".format(noises[tag]*s.timeseries.photonsfromstar*s.timeseries.nsubexposures, noises[tag]/s.timeseries.exposurenoise, noises[tag]*1e6), fontsize=6, color=s.plotting[tag], horizontalalignment='center', alpha=0.7)
        axhistogram.text(np.exp(span*positions[tag] + left), ylevel + ywidth, tag, fontsize=4, color=s.plotting[tag], horizontalalignment='center', alpha=0.7)
        if tag == 'nocosmics':
            scale = np.std(residuals)*5
            axresiduals.set_ylim(-scale, scale)

    axresiduals.set_xlabel('Time from Mid-transit (days)')
    axresiduals.set_ylabel('O-C')
    nsigma = 4
    axtransit.set_ylabel('Relative Flux')
    axtransit.set_title('{0}; {4} subexposures\non a {1:.1f} magnitude star with {2:.0f}% contained flux\nwith a {3}'.format(s.name, s.timeseries.mag, s.timeseries.containedflux*100, label, s.timeseries.nsubexposures))


    axresiduals.set_xlim(-0.05*period, 0.05*period)
    plt.setp(axhistogram.get_yticklabels(), visible=False)
    plt.draw()
    filename = (s.strategyprefix() + '{0:.1f}mag_{1}_{2}exp'.format(s.timeseries.mag, label, s.timeseries.nsubexposures).replace('.','p').replace(' ','')).replace('crdemo','transit') + ".pdf"
    plt.savefig(filename)


def testTransit(s=mean(), mag=10.0, nsubexposures=60, label='Hot Jupiter around Sun'):

    if label == 'Hot Jupiter around Sun':
        rs_over_a = 1.0/10.0
        period = 3.631246
        k = 0.1
    if label == 'Habitable Planet around WD':
        rs_over_a = 1.0/100.0
        period = 0.9351425135
        k = 1.0
    nexposures = 1324*900/nsubexposures

    t = Timeseries(nexposures=nexposures, nsubexposures=nsubexposures, mag=mag)
    t.createTransit(period=period, k=k, rs_over_a=rs_over_a)
    s.calculate(t)
    plottransit(s, label)
