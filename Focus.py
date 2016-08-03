"""Keep track of camera focus (to make cartoon drifts)."""
import numpy as np
import matplotlib.pylab as plt


class Focus(object):
    def __init__(self, camera=None, span=(0.0, 10.0)):

        # store the input camera
        self.camera = camera
        if self.camera.variablefocus:
            self.model = self.variablemodel
        else:
            self.model = self.constantmodel

        # what do you want the absolute focus range to be?
        self.span = span
        self.orbit = 13.7
        self.sigma = 0.1
        self.jump = 0.5
        self.best, self.worst = self.span
        self.scale = self.worst - self.best

    def constantmodel(self, counter):
        try:
            counter.shape
            return np.zeros_like(counter)
        except AttributeError:
            return 0.0

    def variablemodel(self, counter):
        bjd = self.camera.counterToBJD(counter)
        phase = (((bjd - self.camera.bjd0) / self.orbit + 0.5) % 1) - 0.5
        best, worst = self.span

        return np.exp(-0.5 * phase ** 2 / self.sigma ** 2) * self.scale + best

    def writeModel(self, outfile):
        time = np.linspace(0, 2 * self.orbit, 1000)
        counter = 24 * 60 * 60 / self.camera.cadence * time
        plt.figure('focus timeseries')
        assert (time.shape == counter.shape)
        plt.plot(time, self.model(counter), linewidth=2)
        plt.xlabel('Time from Observation Start (days)')
        plt.ylabel('Focus (um)')
        plt.xlim(np.min(time), np.max(time))
        plt.ylim(*self.span)
        plt.draw()
        plt.savefig(outfile.replace('.txt', '.pdf'))

        with open(outfile, 'w') as f:
            if self.camera.variablefocus:
                f.write('x = (((bjd - {})/{} + 0.5) % 1) - 0.5\n'.format(self.camera.bjd0, self.orbit))
                f.write('focus = {best} + {scale}*exp(-0.5*(x/{sigma})**2)'.format(**self.__dict__))
            else:
                f.write('fixed to 0.0 microns')


def plothist2d(hist, title=None, log=False,
               xtitle=None, ytitle=None, filename=None):
    """Plot a 2D histogram."""
    map = hist[0]
    x = hist[1][1:] + (hist[1][0] - hist[1][1]) / 2.0
    y = hist[2][1:] + (hist[2][0] - hist[2][1]) / 2.0
    fig = plt.figure(figsize=(10, 10))
    plt.clf()
    plt.subplots_adjust(hspace=0, wspace=0)
    ax_map = fig.add_subplot(2, 2, 3)
    ax_vert = fig.add_subplot(2, 2, 4, sharey=ax_map)
    ax_hori = fig.add_subplot(2, 2, 1, sharex=ax_map)

    ax_hori.plot(x, np.sum(map, 0) / np.sum(map), marker='o', color='black', linewidth=3)
    ax_vert.plot(np.sum(map, 1) / np.sum(map), y, marker='o', color='black', linewidth=3)
    if log:
        ax_vert.semilogx()
        ax_hori.semilogy()
    if log:
        bottom = np.min(map[map > 0]) / np.maximum(np.sum(map, 0).max(), np.sum(map, 1).max())
    else:
        bottom = 0
    top = 1
    ax_hori.set_ylim(bottom, top)
    ax_vert.set_xlim(bottom, top)

    ax_vert.tick_params(labelleft=False)
    ax_hori.tick_params(labelbottom=False)
    if title is not None:
        ax_hori.set_title(title)
    if xtitle is not None:
        ax_map.set_xlabel(xtitle)
    if ytitle is not None:
        ax_map.set_ylabel(ytitle)

    xhalf, yhalf = (x[1] - x[0]) / 2.0, (y[1] - y[0]) / 2.0
    kw = dict(cmap='gray_r',
              extent=[x.min() - xhalf, x.max() + xhalf,
                      y.min() - yhalf, y.max() + yhalf],
              interpolation='nearest')
    if log:
        y = np.log(map)
    else:
        y = map
    ax_map.imshow(y, **kw)
    if filename is not None:
        fig.savefig(filename)
