import seaborn as sns
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.font_manager as font_manager
import marsilea as ma
import marsilea.plotter as mp
# plt.rcParams['figure.dpi'] = 150
font_dirs = ["/public/home/mowp/test/fonts/"]
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
for font_file in font_files:
    font_manager.fontManager.addfont(font_file)
plt.rcParams["font.family"] = "Arial"
sns.despine(top=True, right=True)
from itertools import product
from functools import reduce
import patchworklib as pw
from jpy_tools.otherTools import pwRecoverSeaborn, pwStack, F, pwShow
pw.show = pwShow
fc_recoverSns = pwRecoverSeaborn()
import seaborn.objects as so
from cycler import cycler
dt_snsStyle = {**sns.axes_style('ticks'), "legend.frameon": False, "axes.spines.top":False, "axes.spines.right":False, "axes.prop_cycle":cycler(color=['grey'])}
import pandas as pd
import numpy as np
from typing import Literal
import sklearn


def heatmapMultiGeneInDifferentTime(ad, layer, ls_genes, ls_names, cbLabel, sort=False, scale=False, cbSide='right', width=3.3, height=6, cmap='Reds',**kwargs):

    _df = ad[:, ls_genes].to_df(layer).T
    _df.index = ls_names
    if sort:
        _df = _df.assign(temp=_df.idxmax(1)).sort_values('temp').drop(columns='temp')
    if scale:
        _ar = sklearn.preprocessing.minmax_scale(_df, axis=1)
    else:
        _ar = _df.values

    _df = pd.DataFrame(_ar, index=_df.index, columns=_df.columns)

    # import pdb; pdb.set_trace()
    h1 = ma.Heatmap(
        _ar, cmap=cmap, label=cbLabel, width=width, height=height, **kwargs
    )
    h1.hsplit(labels=_df.index, order=_df.index, spacing=0.005)
    # h1.add_dendrogram('right')
    h1.add_left(ma.plotter.Labels(_df.index, fontstyle='italic'), pad=0.1)
    h1.add_top(ma.plotter.Labels(_df.columns, rotation=0, align="center"), pad=0.01)
    h1.add_legends(
        pad=0.1,
    )
    h1.add_top(ma.plotter.Title("CT (Hour)"))
    # h1.add_dendrogram('left', method='ward')

    h1.add_legends(side=cbSide)
    h1.render()
    return h1

def heatmapSingleGeneMultiCluster(ad, layer, gene, geneName, clusterKey, cbLabel, cbSide, sortby: Literal["cluster", "peak", "phase"] = "phase", scale=False, cmap='Reds', inputIsInObs=False,**kwargs):
    def getPhase(x, y):
        def _cos(x, a, b, c):
            return a * np.cos((np.pi / 12) * x - (np.pi / 12) * b) + c

        from scipy.optimize import curve_fit

        popt, pcov = curve_fit(_cos, x, y, bounds=([0, -12, 0], [np.inf, 12, np.inf]))
        std = np.sqrt(np.diag(pcov))[1]
        return popt[1]

    if not inputIsInObs:
        _df = (
            ad[:, gene]
            .to_df(layer)
            .merge(ad.obs[["CT", clusterKey]], left_index=True, right_index=True)
            .pivot_table(index=clusterKey, columns="CT", values=gene)
        )
    else:
        _df = ad.obs[[gene, "CT", clusterKey]].pivot_table(index=clusterKey, columns="CT", values=gene)

    if scale:
        _ar = sklearn.preprocessing.minmax_scale(_df, axis=1)
        _df = pd.DataFrame(_ar, index=_df.index, columns=_df.columns)

    if sortby == "cluster":
        _df = _df.sort_index(key=lambda x: x.str.split(',').str[0].astype(int))
    elif sortby == "peak":
        _df = _df.loc[_df.idxmax(axis=1).sort_values().index]
    elif sortby == "phase":
        _df = (
            _df.assign(
                phase=[
                    getPhase(np.linspace(0, 22, 12), x[1].values)
                    for x in _df.iterrows()
                ]
            )
            .sort_values("phase")
            .drop(columns=["phase"])
        )

    h = ma.Heatmap(
            _df.values, cmap=cmap, label=cbLabel, width=3.3, height=6, **kwargs
        )
    # h.hsplit(labels=_df.index, order=_df.index, spacing=0.005)

    h.add_top(ma.plotter.Labels(_df.columns, rotation=0, align="center"), pad=0.01)
    h.add_legends(
        pad=0.1,
    )
    h.add_top(ma.plotter.Title("CT (Hour)", fontsize=10))
    h.add_top(ma.plotter.Title(geneName, fontsize=16))
    h.add_left(ma.plotter.Labels(_df.index), pad=0.1, name='clusterLabel')
    h.add_legends(side=cbSide)
    h.render()
    return h

class SoLineForMa(mp.base.StatsBase):
    def __init__(self, data, df, label, gene, name):
        self.set_data(data)
        self.df = df
        self.label = label
        self.gene = gene
        self.name = name

    def render_ax(self, spec):
        ax = spec.ax
        data = self.df
        gene = self.gene
        name = self.name
        g = (
                so.Plot(data.query("gene==@gene"), x="CT", y="CPM")
                .add(so.Dot(color="#4C72B0"))
                .add(so.Line(color="#4C72B0"))
                .theme(dt_snsStyle)
                .label(color="", y=f"{name}\n(log$_{2}$CPM)", x="")
                .on(ax)
                .plot()
            )
        # ax.get_xaxis().set_visible(True)
        g._repr_png_()
        # plt.close()

import matplotlib.path as mpath
star = mpath.Path.unit_regular_star(6)
class SoMarkForJtkAnno(mp.base.StatsBase):
    def __init__(self, data, label, ls_circadianCluster, marksize=None):
        self.set_data(data)
        self.ls_circadianCluster = ls_circadianCluster
        # self.ls_orderedAllCluster = ls_orderedAllCluster
        self.label = label
        self.marksize=marksize

    def render_ax(self, spec):
        ax = spec.ax
        ls_circadianCluster = self.ls_circadianCluster
        ls_orderedAllCluster = spec.data.tolist()[0]
        for i, cluster in enumerate(ls_orderedAllCluster):
            if cluster in ls_circadianCluster:
                ax.scatter(0, i, marker=star, c='#C2383A', s=self.marksize)
        print(ax.get_yticklabels())        
        ax.set_ylim(len(ls_orderedAllCluster)-0.5, -0.5)
        ax.axis("off")