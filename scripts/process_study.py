from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

root = Path(__file__).resolve().parents[1]
data = root / 'results/study'
figures = root / 'figures'
figures.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.serif':['DejaVu Serif'],'mathtext.fontset':'stix','font.size':10,'axes.labelsize':10,'axes.linewidth':.6,'lines.linewidth':1.3,'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':8,'savefig.bbox':'tight','pdf.fonttype':42})
colors = ['#2166ac','#4d9221','#b2182b']
gaussian = pd.read_csv(data/'gaussian_records.csv')
transformation = pd.read_csv(data/'transformation_records.csv')
eeg = pd.read_csv(data/'eeg_windows.csv')
spectra = pd.read_csv(data/'eeg_spectra.csv')
spectra.groupby(['condition','n','tau','f'])[['qser','crossing']].mean().to_csv(data/'eeg_mean_spectra.csv')
for name, frame, groups in [('gaussian',gaussian,['model','n','tau']),('transformation',transformation,['margin','tau']),('eeg',eeg,['condition','n','tau'])]:
    values = [c for c in frame.select_dtypes('number').columns if c not in groups+['rep','window']]
    summary = frame.groupby(groups)[values].agg(['mean','sem','count'])
    summary.columns = ['_'.join(c) for c in summary.columns]
    summary.to_csv(data/(name+'_summary.csv'))

fig, axes = plt.subplots(1,2,figsize=(7.1,2.9),sharey=True)
for ax, model, title in zip(axes,[1,2],['Independent','AR(1), coefficient 0.6']):
    for tau,color in zip([.1,.5,.9],colors):
        d=gaussian[(gaussian.model==model)&(gaussian.tau==tau)].groupby('n').oracle.agg(['mean','sem'])
        ax.errorbar(d.index,d['mean'],yerr=1.96*d['sem'],marker='o',ms=4,capsize=2,color=color,label=fr'$\tau={tau:g}$')
    ax.set(xscale='log',yscale='log',xlabel='Record length',title=title)
    ax.set_xticks([128,256,512],['128','256','512']);ax.xaxis.set_minor_formatter(NullFormatter());ax.yaxis.set_minor_formatter(NullFormatter());ax.set_yticks([.1,.2,.3,.4],['0.1','0.2','0.3','0.4']);ax.grid(alpha=.15)
axes[0].set_ylabel('Relative reconstruction loss');axes[1].legend(frameon=False)
fig.tight_layout();fig.savefig(figures/'qser_approximation.pdf');plt.close(fig)

fig, axes = plt.subplots(1,2,figsize=(7.1,2.9))
for margin,color,label in zip([1,2,3],colors,['Normal','Lognormal',r'Student $t_3$']):
    d=transformation[transformation.margin==margin].groupby('tau').mean(numeric_only=True)
    for ax,value in zip(axes,['variance','loss']):ax.plot(d.index,d[value],marker='o',ms=4,color=color,label=label)
axes[0].set(yscale='log',ylabel='QSER variance',title='Marginal amplitude')
axes[1].set(ylabel='Relative reconstruction loss',title='Density-calibrated comparison')
for ax in axes:ax.set_xlabel(r'Quantile level $\tau$');ax.set_xticks([.1,.5,.9]);ax.grid(alpha=.15)
axes[0].legend(frameon=False);fig.tight_layout(w_pad=2);fig.savefig(figures/'qser_transformation.pdf');plt.close(fig)

fig, axes = plt.subplots(2,3,figsize=(7.1,4.4),sharex=True,sharey=True)
for row,condition in enumerate(['open','closed']):
    for col,tau in enumerate([.1,.5,.9]):
        ax=axes[row,col]
        d=spectra[(spectra.condition==condition)&(spectra.tau==tau)].groupby('f')[['qser','crossing']].mean()
        ax.axvspan(8,13,color='#d4d7db',alpha=.5)
        ax.plot(160*d.index,d.qser,color=colors[0],label='QSER')
        ax.plot(160*d.index,d.crossing,color='#b35806',ls='--',label='Crossing')
        ax.set(xlim=(0,40),ylim=(0,20.5));ax.set_xticks([0,10,20,30,40]);ax.grid(alpha=.15)
        if row==0:ax.set_title(fr'$\tau={tau:g}$')
        if col==0:ax.set_ylabel(('Eyes open' if row==0 else 'Eyes closed')+'\nNormalized spectrum')
        if row==1:ax.set_xlabel('Frequency (Hz)')
axes[0,2].legend(frameon=False);fig.tight_layout(w_pad=.8,h_pad=1);fig.savefig(figures/'qser_eeg.pdf');plt.close(fig)
