from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

W=Path(__file__).resolve().parents[1]; E=W/'results/ensemble'; F=W/'figures'
F.mkdir(parents=True,exist_ok=True)
mg=pd.read_csv(E/'marginal.csv'); sig=np.array([1,3])/np.sqrt(5)
data={}; risks=[]; means=[]; summary=[]; checks=[]
for n in [128,512]:
    freq=np.arange(1,n//2)/n
    for model,r in [(1,0),(2,.8)]:
        ordinary=pd.read_csv(E/f'Ordinary_{n}_{model}.csv').pivot(index='rep',columns='f',values='value').to_numpy()
        data[n,model,'ordinary']=ordinary
        ch=pd.read_csv(E/f'checks_{n}_{model}.csv')
        checks.append({'n':n,'model':model,'fft_error':float(ch.fft.max()),'cross_scale_error':float(ch.cross_scale.max()),'qp_scale_relative_error':float(ch.qp_scale.max())})
        for row in mg.itertuples():
            tau,q,g=row.tau,row.q,row.g
            A=(norm.cdf(q/sig[0])-norm.cdf(q/sig[1]))**2/4
            SC=tau*(1-tau)-A+A*(1-r*r)/(1-2*r*np.cos(2*np.pi*freq)+r*r)
            SU=SC/g**2
            h=np.arange(1,n)
            SUfejer=(tau*(1-tau)+2*A*(np.cos(2*np.pi*freq[:,None]*h)@((1-h/n)*r**h)))/g**2
            data[n,model,tau,'truth']=SU
            data[n,model,tau,'fejer']=SUfejer
            data[n,model,tau,'g']=g
            for method in ['QP','Crossing','Oracle']:
                df=pd.read_csv(E/f'{method}_{n}_{model}.csv')
                P=df[np.isclose(df.tau,tau)].pivot(index='rep',columns='f',values='value').to_numpy()
                if method!='QP': P=P/g**2
                data[n,model,tau,method]=P
                mean=P.mean(axis=0); se=P.std(axis=0,ddof=1)/np.sqrt(len(P))
                for f,m,s,truth,fejer in zip(freq,mean,se,SU,SUfejer):
                    means.append(dict(n=n,model=model,tau=tau,method=method,f=f,mean=m,se=s,truth=truth,oracle_fejer=fejer))
                summary.append(dict(n=n,model=model,tau=tau,method=method,individual_relative_ise=float(np.mean(np.sum((P-SU)**2,axis=1)/np.sum(SU**2))),pooled_relative_ise=float(np.sum((mean-SU)**2)/np.sum(SU**2)),fejer_relative_bias=float(np.sum((SUfejer-SU)**2)/np.sum(SU**2))))
                for R in [1,5,20,100]:
                    Pbar=P.reshape(-1,R,len(freq)).mean(axis=1)
                    ise=np.sum((Pbar-SU)**2,axis=1)/np.sum(SU**2)
                    risks.append(dict(n=n,model=model,tau=tau,method=method,R=R,mean=float(ise.mean()),se=float(ise.std(ddof=1)/np.sqrt(len(ise))),groups=len(ise)))
pd.DataFrame(means).to_csv(E/'mean_spectra.csv',index=False)
pd.DataFrame(risks).to_csv(E/'averaging_risks.csv',index=False)
pd.DataFrame(summary).to_csv(E/'ensemble_summary.csv',index=False)
(E/'checks_summary.json').write_text(json.dumps(checks,indent=2))

plt.rcParams.update({'font.family':'serif','font.serif':['DejaVu Serif'],'mathtext.fontset':'stix','font.size':9,'axes.labelsize':9,'axes.linewidth':.6,'lines.linewidth':1.15,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':7.5,'savefig.bbox':'tight','pdf.fonttype':42})
blue='#2166ac';orange='#b35806';green='#4d9221';red='#b2182b'
n=512;f=np.arange(1,n//2)/n
fig,ax=plt.subplots(1,3,figsize=(7.1,2.55))
for model,color,label in [(1,green,'Independent mixture'),(2,blue,'Persistent volatility')]:
    ax[0].plot(f,data[n,model,'ordinary'].mean(axis=0),color=color,lw=.8,label=label)
ax[0].axhline(1,color='black',ls=':',lw=1.2,label='Population spectrum')
ax[0].set_ylim(.72,1.35);ax[0].set_title('Ordinary periodograms',fontsize=9);ax[0].set_ylabel('Spectrum')
ax[0].legend(frameon=False,loc='upper right',fontsize=6.2)
for a,tau in zip(ax[1:],[.1,.5]):
    g=data[n,2,tau,'g'];v=tau*(1-tau)/g**2
    a.plot(f,data[n,2,tau,'QP'].mean(axis=0)/v,color=blue,label='Mean QP (QSER)')
    a.plot(f,data[n,2,tau,'Crossing'].mean(axis=0)/v,color=orange,ls='--',label='Mean crossing periodogram')
    a.plot(f,data[n,2,tau,'truth']/v,color='black',ls=':',lw=1.3,label='Population shape')
    a.set_title(fr'Persistent volatility, $\tau={tau:g}$',fontsize=9)
    if tau==.1:a.legend(frameon=False,loc='upper right',fontsize=6.2)
    a.set_ylabel('Normalized spectrum');a.set_ylim(0,2.1 if tau==.1 else 1.4)
for a in ax:
    a.set_xlim(0,.5);a.set_xticks([0,.25,.5]);a.set_xlabel(r'$f=\omega/(2\pi)$');a.grid(alpha=.15)
fig.tight_layout(w_pad=1.0);fig.savefig(F/'qser_ensemble_spectra.pdf');plt.close(fig)

fig,ax=plt.subplots(1,2,figsize=(6.6,2.65))
tau=.1;P=data[n,2,tau,'QP'];T=data[n,2,tau,'truth']
for b,color,label in [(1,blue,r'$Y_t$'),(2,red,r'$2Y_t$')]:
    ax[0].plot(f,b*b*P.mean(axis=0),color=color,label=label)
    ax[0].plot(f,b*b*T,color=color,ls=':',lw=1.2)
ax[0].set(xlim=(0,.5),xlabel=r'$f=\omega/(2\pi)$',ylabel='Unnormalized quantile spectrum',title=r'Same crossings, different amplitude')
ax[0].set_ylim(0,4*T.max()*1.07);ax[0].set_xticks([0,.25,.5])
ax[0].legend(frameon=False);ax[0].grid(alpha=.15)
rr=pd.DataFrame(risks)
for tau,color in [(.1,blue),(.5,green),(.9,red)]:
    for method,ls,mk in [('QP','-','o'),('Crossing','--','s')]:
        z=rr[(rr.n==512)&(rr.model==2)&np.isclose(rr.tau,tau)&(rr.method==method)]
        ax[1].plot(z.R,z['mean'],color=color,ls=ls,marker=mk,ms=3,label=fr'{method}, $\tau={tau:g}$')
ax[1].set(xscale='log',yscale='log',xlabel='Records averaged',ylabel='Relative integrated squared error',title='Averaging independent records')
ax[1].set_xticks([1,5,20,100],['1','5','20','100']);ax[1].legend(frameon=False,fontsize=6.3,ncol=2);ax[1].grid(alpha=.15)
fig.tight_layout(w_pad=1.3);fig.savefig(F/'qser_scale_averaging.pdf');plt.close(fig)
