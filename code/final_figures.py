"""Regenerate final-review comparison data and figures from recorded experiments."""
from pathlib import Path
import json, csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import beta
from certificates import AnytimeCertificate
from monitoring_audit import cg23_beta
ROOT=Path(__file__).resolve().parents[1]

def save(fig,name):
    fig.tight_layout(pad=.45)
    fig.savefig(ROOT/'figures'/f'{name}.pdf',bbox_inches='tight')
    fig.savefig(ROOT/'figures'/f'{name}.png',dpi=220,bbox_inches='tight')
    plt.close(fig)

def main():
    rows=[];cert=AnytimeCertificate(.05)
    for k in [2,5,23,100,1000]:
        n=100000;q=float(beta.ppf(.95,k,n-k+1));cg=cg23_beta(n,k,.05);u=cert(n,k)
        rows.append(dict(n=n,k=k,beta=q,CG23=cg,anytime=u,n_beta=n*q,n_CG23=n*cg,n_anytime=n*u,any_CG=u/cg,any_beta=u/q))
    pd.DataFrame(rows).to_csv(ROOT/'data/final_fixed_time_comparison.csv',index=False)
    with open(ROOT/'tables/final_comparison.tex','w') as f:
        for r in rows:f.write(f"{r['k']} & {r['n_beta']:.2f} & {r['n_CG23']:.2f} & {r['n_anytime']:.2f} & {r['any_CG']:.2f} & {r['any_beta']:.2f}\\\\\n")
    x=pd.read_csv(ROOT/'data/recrossing_example.csv');n=x.n.to_numpy();current=x.U.to_numpy();arch=np.minimum.accumulate(current)
    fig,ax=plt.subplots(figsize=(3.38,2.32));ax.plot(n,current,label='Current controller')
    ax.plot(n,arch,'--',label='Best-certified stored controller');ax.axhline(.05,ls=':',label='Target 0.05')
    ax.set(xlim=(245,395),ylim=(.036,.073),xlabel='Design sample size $n$',ylabel='Risk certificate')
    ax.axvline(351,ls=':',lw=.6);ax.axvline(362,ls=':',lw=.6)
    ax.annotate('Latest controller unqualified\n' + r'$351\leq n\leq361$',xy=(356,.051),xytext=(310,.066),fontsize=7,arrowprops={'arrowstyle':'->','lw':.6})
    ax.legend(fontsize=6.8,loc='lower left');ax.tick_params(labelsize=8);ax.xaxis.label.set_size(8.5);ax.yaxis.label.set_size(8.5)
    save(fig,'certificate_archive')
    m=json.loads((ROOT/'data/final_monitoring_results.json').read_text())
    fig,ax=plt.subplots(figsize=(3.38,2.3))
    for method in ['CG23','Anytime']:
        rr=[r for r in m['coverage'] if r['method']==method];xx=np.array([r['horizon'] for r in rr]);yy=100*np.array([r['frequency'] for r in rr]);er=np.array([[100*(r['frequency']-r['CP95'][0]) for r in rr],[100*(r['CP95'][1]-r['frequency']) for r in rr]])
        ax.errorbar(xx,yy,yerr=er,marker='o',capsize=3,label=method)
    ax.axhline(5,ls=':',label='Nominal 5%');ax.set(xscale='log',xlabel='Monitoring horizon',ylabel='Paths with a bound crossing (%)',ylim=(-.4,8.3));ax.set_xticks([800,8000],labels=['800','8000']);ax.legend(fontsize=7);ax.tick_params(labelsize=8);ax.xaxis.label.set_size(8.5);ax.yaxis.label.set_size(8.5)
    save(fig,'monitoring_coverage')
    fig,ax=plt.subplots(figsize=(3.38,2.3));ks=[r['k'] for r in rows]
    ax.plot(ks,[r['any_beta'] for r in rows],'-o',label='Anytime / beta benchmark')
    ax.plot(ks,[r['n_CG23']/r['n_beta'] for r in rows],'--s',label='CG23 / beta benchmark')
    ax.plot(ks,[r['any_CG'] for r in rows],':^',label='Anytime / CG23')
    ax.set(xscale='log',xlabel='Boundary complexity $k$',ylabel='Ratio of risk certificates',ylim=(.95,4.8));ax.legend(fontsize=6.8);ax.tick_params(labelsize=8);ax.xaxis.label.set_size(8.5);ax.yaxis.label.set_size(8.5)
    save(fig,'final_complexity_price')
    print(json.dumps(rows,indent=2))
if __name__=='__main__':
    (ROOT/'tables').mkdir(exist_ok=True);main()
