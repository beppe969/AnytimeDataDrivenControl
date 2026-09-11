"""Figures and tables for the third-review revision, from supplied output files."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from experiments import nominal_mpc
ROOT=Path(__file__).resolve().parents[1]


def save(fig,name):
    for ax in fig.axes:
        ax.tick_params(labelsize=8)
        ax.xaxis.label.set_size(9);ax.yaxis.label.set_size(9)
    fig.tight_layout(pad=.6)
    fig.savefig(ROOT/'figures'/f'{name}.pdf',bbox_inches='tight')
    fig.savefig(ROOT/'figures'/f'{name}.png',dpi=210,bbox_inches='tight')
    plt.close(fig)


def main():
    x=pd.read_csv(ROOT/'data/horizon_crossover_sweep.csv')
    fig,ax=plt.subplots(figsize=(3.48,2.45))
    for k,ls in zip([2,5,23,100],['-','--','-.',':']):
        r=x[x.k==k]
        ax.plot(r.N,r.any_union,ls,label=f'$k={k}$',lw=1.4)
    ax.axhline(1,ls='--',lw=.7)
    ax.set(xscale='log',xlim=(1e3,1e12),ylim=(.35,1.30),xlabel='Endpoint $n=N$',ylabel='Anytime / finite-horizon bound')
    ax.set_xticks([1e3,1e5,1e7,1e9,1e12])
    ax.legend(fontsize=7.5,ncol=2,loc='lower left',frameon=False)
    save(fig,'horizon_comparison')

    report=json.loads((ROOT/'data/horizon_baseline_results.json').read_text())
    rows=report['monitoring']['coverage']
    rr=[r for r in rows if r['method']=='CG23']
    xx=np.array([r['inspected_through'] for r in rr]);yy=100*np.array([r['frequency'] for r in rr])
    er=np.array([[100*(r['frequency']-r['CP95'][0]) for r in rr],[100*(r['CP95'][1]-r['frequency']) for r in rr]])
    fig,ax=plt.subplots(figsize=(4.85,2.6))
    ax.errorbar(xx,yy,yerr=er,marker='o',lw=1.3,capsize=3,label='Unadjusted CG23')
    ax.plot(xx,np.zeros(len(xx)),'s--',label='FH, summable CG23, and anytime')
    ax.axhline(5,ls=':',lw=1,label='Nominal 5%')
    ax.set(xscale='log',xlim=(90,9000),ylim=(-.35,8),xlabel='Last inspected prefix',ylabel='Streams with a crossing (%)')
    ax.set_xticks([100,400,1600,8000],labels=['100','400','1600','8000'])
    ax.legend(fontsize=8,frameon=False,loc='upper left')
    save(fig,'monitoring_by_horizon')

    d=pd.read_csv(ROOT/'data/horizon_disturbance_path.csv').iloc[9:]
    fig,ax=plt.subplots(figsize=(3.48,2.55))
    for c,label,ls in [('U_hull','Hull: anytime','-'),('finite_hull','Hull: FH ($N=8000$)','--'),('risk_hull','Hull: exact risk','-.'),('U_box','Box: anytime',':')]:
        ax.loglog(d.n,d[c],ls,label=label,lw=1.4)
    ax.axhline(.02,ls='--',lw=.8,label='Target 0.02')
    ax.set(xlabel='Design sample size $n$',ylabel='Risk / upper bound')
    ax.legend(fontsize=7,frameon=False,loc='lower left')
    save(fig,'disturbance_risk')

    # Main table and supplement endpoint table, each generated from exact data.
    ep=pd.read_csv(ROOT/'data/horizon_endpoint_comparison.csv')
    with open(ROOT/'tables/horizon_main_rows.tex','w') as f:
        for r in ep[ep.N==100000].itertuples():
            f.write(f'{r.k} & {r.beta:.2f} & {r.CG23:.2f} & {r.finite_horizon:.2f} & {r.anytime:.2f} & {r.any_union:.3f}'+r'\\'+'\n')
    with open(ROOT/'tables/horizon_endpoints.tex','w') as f:
        f.write(r'\begin{tabular}{rrrrrrrr}'+'\n'+r'\toprule'+'\n')
        f.write(r'$N$ & $k$ & Beta & CG23 & FH & Summable & Any & Any/FH\\'+'\n'+r'\midrule'+'\n')
        for r in ep.itertuples():
            f.write(f'{r.N} & {r.k} & {r.beta:.2f} & {r.CG23:.2f} & {r.finite_horizon:.2f} & {r.summable_CG23:.2f} & {r.anytime:.2f} & {r.any_union:.3f}'+r'\\'+'\n')
        f.write(r'\bottomrule'+'\n'+r'\end{tabular}'+'\n')
    # Direct deterministic feasibility/cost confirmation at the new selected prefixes.
    W=np.load(ROOT/'data/disturbance_samples.npz')['W'];out=[]
    for rr in report['disturbance']:
        n=rr['first_qualified'];points=np.vstack([np.zeros((1,2)),W[:n]])
        if rr['scheme']=='hull':
            points=points[ConvexHull(points).vertices]
        else:
            lo=points.min(axis=0);hi=points.max(axis=0)
            points=np.array([[lo[0],lo[1]],[hi[0],lo[1]],[hi[0],hi[1]],[lo[0],hi[1]]])
        audit=nominal_mpc(points);assert audit['feasible']
        out.append(dict(**rr,mpc=audit))
    (ROOT/'checks/horizon_mpc_verification.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
