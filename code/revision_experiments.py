"""Verified fixed-time comparison, complexity sweep, selection and coverage diagnostics.

The GC upper equation is taken from Garatti--Campi (2022), Eq.(4), and
cross-checked against Campi--Garatti (2023), Theorem 7. Fixed-time curves
are not used as stopping certificates. All seeds and protocols are explicit.
"""
from pathlib import Path
import json, csv, math, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import beta
from certificates import AnytimeCertificate, gc_upper, gcc_refined_stages, fixed_bound
from certificates_original import AnytimeCertificate as OriginalCertificate
from experiments import lmi_path, save_figure
ROOT=Path(__file__).resolve().parents[1]

def main():
    start=time.perf_counter(); c=AnytimeCertificate(); old=OriginalCertificate()
    grid=[(1000,2),(3425,2),(3425,23),(10000,2),(10000,5),(10000,23),(10000,50),
          (100000,2),(100000,5),(100000,23),(100000,50),(100000,100)]
    rows=[]
    for n,k in grid:
        g=gc_upper(n,k,.05); f=c.single_time(n,k); a=c(n,k); o=old(n,k)
        rows.append([n,k,n*g,n*f,n*o,n*a,a/g,a/f])
    with open(ROOT/'data/fixed_time_comparison.csv','w',newline='') as f:
        w=csv.writer(f); w.writerow(['n','k','n_GC22','n_localized_fixed','n_original_anytime','n_revised_anytime','ratio_to_GC22','ratio_to_own_fixed']);w.writerows(rows)
    print('comparison',rows,flush=True)
    sweep=[]; n=100000
    for k in [1,2,4,5,8,16,23,32,50,64,100,128,256]:
        g=gc_upper(n,k,.05); a=c(n,k); o=old(n,k)
        sweep.append([n,k,n*g,n*a,n*o,a/g,o/g])
    np.savetxt(ROOT/'data/complexity_sweep.csv',sweep,delimiter=',',header='n,k,n_GC22,n_revised,n_original,ratio_revised,ratio_original',comments='')
    sw=np.array(sweep)
    fig,ax=plt.subplots(figsize=(3.5,2.45))
    ax.semilogx(sw[:,1],sw[:,5],'-o',markersize=3,label='Revised certificate')
    ax.semilogx(sw[:,1],sw[:,6],'--',label='Original schedule')
    ax.axhline(1,linestyle=':',linewidth=1)
    ax.set_xlabel('Observed complexity $k$'); ax.set_ylabel('Anytime / GC22 bound')
    ax.legend(fontsize=7,frameon=False); ax.grid(True,alpha=.25)
    save_figure(fig,'complexity_price')
    # Independent specification-change stream, all four decisions made after T=4000.
    p=lmi_path(np.random.default_rng(26090772),4000,check=True); cc=AnytimeCertificate(.05,d=2)
    u=np.array([cc(i+1,int(k)) for i,k in enumerate(p['k'])]); spec=[]
    for eps in [.1,.05,.02,.01]:
        hit=np.flatnonzero(u<=eps)
        if not len(hit): raise RuntimeError('Increase acquisition cap for specification example.')
        ix=int(hit[0]);spec.append(dict(eps=eps,prefix=ix+1,k=int(p['k'][ix]),U=float(u[ix]),risk=float(p['risk'][ix])))
    np.savetxt(ROOT/'data/specification_path.csv',np.c_[np.arange(1,4001),p['k'],p['risk'],u],delimiter=',',header='n,k,risk,U',comments='')
    # Sensitivity: retain original aligned evaluation cost as a separate metric.
    rng=np.random.default_rng(26090602); raw=[]; utable=np.ones((800,3))
    for n in range(1,801):
        for k in range(min(n,2)+1): utable[n-1,k]=cc(n,k)
    bn=-np.array([1.1,.7])/np.linalg.norm([1.1,.7])
    for rep in range(150):
        p=lmi_path(rng,800); uu=utable[np.arange(800),p['k']]
        elig=np.flatnonzero(uu<=.05)
        cl=.8+p['gain']@bn
        costs=(1+.1*np.sum(p['gain']**2,axis=1))/(1-cl**2)
        raw.append([rep,float(costs[elig[0]]),float(np.min(costs[elig]))])
    a=np.array(raw); aligned=dict(mean_first=float(a[:,1].mean()),mean_selected=float(a[:,2].mean()),
        relative_mean_reduction=float(1-a[:,2].mean()/a[:,1].mean()))
    np.savetxt(ROOT/'data/selection_sensitivity.csv',raw,delimiter=',',header='replicate,aligned_first_cost,aligned_selected_cost',comments='')
    # Finite-horizon scalar-minimum coverage stress test, using no dimension assertion beyond K<=1.
    stress=[]; N=5000; reps=2000; rng=np.random.default_rng(26090773)
    deltas=[.05,.5]; thresholds={d:np.array([AnytimeCertificate(d,d=1)(n,1) for n in []]) for d in []}
    for delta in deltas:
        cs=AnytimeCertificate(delta,d=1)
        thresholds[delta]=np.array([cs(n,1) for n in range(1,N+1)])
    counts={d:0 for d in deltas}
    for b in range(20):
        v=np.minimum.accumulate(rng.uniform(size=(100,N)),axis=1)
        for d in deltas: counts[d]+=int(np.any(v>thresholds[d][None,:],axis=1).sum())
    for d in deltas:
        x=counts[d]
        stress.append(dict(delta=d,N=N,replicates=reps,crossings=x,frequency=x/reps,
            binomial_upper95=float(beta.ppf(.95,x+1,reps-x)) if x<reps else 1.))
    # Paired cost improvement under the main nominal plant.
    import pandas as pd
    df=pd.read_csv(ROOT/'data/lmi_replicates.csv')
    first=df[df.method=='Anytime'].sort_values('replicate').nominal_cost.to_numpy()
    post=df[df.method=='Posthoc'].sort_values('replicate').nominal_cost.to_numpy()
    out=dict(comparison=rows,complexity_sweep=sweep,specification_change=dict(seed=26090772,T=4000,selections=spec),
        original_cost_sensitivity=aligned,
        selection_gain=dict(relative_means=float(1-post.mean()/first.mean()),mean_paired_relative=float(np.mean(1-post/first)),
            improved_streams=int(np.sum(post<first-1e-10))),
        scalar_coverage=dict(seed=26090773,results=stress),elapsed_seconds=time.perf_counter()-start)
    (ROOT/'data/revision_results.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__': main()
