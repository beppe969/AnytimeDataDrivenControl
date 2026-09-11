"""Finite-horizon and horizon-free CG23 baselines for the third review.

All comparators are evaluated on the unchanged, archived simulation protocols.
The added comparisons were specified after those data existed; no design rule,
seed, horizon, or performance objective is tuned here. A cap N only covers
prefixes <=N. The summable allocation has no terminal cap.

Run: python code/horizon_baselines.py
"""
from pathlib import Path
import csv, json, math, time
from functools import lru_cache
import numpy as np
import pandas as pd
from scipy.special import betainc
from scipy.optimize import brentq
from scipy.stats import beta
from certificates import AnytimeCertificate
from selection_audit import fast_path, cp, cg23_upper

ROOT=Path(__file__).resolve().parents[1]

@lru_cache(maxsize=180000)
def cg23(n: int,k: int,delta: float) -> float:
    """CG23 Theorem 4 root, evaluated in the risk-budget coordinate t=n*u.

    Negative-binomial identity avoids an O(n) sum. Cancellation of log(n)
    before evaluation preserves accuracy at large n and fixed k.
    """
    if n<1 or not 0<=k<=n or not 0<delta<1:
        raise ValueError('Require n>=1, 0<=k<=n, 0<delta<1.')
    if n==k:return 1.
    n=int(n);k=int(k)
    prodlog=math.fsum(math.log1p(-i/n) for i in range(k))
    c=math.log(delta)+math.lgamma(k+1)-prodlog
    def residual(t):
        u=t/n
        if u>=1:return math.inf
        tail=float(betainc(k+1,n-k,u))
        if tail<=0:return -math.inf
        return c+math.log(tail)-(k+1)*math.log(t)-(n-k)*math.log1p(-u)
    ell=-math.log(delta)
    lo=max(float(k),1e-12)
    hi=min(float(np.nextafter(float(n),0.)), k+2*math.sqrt(k+1)*(math.sqrt(math.log(k+1))+4+math.sqrt(ell))+ell)
    if residual(hi)<0: return 1.
    t=brentq(residual,lo,hi,xtol=2e-12,rtol=2e-14)
    return min(1.,float(np.nextafter(t/n,1.)))


def time_weight(n: int) -> float:
    """Positive weights summing to one over n>=1; stable telescoping form."""
    if n<1: raise ValueError('n must be positive.')
    x=n+math.e-1
    return math.log1p(1/x)/(math.log(x)*math.log(x+1))


def precompute(N):
    out={name:np.ones((N,3)) for name in ['CG23','Anytime','Finite horizon','Summable CG23']}
    cert=AnytimeCertificate(.05,d=2)
    for n in range(1,N+1):
        for k in range(min(2,n)+1):
            out['CG23'][n-1,k]=cg23(n,k,.05)
            out['Anytime'][n-1,k]=cert(n,k)
            out['Finite horizon'][n-1,k]=cg23(n,k,.05/N)
            out['Summable CG23'][n-1,k]=cg23(n,k,.05*time_weight(n))
    return out


def deterministic():
    start=time.perf_counter();checks=[]
    for n in [3,10,50,250,1000,8000]:
        for k in sorted(set([0,1,2,min(7,n-1)])):
            if k>=n:continue
            for eta in [.5,.05,.05/8000]:
                a=cg23(n,k,eta);b=cg23_upper(n,k,eta)
                checks.append(dict(n=n,k=k,delta=eta,fast=a,direct=b,error=abs(a-b)))
    assert max(r['error'] for r in checks)<5e-12
    pd.DataFrame(checks).to_csv(ROOT/'checks/horizon_root_checks.csv',index=False)
    cert=AnytimeCertificate(.05);rows=[]
    for N in [800,8000,100000]:
        for k in [2,5,23,100,1000]:
            if k>=N:continue
            cg=cg23(N,k,.05);fh=cg23(N,k,.05/N);u=cert(N,k)
            ss=cg23(N,k,.05*time_weight(N));q=beta.ppf(.95,k,N-k+1)
            rows.append(dict(N=N,k=k,beta=N*q,CG23=N*cg,finite_horizon=N*fh,anytime=N*u,summable_CG23=N*ss,per_time=N*cert.per_time(N,k),any_union=u/fh))
    pd.DataFrame(rows).to_csv(ROOT/'data/horizon_endpoint_comparison.csv',index=False)
    # This log-spaced sweep diagnoses endpoint comparisons, not random stopping.
    sweep=[]
    for k in [2,5,23,100]:
        for N in sorted(set(int(round(10.**x)) for x in np.linspace(3,12,181))):
            u=cert(N,k);fh=cg23(N,k,.05/N)
            sweep.append(dict(N=N,k=k,anytime=N*u,finite_horizon=N*fh,any_union=u/fh))
    pd.DataFrame(sweep).to_csv(ROOT/'data/horizon_crossover_sweep.csv',index=False)
    crossed={}
    for k in [2,5,23,100]:
        ss=[r for r in sweep if r['k']==k]
        first=next((i for i,r in enumerate(ss) if r['any_union']<=1),None)
        sustained=next((i for i,r in enumerate(ss) if all(x['any_union']<=1 for x in ss[i:])),None)
        crossed[k]=dict(first_grid_N=None if first is None else ss[first]['N'],first_grid_bracket=None if first in [None,0] else [ss[first-1]['N'],ss[first]['N']],below_at_all_later_grid_points_N=None if sustained is None else ss[sustained]['N'])
    return dict(root_checks=len(checks),max_root_error=max(r['error'] for r in checks),endpoint_rows=rows,crossover_grid=crossed,elapsed_seconds=time.perf_counter()-start)


def simulate(N,seed,short=False):
    t=time.perf_counter();reps=10000;grid=precompute(N);names=list(grid)
    rng=np.random.default_rng(seed);idx=np.arange(N);rows=[];stops=[];maxfeas=0.;updates=0;min_gap={name:math.inf for name in names};ambiguous={name:0 for name in names}
    for rep in range(reps):
        angles=rng.uniform(0,2*np.pi,N);a=rng.uniform(-.7,.7,N)
        gain,k,v,up,res=fast_path(angles,a);updates+=up;maxfeas=max(maxfeas,res)
        row=dict(replicate=rep,terminal_k=int(k[-1]),terminal_risk=float(v[-1]))
        for name in names:
            bounds=grid[name][idx,k];gap=np.abs(v-bounds);min_gap[name]=min(min_gap[name],float(gap.min()));ambiguous[name]+=int(np.sum(gap<=1e-12));cross=np.flatnonzero(v>bounds+1e-12)
            row[name]=0 if len(cross)==0 else int(cross[0]+1)
            if short:
                for eps in [.05,.15]:
                    ix=np.flatnonzero(bounds<=eps)
                    n=None if len(ix)==0 else int(ix[0]); deployed=n is not None
                    stops.append(dict(replicate=rep,eps=eps,method=name,deployed=int(deployed),acquired=N if n is None else n+1,k=-1 if n is None else int(k[n]),risk=np.nan if n is None else float(v[n]),risk_exceeds_eps=int(deployed and v[n]>eps)))
        rows.append(row)
        if (rep+1)%2000==0:print('cohort',N,'rep',rep+1,'elapsed',time.perf_counter()-t,flush=True)
    df=pd.DataFrame(rows);df.to_csv(ROOT/f'data/horizon_monitoring_{N}.csv',index=False)
    # Exact replay agreement confirms the original streams were retained.
    if short:
        old=pd.read_csv(ROOT/'data/stopping_audit_replicates.csv')
        new=pd.DataFrame(stops)
        for name,oldname in [('CG23','CG23-repeated'),('Anytime','Anytime')]:
            for eps in [.05,.15]:
                a=old[(old.method==oldname)&(old.eps==eps)].sort_values('replicate');b=new[(new.method==name)&(new.eps==eps)].sort_values('replicate')
                assert np.array_equal(a.acquired,b.acquired)
                assert np.max(abs(a.risk.to_numpy()-b.risk.to_numpy()))<1e-12
        new.to_csv(ROOT/'data/horizon_stopping_replicates.csv',index=False)
    else:
        old=pd.read_csv(ROOT/'data/monitoring_replicates.csv')
        assert np.array_equal(old.first_CG23_failure,df['CG23'])
        assert np.array_equal(old.first_anytime_failure,df['Anytime'])
        assert np.max(abs(old.terminal_risk-df.terminal_risk))<1e-12
    coverage=[]
    horizons=sorted(set([100,200,400,800,N]+([1600,3200,6400] if N>800 else [])))
    for h in horizons:
        for name in names:
            x=int(((df[name]>0)&(df[name]<=h)).sum())
            coverage.append(dict(cap=N,inspected_through=h,method=name,crossings=x,frequency=x/reps,CP95=cp(x,reps)))
    summary=[]
    if short:
        sdf=pd.DataFrame(stops)
        for (eps,name),d in sdf.groupby(['eps','method']):
            x=int(d.risk_exceeds_eps.sum())
            summary.append(dict(eps=float(eps),method=name,deployed=int(d.deployed.sum()),mean_n=float(d.acquired.mean()),sd_n=float(d.acquired.std()),mean_risk=float(d.risk.mean()),failures=x,frequency=x/reps,CP95=cp(x,reps)))
    return dict(protocol=dict(seed=seed,cap=N,replicates=reps,delta=.05,comparison_added_in_third_review=True),coverage=coverage,stopping=summary,replay_matches_archived_data=True,boundary_updates=updates,max_sample_residual=maxfeas,minimum_absolute_risk_bound_gap=min_gap,prefix_comparisons_within_tolerance=ambiguous,crossing_tolerance=1e-12,elapsed_seconds=time.perf_counter()-t)


def hull():
    df=pd.read_csv(ROOT/'data/disturbance_path.csv');N=len(df);extra={};summary=[]
    for shape in ['hull','box']:
        ks=df['k_'+shape].to_numpy(int)
        vals=np.array([cg23(n,int(k),.05/N) for n,k in zip(range(1,N+1),ks)])
        extra['finite_'+shape]=vals
        ns=np.flatnonzero(vals<=.02);n=int(ns[0]+1) if len(ns) else None
        summary.append(dict(scheme=shape,cap=N,first_qualified=n,k=None if n is None else int(ks[n-1]),risk=None if n is None else float(df['risk_'+shape].iloc[n-1]),certificate=None if n is None else float(vals[n-1])))
    out=pd.concat([df,pd.DataFrame(extra)],axis=1)
    out.to_csv(ROOT/'data/horizon_disturbance_path.csv',index=False)
    return summary


def main():
    start=time.perf_counter();out={'deterministic':deterministic()}
    print('deterministic completed',flush=True)
    out['short']=simulate(800,26090781,short=True)
    out['monitoring']=simulate(8000,26090782)
    out['disturbance']=hull();out['elapsed_seconds']=time.perf_counter()-start
    (ROOT/'data/horizon_baseline_results.json').write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items() if k!='deterministic'},indent=2),flush=True)
if __name__=='__main__':main()
