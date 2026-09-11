"""Independent long-horizon audit of repeated certificate inspection.

Protocol fixed before this run: seed 26090782; 10,000 streams; 8,000 observations;
delta=.05; same scenario QP. Comparison is simultaneous path coverage, separate
from first certification at one fixed target. No stopping uses exact risks.
"""
from pathlib import Path
import json,time,math,csv
import numpy as np
from scipy.special import betainc
from scipy.optimize import brentq
from certificates import AnytimeCertificate,log_choose
from selection_audit import fast_path,cp,cg23_upper as cg23_direct
ROOT=Path(__file__).resolve().parents[1]


def cg23_beta(n,k,delta):
    """Exact finite sum via the negative-binomial/binomial identity."""
    if k==n:return 1.
    lc=log_choose(n,k)
    def f(u):
        tail=float(betainc(k+1,n-k,u))
        if tail<=0:return -math.inf
        return math.log(delta/n)+math.log(tail)-(k+1)*math.log(u)-lc-(n-k)*math.log1p(-u)
    lo=max(1e-15,k/n)
    hi=min(1-1e-15,(k+2*math.sqrt(k+1)*(math.sqrt(math.log(k+1))+4+math.sqrt(-math.log(delta)))-math.log(delta))/n)
    while f(hi)<0:hi=1-(1-hi)/2
    root=brentq(f,lo,hi,xtol=2e-15,rtol=1e-14)
    return float(np.nextafter(root,1.))


def main():
    start=time.perf_counter();delta=.05;reps=10000;N=8000;seed=26090782
    maxerr=0.;checks=0
    for n in [5,10,50,250,1000,8000]:
        for k in [0,1,2,min(7,n-1)]:
            for d in [.05,.5]:
                x=cg23_direct(n,k,d);y=cg23_beta(n,k,d);maxerr=max(maxerr,abs(x-y));checks+=1
    assert maxerr<3e-12
    cert=AnytimeCertificate(delta,d=2);ut=np.ones((N,3));gt=np.ones_like(ut)
    for n in range(1,N+1):
        for k in range(min(n,2)+1):ut[n-1,k]=cert(n,k);gt[n-1,k]=cg23_beta(n,k,delta)
    rng=np.random.default_rng(seed);rows=[];counts={(h,m):0 for h in [800,8000] for m in ['CG23','Anytime']};example=None;maxfeas=0;updates=0
    for rep in range(reps):
        aa=rng.uniform(0,2*np.pi,N);a=rng.uniform(-.7,.7,N)
        gain,k,v,up,feas=fast_path(aa,a);maxfeas=max(maxfeas,feas);updates+=up
        u=ut[np.arange(N),k];g=gt[np.arange(N),k]
        cg=np.flatnonzero(v>g+1e-12);anyt=np.flatnonzero(v>u+1e-12)
        firstcg=int(cg[0]+1) if len(cg) else 0;firstany=int(anyt[0]+1) if len(anyt) else 0
        rows.append([rep,firstcg,firstany,int(k[-1]),float(v[-1])])
        for h in [800,8000]:
            counts[h,'CG23']+=int(firstcg>0 and firstcg<=h);counts[h,'Anytime']+=int(firstany>0 and firstany<=h)
        if example is None and firstcg>800:
            example=dict(replicate=rep,first_CG23_crossing=firstcg,k=int(k[firstcg-1]),risk=float(v[firstcg-1]),CG23=float(g[firstcg-1]),anytime=float(u[firstcg-1]))
            np.savetxt(ROOT/'data/monitoring_example.csv',np.c_[np.arange(1,N+1),k,v,g,u],delimiter=',',header='n,k,risk,CG23,U',comments='')
        if (rep+1)%2000==0:print('replicates',rep+1,'elapsed',time.perf_counter()-start,flush=True)
    with open(ROOT/'data/monitoring_replicates.csv','w',newline='') as f:
        w=csv.writer(f);w.writerow(['replicate','first_CG23_failure','first_anytime_failure','terminal_k','terminal_risk']);w.writerows(rows)
    summary=[dict(horizon=h,method=m,crossings=x,frequency=x/reps,CP95=cp(x,reps)) for (h,m),x in counts.items()]
    out=dict(protocol=dict(seed=seed,replicates=reps,horizon=N,delta=delta,source='unchanged QP family'),coverage=summary,validation=dict(direct_root_checks=checks,max_difference=maxerr,max_sample_residual=maxfeas,boundary_updates=updates),first_extended_crossing_example=example,elapsed_seconds=time.perf_counter()-start)
    (ROOT/'data/final_monitoring_results.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':main()
