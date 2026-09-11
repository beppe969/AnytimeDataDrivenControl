"""Final-review experiments, fixed in advance: seed 26090781, 10,000 IID streams.

Same scalar-state/two-input scenario QP as experiments.py. The data law and design
objective are unchanged. Risk is used only for evaluation. All stopping rules use
(n,K_n) and their stated target; budgets are delta=.05, eps=.05 or .15.
The Numba implementation is cross-checked against the original polygon routine.
"""
from pathlib import Path
import csv,json,time,math
import numpy as np
from numba import njit
from scipy.special import gammaln,logsumexp
from scipy.stats import beta
from certificates import AnytimeCertificate,classical_sample_size,gcc_refined_stages
from experiments import lmi_path
ROOT=Path(__file__).resolve().parents[1]


def cg23_upper(n,k,delta):
    """CG23 Theorem 4, Eq.(1): (delta/n)sum_{m=k}^{n-1} C(m,k)/C(n,k)(1-u)^{m-n}=1."""
    if n<1 or not 0<=k<=n or not 0<delta<1: raise ValueError('Invalid arguments')
    if k==n:return 1.
    ms=np.arange(k,n,dtype=float)
    lr=np.zeros_like(ms)
    for i in range(k):lr+=np.log((ms-i)/(n-i))
    def f(u):return math.log(delta/n)+float(logsumexp(lr+(ms-n)*math.log1p(-u)))
    lo=0.;hi=min(1-1e-15,(k+2*math.sqrt(k+1)*(math.sqrt(math.log(k+1))+4+math.sqrt(-math.log(delta)))-math.log(delta))/n)
    while f(hi)<0:hi=1-(1-hi)/2
    for _ in range(48):
        mid=(lo+hi)/2
        if f(mid)<0:lo=mid
        else:hi=mid
    return float(np.nextafter(hi,1.))

@njit(cache=True)
def clip_fast(poly,m,b0,b1,rhs,out):
    count=0
    for j in range(m):
        h=(j+1)%m
        p0,p1=poly[j,0],poly[j,1];q0,q1=poly[h,0],poly[h,1]
        vp=p0*b0+p1*b1-rhs;vq=q0*b0+q1*b1-rhs
        if vp<=0:
            out[count,0]=p0;out[count,1]=p1;count+=1
        if (vp<=0)!=(vq<=0):
            t=vp/(vp-vq);out[count,0]=p0+t*(q0-p0);out[count,1]=p1+t*(q1-p1);count+=1
    return count

@njit(cache=True)
def fast_path(angles,a):
    n=len(a);normals=np.empty((n,2));gains=np.empty((n,2));ks=np.zeros(n,np.int64);risks=np.empty(n)
    for i in range(n):normals[i,0]=math.cos(angles[i]);normals[i,1]=math.sin(angles[i])
    p=np.empty((2*n+12,2));tmp=np.empty_like(p)
    p[0,0]=-2;p[0,1]=-2;p[1,0]=2;p[1,1]=-2;p[2,0]=2;p[2,1]=2;p[3,0]=-2;p[3,1]=2;m=4
    g0,g1=1.1,.7;k=0;updates=0;max_feas=0.
    for i in range(n):
        b0,b1=normals[i,0],normals[i,1]
        changes=abs(a[i]+b0*g0+b1*g1)>.9+2e-14
        m=clip_fast(p,m,b0,b1,.9-a[i],tmp);m=clip_fast(tmp,m,-b0,-b1,.9+a[i],p)
        if m<3:raise ValueError('Degenerate or infeasible polygon')
        if changes:
            best=1e100
            for j in range(m):
                h=(j+1)%m;v0=p[h,0]-p[j,0];v1=p[h,1]-p[j,1];den=v0*v0+v1*v1
                if den<=1e-28:continue
                t=max(0.,min(1.,((1.1-p[j,0])*v0+(.7-p[j,1])*v1)/den))
                x0=p[j,0]+t*v0;x1=p[j,1]+t*v1
                cost=(x0-1.1)**2+(x1-.7)**2
                if cost<best:best=cost;g0=x0;g1=x1
            k=0
            for j in range(i+1):
                res=abs(a[j]+normals[j,0]*g0+normals[j,1]*g1)-.9
                max_feas=max(max_feas,res)
                if abs(res)<2e-9:k+=1
            if k>2:
                k=0
                for j in range(i+1):
                    res=abs(a[j]+normals[j,0]*g0+normals[j,1]*g1)-.9
                    if abs(res)<2e-11:k+=1
            if k<1 or k>2:raise ValueError('Ambiguous active set')
            updates+=1
        gains[i,0]=g0;gains[i,1]=g1;ks[i]=k
        r=math.sqrt(g0*g0+g1*g1)
        if r<=.2:risks[i]=0.
        else:
            theta=math.acos(.2/r);risks[i]=(r*math.sin(theta)-.2*theta)/(.7*math.pi)
    return gains,ks,risks,updates,max_feas


def intervals(mask):
    starts=np.flatnonzero(mask&~np.r_[False,mask[:-1]])
    ends=np.flatnonzero(mask&~np.r_[mask[1:],False])
    return [[int(s+1),int(t+1)] for s,t in zip(starts,ends)]


def cp(x,n):
    return [0. if x==0 else float(beta.ppf(.025,x,n-x+1)),
            1. if x==n else float(beta.ppf(.975,x+1,n-x))]


def main(reps=10000,nmax=800):
    start=time.perf_counter();delta=.05;targets=[.05,.15];seed=26090781
    # Independent validation of fast solver on identical inputs.
    maxgain=0.;maxrisk=0.;ksame=True
    for s in range(30):
        rng=np.random.default_rng(26090800+s);aa=rng.uniform(0,2*np.pi,nmax);a=rng.uniform(-.7,.7,nmax)
        fast=fast_path(aa,a);slow=lmi_path(np.random.default_rng(26090800+s),nmax,check=True)
        maxgain=max(maxgain,float(np.max(abs(fast[0]-slow['gain']))));maxrisk=max(maxrisk,float(np.max(abs(fast[2]-slow['risk']))));ksame&=bool(np.array_equal(fast[1],slow['k']))
    assert maxgain<2e-10 and maxrisk<2e-10 and ksame
    print('fast solver verified',maxgain,maxrisk,ksame,flush=True)
    cert=AnytimeCertificate(delta,d=2);ut=np.ones((nmax,3));gt=np.ones_like(ut)
    for n in range(1,nmax+1):
        for k in range(min(n,2)+1):ut[n-1,k]=cert(n,k);gt[n-1,k]=cg23_upper(n,k,delta)
    plans={e:dict(N=classical_sample_size(2,e,delta),stages=list(map(int,gcc_refined_stages(2,e,delta)))) for e in targets}
    print('plans',plans,flush=True)
    rows=[];archive=[];pathcounts=dict(CG23=0,Anytime=0);updates=0;maxfeas=0.;nonmonotone=0;example=None
    rng=np.random.default_rng(seed)
    for rep in range(reps):
        aa=rng.uniform(0,2*np.pi,nmax);a=rng.uniform(-.7,.7,nmax)
        gain,k,v,changes,feas=fast_path(aa,a);updates+=changes;maxfeas=max(maxfeas,feas)
        u=ut[np.arange(nmax),k];g=gt[np.arange(nmax),k]
        pathcounts['CG23']+=int(np.any(v>g+1e-12));pathcounts['Anytime']+=int(np.any(v>u+1e-12));nonmonotone+=int(np.any(np.diff(v)>1e-10))
        for eps in targets:
            first={}
            for name,bound in [('CG23-repeated',g),('Anytime',u)]:
                h=np.flatnonzero(bound<=eps);first[name]=int(h[0]) if len(h) else None
            first['Fixed N']=plans[eps]['N']-1
            stages=plans[eps]['stages'];first['GCC refined']=stages[-1]-1
            for j,n in enumerate(stages):
                if k[n-1]<=j:first['GCC refined']=n-1;break
            for name,ix in first.items():
                if ix is None:rows.append([rep,eps,name,0,nmax,np.nan,np.nan,0]);continue
                rows.append([rep,eps,name,1,ix+1,int(k[ix]),float(v[ix]),int(v[ix]>eps)])
            tau=first['Anytime']
            if tau is not None:
                bad=u>eps;bad[:tau+1]=False
                runs=intervals(bad)
                # latest-only service is unavailable at these indices. Archive retains prior controller.
                ar=[rep,eps,tau+1,len(runs),int(bad.sum()),max([b-a+1 for a,b in runs],default=0)]
                archive.append(ar)
                if eps==.05 and len(runs) and example is None:
                    example=dict(replicate=rep,first=tau+1,intervals=runs,first_k=int(k[tau]),first_risk=float(v[tau]),first_U=float(u[tau]))
                    np.savetxt(ROOT/'data/recrossing_example.csv',np.c_[np.arange(1,nmax+1),k,v,u,g,gain],delimiter=',',header='n,k,risk,U,CG23,gain1,gain2',comments='')
        if (rep+1)%2000==0:print('replicates',rep+1,'elapsed',round(time.perf_counter()-start,2),flush=True)
    with open(ROOT/'data/stopping_audit_replicates.csv','w',newline='') as f:
        w=csv.writer(f);w.writerow(['replicate','eps','method','deployed','acquired','k','risk','risk_exceeds_eps']);w.writerows(rows)
    with open(ROOT/'data/archive_audit_replicates.csv','w',newline='') as f:
        w=csv.writer(f);w.writerow(['replicate','eps','first','recrossing_intervals','unavailable_prefixes','longest_interval']);w.writerows(archive)
    summary=[]
    for eps in targets:
        for name in ['Fixed N','CG23-repeated','GCC refined','Anytime']:
            rr=[r for r in rows if r[1]==eps and r[2]==name];x=sum(r[-1] for r in rr)
            ns=np.array([r[4] for r in rr]);vs=[r[6] for r in rr if r[3]]
            summary.append(dict(eps=eps,method=name,replicates=reps,deployed=sum(r[3] for r in rr),mean_n=float(ns.mean()),sd_n=float(ns.std(ddof=1)),mean_risk=float(np.mean(vs)),failures=x,frequency=x/reps,CP95=cp(x,reps)))
    ars=[]
    for eps in targets:
        ar=np.array([r for r in archive if r[1]==eps]);x=int(np.sum(ar[:,3]>0))
        ars.append(dict(eps=eps,streams_with_recrossings=x,frequency=x/reps,CP95=cp(x,reps),total_intervals=int(ar[:,3].sum()),total_unavailable_prefixes=int(ar[:,4].sum()),mean_unavailable=float(ar[:,4].mean()),mean_unavailable_given_recrossing=float(ar[ar[:,3]>0,4].mean()) if x else 0.,max_interval=int(ar[:,5].max()),mean_first=float(ar[:,2].mean()),latest_unavailability_uniform_post_first=float(np.mean(ar[:,4]/(nmax-ar[:,2])))))
    import pandas as pd
    hull=pd.read_csv(ROOT/'data/disturbance_path.csv');uh=hull.U_hull.to_numpy();h=int(np.flatnonzero(uh<=.02)[0]);bad=uh>.02;bad[:h+1]=False
    hullout=dict(first=h+1,intervals=intervals(bad),unavailable_prefixes=int(bad.sum()))
    out=dict(protocol=dict(seed=seed,nmax=nmax,replicates=reps,delta=delta,eps=targets,source='same QP as experiments.py',plans=plans),validation=dict(max_gain_difference=maxgain,max_risk_difference=maxrisk,complexities_equal=ksame,prefixes_compared=30*nmax),stopping=summary,archive=ars,first_recrossing_example=example,hull_recrossings=hullout,all_prefix_coverage=[dict(method=m,crossings=x,frequency=x/reps,CP95=cp(x,reps)) for m,x in pathcounts.items()],boundary_updates=updates,max_sample_residual=maxfeas,nonmonotone_streams=nonmonotone,elapsed_seconds=time.perf_counter()-start)
    (ROOT/'data/final_selection_results.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':main()
