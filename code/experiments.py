"""Reproducible studies for A Universal Iterated-Logarithm Law for Data-Driven Control Design.

Run from any directory: python code/experiments.py
No training/selection rule uses the evaluation distribution or its exact risks.
Seeds, parameters, machine-independent summaries, and raw paths are exported.
"""
from __future__ import annotations
from pathlib import Path
import json, math, time
import numpy as np
from scipy.spatial import ConvexHull
from scipy.linalg import solve_discrete_are, solve_discrete_lyapunov
from scipy.optimize import minimize, LinearConstraint, Bounds
from shapely.geometry import Polygon, box
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from certificates import AnytimeCertificate, classical_sample_size, gcc_stages, gcc_refined_stages

ROOT=Path(__file__).resolve().parents[1]
for name in ['figures','data','checks']:
    (ROOT/name).mkdir(exist_ok=True)


def save_figure(fig,name):
    for ax in fig.axes:
        ax.tick_params(labelsize=7)
        ax.xaxis.label.set_size(8)
        ax.yaxis.label.set_size(8)
    fig.tight_layout()
    fig.savefig(ROOT/'figures'/f'{name}.pdf',bbox_inches='tight')
    fig.savefig(ROOT/'figures'/f'{name}.png',dpi=170,bbox_inches='tight')
    plt.close(fig)


class EmpiricalHull:
    def __init__(self):
        self.points=np.zeros((1,2)); self.ids=np.array([-1]); self.hull=None
    def add(self,w,index):
        if self.hull is not None:
            e=self.hull.equations
            if np.max(e[:,:2]@w+e[:,2])<=0:
                return False
        self.points=np.vstack([self.points,w]); self.ids=np.r_[self.ids,index]
        if len(self.points)>=3:
            self.hull=ConvexHull(self.points)
            v=self.hull.vertices
            self.points=self.points[v]; self.ids=self.ids[v]
            self.hull=ConvexHull(self.points)
        return True
    @property
    def k(self): return int(np.sum(self.ids>=0))
    @property
    def area(self): return float(self.hull.volume) if self.hull is not None else 0.


def clip(poly,normal,rhs):
    if len(poly)==0: raise RuntimeError('Empty design polygon.')
    val=poly@normal-rhs
    if np.max(val)<=0: return poly
    out=[]
    for i,p in enumerate(poly):
        q=poly[(i+1)%len(poly)]
        a,b=val[i],val[(i+1)%len(poly)]
        if a<=0: out.append(p)
        if (a<=0) != (b<=0):
            out.append(p+(q-p)*(a/(a-b)))
    if not out: raise RuntimeError('Unexpected infeasibility.')
    return np.asarray(out)


def project(poly,desired):
    best=None; bestcost=math.inf
    for i,p in enumerate(poly):
        q=poly[(i+1)%len(poly)]; v=q-p
        if v@v<=1e-28: continue
        t=np.clip((desired-p)@v/(v@v),0,1)
        x=p+t*v; cost=np.sum((x-desired)**2)
        if cost<bestcost: best,bestcost=x,cost
    if best is None: raise RuntimeError('Projection failed.')
    return best


def exact_lmi_risk(gain,rho=.9,amax=.7):
    r=float(np.linalg.norm(gain)); h=rho-amax
    if r<=h: return 0.
    if r>=rho+amax: raise ValueError('Closed-form branch requires norm < rho+amax.')
    theta=math.acos(h/r)
    return (r*math.sin(theta)-h*theta)/(math.pi*amax)


def lmi_path(rng,nmax,check=False):
    angles=rng.uniform(0,2*np.pi,nmax)
    normals=np.c_[np.cos(angles),np.sin(angles)]
    a=rng.uniform(-.7,.7,nmax)
    desired=np.array([1.1,.7]); gain=desired.copy()
    poly=np.array([[-2.,-2.],[2.,-2.],[2.,2.],[-2.,2.]])
    gains=np.zeros((nmax,2)); ks=np.zeros(nmax,dtype=int); risks=np.zeros(nmax)
    support=[]; updates=[]; audit=[]
    for i in range(nmax):
        b=normals[i]
        changes=abs(a[i]+b@gain)>.9+2e-14
        poly=clip(poly,b,.9-a[i]); poly=clip(poly,-b,.9+a[i])
        if changes:
            gain=project(poly,desired)
            res=np.abs(a[:i+1]+normals[:i+1]@gain)-.9
            support=np.flatnonzero(np.abs(res)<2e-9).tolist()
            if len(support)>2: raise RuntimeError('Ambiguous active set.')
            updates.append(i+1)
            if check:
                pp=np.array([[-2.,-2.],[2.,-2.],[2.,2.],[-2.,2.]])
                for t in support:
                    pp=clip(pp,normals[t],.9-a[t]); pp=clip(pp,-normals[t],.9+a[t])
                reconstructed=project(pp,desired)
                error=float(np.max(np.abs(reconstructed-gain)))
                if error>1e-8: raise RuntimeError(f'Reconstruction error {error}')
                signs=np.sign(a[support]+normals[support]@gain)
                G=normals[support]*signs[:,None]
                lam=np.linalg.lstsq(G.T,desired-gain,rcond=None)[0]
                if np.min(lam)<-1e-8: raise RuntimeError('Negative multiplier.')
                audit.append((i+1,error,float(np.min(lam)),float(np.max(res))))
        gains[i]=gain; ks[i]=len(support); risks[i]=exact_lmi_risk(gain)
    return dict(gain=gains,k=ks,risk=risks,updates=updates,a=a,b=normals,audit=audit)


def tube_support(vertices,A,B,Q,R,J=160):
    P=solve_discrete_are(A,B,Q,R)
    L=-np.linalg.solve(R+B.T@P@B,B.T@P@A)
    F=A+B@L
    S=solve_discrete_lyapunov(F.T,np.eye(2)); Sinv=np.linalg.inv(S)
    gamma=math.sqrt(1-1/np.linalg.eigvalsh(S)[-1])
    rS=math.sqrt(np.max(np.einsum('ni,ij,nj->n',vertices,S,vertices)))
    directions=np.vstack([np.eye(2),-np.eye(2),L,-L])
    C=directions.copy(); sums=np.zeros(6)
    for _ in range(J):
        sums+=np.max(C@vertices.T,axis=1)
        C=C@F
    tails=rS*np.sqrt(np.einsum('ni,ij,nj->n',C,Sinv,C))/(1-gamma)
    return sums+tails,L,P,F,float(np.max(tails))


def nominal_mpc(vertices):
    A=np.array([[1.,.2],[0.,1.]]); B=np.array([[.1],[.25]])
    Q=np.diag([1.,.3]); R=np.array([[.15]])
    h,L,P,F,tail=tube_support(vertices,A,B,Q,R)
    xmax=np.array([1.2,.7]); umax=.65; N=20; x0=np.array([.8,-.2])
    xhi=xmax-h[:2]; xlo=-xmax+h[2:4]
    uhi=umax-h[4]; ulo=-umax+h[5]
    if np.any(x0>xhi) or np.any(x0<xlo) or ulo>=uhi:
        return dict(feasible=False,tube=h.tolist(),tail=tail)
    M=[]; c=[]
    for j in range(1,N+1):
        row=np.zeros((2,N))
        for t in range(j): row[:,t]=(np.linalg.matrix_power(A,j-1-t)@B).ravel()
        M.append(row); c.append(np.linalg.matrix_power(A,j)@x0)
    M=np.vstack(M); c=np.concatenate(c)
    Qbig=np.kron(np.eye(N),Q); Qbig[-2:,-2:]=P
    Hess=2*(M.T@Qbig@M+.15*np.eye(N)); lin=2*M.T@Qbig@c
    objective=lambda u: .5*u@Hess@u+lin@u+c@Qbig@c+x0@Q@x0
    grad=lambda u: Hess@u+lin
    cons=[LinearConstraint(M[:-2],np.tile(xlo,N-1)-c[:-2],np.tile(xhi,N-1)-c[:-2]),
          LinearConstraint(M[-2:],-c[-2:],-c[-2:])]
    ans=minimize(objective,np.zeros(N),jac=grad,method='SLSQP',bounds=Bounds(ulo,uhi),
                 constraints=cons,options={'ftol':1e-11,'maxiter':1000})
    if not ans.success: raise RuntimeError('MPC solve: '+ans.message)
    return dict(feasible=True,cost=float(ans.fun),v0=float(ans.x[0]),tube=h.tolist(),
                tail=tail,L=L.ravel().tolist(),pole_moduli=np.abs(np.linalg.eigvals(F)).tolist(),
                terminal_residual=float(np.max(np.abs(M[-2:]@ans.x+c[-2:]))))


def disturbance_study():
    nmax=8000; delta=.05; eps=.02; seed=26090601
    rng=np.random.default_rng(seed)
    G=np.array([[.03,.02],[0.,.015]])
    W=rng.uniform(-1,1,(nmax,2))@G.T
    corners=np.array([[-1,-1],[1,-1],[1,1],[-1,1]])@G.T
    domain=Polygon(corners); area=domain.area
    hull=EmpiricalHull(); lo=np.zeros(2); hi=np.zeros(2); ilo=np.full(2,-1); ihi=ilo.copy()
    cert=AnytimeCertificate(delta); boxcert=AnytimeCertificate(delta,d=4)
    rows=[]; selected={}; bounds=[]; snapshots={}; changes=0
    for i,w in enumerate(W):
        n=i+1; changes+=hull.add(w,i)
        for t in range(2):
            if w[t]<lo[t]: lo[t]=w[t]; ilo[t]=i
            if w[t]>hi[t]: hi[t]=w[t]; ihi[t]=i
        kb=len(set(np.r_[ilo,ihi][np.r_[ilo,ihi]>=0].tolist()))
        uh=cert(n,hull.k); ub=boxcert(n,kb)
        vh=max(0.,1-hull.area/area)
        vb=max(0.,1-domain.intersection(box(lo[0],lo[1],hi[0],hi[1])).area/area)
        un=cert.per_time(n,hull.k)
        rows.append([n,hull.k,vh,uh,un,kb,vb,ub])
        for key,u,k,v,verts in [('Hull',uh,hull.k,vh,hull.points),
                              ('Box',ub,kb,vb,np.array([[lo[0],lo[1]],[hi[0],lo[1]],[hi[0],hi[1]],[lo[0],hi[1]]])),
                              ('Hull_per_time',un,hull.k,vh,hull.points)]:
            if key not in selected and u<=eps:
                mpc=nominal_mpc(verts)
                if mpc['feasible']:
                    selected[key]=dict(n=n,k=k,risk=v,bound=u,horizon_bound=1-(1-u)**10,
                                       exact_escape=1-(1-v)**10,mpc=mpc)
                    snapshots[key]=verts.copy()
        if n in [100,500,1000,2000,4000,8000]:
            bounds.append(dict(n=n,k=hull.k,risk=vh,bound=uh,per_time=un))
    arr=np.asarray(rows)
    np.savetxt(ROOT/'data'/'disturbance_path.csv',arr,delimiter=',',
               header='n,k_hull,risk_hull,U_hull,U_per_time,k_box,risk_box,U_box',comments='')
    np.savez_compressed(ROOT/'data'/'disturbance_samples.npz',W=W,G=G,**snapshots)
    fig,ax=plt.subplots(figsize=(3.5,3.0))
    for col,label,ls in [(3,'Hull: anytime','-'),(4,'Hull: per-time','--'),
                         (2,'Hull: exact risk','-.'),(7,'Box: anytime',':')]:
        ax.loglog(arr[9:,0],arr[9:,col],linestyle=ls,label=label,linewidth=1.5)
    ax.axhline(eps,linestyle='--',linewidth=1,label=r'Target $\varepsilon=0.02$')
    ax.set_xlabel('Design sample size $n$'); ax.set_ylabel('Violation probability / upper bound')
    ax.legend(fontsize=7,loc='lower center',bbox_to_anchor=(.5,1.02),ncol=2,frameon=False,columnspacing=.8); ax.grid(True,which='major',alpha=.25)
    save_figure(fig,'disturbance_risk')
    fig,ax=plt.subplots(figsize=(3.5,2.0))
    cc=np.vstack([corners,corners[0]])
    ax.plot(cc[:,0],cc[:,1],'--',label='Evaluation support',linewidth=1.1)
    for key,label in [('Hull','Hull'),('Box','Box')]:
        pp=snapshots[key]; pp=np.vstack([pp,pp[0]])
        ax.plot(pp[:,0],pp[:,1],label=label,linewidth=1.5)
    ax.set_xlabel('$w_1$'); ax.set_ylabel('$w_2$'); ax.set_aspect('equal',adjustable='box')
    ax.legend(fontsize=7,loc='lower center',bbox_to_anchor=(.5,1.02),ncol=3,frameon=False,columnspacing=.6,handletextpad=.3); ax.grid(True,alpha=.25)
    save_figure(fig,'disturbance_geometry')
    summary=dict(seed=seed,nmax=nmax,delta=delta,eps=eps,G=G.tolist(),selected=selected,
                 checkpoints=bounds,boundary_changes=int(changes),
                 any_certificate_failure=bool(np.any(arr[:,2]>arr[:,3]+1e-10)),
                 maximum_boundary=int(arr[:,1].max()),
                 exact_risk='Polygon area, using the evaluation-only uniform parallelogram law.')
    return summary


def lmi_study():
    delta=.05; eps=.05; nmax=800; reps=150; seed=26090602
    cert=AnytimeCertificate(delta,d=2)
    us=np.ones((nmax,3)); un=us.copy()
    for n in range(1,nmax+1):
        for k in range(min(2,n)+1):
            us[n-1,k]=cert(n,k); un[n-1,k]=cert.per_time(n,k)
    stages=gcc_stages(2,eps,delta); refined=gcc_refined_stages(2,eps,delta); fixed=classical_sample_size(2,eps,delta)
    rng=np.random.default_rng(seed); rows=[]; representative=None; path_failures=0
    gain_audit=[]; all_audits=[]; nonmonotone_paths=0; global_rise=None
    for rep in range(reps):
        p=lmi_path(rng,nmax,check=True)
        all_audits.extend([(rep,*x) for x in p['audit']])
        ri=np.flatnonzero(np.diff(p['risk'])>1e-10)
        nonmonotone_paths+=int(len(ri)>0)
        if len(ri):
            t=int(ri[np.argmax(np.diff(p['risk'])[ri])])+2
            jump=float(p['risk'][t-1]-p['risk'][t-2])
            if global_rise is None or jump>global_rise['jump']:
                global_rise=dict(replicate=rep,n=t,before=float(p['risk'][t-2]),after=float(p['risk'][t-1]),jump=jump)
        ks=p['k']; risks=p['risk']; uu=us[np.arange(nmax),ks]; nn=un[np.arange(nmax),ks]
        path_failures+=int(np.any(risks>uu+1e-10))
        tau=int(np.flatnonzero(uu<=eps)[0])+1
        taun=int(np.flatnonzero(nn<=eps)[0])+1
        taug=next(N for j,N in enumerate(stages) if ks[N-1]<=j)
        taugr=next(N for j,N in enumerate(refined) if ks[N-1]<=j)
        desired=np.array([1.1,.7]); bnom=np.array([.7,-1.1])/np.linalg.norm(desired)
        def cost(i):
            g=p['gain'][i-1]; cl=.6+bnom@g
            if abs(cl)>=1: return float("inf")
            return float((1+.1*(g@g))/(1-cl*cl))
        eligible=np.flatnonzero(uu<=eps)+1
        chosen=int(min(eligible,key=cost))
        for key,n in [('Anytime',tau),('Per_time',taun),('GCC_T1',taug),('GCC_T4',taugr),('Fixed_N',fixed),('Posthoc',chosen)]:
            rows.append([rep,key,n,int(ks[n-1]),float(risks[n-1]),cost(n),float(uu[n-1])])
        if rep==28:
            representative=(p,uu,nn); gain_audit=p['audit']
    p,uu,nn=representative
    np.savetxt(ROOT/'data'/'lmi_path.csv',np.c_[np.arange(1,nmax+1),p['k'],p['risk'],uu,nn,p['gain']],
               delimiter=',',header='n,k,risk,U_anytime,U_per_time,L1,L2',comments='')
    import csv
    with open(ROOT/'data'/'lmi_replicates.csv','w',newline='') as f:
        w=csv.writer(f); w.writerow(['replicate','method','n','k','risk','nominal_cost','anytime_U']); w.writerows(rows)
    with open(ROOT/'checks'/'lmi_support_audit.csv','w',newline='') as f:
        w=csv.writer(f); w.writerow(['replicate','n','reconstruction_error','minimum_multiplier','maximum_constraint_residual']); w.writerows(all_audits)
    records={}
    for key in ['Anytime','Per_time','GCC_T1','GCC_T4','Fixed_N','Posthoc']:
        a=np.asarray([[r[i] for i in [2,3,4,5]] for r in rows if r[1]==key],float)
        records[key]=dict(mean_n=float(a[:,0].mean()),sd_n=float(a[:,0].std(ddof=1)),
                          min_n=int(a[:,0].min()),max_n=int(a[:,0].max()),mean_risk=float(a[:,2].mean()),
                          max_risk=float(a[:,2].max()),mean_nominal_cost=float(a[:,3].mean()),
                          risk_failures=int(np.sum(a[:,2]>eps)))
    fig,ax=plt.subplots(figsize=(3.5,3.0))
    n=np.arange(1,nmax+1)
    ax.step(n,p['risk'],where='post',label='Exact controller risk',linewidth=1.4)
    ax.set_xscale('log'); ax.set_yscale('log')
    upd=np.array(p['updates'],dtype=int)
    ax.plot(upd,p['risk'][upd-1],'o',markersize=3,label='Boundary changes')
    ax.loglog(n,uu,'--',label='Anytime certificate',linewidth=1.4)
    ax.loglog(n,nn,':',label='Per-time allocation',linewidth=1.4)
    ax.axhline(eps,linestyle='-.',linewidth=1,label=r'Target $\varepsilon=0.05$')
    ax.set_xlabel('Design sample size $n$'); ax.set_ylabel('Violation probability / upper bound')
    ax.legend(fontsize=7,loc='lower center',bbox_to_anchor=(.5,1.02),ncol=2,frameon=False,columnspacing=.8); ax.grid(True,which='major',alpha=.25)
    save_figure(fig,'lmi_risk')
    fig,ax=plt.subplots(figsize=(3.5,3.0))
    ax.plot(n,p['gain'][:,0],label='$L_{n,1}$',linewidth=1.3)
    ax.plot(n,p['gain'][:,1],label='$L_{n,2}$',linewidth=1.3)
    ax.set_xscale('log'); ax.set_xlabel('Design sample size $n$'); ax.set_ylabel('Feedback gain')
    ax.legend(fontsize=7); ax.grid(True,alpha=.25)
    save_figure(fig,'lmi_gain')
    increments=np.diff(p['risk']); rises=np.flatnonzero(increments>1e-10)
    rise=int(rises[np.argmax(increments[rises])])+2 if len(rises) else None
    summary=dict(seed=seed,nmax=nmax,replicates=reps,delta=delta,eps=eps,GCC_T1_stages=stages,GCC_T4_stages=[int(x) for x in refined],representative_replicate=28,
                 classical_N=fixed,methods=records,path_failures=path_failures,
                 representative_updates=len(p['updates']),representative_risk_increases=len(rises),
                 largest_rise_n=rise,largest_rise_before=float(p['risk'][rise-2]) if rise else None,
                 largest_rise_after=float(p['risk'][rise-1]) if rise else None,
                 nonmonotone_paths=nonmonotone_paths,largest_risk_increase=global_rise,
                 audited_boundaries=len(all_audits),
                 max_reconstruction_error=max(x[2] for x in all_audits),
                 min_active_multiplier=min(x[3] for x in all_audits))
    return summary


def main():
    t=time.perf_counter()
    print('Running disturbance study...',flush=True)
    d=disturbance_study(); print(json.dumps(d,indent=2),flush=True)
    print('Running LMI study...',flush=True)
    l=lmi_study(); print(json.dumps(l,indent=2),flush=True)
    out=dict(disturbance=d,lmi=l,elapsed_seconds=time.perf_counter()-t,
             software=dict(numpy=np.__version__,scipy=__import__('scipy').__version__,
                           matplotlib=matplotlib.__version__,shapely=__import__('shapely').__version__))
    (ROOT/'data'/'results.json').write_text(json.dumps(out,indent=2))
    print('Completed in',out['elapsed_seconds'],'seconds',flush=True)

if __name__=='__main__': main()
