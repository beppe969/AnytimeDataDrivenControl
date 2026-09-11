"""Six-state, three-actuator mechanical scenario SDP with a fixed Lyapunov matrix.

The spectral-norm form is equivalent to the 12-by-12 affine LMI in the paper.
SLSQP solves the convex program; KKT and support-only reconstruction are audited.
The 64 uncertainty vertices certify a common strictly feasible feedback analytically
by convexity. Out-of-sample evaluation uses an independent Monte Carlo sample.
The reported anytime bound has a distribution-free interpretation conditional on
almost-sure support reconstruction for the sampling rule. The finite numerical
audits assess solved instances. Evaluation failures count contraction/LMI
violations, which may include plants whose closed-loop matrices remain stable.
"""
from pathlib import Path
import json, time, itertools, sys
import numpy as np
from scipy.linalg import solve_discrete_are
from scipy.optimize import minimize, nnls
from scipy.stats import beta
from certificates import AnytimeCertificate
ROOT=Path(__file__).resolve().parents[1]

def run():
    start=time.perf_counter(); d=6; dt=.05
    vec=np.array([[1,0,0],[1,-1,0],[0,1,-1],[0,0,1.]])
    base=vec.T@vec
    B=np.r_[np.zeros((3,3)),dt*np.eye(3)]
    def matrices(q):
        q=np.atleast_2d(q)
        Ks=base[None]+.20*np.einsum('ni,ij,ik->njk',q[:,:3],vec[:3],vec[:3])
        Ds=np.zeros((len(q),3,3)); Ds[:,np.arange(3),np.arange(3)]=.12+.03*q[:,3:]
        As=np.tile(np.eye(6),(len(q),1,1)); As[:,:3,3:]=dt*np.eye(3)
        As[:,3:,:3]=-dt*Ks; As[:,3:,3:]-=dt*Ds
        return As
    A0=matrices(np.zeros((1,6)))[0]
    P=solve_discrete_are(A0,B,np.diag([1.,1.,1.,.1,.1,.1]),.1*np.eye(3))
    L0=-np.linalg.solve(.1*np.eye(3)+B.T@P@B,B.T@P@A0)
    ev,Q=np.linalg.eigh(np.linalg.inv(P)); X=Q@np.diag(ev/ev.sum())@Q.T
    R=Q@np.diag(np.sqrt(ev/ev.sum()))@Q.T; Ri=np.linalg.inv(R); Br=Ri@B
    verts=np.array(list(itertools.product([-1.,1.],repeat=6)))
    Av=matrices(verts)
    bound=float(np.linalg.svd(Ri@(Av+B@L0)@R,compute_uv=False)[:,0].max())
    if bound>=1: raise RuntimeError(f'No strict common contraction: {bound}')
    rho=(1+bound)/2
    seed=26090761; rng=np.random.default_rng(seed); nmax=800
    q=rng.uniform(-1,1,(nmax,6)); As=matrices(q); Ts=Ri@As@R
    desired=(L0+4.*np.outer(np.array([1.,.3,-.2]),np.array([1.,.2,-.1,.3,.1,-.2]))).flatten()
    def valsgrad(l,inds):
        L=l.reshape(3,6); T=Ts[inds]+Br@L@R
        us,ss,vhs=np.linalg.svd(T,full_matrices=False)
        u=us[:,:,0]; v=vhs[:,0,:]
        grad=np.einsum('ni,nj->nij',u@Br,v@R.T).reshape(-1,18)
        return rho-ss[:,0],-grad
    solves=0; audits=[]; path=[]; l=L0.flatten().copy(); ks=0; supports=[]
    cert=AnytimeCertificate(.05,d=18); selected=None
    for n in range(1,nmax+1):
        # The first design is also optimized, starting from the common feasible gain.
        change=n==1 or valsgrad(l,np.array([n-1]))[0][0]<-1e-10
        if change:
            ids=np.arange(n)
            working=np.unique(np.r_[supports,n-1]).astype(int)
            for cutting in range(n+1):
                res=minimize(lambda x:.5*np.sum((x-desired)**2),l,jac=lambda x:x-desired,method='SLSQP',
                    constraints={'type':'ineq','fun':lambda x:valsgrad(x,working)[0],
                                 'jac':lambda x:valsgrad(x,working)[1]},
                    options={'ftol':2e-12,'maxiter':300})
                if not res.success: raise RuntimeError((n,res.message))
                allvals=valsgrad(res.x,ids)[0]
                if allvals.min()>=-1e-9: break
                working=np.unique(np.r_[working,np.argmin(allvals)]).astype(int)
                l=res.x
            l=res.x; vals,jac=valsgrad(l,ids)
            active=np.flatnonzero(vals<2e-7)
            lam=np.linalg.lstsq((-jac[active]).T,desired-l,rcond=None)[0] if len(active) else np.empty(0)
            if len(active) and lam.min()<-1e-6:
                raise RuntimeError((n,"negative dual",lam.tolist()))
            supports=active[lam>1e-7]; dual=lam[lam>1e-7]
            if len(supports)==0:
                raise RuntimeError((n,'Desired design feasible: add zero-support handling'))
            # Reconstruct from the recorded support constraints only.
            rec=minimize(lambda x:.5*np.sum((x-desired)**2),L0.flatten(),jac=lambda x:x-desired,method='SLSQP',
                constraints={'type':'ineq','fun':lambda x:valsgrad(x,supports)[0],
                             'jac':lambda x:valsgrad(x,supports)[1]},
                options={'ftol':2e-12,'maxiter':300})
            err=float(np.max(np.abs(rec.x-l)))
            G=-jac[supports]; s=np.linalg.svd(G,compute_uv=False)
            kkt=float(np.linalg.norm(l-desired+G.T@dual,np.inf))
            if vals.min()<-1e-7 or err>5e-5 or kkt>1e-4 or len(supports)>18:
                raise RuntimeError((n,'audit',vals.min(),err,kkt,len(supports)))
            ks=len(supports); solves+=1
            print("update",n,ks,"err",err,"KKT",kkt,flush=True)
            audits.append([n,ks,float(vals.min()),err,kkt,float(s[-1]),float(dual.min())])
        U=cert(n,ks) if ks<=n else 1.
        path.append([n,ks,U,.5*np.sum((l-desired)**2),*l.tolist()])
        if selected is None and U<=.05:
            selected=(n,ks,l.copy(),U)
        if n%200==0: print('mechanical',n,'k',ks,'U',U,'solves',solves,flush=True)
    if selected is None: selected=(nmax,ks,l.copy(),U)
    n,ks,l,U=selected
    test=np.random.default_rng(seed+1).uniform(-1,1,(50000,6))
    margin=np.linalg.svd(Ri@(matrices(test)+B@l.reshape(3,6))@R,compute_uv=False)[:,0]-rho
    failures=int(np.sum(margin>0)); N=len(test)
    ci=[float(beta.ppf(.025,failures,N-failures+1)) if failures else 0,
        float(beta.ppf(.975,failures+1,N-failures)) if failures<N else 1]
    out=dict(certificate_status="conditional_on_almost_sure_support_reconstruction",
        support_counts="numerically_identified_at_solved_instances",
        seed=seed,states=6,inputs=3,decision_variables=18,dt=dt,rho=rho,
        common_gain_vertex_norm=bound,vertex_strict_margin=rho-bound,
        nmax=nmax,selected_n=n,selected_k=ks,certificate=U,target=.05,
        selected_gain=l.reshape(3,6).tolist(),X=X.tolist(),common_gain=L0.tolist(),desired_gain=desired.reshape(3,6).tolist(),
        evaluation=dict(event="contraction_LMI_violation",validation_seed=seed+1,
                        N=N,failures=failures,estimate=failures/N,CP95=ci),
        boundary_updates=solves,max_k=max(a[1] for a in audits),
        worst_reconstruction_error=max(a[3] for a in audits),
        worst_KKT_residual=max(a[4] for a in audits),
        min_active_gradient_singular_value=min(a[5] for a in audits),
        min_positive_multiplier=min(a[6] for a in audits),
        elapsed_seconds=time.perf_counter()-start)
    np.savetxt(ROOT/'data/mechanical_path.csv',path,delimiter=',',header='n,k,U,objective,'+','.join(f'L{i}' for i in range(18)),comments='')
    np.savetxt(ROOT/'checks/mechanical_audit.csv',audits,delimiter=',',header='n,k,min_slack,reconstruction_error,KKT_residual,min_singular_value,min_multiplier',comments='')
    (ROOT/'data/mechanical_results.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__': run()
