"""Additional deterministic checks for the revision; no empirical proof claims.

Run from any directory: python code/verify_revision.py
"""
from __future__ import annotations
from pathlib import Path
from itertools import combinations
import math
import json
import numpy as np
from scipy.special import comb
from scipy.stats import binom
from certificates import (AnytimeCertificate, log_envelope, gc_upper,
                          gcc_refined_stages, classical_sample_size, fixed_bound)
from verify import submasks
ROOT=Path(__file__).resolve().parents[1]


def explicit_bound_checks() -> dict:
    cert=AnytimeCertificate(.05); count=0; max_fraction=0.
    for n in (2,5,10,50,500,1000,10000,100000,10**7,10**10):
        for k in sorted(set([0,1,min(2,n-1),min(5,n-1),min(20,n-1),min(64,n-1)])):
            j,m=cert.epoch(n,k)
            factor=1+k*math.log1p(1/(k+1+math.log(j+math.e)))
            eta=.05*cert.weight(k)*cert.epoch_weight(j)/factor
            ell=math.log(20*math.sqrt(k+1)/(3*eta))
            r=math.log(n/(k+1)); H=math.e+2*(k+2)*r; R=2*H*math.log(2*H)
            lam=(math.log(40*math.sqrt(k+1)/(3*.05*cert.weight(k)))
                 +math.log(R+1)+2*math.log(math.log(R+1)))
            assert ell <= lam+1e-12
            u=cert(n,k)
            finite=min(1.,(k+1+(1+1/(k+1))*(math.sqrt(2*k*ell)+2*ell))/n)
            assert u <= finite+2e-14
            assert n*u >= k-1e-12
            assert log_envelope(m,k,u) <= math.log(eta)+1e-10
            max_fraction=max(max_fraction,ell/lam);count+=1
    return {'cases':count,'ell_at_most_Lambda':True,'finite_envelope':True,
            'maximum_ell_over_Lambda':max_fraction}


def exact_order_statistic_checks() -> dict:
    z=np.random.default_rng(26090784).uniform(size=8)
    pairs=0
    for k in (1,2,3,4):
        B={};level={};accepted={}
        for mask in range(1<<len(z)):
            ids=[i for i in range(len(z)) if mask>>i&1]
            order=sorted(ids,key=lambda i:z[i],reverse=True)
            B[mask]=set(order[:k]);level[mask]=(z[order[k-1]] if len(order)>=k else -math.inf)
            accepted[mask]=lambda value, mm=mask: (value<=level[mm] or any(value==z[ii] for ii in B[mm]))
            assert all(accepted[mask](z[ii]) for ii in ids)
        for full in range(1<<len(z)):
            for retained in submasks(full):
                ids={i for i in range(len(z))if retained>>i&1}
                omitted=[i for i in range(len(z))if full>>i&1 and not retained>>i&1]
                assert all(accepted[retained](z[i]) for i in omitted)==B[full].issubset(ids)
                if B[full].issubset(ids):assert B[retained]==B[full]
                pairs+=1
    # Unique optimizer with duplicate constraints: support reconstruction fails.
    opt=lambda xs: max([0.]+xs)
    full=[1.,1.]
    support=[i for i in range(2)if opt(full[:i]+full[i+1:])!=opt(full)]
    assert not support and opt(full)!=opt([])
    return {'retained_omitted_pairs':pairs,'all_k_witnesses_projective':True,'sample_consistency':True,
            'duplicate_constraint_counterexample':True}


def refined_reference(d:int,eps:float,delta:float) -> list[int]:
    """Independent direct-sum implementation, suitable only for small d."""
    bars=[classical_sample_size(max(1,j),eps,delta)for j in range(d+1)]
    def h(m,k):return comb(m,k)*(1-eps)**(m-k)
    def total(k,a,b):return sum(h(m,k)for m in range(max(k,a),b+1))
    def pick(k,lam):
        n=bars[k];rhs=lam*total(k,k,bars[k])
        while h(n,k)>rhs:n+=1
        return n
    N=[0]*(d+1);lam0=delta/(bars[d]+1);N[d]=pick(d,lam0)
    for k in range(d-1,-1,-1):
        lam=lam0
        for j in range(d,k,-1):
            mu=max(0.,(h(N[j],k)-lam*total(k,bars[j-1]+1,bars[j]))
                   /(lam*total(k,k,bars[j-1])))
            if mu>=1:raise RuntimeError('Published algorithm abort condition')
            lam*=1-mu
        N[k]=pick(k,lam)
    return list(np.maximum.accumulate(N).astype(int))


def benchmark_checks() -> dict:
    tested=[]
    for d in (1,2,4,8):
        for eps,delta in ((.05,.05),(.10,.01)):
            a=gcc_refined_stages(d,eps,delta);b=refined_reference(d,eps,delta)
            assert a==b
            tested.append({'d':d,'eps':eps,'delta':delta,'stages':[int(x)for x in a]})
    assert gcc_refined_stages(2,.05,.05)==[92,132,162]
    expected={(1000,2):9.821819537963712,(3425,23):41.3887923217847,
              (10000,23):41.46143823462574}
    gc=[]
    for (n,k),value in expected.items():
        got=n*gc_upper(n,k,.05);assert abs(got-value)<2e-9
        # Independent direct polynomial residual (small-to-medium sample sizes).
        u=got/n; t=1-u
        s1=sum(comb(m,k)/comb(n,k)*t**(m-n) for m in range(k,n))*.05/(2*n)
        s2=sum(comb(m,k)/comb(n,k)*t**(m-n) for m in range(n+1,4*n+1))*.05/(6*n)
        assert abs(s1+s2-1)<2e-9
        gc.append({'n':n,'k':k,'n_upper':got,'root_residual':s1+s2-1})
    return {'refined_stage_cases':tested,'GC22_root_checks':gc}


if __name__=='__main__':
    out={'explicit_bounds':explicit_bound_checks(),
         'order_statistics':exact_order_statistic_checks(),
         'benchmarks':benchmark_checks()}
    (ROOT/'checks'/'revision_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
