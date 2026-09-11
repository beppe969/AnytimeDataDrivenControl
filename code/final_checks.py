"""Final numerical consistency checks, separate from mathematical proofs."""
from pathlib import Path
import json, math, re
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import beta
from certificates import AnytimeCertificate,log_envelope
from prior_refinement import RefinedPriorCertificate
from monitoring_audit import cg23_beta
ROOT=Path(__file__).resolve().parents[1]

def main():
    checks=[]
    for prior,gamma in [('power',.25),('power',.5),('power',1.),('critical',.5)]:
        c=RefinedPriorCertificate(prior=prior,gamma=gamma)
        J=100000;w=np.array([c.epoch_weight(j) for j in range(J)])
        err=abs(w.sum()+c.tail(J)-1.)
        assert w.min()>0 and err<5e-14
        for n in [10,1000,10**6]:
            for k in [0,1,2,5]:
                j,m=c.epoch(n,k);a=c.edges[k];fac=1+k*math.log(a[j+1]/a[j])
                eta=c.delta*c.weight(k)*c.epoch_weight(j)/fac
                assert log_envelope(m,k,c(n,k))<=math.log(eta)+1e-9
        checks.append({'prior':prior,'gamma':gamma,'telescoping_residual':err})
    # Default and power gamma=1 must be identical up to floating-point evaluation.
    d=AnytimeCertificate();g=RefinedPriorCertificate(prior='power',gamma=1.)
    assert max(abs(d(n,k)-g(n,k)) for n in [20,100,10000] for k in [0,1,2,5])<1e-13
    root=[]
    for k in [1,2,5]:
        f=lambda x:.05*quad(lambda t:t**k*math.exp(x*(1-t)),0.,1.,epsabs=1e-11)[0]-1
        x=brentq(f,0,100)
        nn=[1000,10000,100000];values=[n*cg23_beta(n,k,.05) for n in nn]
        assert abs(values[-1]-x)<abs(values[0]-x)
        root.append(dict(k=k,limit=x,n=nn,n_bound=values))
    # Independently reconstruct reported counts from per-stream records.
    sel=pd.read_csv(ROOT/'data/stopping_audit_replicates.csv');summary=json.loads((ROOT/'data/final_selection_results.json').read_text())
    for r in summary['stopping']:
        z=sel[(sel['eps']==r['eps'])&(sel.method==r['method'])]
        assert len(z)==10000 and int(z.risk_exceeds_eps.sum())==r['failures']
        assert abs(z.acquired.mean()-r['mean_n'])<1e-12
    mon=pd.read_csv(ROOT/'data/monitoring_replicates.csv');ms=json.loads((ROOT/'data/final_monitoring_results.json').read_text())
    for r in ms['coverage']:
        z=mon['first_CG23_failure' if r['method']=='CG23' else 'first_anytime_failure']
        assert int(((z>0)&(z<=r['horizon'])).sum())==r['crossings']
    c=AnytimeCertificate(.05,d=2)
    default_monotonic=[]
    for k in range(3):
        u=np.array([c(n,k) for n in range(max(k,1),8001)])
        assert np.max(np.diff(u))<=1e-13
        default_monotonic.append(k)
    ar=pd.read_csv(ROOT/'data/archive_audit_replicates.csv');a=ar[ar.eps==.05]
    assert int((a.recrossing_intervals>0).sum())==498 and int(a.unavailable_prefixes.sum())==23641
    x=pd.read_csv(ROOT/'data/recrossing_example.csv')
    assert list(x[(x.n>=351)&(x.n<=361)].k.unique())==[2]
    assert x.loc[x.n==350,'k'].item()==1
    out={'alternative_priors':checks,'CG23_fixed_k_limits':root,'reported_counts_match_raw_outputs':True,'fixed_k_certificate_nonincreasing_through_8000':default_monotonic,'archive_sample_counts':{'streams':498,'intervals':499,'prefixes':23641}}
    (ROOT/'checks/final_math_and_data_checks.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
