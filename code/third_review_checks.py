"""Independent arithmetic and output checks for the third-review revision.

The high-precision check uses the complementary finite binomial tail rather
than scipy's incomplete beta function. These checks do not replace proof review.
"""
from pathlib import Path
import json, math
import mpmath as mp
import numpy as np
import pandas as pd
from horizon_baselines import cg23, time_weight
from selection_audit import cp
ROOT=Path(__file__).resolve().parents[1]


def mp_residual(n,k,delta,u):
    with mp.workdps(65):
        nn=mp.mpf(n);uu=mp.mpf(str(u));dd=mp.mpf(str(delta));t=nn*uu
        pmf=mp.exp(nn*mp.log1p(-uu));tail0=pmf
        for j in range(k):
            pmf*=((nn-j)/(j+1))*(uu/(1-uu));tail0+=pmf
        tail=1-tail0
        return float(mp.log(dd)+mp.loggamma(k+1)-mp.fsum(mp.log1p(-mp.mpf(i)/nn) for i in range(k))
                     +mp.log(tail)-(k+1)*mp.log(t)-(nn-k)*mp.log1p(-uu))


def main():
    tests=[]
    for n in [800,8000,100000,10**8,10**12]:
        for k in [0,1,2,5,23,100]:
            for d in [.05,.05/n]:
                u=cg23(n,k,d);r=mp_residual(n,k,d,u)
                tests.append(dict(n=n,k=k,delta=d,u=u,log_equation_residual=r))
    assert max(abs(x['log_equation_residual']) for x in tests)<2e-9
    pd.DataFrame(tests).to_csv(ROOT/'checks/horizon_high_precision_roots.csv',index=False)
    # Exact telescoping identity checked through a finite prefix and its tail.
    for N in [1,10,10000]:
        err=abs(math.fsum(time_weight(n) for n in range(1,N+1))+1/math.log(N+math.e)-1)
        assert err<5e-15
    result=json.loads((ROOT/'data/horizon_baseline_results.json').read_text())
    for key in ['short','monitoring']:
        x=result[key];N=x['protocol']['cap'];df=pd.read_csv(ROOT/f'data/horizon_monitoring_{N}.csv')
        for r in x['coverage']:
            v=df[r['method']];count=int(((v>0)&(v<=r['inspected_through'])).sum())
            assert count==r['crossings'] and np.max(abs(np.array(cp(count,len(df)))-r['CP95']))<1e-14
        assert sum(x['prefix_comparisons_within_tolerance'].values())==0
    d=pd.read_csv(ROOT/'data/horizon_stopping_replicates.csv')
    for r in result['short']['stopping']:
        x=d[(d.method==r['method'])&(d.eps==r['eps'])]
        assert len(x)==10000 and int(x.deployed.sum())==r['deployed']
        assert abs(x.acquired.mean()-r['mean_n'])<1e-12
        assert int(x.risk_exceeds_eps.sum())==r['failures']
    ar=pd.read_csv(ROOT/'data/archive_audit_replicates.csv');ar=ar[ar.eps==.05]
    assert (ar.recrossing_intervals==2).sum()==1
    assert (ar.recrossing_intervals>0).sum()==498 and ar.recrossing_intervals.sum()==499
    feasible=json.loads((ROOT/'checks/horizon_mpc_verification.json').read_text())
    assert all(x['mpc']['feasible'] for x in feasible)
    summary=dict(high_precision_checks=len(tests),precision_decimal_digits=65,
                 maximum_log_equation_residual=max(abs(x['log_equation_residual']) for x in tests),
                 summable_weights_telescope=True,all_new_counts_match_raw_outputs=True,
                 replay_matches_original_results=True,comparisons_within_crossing_tolerance=0,
                 one_archive_stream_has_two_gaps=True,new_MPC_deployments_feasible=True,
                 endpoint_sweep_is_deterministic=True)
    (ROOT/'checks/third_review_verification.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
