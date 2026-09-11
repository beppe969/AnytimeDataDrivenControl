#!/usr/bin/env python3
"""Small deterministic execution check; it does not reproduce cohort statistics."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import sys


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace',type=Path,required=True)
    args=parser.parse_args();workspace=args.workspace.resolve()
    sys.path.insert(0,str(workspace/'code'))
    import numpy as np
    from certificates import AnytimeCertificate, fixed_bound, log_envelope, classical_sample_size, gcc_refined_stages
    from selection_audit import fast_path, cg23_upper
    from horizon_baselines import cg23
    from experiments import lmi_path
    roots=0
    for n in [10,100,1000]:
        for k in [0,1,2]:
            u=fixed_bound(n,k,.01)
            assert k/n<u<1
            assert log_envelope(n,k,u)<=math.log(.01)+1e-10
            assert abs(cg23(n,k,.05)-cg23_upper(n,k,.05))<5e-12
            roots+=1
    assert classical_sample_size(2,.05,.05)==93
    assert gcc_refined_stages(2,.05,.05)==[92,132,162]
    cert=AnytimeCertificate(.05,d=2)
    assert cert(1,1)==1.
    assert abs(100000*AnytimeCertificate(.05)(100000,2)-21.38364321847749)<1e-10
    n=256;seed=26090800;rng=np.random.default_rng(seed)
    angles=rng.uniform(0,2*np.pi,n);offsets=rng.uniform(-.7,.7,n)
    fast=fast_path(angles,offsets)
    slow=lmi_path(np.random.default_rng(seed),n,check=True)
    gain_error=float(np.max(np.abs(fast[0]-slow['gain'])))
    risk_error=float(np.max(np.abs(fast[2]-slow['risk'])))
    same=bool(np.array_equal(fast[1],slow['k']))
    assert gain_error<2e-10 and risk_error<2e-10 and same
    result={'status':'passed','scope':'small execution check; paper cohorts were not regenerated',
            'root_cases':roots,'qp_prefixes':n,'seed':seed,'maximum_gain_difference':gain_error,
            'maximum_risk_difference':risk_error,'support_counts_identical':same,
            'fixed_sample_size':93,'gcc_refined_stages':[92,132,162]}
    (workspace/'checks/smoke.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
