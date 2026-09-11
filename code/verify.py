"""Deterministic implementation checks and finite structural audits.

These tests supplement the mathematical proofs; finite tests do not establish
an all-sample probabilistic theorem. Run: python code/verify.py
"""
from __future__ import annotations
from pathlib import Path
import json
import math
import numpy as np
from scipy.stats import beta
from shapely.geometry import MultiPoint, Point
from certificates import log_choose, log_envelope, fixed_bound, AnytimeCertificate
from experiments import clip, project

ROOT = Path(__file__).resolve().parents[1]


def scalar_checks() -> dict:
    rng = np.random.default_rng(20260906)
    count = 0
    maximum_error = 0.0
    for n in range(2, 45):
        for k in range(n):
            for u in rng.uniform(0.0001, 0.9999, 6):
                reference = min(log_choose(n, k) - log_choose(m, k)
                                + (n-m)*math.log1p(-u) for m in range(k, n+1))
                error = abs(reference-log_envelope(n, k, float(u)))
                maximum_error = max(maximum_error, error)
                assert error < 1e-11
                count += 1
    roots = 0
    for n in (10, 100, 10000, 10**9, 10**12):
        for k in (0, 1, 2, 5):
            for eta in (0.1, 0.001, 1e-8):
                u = fixed_bound(n, k, eta)
                assert log_envelope(n, k, u) <= math.log(eta)+1e-10
                assert u > k/n
                ell = math.log(20*math.sqrt(k+1)/(3*eta))
                relaxed = min(1, (k+math.sqrt(2*k*ell)+2*ell)/n)
                assert u <= relaxed*(1+2e-12)
                roots += 1
    beta_checks = 0
    for k in range(1, 8):
        for n in (k, k+1, k+10, 200):
            for u in np.linspace(.001, .999, 30):
                assert beta.sf(u, k, n-k+1) <= math.exp(log_envelope(n,k,u))+1e-12
                beta_checks += 1
    for k in range(6):
        for u in (.01, .1, .4, .9):
            vals = [log_envelope(n,k,u) for n in range(max(1,k), 200)]
            assert np.all(np.diff(vals) <= 1e-12)
    certificate = AnytimeCertificate(.05)
    for n in (10, 100, 1000, 10000):
        for k in range(min(n, 5)):
            j, m = certificate.epoch(n,k)
            assert m <= n
            assert certificate.edges[k][j] <= n < certificate.edges[k][j+1]
            factor = 1+k*math.log1p(1/(k+1+math.log(j+math.e)))
            eta = .05*certificate.weight(k)*certificate.epoch_weight(j)/factor
            assert log_envelope(m,k,certificate(n,k)) <= math.log(eta)+1e-10
    return dict(brute_force_comparisons=count, maximum_log_envelope_error=maximum_error,
                root_checks=roots, beta_tail_checks=beta_checks,
                envelope_monotonicity=True, epoch_checks=True)


def submasks(mask: int):
    sub = mask
    while True:
        yield sub
        if sub == 0:
            break
        sub = (sub-1) & mask


def geometry_checks() -> dict:
    rng = np.random.default_rng(26090603)
    w = rng.normal(size=(9,2))
    n = len(w)
    hulls = {}; hull_boundaries = {}; boxes = {}; box_boundaries = {}
    ids_by_mask = {}
    for mask in range(1 << n):
        ids = [i for i in range(n) if mask & (1 << i)]
        ids_by_mask[mask] = ids
        points = np.vstack([np.zeros((1,2)), w[ids]])
        hulls[mask] = MultiPoint(points).convex_hull
        if hulls[mask].geom_type == 'Polygon':
            vertices = np.asarray(hulls[mask].exterior.coords)
        else:
            vertices = np.asarray(hulls[mask].coords)
        hull_boundaries[mask] = {i for i in ids if np.min(np.linalg.norm(vertices-w[i],axis=1))<1e-12}
        lo, hi = points.min(axis=0), points.max(axis=0)
        boxes[mask] = (lo, hi)
        box_boundaries[mask] = {i for i in ids if np.any(np.abs(w[i]-lo)<1e-12)
                               or np.any(np.abs(w[i]-hi)<1e-12)}
    pairs = 0
    for full in range(1 << n):
        for retained in submasks(full):
            omitted = set(ids_by_mask[full])-set(ids_by_mask[retained])
            retained_ids = set(ids_by_mask[retained])
            accepted = all(hulls[retained].distance(Point(w[i]))<1e-10 for i in omitted)
            keeps = hull_boundaries[full].issubset(retained_ids)
            assert accepted == keeps
            if keeps:
                assert hull_boundaries[retained] == hull_boundaries[full]
            lo,hi = boxes[retained]
            accepted_box = all(np.all(w[i]>=lo-1e-12) and np.all(w[i]<=hi+1e-12)
                               for i in omitted)
            keeps_box = box_boundaries[full].issubset(retained_ids)
            assert accepted_box == keeps_box
            if keeps_box:
                assert box_boundaries[retained] == box_boundaries[full]
            pairs += 1
    return dict(observations=n, retained_omitted_pairs_per_scheme=pairs,
                hull_projectivity=True, box_projectivity=True)


def lmi_checks() -> dict:
    rng = np.random.default_rng(26090604)
    n = 8
    a = rng.uniform(-.7,.7,n); angles=rng.uniform(0,2*np.pi,n)
    b = np.c_[np.cos(angles),np.sin(angles)]
    target = np.array([1.1,.7]); gains={}; ids_by_mask={}; boundaries={}
    for mask in range(1 << n):
        ids = [i for i in range(n) if mask & (1 << i)]
        ids_by_mask[mask] = ids
        if not ids or np.max(np.abs(a[ids]+b[ids]@target))<=.9:
            gains[mask]=target.copy()
        else:
            poly = np.array([[-2.,-2.],[2.,-2.],[2.,2.],[-2.,2.]])
            for i in ids:
                poly=clip(poly,b[i],.9-a[i]); poly=clip(poly,-b[i],.9+a[i])
            gains[mask]=project(poly,target)
    for mask in range(1 << n):
        boundaries[mask]={i for i in ids_by_mask[mask]
                          if np.max(np.abs(gains[mask]-gains[mask^(1 << i)]))>1e-9}
    pairs=0; error=0.
    for full in range(1 << n):
        bmask=sum(1 << i for i in boundaries[full])
        error=max(error,float(np.max(np.abs(gains[bmask]-gains[full]))))
        assert np.allclose(gains[bmask],gains[full],rtol=0,atol=1e-9)
        for retained in submasks(full):
            omitted=sorted(set(ids_by_mask[full])-set(ids_by_mask[retained]))
            accepted=(not omitted) or bool(np.all(np.abs(a[omitted]+b[omitted]@gains[retained])<=.9+1e-10))
            keeps=boundaries[full].issubset(ids_by_mask[retained])
            assert accepted == keeps
            if keeps: assert boundaries[retained] == boundaries[full]
            pairs+=1
    return dict(observations=n, retained_omitted_pairs=pairs, projectivity=True,
                maximum_reconstruction_error=error,
                maximum_support_size=max(map(len,boundaries.values())))


if __name__ == '__main__':
    results=dict(scalar=scalar_checks(), geometry=geometry_checks(), lmi=lmi_checks())
    output=ROOT/'checks'/'verification.json'
    output.write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))
