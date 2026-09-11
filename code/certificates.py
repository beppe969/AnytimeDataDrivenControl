"""Anytime projective-boundary certificates. Python 3.10+, NumPy/SciPy.

All schedules and complexity weights are fixed before data inspection.
The root returned by fixed_bound is the upper endpoint of a bisection bracket.
"""
from __future__ import annotations
from bisect import bisect_right
from functools import lru_cache
import math
import numpy as np
from scipy.special import gammaln, logsumexp
from scipy.stats import binom


def log_choose(n: int, k: int) -> float:
    if not 0 <= k <= n:
        return -math.inf
    k = min(k, n-k)
    if k == 0:
        return 0.0
    return math.fsum(math.log(n-i) for i in range(k)) - math.lgamma(k+1)


def log_envelope(n: int, k: int, u: float) -> float:
    """Log Q_{n,k}(u); for n>k its nontrivial part is strictly decreasing."""
    if n < 1 or not 0 <= k <= n or not 0 <= u <= 1:
        raise ValueError('Require n>=1, 0<=k<=n, and 0<=u<=1.')
    if k == n or u <= k/n:
        return 0.0
    if u == 1:
        return -math.inf
    if k == 0:
        return n*math.log1p(-u)
    m = min(n-1, max(k, int(k/u)))
    ratio = math.fsum(math.log((n-i)/(m-i)) for i in range(k))
    return min(0.0, ratio+(n-m)*math.log1p(-u))


@lru_cache(maxsize=300000)
def fixed_bound(n: int, k: int, eta: float) -> float:
    """Root Q_{n,k}(u)=eta, with the convention b(n,n,eta)=1."""
    if n < 1 or not 0 <= k <= n or not 0 < eta < 1:
        raise ValueError('Invalid n, k, or eta.')
    if k == n:
        return 1.0
    if k == 0:
        return -math.expm1(math.log(eta)/n)
    target = math.log(eta)
    ell = -target + math.log((20/3)*math.sqrt(k+1))
    lo = k/n
    hi = min(1.0, (k+math.sqrt(2*k*ell)+2*ell)/n)
    while log_envelope(n,k,hi) > target:
        hi = min(1.0,2*hi)
    for _ in range(47):
        mid = (lo+hi)/2
        if log_envelope(n,k,mid) > target:
            lo = mid
        else:
            hi = mid
    return min(1.0, np.nextafter(hi, 1.0))


class AnytimeCertificate:
    """Revised shrinking epochs with exact epoch cost; optional ceiling d.

    d=None: pi_k=6/[pi^2(k+1)^2].  Finite d: pi_k=1/(d+1).
    """
    def __init__(self, delta: float = .05, d: int | None = None):
        if not 0 < delta < 1 or (d is not None and d < 0):
            raise ValueError('Invalid confidence or dimension.')
        self.delta, self.d = float(delta), d
        self.edges: dict[int,list[float]] = {}
        self.values: dict[tuple[int,int],float] = {}

    def weight(self,k: int) -> float:
        if self.d is not None:
            if not 0 <= k <= self.d:
                raise ValueError('Observed k exceeds the asserted ceiling.')
            return 1/(self.d+1)
        return 6/(math.pi**2*(k+1)**2)

    @staticmethod
    def epoch_weight(j: int) -> float:
        x = j+math.e
        return math.log1p(1/x)/(math.log(x)*math.log(x+1))

    def epoch(self,n: int,k: int) -> tuple[int,int]:
        if n <= k:
            return -1,n
        a = self.edges.setdefault(k,[float(k+1)])
        while a[-1] <= n:
            j=len(a)-1
            a.append(a[-1]*(1+1/(k+1+math.log(j+math.e))))
        j=bisect_right(a,n)-1
        return j,math.ceil(a[j])

    def __call__(self,n: int,k: int) -> float:
        if n < 1 or not 0 <= k <= n:
            raise ValueError('Invalid n or k.')
        if n == k:
            return 1.0
        j,m=self.epoch(n,k)
        key=(j,k)
        if key not in self.values:
            factor=1+k*math.log1p(1/(k+1+math.log(j+math.e)))
            eta=self.delta*self.weight(k)*self.epoch_weight(j)/factor
            self.values[key]=fixed_bound(m,k,eta)
        return self.values[key]

    def per_time(self,n: int,k: int) -> float:
        eta=self.delta*self.weight(k)*6/(math.pi**2*n*n)
        return fixed_bound(n,k,eta)

    def single_time(self,n: int,k: int) -> float:
        return fixed_bound(n,k,self.delta*self.weight(k))


def classical_sample_size(d: int, eps: float, delta: float) -> int:
    """Smallest n>=d with P[Bin(n,eps)<=d-1]<=delta."""
    if d < 1 or not 0 < eps < 1 or not 0 < delta < 1:
        raise ValueError('Invalid parameters.')
    lo,hi=d-1,max(d,2*d)
    while binom.cdf(d-1,hi,eps)>delta:
        hi*=2
    while hi-lo>1:
        mid=(lo+hi)//2
        if binom.cdf(d-1,mid,eps)>delta:
            lo=mid
        else:
            hi=mid
    return hi


def gcc_stages(d: int, eps: float, delta: float) -> list[int]:
    """Garatti--Care--Campi (2023), Theorem 1, equation (6).

    This is the explicitly identified basic stage construction, not their
    sharper Theorem 4 construction.
    """
    out=[]
    for j in range(d+1):
        m0=classical_sample_size(max(1,j),eps,delta)
        ms=np.arange(j,m0+1)
        terms=gammaln(ms+1)-gammaln(j+1)-gammaln(ms-j+1)
        terms+=(ms-j)*math.log1p(-eps)
        target=math.log(delta)-math.log(d+1)-math.log(m0+1)+logsumexp(terms)
        def f(n):
            return log_choose(n,j)+(n-j)*math.log1p(-eps)
        lo,hi=m0-1,m0
        while f(hi)>target:
            hi*=2
        while hi-lo>1:
            mid=(lo+hi)//2
            if f(mid)>target:
                lo=mid
            else:
                hi=mid
        out.append(hi)
    if any(out[i]>out[i+1] for i in range(d)):
        raise RuntimeError('Stage sizes unexpectedly nonmonotone.')
    return out


if __name__=='__main__':
    rng=np.random.default_rng(7351)
    for n in range(2,40):
        for k in range(n):
            for u in rng.uniform(.001,.999,8):
                brute=min(log_choose(n,k)-log_choose(m,k)+(n-m)*math.log1p(-u)
                          for m in range(k,n+1))
                assert abs(brute-log_envelope(n,k,u))<3e-12
    for n in [10,100,10000]:
        for k in [0,1,min(5,n-1)]:
            for eta in [.1,.001,1e-8]:
                u=fixed_bound(n,k,eta)
                assert log_envelope(n,k,u)<=math.log(eta)+2e-10
    c=AnytimeCertificate()
    print('Certificate tests passed; GCC stages:',gcc_stages(2,.05,.05))
    print('n,k,anytime,per-time:',[(n,2,c(n,2),c.per_time(n,2)) for n in [100,1000,10000]])


def gcc_refined_stages(d: int, eps: float, delta: float) -> list[int]:
    """GCC (2023), Theorem 4 / Algorithm 2, exactly its downward recursion.

    An error is raised on the paper's mu>=1 abort condition. The basic stages
    remain available as an explicitly labelled fallback; never silently used.
    """
    if d<1 or not 0<eps<1 or not 0<delta<1:
        raise ValueError('Require d>=1 and eps,delta in (0,1).')
    bars=[classical_sample_size(max(1,j),eps,delta) for j in range(d+1)]
    logsurv=math.log1p(-eps)
    def log_sum(k,lo,hi):
        if hi<lo: return -math.inf
        m=np.arange(max(k,lo),hi+1)
        return float(logsumexp(gammaln(m+1)-gammaln(k+1)-gammaln(m-k+1)+(m-k)*logsurv))
    def choose_stage(k,lam):
        target=math.log(lam)+log_sum(k,k,bars[k])
        def f(n): return log_choose(n,k)+(n-k)*logsurv
        lo,hi=bars[k]-1,bars[k]
        while f(hi)>target: hi*=2
        while hi-lo>1:
            mid=(lo+hi)//2
            if f(mid)>target: lo=mid
            else: hi=mid
        return hi
    out=[0]*(d+1)
    base=delta/(bars[d]+1)
    out[d]=choose_stage(d,base)
    for k in range(d-1,-1,-1):
        lam=base
        for j in range(d,k,-1):
            a=log_choose(out[j],k)+(out[j]-k)*logsurv
            b=math.log(lam)+log_sum(k,bars[j-1]+1,bars[j])
            if a<=b: mu=0.0
            else:
                lognum=a+math.log(-math.expm1(b-a))
                logden=math.log(lam)+log_sum(k,k,bars[j-1])
                mu=math.exp(lognum-logden)
            if mu>=1:
                raise RuntimeError(f'GCC Algorithm 2 abort: mu[{k},{j}]={mu}')
            lam*=1-mu
        out[k]=choose_stage(k,lam)
    return list(np.maximum.accumulate(out).astype(int))


def gc_upper(n: int,k: int,delta: float) -> float:
    """Garatti--Campi (2022), Theorem 1 Eq.(4), upper risk root.

    Also the upper endpoint of Campi--Garatti JMLR (2023), Theorem 7.
    Uses the published polynomial normalized by binom(n,k)*(1-u)^(n-k).
    Returns the conservative (upper) bisection endpoint. The bound is fixed-n.
    """
    if not n>=1 or not 0<=k<=n or not 0<delta<1:
        raise ValueError('Invalid n,k,delta.')
    if k==n: return 1.0
    m1=np.arange(k,n,dtype=float); m2=np.arange(n+1,4*n+1,dtype=float)
    def ratios(ms):
        # Products avoid cancellation in log-gamma at large n and small k.
        out=np.zeros_like(ms)
        for i in range(k): out+=np.log((ms-i)/(n-i))
        return out
    l1=ratios(m1); l2=ratios(m2)
    def residual(u):
        if u>=1: return math.inf
        lt=math.log1p(-u)
        return float(np.logaddexp(math.log(delta/(2*n))+logsumexp(l1+(m1-n)*lt),
                       math.log(delta/(6*n))+logsumexp(l2+(m2-n)*lt)))
    lo=k/n
    if residual(lo)>=0:
        raise RuntimeError('Expected upper root above k/n.')
    hi=min(1.0,(k+2*math.sqrt(k+1)*(math.sqrt(math.log(k+1))+4+math.sqrt(-math.log(delta)))-math.log(delta))/n)
    while residual(hi)<0: hi=min(1.,2*hi)
    for _ in range(48):
        mid=(lo+hi)/2
        if residual(mid)<0: lo=mid
        else: hi=mid
    return float(np.nextafter(hi,1.))
