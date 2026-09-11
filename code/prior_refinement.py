"""Alternative prespecified epoch priors; all published main experiments use default.

The critical prior gives the fixed-complexity second-order expansion in Remark 9.
Asymptotic improvements do not imply pointwise improvement at finite sample sizes.
"""
from __future__ import annotations
import math
from certificates import AnytimeCertificate

class RefinedPriorCertificate(AnytimeCertificate):
    def __init__(self, delta: float=.05, d: int|None=None,
                 prior: str='critical', gamma: float=.5):
        super().__init__(delta,d)
        if prior not in ('critical','power') or gamma<=0:
            raise ValueError('Use critical or power prior and positive gamma.')
        self.prior=prior;self.gamma=float(gamma)

    def tail(self,j: int) -> float:
        if j<0:raise ValueError('Nonnegative epoch index required.')
        if self.prior=='power':return math.log(j+math.e)**(-self.gamma)
        return 1/math.log(math.log(j+math.exp(math.e)))

    def epoch_weight(self,j: int) -> float:
        if j<0:raise ValueError('Nonnegative epoch index required.')
        if self.prior=='power':
            x=j+math.e;l=math.log(x);d=math.log1p(1/x)
            return l**(-self.gamma)*(-math.expm1(-self.gamma*math.log1p(d/l)))
        x=j+math.exp(math.e);l=math.log(x);v=math.log(l)
        dv=math.log1p(math.log1p(1/x)/l)
        return dv/(v*(v+dv))
