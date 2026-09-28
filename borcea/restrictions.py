"""Proposition 5.1, restricted to items (1), (2), (3), (4) and (5).

Notation (Section 5): m = n - 1, k = (m - 1)/m, A = a^2, r = |mu|^2,
v = s - r, b = 2 a Re(mu) and

    K = sqrt(A) (s (n - A)/m)^{m/2},   L = min{ (1 - s) sqrt(m A),  sqrt((1 - s)(A - 1/m)) }.

The only facts about a counterexample used by this module are

    (1)  0 <= v <= s < 1,
    (2)  n s (1 + A) >= n^{2/m} m,
    (3)  |sqrt(A) - sqrt(r)| <= L,
    (4)  |Im mu| <= L,
    (5)  max{A + r - L^2, 1 + A r - k v} <= b <= 2 sqrt(A r),

together with A <= m (Lemma 2.1, the definition of the search space) and
the product bounds of Lemma 2.3, which give Pi |p(0)| <= K <= s^{m/2}.  All
functions accept a block of degrees nlo <= n <= nhi; every bound is made
uniform in n by the monotonicity facts of Section 6.4.
"""
from fractions import Fraction as F
from functools import lru_cache
from .interval import I, imin, imax, SCALE


@lru_cache(None)
def derivative_constant(n):
    """Enclosure of c_n = (m/n) n^{2/m}; its m-th power is found by integer bisection."""
    m = n - 1
    target = n * n * SCALE ** m
    lo, hi = SCALE, n * n * SCALE
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid ** m <= target:
            lo = mid
        else:
            hi = mid
    return I.raw(lo, hi) * F(m, n)


def violates_2(nlo, nhi, A, s):
    """True if (2) fails at every point of the box and every degree of the block.

    (2) reads (n s (1+A)/m)^m >= n^2.  The base is at most nlo s (1+A)/(nlo - 1);
    if it is <= 1 the inequality fails, otherwise the base to the power nhi - 1
    bounds the left side for all degrees of the block (Section 6.4)."""
    base = s * (1 + A) * F(nlo, nlo - 1)
    return base.hi <= SCALE or (base.upper() ** (nhi - 1)).hi < nlo * nlo * SCALE


def contract_s_by_2(n, A, s):
    """Lower bound s >= c_n/(1 + A) from (2), for a single degree n.
    Returns the contracted s, or None when the box is empty."""
    s_min = (derivative_constant(n) / (1 + A)).lo
    if s_min > s.hi:
        return None
    return I.raw(max(s.lo, s_min), s.hi)


def deficit_radius(A, s, e, m):
    """Upper bound for L.  Both entries increase with m, so a block uses m = nhi - 1.
    Here e = sqrt(1 - s)."""
    eps = 1 - s
    return imin(eps * (m * A).sqrt(), e * (A - I(F(1, m))).positive().sqrt())


class MeanBounds:
    """Consequences of (1), (3), (5) for mu on a box (see mean_bounds)."""
    __slots__ = ('a', 'L', 'r', 'U', 'b')

    def __init__(self, a, L, r, U, b):
        self.a, self.L, self.r, self.U, self.b = a, L, r, U, b


def mean_bounds(A, s, e, k, m, r_box=None):
    """Enclose r = |mu|^2 and bound b = 2 a Re(mu) from below using (1), (3), (5).

    Admissible values of x = sqrt(r):
      (1)  x <= sqrt(s);
      (3)  a - L <= x <= a + L;
      (5)  1 + A x^2 - k(s - x^2) <= 2 a x, i.e. (A + k) x^2 - 2 a x + 1 - k s <= 0,
           so x lies between the roots x_-+ = (a -+ sqrt(Delta))/(A + k),
           Delta = k((A + k) s - 1).  Both roots are positive since 1 - k s > 0.
           (The first lower bound in (5) against its upper bound is equivalent to (3).)
    If r_box is given (three-parameter checker), r is also intersected with it.

    Lower bounds for b.  Both lower bounds in (5) increase with r, so they may be
    evaluated at the smallest admissible r.  The following evaluations are used:
      A + r - L^2                    first bound of (5);
      A s + r + (1 - s)/m            the same with L^2 <= (1 - s)(A - 1/m), expanded;
      1 + A r - k(s - r)             second bound of (5), with v = s - r;
      2 A (1 - L/a)                  first bound at r = (a - L)^2, from (3); written
                                     with L/a = min{(1 - s) sqrt(m), e sqrt(1 - 1/(mA))};
      2 a x_-                        second bound at x = x_-, where it equals 2 a x_-.
    Returns None if no admissible r exists, otherwise MeanBounds with the interval
    r, the upper bound U >= |mu| and the lower bound b (a point interval, b >= 0)."""
    a = A.sqrt()
    L = deficit_radius(A, s, e, m)
    delta = k * ((A + k) * s - 1)
    if delta.hi < 0:
        return None
    root = delta.positive().sqrt()
    x_minus = (a - root) / (A + k)
    x_plus = (a + root) / (A + k)
    lo = max((a - L).positive().lo, x_minus.positive().lo)
    hi = min(s.sqrt().hi, (a + L).hi, x_plus.hi)
    if lo > hi:
        return None
    r = I.raw((I.raw(lo, lo) * I.raw(lo, lo)).lo, (I.raw(hi, hi) * I.raw(hi, hi)).hi)
    if r_box is not None:
        rlo, rhi = max(r.lo, r_box.lo), min(r.hi, r_box.hi)
        if rlo > rhi:
            return None
        r = I.raw(rlo, rhi)
        hi = min(hi, r.sqrt().hi)
    eps = 1 - s
    L_rel = imin(eps * I(m).sqrt(), e * (1 - 1 / (m * A)).positive().sqrt())
    b = imax((A + r - L * L).lower(),
             (A * s + r + eps * F(1, m)).lower(),
             (1 + A * r - k * (s - r)).lower(),
             (2 * A * (1 - L_rel)).lower(),
             (2 * a * x_minus).lower()).positive().lower()
    return MeanBounds(a, L, r, I.raw(hi, hi), b)


def product_bounds(nlo, nhi, alo, ahi, s):
    """Enclosures of K = sqrt(A) (s (n - A)/m)^{m/2} = g s^{m/2}, where
    g = sqrt(A)((n - A)/m)^{m/2} <= 1 by Lemma 2.3, so that Pi |p(0)| <= K <= s^{m/2}.

    g increases with n, and for fixed n increases on [0, 1] and decreases after
    (Section 6.4).  Returns (K, K_low): an upper bound for K, and a lower bound
    for the explicit majorant K used only in the centered boundary criteria
    (Proposition 5.4 and (6.1)).  alo, ahi are the exact rational endpoints of
    the A-interval."""
    ml, mh = nlo - 1, nhi - 1
    ap = F(1) if alo <= 1 <= ahi else (ahi if ahi < 1 else alo)
    g = imin(I(ap).sqrt() * I((F(nhi) - ap) / mh).halfpower(mh), 1).upper()
    K = (g * s.upper().halfpower(ml)).upper()
    if ahi >= nlo:
        g_low = I(0)
    else:
        g_low = imin(*[I(x).sqrt() * I((F(nlo) - x) / ml).halfpower(ml) for x in (alo, ahi)]).lower()
    K_low = (g_low * s.lower().halfpower(mh)).lower()
    return K, K_low
