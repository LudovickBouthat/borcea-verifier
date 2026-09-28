"""Tests from the endpoint identity (Proposition 5.8 and Corollary 5.9).

Notation (Section 5.5): D(t) = 1 - a t mu and omega_j = a(mu - q_j) as in (3.6),
e_j(omega) the elementary symmetric polynomials of the omega_j, with
|e_j(omega)| <= alpha_j (A v)^{j/2} (Lemma 5.6), J_j(d) as defined before
Lemma 5.6, and, as defined before Proposition 5.8,

    C_T = n/m + A((1 - T)^2 + (1 - 2T)/m),
    B_j = (-1)^j (a mu)^{m-1-j} (m - j - m a mu)   (0 <= j < m),   B_m = (-1)^{m+1} m.

The dual test encloses the product w = a mu directly (w is a local name only):
Re(w) = b/2, |w| = sqrt(Ar), and |Im(w)| = a |Im(mu)| <= a L by (4).  This is
the enclosure of mu in Remark 5.2 multiplied by a = sqrt(A).
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb
from .interval import I, C, imin, square, SCALE


@lru_cache(None)
def alpha(m):
    """Coefficient constants alpha_j, 2 <= j <= m, of (5.8) and Lemma 5.6 (rounded up)."""
    out = {2: I(F(m, 2))}
    for j in range(3, m + 1):
        out[j] = I(1) if j == m else (I(F(m, j)).halfpower(j) * I(F(m, m - j)).halfpower(m - j)).upper()
    return out


SHIFTED_ENDPOINTS = (F(1, 2), F(5, 8), F(3, 8), F(3, 4), F(7, 8), F(1, 4))


def shifted_test(n, A, s, r, v, b):
    """Corollary 5.9, inequality (5.14): a region is excluded if, for some T,
        sqrt(Ar) T (s C_T)^{m/2} + d_T^n + n sqrt(Ar) sum_j alpha_j (Av)^{j/2} T^{j+1} J_j(d_T) < 1,
    with d_T = min{ sqrt(1 - bT + A r T^2), 1 - T + T sqrt(k v) } >= |D(T)|.
    Returns the endpoint T used, or None."""
    m = n - 1
    k = I(F(m - 1, m))
    alphas = alpha(m)
    root = (A * r).sqrt().upper()
    sqav = (A * v).upper().sqrt()
    hp = [I(1)]
    for _ in range(m):
        hp.append(hp[-1] * sqav)
    for t in SHIFTED_ENDPOINTS:
        CT = I(F(n, m)) + A * ((1 - t) ** 2 + F(1 - 2 * t, m))
        main = root * t * (s * CT.positive()).upper().halfpower(m)
        d = imin((1 - b * t + A * r * t * t).positive().sqrt(),
                 1 - t + t * (k * v).positive().sqrt()).upper()
        dp = [I(1)]
        for _ in range(m):
            dp.append(dp[-1] * d)
        rem = I(0)
        for j in range(2, m + 1):
            J = sum((comb(j + l, j) * dp[l] for l in range(m - j + 1)), I(0)) / (n * comb(m, j))
            rem += alphas[j] * hp[j] * (t ** (j + 1)) * J
        if (main + d ** n + n * root * rem).hi < I(1).lo:
            return str(t)
    return None


def dual_test(n, A, s, r, v, b, imag_cap):
    """Proposition 5.8, inequality (5.11): a region is excluded if, for some T and
    some complex multiplier lambda,
        |T J_0(D(T)) + lambda B_0| > T (s C_T)^{m/2}/n
                                     + sum_{j>=2} alpha_j (Av)^{j/2} |T^{j+1} J_j(D(T)) + lambda B_j|.
    w = a mu is enclosed by Re(w) in [b/2, sqrt(Ar)] and
    |Im(w)| <= min{ sqrt(Ar - Re(w)^2), imag_cap }, where imag_cap = a L comes from (4).
    Multipliers are proposed in floating point and then used as exact binary
    rationals; since (5.11) holds for every lambda, their provenance is irrelevant.
    Returns the endpoint T used, or None."""
    m = n - 1
    alphas = alpha(m)
    ar = A * r
    w_re = I.raw(b.lo // 2, ar.sqrt().hi)
    w_im = (ar - square(w_re)).positive().sqrt()
    w_im = I.raw(0, min(w_im.hi, imag_cap.hi))
    w = C(w_re, w_im)
    wp = [C(1)]
    for _ in range(m):
        wp.append(wp[-1] * w)
    B = [((-1) ** j) * ((m - j) * wp[m - 1 - j] - m * wp[m - j]) for j in range(m)]
    B.append(C((-1) ** (m + 1) * m))
    h = (A * v).upper().sqrt()
    hp = [I(1)]
    for _ in range(m):
        hp.append(hp[-1] * h)
    midw = complex((w_re.lo + w_re.hi) / (2 * SCALE), (w_im.lo + w_im.hi) / (2 * SCALE))
    midB = [((-1) ** j) * ((m - j) * midw ** (m - 1 - j) - m * midw ** (m - j)) for j in range(m)]
    midB.append(complex((-1) ** (m + 1) * m))
    if n == 4:
        endpoints = (F(1),)
    elif n in (6, 7):
        endpoints = (F(1), F(3, 4), F(7, 8), F(1, 2))
    else:
        endpoints = (F(1), F(3, 4), F(1, 2))
    for t in endpoints:
        d = 1 - C(I(t)) * w                    # D(T)
        dp = [C(1)]
        for _ in range(m):
            dp.append(dp[-1] * d)
        J = [sum((comb(j + l, j) * dp[l] for l in range(m - j + 1)), C(0)) * C(I(t ** (j + 1) / F(n * comb(m, j))))
             for j in range(m + 1)]
        CT = I(F(n, m)) + A * ((1 - t) ** 2 + F(1 - 2 * t, m))
        rhs = I(t) * ((s * CT.positive()).upper().halfpower(m)) / n
        midd = 1 - float(t) * midw
        choices = [0j]
        indices = range(2, m + 1) if n in (6, 7) else dict.fromkeys((2, 3, m))
        for j in indices:
            if abs(midB[j]) > 1e-14:
                midJ = float(t) ** (j + 1) * sum(comb(j + l, j) * midd ** l for l in range(m - j + 1)) / (n * comb(m, j))
                choices.append(-midJ / midB[j])
        for lam0 in choices:
            lam = C(lam0)
            lower = (J[0] + lam * B[0]).norm().lo
            upper = rhs + sum((alphas[j] * hp[j] * (J[j] + lam * B[j]).norm() for j in range(2, m + 1)), I(0))
            if lower > upper.hi:
                return str(t)
    return None
