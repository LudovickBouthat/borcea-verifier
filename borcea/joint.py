"""Three-parameter exclusions for a single degree 4 <= n <= 26 (Sections 5 and 6).

A box is given in the coordinates (A, e, f) with e = sqrt(1 - s) and
f = sqrt(v/s), so that s = 1 - e^2, r = s(1 - f^2), v = s f^2.  The search
space is [0, n - 1] x [0, 1] x [0, 1].

check(n, box, parts) returns the name of an exclusion that holds on the whole
box, or None.  Names are for information only; the verifier never trusts them.
"""
from fractions import Fraction as F
from math import comb
from .interval import I, imin, SCALE
from .quadrature import uniform_mesh
from .restrictions import (violates_2, contract_s_by_2, deficit_radius,
                           mean_bounds, product_bounds)
from .endpoint import alpha, shifted_test, dual_test

QUARTIC_LOCAL = 2401   # Section 5.6, (5.18): in degree four, 2401 v <= s is impossible.


def check(n, box, parts=12):
    parts = max(parts, 12)
    al, ah, el, eh, fl, fh = box
    m = n - 1
    k = I(F(m - 1, m))

    # Quartic local lemma (Proposition 5.10), recognized through f^2 = v/s.
    if n == 4 and QUARTIC_LOCAL * fh * fh <= 1:
        return 'quartic-local'

    # Search space: A <= m (Lemma 2.1).
    if al > m:
        return 'domain'
    ah = min(ah, F(m))
    A = I(al, ah)
    e = I(el, eh)
    s = 1 - e * e

    # Restriction (2), then the lower bound it gives for s.
    if violates_2(n, n, A, s):
        return 'derivative-product'
    if al <= 0:
        return None
    contracted = contract_s_by_2(n, A, s)
    if contracted is None:
        return 'derivative-product'
    if contracted.lo != s.lo:
        s = contracted
        e = (1 - s).positive().sqrt()

    # Restriction (1) is built into the coordinates: 0 <= v <= s <= 1.
    f = I(fl, fh)
    v = s * f * f
    r = s * (1 - f * f)
    L = deficit_radius(A, s, e, m)

    # Restriction (3): contract A to [(sqrt r - L)^2, (sqrt r + L)^2].
    z = r.sqrt()
    lo = max(A.lo, ((z - L).positive().lower() ** 2).lo)
    hi = min(A.hi, ((z + L).upper() ** 2).hi)
    if lo > hi:
        return 'restriction-3'
    if lo > A.lo or hi < A.hi:
        al, ah = F(lo, SCALE), F(hi, SCALE)
        A = I(al, ah)
        if violates_2(n, n, A, s):
            return 'derivative-product'

    # Restrictions (1), (3), (5): contract r, and v = s - r with it; lower bound for b.
    mb = mean_bounds(A, s, e, k, m, r_box=r)
    if mb is None:
        return 'restriction-5'
    a, L, r, b = mb.a, mb.L, mb.r, mb.b
    vl, vh = max(v.lo, s.lo - r.hi), min(v.hi, s.hi - r.lo)
    if vl > vh:
        return 'restriction-1'
    v = I.raw(vl, vh)
    z = r.sqrt()
    ar, av = A * r, A * v
    root = ar.sqrt()
    # Upper bound in (5): b <= 2 sqrt(Ar).
    if b.lo > (2 * root).hi:
        return 'restriction-5'

    # Restriction (4), |Im mu| <= L, enters only through the dual test, as |Im(a mu)| <= a L.
    imag_cap = a * L
    if n == 4 and dual_test(n, A, s, r, v, b, imag_cap) is not None:
        return 'dual:1'

    # d >= |D(1)|: |D(1)|^2 = 1 - b + Ar = beta(1) - Av <= 1 - Av (Lemma 3.1).
    d2 = imin(1 - b + ar, 1 - av)
    if d2.hi < 0:
        return 'negative-D'
    d = d2.positive().sqrt().upper()
    K, K_low = product_bounds(n, n, al, ah, s)

    # Centered test with coefficient form (Proposition 5.4, Corollary 5.7):
    # 1 <= K sqrt(r) + v T, T = k d^{n-2} + n A sqrt(Ar) sum_j alpha_j (Av)^{(j-2)/2} J_j(d).
    alphas = alpha(m)
    powers = [I(1)]
    for _ in range(m):
        powers.append(powers[-1] * d)
    sqav = av.upper().sqrt()
    avpowers = [I(1)]
    for _ in range(m - 2):
        avpowers.append(avpowers[-1] * sqav)
    coef = I(0)
    for j in range(2, m + 1):
        J = sum((comb(j + l, j) * powers[l] for l in range(m - j + 1)), I(0)) / (n * comb(m, j))
        coef += alphas[j] * avpowers[j - 2] * J
    T = (k * d.halfpower(2 * (n - 2)) + n * root * A * coef).upper()
    if (K * z + v * T).hi < SCALE:
        return 'centered-coefficient'
    if 2 * T.hi <= K_low.lo:
        return 'centered-coefficient-boundary'

    # Shifted endpoint test (Corollary 5.9).
    t = shifted_test(n, A, s, r, v, b)
    if t is not None:
        return 'shifted:' + t

    # Endpoint duality (Proposition 5.8) in degrees 5 to 13.
    if 5 <= n <= 13:
        t = dual_test(n, A, s, r, v, b, imag_cap)
        if t is not None:
            return 'dual:' + t

    # Direct test (Proposition 5.3), integrated by the chord rule.
    mesh = uniform_mesh(parts)
    As = (A * s).upper()

    def beta_up(t):
        # beta(t) <= min{1, 1 - b t + A s t^2} (Lemma 3.1 and Remark 5.2).
        return imin((1 - b * t + As * t * t).positive(), 1).upper()

    remainder = I(0)
    for l, h, weights in mesh:
        wl, wh, _ = weights[0]
        part = wl * beta_up(l).halfpower(m - 1) + wh * beta_up(h).halfpower(m - 1)
        gl = (s * (A * (1 - l) * (1 - l) + k)).upper()
        gh = (s * (A * (1 - h) * (1 - h) + k)).upper()
        if max(gl.hi, gh.hi) <= SCALE:   # gamma <= 1 on the whole subinterval
            part = imin(part, wl * gl.halfpower(m - 1) + wh * gh.halfpower(m - 1))
        remainder += part
    cm = I(F(m, m - 1)).halfpower(m - 1).upper()
    direct = (imin((k * s).upper().halfpower(m), beta_up(I(1)).halfpower(m))
              + F(m, n) * K * z + cm * m * A * s * remainder)
    if direct.hi < SCALE:
        return 'direct'

    # Centered test with integral form (Corollary 5.5), same chord rule.
    c2 = I(F(m, m - 2)).halfpower(m - 2).upper()
    cache = {}
    integral = I(0)
    for l, h, weights in mesh:
        for t in (l, h):
            if (t.lo, t.hi) not in cache:
                bt = beta_up(t)
                cache[t.lo, t.hi] = (bt.halfpower(m - 1) / 2, bt.halfpower(m - 2) / 2)
        vl_, vh_ = cache[l.lo, l.hi], cache[h.lo, h.hi]
        wl, wh, _ = weights[1]
        second = c2 * m * (m - 1) * A * (wl * vl_[1] + wh * vh_[1])
        den = 1 - root * h        # |D(t)| >= 1 - t sqrt(Ar) >= 1 - h sqrt(Ar)
        if den.lo > 0:
            first = cm * m * A * (wl * vl_[0] + wh * vh_[0]) / den
            integral += imin(first, second)
        else:
            integral += second
    T_int = (k * d.halfpower(2 * (n - 2)) + n * root * integral).upper()
    if (K * z + v * T_int).hi < SCALE:
        return 'centered-integral'
    if 2 * T_int.hi <= K_low.lo:
        return 'centered-integral-boundary'
    return None
