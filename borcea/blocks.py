"""Two-parameter exclusions for a block of degrees nlo <= n <= nhi, 27 <= n <= 2047 (Section 6.4).

A region is given by an integer degree block and a box in the coordinates
(A, e), e = sqrt(1 - s).  The coordinate r = |mu|^2 is not subdivided: it
ranges over the whole interval allowed by restrictions (1), (3), (5), and
v = s - r.  Every bound is uniform in n in the block: powers of bases in
[0, 1] use their smallest exponent, nonnegative prefactors their maximum,
c_m < 5/3 and c_{m,2} < 3, and m/n, k, L and the factor sqrt(A)((n - A)/m)^{m/2}
of K increase with n.

check(nlo, nhi, al, ah, el, eh, parts, mesh) returns the name of an exclusion
valid on the whole region, or None.
"""
from fractions import Fraction as F
from .interval import I, imin, imax, SCALE
from .quadrature import MESHES
from .restrictions import violates_2, mean_bounds, product_bounds


def check(nlo, nhi, al, ah, el, eh, parts=12, mesh='uniform'):
    ml, mh = nlo - 1, nhi - 1

    # Search space: A <= m (Lemma 2.1).
    if al > mh:
        return 'domain'
    ah = min(ah, F(mh))
    A = I(al, ah)
    e = I(el, eh)
    s = 1 - e * e
    k = I(F(ml - 1, ml), F(mh - 1, mh))

    # Restriction (2).
    if violates_2(nlo, nhi, A, s):
        return 'derivative-product'
    if al <= 0:
        return None

    # Restrictions (1), (3), (5): admissible interval of |mu| = sqrt(r) and a
    # lower bound for b = 2 a Re(mu).
    mb = mean_bounds(A, s, e, k, mh)
    if mb is None:
        return 'restriction-5'
    a, U, b = mb.a, mb.U, mb.b
    As = (A * s).upper()
    if (1 - b + As).hi < 0:                 # beta(1) = 1 - b + A s >= 0
        return 'negative-endpoint'
    if b.lo > (2 * a * U).hi:               # upper bound in (5)
        return 'restriction-5'
    K, K_low = product_bounds(nlo, nhi, al, ah, s)

    def beta_up(t):
        value = 1 - b * t + As * t * t
        if value.hi < 0:
            return None
        return imin(value.positive(), 1).upper()

    # Direct test (Proposition 5.3) with sqrt(r) <= U.
    end = beta_up(I(1))
    if end is None:
        return 'negative-endpoint'
    first_terms = (imin((k * s).upper().halfpower(ml), end.halfpower(ml))
                   + F(mh, nhi) * K * U)
    H = (k * end.halfpower(nlo - 2)).upper()      # k d^{n-2} with d^2 = |D(1)|^2 <= beta(1)
    remainder = I(0)
    segments = []
    for l, h, weights in MESHES[mesh](parts):
        bl, bh = beta_up(l), beta_up(h)
        if bl is None or bh is None:
            return 'negative-interior'
        blp, bhp = bl.halfpower(nlo - 2), bh.halfpower(nlo - 2)
        wl, wh, _ = weights[0]
        part = wl * blp + wh * bhp
        gl = (s * (A * (1 - l) * (1 - l) + k)).upper()
        gh = (s * (A * (1 - h) * (1 - h) + k)).upper()
        if gl.hi <= SCALE and gh.hi <= SCALE:
            part = imin(part, wl * gl.halfpower(nlo - 2) + wh * gh.halfpower(nlo - 2))
        remainder += part
        segments.append((l, h, weights, bl, bh, blp, bhp))
    direct = first_terms + F(5, 3) * mh * As * remainder
    if direct.hi < SCALE:
        return 'direct'

    # Centered test, integral form, and the monotonicity criterion (6.1):
    # with z = sqrt(r) in [0, sqrt(s)], the right side of (5.6) is
    # Psi(z) = K z + (s - z^2)(H + C z); if K >= 2 H sqrt(s) + 2 C s then
    # Psi(z) <= Psi(sqrt(s)) = K sqrt(s) <= s^{n/2} < 1 (Lemma 2.3).
    a_up = a.upper()
    integral = I(0)
    for l, h, weights, bl, bh, blp, bhp in segments:
        wl, wh, mass = weights[1]
        if nlo >= 5:
            second_int = wl * bl.halfpower(nlo - 3) + wh * bh.halfpower(nlo - 3)
        else:
            second_int = mass * imax(bl, bh).halfpower(nlo - 3)
        second = F(3, 2) * mh * (mh - 1) * A.upper() * second_int
        den = 1 - a_up * U * h
        if den.lo > 0:
            first = F(5, 6) * mh * A.upper() * (wl * blp + wh * bhp) / den
            integral += imin(first, second)
        else:
            integral += second
    Cc = (nhi * a_up * integral).upper()
    if K_low.lo >= (2 * H * s.sqrt() + 2 * Cc * s).hi:
        return 'centered-monotone'
    return None
