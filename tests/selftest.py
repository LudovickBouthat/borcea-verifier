"""Self-tests run by verify_all.py before any certificate is replayed.

1. Interval arithmetic encloses exact rational results.
2. Quadrature weights are nonnegative and sum to the exact moments.
3. The coefficient identities behind Lemma 5.6 and Proposition 5.8 hold exactly
   for random rational configurations.
4. The constants of the quartic local lemma (Proposition 5.10, radius 1/20) check.
5. The restriction code is sound: at random points satisfying (1), (2), (3), (5),
   every box containing the point is accepted by the restriction layer and the
   enclosures contain the point's r and b.
6. The verifier fails closed on forged, unresolved or incomplete certificates.
7. A box around a known obstruction is not certified.
"""
import math
import os
import sys
from fractions import Fraction as F
from random import Random
from math import comb

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from borcea.interval import I, C, SCALE
from borcea.quadrature import uniform_mesh, geometric_mesh
from borcea.restrictions import violates_2, mean_bounds, deficit_radius
from borcea import certificate


def contains(iv, x):
    assert F(iv.lo, SCALE) <= x <= F(iv.hi, SCALE), (iv.lo, iv.hi, x)


def test_arithmetic():
    rng = Random(401)
    for _ in range(2000):
        a = F(rng.randrange(-100, 101), rng.randrange(1, 100))
        b = F(rng.randrange(1, 101), rng.randrange(1, 100))
        contains(I(a) + I(b), a + b)
        contains(I(a) - I(b), a - b)
        contains(I(a) * I(b), a * b)
        contains(I(a) / I(b), a / b)
        root = I(b).sqrt()
        assert F(root.lo, SCALE) ** 2 <= b <= F(root.hi, SCALE) ** 2
        power = I(b).halfpower(7)
        assert F(power.lo, SCALE) ** 2 <= b ** 7 <= F(power.hi, SCALE) ** 2


def test_quadrature():
    for mesh in (uniform_mesh, geometric_mesh):
        for q in (8, 12, 48, 192):
            totals = [I(0), I(0)]
            for l, h, weights in mesh(q):
                for j in (0, 1):
                    wl, wh, mass = weights[j]
                    assert wl.lo >= 0 and wh.lo >= 0
                    totals[j] += wl + wh
            contains(totals[0], F(1, 2))
            contains(totals[1], F(1, 3))


def _mul(a, b):
    c = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return c


def test_coefficient_identities():
    rng = Random(19)
    for m in range(3, 13):
        n = m + 1
        for _ in range(20):
            pairs = [(F(rng.randrange(6, 9), 10), F(rng.randrange(1, 3), 10)) for _ in range(m // 2)]
            lone = F(4, 5) if m % 2 else None
            Q = (sum(2 * x for x, y in pairs) + (lone or 0)) / m
            a = (sum(2 * x / (x * x + y * y) for x, y in pairs) + (1 / lone if lone else 0)) / m
            s = (sum(2 * (x * x + y * y) for x, y in pairs) + (lone * lone if lone else 0)) / m
            v = s - Q * Q
            w = a * Q                      # a mu
            eq, ey = [F(1)], [F(1)]
            for x, y in pairs:
                eq = _mul(eq, [1, 2 * x, x * x + y * y])
                z = x - Q
                ey = _mul(ey, [1, 2 * z, z * z + y * y])
            if lone:
                eq = _mul(eq, [1, lone])
                ey = _mul(ey, [1, lone - Q])
            assert ey[1] == 0 and eq[m - 1] == m * a * eq[m]
            E = [(-a) ** j * ey[j] for j in range(m + 1)]   # e_j(omega), omega_j = a(mu - q_j)
            B = [(-1) ** j * ((m - j) * w ** (m - 1 - j) - m * w ** (m - j)) for j in range(m)] + [F((-1) ** (m + 1) * m)]
            assert sum(E[j] * B[j] for j in range(m + 1)) == 0
            for T in (F(1), F(1, 2), F(5, 8), F(7, 8)):
                d = 1 - T * w              # D(T)
                J = [sum(F(comb(j + l, j)) * d ** l for l in range(m - j + 1)) / (n * comb(m, j)) for j in range(m + 1)]
                integral = sum((-a) ** j * eq[j] * T ** (j + 1) / (j + 1) for j in range(m + 1))
                assert integral == sum(E[j] * T ** (j + 1) * J[j] for j in range(m + 1))
            for j in range(2, m + 1):
                csq = F(m * m, 4) if j == 2 else (F(m, j) ** j * F(m, m - j) ** (m - j) if j < m else F(1))
                assert ey[j] ** 2 <= csq * v ** j


def test_quartic_constants():
    # Proposition 5.10 with |zeta_j| <= 1/20 after normalization a = 1, b_p = 0.
    rho = F(1, 20)
    tau = 3 * rho * rho
    assert 4 * rho == F(1, 5)                        # |4 z1 z2 z3| <= tau/10
    assert 4 * rho / (1 - rho) == F(4, 19)           # 4|z|^3/(1-|z|) <= (4/19)|z|^2
    assert F(1, 3) - F(4, 57) >= F(1, 4)             # S >= 1 + 2R/3 + tau/4
    assert 1 - tau - tau / 10 > 0                    # the lower bound for V^2 is positive
    # (1 - R - tau/10)(1 + 4R/3 + tau/2) >= 1 for |R| <= tau <= 3/400:
    # = 1 + R/3 + 2tau/5 - 4R^2/3 - 19 R tau/30 - tau^2/20.
    worst = -F(1, 3) + F(2, 5) - F(4, 3) * tau - F(19, 30) * tau - tau / 20
    assert worst > 0
    # Recognition: 2401 v <= s gives 6 v/r <= 1/400, i.e. max|zeta_j|/a <= 1/20.
    assert F(6, 2400) == F(1, 400)


def test_restriction_soundness():
    """Random points satisfying (1), (2), (3), (5) must survive the restriction layer."""
    rng = Random(7)
    tested = 0
    while tested < 1500:
        n = rng.choice([4, 5, 6, 9, 14, 26, 27, 40, 300, 2047])
        m = n - 1
        k = F(m - 1, m)
        A = F(rng.randrange(1, 10 ** 6), 10 ** 6) * rng.choice([1, 2, 4, 8])
        e = F(rng.randrange(1, 10 ** 6), 10 ** 6)
        s = 1 - e * e
        if A > m or not (n * s * (1 + A) / m) ** m >= n * n:
            continue
        f = F(rng.randrange(0, 10 ** 6), 10 ** 6)
        v, r = s * f * f, s * (1 - f * f)
        a, z = math.sqrt(A), math.sqrt(r)
        L = min(float(1 - s) * math.sqrt(m * A), math.sqrt(max(0.0, float((1 - s) * (A - F(1, m))))))
        lower = max(float(A + r) - L * L, float(1 + A * r - k * v))
        upper = 2 * math.sqrt(float(A * r))
        if abs(a - z) > L * (1 - 1e-9) or lower > upper * (1 - 1e-9):
            continue
        b = (lower + upper) / 2
        w = F(1, 2 ** rng.randrange(8, 30))
        Abox, ebox = I(A - w, A + w).positive(), I(max(F(0), e - w), min(F(1), e + w))
        if F(Abox.lo, SCALE) <= 0:
            continue
        sbox = 1 - ebox * ebox
        assert not violates_2(n, n, Abox, sbox)
        kk = I(k)
        mb = mean_bounds(Abox, sbox, ebox, kk, m, r_box=I(r - w, r + w).positive())
        assert mb is not None
        contains(mb.r, r)
        assert F(mb.b.lo, SCALE) <= F(b) + F(1, 10 ** 12)
        assert F(mb.U.hi, SCALE) >= F(z) - F(1, 10 ** 12)
        assert F((mb.a * mb.L).hi, SCALE) >= F(a * L) - F(1, 10 ** 12)
        tested += 1


def test_fail_closed():
    import io
    from contextlib import redirect_stdout
    bad = [
        {'format': certificate.JOINT_FORMAT, 'degrees': [{'n': 4, 'tree': {'unresolved': True}}]},
        {'format': certificate.JOINT_FORMAT, 'degrees': [{'n': n, 'tree': {'test': 'forged', 'parts': 12}} for n in range(4, 27)]},
        {'format': 'other', 'degrees': []},
    ]
    for doc in bad:
        try:
            with redirect_stdout(io.StringIO()):
                certificate.verify_joint(doc)
        except (ValueError, KeyError):
            pass
        else:
            raise AssertionError('invalid joint certificate accepted')
    bad = [
        {'format': certificate.BLOCKS_FORMAT, 'blocks': [{'nlo': 27, 'nhi': 2047, 'tree': {'unresolved': True}}]},
        {'format': certificate.BLOCKS_FORMAT, 'blocks': [{'nlo': 27, 'nhi': 2047, 'tree': {'test': 'direct', 'parts': 12}}]},
        {'format': certificate.BLOCKS_FORMAT, 'blocks': [{'nlo': 27, 'nhi': 100, 'tree': {'test': 'direct', 'parts': 12}}]},
        {'format': certificate.BLOCKS_FORMAT, 'blocks': [{'nlo': 28, 'nhi': 2047, 'tree': {'test': 'direct', 'parts': 12}}]},
    ]
    for doc in bad:
        try:
            with redirect_stdout(io.StringIO()):
                certificate.verify_blocks(doc)
        except (ValueError, KeyError):
            pass
        else:
            raise AssertionError('invalid block certificate accepted')


def test_obstruction():
    """The two-parameter estimates must not certify a box around the degree-6
    obstruction A = 6/5, e = sqrt(1/6) (where the direct/centered bounds exceed 1)."""
    from borcea import blocks
    assert blocks.check(6, 6, F(6, 5), F(6, 5), F(4082, 10000), F(4083, 10000), 64, 'uniform') is None
    assert blocks.check(6, 6, F(6, 5), F(6, 5), F(4082, 10000), F(4083, 10000), 48, 'geometric') is None


def run():
    test_arithmetic()
    test_quadrature()
    test_coefficient_identities()
    test_quartic_constants()
    test_restriction_soundness()
    test_fail_closed()
    test_obstruction()
    return 'PASS: arithmetic, quadrature, coefficient identities, quartic constants, restriction soundness, fail-closed and obstruction checks.'


if __name__ == '__main__':
    print(run())
