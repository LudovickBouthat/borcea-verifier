"""Positive-weight chord rule (Lemma 6.2) and the two rational meshes (Section 6.3).

For consecutive mesh points l < h and j in {1, 2},

    int_l^h t^j g(t) dt <= W_l g(l) + W_h g(h)

for every nonnegative convex g, with the nonnegative rational weights

    W_l = int_l^h t^j (h - t)/(h - l) dt,   W_h = int_l^h t^j (t - l)/(h - l) dt.

A mesh is returned as a tuple of (l, h, ((W_l, W_h, mass) for j = 1, 2)),
all as intervals enclosing the exact rationals.
"""
from fractions import Fraction as F
from functools import lru_cache
from .interval import I


def _weights(points):
    out = []
    for l, h in zip(points, points[1:]):
        weights = []
        for j in (1, 2):
            mass = (h ** (j + 1) - l ** (j + 1)) / (j + 1)
            moment = (h ** (j + 2) - l ** (j + 2)) / (j + 2)
            weights.append((I((h * mass - moment) / (h - l)),
                            I((moment - l * mass) / (h - l)), I(mass)))
        out.append((I(l), I(h), tuple(weights)))
    return tuple(out)


@lru_cache(None)
def uniform_mesh(q):
    """Union of the points j/q and j/(16q), 0 <= j <= q."""
    points = sorted({F(j, q) for j in range(q + 1)} | {F(j, 16 * q) for j in range(q + 1)})
    return _weights(points)


@lru_cache(None)
def geometric_mesh(q):
    """0, 2^-15, ..., 2^-1, 1, each gap cut into max(1, ceil(q/12)) equal parts,
    together with the reflection of every point about 1/2."""
    sub = max(1, (q + 11) // 12)
    breaks = [F(0)] + [F(1, 2 ** j) for j in range(15, 0, -1)] + [F(1)]
    points = {F(0), F(1)}
    for l, h in zip(breaks, breaks[1:]):
        for j in range(sub + 1):
            t = l + (h - l) * F(j, sub)
            points.update((t, 1 - t))
    return _weights(sorted(points))


MESHES = {'uniform': uniform_mesh, 'geometric': geometric_mesh}
