"""Outward-rounded interval arithmetic on integers (Section 6.2).

An interval is a pair of integers [lo, hi] representing [lo/S, hi/S] with
S = 2**96.  Every operation rounds outwards, so the computed interval contains
the exact value at every point of its inputs.  No floating-point number is
ever used to accept a box.
"""
from fractions import Fraction
from math import isqrt

BITS = 96
SCALE = 1 << BITS


def _ceildiv(a, b):
    return -((-a) // b)


class I:
    """Real interval [lo/SCALE, hi/SCALE] with integer endpoints."""
    __slots__ = ('lo', 'hi')

    def __init__(self, value=0, hi=None, _raw=False):
        if _raw:
            self.lo, self.hi = value, hi
        elif isinstance(value, I):
            self.lo, self.hi = value.lo, value.hi
        elif isinstance(value, int) and (hi is None or isinstance(hi, int)):
            self.lo = value * SCALE
            self.hi = (value if hi is None else hi) * SCALE
        else:
            x = Fraction(value)
            y = x if hi is None else Fraction(hi)
            self.lo = x.numerator * SCALE // x.denominator
            self.hi = _ceildiv(y.numerator * SCALE, y.denominator)
        if self.lo > self.hi:
            raise ValueError('reversed interval')

    @staticmethod
    def raw(lo, hi):
        return I(lo, hi, _raw=True)

    def __add__(self, other):
        other = other if isinstance(other, I) else I(other)
        return I.raw(self.lo + other.lo, self.hi + other.hi)
    __radd__ = __add__

    def __neg__(self):
        return I.raw(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -I(other)

    def __rsub__(self, other):
        return I(other) + -self

    def __mul__(self, other):
        other = other if isinstance(other, I) else I(other)
        if self.lo >= 0 and other.lo >= 0:
            return I.raw(self.lo * other.lo // SCALE, _ceildiv(self.hi * other.hi, SCALE))
        p = (self.lo * other.lo, self.lo * other.hi, self.hi * other.lo, self.hi * other.hi)
        return I.raw(min(p) // SCALE, _ceildiv(max(p), SCALE))
    __rmul__ = __mul__

    def __truediv__(self, other):
        other = other if isinstance(other, I) else I(other)
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError('interval contains zero')
        if other.hi < 0:
            return (-self) / (-other)
        if self.lo >= 0:
            return I.raw(self.lo * SCALE // other.hi, _ceildiv(self.hi * SCALE, other.lo))
        pairs = ((self.lo, other.lo), (self.lo, other.hi), (self.hi, other.lo), (self.hi, other.hi))
        return I.raw(min(x * SCALE // y for x, y in pairs),
                     max(_ceildiv(x * SCALE, y) for x, y in pairs))

    def __rtruediv__(self, other):
        return I(other) / self

    def __pow__(self, n):
        if not isinstance(n, int) or n < 0:
            raise ValueError('nonnegative integer powers only')
        if self.lo < 0:
            raise ValueError('power base must be nonnegative')
        out, base = I(1), self
        while n:
            if n & 1:
                out = out * base
            n //= 2
            if n:
                base = base * base
        return out

    def sqrt(self):
        if self.lo < 0:
            raise ValueError('negative square root')
        lo = isqrt(self.lo * SCALE)
        hi = isqrt(self.hi * SCALE)
        return I.raw(lo, hi + (hi * hi < self.hi * SCALE))

    def halfpower(self, numerator):
        """x**(numerator/2) for x >= 0."""
        out = self ** (numerator // 2)
        return out * self.sqrt() if numerator % 2 else out

    def positive(self):
        """Intersection with [0, infinity)."""
        return I.raw(max(0, self.lo), max(0, self.hi))

    def upper(self):
        return I.raw(self.hi, self.hi)

    def lower(self):
        return I.raw(self.lo, self.lo)


def imin(*xs):
    xs = [I(x) for x in xs]
    return I.raw(min(x.lo for x in xs), min(x.hi for x in xs))


def imax(*xs):
    xs = [I(x) for x in xs]
    return I.raw(max(x.lo for x in xs), max(x.hi for x in xs))


def square(x):
    """Enclosure of x**2 for a real interval that may contain 0."""
    if x.lo >= 0:
        return x * x
    if x.hi <= 0:
        return (-x) * (-x)
    return imax(x.lower() * x.lower(), x.upper() * x.upper()).upper() * I(0, 1)


class C:
    """Complex interval: a rectangle re + i im."""

    def __init__(self, re=0, im=0):
        if isinstance(re, C):
            self.re, self.im = re.re, re.im
        elif isinstance(re, complex):
            self.re, self.im = I(re.real), I(re.imag)
        else:
            self.re, self.im = I(re), I(im)

    def __add__(self, o):
        o = C(o)
        return C(self.re + o.re, self.im + o.im)
    __radd__ = __add__

    def __neg__(self):
        return C(-self.re, -self.im)

    def __sub__(self, o):
        return self + -C(o)

    def __rsub__(self, o):
        return C(o) + -self

    def __mul__(self, o):
        o = C(o)
        return C(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    __rmul__ = __mul__

    def norm(self):
        return (square(self.re) + square(self.im)).positive().sqrt()


def lt_one(x):
    """True when the interval lies strictly below 1."""
    return x.hi < SCALE
