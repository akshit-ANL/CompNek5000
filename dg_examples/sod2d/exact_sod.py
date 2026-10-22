#!/usr/bin/env python3
"""Exact Sod solution (Toro, ch. 4) and comparison with prof_*.dat.

    python exact_sod.py            # merge prof_*.dat, print L1/Linf errors
    python exact_sod.py --table    # print the exact solution on a grid

Standard library only. Columns of prof files: x y rho u p.
"""
import glob
import math
import sys

GAMMA = 1.4
LEFT = (1.0, 0.0, 1.0)      # rho, u, p
RIGHT = (0.125, 0.0, 0.1)
X0 = 0.5


def _f(p, rho, pk):
    g = GAMMA
    a = math.sqrt(g * pk / rho)
    if p > pk:  # shock
        A = 2.0 / ((g + 1.0) * rho)
        B = (g - 1.0) / (g + 1.0) * pk
        return (p - pk) * math.sqrt(A / (p + B))
    return 2.0 * a / (g - 1.0) * ((p / pk) ** ((g - 1.0) / (2.0 * g)) - 1.0)


def _df(p, rho, pk):
    g = GAMMA
    a = math.sqrt(g * pk / rho)
    if p > pk:
        A = 2.0 / ((g + 1.0) * rho)
        B = (g - 1.0) / (g + 1.0) * pk
        return math.sqrt(A / (B + p)) * (1.0 - (p - pk) / (2.0 * (B + p)))
    return 1.0 / (rho * a) * (p / pk) ** (-(g + 1.0) / (2.0 * g))


def star_state():
    (rl, ul, pl), (rr, ur, pr) = LEFT, RIGHT
    p = 0.5 * (pl + pr)
    for _ in range(100):
        f = _f(p, rl, pl) + _f(p, rr, pr) + (ur - ul)
        df = _df(p, rl, pl) + _df(p, rr, pr)
        pn = max(p - f / df, 1e-12)
        if abs(pn - p) < 1e-14 * (pn + p):
            p = pn
            break
        p = pn
    u = 0.5 * (ul + ur) + 0.5 * (_f(p, rr, pr) - _f(p, rl, pl))
    return p, u


def exact(x, t):
    """Return (rho, u, p) at position x and time t (left rarefaction,
    right shock configuration, as in the Sod problem)."""
    g = GAMMA
    (rl, ul, pl), (rr, ur, pr) = LEFT, RIGHT
    ps, us = star_state()
    s = (x - X0) / t
    al = math.sqrt(g * pl / rl)
    rsl = rl * (ps / pl) ** (1.0 / g)
    asl = al * (ps / pl) ** ((g - 1.0) / (2.0 * g))
    rsr = rr * ((ps / pr + (g - 1.0) / (g + 1.0)) /
                ((g - 1.0) / (g + 1.0) * ps / pr + 1.0))
    ss = ur + 0.0 + math.sqrt(g * pr / rr) * math.sqrt(
        (g + 1.0) / (2.0 * g) * ps / pr + (g - 1.0) / (2.0 * g))
    if s < ul - al:
        return rl, ul, pl
    if s < us - asl:  # inside the rarefaction fan
        c = 2.0 / (g + 1.0) * (al + (g - 1.0) / 2.0 * (ul - s))
        rho = rl * (c / al) ** (2.0 / (g - 1.0))
        return rho, 2.0 / (g + 1.0) * (al + (g - 1.0) / 2.0 * ul + s), \
            pl * (c / al) ** (2.0 * g / (g - 1.0))
    if s < us:
        return rsl, us, ps
    if s < ss:
        return rsr, us, ps
    return rr, ur, pr


def main():
    t = 0.2
    if "--table" in sys.argv:
        ps, us = star_state()
        print("p* = %.5f  u* = %.5f" % (ps, us))
        for i in range(0, 21):
            x = i / 20.0
            print("%.2f  %.5f %.5f %.5f" % ((x,) + exact(x, t)))
        return
    files = sorted(glob.glob("prof_*.dat"))
    if not files:
        sys.exit("no prof_*.dat files found")
    n, l1, linf, l1u = 0, 0.0, 0.0, 0.0
    for fn in files:
        for line in open(fn):
            x, y, rho, u, p = map(float, line.split())
            re, ue, pe = exact(x, t)
            l1 += abs(rho - re)
            l1u += abs(p - pe)
            linf = max(linf, abs(rho - re))
            n += 1
    print("points: %d" % n)
    print("density  L1 = %.4e   Linf = %.4e" % (l1 / n, linf))
    print("pressure L1 = %.4e" % (l1u / n))


if __name__ == "__main__":
    main()
