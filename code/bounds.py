import gmpy2
from gmpy2 import mpz
from mpmath import iv
from mpmath.libmp import to_rational

P = 320
IV_PREC = 400


def bounds_exact(q, Cmax):
    q = mpz(q)
    J = [0] * (Cmax + 1)
    Bu = [0] * (Cmax + 1)
    Q = mpz(1)
    for c in range(1, Cmax + 1):
        Q *= q
        j = Q.bit_length()
        N = (mpz(1) << j) - Q
        assert N > 0
        num = c * (2 * Q + N)
        den = 2 * q * N
        J[c] = j
        Bu[c] = int(-((-num) // den))
    return J, Bu


def endpoints(x):
    a, b = x._mpi_
    return to_rational(a), to_rational(b)


def floor_div(p, d):
    return int(p) // int(d)


def ceil_div(p, d):
    return -((-int(p)) // int(d))


def log_enclosures(q):
    iv.prec = IV_PREC
    L = iv.log(iv.mpf(q)) / iv.log(iv.mpf(2))
    l2 = iv.log(iv.mpf(2))
    (pa, da), (pb, db) = endpoints(L)
    (ra, ea), _ = endpoints(l2)
    Llo = floor_div(pa << P, da)
    Lhi = ceil_div(pb << P, db)
    l2lo = floor_div(ra << P, ea)
    assert Llo < Lhi and Lhi - Llo <= 4, "enclosure too wide"
    return Llo, Lhi, l2lo


def bounds_interval(q, Cmax):
    Llo, Lhi, l2lo = log_enclosures(q)
    J = [0] * (Cmax + 1)
    Bu = [0] * (Cmax + 1)
    one = 1 << P
    for c in range(1, Cmax + 1):
        flo = (c * Llo) >> P
        fhi = (c * Lhi) >> P
        if flo != fhi:
            raise ArithmeticError(f"q={q}, c={c}: enclosure cannot decide ceil(c log2 q)")
        j = flo + 1
        dlo = j * one - c * Lhi
        assert dlo > 0
        num = c << (2 * P)
        den = q * dlo * l2lo
        J[c] = j
        Bu[c] = -((-num) // den)
    return J, Bu


def certified_bounds(q, Cmax):
    JA, BA = bounds_exact(q, Cmax)
    JB, BB = bounds_interval(q, Cmax)
    for c in range(1, Cmax + 1):
        if JA[c] != JB[c]:
            raise AssertionError(f"q={q} c={c}: j mismatch {JA[c]} vs {JB[c]}")
        if BB[c] > BA[c] + 1 or 12 * BA[c] > 13 * BB[c] + 12:
            raise AssertionError(f"q={q} c={c}: bound mismatch {BA[c]} vs {BB[c]}")
    Bup = [max(a, b) for a, b in zip(BA, BB)]
    worst = 0
    for c, _ in records(Bup):
        d = abs(BA[c] - BB[c])
        if d > 1:
            raise AssertionError(f"q={q} record c={c}: {BA[c]} vs {BB[c]}")
        worst = max(worst, d)
    return JA, Bup, worst


def records(Bup):
    out, best = [], 0
    for c in range(1, len(Bup)):
        if Bup[c] > best:
            best = Bup[c]
            out.append((c, Bup[c]))
    return out


def theorem_table(Bup, recs):
    rows = []
    for (c, b), (cn, _) in zip(recs, recs[1:]):
        rows.append({"C": cn - 1, "B": b, "record_c": c})
    return rows


def staircase(Bup, C):
    out, best = [], 0
    for c in range(C, 0, -1):
        if Bup[c] > best:
            best = Bup[c]
            out.append((c, Bup[c]))
    return out[::-1]


def continued_fraction_log2(q, nterms=60):
    iv.prec = 4000
    x = iv.log(iv.mpf(q)) / iv.log(iv.mpf(2))
    a = []
    for _ in range(nterms):
        (pa, da), (pb, db) = endpoints(x)
        fa, fb = floor_div(pa, da), floor_div(pb, db)
        if fa != fb:
            break
        a.append(fa)
        x = 1 / (x - fa)
    return a


def convergents(a):
    p0, p1, s0, s1 = 1, a[0], 0, 1
    out = [(p1, s1)]
    for t in a[1:]:
        p0, p1 = p1, t * p1 + p0
        s0, s1 = s1, t * s1 + s0
        out.append((p1, s1))
    return out


def check_convergents_exact(q, conv, smax):
    q = mpz(q)
    ok = []
    for i, (p, s) in enumerate(conv):
        if s > smax:
            break
        below = (mpz(1) << p) < q ** s
        if below != (i % 2 == 0):
            raise AssertionError(f"q={q}: convergent {p}/{s} on wrong side")
        ok.append(i)
    return len(ok)
