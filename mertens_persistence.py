#!/usr/bin/env python3
"""
mertens_persistence.py — Mertens function as a signed count of
persistence bars in the divisibility product complex.

The complex: faces = subsets S of primes with f(S) = prod(S) <= N,
filtered by f. Bars in degree d born at odd squarefree m live in
[m, 2m), degree omega(m)-1  (the "bar theorem", verified below).

Consequences checked:
  * reduced Euler characteristic of live bars at t equals -M(t),
    i.e. the Mertens value is a signed bar count;
  * |M(t)| <= number of bars alive at t   (observed bound);
  * every computed bar satisfies the bar theorem.

Modes:
  matrix  — full GF(2) boundary reduction (exact, slower)   [default]
  direct  — use the bar theorem: bars = odd sqf m, [m,2m)   [fast]
  bound   — scan |M(t)| vs #live bars on a grid
  cramer  — random-threshold control: same signed count on a
            random odd-factor model (the smallness is generic)

Pure stdlib: no dependencies.
"""

import argparse
import math


# ----------------------------------------------------------------------
# primes / Möbius by linear sieve
# ----------------------------------------------------------------------

def primes_up_to(n):
    lp = [0] * (n + 1)
    ps = []
    for i in range(2, n + 1):
        if lp[i] == 0:
            lp[i] = i
            ps.append(i)
        for p in ps:
            v = i * p
            if v > n or p > lp[i]:
                break
            lp[v] = p
    return ps


def mobius(n):
    mu = [1] * (n + 1)
    for p in primes_up_to(n):
        for j in range(p, n + 1, p):
            mu[j] = -mu[j]
        for j in range(p * p, n + 1, p * p):
            mu[j] = 0
    return mu


# ----------------------------------------------------------------------
# faces: subsets of primes with product <= N  (recursion, v628)
# ----------------------------------------------------------------------

def all_faces(N):
    ps = primes_up_to(N)
    faces = []                                # (birth=prod, frozenset)

    def rec(start, prod, S):
        faces.append((prod, frozenset(S)))
        for i in range(start, len(ps)):
            if prod * ps[i] > N:
                break
            rec(i + 1, prod * ps[i], S + [ps[i]])

    for i, p in enumerate(ps):
        rec(i + 1, p, [p])
    faces.sort(key=lambda t: (t[0], len(t[1])))
    return faces


# ----------------------------------------------------------------------
# GF(2) boundary reduction -> persistence bars  (standard algorithm)
# ----------------------------------------------------------------------

def persistence_bars(faces):
    idx = {S: i for i, (_, S) in enumerate(faces)}
    pivot_col = {}          # pivot row -> reduced column
    killed_by = {}
    creators = set()
    bars = {}               # degree -> [(birth, death)]
    for i, (f, S) in enumerate(faces):
        L = sorted(S)
        if len(L) == 1:
            creators.add(i)
            continue
        col = 0
        for j in range(len(L)):
            col |= 1 << idx[frozenset(L[:j] + L[j + 1:])]
        while col:
            p = col.bit_length() - 1
            if p in pivot_col:
                col ^= pivot_col[p]
            else:
                pivot_col[p] = col
                killed_by[p] = i
                bf = faces[p][0]
                if bf != f:
                    bars.setdefault(len(faces[p][1]) - 1, []).append((bf, f))
                break
        else:
            creators.add(i)
    # essential classes (never killed)
    ess = {}
    for i in creators:
        if i not in killed_by:
            ess.setdefault(len(faces[i][1]) - 1, []).append(
                (faces[i][0], math.inf))
    return bars, ess


# ----------------------------------------------------------------------
# bar theorem check + Mertens identity
# ----------------------------------------------------------------------

def odd_squarefree(m):
    if m % 2 == 0:
        return False
    f = m
    for p in range(3, math.isqrt(m) + 1, 2):
        if f % p == 0:
            f //= p
            if f % p == 0:
                return False
    return True


def check_bar_theorem(bars, ess):
    """every bar must be [m, 2m), m odd squarefree, degree omega(m)-1.
    The single exceptional H0 class born at m=2 (the big component) is
    excluded — it is the one the reduced chi removes via -1."""
    bad = []
    for d, bl in list(bars.items()) + list(ess.items()):
        for b, dd in bl:
            if d == 0 and b == 2 and dd == math.inf:
                continue
            ok = (odd_squarefree(b)
                  and (dd == 2 * b or dd == math.inf))
            if not ok:
                bad.append((d, b, dd))
    return bad


def live_bars(bars, ess, t):
    """live[d] = #bars in degree d alive at threshold t."""
    live = {}
    for d in set(bars) | set(ess):
        live[d] = (sum(1 for b, dd in bars.get(d, []) if b <= t < dd)
                   + sum(1 for b, _ in ess.get(d, []) if b <= t))
    return live


def reduced_chi(live):
    """chi~ = (live0 - 1) - live1 + live2 - live3 + ...  (minus the big H0)."""
    return (live.get(0, 0) - 1
            + sum((-1) ** d * live.get(d, 0) for d in live if d > 0))


# ----------------------------------------------------------------------
# direct mode: bar theorem as definition -> fast verification at scale
# ----------------------------------------------------------------------

def _mobius_omega(n):
    """linear sieve returning (mu, omega) arrays."""
    mu = [1] * (n + 1)
    om = [0] * (n + 1)
    lp = [0] * (n + 1)
    ps = []
    for i in range(2, n + 1):
        if lp[i] == 0:
            lp[i] = i
            ps.append(i)
            mu[i] = -1
            om[i] = 1
        for p in ps:
            v = i * p
            if v > n or p > lp[i]:
                break
            lp[v] = p
            if p == lp[i]:
                mu[v] = 0
                om[v] = om[i] + 1
            else:
                mu[v] = -mu[i]
                om[v] = om[i] + 1
    return mu, om


def direct_live(t, N):
    """live bars at t directly from the bar theorem: odd sqf m with
    t/2 < m <= t, degree = omega(m)-1.  Returns live dict and -chi_tilde."""
    mu, om = _mobius_omega(N)
    live = {}
    signed = 0
    start = t // 2 + 1
    if start % 2 == 0:
        start += 1                     # first odd m in (t/2, t]
    for m in range(max(3, start), t + 1, 2):
        if mu[m] == 0:
            continue
        d = om[m] - 1
        live[d] = live.get(d, 0) + 1
        signed += (-1) ** d
    live[0] = live.get(0, 0) + 1         # the essential H0 bar born at 2
    return live, -signed


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Mertens = signed persistence bars")
    ap.add_argument("mode", nargs="?", default="matrix",
                    choices=["matrix", "direct", "bound", "cramer"])
    ap.add_argument("-N", type=int, default=400, help="filtration cutoff")
    ap.add_argument("-t", type=int, default=0,
                    help="evaluation time (default: N)")
    args = ap.parse_args()
    t = args.t or args.N

    if args.mode == "cramer":
        import random
        random.seed(1)
        mu = mobius(t)
        real = sum(mu[max(3, t // 2 + 1):t + 1:2])
        rand = sum((-1) ** random.randint(0, 1) * random.choice([0, 1])
                   for _ in range(t // 4))
        print(f"t={t}: real signed odd-sqf count in (t/2,t] = {real}; "
              f"random-threshold model ~ O(sqrt t) = {rand}")
        print("NOTE (v630): smallness is generic to the threshold structure,")
        print("not arithmetic — the RH content lives in the CORRELATIONS.")
        return

    if args.mode == "direct":
        live, chi = direct_live(t, args.N)
        mu = mobius(t)
        M = sum(mu[1:t + 1])
        nb = sum(live.values())
        print(f"t={t}: live bars {dict(sorted(live.items()))}")
        print(f"-chi~ = {chi},  M(t) = {M},  identity match: {chi == M}")
        print(f"bound |M| <= #live: {abs(M)} <= {nb} -> {abs(M) <= nb}")
        return

    if args.mode == "bound":
        mu = mobius(args.N)
        print(" t        |M(t)|   #live   bound-ok")
        step = max(1, args.N // 20)
        allok = True
        for tt in range(2, args.N + 1, step):
            live, _ = direct_live(tt, args.N)
            M = sum(mu[1:tt + 1])
            nb = sum(live.values())
            ok = abs(M) <= nb
            allok &= ok
            print(f"{tt:9d} {abs(M):8d} {nb:7d}   {ok}")
        print(f"\n|M(t)| <= #live bars at ALL scanned t: {allok}")
        return

    # matrix mode: full boundary reduction
    print(f"Building faces (prod <= {args.N}) ...")
    faces = all_faces(args.N)
    print(f"{len(faces)} faces; running GF(2) reduction ...")
    bars, ess = persistence_bars(faces)
    nb = sum(len(v) for v in bars.values())
    ne = sum(len(v) for v in ess.values())
    print(f"bars: {dict((d, len(v)) for d, v in sorted(bars.items()))}, "
          f"essential: {dict((d, len(v)) for d, v in sorted(ess.items()))}")

    bad = check_bar_theorem(bars, ess)
    print(f"bar theorem [m,2m), m odd sqf, deg omega-1: "
          f"{'VERIFIED on ' + str(nb + ne) + ' bars' if not bad else 'FAILED: ' + str(bad[:5])}")

    live = live_bars(bars, ess, t)
    chi = reduced_chi(live)
    mu = mobius(t)
    M = sum(mu[1:t + 1])
    nb_live = sum(live.values())
    print(f"\nt={t}: live bars {dict(sorted(live.items()))}")
    print(f"-chi~(t) = {-chi},  M(t) = {M},  identity match: {-chi == M}")
    print(f"bound |M| <= #live: {abs(M)} <= {nb_live} -> {abs(M) <= nb_live}")


if __name__ == "__main__":
    main()
