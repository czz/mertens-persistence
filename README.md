# mertens-persistence

The Mertens function M(t) = Σ_{n≤t} μ(n) as a **signed count of
persistence bars** in the divisibility product complex — a new
combinatorial-topological object connecting Mertens to TDA.

Single file, pure Python stdlib (no dependencies).

## The object

Faces = subsets S of primes with product f(S) = Π S ≤ N, filtered by f.
Persistent homology (GF(2) boundary reduction) gives bars with an exact
structure — the **bar theorem**:

    every bar is [m, 2m) with m odd squarefree, in degree ω(m)−1.

Then the reduced Euler characteristic of bars alive at threshold t
satisfies the identity

    −χ̃(t) = M(t)

i.e. the Mertens value is a signed count of live persistence bars —
and the observed bound **|M(t)| ≤ #bars alive at t** holds at every
scanned t (the "many bars, small signed sum" window into RH-scale
cancellation).

## Modes

    python3 mertens_persistence.py matrix -N 400    # full GF(2) reduction:
                                                    # bars + bar-theorem check
    python3 mertens_persistence.py direct -N 100000 # bar theorem as formula,
                                                    # verifies -chi~ = M(t)
    python3 mertens_persistence.py bound -N 100000  # |M(t)| <= #live scan
    python3 mertens_persistence.py cramer -N 10000  # random-threshold control

## Verified

- bar theorem on all computed bars (matrix mode, N=400: 162/162)
- identity −χ̃(t) = M(t) exact at t = 100 … 50000 (matrix + direct)
- bound |M(t)| ≤ #live bars at every scanned t ≤ 100000
- honest negative note: the signed-count smallness is generic to the
  threshold structure (random model gives ~√t too) — the arithmetic
  content lives in the correlations, not the topology alone.

## Context

Extracted from the RH research project (scripts v608–v635, documented in
RH_final_consolidated_v472.md §"Mertens = conteggio firmato di barre di
persistenza"; the order-complex χ̃ = M−1 side is Philip Hall's theorem).
