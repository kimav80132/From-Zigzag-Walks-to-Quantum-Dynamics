"""The Pauli test of the particle definition (III.15a: an event = emitter + absorber, weight = (sum of out-paths) x (sum of return paths)).
Two identical walkers, emitters A1, A2, absorbers X1, X2. Both pairings count. Does a relative sign come out of the walk by itself?
Lattice: 1+1 t-reading, column over (site, sense), one tick = shift (+ right, - left) then the coin R = cos(a dt) 1 + sin(a dt) T. U = R S over N ticks, orthogonal.
1. single walker, III.15a: closure sum_X U_{XA}^2 = 1 and the loop weight (U^T U)_{AA} = +1: a closed out-and-back loop of one walker weighs +1, not -1.
2. two independent walkers: the joint step is U (x) U, which commutes with the exchange of the two walkers: the step puts NO sign on the exchanged pairing.
   With III.15a's rule and both pairings, the event weight is (D + E)^2, D = U_{X1A1}U_{X2A2}, E = U_{X2A1}U_{X1A2}; at X1 = X2 it is 4D^2, not 0.
   Exclusion would need (D - E)^2, i.e. a -1 on the exchanged pairing = a -1 per closed loop (DE is the one-loop configuration through all four points).
3. the x-clock reading: one path visits each clock value x exactly once (III.1: the walker always advances in x). So two crossings of ONE path at a fixed t are never
   at the same x: the several particles the forward-time observer sees from one path (III.1) never coincide. Exclusion for crossings of one path, kinematic.
4. the fermion sign used in the Casimir document (-1/2 Tr|H|: a -1 per closed loop) is an input there; the walk's own loop closes with +1 (item 1)."""
import numpy as np, itertools
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
L, N, a, dt = 24, 8, 1.0, 0.25
T2 = np.array([[0, 1], [-1, 0.]]); R = np.cos(a*dt)*np.eye(2) + np.sin(a*dt)*T2
dim = 2*L
idx = lambda x, e: 2*(x % L) + (0 if e == +1 else 1)
S = np.zeros((dim, dim))                       # shift: + moves right, - moves left (periodic)
for x in range(L):
    S[idx(x+1, +1), idx(x, +1)] = 1; S[idx(x-1, -1), idx(x, -1)] = 1
C = np.kron(np.eye(L), R)                      # the coin at every site
step = C @ S
U = np.linalg.matrix_power(step, N)
ok("0. the step is orthogonal: U^T U = 1", np.allclose(U.T @ U, np.eye(dim)))
A = idx(0, +1)
ok("1. one walker: sum_X U_{XA}^2 = 1 (III.15a (C)4) and the closed loop (U^T U)_{AA} = +1: a loop of one walker weighs plus one",
   abs(np.sum(U[:, A]**2) - 1) < 1e-12 and abs((U.T @ U)[A, A] - 1) < 1e-12)
# 2. two independent walkers
A1, A2 = idx(0, +1), idx(6, +1)
P = np.zeros((dim*dim, dim*dim))               # exchange of the two walkers on the joint column
for i in range(dim):
    for j in range(dim):
        P[j*dim + i, i*dim + j] = 1
UU = np.kron(U, U)
ok("2a. the joint step U (x) U commutes with the exchange P of the two walkers: the step is blind to which walker is which, no sign anywhere",
   np.allclose(P @ UU, UU @ P))
D = lambda X1, X2: U[X1, A1]*U[X2, A2]
E = lambda X1, X2: U[X2, A1]*U[X1, A2]
# event weight by III.15a's rule with both pairings in the out-sum and both in the return-sum (return = transpose = same numbers)
W_rule = lambda X1, X2: (D(X1, X2) + E(X1, X2))**2
W_excl = lambda X1, X2: (D(X1, X2) - E(X1, X2))**2
same = [(X, X) for X in range(dim) if abs(D(X, X)) > 1e-6]
ok(f"2b. same endpoint X1 = X2: the rule gives (D+E)^2 = 4D^2 > 0 at {len(same)} sites (no exclusion); exclusion would need (D-E)^2 = 0 there",
   all(abs(W_rule(X, X) - 4*D(X, X)**2) < 1e-12 and W_rule(X, X) > 0 and W_excl(X, X) < 1e-12 for X, _ in same))
tot_rule = sum(W_rule(X1, X2) for X1 in range(dim) for X2 in range(dim))
tot_excl = sum(W_excl(X1, X2) for X1 in range(dim) for X2 in range(dim))
ok(f"2c. summed over all endpoint pairs: rule {tot_rule:.6f}, exclusive {tot_excl:.6f}; both = 2(1 + or - overlap^2) with overlap = sum_X U_XA1 U_XA2 = {np.dot(U[:, A1], U[:, A2]):.4f}: "
   "the cross term is the one-loop configuration through all four points, and the walk gives it +1. (Over ordered pairs; with the 1/2 of identical walkers each sums to 1.)",
   abs(tot_rule - 2*(1 + np.dot(U[:, A1], U[:, A2])**2)) < 1e-9 and abs(tot_excl - 2*(1 - np.dot(U[:, A1], U[:, A2])**2)) < 1e-9)
# 3. the x-clock: one path, one position per clock value
rng = np.random.default_rng(1)
def xclock_path(n):            # III.1: advance in x every tick, t moves by +-1, sense flips with prob a*dx/c
    t, e, pts = 0, +1, []
    for x in range(n):
        pts.append((x, t)); t += e
        if rng.random() < a*dt: e = -e
    return pts
paths = [xclock_path(60) for _ in range(2000)]
ok("3. x-clock reading: every path visits each x exactly once, so two crossings of one path at a fixed t are at different x (2000 paths, 60 steps): "
   "the several particles the forward-time observer sees from ONE path never coincide — exclusion for one path, kinematic",
   all(len(set(x for x, t in p)) == len(p) for p in paths)
   and all(len(set(x for x, t in p if t == t0)) == sum(1 for x, t in p if t == t0) for p in paths[:200] for t0 in set(t for x, t in p)))
# 4. the noise belongs to the event (II.1, III.15b (A)): in the walk's clock, two walkers at the same (t, sense) at clock x draw the SAME coin there, and forever after
rng = np.random.default_rng(7)
def two_walkers_shared_coin(n, start1, start2):
    (t1, e1), (t2, e2) = start1, start2; same = []
    for x in range(n):
        coins = {}                                           # one coin per event (x, t): the event's, not the walker's
        for key in {(x, t1), (x, t2)}: coins[key] = rng.random() < a*dt
        f1, f2 = coins[(x, t1)], coins[(x, t2)]
        t1 += e1; t2 += e2
        if f1: e1 = -e1
        if f2: e2 = -e2
        same.append((t1, e1) == (t2, e2))
    return same
runs_same = [two_walkers_shared_coin(40, (0, +1), (0, +1)) for _ in range(300)]
runs_opp = [two_walkers_shared_coin(40, (0, +1), (0, -1)) for _ in range(300)]
ok("4. event-attached coins: two walkers starting in the same state (same t, same sense) stay identical at every step in 300/300 runs - one process, twice the weight; "
   "two starting at the same event with opposite senses separate at once in 300/300 runs - two processes. One future per point and direction",
   all(all(r) for r in runs_same) and all(not r[0] for r in runs_opp))
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
