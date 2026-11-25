"""Can the companion rule be purely kinematic? Two identical walkers; free step (shift then coin); if the next shift would put them on the same point they bounce:
both senses reversed, positions unchanged for that tick, no coin charged (elastic collision of two light-like particles = exchange of velocities).
Same lattice parity, so the distance is always even and coincidence is the only way to cross.
1. the bounce shift S_c is a bijection on non-coincident ordered pairs, so the joint step (C x C) S_c is orthogonal: the weight is conserved (the path rule leaked; this does not).
2. no two walkers ever share a point: exclusion, kinematic.
3. the distribution over absorber pairs P(x1,e1; x2,e2), compared with the free fermion (antisymmetric column, U x U) and the free boson (symmetric column):
   does the elastic bounce reproduce the fermion's |D - E|^2 (Girardeau: hard-core = fermion in position), the boson's, or neither?
(Check 3 rebuilt 10 Nov 2026.  The earlier version printed the largest single differences and passed "verdict recorded above" without a test; the numbers
 the text quotes, total variations and the head-on weight, were not computed.  The boson's distribution is now taken over ALL its outcomes, including
 two walkers in one state, which carry 0.113 of its weight; without them it is not a distribution and its distance comes out too small.)"""
import numpy as np
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
L, N, a, dt = 16, 7, 1.0, 0.3
c_, s_ = np.cos(a*dt), np.sin(a*dt)
T2 = np.array([[0, 1], [-1, 0.]]); R = c_*np.eye(2) + s_*T2
dim = 2*L; idx = lambda x, e: 2*(x % L) + (0 if e == +1 else 1); site = lambda i: (i//2, +1 if i % 2 == 0 else -1)
S = np.zeros((dim, dim))
for x in range(L): S[idx(x+1, +1), idx(x, +1)] = 1; S[idx(x-1, -1), idx(x, -1)] = 1
C = np.kron(np.eye(L), R); U = np.linalg.matrix_power(C @ S, N)
# joint space: ordered pairs (i, j) of single-walker states, i != j and different sites
pairs = [(i, j) for i in range(dim) for j in range(dim) if site(i)[0] != site(j)[0]]
pidx = {p: k for k, p in enumerate(pairs)}; J = len(pairs)
# free shift of a single state
def shift(i):
    x, e = site(i); return idx(x + e, e)
# bounce shift on pairs
Sc = np.zeros((J, J))
for (i, j), k in pidx.items():
    i2, j2 = shift(i), shift(j)
    if site(i2)[0] == site(j2)[0]:                        # would coincide: bounce, senses reversed, positions kept
        xi, ei = site(i); xj, ej = site(j)
        i2, j2 = idx(xi, -ei), idx(xj, -ej)
    Sc[pidx[(i2, j2)], k] = 1
ok("1. the bounce shift is a permutation of the non-coincident ordered pairs (each column and each row has one 1)",
   np.allclose(Sc.sum(0), 1) and np.allclose(Sc.sum(1), 1))
# joint coin on pairs (site-local, so it maps pairs to pairs with the same sites)
Cj = np.zeros((J, J))
for (i, j), k in pidx.items():
    xi, ei = site(i); xj, ej = site(j)
    for ei2 in (+1, -1):
        for ej2 in (+1, -1):
            Cj[pidx[(idx(xi, ei2), idx(xj, ej2))], k] = R[0 if ei2 == +1 else 1, 0 if ei == +1 else 1]*R[0 if ej2 == +1 else 1, 0 if ej == +1 else 1]
step_hc = Cj @ Sc
ok("1'. the joint step with the bounce is orthogonal: the two-walker weight is conserved", np.allclose(step_hc.T @ step_hc, np.eye(J)))
A1, A2 = idx(0, +1), idx(4, +1)
v = np.zeros(J); v[pidx[(A1, A2)]] = 1
w = np.linalg.matrix_power(step_hc, N) @ v
ok(f"2. after {N} ticks the weight is on non-coincident pairs only by construction, total {np.sum(w**2):.12f}: no two walkers share a point, nothing lost",
   abs(np.sum(w**2) - 1) < 1e-12)
# unordered distribution over absorber pairs for the three rules
def unordered(prob_pairs):
    acc = {}
    for (i, j), p in prob_pairs.items():
        key = tuple(sorted((i, j))); acc[key] = acc.get(key, 0) + p
    return acc
P_hc = unordered({p: w[k]**2 for p, k in pidx.items()})
D = lambda i, j: U[i, A1]*U[j, A2]; E = lambda i, j: U[j, A1]*U[i, A2]
P_F = unordered({(i, j): 0.5*(D(i, j) - E(i, j))**2 for i in range(dim) for j in range(dim)})          # i = j carries 0
P_B = unordered({(i, j): 0.5*(D(i, j) + E(i, j))**2 for i in range(dim) for j in range(dim)})          # i = j: two walkers in one state
tv = lambda P_, Q_: 0.5*sum(abs(P_.get(k, 0) - Q_.get(k, 0)) for k in set(P_) | set(Q_))
def positions(P):                                         # positions only: summed over the senses
    acc = {}
    for (i, j), p in P.items():
        key = tuple(sorted((site(i)[0], site(j)[0]))); acc[key] = acc.get(key, 0) + p
    return acc
X_hc, X_F, X_B = positions(P_hc), positions(P_F), positions(P_B)
same_site = lambda P: sum(p for (i, j), p in P.items() if site(i)[0] == site(j)[0] and i != j)        # one point, opposite senses: the head-on moment
same_state = lambda P: sum(p for (i, j), p in P.items() if i == j)
print(f"   total variation, absorber pairs with senses:  bounce-fermion {tv(P_hc, P_F):.3f}   bounce-boson {tv(P_hc, P_B):.3f}   fermion-boson {tv(P_F, P_B):.3f}")
print(f"   total variation, positions alone:             bounce-fermion {tv(X_hc, X_F):.3f}   bounce-boson {tv(X_hc, X_B):.3f}   fermion-boson {tv(X_F, X_B):.3f}")
print(f"   weight on one point, opposite senses:  fermion {same_site(P_F):.4f}   boson {same_site(P_B):.4f}   bounce {same_site(P_hc):.4f}")
print(f"   weight on one state:                   fermion {same_state(P_F):.4f}   boson {same_state(P_B):.4f}   bounce {same_state(P_hc):.4f}")
ok("3a. the three are distributions: each sums to 1", all(abs(sum(P.values()) - 1) < 1e-12 for P in (P_hc, P_F, P_B)))
ok(f"3b. the bounce is neither: total variation {tv(P_hc, P_F):.2f} from the fermion's (D-E)^2 and {tv(P_hc, P_B):.2f} from the boson's (D+E)^2; in positions alone "
   f"{tv(X_hc, X_F):.2f} and {tv(X_hc, X_B):.2f}. Hard-core walkers are a third statistics here, also in position (no Girardeau map for this walk)",
   min(tv(P_hc, P_F), tv(P_hc, P_B), tv(X_hc, X_F), tv(X_hc, X_B)) > 0.1)
ok(f"3c. why: the free fermion has weight {same_site(P_F):.3f} on two walkers at one point with opposite senses (the head-on moment) and exactly 0 on one state; "
   "the bounce has 0 on both: it forbids two at one point, the rule forbids two in one state",
   same_site(P_F) > 0.04 and same_state(P_F) < 1e-28 and same_site(P_hc) == 0 and same_state(P_hc) == 0)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
