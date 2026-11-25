"""VI.1, remark after (D) (master v43, 11 Nov 2026): what determinism gives for exclusion and what it does not; N walkers as a volume element.

  1. A deterministic, reversible step on a finite list of states is a permutation.  It keeps an occupation that is 0 or 1 at 0 or 1 - and keeps an
     occupation of 2 at 2 just as well: it keeps an exclusion, it does not make one (Lynden-Bell's incompressible phase fluid, in the smallest setting).
  2. Volumes: for columns u, v of the walk on a ring, u ^ u = 0 and u ^ v = -(v ^ u); the step acts on the two-walker column as C -> U C U^T and keeps
     it antisymmetric; its matrix on such columns is the matrix of the 2 x 2 minors of U.  For N = 3 the same with the 3 x 3 minors.
  3. The step keeps the full volume: det U = 1, and the top minor (all slots filled) is 1.
  4. The pair of counters along the walk's clock: with the loss in (III.2) the area shrinks as e^{-2 a tau}; with the loss out (III.5a) the map is a
     squeeze of determinant 1 (Tr M_k = 0); read in t it is a rotation of determinant 1.
"""
import itertools, sys
import numpy as np, scipy.linalg as sl
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all = ok_all and bool(cond); print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()
rng = np.random.default_rng(5)

# ---- 1. a permutation of states
n = 12; perm = rng.permutation(n)
occ01 = np.zeros(n, int); occ01[rng.choice(n, 5, replace=False)] = 1
occ2 = occ01.copy(); occ2[np.flatnonzero(occ01)[0]] = 2
def step(occ):
    out = np.zeros_like(occ)
    for i in range(n): out[perm[i]] += occ[i]
    return out
a01, a2 = occ01.copy(), occ2.copy(); keep01 = True; keep2 = True
for _ in range(50):
    a01, a2 = step(a01), step(a2); keep01 = keep01 and a01.max() == 1 and a01.sum() == 5; keep2 = keep2 and a2.max() == 2
merge = np.zeros(n, int); i1, i2 = np.flatnonzero(occ01)[:2]
for i in range(n): merge[i1 if i == i2 else i] += occ01[i]                     # a step that sends two states to one: not reversible
ok("1. a reversible deterministic step on 12 states (a permutation), 50 steps: an occupation that is 0 or 1 stays 0 or 1 (%s), and one that starts with a 2 "
   "keeps its 2 (%s): the step keeps an exclusion and does not make one. A step that sends two states to one piles two walkers up (largest occupation %d)"
   % (keep01, keep2, merge.max()), keep01 and keep2 and merge.max() == 2)

# ---- 2. volumes: III.15a's tick on a ring
def tick(L, th):
    D = 2*L; U = np.zeros((D, D)); c, s = np.cos(th), np.sin(th)
    for y in range(L):
        yp, ym = (y - 1) % L, (y + 1) % L
        U[2*y, 2*yp] += c; U[2*y, 2*ym + 1] += s; U[2*y + 1, 2*yp] += -s; U[2*y + 1, 2*ym + 1] += c
    return U
U = tick(3, 0.6); D = 6
wedge2 = lambda u, v: np.outer(u, v) - np.outer(v, u)
u, v, w3 = rng.normal(size=D), rng.normal(size=D), rng.normal(size=D)
same = np.abs(wedge2(u, u)).max(); anti = np.abs(wedge2(u, v) + wedge2(v, u)).max()
stepC = np.abs(wedge2(U@u, U@v) - U@wedge2(u, v)@U.T).max()
pairs = list(itertools.combinations(range(D), 2))
comp2 = np.array([[np.linalg.det(U[np.ix_(r, cl)]) for cl in pairs] for r in pairs])    # the 2 x 2 minors of U
vec = lambda C: np.array([C[i, j] for i, j in pairs])
minors2 = np.abs(vec(wedge2(U@u, U@v)) - comp2@vec(wedge2(u, v))).max()
def wedge3(a, b, c_):
    T = np.einsum('i,j,k->ijk', a, b, c_); out = np.zeros_like(T)
    for p_, sg in (((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)): out += sg*np.transpose(T, p_)
    return out
trip = list(itertools.combinations(range(D), 3))
comp3 = np.array([[np.linalg.det(U[np.ix_(r, cl)]) for cl in trip] for r in trip])
vec3 = lambda W: np.array([W[i, j, k] for i, j, k in trip])
minors3 = np.abs(vec3(wedge3(U@u, U@v, U@w3)) - comp3@vec3(wedge3(u, v, w3))).max(); same3 = np.abs(wedge3(u, v, u)).max()
ok("2. volumes, on a ring of 3 sites (6 slots): u ^ u = 0 (%.0e) and u ^ v = -(v ^ u) (%.0e); the step takes u ^ v to (Uu) ^ (Uv) = U (u ^ v) U^T (%.0e), "
   "which on the 15 pairs of slots is the matrix of the 2 x 2 minors of U (%.0e); three walkers: u ^ v ^ u = 0 (%.0e), and the step is the matrix of the "
   "3 x 3 minors on the 20 triples (%.0e)" % (same, anti, stepC, minors2, same3, minors3),
   same == 0 and anti == 0 and stepC < 1e-13 and minors2 < 1e-13 and same3 < 1e-13 and minors3 < 1e-13)

# ---- 3. the full volume
ok("3. the step keeps volume: U^T U = 1 (%.0e), det U = %.12f, so the one filling with every slot occupied has weight 1"
   % (np.abs(U.T@U - np.eye(D)).max(), np.linalg.det(U)), np.abs(U.T@U - np.eye(D)).max() < 1e-13 and abs(np.linalg.det(U) - 1) < 1e-12)

# ---- 4. the pair of counters, three readings
a, k, tau = 1.0, 0.7, 1.3; w = np.hypot(a, k)
loss_in = sl.expm(tau*a*np.array([[-1., 1.], [1., -1.]]))                     # III.2 at rest, probabilities
loss_out = sl.expm(tau*np.array([[-k, a], [a, k]]))                          # L_k(tau) = e^{tau M_k}
in_t = sl.expm(tau*a*np.array([[0., 1.], [-1., 0.]]))                        # the coin read in t at rest
ev = np.sort(np.linalg.eigvalsh(loss_out))
ok("4. the pair of counters over a clock a tau = 1.3: loss in, the area shrinks, det = %.6f = e^{-2 a tau} (%.0e); loss out, a squeeze by e^{-w tau}, e^{+w tau} "
   "= %.4f, %.4f with det = %.12f (Tr M_k = 0); read in t, a rotation with det = %.12f" % (np.linalg.det(loss_in), abs(np.linalg.det(loss_in) - np.exp(-2*a*tau)),
   ev[0], ev[1], np.linalg.det(loss_out), np.linalg.det(in_t)),
   abs(np.linalg.det(loss_in) - np.exp(-2*a*tau)) < 1e-12 and abs(np.linalg.det(loss_out) - 1) < 1e-12 and abs(np.linalg.det(in_t) - 1) < 1e-12
   and abs(ev[1] - np.exp(w*tau)) < 1e-10 and np.abs(in_t.T@in_t - np.eye(2)).max() < 1e-13)

print("ALL PASS" if ok_all else "SOME FAILED"); raise SystemExit(0 if ok_all else 1)
