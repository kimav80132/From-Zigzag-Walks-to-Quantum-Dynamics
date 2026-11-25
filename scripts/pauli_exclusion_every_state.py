"""Rebuilt 9 Nov 2026 against master v34 (check 4: the pair is prepared one tick before the head-on point).
VI.1 (D), revised after the audit of 2026-11-01. The two-walker column C(I,J), I = (site, sense); free step C -> U C U^T; exchange C -> C^T.
1. Exclusion in the site-sense basis alone does not choose the sign: the symmetric column C = (11^T - 1) (x) 1_2 (two walkers in the same sense,
   never on one site, uniform) is kept by the free step exactly, tick after tick, and has nothing on any same-site pair (the audit's counterexample).
2. Exclusion in every state, v^T C v = 0 for every column v, gives C^T = -C by polarization (checked on random real and complex C: the symmetric part
   of C is exactly what v^T C v reads; it is zero for all v iff the symmetric part vanishes).
3. The free step keeps C^T = -C, and keeps C^T = +C: it does not choose; the weight of each sector is conserved (U orthogonal).
4. The even column of the old argument (two walkers meeting head-on at one site, opposite senses, symmetrized) is fed into same-site same-sense
   components by the coin at that site; the odd one is not: the old argument was right about that column and wrong as a uniqueness claim."""
import sys
import numpy as np
fails = 0
def ok(n, b):
    global fails
    print(("PASS  " if b else "FAIL  ") + n); fails += (not b)
rng = np.random.default_rng(3)
def step(L, th=0.37):
    R = np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]]); S = np.zeros((2*L, 2*L))
    for x in range(L): S[2*((x + 1) % L), 2*x] = 1; S[2*((x - 1) % L) + 1, 2*x + 1] = 1
    return np.kron(np.eye(L), R) @ S
res = []
for L in (3, 5, 8):
    U = step(L); C = np.kron(np.ones((L, L)) - np.eye(L), np.eye(2)); C /= np.linalg.norm(C)
    same_site = max(abs(C[2*x:2*x+2, 2*x:2*x+2]).max() for x in range(L))
    kept = all(np.allclose(np.linalg.matrix_power(U, n) @ C @ np.linalg.matrix_power(U, n).T, C) for n in (1, 2, 7))
    res.append(np.allclose(C, C.T) and same_site == 0 and kept)
ok("1. a symmetric two-walker column with nothing on any same-site pair, kept exactly by the free step (L = 3, 5, 8): exclusion in one basis admits the even sector", all(res))
def polar(C, complex_=False):
    n = C.shape[0]; vs = [rng.normal(size=n) + (1j*rng.normal(size=n) if complex_ else 0) for _ in range(200)]
    reads = max(abs(v @ C @ v) for v in vs); return reads
C = rng.normal(size=(6, 6)); A = C - C.T; Sy = C + C.T
Cc = rng.normal(size=(6, 6)) + 1j*rng.normal(size=(6, 6)); Ac = Cc - Cc.T
ok("2. v^T C v = 0 for all v  <=>  C^T = -C: an antisymmetric C reads 0 for every v (real %.1e, complex %.1e); a symmetric part is read (%.2f); and e_I, e_I + e_J give C_IJ + C_JI = 0" % (
   polar(A), polar(Ac, True), polar(Sy)), polar(A) < 1e-12 and polar(Ac, True) < 1e-12 and polar(Sy) > 0.1)
U = step(5); Cr = rng.normal(size=(10, 10)); Ao, Ae = Cr - Cr.T, Cr + Cr.T
ok("3. the free step keeps both parities, (U C U^T)^T = U C^T U^T, and conserves the weight of each: it does not choose the sign",
   np.allclose((U @ Ao @ U.T).T, -(U @ Ao @ U.T)) and np.allclose((U @ Ae @ U.T).T, U @ Ae @ U.T) and abs(np.linalg.norm(U @ Ao @ U.T) - np.linalg.norm(Ao)) < 1e-12)
# 4. the head-on point.  One tick is the shift, then the coin; so the pair that MEETS at site 2 is (1,+; 3,-): after the shift both
#    are on site 2 with opposite senses, (x,+; x,-), and the coin there acts on both.  (The earlier draft started the pair already on
#    one site; the shift then separates it before any coin acts, and it printed 0.000.)
L = 5; th = 0.37; U = step(L, th); e = lambda x, s: np.eye(2*L)[2*x + s]
Ce = np.outer(e(1, 0), e(3, 1)) + np.outer(e(3, 1), e(1, 0)); Co = np.outer(e(1, 0), e(3, 1)) - np.outer(e(3, 1), e(1, 0))
fed = lambda C: max(abs(C[2*x + s, 2*x + s]) for x in range(L) for s in (0, 1))
Ue, Uo = U @ Ce @ U.T, U @ Co @ U.T
ok("4. head-on at site 2, the pair (1,+; 3,-) one tick before it meets: the even column is carried by the coin into the same-state entries "
   "(2,+; 2,+) and (2,-; 2,-), each of size |2 cos(th) sin(th)| = |sin 2th| = %.4f (found %.4f, %.4f); the odd column into none (%.1e); "
   "and the pair started ON one site, (2,+; 2,-), is separated by the shift before any coin acts (even %.1e, odd %.1e)"
   % (abs(np.sin(2*th)), abs(Ue[4, 4]), abs(Ue[5, 5]), fed(Uo),
      fed(U @ (np.outer(e(2, 0), e(2, 1)) + np.outer(e(2, 1), e(2, 0))) @ U.T), fed(U @ (np.outer(e(2, 0), e(2, 1)) - np.outer(e(2, 1), e(2, 0))) @ U.T)),
   abs(abs(Ue[4, 4]) - abs(np.sin(2*th))) < 1e-12 and abs(abs(Ue[5, 5]) - abs(np.sin(2*th))) < 1e-12 and fed(Uo) < 1e-12)
print("ALL PASS" if not fails else "%d FAILED" % fails); sys.exit(1 if fails else 0)
