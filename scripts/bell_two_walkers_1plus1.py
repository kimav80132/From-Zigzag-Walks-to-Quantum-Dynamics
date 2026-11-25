"""VI.3 (C): in 1+1 the sense cannot do it.   Recreated 9 Nov 2026 against master v34.

The construction of VI.3 (B) with the SENSE as the two-valued thing.  Two identical walkers leave one point, one in each sense:
phi_1 launched forward (+), phi_2 launched backward (-); the free tick (shift, then the coin cos + sin T, angle 0.2) for 20 ticks;
detectors: L = more than 8 sites to the left of the source, R = more than 8 to the right.  A setting is a turn of the sense pair
(V.1's first knob: the coin's own rotation), i.e. the detector sorts the sense pair it absorbs along chi(alpha) = (cos alpha/2, sin alpha/2)
and the direction perpendicular to it.  The weights are III.15a's rule on the two-walker event,
    W(alpha, beta) = sum_{x in L, y in R} [ chi_a^T C(x, y) chi_b ]^2,     C = phi_1 phi_2^T -+ phi_2 phi_1^T  (odd / even),  or one term (product),
and S is maximized over the four angles.
  1. S for the odd state, the even state and the product state: all below 2, and alike - each within 0.005 of the product state's
     closed form 2 |r_L| |r_R|, r the Bloch vector of the sense matrix summed over a detector.
  2. why: the walker that reaches the left detector is in the backward sense with weight 0.96 - the sense IS the direction of motion, and
     the side has already read it; the single-site two-walker sense matrix C(x, y) has rank two (the exchange term is there), but its
     weight in the detected pairs is small.
The parameters are this script's; v34's text quotes S = 1.81 and the weight 0.95 from a run whose parameters were not recorded."""
import sys
import numpy as np
from scipy.optimize import minimize
fails = 0
def ok(name, cond):
    global fails
    print(("PASS  " if cond else "FAIL  ") + name); fails += (not cond)
Ls, N, th, d = 144, 20, 0.2, 8
c, s = np.cos(th), np.sin(th); x = np.arange(Ls) - Ls//2
def evolve(col):
    for _ in range(N):
        p, m = np.roll(col[:, 0], 1), np.roll(col[:, 1], -1); col = np.stack([c*p + s*m, -s*p + c*m], axis=1)
    return col
f = (x == 0).astype(float)
p1 = evolve(np.stack([f, 0*f], axis=1)); p2 = evolve(np.stack([0*f, f], axis=1))
Lm, Rm = x < -d, x > d
M = lambda reg, A, B: A[reg].T @ B[reg]                                # sum over the region of A(x) B(x)^T, a 2x2 matrix on the senses
ML = [[M(Lm, a_, b_) for b_ in (p1, p2)] for a_ in (p1, p2)]; MR = [[M(Rm, a_, b_) for b_ in (p1, p2)] for a_ in (p1, p2)]
q = lambda Mx, v: v @ Mx @ v
def weights(kind, A, B):
    out = {}
    for sa in (1, -1):
        al = np.array([np.cos(A/2), np.sin(A/2)]) if sa > 0 else np.array([-np.sin(A/2), np.cos(A/2)])
        for sb in (1, -1):
            be = np.array([np.cos(B/2), np.sin(B/2)]) if sb > 0 else np.array([-np.sin(B/2), np.cos(B/2)])
            D2 = q(ML[0][0], al)*q(MR[1][1], be)                        # walker 1 at L, walker 2 at R
            E2 = q(ML[1][1], al)*q(MR[0][0], be)                        # walker 2 at L, walker 1 at R
            X = q(ML[0][1], al)*q(MR[1][0], be)                         # the exchange term
            out[(sa, sb)] = {'odd': D2 + E2 - 2*X, 'even': D2 + E2 + 2*X, 'product': E2}[kind]
    return out
def corr(kind, A, B):
    W = weights(kind, A, B); return sum(sa*sb*w for (sa, sb), w in W.items())/sum(W.values())
def Smax(kind):
    S = lambda g: abs(corr(kind, g[0], g[2]) - corr(kind, g[0], g[3]) + corr(kind, g[1], g[2]) + corr(kind, g[1], g[3]))
    best = 0.0
    for seed in range(16):
        r = minimize(lambda g: -S(g), np.random.default_rng(seed).uniform(0, 2*np.pi, 4), method='Nelder-Mead',
                     options=dict(xatol=1e-10, fatol=1e-13, maxiter=20000, maxfev=20000))
        best = max(best, -r.fun)
    return best
def bloch(Mx):
    Mn = Mx/np.trace(Mx); return np.hypot(Mn[0, 0] - Mn[1, 1], 2*Mn[0, 1])
S_odd, S_even, S_prod = Smax('odd'), Smax('even'), Smax('product')
closed = 2*bloch(ML[1][1])*bloch(MR[0][0])
ok("1. S maximized over the four settings: odd %.3f, even %.3f, product %.3f - all below 2; the product state's closed form 2|r_L||r_R| = %.3f "
   "(|r_L| = %.4f, |r_R| = %.4f); the three within %.3f of it" % (S_odd, S_even, S_prod, closed, bloch(ML[1][1]), bloch(MR[0][0]),
   max(abs(S_odd - closed), abs(S_even - closed), abs(S_prod - closed))),
   max(S_odd, S_even, S_prod) < 2 and abs(S_prod - closed) < 1e-6 and max(abs(S_odd - closed), abs(S_even - closed)) < 0.005)
tot = ML[0][0] + ML[1][1]; wback = tot[1, 1]/np.trace(tot)
Wd = sum(weights('odd', 0.3, 1.1).values()); D2 = np.trace(ML[0][0])*np.trace(MR[1][1]); E2 = np.trace(ML[1][1])*np.trace(MR[0][0])
vec = lambda t_: np.array([np.cos(t_/2), np.sin(t_/2)])
Xrel = max(abs(2*q(ML[0][1], vec(t1))*q(MR[1][0], vec(t2))) for t1 in np.linspace(0, 2*np.pi, 73) for t2 in np.linspace(0, 2*np.pi, 73))/(D2 + E2)
wsite = (p1**2).sum(axis=1)*(p2**2).sum(axis=1)                     # a site both walkers reach
xs_ = np.where(Lm)[0][np.argmax(wsite[Lm])]; ys_ = np.where(Rm)[0][np.argmax(wsite[Rm])]
Cxy = np.outer(p1[xs_], p2[ys_]) - np.outer(p2[xs_], p1[ys_]); sv = np.linalg.svd(Cxy)[1]
ok("2. a walker absorbed at the left detector is in the backward sense with weight %.3f; of the detected pairs, 'the walker launched backward "
   "at L, the one launched forward at R' carries %.4f and the other pairing %.4f; the exchange term is at most %.3f of the total at any pair of "
   "settings; the two-walker sense matrix at a pair of sites (%d, %d) has rank %d (singular values %.1e, %.1e): the exchange term is there, and "
   "what it entangles has been read by the side"
   % (wback, E2/(D2 + E2), D2/(D2 + E2), Xrel, x[xs_], x[ys_], np.linalg.matrix_rank(Cxy, tol=1e-14), sv[0], sv[1]),
   wback > 0.9 and E2/(D2 + E2) > 0.99 and Xrel < 0.1 and sv[1] > 1e-6*sv[0] and sv[1] > 0)
print("ALL PASS" if not fails else "%d FAILED" % fails); sys.exit(1 if fails else 0)
