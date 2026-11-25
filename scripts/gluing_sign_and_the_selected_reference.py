"""
Exact checks for drafts/gluing_the_route_to_itself_2026_11_07.md.
(1) VI.1's antisymmetric N-walker configurations, summed over N with the wedge (Laplace expansion), give det(1 + U) for any
    one-walker evolution U: Z = sum_N Tr[Lambda^N U] = det(1 + U). Checked by building every N-walker column explicitly.
(2) log det(1 + U) = sum_n (-1)^{n+1} Tr[U^n]/n: the closed history winding n times carries the cyclic exchange sign (-1)^{n-1}.
    The symmetric (undirected, Part II) counting gives 1/det(1 - U) and all plus signs.
(3) Real time, U = e^{-iHT}: Z(T) = prod_s (1 + e^{-iE_s T}) - at rest 2 + 2cos aT - no term dominates.
    Euclidean parameter (III.16), U = e^{-beta H}: Z(beta) = prod_s (1 + e^{-beta E_s}); as beta grows the configuration with
    every E_s < 0 slot filled dominates and -(1/beta) log Z -> sum_{E_s<0} E_s = -1/2 Tr|H|. The symmetric counting has no limit.
Nothing solved or approximated: finite matrices, exact linear algebra; the limit is displayed at increasing beta.
"""
import itertools, sys, math
from collections import Counter
import numpy as np
import scipy.linalg as sl
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all = ok_all and bool(cond); print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()
rng = np.random.default_rng(3)

def wedge_trace(U, N):                        # Tr over the antisymmetric N-walker space = sum over N-subsets of principal minors
    d = U.shape[0]; tot = 0.0
    for S in itertools.combinations(range(d), N):
        tot += np.linalg.det(U[np.ix_(S, S)])  # <S| Lambda^N U |S> = det of the N x N block (the Laplace expansion of the wedge)
    return tot
def sym_trace(U, N):                          # the symmetric (bosonic) count: permanents over multisets
    d = U.shape[0]; tot = 0.0
    for S in itertools.combinations_with_replacement(range(d), N):
        M = U[np.ix_(S, S)]; per = sum(np.prod([M[i, p[i]] for i in range(N)]) for p in itertools.permutations(range(N)))
        norm = np.prod([math.factorial(cnt) for cnt in Counter(S).values()])
        tot += per/norm
    return tot

# (1) the generating function of closed antisymmetric configurations
d = 4; A = rng.normal(size=(d, d)); H = (A + A.T)/2; T = 0.7
U = sl.expm(-1j*H*T)
Z_conf = sum(wedge_trace(U, N) for N in range(d + 1)); Z_det = np.linalg.det(np.eye(d) + U)
ok("1. sum over every antisymmetric N-walker closed configuration (N = 0..%d, built from the minors of U) = %.10f%+.10fi; det(1 + U) = %.10f%+.10fi"
   % (d, Z_conf.real, Z_conf.imag, Z_det.real, Z_det.imag), abs(Z_conf - Z_det) < 1e-10)

# (2) the winding expansion and its signs
logZ = sum((-1)**(n + 1)*np.trace(np.linalg.matrix_power(U*0.5, n))/n for n in range(1, 200))   # scaled U so the series converges
ok("2. log det(1 + U) = sum_n (-1)^(n+1) Tr[U^n]/n to %.1e (U scaled by 1/2 for convergence): the n-fold closed history carries the cyclic "
   "exchange sign (-1)^(n-1); the symmetric count gives 1/det(1 - U), all plus: log = sum_n Tr[U^n]/n to %.1e"
   % (abs(logZ - np.log(np.linalg.det(np.eye(d) + U*0.5))),
      abs(sum(np.trace(np.linalg.matrix_power(U*0.5, n))/n for n in range(1, 200)) + np.log(np.linalg.det(np.eye(d) - U*0.5)))),
   abs(logZ - np.log(np.linalg.det(np.eye(d) + U*0.5))) < 1e-10
   and abs(sum(np.trace(np.linalg.matrix_power(U*0.5, n))/n for n in range(1, 200)) + np.log(np.linalg.det(np.eye(d) - U*0.5))) < 1e-10)   # 10 Nov 2026: the symmetric half was printed, not tested
# the symmetric generating function from permanents, small case
Us = U*0.3; Zs = sum(sym_trace(Us, N) for N in range(0, 7)); Zs_exact = 1/np.linalg.det(np.eye(d) - Us)
ok("   symmetric configurations summed (N = 0..6, permanents) = %.6f%+.6fi vs 1/det(1 - U) = %.6f%+.6fi (series truncated)"
   % (Zs.real, Zs.imag, Zs_exact.real, Zs_exact.imag), abs(Zs - Zs_exact) < 1e-3)

# (3) the rest mode: real time vs the Euclidean parameter
a = 1.0; Hr = np.array([[0, 1j*a], [-1j*a, 0]])        # the user's box at rest, hbar = 1: eigenvalues +-a
for Tt in (0.4, 1.3):
    Zt = np.linalg.det(np.eye(2) + sl.expm(-1j*Hr*Tt))
    ok("3a. real time T = %.1f: det(1 + e^{-iHT}) = %.6f%+.6fi = 2 + 2cos(aT) = %.6f (every configuration a phase, none dominates)"
       % (Tt, Zt.real, Zt.imag, 2 + 2*np.cos(a*Tt)), abs(Zt - (2 + 2*np.cos(a*Tt))) < 1e-12)
rows = []
for beta in (1.0, 3.0, 10.0, 30.0):
    Zb = np.linalg.det(np.eye(2) + sl.expm(-beta*Hr)).real
    terms = {"N=0": 1.0, "N=1, + slot": np.exp(-beta*a), "N=1, - slot": np.exp(+beta*a), "N=2": 1.0}
    E0 = -np.log(Zb)/beta
    rows.append((beta, Zb, terms["N=1, - slot"]/Zb, E0))
print("      Euclidean parameter: Z(beta) = 1 + e^{-beta a} + e^{+beta a} + 1; share of the '- slot filled' configuration; -(1/beta) log Z")
for beta, Zb, share, E0 in rows: print("        beta = %5.1f   Z = %12.4f   share(- filled) = %.6f   -(1/beta) log Z = %+.6f  (target -a = -1)" % (beta, Zb, share, E0))
ok("3b. as beta grows the configuration with the - slot filled takes the whole weight and -(1/beta) log Z -> -a = -1/2 Tr|H| (%+.6f at beta = 30)"
   % rows[-1][3], abs(rows[-1][3] + a) < 1e-3 and rows[-1][2] > 0.999)
# the symmetric count at rest has no limit
ok("3c. the symmetric (undirected) count at rest: 1/(1 - e^{+beta a}) is negative for every beta > 0 (%.4f at beta = 1): the series over "
   "occupations of the - slot does not converge - no reference exists without the exchange sign" % (1/(1 - np.exp(a*1.0))), 1/(1 - np.exp(a)) < 0)

# (3d) with motion: several k, H_k = [[ck, i a],[-i a, -ck]]
ks = np.array([-1.5, -0.5, 0.0, 0.5, 1.5]); c = 1.0
Hk = sl.block_diag(*[np.array([[c*k, 1j*a], [-1j*a, -c*k]]) for k in ks]); w = np.sqrt(c*c*ks**2 + a*a)
for beta in (10.0, 20.0):
    logZ = 0.0                                   # H is block-diagonal in k: det(1 + e^{-beta H}) is the product of the 2x2 block determinants
    for k in ks:
        Hb = np.array([[c*k, 1j*a], [-1j*a, -c*k]]); Ub = sl.expm(-beta*Hb)
        # det(1 + Ub) = 1 + Tr Ub + det Ub exactly (2x2); written through log1p to keep the large and small terms apart
        tr, dt = np.trace(Ub).real, 1.0                      # det e^{-beta H_b} = e^{-beta Tr H_b} = 1 exactly (Tr H_b = 0)
        logZ += np.log(tr) + np.log1p((1 + dt)/tr)
    E0 = -logZ/beta
ok("3d. five k's: -(1/beta) log det(1 + e^{-beta H}) = %+.6f at beta = 20; -sum_k hbar omega_k = %+.6f; -1/2 Tr|H| = %+.6f"
   % (E0, -w.sum(), -0.5*np.abs(np.linalg.eigvalsh(Hk)).sum()), abs(E0 + w.sum()) < 1e-6)

# (4) the x-clock step for readings constant in t (nu = 0, NOT the k = 0 rest mode) is the symmetric exchange, c dQ/dx = a sigma_x Q;
#     its exponential has eigenvalues e^{+-ax/c}, the same as e^{-(x/c) H(k=0)}: an eigenvalue/determinant identity between two
#     different solutions (draft section 6). The moving statement uses M_k below.
sx = np.array([[0, 1.0], [1.0, 0]])
for x in (3.0, 10.0):
    Ux = sl.expm((a*x/c)*sx)                                  # the x-clock evolution of the loss-out counters over a distance x
    Zx = np.linalg.det(np.eye(2) + Ux); Zb = np.linalg.det(np.eye(2) + sl.expm(-(x/c)*Hr)).real
    ev_x = np.sort(np.linalg.eigvals(Ux).real); ev_b = np.sort(np.linalg.eigvals(sl.expm(-(x/c)*Hr)).real)
    ok("4. x-clock, distance x = %4.1f: e^{(ax/c) sigma_x} has eigenvalues %s = those of e^{-(x/c) H} %s; closed-history count det(1 + U_x) = %.6f "
       "= det(1 + e^{-beta H}) at beta = x/c = %.6f; -(c/x) log = %+.6f" % (x, np.round(ev_x, 4), np.round(ev_b, 4), Zx, Zb, -np.log(Zx)/(x/c)),
       np.allclose(ev_x, ev_b) and abs(Zx - Zb) < 1e-9*Zb)

# (5) the moving continuation nu = i c k on the x-clock matrix: M_k = [[-ck, a],[a, ck]] is real, M_k^2 = omega_k^2, and
#     det(1 + e^{s M_k}) = 2 + 2 cosh(omega_k s) -> -(1/s) log -> -omega_k (the other AI's note, checked here on the matrices)
for k in (0.0, 0.7, 2.0):
    Mk = np.array([[-c*k, a], [a, c*k]]); wk = np.sqrt(a*a + c*c*k*k)
    s_ = 12.0; Zs = 2 + 2*np.cosh(wk*s_); E_ = sl.expm(s_*Mk); Zd = 1 + np.trace(E_) + 1.0   # det(1 + E) = 1 + Tr E + det E, det E = e^{s Tr M_k} = 1
    ok("5. k = %.1f: M_k^2 = omega_k^2 to %.0e; det(1 + e^{s M_k}) = 2 + 2cosh(omega s) to %.0e (relative); -(1/s) log = %+.6f vs -omega_k = %+.6f"
       % (k, abs(Mk @ Mk - wk*wk*np.eye(2)).max(), abs(Zd - Zs)/Zs, -np.log(Zd)/s_, -wk),
       abs(Mk @ Mk - wk*wk*np.eye(2)).max() < 1e-12 and abs(Zd - Zs)/Zs < 1e-9 and abs(-np.log(Zd)/s_ + wk) < 1e-3)

# (6) the full note's three items: D = diag(1, i) (III.9's relabelling) conjugates e^{s M_k} to e^{-s h(k)}, the growing eigenvector
#     of M_k maps to the negative-energy eigenvector of h(k); the small-k rate omega_k - a -> (c^2/2a) k^2 = (hbar/2m) k^2 with a = mc^2/hbar
D = np.diag([1.0, 1j])
worst = 0.0; sector_ok = True
for k in (0.0, 0.7, 2.0):
    Mk = np.array([[-c*k, a], [a, c*k]]); hk = np.array([[c*k, 1j*a], [-1j*a, -c*k]])     # h(k)/hbar in the user's box
    worst = max(worst, abs(sl.expm(2.0*Mk) - np.linalg.inv(D) @ sl.expm(-2.0*hk) @ D).max())
    ev, V = np.linalg.eig(Mk); grow = V[:, np.argmax(ev.real)]                            # growing direction of the x-clock matrix
    mapped = D @ grow; e_neg = (hk @ mapped)/np.linalg.norm(mapped)**2                     # is D grow an eigenvector of h with eigenvalue -omega?
    lam = np.vdot(mapped, hk @ mapped)/np.vdot(mapped, mapped)
    sector_ok &= abs(hk @ mapped - lam*mapped).max() < 1e-12 and lam.real < 0
ok("6. e^{s M_k} = D^{-1} e^{-s h(k)} D with D = diag(1, i) (III.9's relabelling) to %.0e at k = 0, 0.7, 2; the growing eigenvector of M_k maps "
   "under D to an eigenvector of h(k) with negative eigenvalue (the - slot): %s" % (worst, sector_ok), worst < 1e-9 and sector_ok)
ks_small = np.array([0.05, 0.1, 0.2]); rate = np.sqrt(a*a + c*c*ks_small**2) - a; diff = (c*c/(2*a))*ks_small**2
ok("   small k: omega_k - a = %s against (c^2/2a) k^2 = %s (ratio %s -> 1 as k -> 0, the O(k^4) term): Part I's diffusion coefficient hbar/2m = c^2/2a"
   % (np.round(rate, 6), np.round(diff, 6), np.round(rate/diff, 4)), np.allclose(rate/diff, 1, atol=1.1e-2) and abs(rate[0]/diff[0] - 1) < 1e-3)
print("ALL PASS" if ok_all else "SOME FAILED"); raise SystemExit(0 if ok_all else 1)
