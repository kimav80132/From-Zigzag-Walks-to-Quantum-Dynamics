"""VI.3 (A), (D): two identical walkers in one packet, 3+1 - the run.   Recreated 9 Nov 2026 against master v34; the text took its numbers from this
run in v35.  Check 1' added 10 Nov.

Setup (this script's; the run behind v34's table was not recorded with its parameters).  Units hbar = c = 1, m = 1.
  - one Gaussian packet f ~ exp(-r^2/2w^2) at rest; walker 1 with axis +x, walker 2 with axis at theta_12 from it (180 degrees unless said);
    each column is the packet times its axis column in the upper pair of the mass-diagonal basis, beta = diag(1,1,-1,-1),
    alpha_k = [[0, sigma_k],[sigma_k, 0]]; the two-walker column is the odd combination of VI.1 (D);
  - free evolution for the time t = 1.5 m w^2, either exactly (Dirac, per momentum) or by IV.10's three-coin step (x-coin, y-coin, z-coin,
    then cos(m eps) - i beta sin(m eps); lattice spacing = eps, N = t/eps ticks);
  - detectors: L is z < -w, R is z > w; a setting sorts the absorbed column by Sigma along n = (cos a, sin a, 0), Sigma_k = diag(sigma_k, sigma_k);
  - weights by III.15a's rule on the two-walker event, VI.3 (2); S maximized over the four angles.
The rows are labelled by 1/(m w), the packet's wavenumber scale over the mass.
 1. the table: S for odd, even, product; E(0, b) = -kappa cos b; odd = even; product = odd/sqrt2 <= 2; S -> 2 sqrt2 as 1/(mw) -> 0.
 2. the walk's step against its continuum limit at 1/(mw) = 0.5, t = 3.2/m, with m eps = 0.4, 0.2, 0.1 on the same lattices: the distance
    between the columns is first order in the step and the gap in S closes.  (At the table's own times a step of m eps = 0.4 is too
    coarse: the note printed after check 2 gives the number, which is the one the text quotes since v35.)
 3. no signalling: the marginal at L is (1/2, 1/2) for every setting at R.
 4. (A): walkers launched with axes 180, 60 and 135 degrees apart give the same S; only the weight of the pair changes, as sin^2(theta_12/2)."""
import sys
import numpy as np
from scipy.optimize import minimize
fails = 0
def ok(name, cond):
    global fails
    print(("PASS  " if cond else "FAIL  ") + name); fails += (not cond)
s0 = np.eye(2); sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
Z2 = np.zeros((2, 2))
beta = np.diag([1, 1, -1, -1]).astype(complex)
alpha = [np.block([[Z2, s], [s, Z2]]) for s in (sx, sy, sz)]
Sigma = [np.block([[s, Z2], [Z2, s]]) for s in (sx, sy, sz)]
def chi_of(theta12):
    """two axis columns: chi_1 along +x, chi_2 at angle theta12 from it in the x-y plane (eigencolumns of n.sigma)"""
    col = lambda ang: np.array([1, np.exp(1j*ang)])/np.sqrt(2)          # eigencolumn of cos(ang) sx + sin(ang) sy, eigenvalue +1
    return col(0.0), col(theta12)
def grids(n, h):
    p = 2*np.pi*np.fft.fftfreq(n, h); PX, PY, PZ = np.meshgrid(p, p, p, indexing='ij')
    x = (np.arange(n) - n//2)*h; X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    return (PX, PY, PZ), (X, Y, Z)
def packet_hat(n, h, w, XYZ):
    X, Y, Z = XYZ; f = np.exp(-(X**2 + Y**2 + Z**2)/(2*w*w)); f /= np.sqrt((f**2).sum())
    return np.fft.fftn(np.fft.ifftshift(f)), f
def exact(fh, P, chi, t, m=1.0):
    PX, PY, PZ = P; E = np.sqrt(PX**2 + PY**2 + PZ**2 + m*m); ce, se = np.cos(E*t), np.sin(E*t)/E
    up = [(ce - 1j*m*se)*fh*chi[0], (ce - 1j*m*se)*fh*chi[1]]
    sp_chi = [PZ*chi[0] + (PX - 1j*PY)*chi[1], (PX + 1j*PY)*chi[0] - PZ*chi[1]]         # (sigma.p) chi
    lo = [-1j*se*fh*sp_chi[0], -1j*se*fh*sp_chi[1]]
    return np.stack([np.fft.fftshift(np.fft.ifftn(c)) for c in up + lo], axis=-1)
def automaton(fh, P, chi, N, eps, m=1.0):
    """IV.10's step: an x-coin, a y-coin, a z-coin (each exp(-i eps p_j alpha_j) on a lattice of spacing eps), then cos(m eps) - i beta sin(m eps)"""
    PX, PY, PZ = P; shape = PX.shape
    U = np.broadcast_to(np.eye(4, dtype=complex), shape + (4, 4)).copy()
    for Pj, aj in zip((PX, PY, PZ), alpha):
        Cj = np.cos(eps*Pj)[..., None, None]*np.eye(4) - 1j*np.sin(eps*Pj)[..., None, None]*aj
        U = Cj @ U
    U = (np.cos(m*eps)*np.eye(4) - 1j*np.sin(m*eps)*beta) @ U
    R = np.broadcast_to(np.eye(4, dtype=complex), shape + (4, 4)).copy(); B = U; k = N
    while k:
        if k & 1: R = B @ R
        B = B @ B; k >>= 1
    u0 = np.array([chi[0], chi[1], 0, 0])
    out = (R @ u0)*fh[..., None]
    return np.stack([np.fft.fftshift(np.fft.ifftn(out[..., i])) for i in range(4)], axis=-1)
def region_mats(phis, Z, d):
    """sums over the two detector regions, z < -d and z > d.  A grid plane lying exactly on an edge counts half (trapezoid rule); run() also places
    the edge midway between two grid planes when it is free to choose the spacing.  Either way the region sums are second order in the spacing.
    (10 Nov 2026: the first version of this script summed whole grid planes with the edge anywhere in a cell, an error of first order in the
    spacing that moved the third decimal of the table: 1.716, 2.490, 2.786 for what are 1.713, 2.486, 2.785.)"""
    tol = 1e-9*max(1.0, abs(d)); zs = Z[0, 0, :]
    wL = np.where(zs < -d - tol, 1.0, np.where(np.abs(zs + d) <= tol, 0.5, 0.0)); wR = np.where(zs > d + tol, 1.0, np.where(np.abs(zs - d) <= tol, 0.5, 0.0))
    Rg = {}
    for name, wz in (('L', wL), ('R', wR)):
        keep = wz > 0; cols = [ph[:, :, keep, :]*np.sqrt(wz[keep])[None, None, :, None] for ph in phis]
        cols = [c_.reshape(-1, 4) for c_ in cols]
        Rg[name] = [[ci.T @ cj.conj() for cj in cols] for ci in cols]       # R^{ij} = sum phi_i phi_j^dagger
    return Rg
def Pax(ang, sgn):
    return 0.5*(np.eye(4) + sgn*(np.cos(ang)*Sigma[0] + np.sin(ang)*Sigma[1]))
def weights(Rg, kind, a, b):
    out = {}
    for sa in (1, -1):
        Pa = Pax(a, sa)
        for sb in (1, -1):
            Pb = Pax(b, sb); tr = lambda P_, M_: np.trace(P_ @ M_)
            D2 = (tr(Pa, Rg['L'][0][0])*tr(Pb, Rg['R'][1][1])).real; E2 = (tr(Pa, Rg['L'][1][1])*tr(Pb, Rg['R'][0][0])).real
            X = (tr(Pa, Rg['L'][0][1])*tr(Pb, Rg['R'][1][0])).real
            out[(sa, sb)] = {'odd': D2 + E2 - 2*X, 'even': D2 + E2 + 2*X, 'product': D2}[kind]
    return out
def corr(Rg, kind, a, b):
    W = weights(Rg, kind, a, b); return sum(sa*sb*w for (sa, sb), w in W.items())/sum(W.values())
def Smax(Rg, kind, seeds=8):
    S = lambda g: abs(corr(Rg, kind, g[0], g[2]) - corr(Rg, kind, g[0], g[3]) + corr(Rg, kind, g[1], g[2]) + corr(Rg, kind, g[1], g[3]))
    best = 0.0
    for seed in range(seeds):
        r = minimize(lambda g: -S(g), np.random.default_rng(seed).uniform(0, 2*np.pi, 4), method='Nelder-Mead', options=dict(xatol=1e-9, fatol=1e-12, maxiter=6000))
        best = max(best, -r.fun)
    return best

def run(w, n, h=None, theta12=np.pi, step_eps=None, tfac=1.5, d=None):
    t = tfac*w*w
    if step_eps is not None: h = step_eps; N = int(round(t/step_eps)); t = N*step_eps
    dd = w if d is None else d
    if h is None:
        h = 14*np.sqrt(w*w/2 + (tfac*w)**2/2)/n
        if dd > 0: h = dd/(max(0, round(dd/h - 0.5)) + 0.5)               # the detector edge midway between two grid planes
    P, XYZ = grids(n, h); fh, f = packet_hat(n, h, w, XYZ); chis = chi_of(theta12)
    phis = [automaton(fh, P, c, N, step_eps) if step_eps is not None else exact(fh, P, c, t) for c in chis]
    edge = sum((abs(p[np.abs(XYZ[2]) > 0.45*n*h])**2).sum() for p in phis)/2
    return region_mats(phis, XYZ[2], w if d is None else d), dict(t=t, h=h, n=n, edge=edge, norm=sum((abs(p)**2).sum() for p in phis)/2)

# ---- 1. the table (exact evolution)
print("  1/(mw)    odd     even    product   kappa = -E(0,0)   max|E(0,b) + kappa cos b|    (t = 1.5 m w^2, d = w)")
table = {}
conv = 0.0
for w in (2.0, 4.0, 12.5, 25.0):
    Rg64, _ = run(w, 64); Rg, info = run(w, 96)
    conv = max(conv, abs(corr(Rg64, 'odd', 0.0, 0.0) - corr(Rg, 'odd', 0.0, 0.0)))
    So, Se, Sp = Smax(Rg, 'odd'), Smax(Rg, 'even'), Smax(Rg, 'product'); kap = -corr(Rg, 'odd', 0.0, 0.0)
    form = max(abs(corr(Rg, 'odd', 0.0, b) + kap*np.cos(b)) for b in np.linspace(0, 2*np.pi, 13))
    table[w] = (So, Se, Sp, kap, form, info)
    print("   %.2f    %.3f   %.3f   %.3f      %.3f              %.1e          (grid %d^3, spacing %.2f, weight near the box edge %.0e)"
          % (1/w, So, Se, Sp, kap, form, info['n'], info['h'], info['edge']))
ok("1a. the table does not depend on the grid: kappa on 64^3 and on 96^3 differ by at most %.0e" % conv, conv < 1e-3)
ok("1. odd = even in every row (largest difference %.0e); product = odd/sqrt2 (%.0e) and <= 2; S = 2 sqrt2 kappa (%.0e); E(0, b) = -kappa cos b (%.0e); "
   "S rises toward 2 sqrt2 = %.3f as 1/(mw) -> 0: %.3f, %.3f, %.3f, %.3f"
   % (max(abs(r[0] - r[1]) for r in table.values()), max(abs(r[2] - r[0]/np.sqrt(2)) for r in table.values()),
      max(abs(r[0] - 2*np.sqrt(2)*r[3]) for r in table.values()), max(r[4] for r in table.values()), 2*np.sqrt(2), *[table[w][0] for w in (2.0, 4.0, 12.5, 25.0)]),
   max(abs(r[0] - r[1]) for r in table.values()) < 1e-6 and max(abs(r[2] - r[0]/np.sqrt(2)) for r in table.values()) < 1e-6
   and max(abs(r[0] - 2*np.sqrt(2)*r[3]) for r in table.values()) < 1e-6 and max(r[4] for r in table.values()) < 1e-9
   and table[2.0][0] < table[4.0][0] < table[12.5][0] < table[25.0][0] < 2*np.sqrt(2) and table[25.0][0] > 2.8 and table[2.0][0] < 2)

# ---- 1'. (added 10 Nov 2026) the text's "the numbers depend on that time and on where the detectors sit": 1/(mw) = 0.5, detectors at d = 0
kap_t = {}
for tf in (1.0, 1.5, 2.0, 3.0, 4.0, 8.0):
    kap_t[tf] = tuple(-corr(run(2.0, n_, tfac=tf, d=0.0)[0], 'odd', 0.0, 0.0) for n_ in (96, 128))
ok("1'. the setup is part of the table: at 1/(mw) = 0.5 with the detectors at d = 0 the coefficient of cos b is %s at t = 1, 1.5, 2, 3, 4, 8 m w^2 "
   "(96^3 and 128^3 agree to %.0e): it is not monotonic in t; with d = w at t = 1.5 m w^2 it is the table's %.3f"
   % (", ".join("%.3f" % kap_t[tf][1] for tf in (1.0, 1.5, 2.0, 3.0, 4.0, 8.0)), max(abs(v[0] - v[1]) for v in kap_t.values()), table[2.0][3]),
   max(abs(v[0] - v[1]) for v in kap_t.values()) < 2e-3 and kap_t[2.0][1] > kap_t[1.0][1] + 0.1 and kap_t[2.0][1] > kap_t[8.0][1] + 0.03)

# ---- 2. the three-coin step against the exact evolution, on the same lattice (spacing = eps), 1/(mw) = 0.5, t = 3.2, d = w
print("   1/(mw) = 0.5, t = 3.2/m:   m eps    ticks    S, the step    S, exact on the same lattice    |column(step) - column(exact)|")
res = []
for n, eps in ((32, 0.4), (64, 0.2), (128, 0.1)):
    w = 2.0; N = int(round(3.2/eps)); P, XYZ = grids(n, eps); fh, f = packet_hat(n, eps, w, XYZ); chis = chi_of(np.pi)
    pa = [automaton(fh, P, c_, N, eps) for c_ in chis]; pe = [exact(fh, P, c_, N*eps) for c_ in chis]
    err = np.sqrt(sum((abs(x_ - y_)**2).sum() for x_, y_ in zip(pa, pe))/2)
    Sa, Se_ = Smax(region_mats(pa, XYZ[2], w), 'odd'), Smax(region_mats(pe, XYZ[2], w), 'odd')
    res.append((eps, N, Sa, Se_, err, abs(sum((abs(p_)**2).sum() for p_ in pa)/2 - 1)))
    print("                              %.1f      %3d       %.3f              %.3f                         %.3f" % (eps, N, Sa, Se_, err))
    del pa, pe
gaps = [abs(r[2] - r[3]) for r in res]; errs = [r[4] for r in res]
ok("2. the walk's step and its continuum limit: the distance between the two columns halves with the step (%.3f, %.3f, %.3f: first order), and "
   "the gap in S closes: %.3f, %.3f, %.3f at m eps = 0.4, 0.2, 0.1; the step keeps the length to %.0e"
   % (errs[0], errs[1], errs[2], gaps[0], gaps[1], gaps[2], max(r[5] for r in res)),
   abs(errs[0]/errs[1] - 2) < 0.2 and abs(errs[1]/errs[2] - 2) < 0.2 and gaps[0] > gaps[1] > gaps[2] and gaps[2] < 0.03 and max(r[5] for r in res) < 1e-12)
# and at the table's own time the coarse step is far from its limit: say so with the number
Ra, ia = run(4.0, 96, step_eps=0.4); Re_, ie = run(4.0, 96, h=0.4, tfac=ia['t']/16.0)
print("   NOTE  at the table's time, 1/(mw) = 0.25, t = 24/m = 60 ticks of m eps = 0.4: the step gives S = %.3f, the exact evolution on the same "
      "lattice %.3f. The step's error grows with the number of ticks; m eps = 0.4 is too coarse for 60 of them." % (Smax(Ra, 'odd'), Smax(Re_, 'odd')))

# ---- 3. no signalling, on the 1/(mw) = 0.25 row
Rg, info = run(4.0, 96)
marg = []
for a_ in (0.0, 0.9):
    for b_ in (0.0, 0.7, 2.1, 4.0):
        W = weights(Rg, 'odd', a_, b_); tot = sum(W.values()); marg.append(((W[(1, 1)] + W[(1, -1)])/tot, (W[(-1, 1)] + W[(-1, -1)])/tot))
dev = max(abs(m_ - 0.5) for pr in marg for m_ in pr)
ok("3. no signalling: the marginal at L is (0.5, 0.5) for two settings at L and four at R, to %.0e" % dev, dev < 1e-9)

# ---- 4. three launch axes
vals = {}
for deg in (180, 60, 135):
    Rg, info = run(12.5, 96, theta12=np.radians(deg))
    Wt = sum(weights(Rg, 'odd', 0.3, 1.2).values()); vals[deg] = (Smax(Rg, 'odd', seeds=12), Wt)
ratio = {deg: vals[deg][1]/vals[180][1] for deg in vals}
ok("4. launched with axes 180, 60, 135 degrees apart: S = %.4f, %.4f, %.4f - the same; the weight of the pair relative to 180 degrees is %.4f, %.4f "
   "against sin^2(theta_12/2) = %.4f, %.4f" % (vals[180][0], vals[60][0], vals[135][0], ratio[60], ratio[135], np.sin(np.radians(30))**2, np.sin(np.radians(67.5))**2),
   max(abs(vals[d_][0] - vals[180][0]) for d_ in (60, 135)) < 1e-4 and abs(ratio[60] - 0.25) < 1e-9 and abs(ratio[135] - np.sin(np.radians(67.5))**2) < 1e-9)
print("ALL PASS" if not fails else "%d FAILED" % fails); sys.exit(1 if fails else 0)
