"""VI.3 (A) "the state from an interaction" and (B'): entanglement made by a scattering.   Rebuilt 9 Nov 2026 against master v34.

Two identical walkers come in from opposite sides, unentangled, and repel each other through V.2's step: the other walker's
potential is a phase per step on the joint state, a Gaussian of range 6 sites and height V0.  The axis is a spectator (V.4: the
field leaves the coin alone), so the two-walker amplitude is a spatial part F(x1 s1, x2 s2) times the two axis columns, and
VI.1's exchange makes the axis part at a pair of detectors (L, R)
        F(L,R) chi1 chi2^T - F(R,L) chi2 chi1^T :        "walker 1 bounced back to L"  minus  "walker 1 went through to R".
The lattice supplies three numbers, all sums over the detector regions and the senses:
        w_b = sum |F(L,R)|^2 (bounce),   w_p = sum |F(R,L)|^2 (pass),   X = sum F(L,R) conj F(R,L)  (complex).
For detector columns alpha (at L), beta (at R) the weight of the outcome is, by III.15a's rule on the two-walker event,
        W = w_b |<alpha|chi1>|^2 |<beta|chi2>|^2 + w_p |<alpha|chi2>|^2 |<beta|chi1>|^2
            - 2 Re[ X <alpha|chi1> <beta|chi2> conj(<alpha|chi2> <beta|chi1>) ].
  1. free walkers pass through each other: w_b = 0, S = 2.
  2. with the repulsion: S maximized over the four detector axes anywhere on the sphere, against (B')'s 2 sqrt(1 + C^2),
     C = 2|X|/(w_b + w_p).
  3. arg X = -pi/2 at every V0 (a symmetric barrier), so sorting in ONE plane reads the real part of an imaginary cross term: S = 2.
Lattice: 1+1, two walkers, each a two-sense column with the rotation coin; 320 sites, coin angle a = 0.5 per tick, packets of
width 14 at sites 100 and 220 with wavenumbers +-0.35, 300 ticks; detector regions: more than 20 sites left / right of the middle.
(The earlier draft of this script had a = 1.0, for which the packets hardly move in 300 ticks, and kept only Re X.)"""
import sys
import numpy as np
from scipy.optimize import minimize
fails = 0
def ok(name, cond):
    global fails
    print(("PASS  " if cond else "FAIL  ") + name); fails += (not cond)

L, a, k0, w, x1, x2, N, sig, gap = 320, 0.5, 0.35, 14.0, 100, 220, 300, 6.0, 20
c, s = np.cos(a), np.sin(a); x = np.arange(L)
def packet(x0, k):
    f = np.exp(-(x - x0)**2/(2*w*w))*np.exp(1j*k*x)
    Hm = np.array([[k, 1j*a], [-1j*a, -k]]); wv, vv = np.linalg.eigh(Hm); u = vv[:, np.argmax(wv)]   # positive-energy column, III.9's basis
    psi = f[:, None]*u[None, :]; return psi/np.linalg.norm(psi)
P0 = np.einsum('xs,yt->xsyt', packet(x1, +k0), packet(x2, -k0))
d = np.abs(x[:, None] - x[None, :]); d = np.minimum(d, L - d)
def scatter(V0):
    ph = np.exp(-1j*V0*np.exp(-d**2/(2*sig*sig)))[:, None, :, None]          # V.2 (F): the potential is a phase per step
    P = P0.copy()
    for n in range(N):
        A0 = np.roll(P[:, 0], 1, axis=0); A1 = np.roll(P[:, 1], -1, axis=0)      # walker 1: shift, then the coin
        P = np.stack([c*A0 + s*A1, -s*A0 + c*A1], axis=1)
        B0 = np.roll(P[..., 0], 1, axis=2); B1 = np.roll(P[..., 1], -1, axis=2)  # walker 2
        P = np.stack([c*B0 + s*B1, -s*B0 + c*B1], axis=3)*ph
    mid = (x1 + x2)//2; Lm = x < mid - gap; Rm = x > mid + gap
    F = P[Lm][:, :, Rm, :]                                                    # walker 1 at L, walker 2 at R: walker 1 bounced
    G = np.transpose(P[Rm][:, :, Lm, :], (2, 3, 0, 1))                         # the same configuration with the labels exchanged: walker 1 passed
    return (np.abs(F)**2).sum(), (np.abs(G)**2).sum(), (F*np.conj(G)).sum(), (np.abs(P)**2).sum()

chi1, chi2 = np.array([1.0, 0.0]), np.array([0.0, 1.0])                        # launched axes +x and -x; theta is measured from x
col = lambda th, ph: np.array([np.cos(th/2), np.exp(1j*ph)*np.sin(th/2)])
def corr(wb, wp, X, A, B):
    tot = num = 0.0
    for sa in (1, -1):
        al = col(A[0], A[1]) if sa > 0 else col(np.pi - A[0], A[1] + np.pi)    # the opposite axis
        for sb in (1, -1):
            be = col(B[0], B[1]) if sb > 0 else col(np.pi - B[0], B[1] + np.pi)
            a1, a2, b1, b2 = np.vdot(al, chi1), np.vdot(al, chi2), np.vdot(be, chi1), np.vdot(be, chi2)
            W = wb*abs(a1)**2*abs(b2)**2 + wp*abs(a2)**2*abs(b1)**2 - 2*(X*a1*b2*np.conj(a2*b1)).real
            tot += W; num += sa*sb*W
    return num/tot
def S_max(wb, wp, X, in_plane=False):
    def S(g):
        if in_plane: A, A2, B, B2 = (g[0], 0), (g[1], 0), (g[2], 0), (g[3], 0)
        else: A, A2, B, B2 = (g[0], g[4]), (g[1], g[5]), (g[2], g[6]), (g[3], g[7])
        return abs(corr(wb, wp, X, A, B) - corr(wb, wp, X, A, B2) + corr(wb, wp, X, A2, B) + corr(wb, wp, X, A2, B2))
    best = 0.0
    for seed in range(12):
        g0 = np.random.default_rng(seed).uniform(0, np.pi, 4 if in_plane else 8)
        r = minimize(lambda g: -S(g), g0, method='Nelder-Mead', options=dict(xatol=1e-10, fatol=1e-12, maxiter=20000, maxfev=20000))
        best = max(best, -r.fun)
    return best

print("lattice: %d sites, coin angle a = %.1f, wavenumbers +-%.2f (k/a = %.2f), width %d, packets at %d and %d, %d ticks, range %d"
      % (L, a, k0, k0/a, w, x1, x2, N, sig))
rows = {}
for V0 in (0.0, 0.15, 0.2, 0.3):
    wb, wp, X, tot = scatter(V0); C = 2*abs(X)/(wb + wp)
    Sopt = S_max(wb, wp, X); Spl = S_max(wb, wp, X, in_plane=True)
    rows[V0] = (wb, wp, X, C, Sopt, Spl, tot)
    print("V0 = %.2f : bounce w_b = %.4f  pass w_p = %.4f  |X|/sqrt(w_b w_p) = %.4f  arg X = %+.4f pi  C = %.4f   S over the sphere = %.4f   "
          "2 sqrt(1 + C^2) = %.4f   S in one plane = %.4f   total %.12f"
          % (V0, wb, wp, abs(X)/np.sqrt(wb*wp) if wb*wp > 1e-12 else 0, np.angle(X)/np.pi, C, Sopt, 2*np.sqrt(1 + C*C), Spl, tot))
f = rows[0.0]
ok("1. free walkers pass through each other: bounce %.4f, pass %.4f, S = %.4f: nothing at the detectors to entangle without an interaction"
   % (f[0], f[1], f[4]), f[0] < 1e-3 and abs(f[4] - 2) < 1e-3)
dev = max(abs(rows[V][4] - 2*np.sqrt(1 + rows[V][3]**2)) for V in (0.15, 0.2, 0.3))
ok("2. with the repulsion: S over the sphere = %.3f, %.3f, %.3f at V0 = 0.15, 0.2, 0.3, against 2 sqrt(1 + C^2) = %.3f, %.3f, %.3f (largest "
   "difference %.0e); at V0 = 0.2: bounce %.2f, pass %.2f, |X|/sqrt(w_b w_p) = %.2f - what keeps S below 2 sqrt2 = %.3f"
   % (rows[0.15][4], rows[0.2][4], rows[0.3][4], *[2*np.sqrt(1 + rows[V][3]**2) for V in (0.15, 0.2, 0.3)], dev, rows[0.2][0], rows[0.2][1],
      abs(rows[0.2][2])/np.sqrt(rows[0.2][0]*rows[0.2][1]), 2*np.sqrt(2)), dev < 1e-4 and all(rows[V][4] > 2.1 for V in (0.15, 0.2, 0.3)))
ph_dev = max(abs(np.angle(rows[V][2]) + np.pi/2) for V in (0.15, 0.2, 0.3))
ok("3. arg X = -pi/2 at every V0 (largest deviation %.1e rad): the symmetric barrier; and with all four axes sorted in one plane "
   "S = %.4f, %.4f, %.4f - the real part of an imaginary cross term" % (ph_dev, rows[0.15][5], rows[0.2][5], rows[0.3][5]),
   ph_dev < 5e-3 and all(abs(rows[V][5] - 2) < 2e-3 for V in (0.15, 0.2, 0.3)))
print("ALL PASS" if not fails else "%d FAILED" % fails); sys.exit(1 if fails else 0)
