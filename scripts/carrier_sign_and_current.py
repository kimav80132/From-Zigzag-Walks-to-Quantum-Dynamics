"""VI.2 (E), (F): the sign of the carrier, and the conservation of the source.   Rebuilt 9 Nov 2026 against master v34.

1. The carrier carries a one-form index.  Static interaction energy of two like unit sources a distance r apart:
     scalar carrier,   static Lagrangian -1/2 (grad phi)^2 - rho phi:   E = -1/(4 pi r)   like sources ATTRACT
     one-form carrier, static Lagrangian +1/2 (grad a_0)^2 - rho a_0:   E = +1/(4 pi r)   like sources REPEL
   The sign is eta_00 = -1 on the time component (1c).  1a, 1b: the energies from the field on a grid (Gauss's law by FFT, two
   Gaussian sources, self-energies subtracted).  The grid is a periodic box of side L, so what it must give is the periodic-box
   potential, not 1/(4 pi r) bare:
       E_box(r) = erf(r/2w)/(4 pi r) - xi/(4 pi L) + r^2/(6 L^3) + w^2/L^3 + O(r^4/L^5),   xi = 2.837297... (the cubic box's
   Madelung constant), erf(r/2w)/(4 pi r) being two Gaussian sources of width w in open space; the constant is the uniform
   neutralizing background the k = 0 mode removes, the r^2 term its potential.  The earlier
   draft compared with the bare 1/(4 pi r) at L = 32 and failed by exactly that constant, 0.00706.  Checked here at two box sizes.
2. Maxwell needs d_nu j^nu = 0.  The source is the t-reading's content (IV.4d (B)): rho = sum of squares, j = c (squares of the +
   sense minus squares of the - sense).  2a, symbolic, on V.2 (E)'s four real lines with arbitrary Phi(t,x), A_x(t,x):
   d_t rho + d_x j = 0 identically.  (The earlier draft used complex functions and forgot to declare Phi and A real.)
   2b, lattice: V.2's two-leg tick with a potential step keeps the total length, tick by tick."""
import sys
import numpy as np, sympy as sp
fails = 0
def ok(name, cond):
    global fails
    print(("PASS  " if cond else "FAIL  ") + name); fails += (not cond)

# ---- 1. static energies on a periodic grid
XI = 2.8372974794806
def energies(N, L, rs, w):
    h = L/N; x = (np.arange(N) - N//2)*h; X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    k = 2*np.pi*np.fft.fftfreq(N, h); KX, KY, KZ = np.meshgrid(k, k, k, indexing='ij'); K2 = KX**2 + KY**2 + KZ**2; K2[0, 0, 0] = 1.0
    def src(cx):
        s = np.exp(-((X - cx)**2 + Y**2 + Z**2)/(2*w*w)); return s/(s.sum()*h**3)
    def half_grad2(rho):                                   # 1/2 int (grad phi)^2 with -lap phi = rho
        rk = np.fft.fftn(rho); rk[0, 0, 0] = 0
        return 0.5*(np.abs(rk)**2/K2).sum()*h**3/N**3
    out = {}
    for r in rs:
        r1, r2 = src(-r/2), src(+r/2)
        out[r] = half_grad2(r1 + r2) - half_grad2(r1) - half_grad2(r2)      # the cross term: + for the one-form, - for the scalar
    return out
from math import erf
box = lambda r, L, w: erf(r/(2*w))/(4*np.pi*r) - XI/(4*np.pi*L) + r*r/(6*L**3) + w*w/L**3
e32 = energies(128, 32.0, (4.0, 8.0), 0.8); e64 = energies(128, 64.0, (4.0, 8.0), 1.6)
d32 = max(abs(e32[r] - box(r, 32.0, 0.8)) for r in (4.0, 8.0)); d64 = max(abs(e64[r] - box(r, 64.0, 1.6)) for r in (4.0, 8.0))
g = lambda r, w: erf(r/(2*w))/(4*np.pi*r)                                    # two Gaussian sources of width w in open space
diff32 = e32[4.0] - e32[8.0]; diff64 = e64[4.0] - e64[8.0]; bare32 = g(4.0, 0.8) - g(8.0, 0.8); bare64 = g(4.0, 1.6) - g(8.0, 1.6)
ok("1a. one-form carrier (eta_00 = -1): like sources repel. Box of side 32, sources of width 0.8: E(4) = %+.5f, E(8) = %+.5f against the "
   "periodic-box potential %+.5f, %+.5f (open space: %.5f, %.5f; the box's constant -xi/(4 pi L) = %.5f); side 64, width 1.6: %+.5f, %+.5f "
   "against %+.5f, %+.5f; largest difference %.1e and %.1e"
   % (e32[4.0], e32[8.0], box(4.0, 32.0, 0.8), box(8.0, 32.0, 0.8), g(4.0, 0.8), g(8.0, 0.8), -XI/(4*np.pi*32), e64[4.0], e64[8.0],
      box(4.0, 64.0, 1.6), box(8.0, 64.0, 1.6), d32, d64), e32[4.0] > 0 and e32[8.0] > 0 and d32 < 1e-4 and d64 < 1e-4)
ok("1a'. the box's constant drops out of the difference: E(4) - E(8) = %.5f (side 32), %.5f (side 64), against the open-space %.5f, %.5f; "
   "what is left, %.1e and %.1e, is the box's r^2 term (4^2 - 8^2)/(6 L^3) = %.1e, %.1e, eight times smaller in the larger box"
   % (diff32, diff64, bare32, bare64, diff32 - bare32, diff64 - bare64, -48/(6*32.0**3), -48/(6*64.0**3)),
   abs(diff32 - bare32 + 48/(6*32.0**3)) < 1e-4 and abs(diff64 - bare64 + 48/(6*64.0**3)) < 1e-4 and abs(diff64 - bare64) < abs(diff32 - bare32))
ok("1b. scalar carrier: the same field energy with the other sign of the static Lagrangian: E(4) = %+.5f, E(8) = %+.5f - like sources attract"
   % (-e32[4.0], -e32[8.0]), -e32[4.0] < 0 and -e32[8.0] < 0)
# the sign as eta_00: static -1/4 F^2 for a_0 alone is +1/2 (grad a_0)^2; a scalar's -1/2 eta^{ij} d_i phi d_j phi is -1/2 (grad phi)^2
eta = sp.diag(-1, 1, 1, 1); g1, g2, g3 = sp.symbols('g1 g2 g3', real=True); grad = [g1, g2, g3]
F = sp.zeros(4, 4)
for i in range(1, 4): F[0, i] = -grad[i-1]; F[i, 0] = grad[i-1]            # F_{0i} = d_0 a_i - d_i a_0 = -d_i a_0 (static, a_0 only)
Fup = eta*F*eta; L_vec = -sp.Rational(1, 4)*sum(F[m, n]*Fup[m, n] for m in range(4) for n in range(4))
L_sc = -sp.Rational(1, 2)*sum(gi**2 for gi in grad)
ok("1c. static Lagrangian densities: one-form %s, scalar %s: opposite signs, and the one-form's sign is eta_00 = -1 (with eta_00 = +1 it would be %s)"
   % (sp.simplify(L_vec), L_sc, sp.simplify(-sp.Rational(1, 4)*sum(F[m, n]*(sp.eye(4)*F*sp.eye(4))[m, n] for m in range(4) for n in range(4)))),
   sp.simplify(L_vec + L_sc) == 0)

# ---- 2a. conservation with the field on, V.2 (E)'s real lines
t, x = sp.symbols('t x', real=True); e, hb, c, a = sp.symbols('e hbar c a', positive=True)
Phi = sp.Function('Phi', real=True)(t, x); Ax = sp.Function('A', real=True)(t, x)
Q = {s: sp.Function('Q' + n, real=True)(t, x) for s, n in ((1, 'p'), (-1, 'm'))}
Qp = {s: sp.Function('Qprime' + n, real=True)(t, x) for s, n in ((1, 'p'), (-1, 'm'))}
b = lambda s: (e/hb)*(Phi - s*c*Ax)
dQ = {s: -s*c*sp.diff(Q[s], x) + s*a*Q[-s] + b(s)*Qp[s] for s in (1, -1)}        # V.2 (E), unprimed line
dQp = {s: -s*c*sp.diff(Qp[s], x) + s*a*Qp[-s] - b(s)*Q[s] for s in (1, -1)}      # V.2 (E), primed line
rho = sum(Q[s]**2 + Qp[s]**2 for s in (1, -1)); j = c*sum(s*(Q[s]**2 + Qp[s]**2) for s in (1, -1))
dt_rho = sum(2*Q[s]*dQ[s] + 2*Qp[s]*dQp[s] for s in (1, -1))
ok("2a. on V.2 (E)'s four real lines with arbitrary Phi(t,x), A_x(t,x): d_t rho + d_x j = 0 identically, rho = sum of the four squares, "
   "j = c (squares of + minus squares of -); the field terms cancel inside each pair, the flips between the senses",
   sp.simplify(sp.expand(dt_rho + sp.diff(j, x))) == 0)

# ---- 2b. the lattice: V.2's two-leg tick (shift, coin, turn of each pair) with a potential step
L, a_, V0, nt = 400, 0.3, 0.8, 300
xs = np.arange(L); ca, sa = np.cos(a_), np.sin(a_)
th = -V0*(xs > 250)                                                             # theta = -(e/hbar) Phi dt, Phi a step of height V0
ct, st = np.cos(th), np.sin(th)
def tick(q):
    qp, qpp, qm, qmp = np.roll(q[0], 1), np.roll(q[1], 1), np.roll(q[2], -1), np.roll(q[3], -1)
    n0, n2 = ca*qp + sa*qm, -sa*qp + ca*qm; n1, n3 = ca*qpp + sa*qmp, -sa*qpp + ca*qmp
    return np.array([ct*n0 - st*n1, st*n0 + ct*n1, ct*n2 - st*n3, st*n2 + ct*n3])
g = np.exp(-(xs - 150.0)**2/(2*12.0**2)); q = np.array([g*np.cos(0.6*xs), g*np.sin(0.6*xs), 0*g, 0*g]); q /= np.sqrt((q**2).sum())
worst = 0.0
for n in range(nt):
    q = tick(q); worst = max(worst, abs((q**2).sum() - 1))
right = (q[:, xs > 250]**2).sum()
ok("2b. lattice: a packet sent at a potential step of height %.1f, %d ticks: the total length sum(Q^2 + Q'^2) stays 1 to %.0e at every tick, "
   "with %.3f of it past the step and %.3f reflected" % (V0, nt, worst, right, 1 - right), worst < 1e-13 and 0.01 < right < 0.99)
print("      III.3's (q, j) = (P+ - P-, c(P+ + P-)) is the walk's own clock reading's conserved pair; the t-reading's is the squares. "
      "Both are conserved; the source of VI.2 (F) is the second.")
print("ALL PASS" if not fails else "%d FAILED" % fails); sys.exit(1 if fails else 0)
