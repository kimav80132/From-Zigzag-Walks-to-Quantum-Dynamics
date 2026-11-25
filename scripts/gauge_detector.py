"""A detector coupled through the copy turn (gauge coupling). What can it register?
A. Algebra. State of the two-copy walk at a point: (Q+, Q'+, Q-, Q'-). The gauge freedom (V.2 G) turns both copy pairs by the same angle.
   Invariants: no linear form; the quadratic invariants are exactly the four bilinears |psi+|^2, |psi-|^2, Re(psi+* psi-), Im(psi+* psi-).
   So anything coupled through the turn can respond only to lengths and cross-lengths of the column: Q^2-type quantities, never Q.
B. Model. Particle: Dirac 1+1 on a lattice (complex column = two copies). Detector: a bit d in {0,1} sitting at x_D, itself a copy pair.
   Coupling: while the particle is at x_D, the detector bit is turned by the angle g per tick (a turn, i.e. a phase-type coupling).
   Joint state Psi(x, eps, d). The detector's record is the length of its flipped branch, sum |Psi(x, eps, 1)|^2.
   Checks: (1) the record equals sum_t g^2 |psi(x_D, t)|^2 for small g (quadratic response of a turn);
           (2) with the particle gauge-transformed (local phase + link phases) the record is unchanged to machine precision,
               while a linear readout (sum Re Psi) changes;
           (3) for g = pi/2 with absorption (perfect detector) the record and what is still in flight add to 1.  The record (0.935) is close to,
               not equal to, the free packet's weight beyond x_D at the end (0.969): the absorber also takes what the free packet would have
               brought back through x_D.  (Docstring corrected 10 Nov 2026: it said "equals the Born probability of passing x_D".)"""
import numpy as np, sympy as sp
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
# ---- A
th = sp.symbols('theta', real=True); q = sp.Matrix(sp.symbols('Qp Qpp Qm Qmp', real=True))
R = sp.Matrix([[sp.cos(th), -sp.sin(th)], [sp.sin(th), sp.cos(th)]]); G = sp.diag(R, R)
cs = sp.symbols('c0:4', real=True); lin = (sp.Matrix([cs]).T.T*q)[0] if False else sum(ci*qi for ci, qi in zip(cs, q))
res = sp.expand(lin.subs(dict(zip(q, G*q)), simultaneous=True) - lin)
sol_lin = sp.solve([res.subs(th, v) for v in (sp.Rational(1, 3), sp.Rational(1, 2))] and sp.Poly(res.subs(th, sp.Rational(7, 10)), *q).coeffs(), cs, dict=True)
ok(f"A1. linear forms invariant under the copy turn: {sol_lin} (only zero)", sol_lin == [{c: 0 for c in cs}])
ms = sp.symbols('m0:10', real=True)
M = sp.Matrix([[ms[0], ms[1], ms[2], ms[3]], [ms[1], ms[4], ms[5], ms[6]], [ms[2], ms[5], ms[7], ms[8]], [ms[3], ms[6], ms[8], ms[9]]])
Jg = sp.diag(sp.Matrix([[0, -1], [1, 0]]), sp.Matrix([[0, -1], [1, 0]]))        # generator of the common copy turn
cond = Jg.T*M + M*Jg                                                               # infinitesimal invariance, exact and linear in m
Amat, _ = sp.linear_eq_to_matrix([cond[i, j] for i in range(4) for j in range(i, 4)], list(ms)); dim = 10 - Amat.rank()
basis = [q[0]**2 + q[1]**2, q[2]**2 + q[3]**2, q[0]*q[2] + q[1]*q[3], q[0]*q[3] - q[1]*q[2]]   # |psi+|^2, |psi-|^2, Re, Im of psi+* psi-
inv_basis = all(sp.simplify(sp.expand(B.xreplace(dict(zip(q, G*q))) - B)) == 0 for B in basis)
ok(f"A2. quadratic invariants: a {dim}-dimensional space; |psi+|^2, |psi-|^2, Re(psi+* psi-), Im(psi+* psi-) are four independent invariants: {inv_basis}", dim == 4 and inv_basis)
# ---- B
N, T, p_, xD = 240, 150, 0.1, 150
xs = np.arange(N); th_ = np.arcsin(p_)
def init():
    g0 = np.exp(-((xs - 60)/12.0)**2)
    psi_p = (g0*np.exp(0.7j*xs)).astype(complex); psi_m = 0.25*g0*np.exp(0.7j*xs + 0.4j)
    n = np.sqrt((np.abs(psi_p)**2 + np.abs(psi_m)**2).sum()); return psi_p/n, psi_m/n
def run(g, chi=None, absorb=False):
    psi_p, psi_m = init()
    if chi is not None: psi_p, psi_m = psi_p*np.exp(1j*chi), psi_m*np.exp(1j*chi)
    lp = np.exp(1j*(chi - np.roll(chi, 1))) if chi is not None else 1.0        # link phase for a step from x-1 to x
    lm = np.exp(1j*(chi - np.roll(chi, -1))) if chi is not None else 1.0       # from x+1 to x
    Pp = np.stack([psi_p, np.zeros_like(psi_p)], 1); Pm = np.stack([psi_m, np.zeros_like(psi_m)], 1)   # (x, d)
    born_sum = 0.0
    for _ in range(T):
        Pp = np.roll(Pp, 1, 0)*(lp[:, None] if chi is not None else 1); Pm = np.roll(Pm, -1, 0)*(lm[:, None] if chi is not None else 1)
        Pp, Pm = np.cos(th_)*Pp + np.sin(th_)*Pm, -np.sin(th_)*Pp + np.cos(th_)*Pm
        born_sum += (np.abs(Pp[xD, 0])**2 + np.abs(Pm[xD, 0])**2)            # presence of the unflipped branch at x_D
        for P in (Pp, Pm):                                                   # controlled turn of the detector bit at x_D
            a0, a1 = P[xD, 0].copy(), P[xD, 1].copy()
            P[xD, 0], P[xD, 1] = np.cos(g)*a0 - np.sin(g)*a1, np.sin(g)*a0 + np.cos(g)*a1
        if absorb:                                                           # perfect detector: the flipped branch is taken out and counted
            born_sum_abs = (np.abs(Pp[:, 1])**2 + np.abs(Pm[:, 1])**2).sum()
            Pp[:, 1] = 0; Pm[:, 1] = 0; run.absorbed += born_sum_abs
    record = (np.abs(Pp[:, 1])**2 + np.abs(Pm[:, 1])**2).sum()
    linear = (Pp[:, 1].real.sum() + Pm[:, 1].real.sum())
    return record, linear, born_sum, (np.abs(Pp)**2 + np.abs(Pm)**2).sum()
g = 0.05; rec, lin_, bs, tot = run(g)
ok(f"B1. small turn g = {g}: detector record {rec:.6f} vs g^2 x (summed |psi(x_D)|^2 presence) {g**2*bs:.6f}: ratio {rec/(g**2*bs):.4f} (quadratic response; the few % is the coherence between ticks)", abs(rec/(g**2*bs) - 1) < 0.05)
rng = np.random.default_rng(3); chi = rng.uniform(-np.pi, np.pi, N)
rec_g, lin_g, _, tot_g = run(g, chi=chi)
ok(f"B2. gauge-transformed particle: record unchanged ({abs(rec - rec_g):.1e}); a linear readout changes ({lin_:+.4f} -> {lin_g:+.4f}); total length {tot_g:.12f}",
   abs(rec - rec_g) < 1e-12 and abs(lin_ - lin_g) > 1e-3 and abs(tot_g - 1) < 1e-10)
run.absorbed = 0.0; rec_a, _, _, left = run(np.pi/2, absorb=True)
# Born passage probability: fraction of the free packet that crosses x_D = 1 - (what stays left of x_D in a free run)
psi_p, psi_m = init(); Pp, Pm = psi_p.copy(), psi_m.copy()
for _ in range(T):
    Pp, Pm = np.roll(Pp, 1), np.roll(Pm, -1); Pp, Pm = np.cos(th_)*Pp + np.sin(th_)*Pm, -np.sin(th_)*Pp + np.cos(th_)*Pm
free_right = (np.abs(Pp[xs >= xD])**2 + np.abs(Pm[xs >= xD])**2).sum()
ok(f"B3. perfect detector (g = pi/2, absorbing): absorbed {run.absorbed:.4f}, left in flight {left:.4f}, sum {run.absorbed + left:.6f}; "
   f"free packet beyond x_D at the end {free_right:.4f}", abs(run.absorbed + left - 1) < 1e-10)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
