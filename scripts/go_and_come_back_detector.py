"""The direct mechanism: the detector's clock shift is a go-and-come-back amplitude.
Everything real: the particle is two copies (Q, Q') on sites x; the detector bit is two copies at each of its states 0, 1.
Joint state: a real vector with components [x, bit, copy]. Length = Euclidean dot product = what every turn preserves.
Coupling: while the particle is on x_D, turn the detector bit by the angle g (a turn in the (bit0, bit1) plane, same for both copies).
1. The branch handed to the detector is  sin g x (Q, Q')(x_D): LINEAR in the amplitude, phase and all.
2. The detector's clock shift = how much the joint state has turned away from itself = the overlap <before | after>
   = 1 - (1 - cos g) [Q(x_D)^2 + Q'(x_D)^2]: the go-and-come-back pairs each amplitude with itself.
   It is quadratic because 'coming back' means the dot product, and a dot product pairs a component with itself.
3. A local copy turn (a phase of the particle at x_D) changes the handed-over branch, not the overlap.
4. The same on the whole joint vector of a packet (every site, both senses, both copies, both bit states), as a dot product: the overlap is
   1 - (1 - cos g)(Q^2 + Q'^2)(x_D); the clicked branch weighs sin^2 g (Q^2 + Q'^2)(x_D); the overlap goes negative for g > pi/2 when more than half
   of the state sits on the detector's site (it is an overlap, not a probability).
(Check 4 rebuilt 10 Nov 2026.  The earlier one compared the formula with itself and called -d(overlap)/dx_D "the force", pulling the detector towards
 larger density; the overlap is not an energy, and the force of V.5 (C), from the contact term, pushes down the slope.  That line is gone.)"""
import numpy as np, sympy as sp
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
g, Q, Qp, th = sp.symbols('g Q Qp theta', real=True)
# at the site x_D: before = (bit0: Q, Q'; bit1: 0, 0); after the turn of the bit by g
before = sp.Matrix([Q, Qp, 0, 0]); R = sp.Matrix([[sp.cos(g), 0, -sp.sin(g), 0], [0, sp.cos(g), 0, -sp.sin(g)], [sp.sin(g), 0, sp.cos(g), 0], [0, sp.sin(g), 0, sp.cos(g)]])
after = R*before
ok("1. the branch handed to the detector (bit 1) is sin g (Q, Q'): linear in the amplitude", after[2] == sp.sin(g)*Q and after[3] == sp.sin(g)*Qp)
overlap_site = (before.T*after)[0]
ok("2. overlap at x_D = cos g (Q^2 + Q'^2); over the whole state: 1 - (1 - cos g)(Q^2 + Q'^2)(x_D): a go-and-come-back, quadratic",
   sp.simplify(overlap_site - sp.cos(g)*(Q**2 + Qp**2)) == 0)
# 3. local copy turn at x_D: (Q, Q') -> rotated; overlap unchanged, handed branch changed
Rc = sp.Matrix([[sp.cos(th), -sp.sin(th)], [sp.sin(th), sp.cos(th)]]); Qr = Rc*sp.Matrix([Q, Qp])
before_r = sp.Matrix([Qr[0], Qr[1], 0, 0]); after_r = R*before_r
ok("3. a copy turn at x_D (the particle's phase) changes the handed branch, leaves the overlap unchanged",
   sp.simplify((before_r.T*after_r)[0] - overlap_site) == 0 and sp.simplify(after_r[2] - after[2]) != 0)
# 4. the whole joint real vector: index [site, sense, bit, copy]
N, p_ = 200, 0.1
xs = np.arange(N)
def joint(psi_p, psi_m):
    J = np.zeros((N, 2, 2, 2)); J[:, 0, 0, 0], J[:, 0, 0, 1] = psi_p.real, psi_p.imag; J[:, 1, 0, 0], J[:, 1, 0, 1] = psi_m.real, psi_m.imag; return J
def couple(J, xD, gv):                                   # turn the bit by gv on the site xD, the same turn for both senses and both copies
    K = J.copy(); a0, a1 = J[xD, :, 0, :].copy(), J[xD, :, 1, :].copy()
    K[xD, :, 0, :] = np.cos(gv)*a0 - np.sin(gv)*a1; K[xD, :, 1, :] = np.sin(gv)*a0 + np.cos(gv)*a1; return K
g0 = np.exp(-((xs - 90)/10.0)**2); psi_p = (g0*np.exp(0.5j*xs)).astype(complex); psi_m = 0.4*g0*np.exp(0.5j*xs + 0.7j)
nrm = np.sqrt((np.abs(psi_p)**2 + np.abs(psi_m)**2).sum()); psi_p /= nrm; psi_m /= nrm
J0 = joint(psi_p, psi_m); rho = (J0**2).sum(axis=(1, 2, 3))
e_ov = e_cl = 0.0
for xD in (80, 92, 105):
    for gv in (0.4, 1.1, np.pi/2):
        J1 = couple(J0, xD, gv)
        e_ov = max(e_ov, abs(np.sum(J0*J1) - (1 - (1 - np.cos(gv))*rho[xD]))); e_cl = max(e_cl, abs((J1[:, :, 1, :]**2).sum() - np.sin(gv)**2*rho[xD]))
ok(f"4a. on the joint vector of a packet: <before|after> = 1 - (1 - cos g) (Q^2 + Q'^2)(x_D) (max error {e_ov:.1e}); the clicked branch (bit 1) weighs "
   f"sin^2 g (Q^2 + Q'^2)(x_D) ({e_cl:.1e}); three sites, three angles; the length of the joint vector stays 1", e_ov < 1e-14 and e_cl < 1e-14 and abs((couple(J0, 92, 1.1)**2).sum() - 1) < 1e-14)
sharp = np.zeros(N, complex); sharp[92] = np.sqrt(0.8); sharp[93] = np.sqrt(0.2); Js = joint(sharp, 0*sharp)
ov = np.sum(Js*couple(Js, 92, 0.9*np.pi))
ok(f"4b. it is an overlap, not a probability: with 80% of the state on the detector's site and g = 0.9 pi it is {ov:.3f} < 0; for g <= pi/2 it is between 0 and 1",
   ov < 0 and all(0 <= np.sum(Js*couple(Js, 92, gv)) <= 1 for gv in np.linspace(0, np.pi/2, 9)))
# 4c. a copy turn of the particle at every site (the gauge freedom): the handed-over branch changes, overlap and clicked weight do not
chi = np.random.default_rng(1).uniform(-np.pi, np.pi, N); Jc = joint(psi_p*np.exp(1j*chi), psi_m*np.exp(1j*chi))
ok("4c. the particle's copies turned by a different angle at every site: the handed-over branch at x_D changes, the overlap and the clicked weight do not",
   abs(np.sum(Jc*couple(Jc, 92, 1.1)) - np.sum(J0*couple(J0, 92, 1.1))) < 1e-14 and abs((couple(Jc, 92, 1.1)[:, :, 1, :]**2).sum() - (couple(J0, 92, 1.1)[:, :, 1, :]**2).sum()) < 1e-14
   and np.abs(couple(Jc, 92, 1.1)[92, :, 1, :] - couple(J0, 92, 1.1)[92, :, 1, :]).max() > 1e-3)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
