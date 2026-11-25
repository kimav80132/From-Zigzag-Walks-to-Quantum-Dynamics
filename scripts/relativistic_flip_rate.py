"""'Rest flip rate + EM flip rate' as 'rest mass + relativistic mass'?
1. free, moving: the turn rate of the column in t is sqrt(a^2 + c^2 k^2) = gamma a: a 'relativistic flip rate', frame-dependent like gamma m.
   The coin itself (sense-flip probability per step) is a dx/c in every frame: the rest value.
2. with Phi: total turn rate (a + e Phi/hbar) at rest. A CONSTANT Phi changes it but nothing observable: the density evolution is identical
   (gauge). So the EM part of the 'flip rate' is like potential energy: gauge-dependent, not a property of the particle alone.
3. the sense-flip (pair) rate is untouched by Phi: in the two-copy step the coin probability is a dx/c with or without Phi."""
import numpy as np, sympy as sp
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
a, c, k, eta = sp.symbols('a c k eta', positive=True)
A = sp.diag(1, -1); T = sp.Matrix([[0, 1], [-1, 0]])
G = -sp.I*c*k*A + a*T
ev = list(G.eigenvals().keys())
ok("1. turn rate in t of a moving free column: eigenvalues +- i sqrt(a^2 + c^2 k^2) = +- i gamma a with k = (a/c) sinh(eta)",
   all(sp.simplify(e_**2 + a**2 + c**2*k**2) == 0 for e_ in ev) and sp.simplify(sp.sqrt(a**2 + c**2*(a/c*sp.sinh(eta))**2) - a*sp.cosh(eta)) == 0)
# 2. constant Phi is unobservable: lattice two-copy walk with Phi = const vs 0
N, Tn, p_, V0 = 200, 60, 0.1, 0.4
xs = np.arange(-N, N + 1); g0 = np.exp(-((xs)/10.0)**2)*np.cos(0.5*xs); g1 = np.exp(-((xs)/10.0)**2)*np.sin(0.5*xs)
def run(Phi):
    psi_p, psi_m = (g0 + 1j*g1).astype(complex), 0.3*(g1 - 1j*g0)
    th = np.arcsin(p_)
    for _ in range(Tn):
        psi_p, psi_m = np.roll(psi_p, 1), np.roll(psi_m, -1)
        psi_p, psi_m = np.cos(th)*psi_p + np.sin(th)*psi_m, -np.sin(th)*psi_p + np.cos(th)*psi_m
        psi_p, psi_m = psi_p*np.exp(-1j*Phi), psi_m*np.exp(-1j*Phi)
    return np.abs(psi_p)**2 + np.abs(psi_m)**2, psi_p
d0, f0 = run(0.0); d1, f1 = run(V0)
ok(f"2. constant Phi: densities identical (max diff {np.abs(d0 - d1).max():.1e}); the column differs only by the overall phase e^{{-i Phi t}} ({np.abs(f1 - f0*np.exp(-1j*V0*Tn)).max():.1e})",
   np.abs(d0 - d1).max() < 1e-12 and np.abs(f1 - f0*np.exp(-1j*V0*Tn)).max() < 1e-12)
# 3 (rebuilt 10 Nov 2026: was printed without a test).  The two-copy step of V.2 (B): R(theta)[(1 - a dx/c) stay + (a dx/c) flipped]; the weight the step
#   moves into a sense from the other one is |R(theta) (a dx/c)(P, P')|, and R is a turn.
th_, dxs, P1, P2 = sp.symbols('theta dx P1 P2', real=True)
Rth = sp.Matrix([[sp.cos(th_), -sp.sin(th_)], [sp.sin(th_), sp.cos(th_)]])
moved = Rth*(a*dxs/c)*sp.Matrix([P1, P2])
ok("3. in the two-copy step the weight flipped into a sense is (a dx/c) times the other sense's length, whatever the turn angle (so whatever Phi and A_x): the field does not touch the coin",
   sp.simplify((moved.T*moved)[0] - (a*dxs/c)**2*(P1**2 + P2**2)) == 0)
# 4. the three readings of the one clock theta = a gamma (t - V x/c^2), and of the field's term (V.4)
tt, xx, V, Phi_, Ax_ = sp.symbols('t x V Phi A_x', real=True); gam = 1/sp.sqrt(1 - V**2/c**2)
clock = a*gam*(tt - V*xx/c**2)
ok("4a. theta = a gamma (t - V x/c^2): at a fixed position it turns at gamma a; along the motion (x = V t) at a/gamma; per proper time (dt = gamma dtau) at a",
   sp.simplify(sp.diff(clock, tt) - gam*a) == 0 and sp.simplify(sp.diff(clock.subs(xx, V*tt), tt) - a/gam) == 0 and sp.simplify(sp.diff(clock.subs(xx, V*tt), tt)*gam - a) == 0)
turn = lambda dt_, dx_: e_*(Phi_*dt_ - Ax_*dx_)                 # the copy turn of a step, times hbar: e(Phi dt - A_x dx)
e_ = sp.symbols('e', positive=True); dt_ = sp.symbols('dt', positive=True)
ok("4b. the field's term the same three ways: e Phi at a fixed position, e (Phi - V A_x) along the motion, e Phi' per proper time with Phi' = gamma (Phi - V A_x), "
   "the potential in the walker's frame",
   sp.simplify(turn(dt_, 0)/dt_ - e_*Phi_) == 0 and sp.simplify(turn(dt_, V*dt_)/dt_ - e_*(Phi_ - V*Ax_)) == 0
   and sp.simplify(turn(dt_, V*dt_)/(dt_/gam) - e_*gam*(Phi_ - V*Ax_)) == 0)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
