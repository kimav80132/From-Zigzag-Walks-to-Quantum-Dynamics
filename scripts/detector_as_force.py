"""V.5 (C): detection as a force from the shift of the detector's effective flip rate (the Casimir mechanism for one particle).
Second model of the detector: a contact term V = kappa delta(x - x_D).  1+1, complex column = the two real copies, |psi|^2 = Q^2 + Q'^2.
 1a. first-order shift kappa |psi(x_D)|^2, independent of the phase.
 1b. second-order shift, by exact diagonalization of a lattice Dirac Hamiltonian: it is the bilinear psi^dag(x_D) G psi(x_D), G = sum over the other
     states of phi_n(x_D) phi_n^dag(x_D)/(E_0 - E_n), a 2x2 matrix at the detector's site; unchanged by a gauge turn of the particle.
 2.  the force as the gradient of the shift (first order).
 3.  lattice: a packet passes a pointer.  3b the push and the pull-back, net momentum and displacement; 3c the gauge turn; 3c' the clock; and
 3d. the force law itself, not assumed: the joint state of the particle and a quantum pointer (mass M, its own wave packet) is evolved, and the
     pointer's momentum is compared, tick by tick, with -kappa <d_X rho> from the particle alone.
 4.  the Casimir document's two signs.
(Rebuilt 10 Nov 2026.  Earlier: 1b and 3c were printed without a test, and 3a compared an expression with itself.  FINDING in 1b: for the two-sense
 column the second-order shift is psi^dag G psi with a 2x2 G, a combination of the four invariants of V.5 (B); it is (Q^2 + Q'^2) times a number of
 the detector alone in the limit 1b' tests, a detector whose own level spacing is large against the particle's energies.)"""
import sympy as sp, numpy as np
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
x, xD, g, hb, L, c = sp.symbols('x x_D g hbar L c', real=True)
R, S = sp.Function('R', positive=True)(x), sp.Function('S', real=True)(x)
psi = R*sp.exp(sp.I*S/hb)
dE1 = g*(sp.conjugate(psi)*psi).subs(x, xD)
ok("1a. first-order shift of the detector energy: kappa |psi(x_D)|^2 = kappa R(x_D)^2, independent of the phase S (gauge invariant)",
   sp.simplify(dE1 - g*R.subs(x, xD)**2) == 0)
# ---- 1b. second order, exactly, on a lattice Dirac Hamiltonian with a smooth confining potential (non-degenerate levels)
Ns = 60; rng = np.random.default_rng(5)
sh = np.roll(np.eye(Ns), 1, axis=0); P = (-1j/2)*(sh - sh.T)                      # Hermitian lattice momentum
sz, sx = np.diag([1., -1.]), np.array([[0, 1.], [1, 0]])
U = rng.uniform(-0.6, 0.6, Ns)                                                    # an irregular potential, the same for both senses: no degeneracies
H0 = np.kron(P, sz) + 0.4*np.kron(np.eye(Ns), sx) + np.kron(np.diag(U), np.eye(2))    # index = 2*site + sense
E, V = np.linalg.eigh(H0)
gaps = np.minimum(np.diff(E)[1:], np.diff(E)[:-1]); n0 = 1 + int(np.argmax(gaps))  # the best isolated level
dens = V[0::2, n0]**2 if np.isrealobj(V) else np.abs(V[0::2, n0])**2 + np.abs(V[1::2, n0])**2
dens = np.abs(V[0::2, n0])**2 + np.abs(V[1::2, n0])**2; xd = int(np.argmax(dens))  # the detector where that state is
Pd = np.zeros((2*Ns, 2*Ns)); Pd[2*xd, 2*xd] = Pd[2*xd + 1, 2*xd + 1] = 1.0
col = V[[2*xd, 2*xd + 1], n0]                                                     # the column at the detector's site
G = sum(np.outer(V[[2*xd, 2*xd + 1], n], V[[2*xd, 2*xd + 1], n].conj())/(E[n0] - E[n]) for n in range(2*Ns) if n != n0)
second_pt = (col.conj() @ G @ col).real; first_pt = (col.conj() @ col).real
def level(kap):                                                                   # the same level followed to coupling kap
    Ek, Vk = np.linalg.eigh(H0 + kap*Pd); return Ek[np.argmax(np.abs(Vk.conj().T @ V[:, n0]))]
kk = 1e-3; second_num = (level(kk) + level(-kk) - 2*E[n0])/(2*kk**2); first_num = (level(kk) - level(-kk))/(2*kk)
scalar_guess = first_pt*np.trace(G).real/2                                         # "(Q^2 + Q'^2) times one number"
ph = np.exp(1j*rng.uniform(-np.pi, np.pi, Ns)); Ug = np.kron(np.diag(ph), np.eye(2))
Hg = Ug @ H0 @ Ug.conj().T; Eg, Vg = np.linalg.eigh(Hg); colg = Vg[[2*xd, 2*xd + 1], n0]
Gg = sum(np.outer(Vg[[2*xd, 2*xd + 1], n], Vg[[2*xd, 2*xd + 1], n].conj())/(Eg[n0] - Eg[n]) for n in range(2*Ns) if n != n0)
ok(f"1b. exact levels of H0 + kappa P_D: first order {first_num:.8f} = |psi(x_D)|^2 = {first_pt:.8f}; second order {second_num:.6f} = psi^dag G psi = {second_pt:.6f}, "
   f"a bilinear in the column at x_D; it is NOT |psi(x_D)|^2 times one number (that would be {scalar_guess:.6f}); "
   f"a random gauge turn of the particle leaves it at {(colg.conj() @ Gg @ colg).real:.6f}",
   abs(first_num - first_pt) < 1e-5*first_pt and abs(second_num - second_pt) < 1e-4*abs(second_pt) and abs(scalar_guess - second_pt) > 0.02*abs(second_pt)
   and abs((colg.conj() @ Gg @ colg).real - second_pt) < 1e-9)
# 1b'. the detector with states of its own: two levels 0 and Delta, coupling kappa P_D x (its flip); the joint level by exact diagonalization.
#      When Delta is large against the particle's energies the sum over the particle's states closes, G -> -(1/Delta) 1, and the shift is
#      |psi(x_D)|^2 times a number of the detector alone; when it is not, the shift is the general bilinear.
rows = []
for Delta in (0.5, 5.0, 50.0, 500.0):
    Hj0 = np.kron(H0, np.eye(2)) + np.kron(np.eye(2*Ns), np.diag([0., Delta])); Vj = np.kron(Pd, sx); v0 = np.kron(V[:, n0], [1., 0.])
    def levelj(kap):
        Ek, Vk = np.linalg.eigh(Hj0 + kap*Vj); return Ek[np.argmax(np.abs(Vk.conj().T @ v0))]
    kq = 1e-3*max(1.0, Delta**0.5); num = (levelj(kq) + levelj(-kq) - 2*E[n0])/(2*kq**2)
    GD = sum(np.outer(V[[2*xd, 2*xd + 1], n], V[[2*xd, 2*xd + 1], n].conj())/(E[n0] - E[n] - Delta) for n in range(2*Ns))
    rows.append((Delta, num, (col.conj() @ GD @ col).real, -first_pt/Delta, np.linalg.eigvalsh(-Delta*GD)))
    print(f"      Delta = {Delta:6.1f}: second-order shift {num:+.6e}; psi^dag G psi {rows[-1][2]:+.6e}; -|psi(x_D)|^2/Delta {rows[-1][3]:+.6e}; eigenvalues of -Delta G: {rows[-1][4][0]:.4f}, {rows[-1][4][1]:.4f}")
ok(f"1b'. a detector with a level spacing Delta: the shift is psi^dag G psi at every Delta; for Delta large against the particle's energies (band {E.min():.2f} to {E.max():.2f}) "
   f"G -> -(1/Delta) 1 and the shift is (Q^2 + Q'^2)(x_D) times a number of the detector alone (within {100*abs(rows[2][1]/rows[2][3] - 1):.1f}% at Delta = 50, "
   f"{100*abs(rows[3][1]/rows[3][3] - 1):.2f}% at 500); for Delta = 0.5 it is {rows[0][1]/rows[0][3]:.4f} of that",
   all(abs(r_[1] - r_[2]) < 1e-4*abs(r_[2]) for r_ in rows) and abs(rows[3][1]/rows[3][3] - 1) < 2e-3 and abs(rows[2][1]/rows[2][3] - 1) < 2e-2
   and abs(rows[0][1]/rows[0][3] - 1) > 0.5 and np.allclose(rows[3][4], 1, atol=5e-3))
F = -sp.diff(dE1, xD)
ok("2.  force on the detector = -d(Delta E)/d x_D = -kappa d_x |psi|^2 at x_D (gradient of the Born density), first order",
   sp.simplify(F + g*sp.diff(R**2, x).subs(x, xD)) == 0)
# ---- 3. lattice: Dirac 1+1 packet passes a pointer at xD
N, T, p_, xDn, gn = 300, 220, 0.1, 170, 0.02
xs = np.arange(N); th = np.arcsin(p_)
def packet():
    g0 = np.exp(-((xs - 60)/14.0)**2); a_ = (g0*np.exp(0.6j*xs)).astype(complex); b_ = 0.3*g0*np.exp(0.6j*xs + 0.3j)
    nrm = np.sqrt((np.abs(a_)**2 + np.abs(b_)**2).sum()); return a_/nrm, b_/nrm
def history(chi=None):
    psi_p, psi_m = packet(); rho_hist = []
    lp = lm = 1.0
    if chi is not None:
        psi_p, psi_m = psi_p*np.exp(1j*chi), psi_m*np.exp(1j*chi); lp = np.exp(1j*(chi - np.roll(chi, 1))); lm = np.exp(1j*(chi - np.roll(chi, -1)))
    for _ in range(T):
        psi_p, psi_m = np.roll(psi_p, 1)*lp, np.roll(psi_m, -1)*lm
        psi_p, psi_m = np.cos(th)*psi_p + np.sin(th)*psi_m, -np.sin(th)*psi_p + np.cos(th)*psi_m
        rho_hist.append(np.abs(psi_p)**2 + np.abs(psi_m)**2)
    return np.array(rho_hist)
rho_t = history(); M = 50.0
force = -gn*(rho_t[:, xDn + 1] - rho_t[:, xDn - 1])/2                              # -kappa d_x rho at x_D, tick by tick
pD = np.cumsum(force); xshift = np.sum(pD)/M; phase = gn*rho_t[:, xDn].sum()
ok(f"3b. a packet that passes pushes the pointer along its own direction, down the slope of |psi|^2, then pulls it back ({force[:np.argmax(rho_t[:, xDn])].sum():+.2e} while it arrives, "
   f"{force[np.argmax(rho_t[:, xDn]):].sum():+.2e} while it leaves), net momentum {pD[-1]:+.1e}, but the pointer is DISPLACED by {xshift:+.5f} lattice units",
   force[:np.argmax(rho_t[:, xDn])].sum() > 0 > force[np.argmax(rho_t[:, xDn]):].sum() and abs(pD[-1]) < 1e-3*np.abs(pD).max() and abs(xshift) > 1e-5)
ok(f"3c'. the detector's clock ran fast while the particle was there: extra phase kappa x (time-integrated |psi(x_D)|^2) = {phase:.5f} = kappa x {rho_t[:, xDn].sum():.4f}",
   abs(rho_t[:, xDn].sum() - 1.0) < 0.05)
chi = np.random.default_rng(3).uniform(-np.pi, np.pi, N)
ok(f"3c. a gauge turn of the copies (a random angle at every site, with its link turns) leaves |psi(x_D, t)|^2 at every tick, and so the force, the displacement "
   f"and the clock, unchanged: max difference {np.abs(history(chi) - rho_t).max():.1e}", np.abs(history(chi) - rho_t).max() < 1e-12)
# ---- 3d. the force law from the joint evolution: particle x quantum pointer (coordinate X on the same sites, mass M, a packet of its own)
NX = 64; Xs = xDn - NX//2 + np.arange(NX); sigX, kap = 5.0, 0.02
chiX = np.exp(-((Xs - xDn)/(2*sigX))**2).astype(complex); chiX /= np.sqrt((np.abs(chiX)**2).sum())
pp, pm = packet(); Jp = pp[:, None]*chiX[None, :]; Jm = pm[:, None]*chiX[None, :]              # joint amplitude [x, X]
PX = 2*np.pi*np.fft.fftfreq(NX); kin = np.exp(-1j*PX**2/(2*M))
contact = np.zeros((N, NX), bool); contact[Xs, np.arange(NX)] = True                           # x = X
p_joint, p_law, acc = [], [], 0.0
w = np.abs(chiX)**2
for k in range(T):
    Jp, Jm = np.roll(Jp, 1, axis=0), np.roll(Jm, -1, axis=0)
    Jp, Jm = np.cos(th)*Jp + np.sin(th)*Jm, -np.sin(th)*Jp + np.cos(th)*Jm
    Jp[contact] *= np.exp(-1j*kap); Jm[contact] *= np.exp(-1j*kap)                             # the contact term for one tick
    Jp = np.fft.ifft(np.fft.fft(Jp, axis=1)*kin[None, :], axis=1); Jm = np.fft.ifft(np.fft.fft(Jm, axis=1)*kin[None, :], axis=1)
    Fp, Fm = np.fft.fft(Jp, axis=1), np.fft.fft(Jm, axis=1)
    p_joint.append(((np.abs(Fp)**2 + np.abs(Fm)**2).sum(0)*PX).sum()/NX)                       # <P> of the pointer
    acc += -kap*np.sum(w*(rho_t[k, Xs + 1] - rho_t[k, Xs - 1])/2); p_law.append(acc)           # -kappa <d_X rho>, from the particle alone
p_joint, p_law = np.array(p_joint), np.array(p_law)
ok(f"3d. the pointer's momentum from the joint evolution follows -kappa <d_X |psi|^2> accumulated tick by tick: largest |<P>| {np.abs(p_joint).max():.3e}, "
   f"largest difference from the law {np.abs(p_joint - p_law).max():.1e} ({100*np.abs(p_joint - p_law).max()/np.abs(p_joint).max():.1f}%); net momentum after the passage {p_joint[-1]:+.1e}",
   np.abs(p_joint - p_law).max() < 0.03*np.abs(p_joint).max() and abs(p_joint[-1]) < 0.03*np.abs(p_joint).max())
n = sp.symbols('n', positive=True)
En = n*sp.pi*hb*c/L; Ecas = -sp.pi*hb*c/(24*L)
ok("4.  one particle in mode n: dE_n/dL = -n pi hbar c/L^2 < 0 (pushes the wall out); the loop sum: d/dL(-pi hbar c/24L) = +pi hbar c/(24 L^2) (pulls in)",
   sp.simplify(sp.diff(En, L) + n*sp.pi*hb*c/L**2) == 0 and sp.simplify(sp.diff(Ecas, L) - sp.pi*hb*c/(24*L**2)) == 0)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
