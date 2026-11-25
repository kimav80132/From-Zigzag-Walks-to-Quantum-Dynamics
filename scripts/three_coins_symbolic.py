"""IV.10 — the three-coin step, symbolically.
Rules: (1) along each axis the column decides a coin, +-c per step, by the projectors (1 +- alpha_j)/2 [III.1 per axis];
       (2) the same rule on all three axes [IV.1];  (3) the coins are taken in turn [IV.2a];  then III.8's turn matrix for the mass.
Claims checked here with sympy, the limit taken as in I.2 (Taylor in the step, divide by the step):
  A. one step = 1 - (i eps/hbar)(c alpha.p + beta m c^2) + O(eps^2)   [momentum space, 4x4]
  B. the O(eps^2) term depends on the order of the coins only through commutators -> no order is an input
  C. in position space: psi(t+eps) = psi - eps c alpha.grad psi - (i eps m c^2/hbar) beta psi + O(eps^2)  => Dirac (IV.5)
  D. massless, 2x2: the column's axis s = chi^dag sigma chi obeys ds/dt = (2c/hbar) p x s : precession about p-hat,
     angle 2|p| c eps/hbar per step  (the flips ARE the rotation)
"""
import sympy as sp

eps, hbar, c, m = sp.symbols('epsilon hbar c m', positive=True)
px, py, pz = sp.symbols('p_x p_y p_z', real=True)
I2 = sp.eye(2); Z2 = sp.zeros(2)
sx = sp.Matrix([[0, 1], [1, 0]]); sy = sp.Matrix([[0, -sp.I], [sp.I, 0]]); sz = sp.Matrix([[1, 0], [0, -1]])
sig = [sx, sy, sz]; p = [px, py, pz]
alpha = [sp.Matrix(sp.BlockMatrix([[s, Z2], [Z2, -s]])) for s in sig]           # Weyl representation
beta = sp.Matrix(sp.BlockMatrix([[Z2, I2], [I2, Z2]]))
I4 = sp.eye(4)
ok_all = True
def ok(name, cond):
    global ok_all; ok_all &= bool(cond); print(("PASS  " if cond else "FAIL  ") + name)

# ---- A. momentum space: coin_j = P_j^+ e^{-i p_j c eps/hbar} + P_j^- e^{+i p_j c eps/hbar}
def coin(j):
    ph = sp.exp(-sp.I * p[j] * c * eps / hbar)
    return (I4 + alpha[j]) / 2 * ph + (I4 - alpha[j]) / 2 / ph
def mass():
    th = m * c**2 * eps / hbar
    return sp.cos(th) * I4 - sp.I * sp.sin(th) * beta
def step(order):
    U = I4
    for j in order:
        U = coin(j) * U
    return mass() * U
H = c * sum((p[j] * alpha[j] for j in range(3)), sp.zeros(4)) + m * c**2 * beta
U = step((0, 1, 2))
first = (U - I4).applyfunc(lambda e: sp.series(e, eps, 0, 2).removeO()).applyfunc(sp.simplify)
ok("A. one step = 1 - (i eps/hbar) H_Dirac + O(eps^2)", sp.simplify(first - (-sp.I * eps / hbar) * H) == sp.zeros(4))

# ---- B. second order: the order-dependent part is a commutator (vanishes in the limit)
def second(order):
    Uo = step(order)
    return (Uo - I4).applyfunc(lambda e: sp.series(e, eps, 0, 3).removeO().coeff(eps, 2)).applyfunc(sp.simplify)
d = sp.simplify(second((0, 1, 2)) - second((2, 1, 0)))
comm = lambda A, B: A * B - B * A
expected = -(c / hbar)**2 / 2 * (comm(px * alpha[0], py * alpha[1]) + comm(px * alpha[0], pz * alpha[2]) + comm(py * alpha[1], pz * alpha[2])) * 2
# the difference between the two orders is minus twice the half-commutator sum; check it is built only from commutators:
ok("B. the two orders differ at O(eps^2) by a sum of commutators [alpha_j p_j, alpha_k p_k]; same limit",
   sp.simplify(d - (-(c / hbar)**2) * (comm(px * alpha[0], py * alpha[1]) + comm(px * alpha[0], pz * alpha[2]) + comm(py * alpha[1], pz * alpha[2]))) == sp.zeros(4)
   or sp.simplify(d + (-(c / hbar)**2) * (comm(px * alpha[0], py * alpha[1]) + comm(px * alpha[0], pz * alpha[2]) + comm(py * alpha[1], pz * alpha[2]))) == sp.zeros(4))

# ---- C. position space: the coins as shifts, Taylor in the step (I.2's move)
x, y, z, t = sp.symbols('x y z t', real=True)
f = [sp.Function('psi%d' % k)(x, y, z) for k in range(4)]
psi = sp.Matrix(f)
coords = [x, y, z]
def shift(vec, j, a):
    return vec.applyfunc(lambda e: e.subs(coords[j], coords[j] + a))
cur = psi
for j in range(3):
    cur = (I4 + alpha[j]) / 2 * shift(cur, j, -c * eps) + (I4 - alpha[j]) / 2 * shift(cur, j, +c * eps)
th = m * c**2 * eps / hbar
cur = sp.cos(th) * cur - sp.I * sp.sin(th) * beta * cur
rate = ((cur - psi) / eps).applyfunc(lambda e: sp.series(e, eps, 0, 1).removeO()).applyfunc(sp.simplify)
grad = sum((alpha[j] * psi.applyfunc(lambda e: sp.diff(e, coords[j])) for j in range(3)), sp.zeros(4, 1))
dirac_rhs = -c * grad - sp.I * m * c**2 / hbar * beta * psi
ok("C. (psi(t+eps)-psi)/eps -> -c alpha.grad psi - (i m c^2/hbar) beta psi   [i hbar d_t psi = H_Dirac psi]",
   sp.simplify(rate - dirac_rhs) == sp.zeros(4, 1))

# ---- D. massless 2x2: the axis precesses about p
U2 = I2
for j in range(3):
    ph = sp.exp(-sp.I * p[j] * c * eps / hbar)
    U2 = ((I2 + sig[j]) / 2 * ph + (I2 - sig[j]) / 2 / ph) * U2
a1, a2, b1, b2 = sp.symbols('a1 a2 b1 b2', real=True)
chi = sp.Matrix([a1 + sp.I * a2, b1 + sp.I * b2])
sv = [sp.simplify((chi.H * s * chi)[0]) for s in sig]                       # s = chi^dag sigma chi (unnormalized is fine: linear)
chi1 = U2 * chi
sv1 = [sp.simplify(sp.series((chi1.H * s * chi1)[0], eps, 0, 2).removeO()) for s in sig]
pv = sp.Matrix(p); svec = sp.Matrix(sv)
cross = pv.cross(svec)
ok("D. axis after one step = s + eps (2c/hbar) p x s + O(eps^2): precession about p-hat by 2|p|c eps/hbar",
   all(sp.simplify(sv1[k] - sv[k] - eps * 2 * c / hbar * cross[k]) == 0 for k in range(3)))

# ---- E-H (added 10 Nov 2026): the lines of IV.10 the text marks "checked" and cites no script for
import mpmath as mp, numpy as np
from scipy.integrate import quad
from scipy.special import j1 as J1
mp.mp.dps = 30
dK = lambda sv, mv: mp.diff(lambda u: mp.besselk(1, mv*u)/u, sv)
ok("E. d/ds [K_1(ms)/s] = -(m/s) K_2(ms), so d_mu Delta = -(m^2/4 pi^2)(x_mu/s^2) K_2(ms) for Delta = m K_1(ms)/(4 pi^2 s)",
   all(abs(dK(mp.mpf(sv), mp.mpf(mv)) + mp.mpf(mv)/mp.mpf(sv)*mp.besselk(2, mp.mpf(mv)*mp.mpf(sv))) < mp.mpf('1e-18') for sv, mv in (('0.7', '1.3'), ('2.2', '0.6'))))
# F. the factorization with gammas squaring to (+,-,-,-): gamma^0 = beta, gamma^k = beta alpha_k; d = (d_t, d_x, d_y, d_z) as commuting symbols
d0, d1, d2, d3, mm = sp.symbols('d0 d1 d2 d3 m_', real=True)
gam = [beta] + [beta*alpha[k] for k in range(3)]
gd = d0*gam[0] + d1*gam[1] + d2*gam[2] + d3*gam[3]
ok("F. (gamma^0)^2 = 1, (gamma^k)^2 = -1, all anticommute; (gamma.d)^2 = d_t^2 - grad^2 = -Box, and (i gamma.d + m)(i gamma.d - m) = Box - m^2 with Box = -d_t^2 + grad^2",
   gam[0]*gam[0] == I4 and all(gam[k]*gam[k] == -I4 for k in (1, 2, 3)) and all(gam[i_]*gam[j_] + gam[j_]*gam[i_] == sp.zeros(4) for i_ in range(4) for j_ in range(i_))
   and sp.expand(gd*gd - (d0**2 - d1**2 - d2**2 - d3**2)*I4) == sp.zeros(4)
   and sp.expand((sp.I*gd + mm*I4)*(sp.I*gd - mm*I4) - ((-d0**2 + d1**2 + d2**2 + d3**2) - mm**2)*I4) == sp.zeros(4))
ok("G. two different axes composed: P_x^+ P_y^+ = (1 + sigma_x + sigma_y + i sigma_z)/4, the phase appears with the second axis",
   (I2 + sx)/2*(I2 + sy)/2 == (I2 + sx + sy + sp.I*sz)/4)
# H. the bulk of D(t, x) = int d^3p/(2 pi)^3 e^{ip.x} sin(Et)/E inside the cone: -(m/4 pi s) J_1(ms), by the radial integral with a soft cutoff
def D_radial(tt, r, mv=1.0, cut=0.03):
    f = lambda q: q*np.sin(q*r)*np.sin(np.sqrt(q*q + mv*mv)*tt)/np.sqrt(q*q + mv*mv)*np.exp(-cut*cut*q*q)
    return quad(f, 0, 150, limit=6000)[0]/(2*np.pi**2*r)
errH = 0
for tt, r in ((3.2, 0.8), (3.2, 1.6), (2.4, 0.5)):
    s_ = np.sqrt(tt*tt - r*r); errH = max(errH, abs(D_radial(tt, r) + J1(s_)/(4*np.pi*s_))/abs(J1(s_)/(4*np.pi*s_)))
ok(f"H. inside the cone D = -(m/4 pi s) J_1(ms) = (1/2 pi)[-(m/2s) J_1(ms)]: radial integral at three points, largest relative difference {errH:.1e}", errH < 5e-3)
# D', numerically: the axis after one step, to O(eps^2)
def axis_after(ev, pvec, chi_):
    U2n = np.eye(2, dtype=complex); S = [np.array(sig[j].tolist(), dtype=complex) for j in range(3)]
    for j in range(3):
        ph = np.exp(-1j*pvec[j]*ev); U2n = ((np.eye(2) + S[j])/2*ph + (np.eye(2) - S[j])/2/ph) @ U2n
    c1 = U2n @ chi_; return np.array([np.vdot(c1, S[k] @ c1).real for k in range(3)]), np.array([np.vdot(chi_, S[k] @ chi_).real for k in range(3)])
rng = np.random.default_rng(2); pv_ = rng.normal(size=3); chi_ = rng.normal(size=2) + 1j*rng.normal(size=2); chi_ /= np.linalg.norm(chi_)
errs = []
for ev in (1e-2, 5e-3, 2.5e-3):
    s1, s0 = axis_after(ev, pv_, chi_); errs.append(np.linalg.norm(s1 - s0 - ev*2*np.cross(pv_, s0)))
ok(f"D'. the same numerically (hbar = c = 1): |s(t+eps) - s - 2 eps p x s| = {errs[0]:.1e}, {errs[1]:.1e}, {errs[2]:.1e} at eps = 0.01, 0.005, 0.0025: a quarter per halving, O(eps^2)",
   all(0.2 < errs[i_ + 1]/errs[i_] < 0.3 for i_ in range(2)))
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
