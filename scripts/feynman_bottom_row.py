"""Bottom row of the 3+1 table: the Feynman (source) propagator inside the cone from the three-coin automaton.
Feynman's rule: positive-energy part forward in time.  For t>0:
   K_F(t,x) = int d^3p/(2pi)^3 e^{ipx} e^{-iEt} Pi_+(p),   Pi_+ = (1 + H/E)/2,
and K_F = (i d_t + H) Delta_F  with the scalar Feynman function
   Delta_F(t,x) = int d^3p/(2pi)^3 e^{ipx} e^{-iE|t|}/(2E)
   = m K1(m s')/(4 pi^2 s') outside the cone (s'^2 = r^2 - t^2; the killed-Brownian result),
   = C * H1^{(2)}(m s)/s   inside (s^2 = t^2 - r^2), C fixed below from a direct radial integral.
Lower block of K_F from a start chi0 in the upper block is  beta m Delta_F chi0  ->  component 3 = m Delta_F.
The automaton's own positive-energy projector is built from the eigenphases of its unitary step U(p).
"""
import numpy as np
from scipy.integrate import quad
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
# (Rebuilt 10 Nov 2026: the earlier version printed its tables and tested nothing.)
from scipy.special import hankel2, k1
m, t = 1.0, 3.2
# ---- 1. the constant C in Delta_F = C H1^(2)(ms)/s inside the cone, by a direct radial integral
def DeltaF_radial(tt, r, eta=0.02):
    f = lambda p: p * np.sin(p * r) * np.exp(-1j * np.sqrt(p * p + m * m) * abs(tt)) / (2 * np.sqrt(p * p + m * m)) * np.exp(-eta * eta * p * p)
    re = quad(lambda p: f(p).real, 0, 200, limit=4000)[0]; im = quad(lambda p: f(p).imag, 0, 200, limit=4000)[0]
    return (re + 1j * im) / (2 * np.pi**2 * r)
ratios = []
for (tt, r) in [(3.2, 0.8), (3.2, 1.6), (2.4, 0.5)]:
    s = np.sqrt(tt * tt - r * r); D = DeltaF_radial(tt, r); ratios.append(D * s / hankel2(1, m * s))
    print(f"   t={tt} r={r} s={s:.3f}:  Delta_F s / H1^(2)(ms) = {ratios[-1]:+.5f}      (i m/(8 pi) = {1j*m/(8*np.pi):+.5f})")
ok(f"1a. inside the cone Delta_F = (i m/8 pi) H1^(2)(ms)/s: the ratio is the same constant at three unrelated points, within {max(abs(q - 1j*m/(8*np.pi)) for q in ratios)/(m/(8*np.pi)):.1e} of i m/8 pi",
   max(abs(q - 1j*m/(8*np.pi)) for q in ratios) < 2e-3*m/(8*np.pi))
r0 = (DeltaF_radial(0.0, 1.2) / (m * k1(m * 1.2) / (4 * np.pi**2 * 1.2)))
ok(f"1b. at t = 0 (outside the cone) Delta_F is the killed-Brownian m K1(mr)/(4 pi^2 r): ratio {r0.real:.4f} (the regulator e^(-eta^2 p^2), eta = 0.02, accounts for the 0.1%)",
   abs(r0 - 1) < 3e-3)
C = 1j * m / (8 * np.pi)

# ---- 2. automaton with the Feynman projector, vs continuum, vs analytic
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
Z2 = np.zeros((2, 2), complex); I2 = np.eye(2)
alpha = [np.block([[s_, Z2], [Z2, -s_]]) for s_ in (sx, sy, sz)]
beta = np.block([[Z2, I2], [I2, Z2]])
chi0 = np.array([1, 0, 0, 0], complex)

def run(N, w, pad=6, lattice=True):
    eps = t / N; L = 2 * N + pad
    k = 2 * np.pi * np.fft.fftfreq(L, d=eps)
    px, py, pz = np.meshgrid(k, k, k, indexing='ij')
    P = np.stack([px, py, pz], -1)
    g = np.exp(-(px**2 + py**2 + pz**2) * w**2 / 2)
    E = np.sqrt(px**2 + py**2 + pz**2 + m**2)
    H = sum(P[..., j, None, None] * alpha[j] for j in range(3)) + m * beta
    Pi_c = (np.eye(4) + H / E[..., None, None]) / 2
    Uc = np.exp(-1j * E * t)[..., None, None] * Pi_c                      # continuum e^{-iEt} Pi_+
    psi_c = np.fft.ifftn((Uc @ chi0) * g[..., None], axes=(0, 1, 2)) / eps**3
    if not lattice:
        return eps, None, psi_c
    U = np.cos(m * eps) * np.eye(4) - 1j * np.sin(m * eps) * beta
    U = np.broadcast_to(U, (L, L, L, 4, 4)).copy()
    for j in (2, 1, 0):
        Pp = (np.eye(4) + alpha[j]) / 2; Pm = (np.eye(4) - alpha[j]) / 2
        ph = np.exp(-1j * P[..., j] * eps)[..., None, None]
        U = U @ (Pp * ph + Pm * np.conj(ph))
    lam, V = np.linalg.eig(U)                                             # U = V diag(lam) V^-1, |lam|=1
    pos = (np.angle(lam) < 0)                                             # lam = e^{-i E_lat eps}, E_lat > 0
    Vinv = np.linalg.inv(V)
    UN_plus = V @ ((lam**N * pos)[..., :, None] * Vinv)                    # U^N Pi_+^lattice
    psi_l = np.fft.ifftn((UN_plus @ chi0) * g[..., None], axes=(0, 1, 2)) / eps**3
    return eps, psi_l, psi_c

def pt(psi, eps, x):
    i = [int(round(v / eps)) % psi.shape[0] for v in x]
    return psi[i[0], i[1], i[2]]

pts = [(0.5, 0, 0), (0, 1.0, 0), (0, 0, 1.5), (1.0, 0.5, 0.5), (-1.25, 0.25, 0.5), (1.5, 1.0, 0.5)]
w = 0.35
vals = {}
for key, N, ww, lat in (("a32", 32, w, True), ("a64", 64, w, True), ("c15", 64, 0.15, False)):
    e_, l_, c_ = run(N, ww, lattice=lat)
    if lat: vals[key] = {x: pt(l_, e_, x)[2] for x in pts}; vals["c" + key[1:]] = {x: pt(c_, e_, x)[2] for x in pts}
    else: vals[key] = {x: pt(c_, e_, x)[2] for x in pts}
    del l_, c_
print(f"\nLower block, component 3 = m Delta_F.  Target inside the cone: (i m^2/8 pi) H1^(2)(ms)/s")
print("point            s    | automaton+Pi N=32 | N=64        | continuum Pi_+ (w=.35) | continuum (w=.15) | analytic, unsmoothed")
g32, g64, rel15, rel35 = [], [], [], []
for x in pts:
    s = np.sqrt(t * t - np.dot(x, x)); an = m * C * hankel2(1, m * s) / s
    print(f"{str(x):18s}{s:.2f} | {vals['a32'][x]:+.4f} | {vals['a64'][x]:+.4f} | {vals['c64'][x]:+.4f} | {vals['c15'][x]:+.4f} | {an:+.4f}")
    g32.append(abs(vals['a32'][x] - vals['c32'][x])); g64.append(abs(vals['a64'][x] - vals['c64'][x]))
    rel15.append(abs(vals['c15'][x] - an)/abs(an)); rel35.append(abs(vals['c64'][x] - an)/abs(an))
g32, g64 = np.array(g32), np.array(g64)
ok(f"2. the automaton with its own positive-energy projector converges onto the continuum e^(-iEt) Pi_+ column: summed gap {g32.sum():.4f} at N = 32, {g64.sum():.4f} at N = 64 "
   f"(a factor {g32.sum()/g64.sum():.1f}; no entry grows by more than 1e-4)", g32.sum()/g64.sum() > 2 and np.all(g64 < g32 + 1e-4))
ok(f"3. the continuum column lands on (i m^2/8 pi) H1^(2)(ms)/s as the smoothing narrows: within {100*max(rel35):.0f}% at w = 0.35, within {100*max(rel15):.1f}% at w = 0.15, six interior points",
   max(rel15) < 0.03 and max(rel15) < 0.5*max(rel35))
from scipy.special import j1, y1
xx = pts[0]; s = np.sqrt(t*t - np.dot(xx, xx))
ok("4. against the evolution kernel's lower block, (i m^2/4 pi s) J1: Feynman's rule keeps half of the J1 and adds the Y1, H1^(2) = J1 - i Y1",
   abs(m*C*hankel2(1, m*s)/s - (0.5*1j*m**2*j1(m*s)/(4*np.pi*s) + m**2*y1(m*s)/(8*np.pi*s))) < 1e-15)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
