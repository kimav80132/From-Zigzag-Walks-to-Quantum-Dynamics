"""IV.10: the three-coin automaton against the continuum Dirac column, exactly in momentum space on the same grid.
Lattice step (per p):  U(p) = Umass * Prod_j [P_j^+ e^{-i p_j eps} + P_j^- e^{+i p_j eps}]   (unitary)
Continuum:             e^{-i H(p) t},  H = alpha.p + beta m
Both applied to the same smooth start  g(p) chi0  (Gaussian of width w), then inverse FFT.  t = 3.2/m, m = 1.
Tested (rebuilt 10 Nov 2026: the earlier version printed the table and tested nothing):
 1. the automaton converges onto the continuum column: the gap at N = 64 is smaller than at N = 32, point by point and component by component;
 2. the continuum column is the bulk formula of IV.10 smoothed: formula + (w^2/2) Laplacian(formula), at the six points away from the cone;
 3. at (0, 0, 2.5), 0.7/m inside the cone, the smoothing reaches the skin: the smoothed column has the opposite sign from the bulk formula,
    and the automaton follows the smoothed column.
Needs about 3 GB of memory and about a minute."""
import numpy as np
from scipy.special import j1, jv
m, t = 1.0, 3.2
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
Z2 = np.zeros((2, 2), complex); I2 = np.eye(2)
alpha = [np.block([[s, Z2], [Z2, -s]]) for s in (sx, sy, sz)]
beta = np.block([[Z2, I2], [I2, Z2]])
chi0 = np.array([1, 0, 0, 0], complex)

def run(N, w, pad=6):
    eps = t / N
    L = 2 * N + pad
    k = 2 * np.pi * np.fft.fftfreq(L, d=eps)
    px, py, pz = np.meshgrid(k, k, k, indexing='ij')
    P = np.stack([px, py, pz], -1)                       # (L,L,L,3)
    g = np.exp(-(px**2 + py**2 + pz**2) * w**2 / 2)      # Gaussian envelope, unit total weight
    # continuum: e^{-iHt} = cos(Et) - i sin(Et) H/E
    E = np.sqrt(px**2 + py**2 + pz**2 + m**2)
    H = sum(P[..., j, None, None] * alpha[j] for j in range(3)) + m * beta
    Uc = np.cos(E * t)[..., None, None] * np.eye(4) - 1j * (np.sin(E * t) / E)[..., None, None] * H
    # lattice step
    U = np.cos(m * eps) * np.eye(4) - 1j * np.sin(m * eps) * beta
    U = np.broadcast_to(U, (L, L, L, 4, 4)).copy()
    for j in (2, 1, 0):                                  # applied in order x, y, z (rightmost first)
        Pp = (np.eye(4) + alpha[j]) / 2; Pm = (np.eye(4) - alpha[j]) / 2
        ph = np.exp(-1j * P[..., j] * eps)[..., None, None]
        Sj = Pp * ph + Pm * np.conj(ph)
        U = U @ Sj
    # U^N by repeated squaring on the exponent
    UN = np.broadcast_to(np.eye(4), U.shape).copy(); B = U; n = N
    while n:
        if n & 1: UN = UN @ B
        B = B @ B; n >>= 1
    psi_l = np.fft.ifftn((UN @ chi0) * g[..., None], axes=(0, 1, 2)) / eps**3
    psi_c = np.fft.ifftn((Uc @ chi0) * g[..., None], axes=(0, 1, 2)) / eps**3
    return eps, psi_l, psi_c

def pt(psi, eps, x):
    i = [int(round(v / eps)) % psi.shape[0] for v in x]
    return psi[i[0], i[1], i[2]]

ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
pts = [(0.5, 0, 0), (0, 1.0, 0), (0, 0, 1.5), (1.0, 0.5, 0.5), (-1.25, 0.25, 0.5), (1.5, 1.0, 0.5), (0, 0, 2.5)]
w = 0.35
res = {}
for N in (32, 64):
    e_, l_, c_ = run(N, w); res[N] = (e_, {x: pt(l_, e_, x).copy() for x in pts}, {x: pt(c_, e_, x).copy() for x in pts}); del l_, c_
# the bulk of IV.10 for the start chi0 = (1,0,0,0): K = (d_t - alpha.grad - i m beta) D, D = -(m/4 pi s) J1(ms) inside the cone
fL = lambda x: 1j*m**2*j1(m*np.sqrt(t*t - np.dot(x, x)))/(4*np.pi*np.sqrt(t*t - np.dot(x, x)))                 # component 3
fR = lambda x: m**2*jv(2, m*np.sqrt(t*t - np.dot(x, x)))/(4*np.pi*(t*t - np.dot(x, x)))*(t + x[2])              # component 1
def lap(f, x, h=1e-3):
    x = np.array(x, float); out = -6*f(x)
    for j in range(3):
        d = np.zeros(3); d[j] = h; out += f(x + d) + f(x - d)
    return out/h**2
print(f"smoothing width w = {w}.  Columns: automaton N=32 | N=64 | continuum (same smoothing) | bulk formula | bulk + (w^2/2) Laplacian")
gap = {32: [], 64: []}; rel_s = []; rel_raw = []
for name, comp, f in (("lower block, component 3, i m^2 J1(ms)/(4 pi s)", 2, fL), ("upper block, component 1, m^2 J2(ms)(t+z)/(4 pi s^2)", 0, fR)):
    print(name)
    for x in pts:
        s_ = np.sqrt(t*t - np.dot(x, x)); a32, a64, cc = res[32][1][x][comp], res[64][1][x][comp], res[64][2][x][comp]
        sm = f(x) + w**2/2*lap(f, x)
        print(f"  x={str(x):20s} s={s_:.2f}:  {a32:+.4f} | {a64:+.4f} | {cc:+.4f} | {f(x):+.4f} | {sm:+.4f}")
        gap[32].append(abs(a32 - res[32][2][x][comp])); gap[64].append(abs(a64 - cc))
        if x != (0, 0, 2.5): rel_s.append(abs(cc - sm)/abs(cc)); rel_raw.append(abs(cc - f(x))/abs(cc))
g32, g64 = np.array(gap[32]), np.array(gap[64])
big = g32 > 2e-4                                           # where the N = 32 gap is above the noise of the comparison
ok(f"1. automaton -> continuum: summed over the 14 entries the gap falls from {g32.sum():.4f} (N = 32) to {g64.sum():.4f} (N = 64), a factor {g32.sum()/g64.sum():.1f}; "
   f"entry by entry, where the N = 32 gap exceeds 2e-4, it shrinks by {np.min(g32[big]/g64[big]):.1f} to {np.max(g32[big]/g64[big]):.1f} times",
   g32.sum()/g64.sum() > 2.5 and np.all(g64[big] < g32[big]))
ok(f"2. the continuum column is the bulk formula smoothed: it differs from the bare formula by {100*min(rel_raw):.0f}% to {100*max(rel_raw):.0f}% at the six points "
   f"away from the cone, and from formula + (w^2/2) Laplacian by at most {100*max(rel_s):.1f}%", max(rel_s) < 0.03 and max(rel_raw) > 0.05)
x7 = (0, 0, 2.5)
ok(f"3. at (0,0,2.5), 0.7/m inside the cone: upper block, bulk formula {fR(x7):+.3f}, smoothed continuum {res[64][2][x7][0].real:+.3f}, automaton {res[64][1][x7][0].real:+.3f}: "
   "the smoothing reaches the skin on the cone, the sign is opposite, and the automaton follows the smoothed column",
   fR(x7) > 0.03 and res[64][2][x7][0].real < -0.02 and abs(res[64][1][x7][0] - res[64][2][x7][0]) < 0.2*abs(res[64][2][x7][0]))
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
