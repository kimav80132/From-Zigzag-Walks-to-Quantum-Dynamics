"""III.16: the source answer from Part II's walk with tau summed, as a positive Monte Carlo.
  G(s) = int_0^inf dtau e^{-m^2 tau} (4 pi tau)^{-d/2} e^{-s^2/4tau}     (Euclidean signature, D = 1)
       = K_0(ms)/(2 pi) in 1+1,   m K_1(ms)/(4 pi^2 s) in 3+1,
       = the expected time a Brownian particle (<dx_i^2> = 2 dtau) killed at the rate m^2 spends per unit volume at distance s.
1. the two integrals, by quadrature.
2. the Monte Carlo in four dimensions: time spent in each radial shell / shell volume, against the exact density AVERAGED OVER THE SAME SHELL,
   with the sampling error from ten independent batches.
(Rebuilt 10 Nov 2026: the earlier version compared the shell average with the value at the shell's centre and printed a table with no test;
 the density is steep, so that comparison is off by up to 2% in the inner shells for a reason that has nothing to do with the walk.)"""
import numpy as np
from scipy.special import k0, k1
from scipy.integrate import quad
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
m = 1.0
# 1. the integrals
e1 = max(abs(quad(lambda tau: np.exp(-m**2*tau - s**2/(4*tau))/(4*np.pi*tau), 0, np.inf)[0] - k0(m*s)/(2*np.pi)) for s in (0.3, 1.0, 2.5))
e2 = max(abs(quad(lambda tau: np.exp(-m**2*tau - s**2/(4*tau))/(4*np.pi*tau)**2, 0, np.inf)[0] - m*k1(m*s)/(4*np.pi**2*s)) for s in (0.3, 1.0, 2.5))
ok(f"1. int dtau e^(-m^2 tau) (4 pi tau)^(-d/2) e^(-s^2/4tau) = K_0(ms)/(2 pi) for d = 2 (max error {e1:.1e}) and m K_1(ms)/(4 pi^2 s) for d = 4 ({e2:.1e})",
   e1 < 1e-9 and e2 < 1e-9)
# 2. the Monte Carlo
rng = np.random.default_rng(11)
n, nb, dtau = 200_000, 10, 2e-3
edges = np.linspace(0.3, 3.0, 10)
occ = np.zeros((nb, len(edges) - 1))
x = np.zeros((n, 4)); batch = np.arange(n) % nb
life = rng.exponential(1/m**2, size=n)                      # killing time
tau = 0.0; alive = np.ones(n, bool)
while alive.any():
    idx = np.nonzero(alive)[0]
    x[idx] += rng.normal(scale=np.sqrt(2*dtau), size=(idx.size, 4))
    tau += dtau
    s = np.linalg.norm(x[idx], axis=1)
    which = np.digitize(s, edges) - 1; good = (which >= 0) & (which < len(edges) - 1)
    np.add.at(occ, (batch[idx][good], which[good]), dtau)
    alive[idx[life[idx] < tau]] = False
shell_vol = 2*np.pi**2*(edges[1:]**4 - edges[:-1]**4)/4      # 4D: dV = 2 pi^2 s^3 ds
G_b = occ/(n/nb)/shell_vol                                   # one estimate per batch
G_mc = G_b.mean(0); G_se = G_b.std(0, ddof=1)/np.sqrt(nb)
G = lambda s: m*k1(m*s)/(4*np.pi**2*s)
G_shell = np.array([quad(lambda s: G(s)*2*np.pi**2*s**3, lo, hi)[0] for lo, hi in zip(edges[:-1], edges[1:])])/shell_vol
sc = 0.5*(edges[1:] + edges[:-1])
print("   shell          MC density   +-        exact, shell average   (MC - exact)/+-    exact at the centre")
for lo, hi, g_, e_, ex, cen in zip(edges[:-1], edges[1:], G_mc, G_se, G_shell, G(sc)):
    print(f"   {lo:.1f} - {hi:.1f}      {g_:.5f}   {e_:.5f}        {ex:.5f}              {(g_ - ex)/e_:+.1f}              {cen:.5f}")
dev = (G_mc - G_shell)/G_se; rel = np.abs(G_mc/G_shell - 1)
ok(f"2. occupation density of the killed walk = m K_1(ms)/(4 pi^2 s), shell by shell: largest deviation {np.max(np.abs(dev)):.1f} standard errors, "
   f"largest relative difference {100*rel.max():.1f}% (inner seven shells, s < 2.4/m: {100*rel[:7].max():.1f}%)", np.max(np.abs(dev)) < 3.5 and rel[:7].max() < 0.015)
ok(f"3. every contribution is a time spent, >= 0: no signs and no phases in the sum (smallest batch entry {occ.min():.3g})", occ.min() >= 0)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
