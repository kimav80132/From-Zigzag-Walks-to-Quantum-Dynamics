"""
Checks of the draft 'Plates made of charged walkers' (4 Nov 2026).
  1. TE/TM matching: solving the 2x2 systems of section 3 gives r_TE = (kz - kmz)/(kz + kmz), r_TM = (eps kz - kmz)/(eps kz + kmz),
     continued to i xi: (kappa - kappa_m)/(kappa + kappa_m), (eps kappa - kappa_m)/(eps kappa + kappa_m)
  2. slab: summing the internal round trips gives r (1 - e^{-2 km d}) / (1 - r^2 e^{-2 km d})
  3. strong response: r_TE -> -1, r_TM -> +1; ideal plates give exactly -pi^2/240
  4. the table of section 6: force / ideal force for omega_p L / c = 0.1 .. 200
  5. section 8: the fractional-sum series for sum 1/(n+x)^2, and B(x) = 2 zeta(2) + 6 zeta(4) x^2 + 10 zeta(6) x^4 + ...
"""
import numpy as np, sympy as sp, mpmath as mp
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all = ok_all and bool(cond); print(("PASS  " if cond else "FAIL  ") + msg)
# 1
kz, kmz, eps = sp.symbols('k_z k_mz epsilon')
rTE, tTE = sp.symbols('r t')
sTE = sp.solve([sp.Eq(1 + rTE, tTE), sp.Eq(kz*(1 - rTE), kmz*tTE)], [rTE, tTE])
sTM = sp.solve([sp.Eq(1 + rTE, tTE), sp.Eq(kz*(1 - rTE), kmz/eps*tTE)], [rTE, tTE])
c1 = sp.simplify(sTE[rTE] - (kz - kmz)/(kz + kmz)) == 0 and sp.simplify(sTM[rTE] - (eps*kz - kmz)/(eps*kz + kmz)) == 0
ok("1. matching: r_TE = %s, r_TM = %s (k_z -> i kappa, k_mz -> i kappa_m leaves both ratios as stated)" % (sp.simplify(sTE[rTE]), sp.simplify(sTM[rTE])), c1)
# 2
r, E = sp.symbols('r E')            # E = e^{-2 km d}; front: r; through front (1+r)... use t t' = 1 - r^2 and internal reflection -r
series = r + (1 - r**2)*(-r)*E/(1 - r**2*E)
ok("2. slab: r + t t' (-r) E sum_n (r^2 E)^n = %s" % sp.factor(sp.simplify(series)), sp.simplify(series - r*(1 - E)/(1 - r**2*E)) == 0)
# 3, 4
xg, wg = np.polynomial.legendre.leggauss(240)
def force_ratio(alpha, X=40.0, nodes=240):
    xg, wg = np.polynomial.legendre.leggauss(nodes)
    x = (xg + 1)/2*X; wx = wg*X/2; mu = (xg + 1)/2; wmu = wg/2
    XX, MM = np.meshgrid(x, mu, indexing='ij'); W = np.outer(wx, wmu)
    km = np.sqrt(XX**2 + alpha**2); epsr = 1 + alpha**2/(XX*MM)**2
    rte = (XX - km)/(XX + km); rtm = (epsr*XX - km)/(epsr*XX + km)
    s = sum(rr**2*np.exp(-2*XX)/(1 - rr**2*np.exp(-2*XX)) for rr in (rte, rtm))
    F = -1/(2*np.pi**2)*(W*XX**3*s).sum()
    return F/(-np.pi**2/240)
# ideal
x = (xg + 1)/2*40; wx = wg*20
Fid = -1/(2*np.pi**2)*(wx*x**3*2*np.exp(-2*x)/(1 - np.exp(-2*x))).sum()
ok("3. ideal plates (r^2 = 1 both polarizations): F L^4/(A hbar c) = %.10f vs -pi^2/240 = %.10f; strong response at alpha = 1e4: ratio %.4f" % (Fid, -np.pi**2/240, force_ratio(1e4)), abs(Fid + np.pi**2/240) < 1e-10 and abs(force_ratio(1e4) - 1) < 1e-3 and force_ratio(1e3) < force_ratio(1e4) < 1)   # 10 Nov 2026: the strong-response ratio was printed, not tested
table = {0.1: 0.0188044, 1: 0.1538280, 5: 0.4655710, 10: 0.6396879, 50: 0.9022738, 200: 0.9739226}
got = {a: force_ratio(a) for a in table}
ok("4. force / ideal force: " + ", ".join("%g: %.7f (draft %.7f)" % (a, got[a], table[a]) for a in table), all(abs(got[a] - table[a]) < 2e-6 for a in table))
got120 = {a: force_ratio(a, nodes=120) for a in table}; d120 = max(abs(got120[a] - got[a]) for a in table)
gotX = {a: force_ratio(a, X=60.0, nodes=360) for a in table}; dX = max(abs(gotX[a] - got[a]) for a in table)
mono = all(got[a] < got[b] for a, b in zip(list(table)[:-1], list(table)[1:]))
ok("4'. the table does not depend on the quadrature: 120 nodes against 240 differ by at most %.0e; cutting x at 60 with 360 nodes instead of 40 with 240, %.0e; "
   "the ratio rises with omega_p L/c at every step" % (d120, dX), d120 < 5e-8 and dX < 5e-8 and mono)   # 10 Nov 2026: the document's '120 and 240 nodes agree to the digits shown' was not in the script
# 5
mp.mp.dps = 30
xv = mp.mpf('0.3')
lhs = mp.nsum(lambda n: 1/(n + xv)**2, [1, mp.inf])
ser = mp.zeta(2) + mp.nsum(lambda m: (-1)**m*(m + 1)*mp.zeta(m + 2)*xv**m, [1, mp.inf])
B = mp.nsum(lambda n: 1/(n + xv)**2 + 1/(n - xv)**2, [1, mp.inf])
Bser = 2*mp.zeta(2) + mp.nsum(lambda k: 2*(2*k + 1)*mp.zeta(2*k + 2)*xv**(2*k), [1, mp.inf])
ok("5. at x = 0.3: sum 1/(n+x)^2 = %s, series zeta(2) - 2zeta(3)x + 3zeta(4)x^2 - ... = %s; B(x) = %s vs 2zeta(2) + 6zeta(4)x^2 + 10zeta(6)x^4 + ... = %s"
   % (mp.nstr(lhs, 15), mp.nstr(ser, 15), mp.nstr(B, 15), mp.nstr(Bser, 15)), abs(lhs - ser) < mp.mpf(10)**-20 and abs(B - Bser) < mp.mpf(10)**-20)
print("ALL PASS" if ok_all else "SOME FAILED"); raise SystemExit(0 if ok_all else 1)
