"""
The coin's wall table is one formula, and the "coincidence" with the scalar is central charge c = 1.

Casimir doc, "The coin's version": E(L) = (hbar c / 2 pi L) sum_{n>=1} cos(n phi)/n^2 for a wall factor e^{i phi} per round trip,
with the table  phi = 0: +0.261799,  phi = pi: -0.130900,  phi = pi/2: -0.032725  (L = 1, hbar = c = 1).

Checks (hbar = c = 1):
  1. the round-trip sum equals (pi/2L) B_2(phi/2pi), B_2(x) = x^2 - x + 1/6 (Bernoulli polynomial), at 200 random phi
  2. the same number is  E L / pi = h - 1/24,  h = (1/2)(phi/2pi - 1/2)^2 :  the Casimir energy of a conformal field of
     central charge c = 1 on a circle of length 2L (E = -(2 pi / 2L)(c/24 - h), the weight h of the twisted sector),
     identically in phi; and the three table rows: h = 1/8 (phi = 0), h = 0 (phi = pi), h = 1/32 (phi = pi/2)
  3. the scalar between Dirichlet walls, E = -pi/24L, is the same formula at h = 0: the smooth-cutoff mode sum
     (1/2) sum_n (n pi/L) e^{-eps n pi/L} minus its 1/eps^2 bulk tends to -pi/24L
  4. the coin's turn-matrix wall (phi = pi, half-integer modes) gives the same finite part by its own mode sum:
     -(1/2) sum_{n>=0} 2 (n+1/2)(pi/L) e^{-eps (n+1/2) pi/L} (both signs of energy, the sea's -1/2 Tr|H|) minus its bulk -> -pi/24L
  so the two numbers agree because both fields have c = 1 (a Dirac fermion in 1+1 bosonizes to one scalar), not by accident.
"""
import numpy as np, sys
from mpmath import mp, mpf, nsum, cos, pi, inf, exp, polylog, mpc, re as mre
mp.dps = 30
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all = ok_all and bool(cond); print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()
L = mpf(1)
def E_sum(phi): return (1/(2*pi*L))*mre(polylog(2, exp(mpc(0, 1)*phi)))   # sum_n cos(n phi)/n^2 = Re Li_2(e^{i phi})
def E_B2(phi): x = phi/(2*pi); return (pi/(2*L))*(x**2 - x + mpf(1)/6)
def E_cft(phi): h = (phi/(2*pi) - mpf(1)/2)**2/2; return (pi/L)*(h - mpf(1)/24)
rng = np.random.default_rng(3); phis = [mpf(float(p)) for p in rng.uniform(0, 2*np.pi, 200)]
d1 = max(abs(E_sum(p) - E_B2(p)) for p in phis)
ok("1. sum_n cos(n phi)/n^2 / (2 pi L) = (pi/2L) B_2(phi/2pi) at 200 random phi: max difference %.1e" % float(d1), d1 < mpf(10)**-20)
d2 = max(abs(E_B2(p) - E_cft(p)) for p in phis)
rows = [(0, mpf(1)/8, 0.261799), (pi, mpf(0), -0.130900), (pi/2, mpf(1)/32, -0.032725)]
tab = all(abs(E_cft(p) - mpf(v)) < 1e-6 and abs((p/(2*pi) - mpf(1)/2)**2/2 - h) < 1e-25 for p, h, v in rows)
ok("2. E L / pi = h - 1/24 with h = (1/2)(phi/2pi - 1/2)^2, identically (max difference %.1e); the table: phi = 0 -> h = 1/8, +0.261799; "
   "phi = pi -> h = 0, -0.130900; phi = pi/2 -> h = 1/32, -0.032725: the c = 1 formula E = -(2pi/2L)(c/24 - h) on the unfolded circle of length 2L" % float(d2), d2 < mpf(10)**-25 and tab)
def finite_part(f, eps):   # f(eps) = bulk/eps^2 + finite + O(eps^2): Richardson in eps^2
    a, b = f(eps), f(eps/2); return (4*b - a)/3
a_s = pi/L
scal = lambda eps: (mpf(1)/2)*nsum(lambda n: (n*a_s)*exp(-eps*n*a_s), [1, inf]) - 1/(2*eps**2*a_s)
v3 = finite_part(scal, mpf('0.01'))
ok("3. Dirichlet scalar, smooth cutoff: (1/2) sum n pi/L e^{-eps n pi/L} - bulk -> %.6f = -pi/24L = %.6f (h = 0 row)" % (float(v3), float(-pi/24/L)), abs(v3 + pi/(24*L)) < 1e-6)
coin = lambda eps: -(mpf(1)/2)*nsum(lambda n: 2*((n + mpf(1)/2)*a_s)*exp(-eps*(n + mpf(1)/2)*a_s), [0, inf]) + 1/(eps**2*a_s)
v4 = finite_part(coin, mpf('0.01'))
ok("4. coin, turn-matrix wall (half-integer modes, both signs of energy, -1/2 Tr|H|): finite part %.6f = -pi/24L: the same c = 1 number, bosonization" % float(v4), abs(v4 + pi/(24*L)) < 1e-6)
print("ALL PASS" if ok_all else "SOME FAILED"); raise SystemExit(0 if ok_all else 1)
