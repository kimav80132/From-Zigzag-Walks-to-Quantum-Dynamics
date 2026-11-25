"""VI.3 (A'), the scattered pair in closed form, in the walk's objects: the pair of columns at the detectors M = r chi1 chi2^T - t chi2 chi1^T,
its split by (0) into a symmetric part and a T part, M^T T M = (det M) T, the concurrence C = 2|det M|/|M|^2, the correlation on the sphere
with IV.3's columns chi(theta, phi) = (cos th/2, e^{i phi} sin th/2), and the CHSH value 2 sqrt(1 + C^2) at the stated angles."""
import sympy as sp, numpy as np
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
rr, tt, d, th12 = sp.symbols('r t delta theta12', positive=True); r = rr; t = tt*sp.exp(sp.I*d)
T = sp.Matrix([[0, 1], [-1, 0]]); chi1 = sp.Matrix([1, 0]); chi2 = sp.Matrix([sp.cos(th12/2), sp.sin(th12/2)])
M = r*chi1*chi2.T - t*chi2*chi1.T
Ssym = (chi1*chi2.T + chi2*chi1.T)/2; A = (chi1*chi2.T - chi2*chi1.T)/2
ok("9.  M = (r - t) S + (r + t) A, with A = (chi1^T T chi2) T/2 by (0): the pair's symmetric part and its T part",
   sp.simplify(M - ((r - t)*Ssym + (r + t)*A)) == sp.zeros(2) and sp.simplify(A - (chi1.T*T*chi2)[0]*T/2) == sp.zeros(2))
detM = sp.simplify(M.det())
ok("10. M^T T M = (det M) T, det M = r t e^{i delta} sin^2(theta12/2): the pair turned by T between its columns and read against itself is T again, times det M",
   sp.simplify(M.T*T*M - detM*T) == sp.zeros(2) and sp.simplify(detM - rr*tt*sp.exp(sp.I*d)*sp.sin(th12/2)**2) == 0)
normM = sp.simplify(sp.expand(sum(sp.expand_complex(M[i, j]*sp.conjugate(M[i, j])) for i in range(2) for j in range(2))))
absdet = rr*tt*sp.sin(th12/2)**2                                   # |det M|, from 10
Cexpr = 2*rr*tt*sp.sin(th12/2)**2/(rr**2 + tt**2 - 2*rr*tt*sp.cos(d)*sp.cos(th12/2)**2)
ok("11. |M|^2 = |r|^2 + |t|^2 - 2|rt| cos delta cos^2(th12/2), so C = 2|det M|/|M|^2 = 2|rt| sin^2(th12/2) / (|r|^2 + |t|^2 - 2|rt| cos delta cos^2(th12/2))",
   sp.simplify(normM - (rr**2 + tt**2 - 2*rr*tt*sp.cos(d)*sp.cos(th12/2)**2)) == 0 and sp.simplify(2*absdet/normM - Cexpr) == 0)
ok("12. opposite axes: C = 2|rt|/(|r|^2+|t|^2); symmetric barrier (delta = pi/2): C = 2|rt| sin^2(th12/2)/(|r|^2+|t|^2)",
   sp.simplify(Cexpr.subs(th12, sp.pi) - 2*rr*tt/(rr**2 + tt**2)) == 0 and sp.simplify(Cexpr.subs(d, sp.pi/2) - 2*rr*tt*sp.sin(th12/2)**2/(rr**2 + tt**2)) == 0)
# 13. the correlation on the sphere, opposite axes (numerical, 300 random cases)
rng = np.random.default_rng(0)
chi = lambda th, ph: np.array([np.cos(th/2), np.exp(1j*ph)*np.sin(th/2)]); P = lambda th, ph: np.outer(chi(th, ph), chi(th, ph).conj())
def E_of(Mn, a, b):
    Pa, Qa = P(*a), np.eye(2) - P(*a); Pb, Qb = P(*b), np.eye(2) - P(*b); W = lambda X, Y: np.linalg.norm(X @ Mn @ Y.T)**2
    return (W(Pa, Pb) + W(Qa, Qb) - W(Pa, Qb) - W(Qa, Pb))/(W(Pa, Pb) + W(Qa, Qb) + W(Pa, Qb) + W(Qa, Pb))
err = 0
for _ in range(300):
    rv, tv, dv = rng.uniform(0.2, 1), rng.uniform(0.2, 1), rng.uniform(0, 2*np.pi); Cv = 2*rv*tv/(rv*rv + tv*tv)
    Mn = rv*np.outer([1, 0], [0, 1]) - tv*np.exp(1j*dv)*np.outer([0, 1], [1, 0]); a = rng.uniform(0, np.pi, 2); b = rng.uniform(0, np.pi, 2)
    err = max(err, abs(E_of(Mn, a, b) - (-np.cos(a[0])*np.cos(b[0]) - Cv*np.sin(a[0])*np.sin(b[0])*np.cos(a[1] - b[1] - dv))))
ok("13. E(a,b) = -cos th_a cos th_b - C sin th_a sin th_b cos(phi_a - phi_b - delta): max error %.1e over 300 random (r, t, delta, axes)" % err, err < 1e-12)
# 14. the CHSH value at the stated angles
Cs = sp.symbols('C', positive=True); thb = sp.atan(Cs)
E = lambda ta, tb: -sp.cos(ta)*sp.cos(tb) - Cs*sp.sin(ta)*sp.sin(tb)
S = E(0, thb) - E(0, sp.pi - thb) + E(sp.pi/2, thb) + E(sp.pi/2, sp.pi - thb)
ok("14. with phi_a - phi_b = delta, th_a = 0, th_a' = pi/2, tan th_b = C, th_b' = pi - th_b:  |S| = 2 sqrt(1 + C^2)", sp.simplify(sp.Abs(S) - 2*sp.sqrt(1 + Cs**2)) == 0)
# 14'. (added 10 Nov 2026) the maximum over all axes, by a search, for three values of C: it does not exceed 2 sqrt(1 + C^2) and reaches it
from scipy.optimize import minimize
def Ssearch(Cv, dv=0.6, tries=40):
    Mn = 1.0*np.outer([1, 0], [0, 1]) - ((1 - np.sqrt(1 - Cv**2))/Cv)*np.exp(1j*dv)*np.outer([0, 1], [1, 0])     # |t|/|r| chosen so that 2|rt|/(r^2+t^2) = C
    f = lambda v: -abs(E_of(Mn, v[0:2], v[4:6]) - E_of(Mn, v[0:2], v[6:8]) + E_of(Mn, v[2:4], v[4:6]) + E_of(Mn, v[2:4], v[6:8]))
    return max(-minimize(f, rng.uniform(0, 2*np.pi, 8), method='Nelder-Mead', options={'xatol': 1e-9, 'fatol': 1e-12, 'maxiter': 20000}).fun for _ in range(tries))
found = {Cv: Ssearch(Cv) for Cv in (0.3, 0.7, 1.0)}
ok("14'. a search over all eight detector angles finds max |S| = " + ", ".join(f"{found[Cv]:.4f}" for Cv in found) + " at C = 0.3, 0.7, 1.0, against 2 sqrt(1 + C^2) = "
   + ", ".join(f"{2*np.sqrt(1 + Cv**2):.4f}" for Cv in found) + ": reached, not exceeded (Horodecki's bound)",
   all(abs(found[Cv] - 2*np.sqrt(1 + Cv**2)) < 1e-4 for Cv in found))
print("    (on the lattice, bell_from_scattering.py's optimizer finds 2.194, 2.489, 2.150 at V0 = 0.15, 0.2, 0.3.)")
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
