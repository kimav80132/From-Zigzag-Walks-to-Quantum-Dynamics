"""VI.3 (B): the two-walker event weight at two detectors, computed in the walk's own objects -- III.8a's half-angle columns chi(n) (real, in the
plane of the detector axes), IV.2a's sorting P(n) = chi(n) chi(n)^T, VI.1's pairings D and E, III.15a's square -- for two slow walkers sharing
one packet with axes +x and -x.  (4) the overlaps; (5) W(+a,+b) = (1/2)[1 - cos(a-b)] rho_L rho_R, the exchange term being -sin a sin b;
(6) E(a,b) = -cos(a-b); (7) |S| = 2 sqrt 2; (8) without the exchange term E = -cos a cos b and |S| <= 2."""
import sympy as sp
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
a, b = sp.symbols('a b', real=True)
chi = lambda th: sp.Matrix([sp.cos(th/2), sp.sin(th/2)])       # III.8a: the column whose axis is at angle th from x-hat, half the angle
P = lambda th: chi(th)*chi(th).T                                 # IV.2a: sort by the axis n: P(n) = (1 + M(n))/2
cp, cm = chi(0), chi(sp.pi)                                      # axes +x and -x: (1,0) and (0,1)
ok("4a. |P_a chi(+x)|^2 = cos^2(a/2) = (1 + cos a)/2 (IV.7's overlap), |P_a chi(-x)|^2 = sin^2(a/2)",
   sp.simplify((cp.T*P(a)*cp)[0] - sp.cos(a/2)**2) == 0 and sp.simplify((cm.T*P(a)*cm)[0] - sp.sin(a/2)**2) == 0)
ok("4b. chi(+x)^T P_a chi(-x) = cos(a/2) sin(a/2) = sin(a)/2: the cross overlap, real", sp.simplify((cp.T*P(a)*cm)[0] - sp.sin(a)/2) == 0)
# (0) and (3), added 10 Nov 2026: the two lines (5) rests on, which the earlier version used without testing
th1, th2 = sp.symbols('theta1 theta2', real=True); Tm = sp.Matrix([[0, 1], [-1, 0]])
c1, c2 = chi(th1), chi(th2)
ok("0. chi1 chi2^T - chi2 chi1^T = (chi1^T T chi2) T for any two axis columns, and chi1^T T chi2 = sin((theta2 - theta1)/2): the odd pair's axis part is the turn T",
   sp.simplify(c1*c2.T - c2*c1.T - (c1.T*Tm*c2)[0]*Tm) == sp.zeros(2) and sp.simplify((c1.T*Tm*c2)[0] - sp.sin((th2 - th1)/2)) == 0)
frob = lambda Mx: sum(Mx[i, j]**2 for i in range(2) for j in range(2))
Dm = (P(a)*cp)*(P(b)*cm).T; Em = (P(a)*cm)*(P(b)*cp).T
ok("3. |D - E|^2 as a sum of squares of entries = |P_a phi1|^2 |P_b phi2|^2 + |P_a phi2|^2 |P_b phi1|^2 - 2 (P_a phi1 . P_a phi2)(P_b phi2 . P_b phi1)",
   sp.simplify(frob(Dm - Em) - ((cp.T*P(a)*cp)[0]*(cm.T*P(b)*cm)[0] + (cm.T*P(a)*cm)[0]*(cp.T*P(b)*cp)[0] - 2*(cp.T*P(a)*cm)[0]*(cm.T*P(b)*cp)[0])) == 0)
def W(PA, PB, exchange=True):
    A11 = (cp.T*PA*cp)[0]; A22 = (cm.T*PA*cm)[0]; A12 = (cp.T*PA*cm)[0]
    B11 = (cp.T*PB*cp)[0]; B22 = (cm.T*PB*cm)[0]; B21 = (cm.T*PB*cp)[0]
    return sp.simplify(A11*B22 + A22*B11 - (2*A12*B21 if exchange else 0))
Pa, Qa, Pb, Qb = P(a), sp.eye(2) - P(a), P(b), sp.eye(2) - P(b)
Wpp, Wpm, Wmp, Wmm = W(Pa, Pb), W(Pa, Qb), W(Qa, Pb), W(Qa, Qb)
ok("5. W(+a,+b) = (1 - cos(a-b))/2, W(+a,-b) = (1 + cos(a-b))/2, the four add to 2 (times rho_L rho_R)",
   sp.simplify(Wpp - (1 - sp.cos(a - b))/2) == 0 and sp.simplify(Wpm - (1 + sp.cos(a - b))/2) == 0 and sp.simplify(Wpp + Wpm + Wmp + Wmm - 2) == 0)
E = sp.simplify((Wpp + Wmm - Wpm - Wmp)/(Wpp + Wpm + Wmp + Wmm))
ok("6. E(a,b) = -cos(a-b)", sp.simplify(E + sp.cos(a - b)) == 0)
S = E.subs({a: 0, b: sp.pi/4}) - E.subs({a: 0, b: 3*sp.pi/4}) + E.subs({a: sp.pi/2, b: sp.pi/4}) + E.subs({a: sp.pi/2, b: 3*sp.pi/4})
ok("7. |S| at (0, pi/2; pi/4, 3pi/4) = 2 sqrt 2 = %.4f" % abs(float(S)), sp.simplify(abs(S) - 2*sp.sqrt(2)) == 0)
Wd = [W(PA, PB, exchange=False) for PA, PB in ((Pa, Pb), (Pa, Qb), (Qa, Pb), (Qa, Qb))]
Ed = sp.simplify((Wd[0] + Wd[3] - Wd[1] - Wd[2])/sum(Wd))
Sd = lambda a1, a2, b1, b2: Ed.subs({a: a1, b: b1}) - Ed.subs({a: a1, b: b2}) + Ed.subs({a: a2, b: b1}) + Ed.subs({a: a2, b: b2})
import numpy as np
rng = np.random.default_rng(0); Smax = max(abs(float(Sd(*rng.uniform(0, 2*np.pi, 4)))) for _ in range(2000))
ok("8. without the exchange term: E = -cos a cos b, and |S| <= 2 (max over 2000 random angle sets: %.4f)" % Smax, sp.simplify(Ed + sp.cos(a)*sp.cos(b)) == 0 and Smax <= 2 + 1e-9)
# 8'. the even combination D + E (added 10 Nov 2026): the text's parenthesis in (8)
def Weven(PA, PB):
    A11 = (cp.T*PA*cp)[0]; A22 = (cm.T*PA*cm)[0]; A12 = (cp.T*PA*cm)[0]; B11 = (cp.T*PB*cp)[0]; B22 = (cm.T*PB*cm)[0]; B21 = (cm.T*PB*cp)[0]
    return sp.simplify(A11*B22 + A22*B11 + 2*A12*B21)
We = [Weven(PA, PB) for PA, PB in ((Pa, Pb), (Pa, Qb), (Qa, Pb), (Qa, Qb))]
Ee = sp.simplify((We[0] + We[3] - We[1] - We[2])/sum(We))
Se = Ee.subs({a: 0, b: -sp.pi/4}) - Ee.subs({a: 0, b: -3*sp.pi/4}) + Ee.subs({a: sp.pi/2, b: -sp.pi/4}) + Ee.subs({a: sp.pi/2, b: -3*sp.pi/4})
ok("8'. the even combination D + E has the exchange term with the other sign, E = -cos a cos b + sin a sin b = -cos(a + b), and the same |S| = 2 sqrt 2 (at b -> -b): "
   "the exchange term is the entanglement, its sign is exclusion", sp.simplify(Ee + sp.cos(a + b)) == 0 and sp.simplify(abs(Se) - 2*sp.sqrt(2)) == 0)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
