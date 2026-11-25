"""IV.8 (C), "the mechanism, as a sum over arrivals": the spin current as a sum over the six arrival directions, with the walk's own weights.
  coin: +1/2 from one side, -1/2 from the other (IV.4c (E));  step: lambda = hbar/2mc (the orbit radius);  speed c;
  quarter-turn about the axis: a walker moving along e_j leaves along s x e_j.
The text's formula:   j_i = -c sum_j [ (1/2) rho(x - lambda e_j) - (1/2) rho(x + lambda e_j) ] eps_ijk s_k
                          = +c sum_j [ (1/2) rho(x - lambda e_j) - (1/2) rho(x + lambda e_j) ] (s x e_j)_i .
Claim: = +(hbar/2m) [grad rho x s]_i = +(hbar/2m) [curl(rho s)]_i for a uniform axis, up to O(lambda^3).
(Rebuilt 10 Nov 2026.  The earlier version used the turn e_j -> e_j x s and so verified MINUS the Gordon current, with the remark "sign = the sense
 of the turn"; the text meanwhile fixes the sense, s x e_j, by the current of (C).  This version uses the text's rule, shows the other sense gives the
 opposite sign, and takes the sign of (C) from the slow limit of the Dirac current itself.)"""
import sympy as sp
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
x, y, z, lam = sp.symbols('x y z lambda', real=True); hbar, m, c = sp.symbols('hbar m c', positive=True)
X = sp.Matrix([x, y, z]); E = [sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 0]), sp.Matrix([0, 0, 1])]
rho = sp.Function('rho')(x, y, z)
s = sp.Matrix(sp.symbols('s1 s2 s3', real=True))          # uniform axis
eps = sp.LeviCivita
shift = lambda j, side: rho.subs({x: x - side*lam*E[j][0], y: y - side*lam*E[j][1], z: z - side*lam*E[j][2]}, simultaneous=True)
def j_text(i):            # the formula as printed
    return -c*sum((shift(j, +1)/2 - shift(j, -1)/2)*sum(eps(i, j, k)*s[k] for k in range(3)) for j in range(3))
def j_turn(i, turn):      # the same sum written with the turned direction: turn(e_j) = s x e_j or e_j x s
    return c*sum((shift(j, +1)/2 - shift(j, -1)/2)*turn(E[j])[i] for j in range(3))
curl_rs = [sum(eps(i, j, k)*sp.diff(rho*s[k], X[j]) for j in range(3) for k in range(3)) for i in range(3)]
gradxs = sp.Matrix([sp.diff(rho, v) for v in X]).cross(s)
lam_val = hbar/(2*m*c)
order = lambda e: sp.series(e, lam, 0, 3).removeO()
ok("1. the printed sum, -c sum_j [rho(x - lambda e_j)/2 - rho(x + lambda e_j)/2] eps_ijk s_k, is the sum over arrivals with the turn e_j -> s x e_j",
   all(sp.expand(j_text(i) - j_turn(i, lambda e: s.cross(e))) == 0 for i in range(3)))
ok("2. it equals +c lambda (grad rho x s)_i = +(hbar/2m)[curl(rho s)]_i with lambda = hbar/2mc; the even orders cancel between the two sides, the next term is O(lambda^3)",
   all(sp.simplify(order(j_text(i)).subs(lam, lam_val) - hbar/(2*m)*curl_rs[i]) == 0 for i in range(3))
   and all(sp.simplify(hbar/(2*m)*curl_rs[i] - hbar/(2*m)*gradxs[i]) == 0 for i in range(3))
   and all(sp.series(j_text(i), lam, 0, 3).removeO().coeff(lam, 2) == 0 and sp.series(j_text(i), lam, 0, 3).removeO().coeff(lam, 0) == 0 for i in range(3)))
ok("3. the other sense of the quarter-turn, e_j -> e_j x s, gives minus this: the sense is not free",
   all(sp.simplify(order(j_turn(i, lambda e: e.cross(s))).subs(lam, lam_val) + hbar/(2*m)*curl_rs[i]) == 0 for i in range(3)))
# 4. which sense the walk has: the slow limit of the Dirac current, uniform axis along z, column at rest (S = 0): j = 2c Re(phi^dag sigma eta), eta = sigma.p phi/2mc
i_ = sp.I; sig = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -i_], [i_, 0]]), sp.diag(1, -1)]
rr = sp.Function('r', positive=True)(x, y, z)            # sqrt(rho)
phi = rr*sp.Matrix([1, 0])                               # axis +z
eta = sum((sig[k]*(-i_*hbar)*sp.diff(phi, X[k]) for k in range(3)), sp.zeros(2, 1))/(2*m*c)
jD = [sp.expand(c*((phi.T*sig[k]*eta)[0] + (phi.T*sig[k]*eta)[0].subs(i_, -i_))) for k in range(3)]      # 2c Re(...): r is real, conjugation is i -> -i
target = (hbar/(2*m)*sp.Matrix([sp.diff(rr**2, v) for v in X]).cross(sp.Matrix([0, 0, 1])))
ok("4. the sign is the walk's: the slow-limit Dirac current of a column at rest with axis z is +(hbar/2m) grad rho x s, the sense of 1.-2.",
   all(sp.simplify(jD[k] - target[k]) == 0 for k in range(3)))
# 5. a sphere of arrival directions instead of six
th, ph = sp.symbols('theta phi', real=True)
n = sp.Matrix([sp.sin(th)*sp.cos(ph), sp.sin(th)*sp.sin(ph), sp.cos(th)])
g = sp.Matrix(sp.symbols('g1 g2 g3', real=True))           # stands for grad rho
avg = lambda f: sp.integrate(sp.integrate(f*sp.sin(th), (th, 0, sp.pi)), (ph, 0, 2*sp.pi))/(4*sp.pi)
sphere = sp.Matrix([sp.simplify(avg(3*n.dot(g)*n.cross(s)[i])) for i in range(3)])
ok("5. int dOmega/4pi 3 (n.grad rho)(n x s) = grad rho x s: the six directions and the sphere agree", sp.simplify(sphere - g.cross(s)) == sp.zeros(3, 1))
a_ = sp.symbols('a', positive=True); D = hbar/(2*m)
ok("6. c lambda = hbar/2m = D, and lambda = c/2a with a = m c^2/hbar: half the mean distance between flips", sp.simplify(c*lam_val - D) == 0
   and sp.simplify(lam_val - (c/(2*a_)).subs(a_, m*c**2/hbar)) == 0)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
