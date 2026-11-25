"""IV.2a: the statements the text marks as checked and cites no script for.  (New, 10 Nov 2026.)
A step of the massless column:  psi(x, t + e) = int dOmega/4pi  W(n) psi(x - e c n, t),  to reproduce d_t psi = -c sigma.grad psi.
1. the averages <n_i n_j> = delta_ij/3 over the sphere; the overlap of two axis columns is cos^2(theta/2).
2. the two conditions int W = 1, int W n_j = sigma_j with W = w0 + lambda sigma.n give w0 = 1, lambda = 3; eigenvalues 4 and -2;
   W = 4 P(n) - 2 P(-n), P(n) = (1 + sigma.n)/2.
3. the inequality: for a positive W (w0 >= |w|) int w.n <= int w0 = 1, but the conditions need 3.  In one dimension the same chain reads 1 <= 1.
4. the holonomy: the product of the projectors round a closed loop of three directions is a positive magnitude times e^{-i Omega/2},
   Omega the solid angle of the spherical triangle (random triangles)."""
import numpy as np, sympy as sp
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
th, ph = sp.symbols('theta phi', real=True)
nv = sp.Matrix([sp.sin(th)*sp.cos(ph), sp.sin(th)*sp.sin(ph), sp.cos(th)])
avg = lambda f: sp.integrate(sp.integrate(f*sp.sin(th), (th, 0, sp.pi)), (ph, 0, 2*sp.pi))/(4*sp.pi)
ok("1a. <n_i n_j> = delta_ij / 3 on the sphere", all(sp.simplify(avg(nv[i]*nv[j]) - (sp.Rational(1, 3) if i == j else 0)) == 0 for i in range(3) for j in range(3)))
sig = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.diag(1, -1)]
chi = lambda t_, p_: sp.Matrix([sp.cos(t_/2), sp.exp(sp.I*p_)*sp.sin(t_/2)])           # the column whose axis is n(theta, phi)
sn = sum((nv[k]*sig[k] for k in range(3)), sp.zeros(2))
res1b = (sn*chi(th, ph) - chi(th, ph)).applyfunc(lambda f: sp.simplify(sp.expand(sp.expand_trig(f.rewrite(sp.exp)))))
ov = (chi(0, 0).H*chi(th, ph))[0]
ok("1b. the column of axis n is the +1 column of sigma.n, and two axis columns overlap by cos^2(theta/2): never zero unless the axes are opposite",
   res1b == sp.zeros(2, 1) and sp.simplify(sp.expand(ov*sp.conjugate(ov)) - sp.cos(th/2)**2) == 0)
w0, lam = sp.symbols('w0 lambda', real=True)
W = w0*sp.eye(2) + lam*sum((nv[k]*sig[k] for k in range(3)), sp.zeros(2))
cond = [sp.simplify(e_) for e_ in list(W.applyfunc(avg) - sp.eye(2))] + [sp.simplify(e_) for j in range(3) for e_ in list((W*nv[j]).applyfunc(avg) - sig[j])]
sol = sp.solve([c_ for c_ in cond if c_ != 0], [w0, lam], dict=True)
Wsol = W.subs(sol[0]); P = lambda sgn: (sp.eye(2) + sgn*sum((nv[k]*sig[k] for k in range(3)), sp.zeros(2)))/2
ok(f"2. int W = 1 and int W n_j = sigma_j with W = w0 + lambda sigma.n give {sol[0]}: W = 1 + 3 sigma.n, eigenvalues 4 and -2, = 4 P(n) - 2 P(-n)",
   sol == [{w0: 1, lam: 3}] and set(sp.simplify(e_) for e_ in Wsol.eigenvals()) == {4, -2} and sp.simplify(Wsol - 4*P(1) + 2*P(-1)) == sp.zeros(2))
# 3. the inequality, numerically on positive weights: W(n) = w0(n) + w(n).sigma with w0 >= |w| drawn at random on a grid of directions, normalized to int W = 1
rng = np.random.default_rng(1); best = 0.0
t_ = np.arccos(1 - 2*(np.arange(400) + 0.5)/400); p_ = np.pi*(1 + 5**0.5)*np.arange(400)                 # 400 directions, equal weights
N3 = np.stack([np.sin(t_)*np.cos(p_), np.sin(t_)*np.sin(p_), np.cos(t_)], 1)
for trial in range(2000):
    w0v = rng.random(400); wv = rng.normal(size=(400, 3)); wv *= (rng.random(400)*w0v/np.linalg.norm(wv, axis=1))[:, None]   # |w| <= w0
    if trial % 2: wv = w0v[:, None]*N3                                                                    # the best case: w along n, |w| = w0
    w0v, wv = w0v/w0v.mean(), wv/w0v.mean()
    best = max(best, np.mean(np.sum(wv*N3, 1)))
ok(f"3. positive weights (w0 >= |w|, int w0 = 1): the largest int w.n reached in 2000 trials is {best:.3f} <= 1; the conditions need 3. In one dimension "
   "the two projectors (1 +- sigma_z)/2 give exactly 1 = 1", best <= 1 + 1e-12 and best > 0.99)
# 4. holonomy on random spherical triangles
def col(n):
    t0, p0 = np.arccos(n[2]), np.arctan2(n[1], n[0]); return np.array([np.cos(t0/2), np.exp(1j*p0)*np.sin(t0/2)])
worst = 0.0
for trial in range(200):
    v = rng.normal(size=(3, 3)); v /= np.linalg.norm(v, axis=1)[:, None]; a_, b_, c_ = v
    Om = 2*np.arctan2(np.dot(a_, np.cross(b_, c_)), 1 + a_@b_ + b_@c_ + c_@a_)             # signed solid angle (Van Oosterom-Strackee)
    ca, cb, cc = col(a_), col(b_), col(c_)
    prod = np.vdot(ca, cc)*np.vdot(cc, cb)*np.vdot(cb, ca)                                 # a -> b -> c -> a, the product of the three overlaps
    worst = max(worst, abs(np.angle(prod*np.exp(1j*Om/2))))
ok(f"4. round a closed loop of three directions the product of the projectors is a positive magnitude times e^(-i Omega/2), Omega the signed solid angle "
   f"enclosed: 200 random triangles, largest phase error {worst:.1e}", worst < 1e-12)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
