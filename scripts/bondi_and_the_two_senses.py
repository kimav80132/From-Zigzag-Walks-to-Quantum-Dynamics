"""Bondi's radar and the walk, at the root. 1+1, c = 1.
1. In the walk's own clock x, an eps = +1 segment keeps u = t - x fixed and an eps = -1 segment keeps v = t + x fixed (III.1: dt = eps dx).
   Bondi's out-leg is a u-line and his return leg a v-line, t1 = u(X), t2 = v(X): the two radar legs ARE the two time senses, read in the clock.
2. The sense swap of III.15a (B) (t -> -t, eps -> -eps) is u <-> -v... precisely: it maps u-lines to v-lines: the operation that makes the return leg
   from the out-leg is the operation that makes Bondi's return signal from his out signal.
3. k-calculus: a boost is (u, v) -> (k u, v/k) on light-cone coordinates and on any vector's light-cone parts; on the readings it is
   Q+ -> Q+/sqrt k, Q- -> sqrt k Q- (III.15c (6)); the content's parts j^u = 2Q-^2, j^v = 2Q+^2 go as (k, 1/k). t1 t2 and Q+ Q- are the invariants. So a record that transforms as
   radar times do is the square of the readings - and no linear map does it, since k != sqrt k (frame_independence_forces_the_square.py 1b).
4. WHERE the rotation reading is (correcting the previous message): a walk read in its own clock is a stochastic step (with loss) whose inverse has
   negative weights - for S in x AND for the telegraph in t. Read on the cut transverse to its clock, with the common loss taken out (III.5a (C)) and the
   exchange multiplied by the direction matrix (III.5a (J)), the step is a rotation - for S in t AND for the telegraph in x. So every walk has a return
   leg on the cut transverse to its clock. III.15b says the physical cut is t. Hence: the return leg is on the physical cut iff the clock is x iff the
   coin is in t (III.6, II.1's sign). The telegraph's return leg exists - on a timelike cut, where no record lives.
5. The noise IS the backward segments: with one sense there is no spread at all (c d_x Q+ = -d_t Q+); with flips, Var(t) along the path = x/(a c)
   (III.5's Kac limit, D = 1/2ac), which is (x/c)^2 divided by the mean number of flips a x/c. Monte Carlo of III.1 against it."""
import numpy as np, sympy as sp
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
# 1-2 symbolic
t, x, dx, eps, k = sp.symbols('t x dx epsilon k', real=True)
u, v = t - x, t + x
du = (t + eps*dx - x - dx) - u; dv = (t + eps*dx + x + dx) - v
ok("1. along a step dt = eps dx: du = (eps - 1) dx, dv = (eps + 1) dx: eps = +1 keeps u, eps = -1 keeps v. Bondi's out-leg (u = t1) is a + segment, his return leg (v = t2) a - segment",
   sp.simplify(du.subs(eps, 1)) == 0 and sp.simplify(dv.subs(eps, -1)) == 0 and sp.simplify(du.subs(eps, -1) + 2*dx) == 0)
ok("2. the sense swap (t -> -t, eps -> -eps) sends u = t - x to -(t + x) = -v: u-lines to v-lines. Return leg = out-leg with the legs' roles exchanged, as in radar",
   sp.simplify(u.subs(t, -t) + v) == 0)
# 3. III.15c (4), from the Lorentz transformation itself, and (6)-(7): which reading goes with which track.  (Rebuilt 10 Nov 2026: the earlier check 3
#    squared sqrt(k) and paired the readings with the tracks the wrong way round; it did not use the boost.)
eta = sp.symbols('eta', real=True); g, b = sp.cosh(eta), sp.tanh(eta); kk = sp.exp(eta)
tp, xp = g*(t - b*x), g*(x - b*t)
r = lambda ex: sp.simplify(sp.expand(ex.rewrite(sp.exp)))
ok("3a. boost t' = gamma(t - beta x), x' = gamma(x - beta t): u' = gamma(1 + beta) u = k u, v' = gamma(1 - beta) v = v/k with k = e^eta = sqrt((1+beta)/(1-beta)); t1' t2' = t1 t2",
   r((tp - xp) - kk*u) == 0 and r((tp + xp) - v/kk) == 0 and r(g*(1 + b) - kk) == 0 and r(g*(1 - b) - 1/kk) == 0
   and r(kk**2 - (1 + b)/(1 - b)) == 0 and r((tp - xp)*(tp + xp) - u*v) == 0)
Qp, Qm, lam, mu = sp.symbols('Q_p Q_m lambda mu', positive=True)
# the lines in null coordinates (III.15c (3)): 2 d_v Q+ = a Q-, 2 d_u Q- = -a Q+.  With d_v = (1/k) d_v', d_u = k d_u', Q+ = lam Q'+, Q- = mu Q'-:
#   2 d_v' Q'+ = a (k mu/lam) Q'-,  2 d_u' Q'- = -a (lam/(k mu)) Q'+.   Same lines iff lam/mu = k; with lam mu = 1: lam = sqrt k.
sol = sp.solve([kk*mu/lam - 1, lam*mu - 1], [lam, mu], dict=True)
ok("3b. the lines keep their form iff lambda/mu = k, and with lambda mu = 1: lambda = sqrt k, mu = 1/sqrt k, i.e. Q'+ = Q+/sqrt k, Q'- = sqrt k Q- (III.15c (6))",
   len(sol) == 1 and r(sol[0][lam] - sp.sqrt(kk)) == 0 and r(sol[0][mu] - 1/sp.sqrt(kk)) == 0)
Qpp, Qmp = Qp/sp.sqrt(kk), sp.sqrt(kk)*Qm
ok("3c. the content's null parts go as the track coordinates they are densities on: j^u = 2Q-^2 -> k (2Q-^2) like u, j^v = 2Q+^2 -> (2Q+^2)/k like v; "
   "Q+ Q- is invariant as t1 t2 is; and the count int Q+^2 du keeps its value, (Q+^2/k)(k du)  (III.15c (7)-(8))",
   r(2*Qmp**2 - kk*2*Qm**2) == 0 and r(2*Qpp**2 - 2*Qp**2/kk) == 0 and r(Qpp*Qmp - Qp*Qm) == 0 and r(Qpp**2*kk - Qp**2) == 0)
# 4. the four readings on the lattice (loss kept where the clock is the process's own)
L, N, a, d = 24, 8, 1.0, 0.25
idx = lambda n_, e: 2*(n_ % L) + (0 if e == +1 else 1); dim = 2*L
S = np.zeros((dim, dim))
for n_ in range(L): S[idx(n_+1, +1), idx(n_, +1)] = 1; S[idx(n_-1, -1), idx(n_, -1)] = 1
T2 = np.array([[0, 1], [-1, 0.]]); sx = np.array([[0, 1], [1, 0.]])
own_clock = np.kron(np.eye(L), (1 - a*d)*np.eye(2) + a*d*sx) @ S           # any walk in its own clock: advance, keep or flip (stochastic)
transverse = np.kron(np.eye(L), np.cos(a*d)*np.eye(2) + np.sin(a*d)*T2) @ S   # the same lines read transverse to the clock, loss out, exchange x direction matrix
ok("4a. in its own clock the step is stochastic and its inverse coin has a negative flip weight: S read in x and the telegraph read in t are this same matrix (x <-> t)",
   np.allclose(own_clock.sum(0), 1) and np.linalg.inv((1 - a*d)*np.eye(2) + a*d*sx)[0, 1] < 0)
ok("4b. transverse to its clock the step is a rotation, U^T U = 1: S read in t and the telegraph read in x are this same matrix (x <-> t); the return leg lives here",
   np.allclose(transverse.T @ transverse, np.eye(dim)))
print("   so: a return leg exists for every walk on the cut transverse to its clock; III.15b puts records on t; S has x as clock (III.6), the telegraph t.")
# 5. the noise is the backward segments
rng = np.random.default_rng(5); c_ = 1.0; X = 40.0; dx_ = 0.01; nsteps = int(X/dx_); npaths = 20000
e = np.ones(npaths); tt = np.zeros(npaths)
for _ in range(nsteps):
    tt += e*dx_/c_
    flip = rng.random(npaths) < a*dx_/c_; e[flip] = -e[flip]
var_mc = tt.var(); var_kac = X/(a*c_)
var_exact = X/(a*c_) - (1 - np.exp(-2*a*X/c_))/(2*a**2)          # the telegraph's exact variance at finite x; the Kac limit is its leading term
se = var_exact*np.sqrt(2/npaths)                                  # sampling error of a variance
ok(f"5. Var(t) along the path at clock x = {X}: Monte Carlo {var_mc:.2f} +- {se:.2f}; exact x/(ac) - (1 - e^(-2ax/c))/(2a^2) = {var_exact:.2f}; "
   f"Kac limit x/(a c) = {var_kac:.2f}; with one sense it is 0. The spread of t is made of the flips",
   abs(var_mc - var_exact) < 3*se and abs(var_exact/var_kac - 1) < 0.02)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
