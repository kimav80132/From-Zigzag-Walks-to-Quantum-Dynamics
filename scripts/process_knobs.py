"""Pick your poison: every first-order knob on the one-dimensional walk, and what it gives.
Step (III.2): P_e(t; x+dx) = (1 - a dx/c) P_e(t - e dx/c; x) + (a dx/c) P_-e(t; x).
Knobs (each a first-order change of the step, real weights):
  K1 coin, common:      a -> a + v                      K2 coin, by sense:  a -> a + e v (flip out of + at a+v, out of - at a-v)
  K3 speed, by sense:   t -> t + e (1 + k_e) dx/c        K4 survival, by sense: (1 - a dx/c) -> (1 - a dx/c - l_e dx)
Read in t (loss out), each knob adds a term to d_t Q = -c A d_x Q + a T Q.  Classify by: does it keep the length Q+^2 + Q-^2 ?
does it keep the light cone (coefficient of d_x) ? is it the EM term -i(e Phi/hbar) 1 + i c (e A_x/hbar) A ?
Then: the most general length-keeping first-order real term on two counters is mu T (a coin change). On two copies (four counters),
the length-keeping terms that leave sense and flips alone are exactly span{ 1 (x) eps }: ONE new knob, and it is the charge."""
import sympy as sp
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
t, x, dx = sp.symbols('t x dx', real=True); c, a = sp.symbols('c a', positive=True)
v, kp, km, lp, lm = sp.symbols('v k_p k_m l_p l_m', real=True)
P = {1: sp.Function('Pp')(t, x), -1: sp.Function('Pm')(t, x)}
Q = {1: sp.Function('Qp')(t, x), -1: sp.Function('Qm')(t, x)}
first = lambda ex: sp.expand(sp.series(ex, dx, 0, 2).removeO())
def lines_from_step(step):                                   # step(e) -> P_e(t; x+dx); returns t-reading of Q with loss a taken out
    out = {}
    for e in (1, -1):
        ln = sp.simplify((first(step(e)) - P[e])*c/dx)        # c d_x P_e = ln
        ln = sp.expand(sp.simplify(ln.subs({P[1]: sp.exp(-a*x/c)*Q[1], P[-1]: sp.exp(-a*x/c)*Q[-1]}).doit()*sp.exp(a*x/c)))
        eq = sp.Eq(c*Q[e].diff(x) - a*Q[e], ln)                # c d_x P e^{ax/c} = c d_x Q - a Q
        out[e] = sp.solve(eq, Q[e].diff(t))[0]
    return out
base = lambda e: (1 - a*dx/c)*P[e].subs(t, t - e*dx/c) + a*dx/c*P[-e]
free = lines_from_step(base)
def extra(step):
    L = lines_from_step(step); return {e: sp.simplify(L[e] - free[e]) for e in (1, -1)}
K1 = extra(lambda e: (1 - (a + v)*dx/c)*P[e].subs(t, t - e*dx/c) + (a + v)*dx/c*P[-e])
K2 = extra(lambda e: (1 - (a + e*v)*dx/c)*P[e].subs(t, t - e*dx/c) + (a - e*v)*dx/c*P[-e])
K3 = extra(lambda e: (1 - a*dx/c)*P[e].subs(t, t - e*(1 + (kp if e > 0 else km))*dx/c) + a*dx/c*P[-e])
K4 = extra(lambda e: (1 - a*dx/c - (lp if e > 0 else lm)*dx)*P[e].subs(t, t - e*dx/c) + a*dx/c*P[-e])
def length_rate(ex):                                         # d_t(Q+^2 + Q-^2) from the extra terms alone (streaming part conserves by itself)
    return sp.simplify(2*Q[1]*ex[1] + 2*Q[-1]*ex[-1])
print("K1 coin common   :", K1[1], "|", K1[-1], "| length rate:", length_rate(K1))
print("K2 coin by sense :", K2[1], "|", K2[-1], "| length rate:", length_rate(K2))
print("K3 speed by sense:", K3[1], "|", K3[-1], "| length rate: (d_x terms) ", length_rate(K3))
print("K4 survival      :", K4[1], "|", K4[-1], "| length rate:", length_rate(K4))
ok("K1: a -> a + v adds v T Q - v A Q (in t the clock's common extra loss reads -v A Q); with the loss taken out by Q = e^{-vx/c} Qt the lines are the free ones with a + v, "
   "and the length of Qt is kept: a change of mass",
   sp.simplify(K1[1] - (v*Q[-1] - v*Q[1])) == 0 and sp.simplify(K1[-1] - (-v*Q[1] + v*Q[-1])) == 0
   and all(sp.simplify(((free[e_] + K1[e_]).subs({Q[1]: sp.exp(-v*x/c)*sp.Function('Rp')(t, x), Q[-1]: sp.exp(-v*x/c)*sp.Function('Rm')(t, x)}).doit()*sp.exp(v*x/c)
            - (-e_*c*sp.Function('Rp' if e_ > 0 else 'Rm')(t, x).diff(x) + e_*(a + v)*sp.Function('Rm' if e_ > 0 else 'Rp')(t, x)))) == 0 for e_ in (1, -1)))
ok("K2: sense-skewed coin adds -v (Q+ + Q-) to both lines: real symmetric, length NOT kept (Weyl rescaling)",
   sp.simplify(K2[1] + v*(Q[1] + Q[-1])) == 0 and sp.simplify(K2[-1] + v*(Q[1] + Q[-1])) == 0 and length_rate(K2) != 0)
# K3 (rebuilt 10 Nov 2026: the earlier check tested the + line only, and printed the length rate without testing it)
full3 = {e_: sp.simplify(free[e_] + K3[e_]) for e_ in (1, -1)}
ok("K3: the speed knob divides each line by 1 + k_e: d_t Q_e = e(-c d_x Q_e + a Q_-e)/(1 + k_e). Per unit t, sense e streams at c/(1+k_e) AND is fed at a/(1+k_e)",
   sp.simplify(full3[1] - (-c*Q[1].diff(x) + a*Q[-1])/(1 + kp)) == 0 and sp.simplify(full3[-1] - (c*Q[-1].diff(x) - a*Q[1])/(1 + km)) == 0)
wl = sp.expand((1 + kp)*2*Q[1]*full3[1] + (1 + km)*2*Q[-1]*full3[-1] + c*sp.diff(Q[1]**2 - Q[-1]**2, x))
pl = sp.simplify((2*Q[1]*full3[1] + 2*Q[-1]*full3[-1]).subs({Q[1].diff(x): 0, Q[-1].diff(x): 0}))
ok(f"K3': what it keeps is the weighted length (1+k_+) Q+^2 + (1+k_-) Q-^2 (flux c(Q+^2 - Q-^2)); the plain length changes at the rate {pl} at a uniform point, zero only if k_+ = k_-",
   sp.simplify(wl) == 0 and sp.simplify(pl - 2*a*Q[1]*Q[-1]*(1/(1 + kp) - 1/(1 + km))) == 0 and sp.simplify(pl.subs(km, kp)) == 0)
sfun = sp.Function('s', positive=True)(x); Rr = sp.Function('R')(t, x)
ok("K3'': a speed that varies with x keeps sum Q^2 with the symmetrized form -(c/2){s, d_x}: Q (s d_x Q + d_x(s Q)) = d_x(s Q^2)",
   sp.simplify(Rr*(sfun*Rr.diff(x) + (sfun*Rr).diff(x)) - (sfun*Rr**2).diff(x)) == 0)
ok("K4: survival skew adds -c e l_e Q_e (t-reading): real diagonal; a common l gives the Phi slot, a sense-odd l_e = e l gives the A_x slot; both real (Weyl)",
   sp.simplify(K4[1] + c*lp*Q[1]) == 0 and sp.simplify(K4[-1] - c*lm*Q[-1]) == 0)   # read in t the backward line flips sign: -c l_e e Q_e
# general real first-order term on two counters that keeps the length
M = sp.Matrix(2, 2, sp.symbols('m0:4', real=True)); A = sp.diag(1, -1); T = sp.Matrix([[0, 1], [-1, 0]])
q = sp.Matrix(sp.symbols('q1 q2', real=True))
cond = sp.Poly(sp.expand((q.T*(M + M.T)*q)[0]), *q).coeffs()
solM = sp.solve(cond, list(M), dict=True)[0]; Mk = M.subs(solM)
ok(f"two counters: real first-order terms that keep the length: M = {Mk.tolist()} = mu T only (a coin change). None commutes with A.",
   Mk == Mk[0, 1]*T)
K = lambda X_, Y_: sp.kronecker_product(X_, Y_); eps = sp.Matrix([[0, 1], [-1, 0]])
N = sp.Matrix(4, 4, sp.symbols('n0:16', real=True)); A4, T4 = K(A, sp.eye(2)), K(T, sp.eye(2))
conds = list(N + N.T) + list(N*A4 - A4*N) + list(N*T4 - T4*N)
solN = sp.solve(conds, list(N), dict=True)[0]; Nk = N.subs(solN)
ok(f"two copies: length-keeping terms leaving sense and flips alone: {Nk.tolist()} = j (1 (x) eps): exactly ONE knob, the charge",
   Nk == Nk[0, 1]*K(sp.eye(2), eps))
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
