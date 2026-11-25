"""III.15, "the same statement from the lines": every numbered line (1)-(8) of the text, checked.  c = 1.
sigma = a^2 (x^2 - t^2)/4.  The lines with the loss out (III.5a (F)):  (d_x + eps d_t) Q_eps = a Q_-eps.
The check labels are the text's equation numbers.  (Rebuilt 10 Nov 2026: the earlier version checked (2), (4), (5), (6) only.)"""
import sympy as sp, mpmath as mp
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
t, x = sp.symbols('t x', real=True); a = sp.symbols('a', positive=True); sig = sp.symbols('sigma', real=True)
S = a**2*(x**2 - t**2)/4

# (1) the lines applied twice, on generic functions
Qp, Qm = sp.Function('Qp')(t, x), sp.Function('Qm')(t, x)
Dp = lambda f: sp.diff(f, x) + sp.diff(f, t); Dm = lambda f: sp.diff(f, x) - sp.diff(f, t)
# lines: Dp Qp = a Qm, Dm Qm = a Qp.  Dm Dp Qp = a Dm Qm = a^2 Qp, and Dm Dp = d_x^2 - d_t^2
ok("(1) (d_x - d_t)(d_x + d_t) = d_x^2 - d_t^2 on any function, so the two lines give (d_x^2 - d_t^2) Q = a^2 Q for each sense",
   sp.simplify(Dm(Dp(Qp)) - (sp.diff(Qp, x, 2) - sp.diff(Qp, t, 2))) == 0
   and sp.simplify(Dm(a*Qm).subs(sp.diff(Qm, x), sp.diff(Qm, t) + a*Qp) - a**2*Qp) == 0
   and sp.simplify(Dp(a*Qp).subs(sp.diff(Qp, x), -sp.diff(Qp, t) + a*Qm) - a**2*Qm) == 0)
# ... and the same two lines solved for d_t (III.5a (H)): d_t Q+ = -d_x Q+ + a Q-, d_t Q- = d_x Q- - a Q+
ok("(1') solved for d_t the lines are the same two equations: d_t Q+ + d_x Q+ - a Q- and -(d_t Q- - d_x Q- + a Q+) are the two lines again",
   sp.simplify((sp.diff(Qp, t) + sp.diff(Qp, x) - a*Qm) - (Dp(Qp) - a*Qm)) == 0
   and sp.simplify(-(sp.diff(Qm, t) - sp.diff(Qm, x) + a*Qp) - (Dm(Qm) - a*Qp)) == 0)

# (2) the wave operator on a function of sigma, with a generic f
f = sp.Function('f')
wave = sp.diff(f(S), x, 2) - sp.diff(f(S), t, 2)
target = a**2*(S*sp.diff(f(sig), sig, 2) + sp.diff(f(sig), sig)).subs(sig, S)
ok("(2) (d_x^2 - d_t^2) f(sigma) = a^2 (sigma f'' + f') for a generic f", sp.simplify(wave.doit() - target.doit()) == 0)

# (3)-(4) sigma f'' + f' = f: the series and its sum
n = sp.symbols('n', integer=True, nonnegative=True)
cn = lambda m: 1/sp.factorial(m)**2
ok("(3)-(4) sigma f'' + f' on sum c_n sigma^n has coefficient (n+1)^2 c_{n+1}; with c_n = 1/(n!)^2 it equals c_n",
   sp.simplify((n + 1)*n*cn(n + 1) + (n + 1)*cn(n + 1) - (n + 1)**2*cn(n + 1)) == 0 and sp.simplify((n + 1)**2*cn(n + 1) - cn(n)) == 0)
g0s = sum(sig**m/sp.factorial(m)**2 for m in range(14))
res = sp.expand(sig*sp.diff(g0s, sig, 2) + sp.diff(g0s, sig) - g0s)
ok("(4) the truncated series solves (3) up to its last term only", sp.Poly(res, sig).degree() == 13 and len(sp.Poly(res, sig).terms()) == 1)
# the second solution is singular at sigma = 0: K_0(2 sqrt(sigma)) solves (3) and diverges logarithmically
s_ = sp.symbols('s', positive=True)
K = sp.besselk(0, 2*sp.sqrt(s_))
mp.mp.dps = 30
resK = lambda v: mp.diff(lambda y: mp.besselk(0, 2*mp.sqrt(y)), v, 2)*v + mp.diff(lambda y: mp.besselk(0, 2*mp.sqrt(y)), v) - mp.besselk(0, 2*mp.sqrt(v))
ok("(4') the other solution of (3) is K_0(2 sqrt(sigma)) (residual %.1e at sigma = 0.7) and it is singular at sigma = 0: K_0 = %.2f at sigma = 1e-12"
   % (abs(resK(mp.mpf('0.7'))), mp.besselk(0, 2*mp.sqrt(mp.mpf('1e-12')))), abs(resK(mp.mpf('0.7'))) < 1e-20 and mp.besselk(0, 2e-6) > 10)

# (5) the + line on the ansatz (x + t) h(sigma)
h = sp.Function('h')
lhs5 = Dp((x + t)*h(S)).doit()
ok("(5) (d_x + d_t)[(x + t) h(sigma)] = 2 h + 2 sigma h' for a generic h",
   sp.simplify(lhs5 - (2*h(sig) + 2*sig*sp.diff(h(sig), sig)).subs(sig, S).doit()) == 0)

# (6) h + sigma h' = (a c0/2) g0 by the series
dn = lambda m: 1/(sp.factorial(m)*sp.factorial(m + 1))
ok("(6) h + sigma h' on sum d_n sigma^n has coefficient (n+1) d_n; with d_n = 1/(n!(n+1)!) it equals 1/(n!)^2", sp.simplify((n + 1)*dn(n) - cn(n)) == 0)
# consistency: the - line.  (d_x - d_t) g0(sigma) = (a^2/2)(x + t) g0', and a Q+bulk = (a^2 c0/2)(x + t) g1: needs g0' = g1
ok("(6') the - line holds too: (d_x - d_t) c0 g0 = a (a c0/2)(x + t) g1 because g0' = g1, i.e. (n+1) c_{n+1} = d_n",
   sp.simplify(Dm(f(S)).doit() - (a**2/2*(x + t)*sp.diff(f(sig), sig)).subs(sig, S).doit()) == 0 and sp.simplify((n + 1)*cn(n + 1) - dn(n)) == 0)

# (7) c0 from the jump across the light line x = t, where the no-flip delta sits.  - line: (d_x - d_t) Q- = a Q+, Q+ = delta(x - t) + bulk.
#     clock x: Q- = c0 g0 on x > t.  (d_x - d_t) step(x - t) = 2 delta(x - t)  ->  2 c0 = a.
#     t-reading: Q- = c0 g0 on t > x. (d_x - d_t) step(t - x) = -2 delta(t - x) ->  -2 c0 = a.
w = sp.symbols('w', real=True)
jump_out = sp.diff(sp.Heaviside(w), w).subs(w, x - t)*(sp.diff(x - t, x) - sp.diff(x - t, t))
jump_in = sp.diff(sp.Heaviside(w), w).subs(w, t - x)*(sp.diff(t - x, x) - sp.diff(t - x, t))
c0 = sp.symbols('c0', real=True)
sol_out = sp.solve(sp.Eq((c0*jump_out).coeff(sp.DiracDelta(x - t)), a), c0)
sol_in = sp.solve(sp.Eq((c0*jump_in).coeff(sp.DiracDelta(t - x)), a), c0)
ok(f"(7) the delta on the light line fixes c0: support x > t (the walk's clock) gives c0 = {sol_out[0]}; support t > x (the t-reading) gives c0 = {sol_in[0]}",
   sol_out == [a/2] and sol_in == [-a/2])

# (8) the sums: g0 = I_0(2 sqrt sigma) = J_0(2 sqrt(-sigma)), g1 = I_1(2 sqrt sigma)/sqrt sigma = J_1(2 sqrt(-sigma))/sqrt(-sigma), both signs of sigma
g0 = lambda v: mp.nsum(lambda m: mp.mpf(v)**m/mp.factorial(m)**2, [0, mp.inf])
g1 = lambda v: mp.nsum(lambda m: mp.mpf(v)**m/(mp.factorial(m)*mp.factorial(m + 1)), [0, mp.inf])
e8 = 0
for v in (0.3, 2.0, 7.5):
    e8 = max(e8, abs(g0(v) - mp.besseli(0, 2*mp.sqrt(v))), abs(g0(-v) - mp.besselj(0, 2*mp.sqrt(v))),
             abs(g1(v) - mp.besseli(1, 2*mp.sqrt(v))/mp.sqrt(v)), abs(g1(-v) - mp.besselj(1, 2*mp.sqrt(v))/mp.sqrt(v)))
ok(f"(8) g0(sigma) = I_0(2 sqrt sigma) for sigma > 0 and J_0(2 sqrt(-sigma)) for sigma < 0; g1 likewise with I_1, J_1 (max error {float(e8):.1e})", e8 < 1e-25)
# and the two kernels of the table, from (8) with (7): outside (x > |t|) and inside (t > |x|)
av = mp.mpf('1.3'); e9 = 0
for (tv, xv) in ((0.4, 1.7), (-0.9, 1.2)):                         # outside the cone: III.14 with the loss out
    sv = av**2*(xv**2 - tv**2)/4; ze = av*mp.sqrt(xv**2 - tv**2)
    e9 = max(e9, abs((av/2)*g0(sv) - (av/2)*mp.besseli(0, ze)),
             abs(av*(av/2)/2*(xv + tv)*g1(sv) - (av/2)*(xv + tv)/mp.sqrt(xv**2 - tv**2)*mp.besseli(1, ze)))
for (tv, xv) in ((1.7, 0.4), (1.2, -0.9)):                         # inside the cone: III.15 in the real basis
    sv = av**2*(xv**2 - tv**2)/4; ze = av*mp.sqrt(tv**2 - xv**2)
    e9 = max(e9, abs((-av/2)*g0(sv) - (-(av/2)*mp.besselj(0, ze))),
             abs(av*(-av/2)/2*(xv + tv)*g1(sv) - (-(av/2)*(xv + tv)/mp.sqrt(tv**2 - xv**2)*mp.besselj(1, ze))))
ok(f"(8') with (7): Q- = c0 g0, Q+bulk = (a c0/2)(x + t) g1 are III.14's (a/2) I_0, (a/2)(x+t) I_1/sqrt(x^2-t^2) outside and "
   f"III.15's -(a/2) J_0, -(a/2)(t+x) J_1/sqrt(t^2-x^2) inside (max error {float(e9):.1e})", e9 < 1e-25)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
