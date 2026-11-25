"""III.15c (C), item 8: the only conserved local densities without derivatives are the quadratic forms.   Recreated 9 Nov 2026 (master v34).

The lines read in t (III.5a (H)):  d_t Q+ = -c d_x Q+ + a Q-,   d_t Q- = +c d_x Q- - a Q+.
Ask for a density f(Q+, Q-) and a flux g(Q+, Q-), no derivatives of the readings and no explicit x, t, with d_t f + d_x g = 0 on
every solution.  With u = Q+, v = Q-:
      d_t f + d_x g = (g_u - c f_u) u_x + (g_v + c f_v) v_x + a (v f_u - u f_v),
and u, v, u_x, v_x are independent at a point, so the three brackets vanish.
  1. the substitution, symbolically, for arbitrary functions f, g.
  2. the theorem for polynomials of degree <= 4 (15 coefficients each, as symbols): the solution space is
        f = k0 + (k/2)(u^2 + v^2),   g = g0 + (c k/2)(u^2 - v^2).
  2'. the same by the text's three steps for arbitrary smooth f: f_uv = 0, then v F'(u) = u G'(v).
  3. with a = 0 the third bracket is absent: every F(u) + G(v) is conserved (9 free coefficients at degree <= 4).
  4. the quartic u^3 v scales under a boost exactly like u^2 (both by 1/k, III.15c (6)), and is not conserved: the boost alone
     would not exclude it, the theorem does."""
import sys
import sympy as sp
fails = 0
def ok(name, cond):
    global fails
    print(("PASS  " if cond else "FAIL  ") + name); fails += (not cond)
u, v, ux, vx = sp.symbols('u v u_x v_x', real=True); a, c = sp.symbols('a c', positive=True)
def D(f, g, a_=a):                                    # d_t f + d_x g on the lines
    return sp.expand(sp.diff(f, u)*(-c*ux + a_*v) + sp.diff(f, v)*(c*vx - a_*u) + sp.diff(g, u)*ux + sp.diff(g, v)*vx)

# 1. the three brackets, for arbitrary f, g
F_, G_ = sp.Function('f')(u, v), sp.Function('g')(u, v)
expr = D(F_, G_)
b1 = sp.diff(G_, u) - c*sp.diff(F_, u); b2 = sp.diff(G_, v) + c*sp.diff(F_, v); b3 = a*(v*sp.diff(F_, u) - u*sp.diff(F_, v))
ok("1. d_t f + d_x g = (g_u - c f_u) u_x + (g_v + c f_v) v_x + a (v f_u - u f_v) for arbitrary f, g", sp.simplify(expr - (b1*ux + b2*vx + b3)) == 0)

# 2. polynomials of degree <= 4
mons = [u**i*v**j for i in range(5) for j in range(5 - i)]
fc = sp.symbols('f0:%d' % len(mons)); gc = sp.symbols('g0:%d' % len(mons))
f = sum(cf*m for cf, m in zip(fc, mons)); g = sum(cg*m for cg, m in zip(gc, mons))
def solve_space(a_):
    eqs = sp.Poly(D(f, g, a_), u, v, ux, vx).coeffs()
    sol = sp.solve(eqs, list(fc) + list(gc), dict=True)[0]
    free = [s for s in list(fc) + list(gc) if s not in sol]
    return sp.expand(f.subs(sol)), sp.expand(g.subs(sol)), free
fs, gs, free = solve_space(a)
k = sp.symbols('k')
# the solution must be f = k0 + (k/2)(u^2+v^2), g = g0 + (ck/2)(u^2 - v^2): three free numbers
fu2 = fs.coeff(u, 2).coeff(v, 0)
target_f = fs.coeff(u, 0).coeff(v, 0) + fu2*(u**2 + v**2); target_g = gs.coeff(u, 0).coeff(v, 0) + c*fu2*(u**2 - v**2)
ok("2. degree <= 4, a != 0: %d of the 30 coefficients stay free (f's constant, g's constant, one scale k); "
   "f = const + (k/2)(u^2 + v^2), g = const + (c k/2)(u^2 - v^2): f = %s, g = %s" % (len(free), fs, gs),
   len(free) == 3 and sp.expand(fs - target_f) == 0 and sp.expand(gs - target_g) == 0)

# 2'. the text's steps for arbitrary smooth f
# bracket 1 and 2: g_u = c f_u, g_v = -c f_v; equal mixed derivatives of g: c f_uv = -c f_uv, so f_uv = 0, f = F(u) + G(v)
Fu, Gv = sp.Function('F')(u), sp.Function('G')(v)
third = v*sp.diff(Fu, u) - u*sp.diff(Gv, v)          # bracket 3 on f = F(u) + G(v)
# v F'(u) = u G'(v)  <=>  F'(u)/u = G'(v)/v = k: check that F = k u^2/2, G = k v^2/2 solves it and that it forces F'' constant
solF, solG = k*u**2/2, k*v**2/2
sep = sp.diff(third/(u*v), u)                         # d/du of [F'(u)/u - G'(v)/v] = (F'(u)/u)'
ok("2'. mixed derivatives of g from brackets 1, 2: d_v(c f_u) - d_u(-c f_v) = 2 c f_uv, so f = F(u) + G(v); bracket 3 then reads "
   "v F'(u) = u G'(v), i.e. F'(u)/u = G'(v)/v, a function of u alone equal to a function of v alone: (F'(u)/u)' = 0, "
   "solved by F = k u^2/2, G = k v^2/2",
   sp.simplify(sp.diff(c*sp.diff(F_, u), v) - sp.diff(-c*sp.diff(F_, v), u) - 2*c*sp.diff(F_, u, v)) == 0
   and sp.simplify(third.subs({Fu: solF, Gv: solG}).doit()) == 0
   and sp.simplify(sep - sp.diff(sp.diff(Fu, u)/u, u)) == 0)

# 3. a = 0
fs0, gs0, free0 = solve_space(0)
ok("3. a = 0: %d coefficients stay free; f has no mixed term (f = F(u) + G(v), any F, G of degree <= 4) and g = c F(u) - c G(v) + const: "
   "it is the mass that pairs the two readings into one square" % len(free0),
   len(free0) == 10 and all(sp.Poly(fs0, u, v).coeff_monomial(u**i*v**j) == 0 for i in range(1, 4) for j in range(1, 5 - i))
   and sp.expand(sp.diff(gs0, u) - c*sp.diff(fs0, u)) == 0 and sp.expand(sp.diff(gs0, v) + c*sp.diff(fs0, v)) == 0)

# 4. the quartic
kb = sp.symbols('k_b', positive=True)
scale = lambda e: sp.simplify(e.subs({u: u/sp.sqrt(kb), v: sp.sqrt(kb)*v}, simultaneous=True)/e)
q = u**3*v
gq = sp.symbols('gq0:%d' % len(mons)); gqp = sum(cg*m for cg, m in zip(gq, mons))
eqs = sp.Poly(D(q, gqp), u, v, ux, vx).coeffs()
ok("4. u^3 v scales by %s under the boost of III.15c (6), u^2 by %s: the same; and no flux g of degree <= 4 conserves u^3 v on the lines "
   "(the equations for g have no solution): the boost does not exclude it, conservation does" % (scale(q), scale(u**2)),
   scale(q) == 1/kb and scale(u**2) == 1/kb and sp.solve(eqs, list(gq), dict=True) == [])
print("ALL PASS" if not fails else "%d FAILED" % fails); sys.exit(1 if fails else 0)
