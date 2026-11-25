"""1+1 -> 2+1 by universality + additivity (the dust-ball -> Einstein move).
Solved case (1+1): along one axis, d_t Q = -c A d_s Q + a T Q with A^2 = 1 (movers at +-c), T^2 = -1 (turn), AT + TA = 0.
P1 universality: the same 1+1 law holds along EVERY direction n (same column, same T).
P2 additivity: the directions combine linearly: a wave along n sees A(n) = n_x A_x + n_z A_z.
Checks:
 1. A(n)^2 = 1 for all unit n  <=>  A_x^2 = A_z^2 = 1 and A_x A_z + A_z A_x = 0
 2. with A_z = diag(1,-1) (the 1+1 counters), the real symmetric A_x is forced to be +-sigma_x
 3. A(n) T + T A(n) = 0 for all n  forces T = +-[[0,1],[-1,0]]
 4. diagonal A_x (every counter carries its own x-direction) is impossible
 5. a coin (1 +- A)/2 has no negative weight  <=>  A is diagonal: so the x-coin must carry a sign
 6. the result: d_t^2 Q = c^2 (d_x^2 + d_z^2) Q - a^2 Q, every component
"""
import sympy as sp
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
nx, nz = sp.symbols('n_x n_z', real=True)
p, q, r, s = sp.symbols('p q r s', real=True)
Az = sp.diag(1, -1)
Ax = sp.Matrix([[p, q], [q, r]])                       # general real symmetric
An = nx*Ax + nz*Az
sq = sp.expand(An*An)
# 1. expansion
ok("1. A(n)^2 = n_x^2 A_x^2 + n_z^2 A_z^2 + n_x n_z (A_x A_z + A_z A_x)",
   sp.simplify(sq - (nx**2*Ax*Ax + nz**2*Az*Az + nx*nz*(Ax*Az + Az*Ax))) == sp.zeros(2))
# require = (nx^2 + nz^2) 1 for all nx, nz: coefficients of nx^2, nz^2, nx nz
eqs = []
for i in range(2):
    for j in range(2):
        e = sp.Poly(sq[i, j] - (nx**2 + nz**2)*(1 if i == j else 0), nx, nz)
        eqs += e.coeffs()
sols = sp.solve(eqs, [p, q, r], dict=True)
ok(f"2. with A_z = diag(1,-1), universality forces A_x = [[0,q],[q,0]] with q = +-1: solutions {sols}",
   all(s_[p] == 0 and s_[r] == 0 and abs(s_[q]) == 1 for s_ in sols) and len(sols) == 2)
Ax0 = sp.Matrix([[0, 1], [1, 0]])
T = sp.Matrix([[p, q], [r, s]])
c1 = list(T*Az + Az*T) + list(T*Ax0 + Ax0*T)
solT = sp.solve(c1, [p, q, r, s], dict=True)[0]
Tg = T.subs(solT)
lam = [x for x in Tg.free_symbols][0]
solT2 = sp.solve(list(Tg*Tg + sp.eye(2)), lam)
ok(f"3. T anticommuting with A_z and A_x, T^2 = -1: T = {Tg} with {lam} in {solT2}",
   (Tg == sp.Matrix([[0, lam], [-lam, 0]]) or Tg == sp.Matrix([[0, -lam], [lam, 0]])) and sorted(sp.solve((Tg*Tg + sp.eye(2))[0, 0], lam)) == [-1, 1])
d1, d2 = sp.symbols('d1 d2', real=True)
D = sp.diag(d1, d2)
ok("4. a diagonal A_x cannot anticommute with A_z unless it is zero (then no x-motion)",
   sp.solve(list(D*Az + Az*D), [d1, d2], dict=True) == [{d1: 0, d2: 0}])
a12 = sp.symbols('a12', real=True)
coin_p = (sp.eye(2) + sp.Matrix([[0, a12], [a12, 0]]))/2; coin_m = (sp.eye(2) - sp.Matrix([[0, a12], [a12, 0]]))/2
ok("5. the two coins (1+A)/2 and (1-A)/2 have off-diagonal entries +a12/2 and -a12/2: one is negative unless a12 = 0",
   coin_p[0, 1] == a12/2 and coin_m[0, 1] == -a12/2)
t, x, z, c, a = sp.symbols('t x z c a', real=True)
Q = sp.Matrix([sp.Function('Q1')(t, x, z), sp.Function('Q2')(t, x, z)])
Tm = sp.Matrix([[0, 1], [-1, 0]])
L = lambda g: -c*(Ax0*g.diff(x) + Az*g.diff(z)) + a*Tm*g
ok("6. d_t^2 Q = c^2 (d_x^2 + d_z^2) Q - a^2 Q for every component",
   sp.simplify(L(L(Q)) - (c**2*(Q.diff(x, 2) + Q.diff(z, 2)) - a**2*Q)) == sp.zeros(2, 1))
# 3+1: three directions need four counters (IV.4a).  (Rebuilt 10 Nov 2026: the earlier check tried only entries -1, 0, 1 of a symmetric matrix.)
g = sp.symbols('g0:4', real=True); Ay = sp.Matrix([[g[0], g[1]], [g[2], g[3]]])              # ANY real 2x2
solA = sp.solve(list(Ay*Ax0 + Ax0*Ay) + list(Ay*Az + Az*Ay), list(g), dict=True)[0]
Ayg = Ay.subs(solA); free = list(Ayg.free_symbols)
ok(f"7. every real 2x2 matrix anticommuting with A_x and A_z is {Ayg.tolist()}: a multiple of the turn, with square -{free[0]}^2 <= 0, never +1. "
   "A third direction does not fit in two real counters (so 3+1 needs four, IV.4a)",
   len(free) == 1 and sp.simplify(Ayg + Ayg.T) == sp.zeros(2) and sp.simplify(Ayg*Ayg + free[0]**2*sp.eye(2)) == sp.zeros(2))
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
