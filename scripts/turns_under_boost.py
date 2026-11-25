"""The mass is a turn in time (III.5/III.9a). What does a boosted observer see?
1. 1+1, rest: Q = (cos at, -sin at) solves d_t Q = -c A d_x Q + a T Q with d_x = 0: a turn of (Q+,Q-) at rate a in t.
2. Boost by eta: the same solution seen moving is a plane wave with omega = a cosh(eta), k = (a/c) sinh(eta): the turn advances in t AND x,
   omega^2 - c^2 k^2 = a^2: the turn per unit proper time is the invariant. The column is turned by S = exp(eta A/2) (IV.4c).
3. Generators (the antisymmetric part of the matrix noise, 3+1): time-space turns A_k (boosts), space-space turns A_j A_k (rotations).
   Under a boost along x, S = exp(eta A_x/2), they mix:  A_z -> cosh A_z + sinh (A_x A_z) ... exactly as (E, B) mix. The set of six is invariant.
4. The mass turn A_z T (exchange, z-clock) is unchanged by the x-boost (A_x commutes with A_z T): the exchange is the same for every runner,
   which was the universality input of IV.4c. The equation keeps its matrices; what mixes is the generators acting on the state."""
import sympy as sp
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
t, x, eta = sp.symbols('t x eta', real=True); a, c = sp.symbols('a c', positive=True)
A = sp.diag(1, -1); T = sp.Matrix([[0, 1], [-1, 0]])
Q0 = sp.Matrix([sp.cos(a*t), -sp.sin(a*t)])
ok("1. rest: Q = (cos at, -sin at) solves d_t Q = a T Q (a turn at rate a)", sp.simplify(Q0.diff(t) - a*T*Q0) == sp.zeros(2, 1))
# boost: t_rest = cosh(eta) t - sinh(eta) x/c ; column S Q
S = sp.cosh(eta/2)*sp.eye(2) + sp.sinh(eta/2)*A
Qb = S*Q0.subs(t, sp.cosh(eta)*t - sp.sinh(eta)*x/c)
res = (Qb.diff(t) + c*A*Qb.diff(x) - a*T*Qb).applyfunc(lambda f: sp.simplify(sp.expand(f.rewrite(sp.exp))))
ok("2a. the boosted solution S Q(t_rest) solves the same lines", res == sp.zeros(2, 1))
om, k = a*sp.cosh(eta), a*sp.sinh(eta)/c
ok("2b. it oscillates as cos(omega t - k x) with omega = a cosh(eta), k = (a/c) sinh(eta); omega^2 - c^2 k^2 = a^2 (the invariant turn rate)",
   sp.simplify(sp.expand((a*(sp.cosh(eta)*t - sp.sinh(eta)*x/c)) - (om*t - k*x))) == 0 and sp.simplify(om**2 - c**2*k**2 - a**2) == 0)
# 3 generators in 3+1
Az = sp.diag(1, -1, -1, 1); Ax = sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]]); Ay = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
T4 = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]])
Sx = sp.cosh(eta/2)*sp.eye(4) + sp.sinh(eta/2)*Ax; Sxi = sp.cosh(eta/2)*sp.eye(4) - sp.sinh(eta/2)*Ax
conj = lambda M: (Sxi*M*Sx).applyfunc(lambda f: sp.simplify(f.rewrite(sp.exp)))
ex = lambda M: M.applyfunc(lambda f: f.rewrite(sp.exp))
ok("3a. boost along x: A_x -> A_x (unchanged); A_z -> cosh(eta) A_z - sinh(eta) A_x A_z (a time-space turn picks up a space-space part)",
   conj(Ax) == Ax and sp.simplify(conj(Az) - ex(sp.cosh(eta)*Az - sp.sinh(eta)*Ax*Az)) == sp.zeros(4))
ex = lambda M: M.applyfunc(lambda f: f.rewrite(sp.exp))
ok("3b. A_x A_z (space-space, in the boost plane) -> cosh A_x A_z - sinh A_z: picks up a time-space part;  A_y A_z (plane perpendicular to the boost) -> unchanged",
   sp.simplify(conj(Ax*Az) - ex(sp.cosh(eta)*Ax*Az - sp.sinh(eta)*Az)) == sp.zeros(4) and conj(Ay*Az) == Ay*Az)
ok("3c'. the pair (A_y, A_x A_y) mixes like (E_y, B_z): A_y -> cosh A_y - sinh A_x A_y,  A_x A_y -> cosh A_x A_y - sinh A_y",
   sp.simplify(conj(Ay) - ex(sp.cosh(eta)*Ay - sp.sinh(eta)*Ax*Ay)) == sp.zeros(4) and sp.simplify(conj(Ax*Ay) - ex(sp.cosh(eta)*Ax*Ay - sp.sinh(eta)*Ay)) == sp.zeros(4))
gens = [Ax, Ay, Az, Ax*Ay, Ay*Az, Az*Ax]
span_ok = all(sp.Matrix.hstack(*[g.reshape(16, 1) for g in gens], conj(G).reshape(16, 1)).rank() == 6 for G in gens)
ok("3c. every boosted generator stays in the span of the six (3 time-space + 3 space-space turns): the set is closed", span_ok)
ok("4. the exchange A_z T is unchanged by the x-boost (same flip exchange for every runner: IV.4c's input)", conj(Az*T4) == Az*T4)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
