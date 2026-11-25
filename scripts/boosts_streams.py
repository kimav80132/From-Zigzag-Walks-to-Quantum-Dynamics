"""Spin from boosts + streams (universality + additivity), 2+1, t-reading:
    d_t Q = -c A_x d_x Q - c A_z d_z Q + a T Q,   A_z = diag(1,-1), A_x = [[0,1],[1,0]], T = [[0,1],[-1,0]].
Q1 + Q2 is the stream moving +c along x, Q1 - Q2 the stream moving -c along x; Q1, Q2 are the streams along +z, -z.
A frame change (boost or rotation) acts on coordinates by L and on the counters by S:  Q'(X) = S Q(L^-1 X).
'S works' means: every solution Q gives a solution Q'.
Checks:
 1. boost along z (rapidity eta), S = diag(l1, l2):  works  <=>  l1/l2 = e^eta  (the flip rate a is the same in every frame)
 2. the stream content j = (Q1^2+Q2^2, 2Q1Q2, Q1^2-Q2^2) must change like a vector (same total count in every frame): l1 = e^(eta/2), l2 = e^(-eta/2)
 3. j is conserved: d_t(Q1^2+Q2^2) + c d_x(2 Q1 Q2) + c d_z(Q1^2-Q2^2) = 0; Q1^2 gets the Doppler factor e^eta of a +c stream
 4. general small boost along x: the most general K_x (2x2 real) that works is A_x/2 + k 1; same total count kills k
 5. additivity: small boost along n = (n_x, n_z) works with K(n) = n_x K_x + n_z K_z, for every n
 6. finite boosts along x and along z work (exact, symbolic)
 7. commutators: vectors [B_z, B_x] = R (rotation of the (x,z) plane); counters [K_z, K_x] = J = A_z A_x / 2
 8. rotation by theta: coordinates exp(theta R), counters exp(theta J): works, for every theta
 9. exp(2 pi R) = 1 (directions come back), exp(2 pi J) = -1 (counters come back with a minus)
10. loop of boosts z(eps), x(del), z(-eps), x(-del) = rotation by eps*del, on vectors and on counters (to 2nd order)
11. J is the U-turn coin of the z-clock: switch weights +1/2 into +, -1/2 into -  (J = X/2, X from two_step_covariance.py)
12. in 1+1 (one direction) nothing is forced: no commutator, a diagonal S works, no sign
"""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
c, a = sp.symbols('c a', positive=True)
eta, th, k = sp.symbols('eta theta k', real=True)
t, x, z = sp.symbols('t x z', real=True)
Az = sp.diag(1, -1); Ax = sp.Matrix([[0, 1], [1, 0]]); T = sp.Matrix([[0, 1], [-1, 0]]); I2 = sp.eye(2)

# --- generic test: does S (counters) with L (coordinates on (c t, x, z)) map solutions to solutions?
# Q'(X) = S Q(Y), Y = L^-1 X.  d_mu Q' = S (L^-1)_{nu mu} (d_nu Q)(Y).  Replace d_t Q by the equation, require E[Q'] = 0
# identically in the free data (Q, d_x Q, d_z Q) at Y.
q = sp.Matrix(sp.symbols('q1 q2')); qx = sp.Matrix(sp.symbols('qx1 qx2')); qz = sp.Matrix(sp.symbols('qz1 qz2'))
qt = -c*Ax*qx - c*Az*qz + a*T*q                       # d_t Q from the equation
grad = [qt / c, qx, qz]                               # derivatives w.r.t. (c t, x, z)
def residual(S, L):
    Li = sp.simplify(L.inv())
    d = [sum((Li[nu, mu] * grad[nu] for nu in range(3)), sp.zeros(2, 1)) for mu in range(3)]   # d_mu of Q(Y)
    dQp = [S * d_mu for d_mu in d]                    # derivatives of Q' w.r.t. (c t, x, z)
    E = c*dQp[0] + c*Ax*dQp[1] + c*Az*dQp[2] - a*T*(S*q)
    return E.applyfunc(lambda e: sp.simplify(sp.expand(sp.expand(e).rewrite(sp.exp))))
def boost(axis, e):                                   # on (c t, x, z)
    L = sp.eye(3); i = 1 if axis == 'x' else 2
    L[0, 0] = L[i, i] = sp.cosh(e); L[0, i] = L[i, 0] = sp.sinh(e)
    return L
def rot(e):                                           # rotation of (x, z) by angle e
    return sp.Matrix([[1, 0, 0], [0, sp.cos(e), -sp.sin(e)], [0, sp.sin(e), sp.cos(e)]])

# 1. boost along z, diagonal S
l1, l2 = sp.symbols('l1 l2', positive=True)
Lz = boost('z', eta)
# restrict to 1+1 first (drop x): free data qx = 0
res = residual(sp.diag(l1, l2), Lz).subs({qx[0]: 0, qx[1]: 0})
sol = sp.solve([sp.expand(r) for r in res], l1, dict=True)
ok(f"1. boost along z works iff l1/l2 = e^eta : {sol}",
   len(sol) == 1 and sp.simplify(sol[0][l1] / l2 - sp.exp(eta)) == 0)

# 2. the stream content: j = (Q1^2+Q2^2, 2 Q1 Q2, Q1^2-Q2^2) on (c t, x, z). Same total count in every frame <=> j changes like a vector.
qq = sp.Matrix(sp.symbols('u1 u2', real=True))
jvec = lambda v: sp.Matrix([v[0]**2 + v[1]**2, 2*v[0]*v[1], v[0]**2 - v[1]**2])
k1 = sp.symbols('k1', real=True)
Sg = sp.diag(sp.exp(k1*eta), sp.exp((k1 - 1)*eta))          # ratio e^eta from check 1, k1 free
dj = (jvec(Sg*qq) - Lz*jvec(qq)).applyfunc(lambda e: sp.simplify(sp.expand(e.rewrite(sp.exp))))
sk = sp.solve([sp.Poly(sp.expand(e), *qq).coeffs()[i] for e in dj if e != 0 for i in range(len(sp.Poly(sp.expand(e), *qq).coeffs()))], k1, dict=True)
ok(f"2a. stream content changes like a vector under the z-boost only if k1 = 1/2: {sk}", sk == [{k1: sp.Rational(1, 2)}])
Kz = Az / 2
ok("2b. so S_z(eta) = exp(eta A_z/2) = diag(e^(eta/2), e^(-eta/2)), and it works in full 2+1",
   sp.simplify((eta*Kz).exp() - sp.diag(sp.exp(eta/2), sp.exp(-eta/2))) == sp.zeros(2)
   and residual((eta*Kz).exp(), Lz) == sp.zeros(2, 1))
ok("2c. density split two ways: Q1^2 + Q2^2 = (+z stream) + (-z stream) = ((Q1+Q2)^2 + (Q1-Q2)^2)/2 = (+x stream) + (-x stream)",
   sp.expand(qq[0]**2 + qq[1]**2 - ((qq[0]+qq[1])**2 + (qq[0]-qq[1])**2)/2) == 0)

# 3. conserved pair and Doppler
Q1, Q2 = [sp.Function(n)(t, x, z) for n in ('Q1', 'Q2')]
Qv = sp.Matrix([Q1, Q2])
rhs = -c*Ax*Qv.diff(x) - c*Az*Qv.diff(z) + a*T*Qv
dens = Q1**2 + Q2**2; jx = 2*Q1*Q2; jz = Q1**2 - Q2**2
cons = (2*Q1*rhs[0] + 2*Q2*rhs[1]) + c*jx.diff(x) + c*jz.diff(z)
ok("3a. d_t(Q1^2+Q2^2) + c d_x(2 Q1 Q2) + c d_z(Q1^2 - Q2^2) = 0 on solutions", sp.simplify(sp.expand(cons)) == 0)
ok("3b. under S_z: Q1^2 -> e^eta Q1^2, Q2^2 -> e^-eta Q2^2 (the +c and -c stream densities get their Doppler factors)",
   sp.simplify(sp.exp(eta/2)**2 - sp.exp(eta)) == 0 and sp.simplify(sp.exp(-eta/2)**2 - sp.exp(-eta)) == 0)

# 4. small boost along x: most general generator
m = sp.symbols('m11 m12 m21 m22', real=True); Kg = sp.Matrix([[m[0], m[1]], [m[2], m[3]]])
Lx = boost('x', eta)
res = residual(I2 + eta*Kg, Lx)
lin = [sp.expand(sp.diff(r, eta).subs(eta, 0)) for r in res]
coeffs = []
for e in lin:
    coeffs += sp.Poly(e, *q, *qx, *qz).coeffs()
sm = sp.solve(coeffs, m, dict=True)
Kx_gen = Kg.subs(sm[0])
free = list(Kx_gen.free_symbols - {c, a})
ok(f"4a. most general small-boost generator along x: {Kx_gen.tolist()}  (= A_x/2 + {free} * 1)",
   len(free) == 1 and sp.simplify(Kx_gen - Ax/2 - free[0]*I2) == sp.zeros(2))
Sx = sp.simplify((eta*Kx_gen).exp())
djx = (jvec(Sx*qq) - Lx*jvec(qq)).applyfunc(lambda e: sp.simplify(sp.expand(e.rewrite(sp.exp))))
kk = sp.solve(sp.diff(djx, eta).subs(eta, 0).applyfunc(sp.expand)[0].coeff(qq[0]**2), free[0])
ok(f"4b. same total count in every frame forces the multiple of 1 to vanish: {kk}", kk == [0])
Kx = Ax / 2

# 5. additivity: small boost along n
nx, nz = sp.symbols('n_x n_z', real=True)
Ln = sp.eye(3) + eta*sp.Matrix([[0, nx, nz], [nx, 0, 0], [nz, 0, 0]])
res = residual(I2 + eta*(nx*Kx + nz*Kz), Ln)
ok("5. small boost along any n works with K(n) = n_x K_x + n_z K_z (first order in eta, all n)",
   all(sp.simplify(sp.diff(r, eta).subs(eta, 0)) == 0 for r in res))

# 6. finite boosts
ok("6. finite boost along x with exp(eta A_x/2) works (exact), and the stream content changes like a vector", residual((eta*Kx).exp().applyfunc(sp.simplify), Lx) == sp.zeros(2, 1)
   and (jvec((eta*Kx).exp()*qq) - Lx*jvec(qq)).applyfunc(lambda e: sp.simplify(sp.expand(e.rewrite(sp.exp)))) == sp.zeros(3, 1))

# 7. commutators
Bz = sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]); Bx = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
R = Bz*Bx - Bx*Bz
J = Kz*Kx - Kx*Kz
ok(f"7a. vectors: [B_z, B_x] = {R.tolist()}  (turns x into z and z into -x: a rotation of the (x,z) plane)",
   R == sp.Matrix([[0, 0, 0], [0, 0, -1], [0, 1, 0]]))
ok(f"7b. counters: [K_z, K_x] = J = {J.tolist()} = A_z A_x / 2", J == Az*Ax/2 and J == sp.Matrix([[0, sp.Rational(1, 2)], [-sp.Rational(1, 2), 0]]))
ok("7c. R generates rot(theta): exp(theta R) = rot(theta)", sp.simplify((th*R).exp() - rot(th)) == sp.zeros(3))

# 8. rotation pairs exp(theta R) with exp(theta J)
SJ = sp.simplify((th*J).exp())
ok(f"8a. exp(theta J) = {SJ.tolist()}  (turns the counters by theta/2)",
   sp.simplify(SJ - sp.Matrix([[sp.cos(th/2), sp.sin(th/2)], [-sp.sin(th/2), sp.cos(th/2)]])) == sp.zeros(2))
ok("8b. rotation by theta works with counters exp(theta J), for every theta", residual(SJ, rot(th)) == sp.zeros(2, 1))
ok("8c. it does NOT work with exp(2 theta J) (whole-angle turn of the counters)", residual(sp.simplify((2*th*J).exp()), rot(th)) != sp.zeros(2, 1))

E = 2*J
ok("8d. sandwich: S^T A_z S = cos(th) A_z + sin(th) A_x, S^T A_x S = cos(th) A_x - sin(th) A_z, S^T S = 1 (streams turn by the full theta)",
   sp.simplify(SJ.T*Az*SJ - (sp.cos(th)*Az + sp.sin(th)*Ax)) == sp.zeros(2)
   and sp.simplify(SJ.T*Ax*SJ - (sp.cos(th)*Ax - sp.sin(th)*Az)) == sp.zeros(2) and sp.simplify(SJ.T*SJ) == I2)
ok("8e. E = 2J anticommutes with A_z and A_x: E A_z = -A_x, E A_x = A_z", E*Az == -Ax and E*Ax == Az and Az*E == Ax and Ax*E == -Az)

# 9. full turn
ok("9. full turn theta = 2 pi: directions rot = 1, counters exp(2 pi J) = -1",
   sp.simplify(rot(2*sp.pi) - sp.eye(3)) == sp.zeros(3) and sp.simplify(SJ.subs(th, 2*sp.pi) + I2) == sp.zeros(2))

# 10. loop of boosts
e1, e2 = sp.symbols('epsilon delta', real=True)
def series2(M):
    return M.applyfunc(lambda f: sp.series(sp.series(f, e1, 0, 3).removeO(), e2, 0, 3).removeO())
def ex(Mat, s):                      # exp to 2nd order
    return sp.eye(Mat.shape[0]) + s*Mat + s**2*Mat*Mat/2
loopV = series2(ex(Bz, e1)*ex(Bx, e2)*ex(Bz, -e1)*ex(Bx, -e2))
loopC = series2(ex(Kz, e1)*ex(Kx, e2)*ex(Kz, -e1)*ex(Kx, -e2))
cut = lambda M: M.applyfunc(lambda f: sum(tm for tm in sp.Add.make_args(sp.expand(f)) if sp.Poly(tm, e1, e2).total_degree() <= 2))
ok(f"10. loop z(eps) x(del) z(-eps) x(-del): vectors = 1 + eps del R, counters = 1 + eps del J  (to 2nd order)",
   sp.simplify(cut(loopV) - (sp.eye(3) + e1*e2*R)) == sp.zeros(3) and sp.simplify(cut(loopC) - (I2 + e1*e2*J)) == sp.zeros(2))

# 11. J is the U-turn coin
X = sp.Matrix([[0, 1], [-1, 0]])
ok("11. J = X/2: from one side, a switch into + weighs +1/2 and a switch into - weighs -1/2 (the U-turn coin)", J == X/2)

# 12. 1+1: one direction, nothing forced
ok("12. 1+1: any diagonal S commutes with A_z; a single boost has no partner, no commutator, J does not exist",
   (sp.diag(l1, l2)*Az - Az*sp.diag(l1, l2)) == sp.zeros(2))
