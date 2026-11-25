"""Dirac uses the step matrices (the 'matrix noise') at first order: gamma^mu d_mu. Tilt them from point to point:
gamma'^mu(x) = Lambda^mu_nu(x) gamma^nu = S^-1(x) gamma^mu S(x).  What coupling comes out?
1. gamma'^mu d_mu psi = S^-1 gamma^mu [d_mu - (d_mu S) S^-1] (S psi): a connection term appears, valued in the generators of S.
2. 1+1, S = exp(theta(x) A/2) (a local boost): connection = -(1/2)(d_mu theta) gamma^mu A; gamma^mu A = epsilon^{mu nu} gamma_nu: a REAL
   vector-like coupling along the dual gradient of theta: Weyl-type (real), no i. It rescales, not rotates: the length changes.
3. 3+1, S = exp(theta(x) A_x A_y / 2) (a local rotation about z): connection = (1/2)(d_mu theta) gamma^mu A_x A_y = the spin connection
   omega_mu^{xy} = d_mu theta, valued in Sigma_z: the coupling of a rotating local frame (gravity/tetrad), not EM.
4. What EM needs: a generator commuting with ALL step matrices (so it does not mix directions or senses). Real 2x2 commutant of {A, T}:
   multiples of 1 only (a scaling). Real 4x4 commutant of {A_x, A_y, A_z, T}: multiples of 1 only. So no tilt of the matrix noise
   is EM; the generator that acts as i needs the doubling (second copy).
"""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
t, x = sp.symbols('t x', real=True); th = sp.Function('theta', real=True)(t, x)
A = sp.diag(1, -1); T = sp.Matrix([[0, 1], [-1, 0]]); beta = sp.I*T
g = [beta, beta*A]; X = (t, x)
psi = sp.Matrix([sp.Function('p1')(t, x), sp.Function('p2')(t, x)])
S = sp.cosh(th/2)*sp.eye(2) + sp.sinh(th/2)*A
Sinv = sp.cosh(th/2)*sp.eye(2) - sp.sinh(th/2)*A
gp = [sp.simplify(Sinv*gm*S) for gm in g]
lhs = sum((gp[mu]*psi.diff(X[mu]) for mu in range(2)), sp.zeros(2, 1))
chi = S*psi
conn = [sp.simplify(S.diff(X[mu])*Sinv) for mu in range(2)]
rhs = Sinv*sum((g[mu]*(chi.diff(X[mu]) - conn[mu]*chi) for mu in range(2)), sp.zeros(2, 1))
ok("1. gamma'^mu d_mu psi = S^-1 gamma^mu [d_mu - (d_mu S) S^-1] (S psi)", (lhs - rhs).applyfunc(lambda f: sp.simplify(sp.expand(f.rewrite(sp.exp)))) == sp.zeros(2, 1))
ok("2a. (d_mu S) S^-1 = (1/2)(d_mu theta) A: the connection is valued in the boost generator", all(sp.simplify(conn[mu] - th.diff(X[mu])*A/2) == sp.zeros(2) for mu in range(2)))
ok("2b. gamma^0 A = gamma^1, gamma^1 A = gamma^0: the connection term -(1/2)(d_mu theta) gamma^mu A = -(1/2)(d_0 theta gamma^1 + d_1 theta gamma^0), real coefficients",
   g[0]*A == g[1] and g[1]*A == g[0])
# 2c: real vector coupling rescales the length: compare i hbar d_t chi = ... with a real term b gamma^0 gamma^1 chi vs the EM term e Phi chi
b, ePhi = sp.symbols('b ePhi', real=True)
Hreal = sp.I*b*g[0]*g[1]                        # i hbar d_t chi = ... + (i b gamma^0 gamma^1) chi  from a real connection -b gamma^1 ... (schematic)
ok("2c. the real connection enters the Hamiltonian as an anti-Hermitian piece (length not conserved); e Phi * 1 is Hermitian",
   sp.simplify(Hreal + Hreal.H) == sp.zeros(2) and sp.simplify(Hreal - Hreal.H) != sp.zeros(2) and (ePhi*sp.eye(2)).H == ePhi*sp.eye(2))
# 3. 3+1 local rotation
Az = sp.diag(1, -1, -1, 1); Ax = sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]]); Ay = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
T4 = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]]); b4 = sp.I*T4
J = Ax*Ay/2; th3 = sp.Function('theta', real=True)(t, x)
S3 = sp.cos(th3/2)*sp.eye(4) + sp.sin(th3/2)*(2*J)        # exp(theta J), J^2 = -1/4
S3inv = sp.cos(th3/2)*sp.eye(4) - sp.sin(th3/2)*(2*J)
conn3 = sp.simplify(S3.diff(x)*S3inv)
ok("3. 3+1 local rotation about z: (d_x S) S^-1 = (1/2)(d_x theta) A_x A_y = omega_x^{xy} Sigma-type: the spin connection of a rotating frame (Sigma_z = -i A_x A_y)",
   sp.simplify(conn3 - th3.diff(x)*Ax*Ay/2) == sp.zeros(4) and sp.simplify(S3*S3inv) == sp.eye(4))
# 4. commutants
m = sp.Matrix(2, 2, sp.symbols('m0:4', real=True))
sol2 = sp.solve(list(m*A - A*m) + list(m*T - T*m), list(m), dict=True)
M4 = sp.Matrix(4, 4, sp.symbols('n0:16', real=True))
sol4 = sp.solve(list(M4*Ax - Ax*M4) + list(M4*Ay - Ay*M4) + list(M4*Az - Az*M4) + list(M4*T4 - T4*M4), list(M4), dict=True)
m2 = m.subs(sol2[0]); m4 = M4.subs(sol4[0])
ok(f"4. real matrices commuting with all step matrices: 1+1 -> {m2.tolist()}; 3+1 -> multiple of identity: {m4 == m4[0, 0]*sp.eye(4)}",
   m2 == m2[0, 0]*sp.eye(2) and m4 == m4[0, 0]*sp.eye(4))
