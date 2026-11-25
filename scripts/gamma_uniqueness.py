"""IV.4a — the extension to three dimensions is forced.
From the process: each spatial coin is a reflection (squares to +1); the time-sense turn squares to -1 (III.6, III.8).
Conditions for three dimensions (IV.1): one speed c in EVERY direction, and Klein-Gordon with the right sign in every direction.
Shown here:
  1. (n.A)^2 = 1 for all unit n  <=>  A_j^2 = 1 and A_j A_k + A_k A_j = 0 (j != k): isotropy forces perpendicular mirrors.
  2. the turn T must anticommute with every A_j and square to -1 (sign of the mass term in every direction).
  3. no solution in odd dimension (determinant), none in dimension 2 (the anticommutant of sigma_x,y,z is {0}).
  4. a solution exists in dimension 4 (two are written: Weyl and Dirac bases), and they are the same up to a change of basis
     (the similarity S is found by solving the linear equations S gamma_W = gamma_D S); Pauli's theorem says this holds for
     every solution in dimension 4.
"""
import numpy as np, itertools
import sympy as sp
rng = np.random.default_rng(0)
ok_all = True
def ok(name, cond):
    global ok_all; ok_all &= bool(cond); print(("PASS  " if cond else "FAIL  ") + name)

# ---- 1. isotropy <=> perpendicular mirrors (symbolic in n, generic noncommuting A)
A1, A2, A3 = sp.symbols('A1 A2 A3', commutative=False)
n1, n2, n3 = sp.symbols('n1 n2 n3', real=True)
sq = sp.expand((n1 * A1 + n2 * A2 + n3 * A3)**2)
# collect the coefficient of each monomial in n: the identity (n.A)^2 = n.n for all n requires
want = {n1**2: A1**2, n2**2: A2**2, n3**2: A3**2, n1 * n2: A1 * A2 + A2 * A1, n1 * n3: A1 * A3 + A3 * A1, n2 * n3: A2 * A3 + A3 * A2}
# extract coefficients by differentiating in the (commuting) n's
coeffs = {n1**2: sp.diff(sq, n1, 2) / 2, n2**2: sp.diff(sq, n2, 2) / 2, n3**2: sp.diff(sq, n3, 2) / 2,
          n1 * n2: sp.diff(sq, n1, n2), n1 * n3: sp.diff(sq, n1, n3), n2 * n3: sp.diff(sq, n2, n3)}
ok("1. (n.A)^2 = sum n_j^2 A_j^2 + sum_{j<k} n_j n_k (A_j A_k + A_k A_j): so '= n.n for all n' forces A_j^2 = 1, {A_j,A_k} = 0",
   all(sp.expand(coeffs[k] - v) == 0 for k, v in want.items()))

# ---- 2. the turn: apply (c d_x- ... ) as in III.6 in direction n: (n.A d + T m)(... ) cross terms vanish iff T anticommutes
# (checked as an algebraic identity on matrices below, together with 4.)

# ---- 3. dimension: odd impossible; 2 impossible
# odd: if A T = -T A with A, T invertible then det(A)det(T) = (-1)^d det(T)det(A) -> d even.
detsign = all(sp.expand((-sp.Matrix(d, d, sp.symbols(f'm0:{d*d}'))).det() - (-1)**d*sp.Matrix(d, d, sp.symbols(f'm0:{d*d}')).det()) == 0 for d in (1, 2, 3, 4, 5))
# and directly for d = 3: with A a mirror (A^2 = 1) of any signature, the matrices T with AT = -TA are all singular
d3 = True
for sig in ((1, 1, 1), (1, 1, -1), (1, -1, -1), (-1, -1, -1)):
    A3d = sp.diag(*sig); Tm = sp.Matrix(3, 3, sp.symbols('t0:9'))
    Ts = Tm.subs(sp.solve(list(A3d*Tm + Tm*A3d), list(Tm), dict=True)[0])
    d3 = d3 and sp.expand(Ts.det()) == 0
ok("3a. odd dimension impossible: det(-M) = (-1)^d det(M) (d = 1..5, symbolic), so AT = -TA with both invertible forces d even; "
   "and directly in d = 3: for a mirror A of any signature every T with AT = -TA has determinant 0", detsign and d3)
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1, -1]).astype(complex)
# the set of 2x2 X with X s + s X = 0 for s in {sx,sy,sz}: linear system in the 4 entries of X
rows = []
for s in (sx, sy, sz):
    M = np.kron(np.eye(2), s) + np.kron(s.T, np.eye(2))      # vec(sX + Xs) = (I⊗s + s^T⊗I) vec(X)
    rows.append(M)
M = np.vstack(rows)
ok("3b. dimension 2: nothing anticommutes with all three of sigma_x, sigma_y, sigma_z  (null space of the linear map is {0})",
   np.linalg.matrix_rank(M) == 4)
print("      (so a fourth perpendicular mirror needs more than two components, and by 3a the next candidate is 4)")

# ---- 4. dimension 4: two solutions and their similarity
Z = np.zeros((2, 2)); I2 = np.eye(2)
def blk(a, b, c, d): return np.block([[a, b], [c, d]])
A_W = [blk(s, Z, Z, -s) for s in (sx, sy, sz)]; T_W = blk(Z, I2, I2, Z) * 1j        # Weyl: T = i beta, T^2 = -1
A_D = [blk(Z, s, s, Z) for s in (sx, sy, sz)]; T_D = blk(I2, Z, Z, -I2) * 1j        # Dirac: T = i beta
def clifford_ok(A, T):
    c1 = all(np.allclose(a @ a, np.eye(4)) for a in A)
    c2 = all(np.allclose(A[j] @ A[k] + A[k] @ A[j], 0) for j in range(3) for k in range(3) if j != k)
    c3 = all(np.allclose(a @ T + T @ a, 0) for a in A) and np.allclose(T @ T, -np.eye(4))
    return c1 and c2 and c3
ok("4a. Weyl basis: three perpendicular mirrors + a turn anticommuting with all, T^2 = -1", clifford_ok(A_W, T_W))
ok("4b. Dirac basis: the same relations", clifford_ok(A_D, T_D))
# KG in every direction: (n.A d + T a)^2 = d^2 - a^2 (cross terms vanish) -- check with numbers d, a and random n
n = rng.normal(size=3); n /= np.linalg.norm(n); d, a = 0.7, 1.3
Op = d * sum(n[j] * A_W[j] for j in range(3)) + a * T_W
ok("2.  (n.A d + a T)^2 = (d^2 - a^2) 1 in every direction n: Klein-Gordon with the right sign, no cross term", np.allclose(Op @ Op, (d * d - a * a) * np.eye(4)))
# similarity: solve S G_W = G_D S for all four generators
G_W = A_W + [T_W]; G_D = A_D + [T_D]
rows = [np.kron(g.T, np.eye(4)) - np.kron(np.eye(4), gd) for g, gd in zip(G_W, G_D)]   # vec(S g) - vec(gd S) with vec column-major
Mbig = np.vstack(rows)
u, sv, vh = np.linalg.svd(Mbig)
null = vh[-1].conj()
S = null.reshape(4, 4, order='F')
ok("4c. the two bases are related by an invertible S with S gamma_W S^-1 = gamma_D (null space of the similarity equations is one-dimensional)",
   (sv < 1e-10).sum() == 1 and abs(np.linalg.det(S)) > 1e-8 and all(np.allclose(S @ g @ np.linalg.inv(S), gd) for g, gd in zip(G_W, G_D)))
# 4d. the two sets written on the page: IV.4's (A_k = diag(sigma_k, -sigma_k), T real, the turn between the copies) and IV.4b (I)'s all-real set.
#     (Added 10 Nov 2026: 4a above uses T = i beta, which is not the T printed in IV.4.)
T_IV4 = blk(Z, I2, -I2, Z).astype(complex)
A_real = [np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]], complex), np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]], complex),
          np.diag([1, -1, -1, 1]).astype(complex)]
T_real = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]], complex)
def similar(G1, G2):
    Mb = np.vstack([np.kron(g.T, np.eye(4)) - np.kron(np.eye(4), gd) for g, gd in zip(G1, G2)])
    u_, sv_, vh_ = np.linalg.svd(Mb); S_ = vh_[-1].conj().reshape(4, 4, order='F')
    return (sv_ < 1e-10).sum() == 1 and abs(np.linalg.det(S_)) > 1e-8 and all(np.allclose(S_ @ g @ np.linalg.inv(S_), gd) for g, gd in zip(G1, G2))
ok("4d. IV.4's set as printed (A_k = diag(sigma_k, -sigma_k), T = [[0,1],[-1,0]] between the copies) and IV.4b (I)'s real set obey the same relations, "
   "and each is carried into the standard set by an invertible S, unique up to a factor",
   clifford_ok(A_W, T_IV4) and clifford_ok(A_real, T_real) and similar(A_W + [T_IV4], G_D) and similar(A_real + [T_real], G_D))
# 4e. the standard set as IV.9 writes it: gamma^0 = F, gamma^k = F A_k with F = -i T: (gamma^0)^2 = 1, (gamma^k)^2 = -1, all anticommute
F = -1j*T_D; gam = [F] + [F @ a_ for a_ in A_D]
ok("4e. gamma^0 = F, gamma^k = F A_k (F = -iT, F^2 = 1): (gamma^0)^2 = 1, (gamma^k)^2 = -1, all four anticommute (IV.9)",
   np.allclose(F @ F, np.eye(4)) and all(np.allclose(gam[k] @ gam[k], -np.eye(4)) for k in (1, 2, 3))
   and all(np.allclose(gam[i] @ gam[j] + gam[j] @ gam[i], 0) for i in range(4) for j in range(i)))
print("    (Pauli 1936: every set of four 4x4 matrices with these relations is related to this one by such an S.)")
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
