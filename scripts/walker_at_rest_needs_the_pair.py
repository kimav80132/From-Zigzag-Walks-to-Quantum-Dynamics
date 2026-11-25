"""Exploration, 11 Nov 2026 (not cited in the master): a test of the user's hypothesis that the second list (the twin) need not be posited.

The master's generalization to 3+1 (IV.4a-IV.4c) asks for the SMALLEST set of counters on which "speed c in every direction" and "the right-sign mass
term" hold, and finds four real counters (one list).  Ask instead for a walker that can STAY AT REST: a state whose flow is zero and stays zero.
Matrices: IV.4b's real A_x, A_y, A_z, T (3+1); the 2 x 2 blocks of III.8 / IV.4c (E) in 1+1 and 2+1.  The lines read in t:  d_t q = -c A_k d_k q + a T q.
Content of a list q:  density q.q,  flow  c q.A_k q  (IV.4d (B)).  A pair (u, v): the same, added over the two lists.

  1. The sixteen products of A_x, A_y, A_z, T fill all 4 x 4 real matrices (in 1+1: A, T fill all 2 x 2).  Then any set of counters obeying the same
     relations is, in some basis, several lists of four (of two): Wedderburn's theorem, cited, not proved here.  So "more counters" = "more lists".
  2. One list, 3+1 and 2+1: the flow is light-like for EVERY list, density^2 = |flow|^2.  It always moves at c.  1+1: density^2 - flow^2 = 4 Q+^2 Q-^2.
  3. One list cannot stay at rest, in 1+1, 2+1 or 3+1: at rest the flow is cos(2at) q.A_k q + sin(2at) q.A_k T q, and both vanish only for q = 0.
  4. A pair: density^2 - |flow|^2 = (2 u.Tv)^2 + (2 u.(T w)v)^2, w = A_x A_y A_z.  The first form is III.8c (8)'s W (route against twin).  A pair has a
     rest frame exactly when these two forms are not both zero; for one list both are zero.  1+1: density^2 - flow^2 = 4 (K^2 + W^2).
  5. A pair stays at rest in exactly two ways: v = T u or v = -T u (the twin a quarter turn ahead, or behind: III.15e's two classes), in 1+1 and 3+1.
  6. The six plane turns (three runs A_k/2, three rotations A_j A_k/2): one list is moved by them in only four independent ways - two combinations
     leave it unchanged; a pair with a rest frame is moved in six; a pair without one (v = (alpha + beta w) u) again in four.
  7. The frame reading: every pair (u, v) = M (u0, T u0) for one and only one M = [1, three rotations, three runs, w] (eight numbers); one list is
     (1 + three rotations) u0 (four numbers).
  8. What the step never reads on the pair: the matrices commuting with the step are four, {1, J, K, JK}; those that keep the weight and square to -1
     are +-J only.  In the frame reading J is the reference turned in one of its own space planes, the one across its spin axis.
  9. In motion: the rest pair run to rapidity eta along z is in the leading class at ck = a sinh(eta), and its flow is steady.
"""
import sys, itertools
import numpy as np, sympy as sp
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all = ok_all and bool(cond); print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()
rng = np.random.default_rng(7)

Ax = sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]])
Ay = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
Az = sp.diag(1, -1, -1, 1); T = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]]); I4 = sp.eye(4)
AA = (Ax, Ay, Az); w = Ax*Ay*Az
n = lambda M: np.array(M.tolist(), float)
A2 = sp.diag(1, -1); X2 = sp.Matrix([[0, 1], [1, 0]]); T2 = sp.Matrix([[0, 1], [-1, 0]])                 # 1+1 and 2+1 blocks

# ---- 0. the relations
rel = all(A*A == I4 for A in AA) and T*T == -I4 and all(A*B == -B*A for A, B in itertools.combinations(AA + (T,), 2)) \
    and all(A.T == A for A in AA) and T.T == -T
ok("0. IV.4b's real matrices: A_k^2 = 1, T^2 = -1, all six pairs anticommute, the A_k symmetric and T antisymmetric", rel)

# ---- 1. the products fill everything
def span_dim(gens, N):
    words = [np.eye(N)]
    for r_ in range(1, len(gens) + 1):
        for comb in itertools.combinations(gens, r_):
            M = np.eye(N)
            for g in comb: M = M@g
            words.append(M)
    return np.linalg.matrix_rank(np.array([W.ravel() for W in words])), len(words)
d4, nw4 = span_dim([n(M) for M in AA + (T,)], 4); d2, nw2 = span_dim([n(A2), n(T2)], 2)
ok("1. the %d products of A_x, A_y, A_z, T span %d = 4 x 4 dimensions: every real 4 x 4 matrix; in 1+1 the %d products of A, T span %d = 2 x 2. So "
   "(Wedderburn, cited) any counters obeying the same relations are several lists of four (of two): more counters means more lists, nothing else"
   % (nw4, d4, nw2, d2), d4 == 16 and d2 == 4)

# ---- 2. one list always moves at c (3+1, 2+1)
q = sp.Matrix(sp.symbols('q1:5', real=True)); dens = (q.T*q)[0]; flow = [(q.T*A*q)[0] for A in AA]
null4 = sp.expand(dens**2 - sum(f**2 for f in flow)) == 0
p_, m_ = sp.symbols('Qp Qm', real=True); q2 = sp.Matrix([p_, m_])
null3 = sp.expand((q2.T*q2)[0]**2 - (q2.T*A2*q2)[0]**2 - (q2.T*X2*q2)[0]**2) == 0
d11 = sp.expand((q2.T*q2)[0]**2 - (q2.T*A2*q2)[0]**2 - 4*p_**2*m_**2) == 0
ok("2. one list: density^2 - |flow|^2 = 0 identically in 3+1 (four counters, three directions) and in 2+1 (two counters, two directions): one list "
   "always moves at c, whatever is in it. In 1+1: density^2 - flow^2 = 4 Q+^2 Q-^2", null4 and null3 and d11)

# ---- 3. one list cannot stay at rest
# at rest q(t) = (cos at + sin at T) q; A_k anticommutes with T, so flow_k(t) = cos(2at) q.A_k q + sin(2at) q.A_k T q
t_, a_ = sp.symbols('t a', real=True)
Et = sp.cos(a_*t_)*I4 + sp.sin(a_*t_)*T
form = all(sp.simplify(((Et*q).T*A*(Et*q))[0] - (sp.cos(2*a_*t_)*(q.T*A*q)[0] + sp.sin(2*a_*t_)*(q.T*A*T*q)[0])) == 0 for A in AA)
sol11 = sp.solve([(q2.T*A2*q2)[0], (q2.T*A2*T2*q2)[0]], [p_, m_], dict=True)
ok("3. one list cannot stay at rest. At rest the flow is cos(2at) q.A_k q + sin(2at) q.A_k T q (checked). 3+1 and 2+1: by 2 the flow has the size of "
   "the density, so it vanishes only for q = 0. 1+1: Q+^2 - Q-^2 = 0 and 2 Q+ Q- = 0 have the only solution %s" % sol11,
   form and sol11 == [{p_: 0, m_: 0}] and sp.expand((q2.T*A2*T2*q2)[0] - 2*p_*m_) == 0)

# ---- 4. a pair: what makes a rest frame
u = sp.Matrix(sp.symbols('u1:5', real=True)); v = sp.Matrix(sp.symbols('v1:5', real=True))
J0 = (u.T*u)[0] + (v.T*v)[0]; Jk = [(u.T*A*u)[0] + (v.T*A*v)[0] for A in AA]
s_form = 2*(u.T*T*v)[0]; p_form = 2*(u.T*T*w*v)[0]
fierz = sp.expand(J0**2 - sum(j**2 for j in Jk) - s_form**2 - p_form**2) == 0
up, um, vp, vm = sp.symbols('up um vp vm', real=True)
K11 = up*um + vp*vm; W11 = up*vm - vp*um
f11 = sp.expand((up**2 + um**2 + vp**2 + vm**2)**2 - (up**2 - um**2 + vp**2 - vm**2)**2 - 4*(K11**2 + W11**2)) == 0
antis = (T*w).T == -(T*w) and T.T == -T
ok("4. a pair (u, v): density^2 - |flow|^2 = (2 u.T v)^2 + (2 u.(T w) v)^2, w = A_x A_y A_z, identically (3+1). Both forms are antisymmetric, so "
   "both vanish when there is one list only; u.T v is III.8c (8)'s W, the route against the twin. A pair has a rest frame exactly when the two "
   "are not both zero. 1+1: density^2 - flow^2 = 4 (K^2 + W^2)", fierz and f11 and antis)

# ---- 5. a pair stays at rest in exactly two ways
J8 = np.block([[np.zeros((4, 4)), -np.eye(4)], [np.eye(4), np.zeros((4, 4))]]); dbl = lambda M: np.kron(np.eye(2), n(M))
A8 = [dbl(A) for A in AA]; T8 = dbl(T)
steady = lambda r: max(max(abs(r@A@r), abs(r@A@T8@r)) for A in A8)             # the two amplitudes of the flow at rest, largest over k
cls_ok = True
for sg in (1, -1):
    cls_ok = cls_ok and all(sp.expand((u.T*A*u)[0] + ((sg*T*u).T*A*(sg*T*u))[0]) == 0 and sp.expand((u.T*A*T*u)[0] + ((sg*T*u).T*A*T*(sg*T*u))[0]) == 0 for A in AA)
# converse: r = (u, Tu) + (x, -Tx); the conditions on x are x . A_k u = 0 and x . (T A_k u) = 0 for k = x, y, z
V = sp.Matrix.hstack(*[A*u for A in AA], *[T*A*u for A in AA])               # 4 x 6
gram = sp.factor(sp.expand((V*V.T).det()))
generic = float(np.mean([steady(rng.normal(size=8)) for _ in range(200)]))
u11 = sp.Matrix([p_, m_]); det11 = sp.expand(sp.Matrix.hstack(A2*u11, T2*A2*u11).det())
ok("5. a pair stays at rest in exactly two ways. Both v = T u and v = -T u have zero flow at all times (identically). Conversely write the pair as "
   "(u, Tu) + (x, -Tx): the conditions on x are x.A_k u = 0, x.(T A_k u) = 0, and these six vectors fill the four counters whenever u is not zero "
   "(Gram determinant = %s), so x = 0 - or u = 0. A random pair has flow amplitude %.2f. In 1+1 the same: det[A u, T A u] = %s. So: twin a quarter "
   "turn ahead, or a quarter turn behind - III.15e's two classes - and nothing else" % (gram, generic, det11),
   cls_ok and sp.expand(gram - 4*(u.T*u)[0]**4) == 0 and generic > 0.5 and sp.expand(det11 + p_**2 + m_**2) == 0)

# ---- 6. the six plane turns
six = [n(A) for A in AA] + [n(Ay*Az), n(Az*Ax), n(Ax*Ay)]                      # runs, rotations (generators up to 1/2)
wn = n(w); Tn = n(T)
r1 = set(); r2 = set(); r3 = set(); smin = 1e9
for _ in range(300):
    uu, vv = rng.normal(size=4), rng.normal(size=4)
    r1.add(int(np.linalg.matrix_rank(np.array([G@uu for G in six]).T, tol=1e-9)))
    r2.add(int(np.linalg.matrix_rank(np.array([np.concatenate([G@uu, G@vv]) for G in six]).T, tol=1e-9)))
    smin = min(smin, (2*uu@Tn@vv)**2 + (2*uu@Tn@wn@vv)**2)
    al, be = rng.normal(size=2); vd = al*uu + be*wn@uu                         # a pair with no rest frame
    r3.add(int(np.linalg.matrix_rank(np.array([np.concatenate([G@uu, G@vd]) for G in six]).T, tol=1e-9)))
ok("6. the six plane turns (three runs, three rotations), 300 random states: they move one list in %s independent ways - two combinations leave it "
   "unchanged; a pair with a rest frame in %s; a pair without one, v = (alpha + beta w) u, in %s again" % (sorted(r1), sorted(r2), sorted(r3)),
   r1 == {4} and r2 == {6} and r3 == {4} and smin > 1e-6)

# ---- 7. the frame reading
even = [np.eye(4)] + six[3:] + six[:3] + [wn]                                  # 1 | yz zx xy | x y z | w
u0 = np.array([1., 0, 0, 0]); v0 = Tn@u0; r0 = np.concatenate([u0, v0])
M8 = np.array([np.concatenate([M@u0, M@v0]) for M in even]).T; M4 = np.array([M@u0 for M in even[:4]]).T
closed = all(np.linalg.matrix_rank(np.vstack([np.array([m.ravel() for m in even]), (p@q_).ravel()])) == 8 for p in even for q_ in even)
ok("7. the frame reading. The eight matrices [1 | rotations yz, zx, xy | runs x, y, z | w] are closed under products; applied to the reference pair "
   "(u0, T u0), u0 = the forward counter, they give every pair once: det = %.0f. One list: (1 + three rotations) u0 gives every list once: det = %.0f. "
   "Eight counters = a size, six plane turns and one more; four counters = a size and three space-plane turns"
   % (np.linalg.det(M8), np.linalg.det(M4)), closed and abs(np.linalg.det(M8)) > 1 and abs(np.linalg.det(M4)) > 0.5)

# ---- 8. what the step never reads
gens8 = A8 + [T8]
rows = []
for i in range(8):
    for j in range(8):
        E_ = np.zeros((8, 8)); E_[i, j] = 1.0
        rows.append(np.concatenate([(E_@G - G@E_).ravel() for G in gens8]))
comm_dim = 64 - np.linalg.matrix_rank(np.array(rows).T, tol=1e-9)
K8 = np.diag([1., 1, 1, 1, -1, -1, -1, -1])
in_comm = all(np.allclose(X@G, G@X) for X in (J8, K8, J8@K8) for G in gens8)
# weight-keeping (orthogonal) with square -1 inside {a + b J + c K + d JK}: on the list index these are 2 x 2 matrices [[a + c, -b + d], [b + d, a - c]]
aa, bb, cc, dd = sp.symbols('aa bb cc dd', real=True)
m2 = sp.Matrix([[aa + cc, -bb + dd], [bb + dd, aa - cc]])
sols = sp.solve(list(m2.T*m2 - sp.eye(2)) + list(m2*m2 + sp.eye(2)), [aa, bb, cc, dd], dict=True)
coef = np.linalg.solve(M8, J8@r0)                                              # J on the reference = which even element?
B = sum(c_*M for c_, M in zip(coef, even))
spin = [float(r0@(-J8@np.kron(np.eye(2), R))@r0) for R in (six[3], six[4], six[5])]   # Sigma_x, Sigma_y, Sigma_z on the reference (VI.7 (8))
lab = ['1', 'rotation yz', 'rotation zx', 'rotation xy', 'run x', 'run y', 'run z', 'w']
which = [lab[i] for i in range(8) if abs(coef[i]) > 1e-12]
ok("8. what the step never reads on the pair: the matrices commuting with all of A_x, A_y, A_z, T on both lists are %d, spanned by 1, J, K, JK; those "
   "that keep the weight and square to -1 are %s, i.e. +-J only. In the frame reading, J (reference) = B (reference) with B = %s x %s, B^2 = -1: "
   "the reference turned in one of its own space planes. The reference's spin is along (%.0f, %.0f, %.0f): the plane is the one across the spin axis"
   % (comm_dim, sols, which, [round(float(c_), 12) for c_ in coef if abs(c_) > 1e-12], *spin),
   comm_dim == 4 and in_comm and len(sols) == 2 and all(s_[aa] == 0 and s_[cc] == 0 and s_[dd] == 0 and abs(s_[bb]) == 1 for s_ in sols)
   and len(which) == 1 and which[0].startswith('rotation') and np.allclose(B@B, -np.eye(4))
   and abs(spin[['rotation yz', 'rotation zx', 'rotation xy'].index(which[0])]) > 0.5 and sum(abs(x) > 1e-12 for x in spin) == 1)

# ---- 9. in motion
import scipy.linalg as sl
dev = 0.0; steady_move = 0.0; found = None
for eta in (0.3, 0.8, 1.7):
    Brun = sl.expm(eta*n(Az)/2); a0 = 1.0; ww = a0*np.cosh(eta)
    for sgn in (1, -1):
        kk = sgn*a0*np.sinh(eta); G = -kk*A8[2]@J8 + a0*T8; Pp = (np.eye(8) + J8@G/ww)/2
        worst = 0.0
        for _ in range(20):
            uu = rng.normal(size=4); rr = np.concatenate([Brun@uu, Brun@Tn@uu]); worst = max(worst, np.abs(Pp@rr - rr).max())
        if worst < 1e-12: found = sgn; dev = max(dev, worst)
    uu = rng.normal(size=4); fl = []
    for th in np.linspace(0, 2*np.pi, 9):
        Eth = np.cos(th)*np.eye(4) + np.sin(th)*Tn; rr = np.concatenate([Brun@Eth@uu, Brun@Eth@Tn@uu]); fl.append([rr@A@rr for A in A8] + [rr@rr])
    steady_move = max(steady_move, np.ptp(np.array(fl), axis=0).max())
ok("9. in motion: the rest pair (u, T u), run to rapidity eta along z (e^{eta A_z/2} on each list), lies in the leading class at ck = %+d a sinh(eta), "
   "w = a cosh(eta) (%.0e, eta = 0.3, 0.8, 1.7); its density and flow do not change as its clock turns (%.0e). The steady moving walker's second "
   "column is the class partner; v(0) = T u(0) is its form at rest" % (found, dev, steady_move), found is not None and dev < 1e-12 and steady_move < 1e-12)

print("ALL PASS" if ok_all else "SOME FAILED"); raise SystemExit(0 if ok_all else 1)
