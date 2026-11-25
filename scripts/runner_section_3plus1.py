"""Checks for the draft section 'IV.4c' (runner step first), in the master's conventions (IV.4b):
   c d_z Q = -A_z d_t Q + a (A_z T) Q + c W_x d_x Q + c W_y d_y Q,   Q = e^{az/c} P,   A_z = diag(eps) = diag(1,-1,-1,1).
(a) runner along x, first order: d_t = d_t' + v d_x',  d_x = d_x' + (v/c^2) d_t'
(b) the runner's lines have the same form  <=>  W_x = K A_z - A_z K,  K W_x - W_x K = A_z,  K (A_z T) = (A_z T) K,  K W_y = W_y K
(c) entry by entry W_ij = K_ij (eps_j - eps_i): zero diagonal
(d) two counters, one runner: K = +-(1/2)[[0,1],[1,0]] (+ k 1), W = -+[[0,-1],[1,0]] ; second runner: no solution
(e) 2+1 walk check: the step with coins from W reproduces the 2+1 P lines (done in walk_1d_to_2d_by_runs.py, W8)
(f) four counters: every solution of both runners' conditions (random starts) has W_x^2 = W_y^2 = -1, W_x W_y = -W_y W_x,
    and W_y W_x = 2 (K_x K_y - K_y K_x)
(g) K_x = A_x/2, K_y = A_y/2 (IV.4b (I)) satisfy everything, and W_x, W_y are exactly the x-turn and y-turn of IV.4b (A)
(h) the four lines (1)-(4) of IV.4b (A) are reproduced
(i) directions: [turn zx, turn zy] = -turn xy; small loop of the two turns = small xy turn; counters obey the same relation
(j) x then y = W_y W_x = A_x A_y = 2 x (xy turn), y then x = -A_x A_y; (A_x A_y)^2 = -1; full xy turn on counters = -1
(k) walker from counter 4: x then y -> 3 -> 1 with +1; y then x -> 2 -> 1 with -1
(l) finite runs along x and along y leave the four lines unchanged
"""
import numpy as np, sympy as sp
from scipy.optimize import least_squares
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
c, a = sp.symbols('c a', positive=True); eps, eta = sp.symbols('epsilon eta', real=True)
Az = sp.diag(1, -1, -1, 1)
Ax = sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]])
Ay = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
T = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]])
Ex = Az*T
# (a)
t, x = sp.symbols('t x'); v = sp.symbols('v')
tp, xp = t + v*x/c**2, x + v*t
J = sp.Matrix([[tp.diff(t), xp.diff(t)], [tp.diff(x), xp.diff(x)]])     # d_t = d tp/dt d_t' + d xp/dt d_x', ...
ok("(a) d_t = d_t' + v d_x',  d_x = (v/c^2) d_t' + d_x'", J == sp.Matrix([[1, v], [v/c**2, 1]]))
# (b) generic 4x4 K, W_x, W_y: first-order residual of the runner's lines
def gen(n, name): return sp.Matrix(n, n, sp.symbols(f'{name}0:{n*n}'))
K, Wx, Wy = gen(4, 'k'), gen(4, 'wx'), gen(4, 'wy')
q, qt, qx, qy = [sp.Matrix(sp.symbols(f'{s}1:5')) for s in ('q', 'qt', 'qx', 'qy')]
line = lambda Q, Qt, Qx, Qy: (-Az*Qt + a*Ex*Q + c*Wx*Qx + c*Wy*Qy)/c       # d_z Q
S = sp.eye(4) + eps*K                                                        # eps = v/c
Qp, Qtp, Qxp, Qyp = S*q, S*(qt - eps*c*qx), S*(qx - eps/c*qt), S*qy          # d_t' = d_t - v d_x, d_x' = d_x - (v/c^2) d_t (first order)
res = (S*line(q, qt, qx, qy) - line(Qp, Qtp, Qxp, Qyp)).applyfunc(lambda e: sp.expand(sp.diff(e, eps).subs(eps, 0)))
grab = lambda vec: sp.Matrix(4, 4, lambda i, j: res[i].coeff(vec[j]))
Ct, Cx, Cy, C0 = grab(qt), grab(qx), grab(qy), grab(q)
ok("(b) first-order runner residual: d_t part = (W_x - (K A_z - A_z K))/c, d_x part = K W_x - W_x K - A_z, d_y part = K W_y - W_y K, flips = a(K E - E K)/c",
   sp.expand(Ct*c - (Wx - (K*Az - Az*K))) == sp.zeros(4) and sp.expand(Cx - (K*Wx - Wx*K - Az)) == sp.zeros(4)
   and sp.expand(Cy - (K*Wy - Wy*K)) == sp.zeros(4) and sp.expand(C0*c/a - (K*Ex - Ex*K)) == sp.zeros(4))
# (c)
epsv = [1, -1, -1, 1]
ok("(c) (K A_z - A_z K)_ij = K_ij (eps_j - eps_i): diagonal always zero",
   all(sp.expand((K*Az - Az*K)[i, j] - K[i, j]*(epsv[j] - epsv[i])) == 0 for i in range(4) for j in range(4)))
# (d) two counters
A2, T2 = sp.diag(1, -1), sp.Matrix([[0, 1], [-1, 0]]); E2 = A2*T2
k2 = gen(2, 'm'); W2 = k2*A2 - A2*k2
sol = sp.solve(list(k2*W2 - W2*k2 - A2) + list(k2*E2 - E2*k2), list(k2), dict=True)
Ws = [W2.subs(s_) for s_ in sol]
ok(f"(d1) one runner in two counters: K = {[k2.subs(s_).tolist() for s_ in sol]}, W = {[w.tolist() for w in Ws]}",
   len(sol) == 2 and all(w == sp.Matrix([[0, -1], [1, 0]]) or w == sp.Matrix([[0, 1], [-1, 0]]) for w in Ws))
Wx2 = sp.Matrix([[0, -1], [1, 0]]); ky = gen(2, 'n'); Wy2 = ky*A2 - A2*ky
sol2 = sp.solve(list(ky*Wy2 - Wy2*ky - A2) + list(ky*E2 - E2*ky) + list(ky*Wx2 - Wx2*ky), list(ky), dict=True)
ok(f"(d2) a second runner (along y) in two counters: solutions {sol2}", sol2 == [])
# (f) four counters, all solutions (numerical, random starts)
Azn, Exn = np.array(Az, float), np.array(Ex, float)
def conds(vv):
    Kx_, Ky_ = vv[:16].reshape(4, 4), vv[16:].reshape(4, 4); Wx_, Wy_ = Kx_@Azn - Azn@Kx_, Ky_@Azn - Azn@Ky_
    r = [Kx_@Wx_ - Wx_@Kx_ - Azn, Kx_@Exn - Exn@Kx_, Kx_@Wy_ - Wy_@Kx_, Ky_@Wy_ - Wy_@Ky_ - Azn, Ky_@Exn - Exn@Ky_, Ky_@Wx_ - Wx_@Ky_]
    return np.concatenate([m_.ravel() for m_ in r])
rng = np.random.default_rng(11); worst = np.zeros(4); n = 0
for _ in range(300):
    s_ = least_squares(conds, rng.normal(size=32), xtol=1e-14, ftol=1e-14, gtol=1e-14)
    if np.linalg.norm(conds(s_.x)) < 1e-9:
        n += 1; Kx_, Ky_ = s_.x[:16].reshape(4, 4), s_.x[16:].reshape(4, 4); Wx_, Wy_ = Kx_@Azn - Azn@Kx_, Ky_@Azn - Azn@Ky_
        worst = np.maximum(worst, [np.abs(Wx_@Wx_ + np.eye(4)).max(), np.abs(Wy_@Wy_ + np.eye(4)).max(),
                                   np.abs(Wx_@Wy_ + Wy_@Wx_).max(), np.abs(Wy_@Wx_ - 2*(Kx_@Ky_ - Ky_@Kx_)).max()])
ok(f"(f) {n} solutions in four counters: max |W_x^2+1|, |W_y^2+1|, |W_xW_y+W_yW_x|, |W_yW_x - 2[K_x,K_y]| = {np.round(worst, 12)}",
   n > 50 and worst.max() < 1e-8)
# (g)
Kx, Ky = Ax/2, Ay/2; WX, WY = Kx*Az - Az*Kx, Ky*Az - Az*Ky
xturn = sp.Matrix([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]])
yturn = sp.Matrix([[0, 0, -1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, -1, 0, 0]])
ok("(g) K_x = A_x/2, K_y = A_y/2: W_x = A_x A_z = x-turn of IV.4b (A), W_y = A_y A_z = y-turn, all conditions hold",
   WX == xturn and WY == yturn and WX == Ax*Az and WY == Ay*Az and Kx*WX - WX*Kx == Az and Ky*WY - WY*Ky == Az
   and Kx*Ex == Ex*Kx and Ky*Ex == Ex*Ky and Kx*WY == WY*Kx and Ky*WX == WX*Ky)
# (h)
tt, xx, yy, zz = sp.symbols('t x y z'); P = sp.Matrix([sp.Function(f'P{i}')(tt, xx, yy, zz) for i in range(1, 5)])
rhs = -Az*P.diff(tt) - a*P + a*Ex*P + c*WX*P.diff(xx) + c*WY*P.diff(yy)
doc = sp.Matrix([-P[0].diff(tt) - a*P[0] + a*P[2] - c*P[1].diff(xx) - c*P[2].diff(yy),
                 P[1].diff(tt) - a*P[1] - a*P[3] + c*P[0].diff(xx) + c*P[3].diff(yy),
                 P[2].diff(tt) - a*P[2] + a*P[0] - c*P[3].diff(xx) + c*P[0].diff(yy),
                 -P[3].diff(tt) - a*P[3] - a*P[1] + c*P[2].diff(xx) - c*P[1].diff(yy)])
ok("(h) the four lines (1)-(4) of IV.4b (A) are reproduced", sp.expand(rhs - doc) == sp.zeros(4, 1))
# (i)
def E_(i, j, n=4):
    m_ = sp.zeros(n); m_[i, j] = 1; return m_
B = {k_: E_(0, i) + E_(i, 0) for k_, i in (('x', 1), ('y', 2), ('z', 3))}
com = lambda P_, Q_: P_*Q_ - Q_*P_
Rzx, Rzy, Rxy = com(B['z'], B['x']), com(B['z'], B['y']), com(B['x'], B['y'])
e1, e2 = sp.symbols('e1 e2'); ex2 = lambda M_, s: sp.eye(4) + s*M_ + s**2*M_*M_/2
loop = (ex2(Rzx, e1)*ex2(Rzy, e2)*ex2(Rzx, -e1)*ex2(Rzy, -e2)).applyfunc(sp.expand)
cut = loop.applyfunc(lambda f: sum(tm for tm in sp.Add.make_args(f) if sp.Poly(tm, e1, e2).total_degree() <= 2))
Kz = Az/2
ok("(i) directions: [turn zx, turn zy] = -turn xy, the small loop of the two turns is 1 - e1 e2 turn xy; counters: same relation",
   com(Rzx, Rzy) == -Rxy and cut == sp.eye(4) - e1*e2*Rxy and com(com(Kz, Kx), com(Kz, Ky)) == -com(Kx, Ky))
# (j)
th = sp.symbols('theta', real=True); Jxy = com(Kx, Ky)
ok("(j) x then y = W_y W_x = A_x A_y = 2 turn(xy);  y then x = W_x W_y = -A_x A_y;  (A_x A_y)^2 = -1;  exp(2 pi turn xy) = -1",
   WY*WX == Ax*Ay and WX*WY == -Ax*Ay and WY*WX == 2*Jxy and (Ax*Ay)**2 == -sp.eye(4)
   and sp.simplify((2*sp.pi*Jxy).exp() + sp.eye(4)) == sp.zeros(4))
ok("(j2) A_x A_y written out", Ax*Ay == sp.Matrix([[0, 0, 0, 1], [0, 0, 1, 0], [0, -1, 0, 0], [-1, 0, 0, 0]]))
# (k)
ok("(k) from counter 4: x-turn -> 3 (entry -1), y-turn -> 1 (entry -1): +1;  y-turn -> 2 (entry +1), x-turn -> 1 (entry -1): -1",
   WX[2, 3] == -1 and WY[0, 2] == -1 and (WY*WX)[0, 3] == 1 and WY[1, 3] == 1 and WX[0, 1] == -1 and (WX*WY)[0, 3] == -1)
# (l) finite runs: coordinates (ct, x, y, z), counters exp(eta A/2) = cosh(eta/2) + sinh(eta/2) A
u = sp.Matrix(sp.symbols('u1:5')); ut, ux, uy = [sp.Matrix(sp.symbols(f'{s}1:5')) for s in ('ut', 'ux', 'uy')]
uz = (-Az*ut + a*Ex*u + c*WX*ux + c*WY*uy)/c                                 # loss taken out
grad = [ut/c, ux, uy, uz]
def finite(Amat, i):
    L = sp.eye(4); L[0, 0] = L[i, i] = sp.cosh(eta); L[0, i] = L[i, 0] = sp.sinh(eta); Li = L.inv()
    Sm = sp.cosh(eta/2)*sp.eye(4) + sp.sinh(eta/2)*Amat
    d = [sum((Li[nu, mu]*grad[nu] for nu in range(4)), sp.zeros(4, 1)) for mu in range(4)]
    dQ = [Sm*e for e in d]
    E = dQ[3] - (-Az*dQ[0] + a/c*Ex*(Sm*u) + WX*dQ[1] + WY*dQ[2])        # d/d(ct) = d_t / c
    return E.applyfunc(lambda e: sp.simplify(sp.expand(sp.expand(e).rewrite(sp.exp))))
ok("(l) finite run along x (counters cosh + sinh A_x) and along y (cosh + sinh A_y) leave the four lines unchanged",
   finite(Ax, 1) == sp.zeros(4, 1) and finite(Ay, 2) == sp.zeros(4, 1))
