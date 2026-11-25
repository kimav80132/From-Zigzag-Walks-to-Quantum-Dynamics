"""A skewed noise covariance C = eta + kappa S (S symmetric, off-diagonal: the ellipsoid/hyperboloid tilted), det kept.
1. In the Fokker-Planck equation it enters as (1/2) C^{mu nu} d_mu d_nu: a change of the SECOND-order operator (the metric), no
   first-derivative A.d term, no A^2 term, no Phi term. A Lorentz Lambda gives Lambda eta Lambda^T = eta: no skew. A symmetric C of the
   same signature is Lambda' eta Lambda'^T for a linear Lambda' that is not Lorentz: a change of coordinates, not a moving observer (V.3).
2. Could a position-dependent skew C(x) give a first-order term?  In Ito form d_mu d_nu (C^{mu nu} rho) it gives the
   'spurious drift' d_nu C^{mu nu}: real, and a gradient-type term; compare with the Lorentz force e F^mu_nu v^nu: that needs
   the velocity v, which a drift from C(x) does not contain -> cannot be the Lorentz force for a general state.
3. What the field DOES do with a symmetric-skew structure: the two-step moment of the Dirac process gamma^mu gamma^nu has
   an antisymmetric part; the field couples to it through [D_mu, D_nu]. Check: a boost of the step matrices
   gamma' = Lambda gamma (vector index) is equivalent to the spinor change S gamma S^-1: the Dirac process is boost-covariant,
   the 'skew' of the matrix noise IS a boost, and it leaves the equation form-invariant (again: a change of observer, not a field).
"""
import sympy as sp
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
eta = sp.diag(-1, 1); th = sp.symbols('theta', real=True)
Lam = sp.Matrix([[sp.cosh(th), sp.sinh(th)], [sp.sinh(th), sp.cosh(th)]])
C = sp.simplify(Lam*eta*Lam.T)
ok(f"1a. a Lorentz Lambda does not skew the covariance at all: Lambda eta Lambda^T = {C.tolist()} = eta (a moving observer sees the same noise matrix)", C == eta)
# 1b (rebuilt 10 Nov 2026: the earlier check only confirmed det = -1 of a matrix built to have it, and 1a's label called the result "off-diagonal").
#    A symmetric C of signature (-,+) is Lambda' eta Lambda'^T for a real linear Lambda' that is NOT Lorentz: a change of coordinates.
import numpy as np
rng = np.random.default_rng(4); okb = True; nonlor = 0
for _ in range(200):
    M_ = rng.normal(size=(2, 2)); Cn = M_ @ np.diag([-1., 1.]) @ M_.T                 # every such C has signature (-,+); recover a Lambda' from C alone
    w_, V_ = np.linalg.eigh(Cn)                                                        # one negative, one positive eigenvalue
    Lp = V_ @ np.diag(np.sqrt(np.abs(w_)))                                             # C = Lp diag(sign w) Lp^T, and eigh orders the negative one first
    okb &= (w_[0] < 0 < w_[1]) and np.allclose(Lp @ np.diag([-1., 1.]) @ Lp.T, Cn)
    nonlor += not np.allclose(Cn, np.diag([-1., 1.]))
ok(f"1b. every symmetric C with one negative and one positive eigenvalue is Lambda' eta Lambda'^T for a real Lambda' built from C (200 random C; {nonlor} of them differ from eta, "
   "so their Lambda' is not a Lorentz transformation): a change of coordinates", okb and nonlor == 200)
s, d1 = sp.symbols('s d1', real=True)
Cs = sp.Matrix([[-(1 + d1), s], [s, (1 - s**2)/(1 + d1)]])
t, x = sp.symbols('t x', real=True); f = sp.Function('f')(t, x)
X = (t, x)
op = lambda M: sum(M[i, j]*sp.diff(f, X[i], X[j]) for i in range(2) for j in range(2))/2
ok("1c. the FP operator with the skewed C has only second derivatives: no term A.d f, no A^2 f, no Phi f (count of first-order terms = 0)",
   all(sp.expand(op(Cs)).coeff(sp.Derivative(f, v)) == 0 for v in X) and sp.expand(op(Cs)).coeff(f) == 0)
# 2 position-dependent skew, Ito form
sx = sp.Function('s', real=True)(t, x); d1x = sp.Function('d', real=True)(t, x)
Cx = sp.Matrix([[-(1 + d1x), sx], [sx, (1 - sx**2)/(1 + d1x)]]); rho = sp.Function('rho', positive=True)(t, x)
ito = sum(sp.diff(Cx[i, j]*rho, X[i], X[j]) for i in range(2) for j in range(2))/2
drift = [sum(sp.diff(Cx[i, j], X[j]) for j in range(2)) for i in range(2)]      # spurious drift d_nu C^{mu nu}
ok("2. Ito form: d_mu d_nu(C rho)/2 = C d d rho/2 + (d_nu C^{mu nu}) d_mu rho + ...: the induced drift is d_nu C^{mu nu}, a function of position only, "
   "no dependence on the velocity v -> it cannot be e F^mu_nu v^nu for all states",
   sp.simplify(sp.expand(ito - (sum(Cx[i, j]*sp.diff(rho, X[i], X[j]) for i in range(2) for j in range(2))/2
                               + sum(drift[i]*sp.diff(rho, X[i]) for i in range(2)) + rho*sum(sp.diff(Cx[i, j], X[i], X[j]) for i in range(2) for j in range(2))/2))) == 0)
# 3 matrix noise: a boost of the step matrices is a spinor change
A = sp.diag(1, -1); T = sp.Matrix([[0, 1], [-1, 0]]); beta = sp.I*T
g0, g1 = beta, beta*A                                   # gamma^0, gamma^1 ; (+,-) signature for these
Lb = sp.Matrix([[sp.cosh(th), sp.sinh(th)], [sp.sinh(th), sp.cosh(th)]])
gp = [Lb[0, 0]*g0 + Lb[0, 1]*g1, Lb[1, 0]*g0 + Lb[1, 1]*g1]            # gamma'^mu = Lambda^mu_nu gamma^nu
Sp = sp.cosh(th/2)*sp.eye(2) + sp.sinh(th/2)*A                          # spinor boost exp(theta A/2)
ok("3. gamma'^mu = Lambda^mu_nu gamma^nu equals S^-1 gamma^mu S with S = exp(theta A/2): the skewed matrix noise is a boost of the column, "
   "the equation keeps its form (an observer change, not a field)",
   all(sp.simplify(gp[mu] - Sp.inv()*[g0, g1][mu]*Sp) == sp.zeros(2) for mu in range(2)))
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
