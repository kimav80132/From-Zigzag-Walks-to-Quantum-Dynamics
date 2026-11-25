"""The other AI's proposal: the field acts on the noise vector as  xi' = (1 + (e/m) F dtau) xi.
1. Lambda = 1 + eps F^mu_nu is an infinitesimal Lorentz transformation: Lambda eta Lambda^T = eta + O(eps^2), det = 1 + O(eps^2).
2. The covariance of the rotated noise:  <xi' xi'> = Lambda eta Lambda^T = eta + eps (F eta + eta F^T) = eta exactly to first order
   (F^mu_nu eta^{nu rho} = F^{mu rho} antisymmetric, so the first-order change is F^{mu rho} + F^{rho mu} = 0).
   So the rotation does NOTHING to the noise matrix at first order: the Fokker-Planck operator is unchanged, no field in the equation.
3. Applied to the mean step instead (the drift, first moment): xi includes the drift b^mu; (1 + eps F) b = b + eps F b: that IS the Lorentz force,
   dv^mu/dtau = (e/m) F^mu_nu v^nu, a change of the first moment.
4. Minimal coupling in the step: the per-step phase e^{i e z.A/hbar}, through the covariance, gives the drift shift (i e/m) A (imaginary_skew.py).
   Its curl is the Lorentz force on the mean velocity: the Madelung mean velocity v = (dS - eA)/m obeys m dv/dtau = e F v on the classical (hbar -> 0) solutions.
"""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
eps, Ex, Ey, Ez, Bx, By, Bz = sp.symbols('epsilon E_x E_y E_z B_x B_y B_z', real=True)
eta = sp.diag(-1, 1, 1, 1)
Fup = sp.Matrix([[0, Ex, Ey, Ez], [-Ex, 0, -Bz, By], [-Ey, Bz, 0, -Bx], [-Ez, -By, Bx, 0]])     # F^{mu nu} (c = 1)
Fmix = Fup*eta                                                                                   # F^mu_nu = F^{mu rho} eta_{rho nu}
Lam = sp.eye(4) + eps*Fmix
ok("1. Lambda eta Lambda^T - eta = O(eps^2) and det Lambda = 1 + O(eps^2)",
   sp.expand(Lam*eta*Lam.T - eta).applyfunc(lambda f: f.coeff(eps, 1)) == sp.zeros(4) and sp.expand(Lam.det()).coeff(eps, 1) == 0)
cov = sp.expand(Lam*eta*Lam.T)
ok("2. covariance of the rotated noise to first order = eta: the field drops out of the noise matrix", cov.applyfunc(lambda f: f.coeff(eps, 1)) == sp.zeros(4))
print("   (the second-order piece is eps^2 F eta F^T, symmetric: it changes the diffusion constant by O(F^2), still no first-order coupling)")
u = sp.Matrix(sp.symbols('u0:4', real=True))
ok("3. on the mean step: (1 + eps F) u - u = eps F^mu_nu u^nu, the Lorentz force on the drift (first moment)", sp.simplify(Lam*u - u - eps*Fmix*u) == sp.zeros(4, 1))
# 4. Madelung velocity with minimal coupling obeys the Lorentz force (classical solutions)
t, x = sp.symbols('t x', real=True); m, e = sp.symbols('m e', positive=True)
S = sp.Function('S', real=True)(t, x); Phi, Ax = sp.Function('Phi', real=True)(t, x), sp.Function('Ax', real=True)(t, x)
X = (t, x); A_low = (-Phi, Ax)                       # A_mu, eta = diag(-1,1), c = 1
pi_low = [sp.diff(S, X[mu]) - e*A_low[mu] for mu in range(2)]              # m v_mu
pi_up = [eta[mu, mu]*pi_low[mu] for mu in range(2)]
HJ = sum(pi_up[mu]*pi_low[mu] for mu in range(2)) + m**2                   # = 0 on classical solutions
# d/dtau (m v_mu) = v^nu d_nu (m v_mu) = v^nu d_nu(d_mu S - e A_mu) = v^nu d_mu d_nu S - e v^nu d_nu A_mu
#   = v^nu d_mu (m v_nu + e A_nu) - e v^nu d_nu A_mu = (1/2) d_mu (m v.v) + e v^nu (d_mu A_nu - d_nu A_mu) = (1/2) d_mu(-m) ... + e F_{mu nu} v^nu
Flow = sp.Matrix(2, 2, lambda mu, nu: sp.diff(A_low[nu], X[mu]) - sp.diff(A_low[mu], X[nu]))
v_up = [pu/m for pu in pi_up]
lhs = [sum(v_up[nu]*sp.diff(pi_low[mu], X[nu]) for nu in range(2)) for mu in range(2)]
rhs = [sp.diff(HJ, X[mu])/(2*m) + e*sum(Flow[mu, nu]*v_up[nu] for nu in range(2)) for mu in range(2)]
ok("4. m dv_mu/dtau = (1/2m) d_mu(m^2 v.v + m^2) + e F_{mu nu} v^nu: on the mass shell the first term vanishes, leaving the Lorentz force",
   all(sp.simplify(sp.expand(lhs[mu] - rhs[mu])) == 0 for mu in range(2)))
