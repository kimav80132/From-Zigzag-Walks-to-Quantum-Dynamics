"""VI.2 (E)-(F): the carrier of the pair-angle field, and what its own law can be.
(E) The carrier. The Casimir document: "the long-range force needs a carrier that never flips". Part II's walker with m = 0 never flips (K = 0 in II.7)
    and has no direction (VI.1 (F): a boson). Its source answer is III.16 at m = 0:
    1. 3+1: int_0^inf dtau/(4 pi tau)^2 exp(-s^2/4 tau) = 1/(4 pi^2 s^2), the m -> 0 limit of (m/4 pi^2 s) K_1(m s).
    2. a static source: integrate over Euclidean time, int dt 1/(4 pi^2 (r^2 + t^2)) = 1/(4 pi r): Coulomb's law from the massless undirected walker.
    3. 1+1: the same at m = 0 is -(1/2 pi) ln s up to a constant; its static field is -sign(x)/2 (a potential |x|/2, confining in one dimension).
    4. 1/(4 pi^2 s^2) is the Fourier transform of 1/k^2: the carrier's kernel IS the propagator of the F^2 action in Feynman gauge, eta_{mu nu}/k^2.
(F) The field's own law. The pair angle per link is a one-form a_mu (the connection of pair_angle_is_a_connection.py). If it is itself a noise of the
    event (the event's clock-zero drawn at random), its law is a local quadratic weight. Two invariances are already in the document: no preferred
    frame (II.1) and the relabelling a -> a + d chi (V.2 (G), derived in VI.2 (A)-(B)). With two derivatives there are three invariant scalars,
    T1 = d_mu a_nu d^mu a^nu, T2 = d_mu a_nu d^nu a^mu, T3 = (d.a)^2, and one mass scalar T0 = a.a.
    5. the relabelling invariance of the equations of motion forces c0 = 0 and c1 + c2 + c3 = 0; T2 and T3 have the same equations of motion
       (they differ by a total derivative). So the weight is unique up to scale: c1 (T1 - T2) = (c1/2) F_{mu nu} F^{mu nu}. Maxwell's free action, forced.
    6. its equations of motion are d_mu F^{mu nu} = 0, and with the walker's turn angle as the source, d_mu F^{mu nu} = j^nu with j^nu = e c psi^dagger alpha^nu psi:
       the source of the field is the walker's own current, the derivative of the walker's energy with respect to a_nu (V.2 (F)).
    What is NOT derived: that the clock-zero is a noise with a local quadratic weight (the statement that the field has its own dynamics), and the scale (e)."""
import sympy as sp, numpy as np
from scipy import integrate, special
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
# ---------- (E) ----------
tau, s, r, t, m = sp.symbols('tau s r t m', positive=True)
K31 = sp.integrate(sp.exp(-s**2/(4*tau))/(4*sp.pi*tau)**2, (tau, 0, sp.oo))
ok("1. 3+1 at m = 0: int dtau (4 pi tau)^-2 exp(-s^2/4tau) = 1/(4 pi^2 s^2); and the massive (m/4 pi^2 s) K_1(ms) tends to it as m -> 0",
   sp.simplify(K31 - 1/(4*sp.pi**2*s**2)) == 0 and abs(float(1e-6/(4*np.pi**2*0.7)*special.kv(1, 1e-6*0.7)) - 1/(4*np.pi**2*0.49)) < 1e-9)
static = sp.integrate(1/(4*sp.pi**2*(r**2 + t**2)), (t, -sp.oo, sp.oo))
ok("2. a static source: int dt 1/(4 pi^2 (r^2 + t^2)) = 1/(4 pi r): Coulomb's law, from Part II's walker at m = 0", sp.simplify(static - 1/(4*sp.pi*r)) == 0)
# 1+1: the tau integral diverges at large tau (the log); the difference between two points is finite
x, L = sp.symbols('x L', positive=True)
# with u = 1/tau the cut integral is E_1(s^2/4L)/(4 pi), the exponential integral; its small-argument form is -gamma - ln z
z = sp.symbols('z', positive=True)
E1_small = sp.series(sp.expint(1, z), z, 0, 2).removeO()
diff11 = sp.limit((E1_small.subs(z, s**2/(4*L)) - E1_small.subs(z, 1/(4*L)))/(4*sp.pi), L, sp.oo)
num = integrate.quad(lambda tt: np.exp(-2.3**2/(4*tt))/(4*np.pi*tt), 0, 1e6)[0] - integrate.quad(lambda tt: np.exp(-1/(4*tt))/(4*np.pi*tt), 0, 1e6)[0]
x_ = 1.3
ok("3. 1+1 at m = 0: G(s) - G(1) = -(1/2 pi) ln s (numerically at s = 2.3 to 1e-6); the static field int dt d_x G = -sign(x)/2, a potential |x|/2 (confining in one dimension)",
   sp.simplify(diff11 + sp.log(s)/(2*sp.pi)) == 0 and abs(num + np.log(2.3)/(2*np.pi)) < 1e-6 and
   abs(integrate.quad(lambda tt: -x_/(2*np.pi*(x_**2 + tt**2)), -np.inf, np.inf)[0] + 0.5) < 1e-9)
# 4. Fourier transform: int d^4k/(2pi)^4 e^{ik.x}/k^2 = 1/(4 pi^2 x^2) in 4 Euclidean dimensions. Check numerically via the radial form:
# (1/(2pi)^4) * (2 pi)^2 ... use the known formula int d^dk e^{ikx} k^{-2} /(2pi)^d = Gamma(d/2-1)/(4 pi^{d/2} x^{d-2}); d = 4: 1/(4 pi^2 x^2)
d = 4
ok("4. 1/(4 pi^2 s^2) is the transform of 1/k^2 in four dimensions (Gamma(d/2 - 1)/(4 pi^{d/2} s^{d-2}) at d = 4): the carrier's kernel is the Feynman-gauge propagator of F^2",
   sp.simplify(sp.gamma(sp.Rational(d, 2) - 1)/(4*sp.pi**sp.Rational(d, 2)*s**(d - 2)) - 1/(4*sp.pi**2*s**2)) == 0)
# ---------- (F) ----------
X = sp.symbols('x0:4', real=True); eta = sp.diag(-1, 1, 1, 1)
a = [sp.Function(f'a{mu}', real=True)(*X) for mu in range(4)]
chi = sp.Function('chi', real=True)(*X)
d = lambda f, mu: sp.diff(f, X[mu])
up = lambda v: [sum(eta[mu, nu]*v[nu] for nu in range(4)) for mu in range(4)]
aU = up(a)
T0 = sum(a[mu]*aU[mu] for mu in range(4))
T1 = sum(eta[mu, mu]*eta[nu, nu]*d(a[nu], mu)**2 for mu in range(4) for nu in range(4))                       # d_mu a_nu d^mu a^nu
T2 = sum(eta[mu, mu]*eta[nu, nu]*d(a[nu], mu)*d(a[mu], nu) for mu in range(4) for nu in range(4))              # d_mu a_nu d^nu a^mu
T3 = sum(eta[mu, mu]*d(a[mu], mu) for mu in range(4))**2                                                       # (d.a)^2
c0, c1, c2, c3 = sp.symbols('c0:4', real=True)
Lag = c0*T0 + c1*T1 + c2*T2 + c3*T3
from sympy.calculus.euler import euler_equations
EL = [eq.lhs for eq in euler_equations(Lag, a, X)]
# relabelling a -> a + d chi: the equations of motion must not change (chi drops out)
sub = {a[mu]: a[mu] + d(chi, mu) for mu in range(4)}
shift = [sp.expand((E.subs(sub) - E).doit()) for E in EL]
# relabelling: put chi = exp(i k.x); the shift of each equation is a polynomial in k times the exponential, which must vanish identically in k
k = sp.symbols('k0:4', real=True)
pw = sp.exp(sp.I*sum(k[mu]*X[mu] for mu in range(4)))
conds = set()
for e_ in shift:
    poly = sp.expand(sp.simplify(e_.subs(chi, pw).doit()/pw))
    conds |= set(sp.Poly(poly, *k).coeffs())
sol = sp.solve(list(conds), [c0, c1, c2, c3], dict=True)
ok(f"5. relabelling invariance of the equations of motion gives {sol}: c0 = 0 and c1 + c2 + c3 = 0", sol == [{c0: 0, c1: -c2 - c3}])
EL2 = [sp.expand(eq.lhs) for eq in euler_equations(T2, a, X)]; EL3 = [sp.expand(eq.lhs) for eq in euler_equations(T3, a, X)]
ok("5a. T2 and T3 have the same equations of motion (a total derivative apart), so the weight is c1 (T1 - T2) up to total derivatives, unique up to scale",
   all(sp.simplify(EL2[mu] - EL3[mu]) == 0 for mu in range(4)))
kk_ = sp.symbols('q0:4', real=True); kU = [sum(eta[mu, nu]*kk_[nu] for nu in range(4)) for mu in range(4)]; ksq = sum(kk_[mu]*kU[mu] for mu in range(4))
b0, b1, b2 = sp.symbols('b0:3', real=True)
Kmat = sp.Matrix(4, 4, lambda mu, nu: (b0 + b1*ksq)*eta[mu, nu] + b2*kU[mu]*kU[nu])
gauge = [sp.Poly(sp.expand(sum(Kmat[mu, nu]*kk_[nu] for nu in range(4))), *kk_).coeffs() for mu in range(4)]
ok("5b. the same in Fourier form: K^{mu nu}(k) = b0 eta + b1 k^2 eta + b2 k k with K^{mu nu} k_nu = 0 for every k forces b0 = 0, b2 = -b1",
   sp.solve([c_ for lst in gauge for c_ in lst], [b0, b1, b2], dict=True) == [{b0: 0, b1: -b2}])
F = [[d(a[nu], mu) - d(a[mu], nu) for nu in range(4)] for mu in range(4)]
F2 = sum(eta[mu, mu]*eta[nu, nu]*F[mu][nu]**2 for mu in range(4) for nu in range(4))
ok("5'. T1 - T2 = F_{mu nu} F^{mu nu} / 2 identically", sp.simplify(sp.expand(T1 - T2 - F2/2)) == 0)
# 6. the source term and its sign.  (Rebuilt 10 Nov 2026: the earlier check used -j.a and accepted either overall sign, so it could not tell
#    div E = rho from div E = -rho; the text says the check is strict, and now it is.)
#    Signature (-,+,+,+).  a_mu (lower index) = (-Phi, A);  j^mu = (rho, j);  j^mu a_mu = -rho Phi + j.A;  E = -grad Phi - d_t A,  B = curl A.
j = [sp.Function(f'j{mu}', real=True)(*X) for mu in range(4)]           # j^mu, upper index
FU = [[eta[mu, mu]*eta[nu, nu]*F[mu][nu] for nu in range(4)] for mu in range(4)]
Phi_ = -a[0]; Avec = [a[1], a[2], a[3]]
Evec = [-d(Phi_, i) - d(Avec[i - 1], 0) for i in (1, 2, 3)]
Bvec = [d(Avec[2], 2) - d(Avec[1], 3), d(Avec[0], 3) - d(Avec[2], 1), d(Avec[1], 1) - d(Avec[0], 2)]
divE = sum(d(Evec[i - 1], i) for i in (1, 2, 3))
curlB = [d(Bvec[2], 2) - d(Bvec[1], 3), d(Bvec[0], 3) - d(Bvec[2], 1), d(Bvec[1], 1) - d(Bvec[0], 2)]
def eom(sign):                                                          # Euler-Lagrange of -F^2/4 + sign * j^mu a_mu, as (expression = 0)
    Ls = -F2/4 + sign*sum(j[mu]*a[mu] for mu in range(4))
    return [sp.expand(eq.lhs) for eq in euler_equations(Ls, a, X)]
plus, minus = eom(+1), eom(-1)
cov = [sp.expand(sum(d(FU[nu][mu], mu) for mu in range(4)) - j[nu]) for nu in range(4)]          # d_mu F^{nu mu} - j^nu
ok("6. with the source term +j^mu a_mu the equations of motion are d_mu F^{nu mu} - j^nu = 0 in every component (the equation's overall factor is free; "
   "the sign of j against F is not, and is tested)",
   all(sp.simplify(plus[nu] - cov[nu]) == 0 or sp.simplify(plus[nu] + cov[nu]) == 0 for nu in range(4))
   and not any(sp.simplify(minus[nu] - cov[nu]) == 0 or sp.simplify(minus[nu] + cov[nu]) == 0 for nu in range(4)))
gauss_plus = sp.solve(plus[0], j[0])[0]; gauss_minus = sp.solve(minus[0], j[0])[0]
ok("6a. in E and B: +j.a gives div E = +rho and curl B - d_t E = +j (like charges repel); F^{0i} = E_i",
   sp.simplify(gauss_plus - divE) == 0 and all(sp.simplify(sp.solve(plus[i], j[i])[0] - (curlB[i - 1] - d(Evec[i - 1], 0))) == 0 for i in (1, 2, 3))
   and all(sp.simplify(FU[0][i] - Evec[i - 1]) == 0 for i in (1, 2, 3)))
ok("6b. the other sign, -j.a, gives div E = -rho: like charges would attract. The sign of the source term is not free once V.2's turn angle is given",
   sp.simplify(gauss_minus + divE) == 0)
ok("6c. j^mu a_mu = -rho Phi + j.A in this signature, the turn angle of V.2 per unit volume and time",
   sp.simplify(sum(j[mu]*a[mu] for mu in range(4)) - (-j[0]*Phi_ + sum(j[i]*Avec[i - 1] for i in (1, 2, 3)))) == 0)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
