"""IV.8: the fluid of the axis, derived in the master instead of cited from the spin/chain documents.
Psi = sqrt(rho) e^{iS/hbar} chi, chi = (cos th/2, e^{i phi} sin th/2), s = chi^dag sigma chi; rho, S, th, phi generic test fields of (x,y,z).
Every identity is checked by exact differentiation of the test fields and evaluation at random points (max |difference| < 1e-9).
A. completed square: (hbar^2/2m)|grad Psi|^2 = (hbar^2/2m)(grad sqrt rho)^2 + (m/2) rho v^2 + (hbar^2/8m) rho |grad s|^2, m v = grad S + hbar sin^2(th/2) grad phi;
   (hbar/m) Im[Psi^dag grad Psi] = rho v;  chi^dag grad chi = i sin^2(th/2) grad phi;  |grad chi|^2 = |grad s|^2/4 + sin^4(th/2)(grad phi)^2.
B. vorticity curl v = (hbar/2m) sin th grad th x grad phi; circulation (h/m)[n + k sin^2(th/2)].
C. slow limit of Dirac: eta = sigma.p phi/2mc; current 2c Re(phi^dag sigma eta) = rho v + (hbar/2m) curl(rho s); rho_R - rho_L; (S_L - S_R)/hbar; spin current = u x s + D curl s.
D. (sigma.pi)^2 = pi^2 - (e hbar/c) sigma.B.
E. angular momenta of a Gaussian cloud with uniform axis: Gordon curl carries hbar; Belinfante half carries hbar/2; orbit hbar/2; moment e hbar/2mc (g = 2)."""
import sympy as sp, numpy as np
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
i = sp.I
x, y, z = sp.symbols('x y z', real=True); X = (x, y, z)
hbar, m, c, e = sp.Rational(7, 5), sp.Rational(3, 2), sp.Rational(11, 4), sp.Rational(2, 3)
rng = np.random.default_rng(3); PTS = rng.uniform(-0.9, 0.9, size=(4, 3))
def zero(ex):
    f = sp.lambdify(X, ex, 'numpy')
    return max(abs(complex(f(*p))) for p in PTS) < 1e-9
sx = sp.Matrix([[0, 1], [1, 0]]); sy = sp.Matrix([[0, -i], [i, 0]]); sz = sp.diag(1, -1); sig = (sx, sy, sz)
rho = sp.exp(x*y/3 - z**2/5 + x/2); S = 2*x*z/3 + y**2/4 - x*y*z/7
th = 1 + x/3 - y*z/4 + z**2/6; ph = y/2 + x*z/3 - x**2/5
chi = sp.Matrix([sp.cos(th/2), sp.exp(i*ph)*sp.sin(th/2)])
Psi = sp.sqrt(rho)*sp.exp(i*S/hbar)*chi
grad = lambda f: sp.Matrix([sp.diff(f, v) for v in X])
div = lambda F: sum(sp.diff(F[k], X[k]) for k in range(3))
curl = lambda F: sp.Matrix([sp.diff(F[2], y) - sp.diff(F[1], z), sp.diff(F[0], z) - sp.diff(F[2], x), sp.diff(F[1], x) - sp.diff(F[0], y)])
conjT = lambda M: M.T.applyfunc(sp.conjugate)          # Hermitian conjugate without assumptions
Re = lambda ex: (ex + sp.conjugate(ex))/2; Im = lambda ex: (ex - sp.conjugate(ex))/(2*i)
s = sp.Matrix([sp.sin(th)*sp.cos(ph), sp.sin(th)*sp.sin(ph), sp.cos(th)])
ok("0. axis: chi^dag sigma chi = (sin th cos phi, sin th sin phi, cos th)", all(zero((conjT(chi)*M*chi)[0] - s[k]) for k, M in enumerate(sig)))
mv = grad(S) + hbar*sp.sin(th/2)**2*grad(ph)
# A
gPsi = [Psi.diff(v) for v in X]
lhs = hbar**2/(2*m)*sum((conjT(gPsi[k])*gPsi[k])[0] for k in range(3))
grad_s_sq = sum(sp.diff(s[j], X[k])**2 for j in range(3) for k in range(3))
rhs = hbar**2/(2*m)*grad(sp.sqrt(rho)).dot(grad(sp.sqrt(rho))) + rho/(2*m)*mv.dot(mv) + hbar**2/(8*m)*rho*grad_s_sq
ok("A1. (hbar^2/2m)|grad Psi|^2 = (hbar^2/2m)(grad sqrt rho)^2 + (m/2) rho v^2 + (hbar^2/8m) rho |grad s|^2, with m v = grad S + hbar sin^2(th/2) grad phi", zero(lhs - rhs))
ok("A2. (hbar/m) Im[Psi^dag grad Psi] = rho v: the convective current is the velocity of the square",
   all(zero(hbar/m*Im((conjT(Psi)*gPsi[k])[0]) - rho*mv[k]/m) for k in range(3)))
ok("A3. chi^dag grad chi = i sin^2(th/2) grad phi, and |grad chi|^2 = |grad s|^2/4 + sin^4(th/2)(grad phi)^2",
   all(zero((conjT(chi)*chi.diff(v))[0] - i*sp.sin(th/2)**2*sp.diff(ph, v)) for v in X)
   and zero(sum((conjT(chi.diff(v))*chi.diff(v))[0] for v in X) - (grad_s_sq/4 + sp.sin(th/2)**4*grad(ph).dot(grad(ph)))))
# B
ok("B1. curl v = (hbar/2m) sin th (grad th x grad phi): irrotational except where the axis twists",
   all(zero(curl(mv/m)[k] - hbar/(2*m)*(sp.sin(th)*grad(th).cross(grad(ph)))[k]) for k in range(3)))
k_, n_, th0 = sp.symbols('k n theta0')
circ = (2*sp.pi*hbar*n_ + hbar*sp.sin(th0/2)**2*2*sp.pi*k_)/m
ok("B2. circulation on a loop with th = th0, phi winding k, S winding n: (h/m)[n + k sin^2(th0/2)]; on the equator the extra is half a quantum, pi hbar/m",
   sp.simplify(circ - 2*sp.pi*hbar/m*(n_ + k_*sp.sin(th0/2)**2)) == 0 and sp.simplify(circ.subs({n_: 0, k_: 1, th0: sp.pi/2}) - sp.pi*hbar/m) == 0)
# C
p_ = lambda col: [sp.Matrix([-i*hbar*sp.diff(col[r], X[k]) for r in range(2)]) for k in range(3)]
sigp = lambda col: sum((sig[k]*p_(col)[k] for k in range(3)), sp.zeros(2, 1))
eta = sigp(Psi)/(2*m*c)
jD = [2*c*Re((conjT(Psi)*sig[k]*eta)[0]) for k in range(3)]
ok("C1. slow-limit Dirac current 2c Re(phi^dag sigma eta), eta = sigma.p phi/2mc, = rho v + (hbar/2m) curl(rho s): the full vector P = rho m v + (hbar/2) curl(rho s) of IV.7",
   all(zero(jD[k] - (rho*mv[k]/m + hbar/(2*m)*curl(rho*s)[k])) for k in range(3)))
ok("C2. rho_R - rho_L = 2 Re(phi^dag eta) = (rho/c)[s.v + (hbar/2m) s.curl s]: III.11's c cos theta = v, projected on the axis, plus the twist",
   zero(2*Re((conjT(Psi)*eta)[0]) - rho/c*(s.dot(mv)/m + hbar/(2*m)*s.dot(curl(s)))))
u = hbar/(2*m)*grad(sp.log(rho))
ok("C3. (S_L - S_R)/hbar = -2 Im(phi^dag eta)/rho = (1/c)[s.u + (hbar/2m) div s]: III.11's c phi = u, projected, plus the divergence of the axis",
   zero(-2*Im((conjT(Psi)*eta)[0])/rho - (s.dot(u) + hbar/(2*m)*div(s))/c))
ok("C4. (hbar/2m rho) curl(rho s) = u x s + (hbar/2m) curl s: for a uniform axis the spin current is the osmotic velocity crossed with the axis",
   all(zero(hbar/(2*m*rho)*curl(rho*s)[k] - (u.cross(s) + hbar/(2*m)*curl(s))[k]) for k in range(3)))
# D
Avec = sp.Matrix([y*z - x**2/3, x*y + z/2, x*z/3 - y**2/4])
col = sp.Matrix([sp.exp(x/3 + y*z/5)*(1 + i*x*y/4), sp.cos(x - y/2) + i*z*x/3])
pi_ = lambda col: [sp.Matrix([-i*hbar*sp.diff(col[r], X[k]) - e/c*Avec[k]*col[r] for r in range(2)]) for k in range(3)]
sigpi = lambda col: sum((sig[k]*pi_(col)[k] for k in range(3)), sp.zeros(2, 1))
pi2 = sum((pi_(pi_(col)[k])[k] for k in range(3)), sp.zeros(2, 1))
B = curl(Avec)
rhsD = pi2 - e*hbar/c*sum((sig[k]*B[k] for k in range(3)), sp.zeros(2, 2))*col
ok("D. (sigma.pi)^2 = pi^2 - (e hbar/c) sigma.B on a two-component column: the Pauli term, g = 2, from the square of the direction matrices",
   all(zero(ex) for ex in (sigpi(sigpi(col)) - rhsD)))
# E
a_ = sp.symbols('a', positive=True); hb, mm, cc, ee = sp.symbols('hbar m c e', positive=True)
rhoG = (sp.pi*a_**2)**sp.Rational(-3, 2)*sp.exp(-(x**2 + y**2 + z**2)/a_**2)
zhat = sp.Matrix([0, 0, 1]); r = sp.Matrix([x, y, z])
Lz = lambda G: sp.integrate(r.cross(G)[2], (x, -sp.oo, sp.oo), (y, -sp.oo, sp.oo), (z, -sp.oo, sp.oo))
ok("E1. int r x (hbar/2) curl(rho s) = hbar s: the Gordon curl term of P carries angular momentum hbar, twice the spin density's hbar/2",
   sp.simplify(Lz(hb/2*curl(rhoG*zhat)) - hb) == 0)
ok("E2. int r x (1/2) curl(rho (hbar/2) s) = (hbar/2) s: the energy-flow (Belinfante) half carries exactly the spin", sp.simplify(Lz(hb/4*curl(rhoG*zhat)) - hb/2) == 0)
ok("E3. the orbit of true_spin.py 3b: radius hbar/2mc at speed c, angular momentum m r c = hbar/2, loop 2 pi r c = (h/m)/2, the half quantum of B2",
   sp.simplify(mm*(hb/(2*mm*cc))*cc - hb/2) == 0 and sp.simplify(2*sp.pi*(hb/(2*mm*cc))*cc - sp.pi*hb/mm) == 0)
ok("E4. magnetic moment of the Gordon current (e/2c) int r x j_spin = (e hbar/2mc) s = (e/mc)(hbar/2) s: g = 2 is the factor between E1 and E2",
   sp.simplify(ee/(2*cc)*Lz(hb/(2*mm)*curl(rhoG*zhat)) - ee*hb/(2*mm*cc)) == 0)
# E2'. the Belinfante half, derived: for the Dirac field T^{0k}_sym gives G = (hbar/2) Im(psi^dag grad psi) + (m c/2) psi^dag alpha psi
# (the second term from the symmetrized time derivative with psi ~ e^{-i m c^2 t/hbar}); with psi = (phi, eta), eta = sigma.p phi/2mc:
G_can_half = [hbar/2*Im((conjT(Psi)*gPsi[k])[0] + (conjT(eta)*eta.diff(X[k]))[0]) for k in range(3)]
G_sym = [G_can_half[k] + m/2*jD[k] for k in range(3)]
ok("E2'. G_sym = (hbar/2) Im(psi^dag grad psi) + (m/2) j = rho m v + (hbar/4) curl(rho s) + (hbar/2) Im(eta^dag grad eta): "
   "half the canonical momentum plus half of m times the Gordon current; the curl part is half the Gordon curl, up to the O(1/c^2) small-component term",
   all(zero(G_sym[k] - (rho*mv[k] + hbar/4*curl(rho*s)[k] + hbar/2*Im((conjT(eta)*eta.diff(X[k]))[0]))) for k in range(3)))
# C5. uniform axis: v = grad S/m stays irrotational; the whole spin current is rho (u x s), built from the osmotic velocity, and it transports no density
thc, phc = sp.Rational(7, 10), sp.Rational(3, 10)
chi_u = sp.Matrix([sp.cos(thc/2), sp.exp(i*phc)*sp.sin(thc/2)]); Psi_u = sp.sqrt(rho)*sp.exp(i*S/hbar)*chi_u
s_u = sp.Matrix([sp.sin(thc)*sp.cos(phc), sp.sin(thc)*sp.sin(phc), sp.cos(thc)])
eta_u = sigp(Psi_u)/(2*m*c); jD_u = [2*c*Re((conjT(Psi_u)*sig[k]*eta_u)[0]) for k in range(3)]
ok("C5. uniform axis: j = rho grad S/m + rho (u x s) exactly; grad S/m is a gradient (irrotational) and div[rho (u x s)] = 0: the circulation is osmotic and carries no density",
   all(zero(jD_u[k] - (rho*grad(S)[k]/m + rho*u.cross(s_u)[k])) for k in range(3)) and zero(div(rho*u.cross(s_u))))
# E2''. (added 10 Nov 2026) the formula E2' starts from, obtained from the symmetric stress tensor itself.  Dirac representation, psi = e^{-i m c^2 t/hbar} (Psi, eta):
#   momentum density G_k = T^{0k}_sym / c = (hbar/2) Im(psi^dag d_k psi) - (hbar/2c) Im(psi^dag alpha_k d_t psi)
tt = sp.symbols('t', real=True)
psi4 = sp.exp(-i*m*c**2*tt/hbar)*sp.Matrix([Psi[0], Psi[1], eta[0], eta[1]])
alpha = [sp.Matrix(sp.BlockMatrix([[sp.zeros(2), M], [M, sp.zeros(2)]])) for M in sig]
psi4d = sp.exp(i*m*c**2*tt/hbar)*sp.Matrix([[sp.conjugate(Psi[0]), sp.conjugate(Psi[1]), sp.conjugate(eta[0]), sp.conjugate(eta[1])]])
G_tensor = [hbar/2*Im((psi4d*psi4.diff(X[k]))[0]) - hbar/(2*c)*Im((psi4d*alpha[k]*psi4.diff(tt))[0]) for k in range(3)]
ok("E2''. the symmetric stress tensor's momentum density, (hbar/2) Im(psi^dag d_k psi) - (hbar/2c) Im(psi^dag alpha_k d_t psi), on psi = e^{-i m c^2 t/hbar}(phi, eta) "
   "is the G_sym of E2': half the canonical momentum plus half of m times the current", all(zero(sp.simplify(G_tensor[k] - G_sym[k]).subs(tt, 0)) for k in range(3)))
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
