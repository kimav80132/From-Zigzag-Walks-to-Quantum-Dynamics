"""Klein-Gordon with the electromagnetic field, factorised; what the two-step moment ('noise') of the Dirac process must be.
A. Part II level (1+1, eta = diag(-1,1)): Madelung of (D^mu D_mu - m^2c^2/hbar^2) psi = 0 with D_mu = d_mu - i e A_mu/hbar:
   the drift becomes (d S - e A)/m, the osmotic term (the noise) is unchanged; a c-number antisymmetric addition to the
   noise, (eta + kappa F)^{mu nu} d_mu d_nu, is identically Box: it does nothing.
B. Factorisation with IV.4b's real matrices: beta = i T, gamma^0 = beta, gamma^k = beta A_k obey the Clifford relations;
   gamma^mu gamma^nu = eta^{mu nu} + (1/2)[gamma^mu, gamma^nu], with (1/2)[gamma^0,gamma^k] = A_k (advances)
   and (1/2)[gamma^j,gamma^k] = -A_j A_k (minus twice the (j,k) turn).
C. Square of Dirac with the field (Hamiltonian form, general Phi(t,x,y,z), A(t,x,y,z)):
   (i hbar d_t - e Phi)^2 psi = [c^2 (p - eA)^2 + m^2c^4 - e hbar c^2 Sigma.B + i e hbar c A.E] psi,
   Sigma_z = -i A_x A_y etc.  Klein-Gordon with the field + the two antisymmetric pieces: turn x B and advance x E.
D. Slow limit: H - e Phi = mc^2 + pi^2/2m - (e hbar/2m) Sigma.B + ...:  g = 2.
"""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
# ---------- A
x0, x1 = sp.symbols('x0 x1', real=True); X = (x0, x1); eta = sp.diag(-1, 1)
e, hb, m, c, kap = sp.symbols('e hbar m c kappa', positive=True)
rho = sp.Function('rho', positive=True)(*X); S = sp.Function('S', real=True)(*X)
A0, A1 = [sp.Function(n, real=True)(*X) for n in ('A0', 'A1')]; A = (A0, A1)
psi = sp.sqrt(rho)*sp.exp(sp.I*S/hb)
D = lambda f, mu: sp.diff(f, X[mu]) - sp.I*e/hb*A[mu]*f
KG = sum(eta[mu, mu]*D(D(psi, mu), mu) for mu in range(2)) - m**2*c**2/hb**2*psi
expr = sp.expand(sp.simplify(KG/psi))                 # all symbols real: split on the explicit i
re, im = sp.simplify(expr.subs(sp.I, 0)), sp.simplify(expr.coeff(sp.I))
pi_ = [sp.diff(S, X[mu]) - e*A[mu] for mu in range(2)]
pi_up = [eta[mu, mu]*pi_[mu] for mu in range(2)]
box_sqrt = sum(eta[mu, mu]*sp.diff(sp.sqrt(rho), X[mu], 2) for mu in range(2))/sp.sqrt(rho)
HJ = box_sqrt - sum(pi_up[mu]*pi_[mu] for mu in range(2))/hb**2 - m**2*c**2/hb**2
cont = sum(sp.diff(rho*pi_up[mu], X[mu]) for mu in range(2))/(hb*rho)
ok("A1. real part: Box sqrt(rho)/sqrt(rho) - (dS - eA).(dS - eA)/hbar^2 - m^2c^2/hbar^2 = 0  (drift shifted, noise term unchanged)",
   sp.simplify(re - HJ) == 0)
ok("A2. imaginary part: d_mu [rho (d^mu S - e A^mu)] = 0  (continuity with the shifted drift)", sp.simplify(im - cont) == 0)
F01 = sp.symbols('F01', real=True); F = sp.Matrix([[0, F01], [-F01, 0]]); f = sp.Function('f')(*X)
noise_op = sum((eta + kap*F)[mu, nu]*sp.diff(f, X[mu], X[nu]) for mu in range(2) for nu in range(2))
ok("A3. c-number antisymmetric noise: (eta + kappa F)^{mu nu} d_mu d_nu f = Box f identically (the field drops out)",
   sp.simplify(noise_op - sum(eta[mu, mu]*sp.diff(f, X[mu], 2) for mu in range(2))) == 0)
# ---------- B
Az = sp.diag(1, -1, -1, 1)
Ax = sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]])
Ay = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
T = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]])
beta = sp.I*T; Ak = [Ax, Ay, Az]
g = [beta] + [beta*Ak_ for Ak_ in Ak]; eta4 = sp.diag(1, -1, -1, -1); I4 = sp.eye(4)   # (+,-,-,-) for these gammas
ok("B1. beta = iT: beta^2 = 1, beta anticommutes with A_x, A_y, A_z; gamma^mu gamma^nu + gamma^nu gamma^mu = 2 eta^{mu nu}",
   beta*beta == I4 and all(beta*Ak_ + Ak_*beta == sp.zeros(4) for Ak_ in Ak)
   and all(g[mu]*g[nu] + g[nu]*g[mu] == 2*eta4[mu, nu]*I4 for mu in range(4) for nu in range(4)))
ok("B2. (1/2)[gamma^0, gamma^k] = A_k (the advances);  (1/2)[gamma^j, gamma^k] = -A_j A_k (minus twice the (j,k) turn)",
   all((g[0]*g[k] - g[k]*g[0])/2 == Ak[k - 1] for k in (1, 2, 3))
   and all((g[j]*g[k] - g[k]*g[j])/2 == -Ak[j - 1]*Ak[k - 1] for j in (1, 2, 3) for k in (1, 2, 3) if j != k))
# ---------- C
t, x, y, z = sp.symbols('t x y z', real=True); R = (x, y, z)
Phi = sp.Function('Phi', real=True)(t, x, y, z); Av = [sp.Function(n, real=True)(t, x, y, z) for n in ('Ax', 'Ay', 'Az')]
psi4 = sp.Matrix([sp.Function(f'p{i}')(t, x, y, z) for i in range(4)])
pi = lambda G, j: -sp.I*hb*G.diff(R[j]) - e*Av[j]*G
H = lambda G: c*sum((Ak[j]*pi(G, j) for j in range(3)), sp.zeros(4, 1)) + m*c**2*beta*G + e*Phi*G
Dt = lambda G: sp.I*hb*G.diff(t) - e*Phi*G
E = [-sp.diff(Phi, R[j]) - sp.diff(Av[j], t) for j in range(3)]
B = [sp.diff(Av[2], y) - sp.diff(Av[1], z), sp.diff(Av[0], z) - sp.diff(Av[2], x), sp.diff(Av[1], x) - sp.diff(Av[0], y)]
Sig = [-sp.I*Ay*Az, -sp.I*Az*Ax, -sp.I*Ax*Ay]
ok("C0. Sigma_k = -i A_j A_l (cyclic) are Hermitian, square to 1, and commute with beta", all(Sg == Sg.H and Sg*Sg == I4 and Sg*beta == beta*Sg for Sg in Sig))
# (i hbar d_t - e Phi) psi = (H - e Phi) psi  ->  apply (i hbar d_t - e Phi) to both sides
lhs = Dt(H(psi4) - e*Phi*psi4)                       # (i hbar d_t - e Phi)^2 psi, using the equation once
pi2 = sum((pi(pi(psi4, j), j) for j in range(3)), sp.zeros(4, 1))
rhs = c**2*pi2 + m**2*c**4*psi4 - e*hb*c**2*sum((Sig[j]*B[j] for j in range(3)), sp.zeros(4))*psi4 \
      + sp.I*e*hb*c*sum((Ak[j]*E[j] for j in range(3)), sp.zeros(4))*psi4
# lhs contains d_t psi; replace it by the equation
sub = {psi4[i].diff(t): ((H(psi4))[i])/(sp.I*hb) for i in range(4)}
lhs_s = lhs.subs(sub).doit()
ok("C1. (i hbar d_t - e Phi)^2 psi = [c^2 (p - eA)^2 + m^2 c^4 - e hbar c^2 Sigma.B + i e hbar c A.E] psi",
   sp.simplify(sp.expand(lhs_s - rhs)) == sp.zeros(4, 1))
# ---------- D
b, p2 = sp.symbols('b p2', real=True)
Hsq = m**2*c**4 + c**2*p2 - e*hb*c**2*b                 # b = Sigma.B eigenvalue, p2 = pi^2
lam = sp.symbols('lambda', positive=True)             # lam = 1/c^2 as the small parameter
ser = sp.series(sp.sqrt(Hsq.subs(c, 1/sp.sqrt(lam))), lam, 0, 1).removeO()
ok(f"D. sqrt(H^2) = m c^2 + p^2/2m - (e hbar/2m) Sigma.B + O(1/c^2): magnetic moment (e hbar/2m) Sigma = (e/m) S, g = 2",
   sp.simplify(ser - (m/lam + p2/(2*m) - e*hb*b/(2*m))) == 0)
