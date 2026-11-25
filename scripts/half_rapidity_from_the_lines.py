"""Where the half comes from: boost the lines of III.5a (H) and ask what the readings must do for the lines to keep their form. c = 1.
Lines: d_t Q = -A d_x Q + a T Q, A = diag(1,-1), T = [[0,1],[-1,0]]. Boost t' = g(t - b x), x' = g(x - b t), g = cosh eta, b = tanh eta.
1. d_t = g(d_t' - b d_x'), d_x = g(-b d_t' + d_x'). Substituting: g(1 - bA) d_t' Q + g(A - b) d_x' Q = a T Q, and A - b = A(1 - bA) since A^2 = 1, so
       g(1 - bA)[d_t' + A d_x'] Q = a T Q,   and g(1 - bA) = e^{-eta A} (= diag(e^{-eta}, e^{eta})): the FULL rapidity on each null line.
2. Multiply by e^{eta A}: [d_t' + A d_x'] Q = a e^{eta A} T Q. The mass line has picked up e^{eta A}. Because A T = -T A, e^{eta A} T = T e^{-eta A}, and the only
   constant rescaling Q = Lambda Q' that commutes with A and restores the form is Lambda = e^{eta A/2}:  e^{-eta A/2} e^{eta A} T e^{eta A/2} = e^{eta A/2} T e^{eta A/2}
       = e^{eta A/2} e^{-eta A/2} T = T.  So Q' = e^{-eta A/2} Q: Q'+ = e^{-eta/2} Q+, Q'- = e^{+eta/2} Q-.  The half is the square root of the null-line factor,
   forced by the mass term, which couples the two lines and must stay form-invariant.
3. In components, line by line: e^{-eta}(d_t' + d_x') Q+ = a Q-,  e^{eta}(d_t' - d_x') Q- = -a Q+. With Q+ = lam Q'+, Q- = mu Q'-: form invariance needs
       e^{-eta} lam = mu  and  e^{eta} mu = lam, i.e. lam/mu = e^{eta} (consistent), and with lam mu = 1 (the content's normalization), lam = e^{eta/2}, mu = e^{-eta/2}.
4. The invariants: Q+ Q- (scalar); the content's null parts 2Q+^2, 2Q-^2 go as (e^{+-eta}); j.j = (j0)^2 - (j1)^2 = 4 Q+^2 Q-^2, so the rest-frame density
   sqrt(j.j) = 2 Q+ Q- is the geometric mean of the two null fluxes, as the proper time sqrt(t1 t2) is the geometric mean of the two null coordinates."""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
eta, a = sp.symbols('eta a', real=True); t, x, tp, xp = sp.symbols('t x tp xp', real=True)
g, b = sp.cosh(eta), sp.tanh(eta)
A = sp.diag(1, -1); T = sp.Matrix([[0, 1], [-1, 0]]); I2 = sp.eye(2)
expA = lambda s: sp.cosh(s)*I2 + sp.sinh(s)*A                     # e^{s A}
ok("1a. g(1 - bA) = e^{-eta A} = diag(e^{-eta}, e^{eta}): the boost puts the full rapidity on each null line",
   (g*(I2 - b*A) - expA(-eta)).applyfunc(lambda f: sp.simplify(f.rewrite(sp.exp))) == sp.zeros(2, 2)
   and (expA(-eta) - sp.diag(sp.exp(-eta), sp.exp(eta))).applyfunc(lambda f: sp.simplify(f.rewrite(sp.exp))) == sp.zeros(2, 2))
ok("1b. A - b = A(1 - bA)", sp.simplify(A - b*I2 - A*(I2 - b*A)) == sp.zeros(2, 2))
# full substitution on generic solutions: Q(t,x) = Q'(t',x') with the primed functions
Qp = sp.Matrix([sp.Function('Qp1')(tp, xp), sp.Function('Qp2')(tp, xp)])
Lam = expA(eta/2)
Q = Lam*Qp                                                            # claim: Q = e^{eta A/2} Q'
# derivatives of Q wrt t, x through t', x'
dt_Q = g*sp.diff(Q, tp) - g*b*sp.diff(Q, xp); dx_Q = -g*b*sp.diff(Q, tp) + g*sp.diff(Q, xp)
residual_unprimed = dt_Q + A*dx_Q - a*T*Q                             # the lines in the unprimed frame, applied to Q = Lam Q'
residual_primed = sp.diff(Qp, tp) + A*sp.diff(Qp, xp) - a*T*Qp       # the same lines in the primed frame, applied to Q'
# claim: residual_unprimed = e^{-eta A/2} residual_primed  (so Q solves the lines iff Q' solves them)
diff = (residual_unprimed - expA(-eta/2)*residual_primed).applyfunc(lambda f: sp.simplify(sp.expand(f.rewrite(sp.exp))))
ok("2. with Q = e^{eta A/2} Q', the unprimed lines on Q equal e^{-eta A/2} times the primed lines on Q': the lines keep their form iff Q'+- = e^{-+eta/2} Q+-",
   diff == sp.zeros(2, 1))
ok("2'. e^{-eta A/2} e^{eta A} T e^{eta A/2} = T (A anticommutes with T)",
   (expA(-eta/2)*expA(eta)*T*expA(eta/2) - T).applyfunc(lambda f: sp.simplify(f.rewrite(sp.exp))) == sp.zeros(2, 2))
lam, mu = sp.symbols('lambda mu', positive=True)
sol = sp.solve([sp.exp(-eta)*lam - mu, sp.exp(eta)*mu - lam, lam*mu - 1], [lam, mu], dict=True)
ok(f"3. the component conditions e^-eta lam = mu, e^eta mu = lam, lam mu = 1 give lam = e^{{eta/2}}, mu = e^{{-eta/2}}: {sol}",
   any(sp.simplify(s_[lam] - sp.exp(eta/2)) == 0 and sp.simplify(s_[mu] - sp.exp(-eta/2)) == 0 for s_ in sol))
Qp_, Qm_ = sp.symbols('Q_p Q_m', positive=True)
j0, j1 = Qp_**2 + Qm_**2, Qp_**2 - Qm_**2
ok("4. j.j = (j0)^2 - (j1)^2 = (j0+j1)(j0-j1) = 4 Q+^2 Q-^2: the rest-frame density 2 Q+ Q- is the geometric mean of the two null fluxes 2Q+^2 and 2Q-^2",
   sp.simplify(j0**2 - j1**2 - 4*Qp_**2*Qm_**2) == 0 and sp.simplify(sp.sqrt((2*Qp_**2)*(2*Qm_**2)) - 2*Qp_*Qm_) == 0)
