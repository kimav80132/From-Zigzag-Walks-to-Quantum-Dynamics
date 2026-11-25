"""Does frame independence (II.1) force the record to be out x return? 1+1, the t-reading of III.5a (H), with and without mass.
A record on a cut is what all frames agree on: a flux of a conserved current through the cut, int j^mu dSigma_mu. So the record must be the time
component of a conserved FOUR-VECTOR. Test what the readings can build.
1. under a run by rapidity eta (IV.4d (D)-(E)) the readings go Q+ -> e^{eta/2} Q+, Q- -> e^{-eta/2} Q-: HALF the boost. The content (Q+^2+Q-^2, Q+^2-Q-^2)
   goes as a four-vector, Lambda(eta) j. The pair (Q+, Q-) does not go as any four-vector: its light-cone parts scale by e^{+-eta/2}, a vector's by e^{+-eta}.
2. with mass there is no conserved current linear in (Q+, Q-) at all: d_t(alpha Q+ + beta Q-) + c d_x(gamma Q+ + delta Q-) = 0 for all solutions has no solution but zero.
   The bilinear one, (Q+^2+Q-^2, c(Q+^2-Q-^2)), is conserved (IV.4d (B)).
3. without mass, Q+ alone is conserved along its light line, (Q+, c Q+) has zero divergence - but it is a null pair that scales by e^{eta/2}, not a four-vector:
   boosting the pair and boosting the solution give different numbers. Only its square, (Q+^2, c Q+^2), is a vector.
So a frame-independent count is bilinear in the readings: one reading on each side (III.8a's vector), and III.15a names the two: the out leg and the return leg."""
import sympy as sp, numpy as np
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
eta = sp.symbols('eta', real=True); Qp, Qm = sp.symbols('Q_p Q_m', real=True)
Lam = sp.Matrix([[sp.cosh(eta), sp.sinh(eta)], [sp.sinh(eta), sp.cosh(eta)]])          # boost on (j^0, j^1), c = 1
j = sp.Matrix([Qp**2 + Qm**2, Qp**2 - Qm**2])
j_run = j.subs({Qp: sp.exp(eta/2)*Qp, Qm: sp.exp(-eta/2)*Qm})
ok("1a. the content of the run readings equals the boosted content: (Q+^2+Q-^2, Q+^2-Q-^2) is a four-vector", (j_run - Lam*j).applyfunc(lambda f: sp.simplify(sp.expand(f.rewrite(sp.exp)))) == sp.zeros(2, 1))
q = sp.Matrix([Qp, Qm]); q_run = q.subs({Qp: sp.exp(eta/2)*Qp, Qm: sp.exp(-eta/2)*Qm})
# is there ANY constant 2x2 M with M q transforming as a vector, Lam (M q) = M q_run, for all Q and eta? Solve.
M = sp.Matrix(2, 2, sp.symbols('m0:4', real=True))
eqs = sp.expand(Lam*M*q - M*q_run)
conds = []
for ex in eqs:
    ser = sp.series(ex, eta, 0, 3).removeO()
    conds += [sp.Poly(ser.coeff(eta, k), Qp, Qm).coeffs() for k in range(3)]
conds = [c_ for lst in conds for c_ in lst]
sol = sp.solve(conds, list(M), dict=True)
ok(f"1b. no linear combination of the readings transforms as a four-vector: the only solution of Lam (M q) = M q_run is M = 0 ({sol})",
   sol == [{m: 0 for m in M}] or all(v == 0 for v in sol[0].values()))
# 2. conserved linear current with mass? Use the equations d_t Q+ = -c d_x Q+ + a Q-, d_t Q- = c d_x Q- - a Q+ on a generic solution (symbolic functions)
t, x = sp.symbols('t x', real=True); c, a = sp.symbols('c a', positive=True)
Qpf, Qmf = sp.Function('Qp')(t, x), sp.Function('Qm')(t, x)
al, be, ga, de = sp.symbols('alpha beta gamma delta', real=True)
dtQp = -c*sp.diff(Qpf, x) + a*Qmf; dtQm = c*sp.diff(Qmf, x) - a*Qpf
div = al*dtQp + be*dtQm + c*sp.diff(ga*Qpf + de*Qmf, x)
div = sp.expand(div)
coeffs = [div.coeff(sp.diff(Qpf, x)), div.coeff(sp.diff(Qmf, x)), div.coeff(Qpf), div.coeff(Qmf)]
sol2 = sp.solve(coeffs, [al, be, ga, de], dict=True)
ok(f"2. with mass, no current linear in the readings is conserved (only alpha = beta = gamma = delta = 0: {sol2}); the bilinear content is",
   all(v == 0 for v in sol2[0].values()) and sp.simplify(sp.expand(2*Qpf*dtQp + 2*Qmf*dtQm + c*sp.diff(Qpf**2 - Qmf**2, x))) == 0)
# 3. massless: (Q+, cQ+) has zero divergence but is not a vector; (Q+^2, c Q+^2) is
null_lin = sp.Matrix([Qp, Qp]); null_sq = sp.Matrix([Qp**2, Qp**2])
ok("3. massless: (Q+, Q+) run gives e^{eta/2}(Q+, Q+), the boost gives e^{eta}(Q+, Q+): the linear null pair is not a vector; its square is",
   sp.simplify(Lam*null_lin - sp.exp(eta)*null_lin) == sp.zeros(2, 1) and sp.simplify(null_lin.subs(Qp, sp.exp(eta/2)*Qp) - sp.exp(eta/2)*null_lin) == sp.zeros(2, 1)
   and sp.simplify(Lam*null_sq - null_sq.subs(Qp, sp.exp(eta/2)*Qp)) == sp.zeros(2, 1))
