"""Does the field reach the coin (flip rate) as an output?
Square of Dirac with the field (em_noise_factorization C1):  (i hbar d_t - e Phi)^2 = c^2 pi^2 + m^2 c^4 - e hbar c^2 Sigma.B + i e hbar c alpha.E.
1. the m^2 c^4 term is untouched: the coin's rate a = m c^2/hbar is not changed by the field (no term  m c^2 x field).
2. the field couples to the TURNS, not to the exchange: [gamma^mu, gamma^nu] never contains beta (the mass matrix).
3. the Zitterbewegung frequency (gap between the two branches at p = 0) in a field: E_+ - E_- = 2 mc^2 for Phi (both shift together),
   and 2 sqrt(m^2c^4 + ...)... B splits by spin: the coin stays at a; what changes is the energy, i.e. the turn rate seen in t.
4. so a 'skew of the coin' is not produced by the field. The field's output lives in the antisymmetric noise (turns) and in the drift.
   Compare: the scalar coupling a -> a + V/hbar does change the coin: m c^2 -> m c^2 + V in the square.
"""
import sympy as sp
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
p, m, c, V, hb, e, B = sp.symbols('p m c V hbar e B', positive=True)
sz, sx = sp.diag(1, -1), sp.Matrix([[0, 1], [1, 0]])
Hem = c*p*sz + m*c**2*sx + V*sp.eye(2)                  # 1+1, V = e Phi
Hsc = c*p*sz + (m*c**2 + V)*sx
sq_em = sp.simplify((Hem - V*sp.eye(2))**2); sq_sc = sp.simplify(Hsc**2)
ok(f"1. EM: (H - e Phi)^2 = {sq_em[0, 0]} x 1: the mass term m^2 c^4 is untouched (coin unchanged)", sp.expand(sq_em[0, 0] - (c**2*p**2 + m**2*c**4)) == 0)
ok(f"2. scalar: H^2 = {sq_sc[0, 0]} x 1: the coin IS changed, m c^2 -> m c^2 + V", sp.expand(sq_sc[0, 0] - (c**2*p**2 + (m*c**2 + V)**2)) == 0)
Az = sp.diag(1, -1, -1, 1); Ax = sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]]); Ay = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
T4 = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]]); beta = sp.I*T4
g = [beta, beta*Ax, beta*Ay, beta*Az]
comms = [(g[i]*g[j] - g[j]*g[i])/2 for i in range(4) for j in range(4) if i < j]
# does any commutator have a component along beta (the mass/exchange matrix)?  beta is Hermitian with tr(beta X) = 0 for the turns
ok("3. none of the six turns [gamma^mu, gamma^nu]/2 has a component along beta (the exchange): tr(beta * turn) = 0 for all six",
   all(sp.simplify((beta*Cm).trace()) == 0 for Cm in comms))
# 4. at rest in B along z (rebuilt 10 Nov 2026: the earlier check squared a square root).  The matrix part of the square, m^2 c^4 - e hbar c^2 Sigma_z B:
Sz = -sp.I*Ax*Ay
Msq = m**2*c**4*sp.eye(4) - e*hb*c**2*B*Sz
ok("4. at rest in B: the square's matrix part m^2 c^4 - e hbar c^2 Sigma_z B has eigenvalues m^2 c^4 -+ e hbar c^2 B, each twice: the energy (the turn rate in t) "
   "splits by spin; the coefficient of 1, m^2 c^4, is the coin's and is not touched; Sigma_z commutes with beta",
   Msq.eigenvals() == {m**2*c**4 - e*hb*c**2*B: 2, m**2*c**4 + e*hb*c**2*B: 2} and Sz*beta == beta*Sz and Sz*Sz == sp.eye(4))
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
