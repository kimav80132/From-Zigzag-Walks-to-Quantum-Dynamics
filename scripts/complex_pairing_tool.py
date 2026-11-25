"""The complex pairing as a tool (draft III.8b). One rule, three uses.
Rule: a real space with a matrix J, J^2 = -1, is a complex space with i = J. Eigenvectors of J give the (z, zbar) basis.
A real matrix M commuting with J is complex-linear in z; one anticommuting with J is antilinear (sends z to a multiple of zbar).
Canonical layer: a canonical pair (q,p), {q,p} = 1, gives z = (q + i p)/sqrt2 with {z, zbar} = -i; Madelung is the polar form with hbar.
Use 1 (I.9): (rho, S) -> psi.  Use 2 (III.9): (Q, Q') -> psi, J the copy turn, commuting with A and T.  Use 3 (Part 5): (x, y) -> z with J = B3."""
import sympy as sp
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
i = sp.I
# --- canonical layer ---
q, p, rho, S, hbar = sp.symbols('q p rho S hbar', real=True); hbar = sp.symbols('hbar', positive=True)
pb = lambda f, g, x, y: sp.diff(f, x)*sp.diff(g, y) - sp.diff(f, y)*sp.diff(g, x)      # {f,g} with {x,y} = 1
z = (q + i*p)/sp.sqrt(2); zb = (q - i*p)/sp.sqrt(2)
ok("C1. {q,p} = 1  ->  z = (q + i p)/sqrt2 has {z, zbar} = -i: the i of quantum mechanics is the bracket of a pair read as one complex number",
   sp.simplify(pb(z, zb, q, p) + i) == 0)
psi = sp.sqrt(rho)*sp.exp(i*S/hbar); psib = sp.sqrt(rho)*sp.exp(-i*S/hbar)
ok("C2. the same with {rho,S} = 1 and the polar form psi = sqrt(rho) e^{iS/hbar}: {psi, psi*} = -i/hbar (I.9). Madelung is use 1 of the rule",
   sp.simplify(pb(psi, psib, rho, S) + i/hbar) == 0)
# --- linear / antilinear theorem, 2x2 ---
B1 = sp.Matrix([[0, 1], [1, 0]]); B2 = sp.diag(1, -1); B3 = sp.Matrix([[0, -1], [1, 0]])
Uinv = sp.Matrix([[1, i], [1, -i]])/sp.sqrt(2); U = Uinv.inv()        # rows of U^{-1}: the two eigenvectors of B3
x, y = sp.symbols('x y', real=True); v = sp.Matrix([x, y])
zz = (Uinv*v)[0]                                                        # = (x + i y)/sqrt2
ok("P1. the pairing: U^{-1} (x,y) = ((x+iy), (x-iy))/sqrt2 = (z, zbar)/sqrt2 up to the factor; J = B3 becomes i on z (and -i on zbar)",
   sp.simplify(zz - (x + i*y)/sp.sqrt(2)) == 0 and sp.simplify(Uinv*B3*U - sp.diag(i, -i)) == sp.zeros(2))
rot = lambda th: sp.cos(th)*sp.eye(2) + sp.sin(th)*B3
th = sp.symbols('theta', real=True)
ok("P2. commuting with J -> complex-linear: e^{theta B3} commutes with B3 and on z it is the number e^{i theta}",
   sp.simplify(rot(th)*B3 - B3*rot(th)) == sp.zeros(2) and sp.simplify(sp.expand((Uinv*rot(th)*v)[0] - sp.exp(i*th)*zz).rewrite(sp.cos)) == 0)
ok("P3. anticommuting with J -> antilinear: B2 z = zbar, B1 z = i zbar (the ancient's 'B1 z = i zbar ... involves conjugation')",
   B2*B3 == -B3*B2 and B1*B3 == -B3*B1
   and sp.simplify((Uinv*B2*v)[0] - sp.conjugate(zz)) == 0 and sp.simplify((Uinv*B1*v)[0] - i*sp.conjugate(zz)) == 0)
# --- the rule for a general real matrix (added 10 Nov 2026: P2, P3 above test it on examples only) ---
M = sp.Matrix(2, 2, sp.symbols('m0:4', real=True)); al, be = sp.symbols('alpha beta', real=True)
com = sp.solve(list(M*B3 - B3*M), list(M), dict=True)[0]; anti = sp.solve(list(M*B3 + B3*M), list(M), dict=True)[0]
Mc, Ma = M.subs(com), M.subs(anti)
ok("P2'. every real M with MJ = JM is alpha 1 + beta J, and on z it is the complex number alpha + i beta: linear",
   sp.simplify(Mc - (Mc[0, 0]*sp.eye(2) + Mc[1, 0]*B3)) == sp.zeros(2) and sp.simplify(sp.expand((Uinv*Mc*v)[0] - (Mc[0, 0] + i*Mc[1, 0])*zz)) == 0)
ok("P3'. every real M with MJ = -JM is alpha B2 + beta B1, and it sends z to (alpha + i beta) zbar: antilinear; and any real M is the sum of one of each",
   sp.simplify(Ma - (Ma[0, 0]*B2 + Ma[0, 1]*B1)) == sp.zeros(2) and sp.simplify(sp.expand((Uinv*Ma*v)[0] - (Ma[0, 0] + i*Ma[0, 1])*sp.conjugate(zz))) == 0
   and sp.simplify((M - B3*M*B3)/2*B3 - B3*(M - B3*M*B3)/2) == sp.zeros(2) and sp.simplify((M + B3*M*B3)/2*B3 + B3*(M + B3*M*B3)/2) == sp.zeros(2))
# --- the doubling is the pairing on the doubled space ---
Jc = sp.Matrix([[0, -1], [1, 0]])
A = sp.diag(1, -1); T = sp.Matrix([[0, 1], [-1, 0]])
kron = sp.kronecker_product
J4 = kron(Jc, sp.eye(2))              # the copy turn of V.2 on (Q, Q')
ok("D1. use 2 is the rule on the doubled space (Q,Q') with J = the copy turn: J^2 = -1 and J commutes with A and T, so A and T stay complex-linear "
   "and the complex column obeys the same Dirac equation (III.9, V.2 (F))",
   J4*J4 == -sp.eye(4) and all(J4*kron(sp.eye(2), M) == kron(sp.eye(2), M)*J4 for M in (A, T)))
ok("D2. on the single real column, T is the only quarter-turn and A anticommutes with it: pairing (Q+,Q-) with i = T makes A antilinear, "
   "so the 1+1 column cannot be read as one complex number with A linear. That is why III.9's i is a second copy and not T",
   A*T == -T*A)
# D2 for every J, not only T (added 10 Nov 2026): what commutes with A is diagonal, and a real diagonal matrix does not square to -1
Jg = sp.Matrix(2, 2, sp.symbols('j0:4', real=True))
comA = sp.solve(list(Jg*A - A*Jg), list(Jg), dict=True)[0]; Jd = Jg.subs(comA)
orth = sp.solve(list(Jg*Jg + sp.eye(2)) + list(Jg.T*Jg - sp.eye(2)), list(Jg), dict=True)
ok("D2'. a real matrix commuting with A = diag(1,-1) is diagonal, and its square has entries j0^2, j3^2 >= 0, never -1: no real J at all makes the "
   "single column complex with A linear. The real J with J^2 = -1 that are also turns (J^T J = 1) are +T and -T, nothing else",
   Jd[0, 1] == 0 and Jd[1, 0] == 0 and (Jd*Jd)[0, 0] == Jd[0, 0]**2 and (Jd*Jd)[1, 1] == Jd[1, 1]**2
   and sorted([tuple(Jg.subs(o_)) for o_ in orth]) == sorted([tuple(T), tuple(-T)]))
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
