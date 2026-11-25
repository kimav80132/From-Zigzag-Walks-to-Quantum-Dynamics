"""Can the second copy be the same walk going back?  Test on the plane wave, where copy 2 must be copy 1 a quarter clock later.
Real walk d_t Q = -cA d_x Q + aT Q.  Reversals that map solutions to solutions, Qr(t,x) = M Q(reversed point):
 (1) lab time reversed, t -> -t:          needs M A = -A M, M T = -T M  ->  M = Ax (swap the senses).
 (2) worldline retraced: the walker's own proper time reversed, tau -> -tau at fixed comoving position sigma:
     this is (1) conjugated by the run B = diag(e^{eta/2}, e^{-eta/2}): M = B Ax B^-1 = B^2 Ax.
Plane wave with rapidity eta: copy1 = B exp(a tau T) e1, required copy2 = B exp(a tau T) T e1.
 (3) at rest, (1) gives -copy2.
 (4) moving, (1) gives -(copy 2 of the wave going the OTHER way): the bridge's failure on moving packets.
 (5) moving, (2) gives -copy2 exactly, for every eta: the second copy IS the walker retracing its own worldline, senses swapped."""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
t, x, a, c, eta = sp.symbols('t x a c eta', real=True)
A = sp.diag(1, -1); T = sp.Matrix([[0, 1], [-1, 0]]); Ax = sp.Matrix([[0, 1], [1, 0]]); I2 = sp.eye(2); e1 = sp.Matrix([1, 0])
m = sp.Matrix(2, 2, sp.symbols('m0:4', real=True))
sol = sp.solve(list(m*A + A*m) + list(m*T + T*m), list(m), dict=True)
ok("1. lab time reversal t -> -t maps solutions to solutions with M = p Ax: the senses swap", m.subs(sol[0]) == m.subs(sol[0])[0, 1]*Ax)
B = sp.diag(sp.exp(eta/2), sp.exp(-eta/2))
tau = t*sp.cosh(eta) - x*sp.sinh(eta)/c; sig = x*sp.cosh(eta) - c*t*sp.sinh(eta)
E = lambda s: sp.cos(a*s)*I2 + sp.sin(a*s)*T
copy1 = B*E(tau)*e1; copy2 = B*E(tau)*T*e1
simp = lambda M: M.applyfunc(lambda f: sp.simplify(f.rewrite(sp.exp)))
# (3) at rest, lab reversal
rest1 = copy1.subs(eta, 0); rest2 = copy2.subs(eta, 0)
ok("3. at rest: Ax * copy1(-t) = -copy2: the returning walker IS the quarter-clock partner", simp(Ax*rest1.subs(t, -t) + rest2) == sp.zeros(2, 1))
# (4) moving, lab reversal: compare with copy2 of the wave with -eta
copy2_back = copy2.subs(eta, -eta)
ok("4. moving: Ax * copy1(-t,x) = -(copy 2 of the wave with rapidity -eta): lab-time reversal pairs the walker with the wave going the other way",
   simp(Ax*copy1.subs(t, -t) + copy2_back) == sp.zeros(2, 1))
# (5) worldline retrace: tau -> -tau, sigma fixed. Find the lab point (t', x') with tau' = -tau, sigma' = sigma.
tp, xp = sp.symbols('tp xp', real=True)
solpt = sp.solve([tp*sp.cosh(eta) - xp*sp.sinh(eta)/c + tau, xp*sp.cosh(eta) - c*tp*sp.sinh(eta) - sig], [tp, xp], dict=True)[0]
M2 = B**2*Ax
retrace = M2*copy1.subs({t: solpt[tp], x: solpt[xp]}, simultaneous=True)
ok("5. moving: (B^2 Ax) * copy1(proper time reversed, comoving position fixed) = -copy2 for every eta: the second copy is the walker retracing its own worldline with senses swapped",
   simp(retrace + copy2) == sp.zeros(2, 1))
# and the retrace map is itself a symmetry of the real walk (conjugate of (1) by the run)
ok("5'. B^2 Ax = B Ax B^-1: the retrace matrix is the sense swap seen by the runner", sp.simplify(B**2*Ax - B*Ax*B.inv()) == sp.zeros(2))
