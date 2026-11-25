"""GR-style chain for the counters, 2+1, t-reading (A_z = diag(1,-1), A_x = [[0,1],[1,0]], T = [[0,1],[-1,0]]).
Element: the +z stream of the 1+1 case, reading psi0 = (1, 0).  Its content (n, j_x, j_z) = (Q1^2+Q2^2, 2Q1Q2, Q1^2-Q2^2) = (1, 0, 1).
Boosts make streams in other directions from it (as boosting u gives every u_A in GR); the counters follow with
S_z = exp(eta A_z/2), S_x = exp(eta A_x/2) (boosts_streams.py).  Universality: each stream obeys the 1+1 law along its own direction.
Additivity: the readings of coexisting streams add (the lines are linear).
Angle phi of a direction measured from +z towards +x:  n = (sin phi, cos phi) = (n_x, n_z).
"""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
eta, e1, e2, phi, th = sp.symbols('eta eta1 eta2 phi theta', real=True)
Az = sp.diag(1, -1); Ax = sp.Matrix([[0, 1], [1, 0]]); I2 = sp.eye(2)
Kz, Kx = Az/2, Ax/2; J = Kz*Kx - Kx*Kz
Sz = lambda e: sp.diag(sp.exp(e/2), sp.exp(-e/2))
Sx = lambda e: sp.Matrix([[sp.cosh(e/2), sp.sinh(e/2)], [sp.sinh(e/2), sp.cosh(e/2)]])
content = lambda Q: sp.Matrix([Q[0]**2 + Q[1]**2, 2*Q[0]*Q[1], Q[0]**2 - Q[1]**2])
def Lb(axis, e):
    L = sp.eye(3); i = 1 if axis == 'x' else 2
    L[0, 0] = L[i, i] = sp.cosh(e); L[0, i] = L[i, 0] = sp.sinh(e); return L
simp = lambda M: M.applyfunc(lambda f: sp.simplify(sp.expand(f.rewrite(sp.exp))))
psi0 = sp.Matrix([1, 0]); k0 = sp.Matrix([1, 0, 1])
psi = lambda p: sp.Matrix([sp.cos(p/2), sp.sin(p/2)])

# 1. boosting the +z stream along x: readings and content
Q = Sx(eta)*psi0
ok(f"1a. S_x(eta) psi0 = {list(Q)}", Q == sp.Matrix([sp.cosh(eta/2), sp.sinh(eta/2)]))
ok("1b. its content = L_x(eta) (1,0,1) = (cosh eta, sinh eta, 1): the boosted stream, tilted towards +x (aberration)",
   simp(content(Q) - Lb('x', eta)*k0) == sp.zeros(3, 1))
# 2. tilt angle: tan phi = sinh eta ; readings proportional to psi(phi), i.e. tan(phi/2) = tanh(eta/2)
t = sp.tanh(eta/2)
ok("2. tan(phi/2) = tanh(eta/2) gives tan(phi) = 2t/(1-t^2) = sinh(eta): the readings are (cos phi/2, sin phi/2) times sqrt(cosh eta)",
   sp.simplify((2*t/(1 - t**2) - sp.sinh(eta)).rewrite(sp.exp)) == 0)
# 3. universality: the stream along phi is the +c mover of A(phi)
A = lambda p: sp.cos(p)*Az + sp.sin(p)*Ax
ok("3a. A(phi) psi(phi) = +psi(phi): psi(phi) moves at +c along phi", sp.simplify(A(phi)*psi(phi) - psi(phi)) == sp.zeros(2, 1))
ok("3b. A(phi) psi(phi+pi) = -psi(phi+pi): the reverse stream moves at -c along phi", sp.simplify(A(phi)*psi(phi + sp.pi) + psi(phi + sp.pi)) == sp.zeros(2, 1))
ok("3c. A(phi) = psi psi^T - psi' psi'^T (built from the two streams), and it equals cos phi A_z + sin phi A_x",
   sp.simplify(psi(phi)*psi(phi).T - psi(phi + sp.pi)*psi(phi + sp.pi).T - A(phi)) == sp.zeros(2))
ok("3d. content of psi(phi) = (1, sin phi, cos phi): one stream of count 1 moving along phi",
   sp.simplify(content(psi(phi)) - sp.Matrix([1, sp.sin(phi), sp.cos(phi)])) == sp.zeros(3, 1))
# 4. full turn
ok("4. psi(phi + 2 pi) = -psi(phi), content unchanged", sp.simplify(psi(phi + 2*sp.pi) + psi(phi)) == sp.zeros(2, 1))
# 5. two boost routes from +z to -z
r_plus = Sz(-e2)*Sx(e1)*psi0; r_minus = Sz(-e2)*Sx(-e1)*psi0
lim = lambda r: sp.Matrix([sp.limit(sp.simplify(r[i]/sp.sqrt(r[0]**2 + r[1]**2)), e2, sp.oo) for i in range(2)])
lp = lim(r_plus.subs(e1, 1)); lm = lim(r_minus.subs(e1, 1))
ok(f"5a. route boost +x then boost -z (repeatedly): reading -> {list(lp)};  route -x then -z: -> {list(lm)}",
   lp == sp.Matrix([0, 1]) and lm == sp.Matrix([0, -1]))
ok("5b. both routes end on the same stream (-z): contents both -> (1, 0, -1)",
   simp(content(lp) - sp.Matrix([1, 0, -1])) == sp.zeros(3, 1) and simp(content(lm) - sp.Matrix([1, 0, -1])) == sp.zeros(3, 1))
# 6. half-turns
Uccw = sp.Matrix([[0, -1], [1, 0]]); Ucw = sp.Matrix([[0, 1], [-1, 0]])
ok("6a. turning phi -> phi + pi (via +x) acts as [[0,-1],[1,0]]; phi -> phi - pi (via -x) acts as [[0,1],[-1,0]]",
   sp.simplify(psi(phi + sp.pi) - Uccw*psi(phi)) == sp.zeros(2, 1) and sp.simplify(psi(phi - sp.pi) - Ucw*psi(phi)) == sp.zeros(2, 1))
ok("6b. the two half-turns differ by -1, and they are -2J and +2J (2J = 2[K_z,K_x])", Uccw == -Ucw and Ucw == 2*J)
ok("6c. turning by any theta: psi(phi - theta) = exp(theta J) psi(phi)", sp.simplify(psi(phi - th) - sp.simplify((th*J).exp())*psi(phi)) == sp.zeros(2, 1))
ok("6d. the order of the boosts sets the sense: [K_z,K_x] = J, [K_x,K_z] = -J", Kx*Kz - Kz*Kx == -J)
# 7. additivity of readings
x, z, ct = sp.symbols('x z ct', real=True); f, g = sp.Function('f'), sp.Function('g')
p1, p2 = sp.symbols('phi1 phi2', real=True)
stream = lambda p, F: psi(p)*F(ct - sp.sin(p)*x - sp.cos(p)*z)
L0 = lambda Qv: Qv.diff(ct) + Ax*Qv.diff(x) + Az*Qv.diff(z)           # the lines with a = 0, c t as time
ok("7a. one stream psi(phi) f(ct - n.r) solves the lines (a = 0), for every phi", sp.simplify(L0(stream(p1, f))) == sp.zeros(2, 1))
ok("7b. two streams added solve them too: readings add", sp.simplify(L0(stream(p1, f) + stream(p2, g))) == sp.zeros(2, 1))
u, w = sp.symbols('u w', real=True)
cross = sp.simplify(content(u*psi(p1) + w*psi(p2)) - content(u*psi(p1)) - content(w*psi(p2)))
ok(f"7c. but contents do not add: cross term {list(cross)}", cross != sp.zeros(3, 1))
ok("7d. cross term of a stream with the same stream turned by 2 pi has the opposite sign",
   sp.simplify(content(u*psi(p1) + w*psi(p1 + 2*sp.pi)) - content((u - w)*psi(p1))) == sp.zeros(3, 1))
