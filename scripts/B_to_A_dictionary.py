"""The link between the ancient Part 5 matrices B1, B2, B3 and the master's A_x, A_y, A_z, T.
Ancient: B1 = [[0,1],[1,0]], B2 = [[1,0],[0,-1]], B3 = [[0,-1],[1,0]]; U^{-1} = (1/sqrt2)[[1,i],[1,-i]].
Master 1+1 (III.8): A = diag(1,-1), T = [[0,1],[-1,0]]; III.8a: m(beta) = cos2b*diag(1,-1) + sin2b*[[0,1],[1,0]], R(alpha).
Master 2+1 (IV.4c/d): A_z = diag(1,-1), A_x = [[0,1],[1,0]], K_z = A_z/2, K_x = A_x/2, J = K_zK_x - K_xK_z.
Master 3+1 real (IV.4b (I)): A_z = diag(1,-1,-1,1), A_x, A_y, T real 4x4.
Every PASS below is one sentence of the chat answer."""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
I2 = sp.eye(2); i = sp.I
B1 = sp.Matrix([[0, 1], [1, 0]]); B2 = sp.diag(1, -1); B3 = sp.Matrix([[0, -1], [1, 0]])
sx, sy, sz = B1, sp.Matrix([[0, -i], [i, 0]]), B2
A1 = sp.diag(1, -1); T = sp.Matrix([[0, 1], [-1, 0]])
m = lambda b: sp.Matrix([[sp.cos(2*b), sp.sin(2*b)], [sp.sin(2*b), -sp.cos(2*b)]])
R = lambda al: sp.Matrix([[sp.cos(al), -sp.sin(al)], [sp.sin(al), sp.cos(al)]])

# ---- 1. in 1+1 and 2+1 the B's ARE the master's matrices ----
ok("1a. B2 = A (III.8 direction matrix) = m(0) = M(x-hat): the half-turn about the in-plane x axis", B2 == A1 and B2 == m(0))
ok("1b. B1 = m(pi/4) = M(y-hat): the half-turn about the in-plane diagonal = III.5a's exchange matrix in the clock = A_x of the 2+1 walk", B1 == m(sp.pi/4))
ok("1c. B3 = -T = B1 B2 = R(+pi/2); T = B2 B1 = m(0) m(pi/4): the quarter-turn is the product of the two half-turns, as III.8a says",
   B3 == -T and B1*B2 == B3 and B2*B1 == T and B3 == R(sp.pi/2))
ok("1d. [B1,B2] = +2 B3 (the old Doc 10 wrote -2B3: that was the slip). {B1,B2,B3} close: [B3,B1] = -2B2, [B3,B2] = 2B1",
   B1*B2 - B2*B1 == 2*B3 and B3*B1 - B1*B3 == -2*B2 and B3*B2 - B2*B3 == 2*B1)
# the 2+1 runner generators of IV.4c/d
Kz, Kx = A1/2, B1/2
Jrun = Kz*Kx - Kx*Kz
ok("1e. IV.4d (F): K_z = B2/2, K_x = B1/2, and the turn they leave, J = K_zK_x - K_xK_z = T/2 = -B3/2. "
   "These are the ancient's (K1, K2, J) = (B2/2, B1/2, B3/2) up to the sign of J: the flatlander's own 2+1 Lorentz algebra so(2,1)",
   Jrun == T/2 and Jrun == -B3/2)
eta, th = sp.symbols('eta theta', real=True)
ok("1f. exp(eta K_z) = diag(e^{eta/2}, e^{-eta/2}): the ancient's 'squeeze of the shadow' e^{eta B2/2} IS the run along z of IV.4d (D), "
   "and exp(eta K_x) = cosh(eta/2) + sinh(eta/2) B1 is the run along x of IV.4d (E). Exponentiating a real B gives a boost, never a rotation",
   sp.simplify((eta*Kz).exp() - sp.diag(sp.exp(eta/2), sp.exp(-eta/2))) == sp.zeros(2)
   and sp.simplify((eta*Kx).exp() - (sp.cosh(eta/2)*I2 + sp.sinh(eta/2)*B1)) == sp.zeros(2))
ok("1g. exp(theta B3) = R(theta): the one real rotation generator; the only real antisymmetric 2x2 (III.8: 'the quarter-turn is forced')",
   sp.simplify((th*B3).exp() - R(th)) == sp.zeros(2))

# ---- 2. the ancient's U is III.9's complex column: a real vector at angle alpha -> (e^{i alpha}, e^{-i alpha})/sqrt2 ----
Uinv = sp.Matrix([[1, i], [1, -i]])/sp.sqrt(2); U = Uinv.inv()
al = sp.symbols('alpha', real=True)
v = sp.Matrix([sp.cos(al), sp.sin(al)])
ok("2a. U^{-1} v(alpha) = (e^{i alpha}, e^{-i alpha})/sqrt2: the real vector becomes the pair (z, zbar), z = x + i y. "
   "The rows of U^{-1} are the eigenvectors of B3 (= -T): this is III.9's complex column e1 + i T e1, on which T is the number -i (polar_matrices.py v)",
   sp.simplify(Uinv*v - sp.Matrix([sp.exp(i*al), sp.exp(-i*al)])/sp.sqrt(2)) == sp.zeros(2, 1)
   and sp.simplify(B3*sp.Matrix([1, -i]) - i*sp.Matrix([1, -i])) == sp.zeros(2, 1))
ok("2b. In that basis the real rotation generator becomes a number: U^{-1} B3 U = i sigma_z = diag(i, -i): B3 z = i z, B3 zbar = -i zbar "
   "(the ancient's 'the i emerged from the coordinate transformation'; III.9's T -> -i on the plane wave)",
   sp.simplify(Uinv*B3*U) == sp.diag(i, -i))
ok("2c. ... and the two half-turns become off-diagonal: U^{-1} B2 U = sigma_x, U^{-1} B1 U = -sigma_y: they send z to a multiple of zbar, "
   "i.e. they are antilinear in z (the ancient's 'B1 z = i zbar ... involves conjugation'). A half-turn anticommutes with the quarter-turn, so it conjugates the i",
   sp.simplify(Uinv*B2*U - sx) == sp.zeros(2) and sp.simplify(Uinv*B1*U + sy) == sp.zeros(2)
   and B2*B3 == -B3*B2 and B1*B3 == -B3*B1)
ok("2d. The ancient's three so(3) generators {B3, iB2, -iB1} = U{i sigma_z, i sigma_x, i sigma_y}U^{-1}: so U just renames the axes, "
   "(x,y,z) -> (z,-x,-y) [U sigma_x U^{-1} = sigma_z, U sigma_y U^{-1} = -sigma_x, U sigma_z U^{-1} = -sigma_y], a rotation of space (det +1)",
   sp.simplify(U*sx*Uinv - sz) == sp.zeros(2) and sp.simplify(U*sy*Uinv + sx) == sp.zeros(2) and sp.simplify(U*sz*Uinv + sy) == sp.zeros(2)
   and sp.Matrix([[0, 0, 1], [-1, 0, 0], [0, -1, 0]]).det() == 1)

# ---- 3. in 3+1 (IV.4b real basis) the ancient's i is the second pair of counters ----
Az = sp.diag(1, -1, -1, 1)
Ax = sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]])
Ay = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
T4 = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]])
kron = sp.kronecker_product
Jb = sp.Matrix([[0, -1], [1, 0]])          # the real 'i' acting between the pairs (1,2) and (3,4)
ok("3a. IV.4b's real matrices in blocks: A_z = diag(B2, -B2), A_x = diag(B1, -B1), A_y = sigma_x (x) 1, T = -J (x) 1 with J = [[0,-1],[1,0]] between the pairs",
   Az == kron(sp.diag(1, -1), B2) and Ax == kron(sp.diag(1, -1), B1) and Ay == kron(sx, I2) and T4 == kron(-Jb, I2))
ok("3b. The three spin generators (IV.4c (K): (1/2)A_iA_j) are the ancient's three, with i -> J between the pairs: "
   "A_zA_x = 1 (x) (-B3) [real, the in-plane quarter-turn on each pair: the ancient's real B3]; "
   "A_yA_z = J (x) B2 [the ancient's i B2]; A_xA_y = (-J) (x) B1 [the ancient's -i B1]",
   Az*Ax == kron(I2, -B3) and Ay*Az == kron(Jb, B2) and Ax*Ay == kron(-Jb, B1))
ok("3c. exp(theta A_yA_z) (v, 0) = (cos theta v, sin theta B2 v): the ancient's e^{theta iB2} v = cos v + i sin B2 v, with 'the imaginary part' = the second pair of counters. "
   "Tilting the coin moves the reading from one pair into the other, reflected",
   (lambda M: sp.simplify(M*sp.Matrix([sp.cos(al), sp.sin(al), 0, 0]) - sp.Matrix([sp.cos(th)*sp.cos(al), sp.cos(th)*sp.sin(al), sp.sin(th)*sp.cos(al), -sp.sin(th)*sp.sin(al)])) == sp.zeros(4, 1))
   ((th*Ay*Az).exp().applyfunc(sp.simplify)))
# ---- 4. which 'i' commutes with what ----
Jw = Ax*Ay*Az
ok("4a. A_xA_yA_z squares to -1 and commutes with A_x, A_y, A_z: it is the i of the massless real 4-column, which is therefore one complex Weyl column (IV.3) exactly",
   sp.simplify(Jw*Jw + sp.eye(4)) == sp.zeros(4) and all(Jw*M == M*Jw for M in (Ax, Ay, Az)))
ok("4b. T anticommutes with that i: the mass is antilinear for the Weyl reading (a Majorana mass couples the column to its conjugate). "
   "So the real 4-column with mass is NOT a complex 2-column; V.2's statement 'real Dirac = neutral (Majorana), charged = complex'",
   T4*Jw == -Jw*T4)
ok("4c. A_xA_y squares to -1 and commutes with A_z and T: it is the i of the massive column moving along z (the 1+1 complex column of III.9 read inside the real 3+1 walker, "
   "i = the turn about the axis of motion; i_is_spin_about_motion.py). It does not commute with A_x, A_y: there is no i for the massive column in all directions at once",
   sp.simplify(Ax*Ay*Ax*Ay + sp.eye(4)) == sp.zeros(4) and Ax*Ay*Az == Az*Ax*Ay and Ax*Ay*T4 == T4*Ax*Ay and Ax*Ay*Ax != Ax*Ax*Ay)
# the commutant of the full set is R: no J commuting with all of A_x,A_y,A_z,T
X = sp.Matrix(4, 4, sp.symbols('x0:16'))
eqs = []
for M in (Ax, Ay, Az, T4):
    eqs += list(X*M - M*X)
sol = sp.solve(eqs, list(X), dict=True)
ok("4d. Nothing real commutes with all four of A_x, A_y, A_z, T except multiples of 1: the charged particle's i (III.9, V.2) must be a second copy, not a matrix on these four counters",
   len(sol) == 1 and (X.subs(sol[0]) - X.subs(sol[0])[0, 0]*sp.eye(4)) == sp.zeros(4))

# ---- 5. the orbit as a vortex: half a quantum of circulation ----
hbar, mm, c = sp.symbols('hbar m c', positive=True)
r = hbar/(2*mm*c); vv = c
ok("5. The one-copy orbit of true_spin.py 3b: radius hbar/2mc at speed c -> angular momentum m r c = hbar/2, circulation 2 pi r c = (h/m)/2: "
   "the walker's own loop is the half-quantum of the old vortex picture",
   sp.simplify(mm*r*vv - hbar/2) == 0 and sp.simplify(2*sp.pi*r*vv - (2*sp.pi*hbar/mm)/2) == 0)
