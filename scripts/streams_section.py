"""Checks for the draft IV.4d (streams), 2+1, t-reading, a = 0 unless said: d_t Q = -c A_z d_z Q - c A_x d_x Q (+ a T Q).
(A) the four streams solve the lines;  (B) the content identity;  (C) the count splits two ways;  (D) half-rapidity from the content;
(E) the x-run from universality; (F) the loop of two runs = the U-turn coin;  (G) one stream along theta = forward + backward along z,
    with the cross term = the x-flow;  (H) with a != 0: the flips move content between the z-streams: d(Q+^2) = 2aQ+Q-, d(Q-^2) = -2aQ+Q-."""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
t, x, z, eta, th = sp.symbols('t x z eta theta', real=True); c, a = sp.symbols('c a', positive=True)
Az = sp.diag(1, -1); Ax = sp.Matrix([[0, 1], [1, 0]]); T = sp.Matrix([[0, 1], [-1, 0]]); f = sp.Function('f')
L0 = lambda Q: Q.diff(t) + c*Az*Q.diff(z) + c*Ax*Q.diff(x)
S = {'+z': (sp.Matrix([1, 0]), t - z/c), '-z': (sp.Matrix([0, 1]), t + z/c), '+x': (sp.Matrix([1, 1]), t - x/c), '-x': (sp.Matrix([1, -1]), t + x/c)}
ok("(A) (1,0)f(t-z/c), (0,1)f(t+z/c), (1,1)f(t-x/c), (1,-1)f(t+x/c) solve the lines", all(sp.simplify(L0(v*f(s))) == sp.zeros(2, 1) for v, s in S.values()))
Qp, Qm = sp.Function('Qp')(t, x, z), sp.Function('Qm')(t, x, z); Q = sp.Matrix([Qp, Qm])
sub = {Qp.diff(t): -c*Qp.diff(z) - c*Qm.diff(x), Qm.diff(t): c*Qm.diff(z) - c*Qp.diff(x)}
ident = sp.expand((sp.diff(Qp**2 + Qm**2, t)).subs(sub) + c*sp.diff(Qp**2 - Qm**2, z) + c*sp.diff(2*Qp*Qm, x))
ok("(B) d_t(Q+^2 + Q-^2) + c d_z(Q+^2 - Q-^2) + c d_x(2 Q+ Q-) = 0", sp.simplify(ident) == 0)
u, w = sp.symbols('u w', real=True)
ok("(C) Q+^2 + Q-^2 = (1/2)(Q+ + Q-)^2 + (1/2)(Q+ - Q-)^2", sp.expand((u + w)**2/2 + (u - w)**2/2 - u**2 - w**2) == 0)
j = lambda q: sp.Matrix([q[0]**2 + q[1]**2, 2*q[0]*q[1], q[0]**2 - q[1]**2])        # (ct, x, z)
Lz = sp.Matrix([[sp.cosh(eta), 0, sp.sinh(eta)], [0, 1, 0], [sp.sinh(eta), 0, sp.cosh(eta)]])
qb = sp.Matrix([sp.exp(eta/2)*u, sp.exp(-eta/2)*w])
ok("(D) readings scaled by e^{+-eta/2} give a content scaled as a vector: j(S Q) = L_z j(Q)", (j(qb) - Lz*j(sp.Matrix([u, w]))).applyfunc(lambda e: sp.simplify(e.rewrite(sp.exp))) == sp.zeros(3, 1))
Sx = sp.Matrix([[sp.cosh(eta/2), sp.sinh(eta/2)], [sp.sinh(eta/2), sp.cosh(eta/2)]])
Lx = sp.Matrix([[sp.cosh(eta), sp.sinh(eta), 0], [sp.sinh(eta), sp.cosh(eta), 0], [0, 0, 1]])
ok("(E) the x-run: (Q+ + Q-) -> e^{eta/2}(Q+ + Q-), (Q+ - Q-) -> e^{-eta/2}(Q+ - Q-), i.e. the cosh/sinh matrix; content scales as a vector along x",
   (j(Sx*sp.Matrix([u, w])) - Lx*j(sp.Matrix([u, w]))).applyfunc(lambda e: sp.simplify(e.rewrite(sp.exp))) == sp.zeros(3, 1))
Kz, Kx = Az/2, Ax/2; J = Kz*Kx - Kx*Kz
ok("(F) [K_z, K_x] = J = [[0,1/2],[-1/2,0]]: 2J is the z-clock U-turn coin, +1/2 into +, -1/2 into -", J == sp.Matrix([[0, sp.Rational(1, 2)], [-sp.Rational(1, 2), 0]]))
half = sp.Matrix([sp.cos(th/2), sp.sin(th/2)]); A_th = sp.cos(th)*Az + sp.sin(th)*Ax
ok("(G) (cos th/2, sin th/2) is the +c stream along theta; sorted along z it is cos^2(th/2) forward + sin^2(th/2) backward, cross term 2Q+Q- = sin(theta) = its x-flow",
   sp.simplify(A_th*half - half) == sp.zeros(2, 1) and sp.simplify(2*half[0]*half[1] - sp.sin(th)) == 0)
subm = {Qp.diff(t): -c*Qp.diff(z) + a*Qm, Qm.diff(t): c*Qm.diff(z) - a*Qp}              # 1+1 with flips
ok("(H) with flips: (d_t + c d_z) Q+^2 = 2a Q+Q-, (d_t - c d_z) Q-^2 = -2a Q+Q-: the flips move content from one z-stream to the other",
   sp.simplify(sp.expand((sp.diff(Qp**2, t) + c*sp.diff(Qp**2, z)).subs(subm)) - 2*a*Qp*Qm) == 0
   and sp.simplify(sp.expand((sp.diff(Qm**2, t) - c*sp.diff(Qm**2, z)).subs(subm)) + 2*a*Qp*Qm) == 0)
