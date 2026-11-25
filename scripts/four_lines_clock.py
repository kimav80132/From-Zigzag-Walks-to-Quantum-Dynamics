"""IV.4b — the four lines in the clock z, line by line, in the labels (P1,P2,P3,P4).
(A)  c dz P = [advance + loss + exchange + x-turn dx + y-turn dy] P
(F)  loss removed, Q = e^{az/c} P
(G-I) solved for dt:  dt Q = -c(A_x dx + A_y dy + A_z dz) Q + a T Q
(J)  left-multiplication by the direction matrix A_z = diag(1,-1,-1,1) sends  exchange -> T,  x-turn -> -A_x,  y-turn -> -A_y
(K)  T^2 = -1, A_j^2 = 1, all six pairs anticommute  => Klein-Gordon with -a^2 in every component
plus: the real A_y is IV.4's imaginary A_y after a unitary change of basis (same algebra, Pauli), and
      with dx = dy = 0 the four lines are III.2 twice.
"""
import numpy as np, sympy as sp
ok_all = True
def ok(name, cond):
    global ok_all; ok_all &= bool(cond); print(("PASS  " if cond else "FAIL  ") + name)
I4 = np.eye(4)
# (A) the five matrices, read off the four lines
adv = np.diag([-1, 1, 1, -1.])                                   # coefficient of d_t
exch = np.array([[0, 0, 1, 0], [0, 0, 0, -1], [1, 0, 0, 0], [0, -1, 0, 0.]])      # times a
xt = np.array([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0.]])        # times c d_x
yt = np.array([[0, 0, -1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, -1, 0, 0.]])        # times c d_y
ok("A. exchange symmetric; x-turn and y-turn antisymmetric", np.allclose(exch, exch.T) and np.allclose(xt, -xt.T) and np.allclose(yt, -yt.T))
ok("A. exchange^2 = 1 (a flip), x-turn^2 = y-turn^2 = -1 (turns), already in the clock",
   np.allclose(exch @ exch, I4) and np.allclose(xt @ xt, -I4) and np.allclose(yt @ yt, -I4))
# the four lines, symbolically, against the user's system
t, x, y, z, c, a = sp.symbols('t x y z c a', real=True)
P = sp.Matrix([sp.Function('P%d' % k, real=True)(t, x, y, z) for k in range(1, 5)])
S = lambda M: sp.Matrix(M.astype(int))
rhsA = S(adv) * P.diff(t) - a * P + a * S(exch) * P + c * S(xt) * P.diff(x) + c * S(yt) * P.diff(y)
P1, P2, P3, P4 = P
lines = sp.Matrix([-P1.diff(t) - a*P1 + a*P3 - c*P2.diff(x) - c*P3.diff(y),
                   +P2.diff(t) - a*P2 - a*P4 + c*P1.diff(x) + c*P4.diff(y),
                   +P3.diff(t) - a*P3 + a*P1 - c*P4.diff(x) + c*P1.diff(y),
                   -P4.diff(t) - a*P4 - a*P2 + c*P3.diff(x) - c*P2.diff(y)])
ok("A. the matrix form reproduces the four lines", sp.simplify(rhsA - lines) == sp.zeros(4, 1))
# (B)-(F): loss removal; (E) needs no d_z inside the bracket -- true by inspection; check the identity (C)
Pc = sp.Matrix([sp.Function('p%d' % k)(z) for k in range(4)])
ok("C. e^{-az/c} c d_z(e^{az/c} P) = a P + c d_z P", sp.simplify(sp.exp(-a*z/c) * c * (sp.exp(a*z/c) * Pc).diff(z) - (a*Pc + c*Pc.diff(z))) == sp.zeros(4, 1))
# (G)-(I): solve for d_t = left-multiply by adv^{-1} = adv
D = -adv                                                         # the direction matrix diag(1,-1,-1,1): minus the coefficient of d_t
Az, T, Ax, Ay = D, D @ exch, -D @ xt, -D @ yt
ok("J. direction matrix times exchange = T (antisymmetric)", np.allclose(T, np.array([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0.]])) and np.allclose(T, -T.T))
ok("J. direction matrix times x-turn = -A_x, times y-turn = -A_y, both symmetric",
   np.allclose(Ax, np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0.]])) and np.allclose(Ay, np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0.]]))
   and np.allclose(Ax, Ax.T) and np.allclose(Ay, Ay.T))
anti = lambda X, Y: X @ Y + Y @ X
ok("K. T^2 = -1; A_x^2 = A_y^2 = A_z^2 = 1; all six pairs anticommute",
   np.allclose(T @ T, -I4) and all(np.allclose(M @ M, I4) for M in (Ax, Ay, Az))
   and all(np.allclose(anti(U, V), 0) for U, V in [(Ax, Ay), (Ax, Az), (Ay, Az), (Ax, T), (Ay, T), (Az, T)]))
# KG symbolically from the t-form
Q = sp.Matrix([sp.Function('Q%d' % k, real=True)(t, x, y, z) for k in range(1, 5)])
L = lambda g: -c * (S(Ax) * g.diff(x) + S(Ay) * g.diff(y) + S(Az) * g.diff(z)) + a * S(T) * g
ok("K. L(L Q) = c^2 Laplacian Q - a^2 Q, every component", sp.simplify(L(L(Q)) - (c**2 * (Q.diff(x, 2) + Q.diff(y, 2) + Q.diff(z, 2)) - a**2 * Q)) == sp.zeros(4, 1))
# the real A_y vs IV.4's imaginary one: same A_x, A_z, T; a unitary U with U Ay_IV4 U^dag = Ay_real commuting with the other three
Ay4 = np.array([[0, -1j, 0, 0], [1j, 0, 0, 0], [0, 0, 0, 1j], [0, 0, -1j, 0]])
C = Ay @ Ay4                                                     # C^2 = -1, commutes with A_x, A_z, T, anticommutes with both A_y's
U = (I4 + C) / np.sqrt(2)
ok("IV.4's imaginary A_y and the real A_y are related by the unitary U = (1 + A_y A_y^{IV.4})/sqrt2, which leaves A_x, A_z, T alone",
   np.allclose(U @ U.conj().T, I4) and np.allclose(U @ Ay4 @ U.conj().T, Ay)
   and all(np.allclose(U @ M @ U.conj().T, M) for M in (Ax, Az, T)))
# with dx = dy = 0: two copies of III.2 (pairs (1,3) and (2,4); the second with the sign P4 -> -P4)
ok("dx = dy = 0: lines (1),(3) are III.2 along z for (P1,P3); lines (2),(4) are III.2 for (P2,-P4)",
   sp.simplify(lines.subs({P2.diff(x): 0, P3.diff(y): 0, P1.diff(x): 0, P4.diff(y): 0, P4.diff(x): 0, P1.diff(y): 0, P3.diff(x): 0, P2.diff(y): 0})
               - sp.Matrix([-P1.diff(t) - a*P1 + a*P3, P2.diff(t) - a*P2 + a*(-P4), P3.diff(t) - a*P3 + a*P1, -P4.diff(t) - a*P4 - a*P2])) == sp.zeros(4, 1))

# ---- the copy index is the framing: sideways step in direction phi = (1/2) keep + (1/2) m(phi/2) on (Q1,Q4) -> (Q2,Q3)
m = lambda b: np.array([[np.cos(2*b), np.sin(2*b)], [np.sin(2*b), -np.cos(2*b)]])
good = True
for phi in np.linspace(0, 2*np.pi, 25):
    Pc = (I4 + np.cos(phi)*Ax + np.sin(phi)*Ay) / 2
    change = np.array([[Pc[1, 0], Pc[1, 3]], [Pc[2, 0], Pc[2, 3]]]); keep = np.array([[Pc[0, 0], Pc[0, 3]], [Pc[3, 0], Pc[3, 3]]])
    good &= np.allclose(change, m(phi/2)/2) and np.allclose(keep, np.eye(2)/2)
ok("framing: sideways coin in direction phi = (1/2)(keep) + (1/2) m(phi/2) on the framing (III.8a half-turn, label along the step)", good)
ok("framing: two reversals along phi1, phi2 compose to R(phi2 - phi1); out +x back -x = -1 (belt trick)",
   np.allclose(m(0.45)@m(0.15), np.array([[np.cos(0.6), -np.sin(0.6)], [np.sin(0.6), np.cos(0.6)]])) and np.allclose(m(np.pi/2)@m(0), -np.eye(2)))

# ---- the screw: the massless complex structure J (commutant of the three A's) pairs (1,4),(2,3) with OPPOSITE orientation
basisA = []
for i in range(4):
    for j in range(i + 1, 4):
        E = np.zeros((4, 4)); E[i, j] = 1; E[j, i] = -1; basisA.append(E)
Mc = np.array([np.concatenate([(E @ A - A @ E).ravel() for A in (Ax, Ay, Az)]) for E in basisA]).T
uu, ss, vv = np.linalg.svd(Mc); nul = vv[np.sum(ss > 1e-9):]
ok("screw: the antisymmetric commutant of A_x, A_y, A_z is one-dimensional (one complex structure J)", nul.shape[0] == 1)
Jm = sum(cf * E for cf, E in zip(nul[0], basisA)); Jm /= np.sqrt(-np.trace(Jm @ Jm) / 4)
if Jm[0, 3] < 0: Jm = -Jm
ok("screw: J pairs 1<->4 and 2<->3 (the framing pairs), with opposite orientation on the two copies; J^2 = -1",
   np.allclose(Jm, np.array([[0, 0, 0, 1], [0, 0, -1, 0], [0, 1, 0, 0], [-1, 0, 0, 0.]])) and np.allclose(Jm @ Jm, -I4))
ok("screw: T anticommutes with J (the flip reverses the screw); T on the face (Q1,Q4)->(Q2,Q3) is the quarter-turn [[0,1],[-1,0]]",
   np.allclose(Jm @ T + T @ Jm, 0) and np.allclose(np.array([[T[1, 0], T[1, 3]], [T[2, 0], T[2, 3]]]), np.array([[0, 1], [-1, 0.]])))
ok("sense = z-direction: coefficient of d_t in (A) is -A_z (states 1,4 forward and +z; 2,3 backward and -z)", np.allclose(adv, -Az))

# ---- the clock form as IV.7a (A) writes it, and the growing mode of "What the P's are"  (added 10 Nov 2026)
ok("IV.7a (A): in the clock the matrices are -A_z on d_t, A_z T on a, A_x A_z on c d_x, A_y A_z on c d_y",
   np.allclose(adv, -Az) and np.allclose(exch, Az @ T) and np.allclose(xt, Ax @ Az) and np.allclose(yt, Ay @ Az))
kx, ky, w, av = 0.7, -0.4, 0.3, 1.3
Mclock = adv*(-1j*w) + av*exch + 1j*kx*xt + 1j*ky*yt                      # the bracket of (F) on e^{i(kx x + ky y - w t)}, c = 1
ev = np.linalg.eigvals(av*exch + 1j*kx*xt + 1j*ky*yt)                     # w = 0
ok(f"the square of the clock operator on a mode is (a^2 + k_perp^2 - w^2) 1; at w = 0 its eigenvalues are +-sqrt(a^2 + k_perp^2), "
   f"so P = e^(-az/c) Q has the growing mode -a + sqrt(a^2 + k_perp^2) = {-av + np.sqrt(av**2 + kx**2 + ky**2):.4f} > 0",
   np.allclose(Mclock @ Mclock, (av**2 + kx**2 + ky**2 - w**2)*I4) and np.allclose(sorted(ev.real), sorted([-np.sqrt(av**2 + kx**2 + ky**2)]*2 + [np.sqrt(av**2 + kx**2 + ky**2)]*2))
   and np.allclose(ev.imag, 0))
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
