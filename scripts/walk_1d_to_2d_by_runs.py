"""Non-circular chain: 1D walk only  ->  runs (boosts)  ->  the 2D lines  ->  the 2D walk with its U-turn signs.
Input is ONLY the 1+1 walk (clock z, counters P+ forward in t, P- backward in t, flips a).
W1  1D walk step -> 1+1 lines
W2  run along z: Q+ -> e^{eta/2} Q+, Q- -> e^{-eta/2} Q-   (same a + mirror; checked in spinor_chain_2plus1.py S4)
W3  directions: loop of runs = turn;  [R, B_z] = -B_x, [R, B_x] = B_z with R = B_z B_x - B_x B_z
W4  counters, x-rule unknown A_x = [[p,q],[r,s]]: the same relations force p = s = 0, q r = 1
W5  a classical x-rule (each counter carries its own x-direction: A_x diagonal) fails
W6  same flip rate along x: T anticommutes with A_x and A_z
W7  z-clock form: the sideways matrix is A_z A_x = [[0,1],[-1,0]]
W8  the 2D walk: time move + sideways move with coins (1 +- A_z A_x)/2 + flips; step -> the 2+1 P lines
"""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
t, x, z = sp.symbols('t x z', real=True); c, a, D = sp.symbols('c a Delta', positive=True)
Pp, Pm = [sp.Function(n) for n in ('Pp', 'Pm')]
first = lambda e: sp.expand(sp.series(e, D, 0, 2).removeO())
# W1
stepP = (1 - a*D/c)*Pp(t - D/c, z) + a*D/c*Pm(t - D/c, z)
stepM = (1 - a*D/c)*Pm(t + D/c, z) + a*D/c*Pp(t + D/c, z)
lp = sp.simplify((first(stepP) - Pp(t, z))/D*c); lm = sp.simplify((first(stepM) - Pm(t, z))/D*c)
ok("W1  1D walk step -> c d_z P+ = -d_t P+ - a P+ + a P-,  c d_z P- = d_t P- - a P- + a P+",
   sp.simplify(lp - (-Pp(t, z).diff(t) - a*Pp(t, z) + a*Pm(t, z))) == 0 and sp.simplify(lm - (Pm(t, z).diff(t) - a*Pm(t, z) + a*Pp(t, z))) == 0)
# W3
Bz = sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]); Bx = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
R = Bz*Bx - Bx*Bz
ok(f"W3  R = {R.tolist()} turns x into z;  R B_z - B_z R = -B_x,  R B_x - B_x R = +B_z", R*Bz - Bz*R == -Bx and R*Bx - Bx*R == Bz)
# W4
p, q, r, s = sp.symbols('p q r s', real=True)
Az = sp.diag(1, -1); Ax = sp.Matrix([[p, q], [r, s]]); Kz, Kx = Az/2, Ax/2; J = Kz*Kx - Kx*Kz
sol = sp.solve(list(J*Kz - Kz*J + Kx) + list(J*Kx - Kx*J - Kz), [p, q, r, s], dict=True)
ok(f"W4a counters must satisfy the same relations: solutions {sol}  (p = s = 0, q r = 1)", sol == [{p: 0, q: 1/r, s: 0}])
ok("W4b with q = r (same unit on both counters): A_x = +-[[0,1],[1,0]];  A_x^2 = 1 and A_x A_z = -A_z A_x",
   all((M*M == sp.eye(2) and M*Az == -Az*M) for M in (sp.Matrix([[0, 1], [1, 0]]), -sp.Matrix([[0, 1], [1, 0]]))))
# W5
d1, d2 = sp.symbols('d1 d2', real=True)
Kd = sp.diag(d1, d2)/2; Jd = Kz*Kd - Kd*Kz
ok("W5  classical x-rule (diagonal): the loop of runs does nothing to the counters (J = 0), so J K_z - K_z J = 0 != -K_x",
   Jd == sp.zeros(2) and sp.solve(list(Jd*Kz - Kz*Jd + Kd), [d1, d2], dict=True) == [{d1: 0, d2: 0}])
# W6
Ax = sp.Matrix([[0, 1], [1, 0]]); T = sp.Matrix([[0, 1], [-1, 0]])
ok("W6  T = [[0,1],[-1,0]] (the 1D flips) anticommutes with A_z and with A_x", T*Az == -Az*T and T*Ax == -Ax*T)
# W7
X = Az*Ax
ok(f"W7  z-clock sideways matrix A_z A_x = {X.tolist()};  flips A_z T = {(Az*T).tolist()}", X == sp.Matrix([[0, 1], [-1, 0]]) and Az*T == sp.Matrix([[0, 1], [1, 0]]))
# W8: the 2D walk; P(t, x, z); step Delta: time move + flips (1D walk), then sideways move ell = Delta with coins
P = [sp.Function('Pp'), sp.Function('Pm')]
def after_time(i, tt, xx):                         # 1D walk part at (tt, xx), arriving on counter i
    if i == 0: return (1 - a*D/c)*P[0](tt - D/c, xx, z) + a*D/c*P[1](tt - D/c, xx, z)
    return (1 - a*D/c)*P[1](tt + D/c, xx, z) + a*D/c*P[0](tt + D/c, xx, z)
Cl = (sp.eye(2) + X)/2; Cr = (sp.eye(2) - X)/2      # from x - ell and from x + ell
ell = D
new = [sum(Cl[i, j]*after_time(j, t, x - ell) + Cr[i, j]*after_time(j, t, x + ell) for j in range(2)) for i in range(2)]
lines = [-P[0](t, x, z).diff(t) - a*P[0](t, x, z) + a*P[1](t, x, z) - c*P[1](t, x, z).diff(x),
         P[1](t, x, z).diff(t) - a*P[1](t, x, z) + a*P[0](t, x, z) + c*P[0](t, x, z).diff(x)]
res = [sp.simplify((first(new[i]) - P[i](t, x, z))/D*c - lines[i]) for i in range(2)]
ok("W8a 2D walk (time move, flips, sideways coins (1 +- A_z A_x)/2) -> exactly the 2+1 P lines", res == [0, 0])
ok("W8b the coins: from x - ell [[1/2, 1/2], [-1/2, 1/2]], from x + ell [[1/2, -1/2], [1/2, 1/2]]: the U-turn weights carry the sign",
   Cl == sp.Matrix([[sp.Rational(1, 2), sp.Rational(1, 2)], [-sp.Rational(1, 2), sp.Rational(1, 2)]]) and
   Cr == sp.Matrix([[sp.Rational(1, 2), -sp.Rational(1, 2)], [sp.Rational(1, 2), sp.Rational(1, 2)]]))
