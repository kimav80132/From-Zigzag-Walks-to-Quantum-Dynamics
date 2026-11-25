"""What the matrices become when Dirac is written in Madelung variables.
One copy: psi = sqrt(rho) e^{iS/hbar} chi, chi = (cos th/2, e^{i phi} sin th/2), axis s = chi^dag sigma chi (IV.6).
1. Every 2x2 matrix M = m0 + m.sigma sandwiched becomes a real function: psi^dag M psi = rho (m0 + m.s). The half-turns become the axis components.
2. The rotation generators e^{i theta sigma_k/2} acting on chi rotate the axis s by the full angle theta about k (half angle on the column, full on the axis).
3. In canonical variables (IV.6) the spin density s_k = (hbar/2) rho s_k generates that same rotation through the bracket: {phi, s_z} = 1 shifts phi.
4. The mass term between two copies: Re[psi_L^dag psi_R] = sqrt(rho_L rho_R) Re[e^{i(S_R-S_L)/hbar} chi_L^dag chi_R], |chi_L^dag chi_R|^2 = (1 + s_L.s_R)/2 (IV.7).
5. The flatlander's real vector v(alpha) is the axis lying in the plane at angle -2 alpha; B3 shifts that angle; the tilt iB2 rotates the axis about x (so the i is gone in s)."""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
i = sp.I
rho, S, th, ph, hbar, theta, m0 = sp.symbols('rho S theta phi hbar vartheta m0', real=True)
rho = sp.symbols('rho', positive=True); hbar = sp.symbols('hbar', positive=True)
sx = sp.Matrix([[0, 1], [1, 0]]); sy = sp.Matrix([[0, -i], [i, 0]]); sz = sp.diag(1, -1); sig = (sx, sy, sz)
chi = sp.Matrix([sp.cos(th/2), sp.exp(i*ph)*sp.sin(th/2)])
psi = sp.sqrt(rho)*sp.exp(i*S/hbar)*chi
simp = lambda e: sp.simplify(sp.expand(e.rewrite(sp.cos)))
s = sp.Matrix([simp((chi.H*M*chi)[0]) for M in sig])
ok("0. axis: chi^dag sigma chi = (sin th cos phi, sin th sin phi, cos th), real, unit length",
   simp(s[0] - sp.sin(th)*sp.cos(ph)) == 0 and simp(s[1] - sp.sin(th)*sp.sin(ph)) == 0 and simp(s[2] - sp.cos(th)) == 0 and simp(s.dot(s) - 1) == 0)
m = sp.Matrix(sp.symbols('m1 m2 m3', real=True))
M = m0*sp.eye(2) + m[0]*sx + m[1]*sy + m[2]*sz
ok("1. any matrix sandwiched: psi^dag (m0 + m.sigma) psi = rho (m0 + m.s). sigma_y, imaginary, gives the real number rho sin th sin phi",
   simp((psi.H*M*psi)[0] - rho*(m0 + m.dot(s))) == 0)
# 2. generators rotate the axis by the full angle
def axis_of(c): return sp.Matrix([simp((c.H*Mk*c)[0]) for Mk in sig])
Rz = sp.Matrix([[sp.cos(theta), -sp.sin(theta), 0], [sp.sin(theta), sp.cos(theta), 0], [0, 0, 1]])
Rx = sp.Matrix([[1, 0, 0], [0, sp.cos(theta), -sp.sin(theta)], [0, sp.sin(theta), sp.cos(theta)]])
Gz = sp.cos(theta/2)*sp.eye(2) - i*sp.sin(theta/2)*sz        # e^{-i theta sigma_z/2}
Gx = sp.cos(theta/2)*sp.eye(2) - i*sp.sin(theta/2)*sx
ok("2. e^{-i theta sigma_z/2} chi has axis R_z(theta) s, and e^{-i theta sigma_x/2} chi has axis R_x(theta) s: half angle on the column, full angle on the axis",
   (axis_of(Gz*chi) - Rz*s).applyfunc(simp) == sp.zeros(3, 1) and (axis_of(Gx*chi) - Rx*s).applyfunc(simp) == sp.zeros(3, 1))
# 3. canonical variables of IV.6: rho1 = rho cos^2, rho2 = rho sin^2, S1 = S, S2 = S + hbar phi; P_phi = (hbar/2)(rho1 - rho2) = (hbar/2) rho cos th
r1, r2, S1, S2 = sp.symbols('rho1 rho2 S1 S2', real=True)
pb = lambda f, g: sum(sp.diff(f, q)*sp.diff(g, p) - sp.diff(f, p)*sp.diff(g, q) for q, p in ((r1, S1), (r2, S2)))   # {rho_a, S_a} = 1
phi_c = (S2 - S1)/hbar; Pphi = hbar/2*(r1 - r2)
rho_c = r1 + r2; costh = (r1 - r2)/(r1 + r2); sinth = 2*sp.sqrt(r1*r2)/(r1 + r2)
s_z = Pphi; s_x = hbar/2*rho_c*sinth*sp.cos(phi_c); s_y = hbar/2*rho_c*sinth*sp.sin(phi_c)
ok("3. {phi, s_z} = 1 and {s_x, s_y} = s_z: the spin density (hbar/2) rho s is the generator; the flow of s_z shifts phi, i.e. rotates the axis about z (IV.6)",
   sp.simplify(pb(phi_c, s_z) - 1) == 0 and sp.simplify(pb(s_x, s_y) - s_z) == 0)
# 4. the mass term between two copies
thL, phL, thR, phR, rL, rR, SL, SR = sp.symbols('theta_L phi_L theta_R phi_R rho_L rho_R S_L S_R', real=True)
rL, rR = sp.symbols('rho_L rho_R', positive=True)
chiL = sp.Matrix([sp.cos(thL/2), sp.exp(i*phL)*sp.sin(thL/2)]); chiR = sp.Matrix([sp.cos(thR/2), sp.exp(i*phR)*sp.sin(thR/2)])
sL = axis_of(chiL); sR = axis_of(chiR)
ov = (chiL.H*chiR)[0]
ok("4. |chi_L^dag chi_R|^2 = (1 + s_L.s_R)/2: the mass (the turn between the copies) becomes the overlap of the two axes times the cosine of the relative phase (IV.7)",
   simp(sp.expand(ov*sp.conjugate(ov)) - (1 + sL.dot(sR))/2) == 0)
# 5. the flatlander's real vector
al, a = sp.symbols('alpha a', real=True)
B2 = sz; B3 = sp.Matrix([[0, -1], [1, 0]])
Uinv = sp.Matrix([[1, i], [1, -i]])/sp.sqrt(2)
v = sp.Matrix([sp.cos(al), sp.sin(al)])
ax = lambda col: axis_of(Uinv*col)
ok("5. the real vector v(alpha) is the axis (cos 2alpha, -sin 2alpha, 0): in the plane, at twice the angle; e^{aB3} v shifts it by 2a in the plane; "
   "e^{a iB2} v rotates it about x by 2a out of the plane — in s the i is gone and the angle is doubled",
   (ax(v) - sp.Matrix([sp.cos(2*al), -sp.sin(2*al), 0])).applyfunc(simp) == sp.zeros(3, 1)
   and (ax((sp.cos(a)*sp.eye(2) + sp.sin(a)*B3)*v) - sp.Matrix([sp.cos(2*al + 2*a), -sp.sin(2*al + 2*a), 0])).applyfunc(simp) == sp.zeros(3, 1)
   and (ax((sp.cos(a)*sp.eye(2) + i*sp.sin(a)*B2)*v) - Rx.subs(theta, -2*a)*sp.Matrix([sp.cos(2*al), -sp.sin(2*al), 0])).applyfunc(simp) == sp.zeros(3, 1))
