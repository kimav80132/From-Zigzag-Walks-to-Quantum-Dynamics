"""The spin current decomposed by direction, from the walk, no model.
Exact: j_k = c psi^dag A_k psi is the net flow of the k-th coin (IV.4d (B): along z it is Q+^2 - Q-^2, along x it is 2 Q+ Q-, the cross term of the two z-streams).
2+1 slow limit (III.9's basis, H = c sigma_z p_z + c sigma_x p_x + i hbar a T, T = [[0,1],[-1,0]], spin axis = the normal y):
  large component along the +1 eigenvector of the mass matrix, small component eta = k phi/2mc; then
  j_z = rho v_z - (hbar/2m) d_x rho,   j_x = rho v_x + (hbar/2m) d_z rho   (= rho (v + u x s) with s = -y-hat, the axis of the positive-energy column),
so the sideways flow of a particle at rest is the osmotic velocity of the perpendicular walk, turned: j_x = +rho u_z, j_z = -rho u_x.
(Docstring corrected 10 Nov 2026: it carried the signs of an earlier draft, s = +y-hat; the checks below already had the text's signs.  Check 0'' added.)
In the walk's words: j_x = 2 Q+ Q- is the interference of the forward and backward z-streams, and III.11 says that interference (the relative phase) is u_z."""
import sympy as sp, numpy as np
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
i = sp.I
x, z = sp.symbols('x z', real=True); hbar, m, c = sp.Rational(7, 5), sp.Rational(3, 2), sp.Rational(11, 4)
sx = sp.Matrix([[0, 1], [1, 0]]); sy = sp.Matrix([[0, -i], [i, 0]]); sz = sp.diag(1, -1); T = sp.Matrix([[0, 1], [-1, 0]])
Mmass = i*T                                     # the mass matrix of III.9, = -sigma_y
ev = Mmass.eigenvects()
e_plus = [v[0] for val, _, v in ev if val == 1][0]; e_minus = [v[0] for val, _, v in ev if val == -1][0]
e_plus = e_plus/sp.sqrt((e_plus.H*e_plus)[0]); e_minus = e_minus/sp.sqrt((e_minus.H*e_minus)[0])
ok("0. mass matrix i T = -sigma_y, eigenvectors (1,-i)/sqrt2 (+1: positive energy at rest) and (1,i)/sqrt2 (-1)",
   Mmass == -sy and sp.simplify(Mmass*e_plus - e_plus) == sp.zeros(2, 1) and sp.simplify(Mmass*e_minus + e_minus) == sp.zeros(2, 1)
   and sp.simplify(sp.Abs((e_plus.H*sp.Matrix([1, -i])/sp.sqrt(2))[0]) - 1) == 0)
ok("0'. the spin axis of the large component is the normal: e_plus^dag sigma e_plus = (0, -1, 0)",
   [sp.simplify((e_plus.H*M*e_plus)[0]) for M in (sx, sy, sz)] == [0, -1, 0])
# the matrix elements the text writes, with its own e_+ = (1,-i)/sqrt2, e_- = (1,i)/sqrt2
tp, tm = sp.Matrix([1, -i])/sp.sqrt(2), sp.Matrix([1, i])/sp.sqrt(2)
ok("0''. e_-^dag sigma_z e_+ = 1, e_-^dag sigma_x e_+ = -i, e_+^dag sigma_z e_- = 1, e_+^dag sigma_x e_- = i, and e_+^dag sigma_z e_+ = e_+^dag sigma_x e_+ = 0",
   [sp.simplify((a_.H*M*b_)[0]) for a_, M, b_ in ((tm, sz, tp), (tm, sx, tp), (tp, sz, tm), (tp, sx, tm), (tp, sz, tp), (tp, sx, tp))] == [1, -i, 1, i, 0, 0])
rho = sp.exp(x*z/3 - z**2/5 + x/2); S = 2*x*z/3 + x**2/4
phi0 = sp.sqrt(rho)*sp.exp(i*S/hbar)
p = lambda f: {'z': -i*hbar*sp.diff(f, z), 'x': -i*hbar*sp.diff(f, x)}
# small component: eta0 = e_minus^dag (sigma_z p_z + sigma_x p_x) e_plus phi0 / 2mc  (the stationary point of the fast pair, III.11)
kz = (e_minus.H*sz*e_plus)[0]; kx = (e_minus.H*sx*e_plus)[0]
eta0 = (kz*p(phi0)['z'] + kx*p(phi0)['x'])/(2*m*c)
psi = phi0*e_plus + eta0*e_minus
jz = c*(psi.H*sz*psi)[0]; jx = c*(psi.H*sx*psi)[0]
u = {'z': hbar/(2*m)*sp.diff(sp.log(rho), z), 'x': hbar/(2*m)*sp.diff(sp.log(rho), x)}
v = {'z': sp.diff(S, z)/m, 'x': sp.diff(S, x)/m}
pts = [(0.3, -0.4), (-0.2, 0.5), (0.6, 0.1)]
zero = lambda ex: max(abs(complex(sp.N(ex.subs({x: a, z: b}), 20))) for a, b in pts) < 1e-9
sgn = -1   # spin axis is -y for e_plus (0'), so s = -y-hat and u x s = (u_z, 0, -u_x): j_x = +rho u_z, j_z = -rho u_x
ok("1. slow limit, 2+1: j_z = rho v_z - rho u_x and j_x = rho v_x + rho u_z, i.e. j = rho (v + u x s) with s = -y-hat: each coin's net flow is Nelson's v along it plus the osmotic velocity of the PERPENDICULAR direction, turned",
   zero(jz - (rho*v['z'] + sgn*rho*u['x'])) and zero(jx - (rho*v['x'] - sgn*rho*u['z'])))
# the walk's words: j_x = c psi^dag sigma_x psi = 2c Re(psi_+^* psi_-): the cross term of the two z-streams; at rest (S = 0) it is -sgn rho u_z
phi0r = sp.sqrt(rho); eta0r = (kz*p(phi0r)['z'] + kx*p(phi0r)['x'])/(2*m*c); psir = phi0r*e_plus + eta0r*e_minus
jx_rest = c*(psir.H*sx*psir)[0]
ok("2. j_x = 2c Re(psi_+^* psi_-): the interference of the forward and backward z-streams (IV.4d (B)); at rest it is rho u_z (up to the sense of the axis): "
   "the relative phase of the two time senses, III.11's c phi = u, read sideways by A_x",
   zero(jx - c*(psi[0].conjugate()*psi[1] + psi[1].conjugate()*psi[0])) and zero(jx_rest + sgn*rho*u['z']))
# the text's intermediate line: 2mc phi^* eta = (rho d_z S - (hbar/2) d_x rho) - i (rho d_x S + (hbar/2) d_z rho), with the text's e_+-
eta_t = (p(phi0)['z'] - i*p(phi0)['x'])/(2*m*c)
ok("3. 2mc phi^* eta = (rho d_z S - (hbar/2) d_x rho) - i (rho d_x S + (hbar/2) d_z rho), eta = (p_z - i p_x) phi/2mc",
   zero(2*m*c*sp.conjugate(phi0)*eta_t - ((rho*sp.diff(S, z) - hbar/2*sp.diff(rho, x)) - i*(rho*sp.diff(S, x) + hbar/2*sp.diff(rho, z)))))
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
