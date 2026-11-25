"""IV.8 (F): the full vortex, worked. A Gaussian cloud rho = (pi a^2)^{-3/2} exp(-r^2/a^2), axis s = z-hat uniform, no orbital phase (S = const).
Then v = grad S/m = 0 and the whole flow is the osmotic circulation u x s (IV8_axis_fluid.py C5).
1. the flow field: u x s = (hbar/m a^2) (-y, x, 0): a rigid rotation about the axis with angular velocity omega = hbar/(m a^2); regular at the centre, no 1/r line.
2. circulation of the Gordon flow on a circle of radius r in the plane: Gamma(r) = 2 pi hbar r^2/(m a^2); = h/m at r = a. The energy-flow (Belinfante) circulation is half: h/2m at r = a.
3. angular momentum about the axis per particle: Gordon hbar; Belinfante hbar/2; the orbit of true_spin 3b hbar/2.
4. the old 1/r vortex v = (hbar/m r) theta-hat is the orbital phase S = hbar*azimuth: angular momentum hbar, singular at the centre, independent of rho. It is not the spin.
5. magnetic moment of the Gordon flow: e hbar/2mc on a spin of hbar/2: g = 2.
6. the old collision integral's missing O(lambda^2) term and Thomas 1/2: the exact flow is rigid (omega uniform), so no model constant survives; the half is Belinfante's, not Thomas's."""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
x, y, z, r, phi_ = sp.symbols('x y z r phi', real=True)
a, hbar, m, c, e = sp.symbols('a hbar m c e', positive=True)
D = hbar/(2*m)
rho = (sp.pi*a**2)**sp.Rational(-3, 2)*sp.exp(-(x**2 + y**2 + z**2)/a**2)
X = (x, y, z)
grad = lambda f: sp.Matrix([sp.diff(f, v) for v in X])
curl = lambda F: sp.Matrix([sp.diff(F[2], y) - sp.diff(F[1], z), sp.diff(F[0], z) - sp.diff(F[2], x), sp.diff(F[1], x) - sp.diff(F[0], y)])
s = sp.Matrix([0, 0, 1]); R = sp.Matrix([x, y, z])
u = D*grad(sp.log(rho))
w = u.cross(s)
omega = hbar/(m*a**2)
ok("1. u x s = (hbar/m a^2)(-y, x, 0) = omega z-hat x r: a rigid rotation of the cloud about its axis, omega = hbar/(m a^2), regular at the centre",
   sp.simplify(w - omega*sp.Matrix([-y, x, 0])) == sp.zeros(3, 1) and sp.simplify(w - s.cross(R)*omega) == sp.zeros(3, 1))
ok("1'. it is divergence-free with rho: div(rho u x s) = 0, and (hbar/2m) curl(rho s) = rho (u x s) (IV8 C4 for a uniform axis)",
   sp.simplify(sum(sp.diff((rho*w)[k], X[k]) for k in range(3))) == 0 and sp.simplify(hbar/(2*m)*curl(rho*s) - rho*w) == sp.zeros(3, 1))
Gamma = 2*sp.pi*r*(omega*r)
ok("2. circulation on a circle of radius r: Gamma(r) = 2 pi hbar r^2/(m a^2); at r = a it is h/m, and the energy-flow half is h/2m there: the half quantum at the cloud's radius",
   sp.simplify(Gamma - 2*sp.pi*hbar*r**2/(m*a**2)) == 0 and sp.simplify(Gamma.subs(r, a) - 2*sp.pi*hbar/m) == 0 and sp.simplify(Gamma.subs(r, a)/2 - sp.pi*hbar/m) == 0)
Lz = lambda G: sp.integrate(R.cross(G)[2], (x, -sp.oo, sp.oo), (y, -sp.oo, sp.oo), (z, -sp.oo, sp.oo))
ok("3. angular momentum per particle about the axis: Gordon flow int r x (m rho u x s) = hbar; energy flow (half the curl) hbar/2; orbit m (hbar/2mc) c = hbar/2",
   sp.simplify(Lz(m*rho*w) - hbar) == 0 and sp.simplify(Lz(m*rho*w/2) - hbar/2) == 0 and sp.simplify(m*(hbar/(2*m*c))*c - hbar/2) == 0)
# 4. the old 1/r vortex: orbital phase S = hbar * azimuth
S_orb = hbar*sp.atan2(y, x)
v_orb = grad(S_orb)/m
ok("4. the old picture's v = (hbar/m r) theta-hat is grad(hbar * azimuth)/m: a phase winding, angular momentum hbar per particle, singular at r = 0, independent of rho. Orbital, not spin",
   sp.simplify(v_orb - hbar/(m*(x**2 + y**2))*sp.Matrix([-y, x, 0])) == sp.zeros(3, 1) and sp.simplify(Lz(m*rho*v_orb) - hbar) == 0)
ok("5. magnetic moment of the Gordon flow (e/2c) int r x (rho u x s) = e hbar/2mc: on the spin hbar/2 that is g = 2",
   sp.simplify(e/(2*c)*Lz(rho*w) - e*hbar/(2*m*c)) == 0)
# 6. the Gaussian is not special: any axisymmetric rho(r_perp) gives u x s = -D (rho'/rho) r_perp theta-hat = (D/r_perp)(-d ln rho/d ln r_perp) theta-hat; circulation per particle integrates to hbar
rp = sp.symbols('r_p', positive=True)
# 6. any axisymmetric profile: azimuthal speed -D d(ln rho)/dr_perp, angular momentum per particle = 2 m D = hbar by one integration by parts; two explicit profiles:
for nm, prof in (("Gaussian", sp.exp(-rp**2/a**2)/(sp.pi*a**2)), ("exponential", sp.exp(-rp/a)/(2*sp.pi*a**2))):
    L = sp.integrate(2*sp.pi*rp*prof*m*rp*(-D*sp.diff(sp.log(prof), rp)), (rp, 0, sp.oo))
    ok(f"6{nm[0]}. {nm} profile: angular momentum per particle of the osmotic circulation = hbar (Gordon), independent of the width a", sp.simplify(L - hbar) == 0)
