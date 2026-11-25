"""Does the spin flow need a hole? Cylindrical radius r, axis z, D = hbar/2m. Azimuthal speed of the osmotic circulation u x s is w = -D d(ln rho)/dr (vortex_worked 6).
1. Gaussian core (no phase winding): rho ~ exp(-r^2/a^2): w = hbar r/(m a^2) -> 0 at the centre; the cloud is fullest there; no hole needed.
2. A hole: an orbital vortex S = hbar*azimuth forces rho -> 0 at the core, rho ~ r^2 exp(-r^2/a^2). Then w = -hbar/(m r) + hbar r/(m a^2): near the core a 1/r vortex of circulation -h/m, opposite to the outer rotation.
3. With the orbital flow v = +hbar/(m r) (same sense as the spin), the total current rho(v + w) is regular at the core: the osmotic circulation cancels the phase vortex's singularity.
4. angular momenta per particle: orbital hbar (phase), Gordon spin flow hbar (any profile), energy-flow spin hbar/2; total energy-flow L + S = 3 hbar/2, charge flow hbar(l + 2s) = 2 hbar: the Lande structure."""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
r, a, hbar, m = sp.symbols('r a hbar m', positive=True); D = hbar/(2*m)
w = lambda rho: -D*sp.diff(sp.log(rho), r)
rhoG = sp.exp(-r**2/a**2)/(sp.pi*a**2); rhoH = r**2*sp.exp(-r**2/a**2)/(sp.pi*a**4)
ok("0. both profiles normalized in the plane", sp.simplify(sp.integrate(2*sp.pi*r*rhoG, (r, 0, sp.oo)) - 1) == 0 and sp.simplify(sp.integrate(2*sp.pi*r*rhoH, (r, 0, sp.oo)) - 1) == 0)
ok("1. no hole: w = hbar r/(m a^2), zero at the centre, rigid; rho maximal at r = 0", sp.simplify(w(rhoG) - hbar*r/(m*a**2)) == 0 and sp.limit(rhoG, r, 0) == 1/(sp.pi*a**2))
ok("2. hole (rho ~ r^2): w = -hbar/(m r) + hbar r/(m a^2): a 1/r vortex at the core of circulation -h/m, opposite in sense to the rigid rotation outside",
   sp.simplify(w(rhoH) - (-hbar/(m*r) + hbar*r/(m*a**2))) == 0 and sp.simplify(sp.limit(2*sp.pi*r*w(rhoH), r, 0) + 2*sp.pi*hbar/m) == 0)
v_orb = hbar/(m*r)
ok("3. orbital vortex v = hbar/(m r) in the same sense: total speed v + w = hbar r/(m a^2), the 1/r parts cancel; the current rho(v + w) is regular at the core",
   sp.simplify(v_orb + w(rhoH) - hbar*r/(m*a**2)) == 0)
L = lambda rho, speed: sp.integrate(2*sp.pi*r*rho*m*r*speed, (r, 0, sp.oo))
ok("4. per particle, hole profile: orbital hbar; Gordon spin flow hbar; energy-flow spin hbar/2. Energy flow total 3hbar/2 = l + s; charge flow 2 hbar = l + 2s",
   sp.simplify(L(rhoH, v_orb) - hbar) == 0 and sp.simplify(L(rhoH, w(rhoH)) - hbar) == 0 and sp.simplify(L(rhoH, w(rhoH)/2) - hbar/2) == 0)
ok("4'. opposite senses (spin down with the orbital vortex up): speeds add to hbar/(m r)(1 + 1) - ... : the 1/r parts double at the core, charge flow 0, energy flow hbar/2",
   sp.simplify(v_orb - w(rhoH) - (2*hbar/(m*r) - hbar*r/(m*a**2))) == 0 and sp.simplify(L(rhoH, v_orb) - L(rhoH, w(rhoH))) == 0 and sp.simplify(L(rhoH, v_orb) - L(rhoH, w(rhoH)/2) - hbar/2) == 0)
