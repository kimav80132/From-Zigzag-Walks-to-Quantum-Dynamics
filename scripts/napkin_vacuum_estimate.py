"""The back-of-envelope numbers of "Next questions", item 6 (master v46). Nothing here is a result of the walk: it is arithmetic on the
cosmological parameters of Planck 2018 (arXiv:1807.06209: H0 = 67.4 km/s/Mpc, Omega_m = 0.315, Omega_b h^2 = 0.0224; flat, so Omega_Lambda = 0.685)."""
import math
G, c, Mpc, m_p = 6.67430e-11, 2.99792458e8, 3.0856775814913673e22, 1.67262192369e-27
H0 = 67.4e3/Mpc; h = 0.674; Om = 0.315; OL = 1 - Om; Ob = 0.0224/h**2
rho_c = 3*H0**2/(8*math.pi*G)                     # kg/m^3
u_c = rho_c*c**2                                   # J/m^3
u_L, u_m, u_b = OL*u_c, Om*u_c, Ob*u_c
n_b = Ob*rho_c/m_p
ok = True
def check(msg, cond):
    global ok; ok = ok and cond; print(("PASS  " if cond else "FAIL  ") + msg)
check("critical density %.2e kg/m^3 = %.2e J/m^3" % (rho_c, u_c), abs(rho_c/8.53e-27 - 1) < 0.01)
check("observed dark energy %.2e J/m^3 (text: 5.3e-10)" % u_L, abs(u_L/5.3e-10 - 1) < 0.02)
check("rest energy of all matter %.2e J/m^3 (text: 2.4e-10)" % u_m, abs(u_m/2.4e-10 - 1) < 0.03)
check("rest energy of baryons %.2e J/m^3 (text: 3.8e-11), %.2f baryons per m^3 (text: 0.25)" % (u_b, n_b), abs(u_b/3.8e-11 - 1) < 0.03 and abs(n_b/0.25 - 1) < 0.03)
check("dark energy / matter = %.1f, dark energy / baryons = %.1f (text: a factor of 2 to 14)" % (u_L/u_m, u_L/u_b), abs(u_L/u_m - 2.2) < 0.1 and abs(u_L/u_b - 13.9) < 0.3)
r1, r3 = Om*2**3/OL, Om*4**3/OL
check("matter / dark energy at redshift 1: %.1f (text: 3.7); at redshift 3: %.1f (text: 29)" % (r1, r3), abs(r1 - 3.7) < 0.05 and abs(r3 - 29.4) < 0.5)
hbar, GeV = 1.054571817e-34, 1.602176634e-10
u_ew = (100*GeV)**4/(hbar*c)**3; u_pl = c**7/(hbar*G**2)
check("a sum cut at 100 GeV: %.1e J/m^3, 10^%.1f times the observed density; cut at the Planck length: %.1e J/m^3, 10^%.1f times (text: about 10^55 to "
      "about 10^123)" % (u_ew, math.log10(u_ew/u_L), u_pl, math.log10(u_pl/u_L)), 54 < math.log10(u_ew/u_L) < 56 and 122.5 < math.log10(u_pl/u_L) < 123.5)
print("ALL PASS" if ok else "SOME FAILED"); raise SystemExit(0 if ok else 1)
