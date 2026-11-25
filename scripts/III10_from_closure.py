"""Does III.10 (the Madelung walk of Dirac, the walk that carries the particle) follow from III.1's loop closure?
Two-copy real walk (no field): copies do not mix, so loops close copy by copy and the closure weight through (x, eps) is
   rho_eps = Q_eps^2 + Q'_eps^2 = |psi_eps|^2   (loop_closure.py per copy, added).
Claim: the closure weights obey III.10's population equations, with III.10's exchange term, as a consequence of III.1's equation:
   (d_t + c d_x) rho_+ = 2a (Q+Q- + Q'+Q'-),  (d_t - c d_x) rho_- = -2a (Q+Q- + Q'+Q'-),
and 2a (Q+Q- + Q'+Q'-) = 2a Re(psi+* psi-) [natural basis] = 2a sqrt(rho+ rho-) sin(phi) [standard basis of III.10].
So III.10 is the forward-in-t bookkeeping of where III.1's loops close. Not assumed: the square."""
import sympy as sp
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
t, x, a, c = sp.symbols('t x a c', real=True)
Qp, Qm, Pp, Pm = [sp.Function(n, real=True)(t, x) for n in ('Qp', 'Qm', 'Pp', 'Pm')]   # copy 1: Q, copy 2: P (primed)
# III.1 read in t, each copy: d_t Q+ = -c d_x Q+ + a Q-,  d_t Q- = +c d_x Q- - a Q+
rules = {Qp.diff(t): -c*Qp.diff(x) + a*Qm, Qm.diff(t): c*Qm.diff(x) - a*Qp, Pp.diff(t): -c*Pp.diff(x) + a*Pm, Pm.diff(t): c*Pm.diff(x) - a*Pp}
rhop, rhom = Qp**2 + Pp**2, Qm**2 + Pm**2
X = Qp*Qm + Pp*Pm
lhs_p = (rhop.diff(t) + c*rhop.diff(x)).subs(rules); lhs_m = (rhom.diff(t) - c*rhom.diff(x)).subs(rules)
ok("1. (d_t + c d_x)(Q+^2 + Q'+^2) = 2a(Q+Q- + Q'+Q'-) and (d_t - c d_x)(Q-^2 + Q'-^2) = -2a(Q+Q- + Q'+Q'-): III.10's population equations from III.1",
   sp.simplify(lhs_p - 2*a*X) == 0 and sp.simplify(lhs_m + 2*a*X) == 0)
# the exchange term in III.10's variables: standard basis psi- -> i psi-
psip, psim = Qp + sp.I*Pp, Qm + sp.I*Pm
ok("2. Q+Q- + Q'+Q'- = Re(psi+* psi-) (natural basis)", sp.simplify(sp.re(sp.expand(sp.conjugate(psip)*psim)) - X) == 0)
r1, r2, s1, s2 = sp.symbols('r1 r2 s1 s2', positive=True)
psip_s, psim_s = sp.sqrt(r1)*sp.exp(sp.I*s1), sp.sqrt(r2)*sp.exp(sp.I*s2)             # III.10's Madelung components, standard basis
psim_nat = psim_s/sp.I                                                                  # psi-_std = i psi-_nat (III.9)
val = sp.re(sp.expand((sp.conjugate(psip_s)*psim_nat).rewrite(sp.cos)))
ok("3. in the standard basis Re(psi+* psi-_nat) = sqrt(rho+ rho-) sin(S- - S+) = sqrt(rho+ rho-) sin(phi): III.10's exchange term exactly",
   sp.simplify(val - sp.sqrt(r1*r2)*sp.sin(s2 - s1)) == 0)
# 4. the phases (III.15a (G), last sentence): theta_eps = atan2(Q'_eps, Q_eps), W = Q+ Q'- - Q'+ Q-
thp, thm = sp.atan2(Pp, Qp), sp.atan2(Pm, Qm); W = Qp*Pm - Pp*Qm
ph_p = (rhop*(thp.diff(t) + c*thp.diff(x))).subs(rules); ph_m = (rhom*(thm.diff(t) - c*thm.diff(x))).subs(rules)
ok("4. rho_+ (d_t + c d_x) theta_+ = a W and rho_- (d_t - c d_x) theta_- = a W, with W = Q+ Q'- - Q'+ Q-: the phases close the same way",
   sp.simplify(ph_p - a*W) == 0 and sp.simplify(ph_m - a*W) == 0)
# 5. the same in III.10's variables: S_+ = hbar theta_+, S_- = hbar (theta_- + pi/2), phi = (S_- - S_+)/hbar: W = -sqrt(rho+ rho-) cos(phi) (III.8c (9))
thA, thB = sp.symbols('thetaA thetaB', real=True)
Wpolar = sp.sqrt(r1)*sp.cos(thA)*sp.sqrt(r2)*sp.sin(thB) - sp.sqrt(r1)*sp.sin(thA)*sp.sqrt(r2)*sp.cos(thB)
Kpolar = sp.sqrt(r1)*sp.cos(thA)*sp.sqrt(r2)*sp.cos(thB) + sp.sqrt(r1)*sp.sin(thA)*sp.sqrt(r2)*sp.sin(thB)
phi = thB + sp.pi/2 - thA
ok("5. with phi = theta_- + pi/2 - theta_+: Q+Q- + Q'+Q'- = sqrt(rho+ rho-) sin(phi) and W = -sqrt(rho+ rho-) cos(phi), "
   "so 4. is III.10's D_+- S_+- = -m c^2 sqrt(rho_-+/rho_+-) cos(phi) with a = m c^2/hbar",
   sp.simplify(Kpolar - sp.sqrt(r1*r2)*sp.sin(phi)) == 0 and sp.simplify(Wpolar + sp.sqrt(r1*r2)*sp.cos(phi)) == 0)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
