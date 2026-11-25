"""Can the field enter as an IMAGINARY skew, carried by the noise?
1. x-clock form of Dirac with Phi and A_x (1+1):  c d_x Q_e + e (d_t + i e Phi/hbar) Q_e - i c (e A_x/hbar) Q_e = a Q_-e.
   The field terms are imaginary and DIAGONAL (loss-like): Phi skewed by the sense (-i eps e Phi/hbar), A_x the same for both senses.
   Neither is an exchange skew.
2. The attachment's exchange skew made imaginary, a -> a +- i delta: adds -i delta (1 + sigma_x) in the t-reading = potential + a
   pseudoscalar term hbar delta sigma_x (Hermitian, probability conserved, but not pure EM).
3. Part II: a phase e^{i e z.A/hbar} on each step. Expanding <e^{-z.d} e^{i e z.A/hbar}>: covariance <z z> = 2 D dtau eta unchanged;
   the drift shift  Delta b^mu = i e <z^mu z^nu> A_nu /(hbar dtau) = (i e/m) A^mu  and the A^2 term both come THROUGH the noise matrix.
4. Part II's own signature (time kick i zeta): the electric part of the phase is a REAL weight on the time step,
   e^{+ e zeta sqrt(2 D dtau) Phi/(hbar c)}.  Weyl's real scaling is London's phase, seen through the imaginary time kick.
5. Part III real walk: the imaginary diagonal needs the second copy: -eps delta J for Phi, +c (e A_x/hbar) J for A_x, J^2 = -1.
   The total length over both copies is conserved; the pair of copies is Dirac with the field.
"""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
t, x = sp.symbols('t x', real=True); c, a, hb, e, m, D, dtau = sp.symbols('c a hbar e m D dtau', positive=True)
Phi, Ax = sp.Function('Phi', real=True)(t, x), sp.Function('Ax', real=True)(t, x)
Qp, Qm = sp.Function('Qp')(t, x), sp.Function('Qm')(t, x)
# 1. Dirac (t-reading, original basis): d_t Q = -c A (d_x - i e Ax/hbar) Q + a T Q - i (e Phi/hbar) Q ; solve for c d_x
dQp = -c*(Qp.diff(x) - sp.I*e*Ax/hb*Qp) + a*Qm - sp.I*e*Phi/hb*Qp
dQm = +c*(Qm.diff(x) - sp.I*e*Ax/hb*Qm) - a*Qp - sp.I*e*Phi/hb*Qm
xclock_p = sp.expand(sp.solve(sp.Eq(Qp.diff(t), dQp), c*Qp.diff(x))[0]) if False else sp.expand(-Qp.diff(t) + a*Qm - sp.I*e*Phi/hb*Qp + sp.I*c*e*Ax/hb*Qp)
xclock_m = sp.expand(+Qm.diff(t) + a*Qp + sp.I*e*Phi/hb*Qm + sp.I*c*e*Ax/hb*Qm)
ok("1. x-clock: c d_x Q+ = -d_t Q+ + a Q- - i(e Phi/hbar) Q+ + i c(e Ax/hbar) Q+ ;  c d_x Q- = +d_t Q- + a Q+ + i(e Phi/hbar) Q- + i c(e Ax/hbar) Q-",
   sp.simplify(sp.expand(xclock_p.subs(Qp.diff(t), dQp) - c*Qp.diff(x))) == 0 and sp.simplify(sp.expand(xclock_m.subs(Qm.diff(t), dQm) - c*Qm.diff(x))) == 0)
print("     -> the field sits on the DIAGONAL, imaginary: Phi with opposite sign for the two senses, Ax with the same sign; no exchange skew.")
# 2. imaginary exchange skew
d = sp.symbols('delta', real=True)
M = -c*sp.diag(1, -1)*0 + a*sp.Matrix([[0, 1], [-1, 0]]) - sp.I*d*sp.Matrix([[1, 1], [1, 1]])     # local part of the t-reading
H = sp.I*hb*M
ok("2. imaginary exchange skew: H_local = i hbar a T + hbar delta (1 + sigma_x) is Hermitian; the extra hbar delta sigma_x anticommutes with sigma_z and sigma_y (a pseudoscalar coupling)",
   sp.simplify(H - H.H) == sp.zeros(2) and (sp.Matrix([[0, 1], [1, 0]])*sp.diag(1, -1) + sp.diag(1, -1)*sp.Matrix([[0, 1], [1, 0]])) == sp.zeros(2))
# 3. Part II: phase on each step, Gaussian step in one direction (the covariance does the work)
A0, b, p = sp.symbols('A b p', real=True)
var = 2*D*dtau; k = -p + sp.I*e*A0/hb                                      # e^{-z p} e^{i e z A/hbar} = e^{z k}
mgf = sp.exp(b*dtau*k + var*k**2/2)                                        # Gaussian step: mean b dtau, variance 2 D dtau
mom = sp.expand(sp.series(mgf, dtau, 0, 2).removeO())                      # first order in dtau
lin = mom.coeff(p, 1); const = mom.coeff(p, 0)
ok(f"3a. first-order (drift) coefficient: {sp.simplify(lin)} = -(b dtau + i e <z^2> A/hbar) : the phase shifts the drift by i e <zz> A/(hbar dtau) = 2 i e D A/hbar",
   sp.simplify(lin - (-(b*dtau + sp.I*e*var*A0/hb))) == 0)
ok(f"3b. zeroth-order term: 1 + i e b dtau A/hbar - e^2 <z^2> A^2/(2 hbar^2): the A^2 term comes from the covariance",
   sp.simplify(const - (1 + sp.I*e*b*dtau*A0/hb - e**2*var*A0**2/(2*hb**2))) == 0)
ok("3c. with D = hbar/2m the drift shift is (i e/m) A", sp.simplify((sp.I*e*var*A0/hb/dtau).subs(D, hb/(2*m)) - sp.I*e*A0/m) == 0)
# 4. imaginary time kick: the electric phase is a real weight
zeta = sp.symbols('zeta', real=True); Phi0 = sp.symbols('Phi0', real=True)
z0 = sp.I*zeta*sp.sqrt(var); A_0 = -Phi0/c                                   # lower index, eta = diag(-1, 1)
ph = sp.exp(sp.I*e*z0*A_0/hb)
ok(f"4. phase on the time kick: exp(i e z^0 A_0/hbar) = {sp.simplify(ph)}: real", sp.im(sp.simplify(ph)) == 0 and sp.simplify(ph - sp.exp(e*zeta*sp.sqrt(var)*Phi0/(hb*c))) == 0)
# 5. two real copies
Q = sp.Matrix([sp.Function(f'q{i}')(t, x) for i in range(4)])
K = lambda P_, R_: sp.kronecker_product(P_, R_); eps2 = sp.Matrix([[0, 1], [-1, 0]])
A4, T4, J4 = K(sp.diag(1, -1), sp.eye(2)), K(sp.Matrix([[0, 1], [-1, 0]]), sp.eye(2)), K(sp.eye(2), -eps2)   # J (q0,q1) = (-q1, q0): J acts as i on q0 + i q1
dQ = -c*A4*Q.diff(x) + a*T4*Q - (e*Phi/hb)*J4*Q + c*(e*Ax/hb)*A4*J4*Q
dens = sp.expand(sum(2*Q[i]*dQ[i] for i in range(4)) + c*sp.diff((Q.T*A4*Q)[0], x))
ok("5a. two real copies with -(e Phi/hbar) J and +c (e Ax/hbar) A J: total length conserved", sp.simplify(dens) == 0)
psi = sp.Matrix([Q[0] + sp.I*Q[1], Q[2] + sp.I*Q[3]])         # one complex column from the two copies
dpsi = sp.Matrix([dQ[0] + sp.I*dQ[1], dQ[2] + sp.I*dQ[3]])
dirac = -c*sp.diag(1, -1)*(psi.diff(x) - sp.I*e*Ax/hb*psi) + a*sp.Matrix([[0, 1], [-1, 0]])*psi - sp.I*e*Phi/hb*psi
ok("5b. the two copies read as one complex column are Dirac with the field: d_t psi = -c A (d_x - i e Ax/hbar) psi + a T psi - i (e Phi/hbar) psi",
   sp.simplify(sp.expand(dpsi - dirac)) == sp.zeros(2, 1))
