"""Do the half-turn matrices give true spin (S_z with eigenvalues +-hbar/2, in its eigenvector basis)?
1. IV.4b's basis is already the eigenvector basis of A_z (A_z is diagonal): readings 1,4 have axis +z, readings 2,3 axis -z.
   The A_i satisfy the Pauli algebra (anticommute, square 1): they are the Pauli matrices in a real 4-dim representation.
2. Dirac's spin operator along z is Sigma_z = -i A_x A_y. On the real readings A_x A_y is a quarter-turn with no real eigenvectors;
   -i A_x A_y (complexified) has eigenvalues +-1, each twice, and commutes with A_z and T: spin along the axis of motion is conserved.
3. What 'spin up/down along z' is on the real readings: the two senses of turning about z. On the plane wave along z, the four
   readings turn about z at the Compton rate a (half-angle: exp(-a s A_xA_y) or exp(+a s A_xA_y)), one sense per spin state;
   the content (flow, along z) does not move. The clock phase e^{-i a tau} of the 1+1 column IS this turning (i_is_spin_about_motion.py)."""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
Az = sp.diag(1, -1, -1, 1)
Ax = sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]])
Ay = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
T4 = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]])
A = (Ax, Ay, Az)
ok("1a. A_z is diagonal: IV.4b's readings are the eigenvector basis of the z direction (+z: readings 1,4; -z: readings 2,3)", Az == sp.diag(1, -1, -1, 1))
ok("1b. Pauli algebra: A_i A_j + A_j A_i = 2 delta_ij", all(sp.simplify(A[i]*A[j] + A[j]*A[i] - (2 if i == j else 0)*sp.eye(4)) == sp.zeros(4) for i in range(3) for j in range(3)))
Sz = -sp.I*Ax*Ay
ok("2a. Sigma_z = -i A_x A_y is Hermitian with eigenvalues {+1: 2, -1: 2}: spin up and spin down along z, each for both time senses",
   Sz.H == Sz and dict(Sz.eigenvals()) == {1: 2, -1: 2})
ok("2b. Sigma_z commutes with A_z and T (motion along z and the clock): spin along the axis of motion is conserved; it anticommutes with A_x, A_y",
   Sz*Az == Az*Sz and Sz*T4 == T4*Sz and Sz*Ax == -Ax*Sz and Sz*Ay == -Ay*Sz)
ok("2c. on the REAL readings A_x A_y has no real eigenvector (eigenvalues +-i): spin up/down is not a direction in the reading space but a sense of turning",
   dict((Ax*Ay).eigenvals()) == {sp.I: 2, -sp.I: 2})
# 3. plane wave along z, from i_is_spin_about_motion.py: Q4 = P (Re psi; -Im psi), psi = (e^{eta/2}, -i e^{-eta/2}) e^{-i a tau}
t, z, a, c, eta, s = sp.symbols('t z a c eta s', real=True)
tau = t*sp.cosh(eta) - z*sp.sinh(eta)/c
e = lambda i: sp.eye(4)[:, i]; P = sp.Matrix.hstack(e(0), e(2), e(3), -e(1))
psi = sp.Matrix([sp.exp(eta/2), -sp.I*sp.exp(-eta/2)])*sp.exp(-sp.I*a*tau)
Q4 = P*sp.Matrix.vstack(psi.applyfunc(sp.re), -psi.applyfunc(sp.im))
J = Ax*Ay
adv = Q4.subs(t, t + s/sp.cosh(eta))                        # clock advanced by s in proper time
rot_plus = (sp.cos(a*s)*sp.eye(4) + sp.sin(a*s)*J)*Q4; rot_minus = (sp.cos(a*s)*sp.eye(4) - sp.sin(a*s)*J)*Q4
simp = lambda M: M.applyfunc(lambda f: sp.simplify(f.rewrite(sp.exp)))
which = "+" if simp(adv - rot_plus) == sp.zeros(4, 1) else ("-" if simp(adv - rot_minus) == sp.zeros(4, 1) else "none")
ok(f"3a. advancing the clock by s turns the four readings about z by the half-angle a s (sense {which}): exp({which}a s A_xA_y) Q4", which != "none")
flow = sp.Matrix([sp.simplify((Q4.T*M*Q4)[0].rewrite(sp.exp)) for M in A])
target = sp.Matrix([-2*sp.cos(2*a*tau), -2*sp.sin(2*a*tau), 2*sp.sinh(eta)])
ok("3b. ONE real copy: its flow is (-2cos 2a tau, -2 sin 2a tau, 2 sinh eta): the segments circle the axis at the full angle 2a tau while the readings turn by a tau. "
   "At rest the walker goes round at c, frequency 2a = 2mc^2/hbar, radius c/2a = hbar/2mc: the Zitterbewegung orbit, as the content of one real column",
   simp(flow - target) == sp.zeros(3, 1))
# two copies a quarter clock apart (the complex column): the transverse flows cancel
Q4b = Q4.subs(t, t + sp.pi/(2*a*sp.cosh(eta)))
flow2 = sp.Matrix([sp.simplify((Q4b.T*M*Q4b)[0].rewrite(sp.exp)) for M in A])
ok("3b'. the partner copy a quarter clock later has the opposite transverse flow; the two together (= the complex column |psi|^2 etc.) flow along z only, steady: "
   "the steady particle is two circling walkers half a circle apart", simp(flow + flow2 - sp.Matrix([0, 0, 4*sp.sinh(eta)])) == sp.zeros(3, 1))
# the other spin state: conjugate psi (i -> -i) gives the opposite sense
Q4c = P*sp.Matrix.vstack(psi.applyfunc(sp.re), psi.applyfunc(sp.im))
advc = Q4c.subs(t, t + s/sp.cosh(eta))
other = "-" if simp(advc - (sp.cos(a*s)*sp.eye(4) - sp.sin(a*s)*J)*Q4c) == sp.zeros(4, 1) else ("+" if simp(advc - (sp.cos(a*s)*sp.eye(4) + sp.sin(a*s)*J)*Q4c) == sp.zeros(4, 1) else "none")
ok(f"3c. the other spin state (i -> -i) turns the other way (sense {other}): spin up/down along the motion = the two senses of turning of the readings about the axis", other != "none" and other != which)
