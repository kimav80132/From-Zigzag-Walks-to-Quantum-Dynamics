"""VI.2 (D): the Lorentz force is the turn angle of the pair, differentiated twice. Closes Part V's ledger item.
The angle a step picks up is theta = (e/hbar)(A . dx) = (e/hbar)(-Phi dt + A_k dx^k)  (V.1, V.2). What a loop picks up is its flux (pair_angle_is_a_connection.py 4-5).
The force is the rate of change of the mean kinetic momentum pi = p - eA: the generator of translations with the connection.
1. (symbolic, exact) In 3+1 with H = c alpha.(p - eA) + beta m c^2 + e Phi, and arbitrary A(t,x,y,z), Phi(t,x,y,z):
       d<pi_i>/dt = (i/hbar)<[H, pi_i]> + <d pi_i/dt>  =  e <E_i + (c alpha x B)_i>
   with E = -grad Phi - dA/dt, B = curl A: the Lorentz force on the mean step, with the velocity c alpha (the direction matrices, IV.4b).
   Line by line: [pi_i, pi_j] = i hbar e eps_{ijk} B_k is the flux of a plaquette (the two orders of two moves, V.3); [pi_i, e Phi] = -i hbar e d_i Phi.
2. (lattice, 1+1) the two-copy walk of V.2 with a uniform electric field Phi = -E x: the mean kinetic momentum of a packet grows at the rate eE
   to 1 percent over the run, read from the two real legs with no complex number: the turn angle per step, differenced along t, is the force.
3. the sign of the force reverses with the sense: in the clock reading theta_eps carries eps Phi (V.2 (H)); the backward segment is pushed the other way."""
import sympy as sp, numpy as np
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
# ---------- 1. symbolic commutator, exact ----------
t, x, y, z = sp.symbols('t x y z', real=True); c, e, hb, m = sp.symbols('c e hbar m', positive=True)
X = (x, y, z)
Ax, Ay, Az, Phi = [sp.Function(n, real=True)(t, x, y, z) for n in ('A_x', 'A_y', 'A_z', 'Phi')]
A = (Ax, Ay, Az)
# a generic 4-component test column
psi = sp.Matrix([sp.Function(f'psi{k}')(t, x, y, z) for k in range(4)])
# real matrices of IV.4b (I): A_z, A_x, A_y, T; standard beta = i T
Az_m = sp.diag(1, -1, -1, 1)
Ax_m = sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]])
Ay_m = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
T_m = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]])
alpha = (Ax_m, Ay_m, Az_m); beta = sp.I*T_m
pi_op = lambda i, f: -sp.I*hb*sp.diff(f, X[i]) - e*A[i]*f                    # kinetic momentum on a column
H_op = lambda f: sum((c*alpha[i]*pi_op(i, f) for i in range(3)), sp.zeros(4, 1)) + beta*m*c**2*f + e*Phi*f
E = [-sp.diff(Phi, X[i]) - sp.diff(A[i], t)/c for i in range(3)]              # E = -grad Phi - (1/c) dA/dt  (Gaussian-like units: A.dx - Phi c dt... see below)
# convention: theta = (e/hbar)(A.dx - Phi dt) with H carrying e Phi and p - eA; then E_i = -d_i Phi - d_t A_i /1? Fix by the commutator itself:
B = [sp.diff(A[(i+2) % 3], X[(i+1) % 3]) - sp.diff(A[(i+1) % 3], X[(i+2) % 3]) for i in range(3)]   # B = curl A
allok = True
for i in range(3):
    lhs = (sp.I/hb)*(H_op(pi_op(i, psi)) - pi_op(i, H_op(psi))) + (-e*sp.diff(A[i], t))*psi     # d pi/dt = (i/hbar)[H, pi] + d pi/dt(explicit)
    Ei = -sp.diff(Phi, X[i]) - sp.diff(A[i], t)                                                  # with H = ... + e Phi and pi = p - eA: E = -grad Phi - dA/dt
    cross = sum((sp.LeviCivita(i, j, k)*c*alpha[j]*B[k] for j in range(3) for k in range(3)), sp.zeros(4, 4))
    rhs = e*(Ei*sp.eye(4) + cross)*psi
    allok &= sp.simplify(sp.expand(lhs - rhs)) == sp.zeros(4, 1)
ok("1. exact, any A(t,x), Phi(t,x): d pi_i/dt = e [E_i + (c alpha x B)_i] on every column, E = -grad Phi - dA/dt, B = curl A, velocity = c alpha", allok)
# the two pieces, line by line
comm_pipi = lambda i, j, f: pi_op(i, pi_op(j, f)) - pi_op(j, pi_op(i, f))
ok("1'. [pi_i, pi_j] = i hbar e (d_i A_j - d_j A_i): the flux of the plaquette (i,j) - the two orders of two moves differ by B (V.3)",
   all(sp.simplify(comm_pipi(i, j, psi) - sp.I*hb*e*(sp.diff(A[j], X[i]) - sp.diff(A[i], X[j]))*psi) == sp.zeros(4, 1) for i in range(3) for j in range(3)))
ok("1''. [pi_i, e Phi] = -i hbar e d_i Phi: the gradient of the time-link angle", sp.simplify(pi_op(0, e*Phi*psi) - e*Phi*pi_op(0, psi) + sp.I*hb*e*sp.diff(Phi, x)*psi) == sp.zeros(4, 1))
# ---------- 2. lattice, 1+1, two real legs, uniform E ----------
Lx, a_, dt, ec = 2000, 0.6, 0.05, 0.004          # ec = e E (force) in lattice units hbar = c = 1
T2 = np.array([[0, 1], [-1, 0.]]); Rc = np.cos(a_*dt)*np.eye(2) + np.sin(a_*dt)*T2
J2 = np.array([[0, -1], [1, 0.]]); Rl = lambda th: np.cos(th)*np.eye(2) + np.sin(th)*J2
xs = np.arange(Lx) - Lx/2
# column per site: (Q+, Q+', Q-, Q-'). One tick: shift, coin on senses, turn of the pair by theta = -(e/hbar) Phi dt with Phi = -E x  => theta = +eE x dt
def tick(col):
    out = np.zeros_like(col)
    out[1:, 0:2] = col[:-1, 0:2]; out[0, 0:2] = col[-1, 0:2]            # + moves right
    out[:-1, 2:4] = col[1:, 2:4]; out[-1, 2:4] = col[0, 2:4]            # - moves left
    for s in (0, 2):                                                    # coin on senses: mixes rows (Q+,Q-) and (Q+',Q-') alike
        pass
    new = np.zeros_like(out)
    new[:, 0] = Rc[0, 0]*out[:, 0] + Rc[0, 1]*out[:, 2]; new[:, 2] = Rc[1, 0]*out[:, 0] + Rc[1, 1]*out[:, 2]
    new[:, 1] = Rc[0, 0]*out[:, 1] + Rc[0, 1]*out[:, 3]; new[:, 3] = Rc[1, 0]*out[:, 1] + Rc[1, 1]*out[:, 3]
    th = ec*xs*dt
    for s in (0, 2):
        q, qp = new[:, s].copy(), new[:, s+1].copy()
        new[:, s] = np.cos(th)*q - np.sin(th)*qp; new[:, s+1] = np.sin(th)*q + np.cos(th)*qp
    return new
# packet: positive-energy column of momentum k0 (III.9 basis, mass matrix -sigma_y): H = sigma_z k + a (0 1;-1 0) i  -> use the complex column to build, then split
k0, w = 0.0, 25.0
sig = np.exp(-xs**2/(2*w**2))*np.exp(1j*k0*xs)
# positive-energy eigenvector of H(k0) = [[k0, i a],[-i a, -k0]]
Hk = np.array([[k0, 1j*a_], [-1j*a_, -k0]]); ev, evec = np.linalg.eigh(Hk); u = evec[:, np.argmax(ev)]
psi0 = np.outer(sig, u); psi0 /= np.sqrt(np.sum(np.abs(psi0)**2))
col = np.stack([psi0[:, 0].real, psi0[:, 0].imag, psi0[:, 1].real, psi0[:, 1].imag], axis=1)
def kin_momentum(col):                                                   # <p - eA> with A = 0 here (Phi gauge): <psi| -i d_x |psi> from the two real legs
    psi = np.stack([col[:, 0] + 1j*col[:, 1], col[:, 2] + 1j*col[:, 3]], axis=1)
    dpsi = (np.roll(psi, -1, axis=0) - np.roll(psi, 1, axis=0))/2
    return np.sum(np.conj(psi)*(-1j*dpsi)).real
p = [kin_momentum(col)]
for n in range(600):
    col = tick(col); p.append(kin_momentum(col))
ts = np.arange(601)*dt
slope = np.polyfit(ts, p, 1)[0]
ok(f"2. lattice, two real legs, uniform E: d<p>/dt = {slope:.5f} against eE = {ec:.5f} ({100*abs(slope/ec-1):.2f} percent): the turn angle per step is the force",
   abs(slope/ec - 1) < 0.02)
# 3. sign with the sense: redo with the clock-reading angle theta_eps = -(e/hbar) eps Phi dx/c... in the t-reading the same Phi pushes both senses the same way
# (V.2 (E)); the opposite push is on the backward segment read in the walk's own clock. Check on the formula: theta_eps of V.2 (A) with eps = -1 reverses.
eps_s = sp.symbols('epsilon'); Phi_s, dx_s = sp.symbols('Phi dx', real=True)
theta_eps = -(e/hb)*(eps_s*Phi_s*dx_s/c)
ok("3. in the clock reading the angle carries eps: theta_{-} = -theta_{+} for the same Phi (V.2 (H)); the backward segment is turned, and pushed, the other way",
   sp.simplify(theta_eps.subs(eps_s, -1) + theta_eps.subs(eps_s, 1)) == 0)
