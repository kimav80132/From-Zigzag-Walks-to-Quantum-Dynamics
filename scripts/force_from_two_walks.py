"""
VI.2 (G): the force, run as one process (3+1, lattice spacing 1, tick 1, c = hbar = 1).

Walker 1 sits at the centre of a 96^3 box; its density is the source. The carrier is Part II's walker at m = 0 (VI.2 (E)):
a simple random walk from the source, one of 6 neighbours per step, absorbed at the box faces. Its total occupation
    G(x) = sum_n P^n delta_0 ,  i.e.  (1 - P) G = delta_0  <=>  -(lattice Laplacian) G = 6 delta_0
(III.16's kernel summed over its clock) is the potential, Phi = e G. The link angle at each site is V.2's time-link angle,
theta_t = -e Phi per tick. Walker 2 is IV.10's three-coin automaton with the mass turn (x-coin, y-coin, z-coin, then
cos(m) - i beta sin(m)), a packet at rest at distance d along x, positive-energy projected with the automaton's own
eigenvectors, turned by theta_t each tick. Its momentum is read from its own lattice modes (FFT).

Checks:
  1. the carrier's occupation, exact (sparse CG) against a Monte Carlo of independent walks, shell by shell
  2. Coulomb from the walk: G(r) - G(16) against (6/4pi)(1/r - 1/16); the gradient against (6/4pi)/r^2
  3. the push on walker 2 is -e <d_x Phi>, tick by tick, and points away from walker 1 (like charges repel)
  4. reversing walker 2's charge reverses the push (toward walker 1)
  5. the walker is a cloud: F d^2 / (e^2 6/4pi) equals the cloud's weight within d of its own centre
     (a point source's field averaged over a round cloud: the mean-value property of a harmonic field)
"""
import numpy as np, time, sys
import scipy.sparse as sps, scipy.sparse.linalg as spla
def ok(msg, cond): print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()

Lb = 96; c0 = Lb//2                                   # box 0..95, source at (48,48,48), faces absorbing
n_in = Lb - 2                                         # interior sites 1..94

# ---------- 1a. exact occupation: (6 - A) G = 6 delta on the interior, G = 0 on the faces ----------
t0 = time.time()
e1 = np.ones(n_in); D1 = sps.diags([-e1[:-1], 2*e1, -e1[:-1]], [-1, 0, 1], format='csr')   # 1d: 2 - (shift up + shift down)
I1 = sps.identity(n_in, format='csr')
K = (sps.kron(sps.kron(D1, I1), I1) + sps.kron(sps.kron(I1, D1), I1) + sps.kron(sps.kron(I1, I1), D1)).tocsr()   # 6 - A
rhs = np.zeros(n_in**3); src = ((c0 - 1)*n_in + (c0 - 1))*n_in + (c0 - 1); rhs[src] = 6.0
Gv, info = spla.cg(K, rhs, rtol=1e-11, maxiter=20000)
G = np.zeros((Lb, Lb, Lb)); G[1:-1, 1:-1, 1:-1] = Gv.reshape(n_in, n_in, n_in)
print("      exact occupation by conjugate gradients: %d unknowns, info %d, %.1f s, G(source) = %.4f (lattice value 6 x 0.2527 = 1.516 for an infinite lattice)" % (n_in**3, info, time.time() - t0, G[c0, c0, c0]))

# ---------- 1b. Monte Carlo of the carrier walk, shell-averaged occupation ----------
rng = np.random.default_rng(7); Nw = 200_000; chunk = 50_000
shell_edges = np.arange(0, 21); occ_mc = np.zeros(len(shell_edges) - 1)
t1 = time.time()
moves = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=np.int16)
for start in range(0, Nw, chunk):
    pos = np.full((chunk, 3), c0, dtype=np.int16)
    while len(pos):
        r = np.sqrt(((pos - c0).astype(np.float32)**2).sum(1))
        sh = r.astype(np.int64); keep = sh < len(occ_mc)
        occ_mc += np.bincount(sh[keep], minlength=len(occ_mc))[:len(occ_mc)]
        pos = pos + moves[rng.integers(0, 6, len(pos))]
        alive = ((pos > 0) & (pos < Lb - 1)).all(1); pos = pos[alive]
occ_mc /= Nw
# exact, summed over the same shells
X, Y, Z = np.meshgrid(*(np.arange(Lb) - c0,)*3, indexing='ij'); R = np.sqrt(X**2 + Y**2 + Z**2)
occ_ex = np.array([G[(R >= lo) & (R < lo + 1)].sum() for lo in shell_edges[:-1]])
rel = np.abs(occ_mc[2:16]/occ_ex[2:16] - 1)
ok("1. the carrier's occupation: exact solve against %d Monte Carlo walks (%.0f s), shells r = 2..15: max deviation %.1f%% (mean %.1f%%)" % (Nw, time.time() - t1, 100*rel.max(), 100*rel.mean()), rel.max() < 0.03)

# ---------- 2. Coulomb from the walk ----------
cC = 6/(4*np.pi)
rs = [4, 6, 8, 12]
ratios = [(G[c0 + r, c0, c0] - G[c0 + 16, c0, c0])/(cC*(1/r - 1/16)) for r in rs]
grad = [(G[c0 + d - 1, c0, c0] - G[c0 + d + 1, c0, c0])/2/(cC/d**2) for d in (8, 12, 16)]
ok("2. Coulomb's law from the walk: [G(r) - G(16)] / [(6/4pi)(1/r - 1/16)] = %s at r = %s (the box shifts G by a nearly uniform amount, which the differences remove); -d_x G / [(6/4pi)/d^2] = %s at d = 8, 12, 16 (closer in, the lattice's own short-distance correction: 6%% at d = 6)"
   % (np.round(ratios, 3).tolist(), rs, np.round(grad, 3).tolist()), all(abs(q - 1) < 0.04 for q in ratios) and all(abs(q - 1) < 0.04 for q in grad))

# ---------- the second walker: IV.10's three-coin automaton ----------
m = 0.6; w = 5.0
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.array([[1, 0], [0, -1]], complex); Z2 = np.zeros((2, 2))
alpha = [np.block([[Z2, s]] + [[s, Z2]]) for s in (s1, s2, s3)]
beta = np.diag([1, 1, -1, -1]).astype(complex)
Pp = [(np.eye(4) + a)/2 for a in alpha]; Pm = [(np.eye(4) - a)/2 for a in alpha]
Mmass = np.cos(m)*np.eye(4) - 1j*np.sin(m)*beta
def tick(psi, phase=None):                             # psi[x, y, z, 4]
    for j in range(3):
        plus = psi @ Pp[j].T; minus = psi @ Pm[j].T
        psi = np.roll(plus, 1, axis=j) + np.roll(minus, -1, axis=j)
    psi = psi @ Mmass.T
    if phase is not None: psi = psi*phase[..., None]
    return psi
kk = 2*np.pi*np.fft.fftfreq(Lb)
def momentum_x(psi):
    ph = np.fft.fftn(psi, axes=(0, 1, 2)); wgt = (np.abs(ph)**2).sum(-1)
    return (wgt*kk[:, None, None]).sum()/wgt.sum()
# positive-energy projector of the automaton, mode by mode: U(p) = M Cz Cy Cx, C_j = P+ e^{-i p_j} + P- e^{+i p_j}
t2 = time.time()
KX, KY, KZ = np.meshgrid(kk, kk, kk, indexing='ij')
Up = np.broadcast_to(np.eye(4, dtype=complex), (Lb, Lb, Lb, 4, 4)).copy()
for j, Kj in enumerate((KX, KY, KZ)):
    Cj = Pp[j][None, None, None]*np.exp(-1j*Kj)[..., None, None] + Pm[j][None, None, None]*np.exp(1j*Kj)[..., None, None]
    Up = Cj @ Up
Up = Mmass[None, None, None] @ Up
lam, V = np.linalg.eig(Up.reshape(-1, 4, 4))
posE = (np.angle(lam) < 0)                            # lambda = e^{-i omega}, omega > 0
Vinv = np.linalg.inv(V)
Ppos = np.einsum('nij,nj,njk->nik', V, posE.astype(complex), Vinv).reshape(Lb, Lb, Lb, 4, 4)
del Up, V, Vinv, lam
print("      automaton's positive-energy projector built mode by mode (%.0f s); %d of %d modes have two positive bands" % (time.time() - t2, (posE.sum(1) == 2).sum(), posE.shape[0]))
def packet(d):
    env = np.exp(-((X - d)**2 + Y**2 + Z**2)/(4*w**2))            # density width w
    psi = np.zeros((Lb, Lb, Lb, 4), complex); psi[..., 0] = env    # spin up, upper (positive-mass) component
    ph = np.fft.fftn(psi, axes=(0, 1, 2)); ph = np.einsum('xyzij,xyzj->xyzi', Ppos, ph)
    psi = np.fft.ifftn(ph, axes=(0, 1, 2)); return psi/np.sqrt((np.abs(psi)**2).sum())
dPhix = (np.roll(G, -1, 0) - np.roll(G, 1, 0))/2                    # d_x G, central difference

# ---------- 3, 4. the push, tick by tick, and its sign ----------
e = 0.3; d = 10; nt = 6
for qs, name in ((+1, "like charges"), (-1, "opposite charges")):
    psi = packet(d); phase = np.exp(-1j*(qs*e)*(e*G)); rows = []
    for n in range(nt):
        rho = (np.abs(psi)**2).sum(-1); expected = -(qs*e)*(e*(rho*dPhix).sum())     # -e2 <d_x Phi>, Phi = e1 G
        p0 = momentum_x(psi); psi = tick(psi, phase); p1 = momentum_x(psi)
        rows.append((p1 - p0, expected))
    rel = max(abs(dp/ex - 1) for dp, ex in rows); sign = np.sign(rows[0][0])
    if qs > 0:
        ok("3. %s: walker 2's momentum gain per tick = -e <d_x Phi> to %.1f%% over %d ticks (gain %.3e, expected %.3e at the first tick); it points %s walker 1"
           % (name, 100*rel, nt, rows[0][0], rows[0][1], "away from" if sign > 0 else "toward"), rel < 0.03 and sign > 0)
    else:
        ok("4. %s: the same push reversed, gain %.3e per tick, toward walker 1 (to %.1f%% of -e <d_x Phi>)" % (name, rows[0][0], 100*rel), rel < 0.03 and sign < 0)

# ---------- 5. the walker is a cloud ----------
out = []
for d in (8, 12, 16):
    psi = packet(d); rho = (np.abs(psi)**2).sum(-1)
    F = -e*(e*(rho*dPhix).sum())                                    # force at t = 0
    Rc = np.sqrt((X - d)**2 + Y**2 + Z**2); inside = (rho*(Rc < d)).sum()
    out.append((d, F*d*d/(e*e*cC), inside))
okc = all(abs(fd/ins - 1) < 0.06 for _, fd, ins in out)
ok("5. the walker is a cloud: F d^2 / (e^2 6/4pi) against the cloud's weight within d of its centre: %s - the field of a point source averaged over a round cloud is the field at the centre on the part inside d (the mean-value property of a harmonic field, which the walk's occupation has)"
   % ", ".join("d = %d: %.3f vs %.3f" % r for r in out), okc)
