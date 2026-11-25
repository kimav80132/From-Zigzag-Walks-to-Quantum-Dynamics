"""
The end-to-end rule, tested: magnetism from two walks (3+1 lattice, c = hbar = 1, tick = spacing = 1).

Rule: a carrier path from a segment of walker 1 (displacement dx1^mu) to a segment of walker 2 (dx2^nu) turns walker 2's
copy pair by   kappa q1 q2 eta_{mu nu} dx1^mu dx2^nu x (carrier path weight).
Summed over all carrier paths from all of walker 1's segments: walker 2's turn per segment is q2 (A_t dt + A_z dz) with
    A_t = kappa q1 eta_tt (G * rho1),   A_z = kappa q1 eta_zz (G * j1),
rho1 = the density of walker 1's time segments (its column's density), j1 = the imbalance of its z-moves (psi^dag alpha_z psi),
G = the carrier walk's occupation (VI.2 (E), (G)): (6 - A) G = 6 source on the lattice.

Setup: steady line currents along z (plane waves in z, Gaussian in x,y), so the static G suffices. Walker 2 is IV.10's
three-coin automaton; its time turn per tick is exp(i q2 A_t) and each z-hop carries exp(+- i q2 A_z) (the rule on its
z-segments). Its transverse momentum gain per tick is measured from its lattice modes.

Checks:
  1. walker 1's z-segments average to its current: j1/rho1 = its group velocity (from the automaton's dispersion)
  2. walker 2's push per tick = q2 < -d_x Phi + alpha_z d_x A_z >  (the Lorentz force, VI.2 (D)), tick by tick
  3. with eta = diag(-1, +1): F(v1, v2) / F(0, 0) = 1 - v1 v2  (parallel currents attract magnetically; antiparallel repel more); tolerance 3%, the lattice's own Lorentz-check level of 2.
  4. with a Euclidean delta in place of eta (sign of the time-time entry flipped, then kappa flipped so like static charges still
     repel): the factor is 1 + v1 v2 - the wrong magnetism. The relative sign of electricity and magnetism is eta's.
"""
import numpy as np, scipy.sparse as sps, scipy.sparse.linalg as spla, sys
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all = ok_all and bool(cond); print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()

Lx = Ly = 72; Lz = 32; m = 0.5; w = 4.0; d = 14                     # box, mass, packet width, separation along x
cx = Lx//2; cy = Ly//2
X, Y = np.meshgrid(np.arange(Lx) - cx, np.arange(Ly) - cy, indexing='ij')

# carrier: z-uniform sources, so the 3d walk's occupation is z-independent: (4 - A_2d) G = 6 s on x,y, faces absorbing
n = Lx - 2; e1 = np.ones(n)
D1 = sps.diags([-e1[:-1], 2*e1, -e1[:-1]], [-1, 0, 1], format='csr'); I1 = sps.identity(n, format='csr')
K2 = (sps.kron(D1, I1) + sps.kron(I1, D1)).tocsr()
def carrier(src2d):
    v, info = spla.cg(K2, 6*src2d[1:-1, 1:-1].ravel(), rtol=1e-12, maxiter=20000); assert info == 0
    G = np.zeros((Lx, Ly)); G[1:-1, 1:-1] = v.reshape(n, n); return G

# the three-coin automaton (Dirac basis), z-hops with link phases
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.array([[1, 0], [0, -1]], complex); Z2 = np.zeros((2, 2))
alpha = [np.block([[Z2, s], [s, Z2]]) for s in (s1, s2, s3)]; beta = np.diag([1, 1, -1, -1]).astype(complex)
Pp = [(np.eye(4) + aa)/2 for aa in alpha]; Pm = [(np.eye(4) - aa)/2 for aa in alpha]
Mm = np.cos(m)*np.eye(4) - 1j*np.sin(m)*beta
def tick(psi, At=None, Az=None):                                      # psi[x, y, z, 4]
    for j in range(3):
        plus = psi @ Pp[j].T; minus = psi @ Pm[j].T
        if j == 2 and Az is not None:
            plus = plus*np.exp(1j*Az)[:, :, None, None]; minus = minus*np.exp(-1j*Az)[:, :, None, None]
        psi = np.roll(plus, 1, axis=j) + np.roll(minus, -1, axis=j)
    psi = psi @ Mm.T
    if At is not None: psi = psi*np.exp(1j*At)[:, :, None, None]
    return psi
kx = 2*np.pi*np.fft.fftfreq(Lx)
def px(psi):
    ph = np.fft.fft(psi, axis=0); wgt = (np.abs(ph)**2).sum(axis=(1, 2, 3)); return (wgt*kx).sum()/wgt.sum()
def U_of(k):                                                          # one tick for a mode (kx, ky, kz)
    U = np.eye(4, dtype=complex)
    for j in range(3): U = (Pp[j]*np.exp(-1j*k[j]) + Pm[j]*np.exp(1j*k[j])) @ U
    return Mm @ U
# positive-energy packet: plane wave e^{i kz z} times a Gaussian in x,y, projected mode by mode on the two positive bands
KX, KY = np.meshgrid(kx, 2*np.pi*np.fft.fftfreq(Ly), indexing='ij')
def packet(x0, nz):
    kz = 2*np.pi*nz/Lz
    env = np.exp(-((X - x0)**2 + Y**2)/(4*w**2))
    spin = np.array([1, 0, 0, 0], complex)
    ph = np.fft.fft2(env)[..., None]*spin                              # (Lx, Ly, 4) in k_perp
    out = np.zeros_like(ph)
    for i in range(Lx):
        for jj in range(Ly):
            if abs(ph[i, jj]).max() < 1e-14: continue
            lam, V = np.linalg.eig(U_of((KX[i, jj], KY[i, jj], kz))); pos = np.angle(lam) < 0
            Vi = np.linalg.inv(V); out[i, jj] = V @ (pos*(Vi @ ph[i, jj]))
    p2 = np.fft.ifft2(out, axes=(0, 1))
    psi = p2[:, :, None, :]*np.exp(1j*kz*np.arange(Lz))[None, None, :, None]
    return psi/np.sqrt((np.abs(psi)**2).sum())
def group_velocity(x0, nz):                                            # d omega / d kz of the positive bands, averaged over the packet's transverse modes
    kz = 2*np.pi*nz/Lz; h = 1e-5
    env = np.exp(-((X - x0)**2 + Y**2)/(4*w**2)); wk = np.abs(np.fft.fft2(env))**2
    def om(k): return np.sort(-np.angle(np.linalg.eig(U_of(k))[0]))[-2:]          # the two positive bands
    tot = 0.0; vs = 0.0
    for i in range(Lx):
        for jj in range(Ly):
            if wk[i, jj] < 1e-12*wk.max(): continue
            dv = (om((KX[i, jj], KY[i, jj], kz + h)) - om((KX[i, jj], KY[i, jj], kz - h))).mean()/(2*h)
            vs += wk[i, jj]*dv; tot += wk[i, jj]
    return vs/tot
az = alpha[2]
def sources(psi):                                                     # per (x,y), summed over z: time-segment density and z-current
    rho = (np.abs(psi)**2).sum(axis=(2, 3))
    jz = np.einsum('xyzi,ij,xyzj->xy', psi.conj(), az, psi).real
    return rho, jz

q1 = q2 = 1.0; kappa = 0.4
def run(n1, n2, metric_tt, kap):
    psi1 = packet(-d//2, n1); psi2 = packet(d - d//2, n2)
    rho1, j1 = sources(psi1)
    Gr = carrier(rho1); Gj = carrier(j1)
    At = kap*q1*metric_tt*Gr; Az = kap*q1*(+1.0)*Gj                    # A_mu = kappa q1 eta_{mu mu} (G * segment density)
    Phi = -At                                                         # walker 2's time turn exp(i q2 A_t) = exp(-i q2 Phi)
    dPhi = (np.roll(Phi, -1, 0) - np.roll(Phi, 1, 0))/2; dAz = (np.roll(Az, -1, 0) - np.roll(Az, 1, 0))/2
    rows = []
    for t in range(3):
        r2 = (np.abs(psi2)**2).sum(axis=(2, 3)); a2 = np.einsum('xyzi,ij,xyzj->xy', psi2.conj(), az, psi2).real
        pred = q2*((-dPhi*r2).sum() + (dAz*a2).sum())
        p0 = px(psi2); psi2 = tick(psi2, q2*At, q2*Az); rows.append((px(psi2) - p0, pred))
    return rows, (j1.sum()/rho1.sum())

nz = 4
v_g = group_velocity(-d//2, nz)
res = {}
for key, (n1, n2) in {"static": (0, 0), "parallel": (nz, nz), "antiparallel": (nz, -nz)}.items():
    res[key] = run(n1, n2, -1.0, kappa)
v1 = res["parallel"][1]
ok("1. walker 1's z-segments average to its current: j1/rho1 = %.4f against the automaton's group velocity d omega/d k_z averaged over the packet's transverse modes = %.4f, k_z = 2 pi %d/%d" % (v1, v_g, nz, Lz), abs(v1 - v_g) < 3e-3)
worst = max(abs(dp/pr - 1) for rows, _ in res.values() for dp, pr in rows)
ok("2. walker 2's push per tick = q2 < -d_x Phi + alpha_z d_x A_z > (the Lorentz force from the rule's potentials), all runs, 3 ticks each: max deviation %.2f%%" % (100*worst), worst < 0.02)
F0 = res["static"][0][0][0]; Fp = res["parallel"][0][0][0]; Fa = res["antiparallel"][0][0][0]
rp, ra = Fp/F0, Fa/F0
ok("3. with eta: like static charges repel (F_x = %+.3e, pushing walker 2 away); parallel currents F/F0 = %.4f vs 1 - v1 v2 = %.4f; antiparallel %.4f vs 1 + v1^2 = %.4f - parallel currents attract magnetically"
   % (F0, rp, 1 - v1*v1, ra, 1 + v1*v1), F0 > 0 and abs(rp - (1 - v1*v1)) < 0.03 and abs(ra - (1 + v1*v1)) < 0.03)
# Euclidean delta: time-time entry +1; to keep like static charges repelling, kappa -> -kappa
E0 = run(0, 0, +1.0, -kappa)[0][0][0]; Ep = run(nz, nz, +1.0, -kappa)[0][0][0]
ok("4. with a Euclidean delta (kappa flipped so like static charges still repel, F_x = %+.3e): parallel currents F/F0 = %.4f vs 1 + v1 v2 = %.4f - the wrong magnetism; the relative sign of electricity and magnetism is eta's"
   % (E0, Ep/E0, 1 + v1*v1), E0 > 0 and abs(Ep/E0 - (1 + v1*v1)) < 0.03)
print("ALL PASS" if ok_all else "SOME FAILED"); raise SystemExit(0 if ok_all else 1)
