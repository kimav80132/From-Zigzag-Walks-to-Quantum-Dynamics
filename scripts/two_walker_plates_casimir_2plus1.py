"""
Two plates made of walkers, coupled through the end-to-end rule: the force from the plates' own fluctuations
(2+1, near field; c = hbar = 1, tick = spacing = 1).

Why 2+1. In 1+1 a neutral plate has no field outside it (Gauss's law in one dimension), so two neutral 1+1 plates do not
feel each other. Here the plates are infinite along y and stacked along z. Per transverse momentum k the carrier's
occupation falls off as e^{-kappa |z|}, and that is what reaches the other plate.

THE PLATE. The 2+1 walk at transverse momentum p (conjugate to y), one tick:
    U(p) = M S C(p),   C(p) = cos p 1 - i sin p sigma_x       the y-step (+1 or -1 by the y-sense)
                       S: sense 1 steps z -> z+1, sense 2 steps z -> z-1
                       M(z) = [[cos m, -sin m], [sin m, cos m]]  the turn; m = 0.30 in the slab (6 sites), pi/4 outside.
For this walk a uniform region's gap is min(m, pi/2 - m) (check 1), so pi/4 is the widest vacuum it allows.
The slab is where turning is weak: walkers there have a lower threshold and are bound in it.
Lowest state: every column with quasi-energy below E_F is filled, at every p. That is the whole sea (all negative
quasi-energies) and a band of conduction walkers bound in the slab, |p| < p_F = 0.05. Each tick moves y and z by one,
so a walker never leaves its sublattice of the checkerboard; the plate holds walkers on both, which shows up as an exact
repeat of the conduction pocket at p ~ pi (U(p + pi) = Pi U(p) Pi, Pi = (-1)^z).

LINES. For transverse momentum k, a line runs from each occupied conduction column o at p to each empty column e at
p + k or p - k. Its weight is the site overlap N(z) = sum_c conj(e_c(z)) o_c(z); its rate is Omega = eps_e - eps_o.
    response      chi(i xi; z, z') = int dp/2pi [ sum_{e at p+k} conj N(z) N(z') / (i xi - Omega) - sum_{e at p-k} N(z) conj N(z') / (i xi + Omega) ]
    fluctuations  C(z, z')        = (1/2) int dp/2pi [ sum_{e at p+k} conj N(z) N(z') + sum_{e at p-k} N(z) conj N(z') ]   (equal time, symmetrized)
The sea's own lines (sea -> empty) are the vacuum's polarization; they are left out (the NOTE gives their size).

THE CARRIER. v(k; z, z') = q^2 g_k(z - z'),  g_k(z) = e^{-kappa |z|} / (2 sinh kappa),  cosh kappa = 2 - cos k:
the carrier walk's occupation per transverse momentum (VI.2 (E)), divided by its 4 neighbours (check 2). It is used
static: the near-field end of the end-to-end rule, good while light crosses the gap fast against the plates' own rate,
a d << 1 (a defined in check 4).

ONE PLATE. Screened response R = chi (1 - v chi)^{-1}. What the plate sends back when the carrier's e^{kappa z} arrives
from the right:  r_R = (q^2 / 2 sinh kappa) sum_{z,z'} e^{kappa (z - zc)} R(z, z') e^{kappa (z' - zc)};  r_L with e^{-kappa (z - zc)}.

TWO PLATES, centre to centre d. No walker crosses the gap, so chi = chi_1 (+) chi_2; only the carrier joins them.
    energy  E(d) = int_{-pi}^{pi} dk/2pi (1/2pi) int_0^inf dxi Re[ ln det(1 - v chi)_both - ln det(1 - v chi)_1 - ln det(1 - v chi)_2 ]
    push    F(d) = - < d V_12 / d d >  =  int dk/2pi Re sum_{z1, z2} kappa v_12(z1, z2) C_12(z1, z2)
            C_12 = the gap-induced correlation of the two plates' charges = -(1/pi) int_0^inf dxi [Hermitian part of chi (1 - v chi)^{-1}]_12

Checks:
  1. the plate: a uniform region's gap is min(m, pi/2 - m); U(p + pi) = Pi U(p) Pi and the mirror hold; every column below
     E_F, at every p, is a slab-bound conduction column; nothing delocalized below 0.78
  2. the carrier per transverse momentum: the 2d carrier walk's occupation, Fourier transformed along y, is 4 g_k
  3. the plate's fluctuations are its response, C = -(1/pi) int dxi Herm chi, because every line is upward: the filled sea
     blocks the conduction walkers' downward lines (alone they would be a lamp)
  4. the plate is a metal: its charge response summed over the slab is the 1D electron gas's, -D k^2/(xi^2 + vF^2 k^2),
     D = 2 vF/pi from the walkers' own band; r -> the plasma sheet's -a k/(xi^2 + a k), a = q^2 D/2; the k -> -k mirror
  5. only round trips cross the gap: det(both)/(det 1 det 2) = 1 - r_R r_L e^{-2 kappa d} (assembled plates, all k, xi, d);
     and from one box holding both slabs' walkers
  6. the force two ways: -dE/dd = the push from the gap-induced correlation; C_12 falls as e^{-kappa d} (one crossing),
     the force per k as e^{-2 kappa d} (a round trip); the first round trip's share
  7. against Lifshitz's near-field formula, the plates' material read from their walkers (nothing fitted): (a) a 1D electron gas
     spread across the slab with the conduction walkers' own density profile; (b) the same gas as a zero-thickness sheet;
     (c) the plasma sheet; plus the quadrature (grids doubled)
NOTE: the size of the omitted vacuum polarization; what the static carrier leaves out.
"""
import numpy as np, scipy.linalg as sl, scipy.sparse as sps, scipy.sparse.linalg as spla, sys, time
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all = ok_all and bool(cond); print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()
def note(msg): print("NOTE  " + msg); sys.stdout.flush()
T0 = time.time()

# ---------------- the walk ----------------
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); I2 = np.eye(2)
def make_walk(Lz, mz):
    dim = 2*Lz; Sh = np.zeros((dim, dim)); Mm = np.zeros((dim, dim))
    for zz in range(Lz):
        Sh[2*((zz + 1) % Lz), 2*zz] = 1; Sh[2*((zz - 1) % Lz) + 1, 2*zz + 1] = 1
        c, s = np.cos(mz[zz]), np.sin(mz[zz]); Mm[2*zz:2*zz + 2, 2*zz:2*zz + 2] = [[c, -s], [s, c]]
    MS = Mm @ Sh
    def U(p): return MS @ np.kron(np.eye(Lz), np.cos(p)*I2 - 1j*np.sin(p)*sx)
    def states(p, allstates=False):                 # orthonormal columns (Schur of a unitary), ascending quasi-energy
        T, Z = sl.schur(U(p), output='complex'); e = -np.angle(np.diag(T))
        if not allstates:
            pos = e > 0; e = e[pos]; Z = Z[:, pos]
        o = np.argsort(e); return e[o], Z[:, o]
    return U, states

Lz = 40; W = 6; zc = Lz//2; z = np.arange(Lz); slab = (z >= zc - W//2) & (z < zc + W//2); zcen = zc - 0.5
m_in, m_out = 0.30, np.pi/4
pF = 0.05; q2 = 0.02
U, states = make_walk(Lz, np.where(slab, m_in, m_out))
near = np.abs(z - zcen) <= W/2 + 3
EF = states(pF)[0][0]; h = 1e-4
vF = (states(pF + h)[0][0] - states(pF - h)[0][0])/(2*h)
inv_mstar = (states(h)[0][0] - 2*states(0.0)[0][0] + states(-h)[0][0])/h**2
D = 2*vF/np.pi; nlin = 2*pF/np.pi; a_pl = q2*D/2

# ---------------- 1. the plate ----------------
P_ = np.linspace(-np.pi, np.pi, 361); K_ = np.linspace(-np.pi, np.pi, 361); PP, KK = np.meshgrid(P_, K_)
gaps = []
for m in (0.30, np.pi/4, 1.20):
    # the walk's own 2x2 tick for a uniform region, U = M S C (det 1), cos eps = Re tr U / 2
    c, s = np.cos(m), np.sin(m)
    Mm2 = np.array([[c, -s], [s, c]])
    # S = diag(e^{-i kz}, e^{i kz}), C = cos p - i sin p sigma_x; tr(M S C) = sum_ij M_ij (S C)_ji
    S0, S1 = np.exp(-1j*KK), np.exp(1j*KK); cp, sp = np.cos(PP), np.sin(PP)
    C = np.array([[cp, -1j*sp], [-1j*sp, cp]])                       # C[i,j] as arrays
    SC = np.array([[S0*C[0, 0], S0*C[0, 1]], [S1*C[1, 0], S1*C[1, 1]]])
    tr = Mm2[0, 0]*SC[0, 0] + Mm2[0, 1]*SC[1, 0] + Mm2[1, 0]*SC[0, 1] + Mm2[1, 1]*SC[1, 1]
    gaps.append((m, np.arccos(np.clip(tr.real/2, -1, 1)).min(), min(m, np.pi/2 - m)))
Pi = np.kron(np.diag((-1.0)**z), I2)
Zm = np.kron(np.eye(Lz)[::-1], sy)                               # z -> 39 - z, senses swapped by sigma_y
p0 = 0.37
res_pi = np.abs(U(p0 + np.pi) - Pi @ U(p0) @ Pi).max(); res_mir = np.abs(Zm @ U(p0) @ Zm.conj().T - U(-p0)).max()
scan_p = np.linspace(0, np.pi, 181); below = []; lowest_free = 9.0; second = 9.0
for p in scan_p:
    e, V = states(p); wgt = (np.abs(V.reshape(Lz, 2, -1))**2).sum(1)[near].sum(0)
    for j in np.where(e < EF)[0]: below.append((p, e[j], wgt[j], j))
    free = wgt < 0.5
    if free.any(): lowest_free = min(lowest_free, e[free].min())
    bound = np.where(wgt >= 0.5)[0]
    if len(bound) > 1: second = min(second, e[bound[1]])
nb_ok = all(j == 0 and w > 0.85 and (pp < pF + 1e-9 or pp > np.pi - pF - 1e-9) for pp, ee, w, j in below)
ok("1. the plate: a uniform region's gap from the walk's own tick is %s (= min(m, pi/2 - m)); U(p+pi) = Pi U(p) Pi to %.0e and the "
   "mirror (z -> 39 - z, senses by sigma_y) takes U(p) to U(-p) to %.0e; scanning p over [0, pi]: the %d columns below E_F = %.5f "
   "are all the lowest slab-bound column (slab weight >= %.3f), at |p| < p_F or |p - pi| < p_F; nothing delocalized below %.4f, "
   "next bound band from %.4f.  vF = %.5f, 1/m* = %.4f, walkers per unit length %.5f (both pockets), D = 2 vF/pi = %.5f"
   % (", ".join("%.4f at m = %.3f" % (g, m) for m, g, _ in gaps), res_pi, res_mir, len(below), EF,
      min(w for _, _, w, _ in below), lowest_free, second, vF, inv_mstar, nlin, D),
   all(abs(g - t) < 2e-3 for _, g, t in gaps) and res_pi < 1e-13 and res_mir < 1e-13 and nb_ok and lowest_free > EF + 0.3 and second > EF + 0.2)

# ---------------- 2. the carrier ----------------
Ly = 256; Lzc = 799; src_z = 399                                  # z unknowns, absorbing 400 sites from the source on both sides
ey = np.ones(Ly); Ky = sps.diags([-ey[:-1], 2*ey, -ey[:-1]], [-1, 0, 1], format='lil'); Ky[0, Ly - 1] = -1; Ky[Ly - 1, 0] = -1
ez = np.ones(Lzc); Kz = sps.diags([-ez[:-1], 2*ez, -ez[:-1]], [-1, 0, 1])
Kc = (sps.kron(Ky.tocsr(), sps.identity(Lzc)) + sps.kron(sps.identity(Ly), Kz)).tocsc()      # 4 - A
rhs = np.zeros(Ly*Lzc); rhs[0*Lzc + src_z] = 4.0
Gc = spla.spsolve(Kc, rhs).reshape(Ly, Lzc)
Gk = np.fft.fft(Gc, axis=0)
def kappa_of(k): return np.arccosh(2 - np.cos(k))
def g_of(k, dz): ka = kappa_of(k); return np.exp(-ka*np.abs(dz))/(2*np.sinh(ka))
worst2 = 0.0; rows2 = []
for n in (2, 8, 32):
    k = 2*np.pi*n/Ly; zz = np.arange(src_z - 30, src_z + 31)
    num = Gk[n, zz].real; ex = 4*g_of(k, zz - src_z); keep = ex > 1e-6*ex.max()       # below that the mode is at roundoff
    worst2 = max(worst2, np.abs(num[keep]/ex[keep] - 1).max()); rows2.append((k, kappa_of(k)))
ok("2. the carrier: the 2d carrier walk's occupation ((4 - A) G = 4 delta, 256 x 799 sites, absorbed 400 sites out), Fourier transformed along y, equals "
   "4 e^{-kappa|z|}/(2 sinh kappa) within 30 sites of the source (down to 1e-6 of its peak) to %.1e at k = %s; kappa = %s (-> k: the 2d Coulomb kernel e^{-k|z|}/2k)"
   % (worst2, ", ".join("%.4f" % k for k, _ in rows2), ", ".join("%.5f" % ka for _, ka in rows2)), worst2 < 1e-6)

# ---------------- lines and the response ----------------
gl_x, gl_w = np.polynomial.legendre.leggauss(12)
def pocket_nodes(k):
    br = sorted(set([-pF, pF] + [b for b in (pF - k, -pF + k) if -pF < b < pF]))
    P, Wt = [], []
    for lo, hi in zip(br[:-1], br[1:]):
        P += list(0.5*(hi - lo)*gl_x + 0.5*(hi + lo)); Wt += list(0.5*(hi - lo)*gl_w)
    return np.array(P), np.array(Wt)
def lines(states_fn, LL, k, npockets=2, center=0.0, with_sea=False):
    """lines from the occupied conduction columns at p to the empty columns at p + k (+1) and p - k (-1)"""
    P, Wt = pocket_nodes(k)
    out = {+1: [[], [], []], -1: [[], [], []]}; sea = {+1: [], -1: []}
    for p, w in zip(P + center, Wt):
        e0, V0 = states_fn(p); occ = np.where(e0 < EF)[0]
        for s in (+1, -1):
            if with_sea:
                e1, V1 = states_fn(p + s*k, allstates=True)
            else:
                e1, V1 = states_fn(p + s*k)
            emp = e1 > EF
            for j in occ:
                o = V0[:, j].reshape(LL, 2)
                E = V1[:, emp].reshape(LL, 2, -1); N = np.einsum('zcj,zc->jz', E.conj(), o)
                Om = e1[emp] - e0[j]; assert (Om > 0).all()
                out[s][0].append(np.full(len(Om), w*npockets/(2*np.pi))); out[s][1].append(Om); out[s][2].append(N)
                if with_sea:
                    sm = e1 < 0; Ns = np.einsum('zcj,zc->jz', V1[:, sm].reshape(LL, 2, -1).conj(), o)
                    sea[s].append((w*npockets/(2*np.pi), e1[sm] - e0[j], Ns))
    pack = {s: (np.concatenate(out[s][0]), np.concatenate(out[s][1]), np.concatenate(out[s][2])) for s in (+1, -1)}
    return (pack, sea) if with_sea else pack
def chi_of(pack, xis, flip=False):
    a, b = (+1, -1) if not flip else (-1, +1)                         # flip: the response at -k
    w, Om, N = pack[a]; A = (w[:, None, None]*N.conj()[:, :, None]*N[:, None, :]).reshape(len(w), -1)
    w2, Om2, N2 = pack[b]; B = (w2[:, None, None]*N2[:, :, None]*N2.conj()[:, None, :]).reshape(len(w2), -1)
    chi = (1.0/(1j*xis[:, None] - Om[None, :])) @ A - (1.0/(1j*xis[:, None] + Om2[None, :])) @ B
    LL = N.shape[1]; return chi.reshape(len(xis), LL, LL)
def corr_of(pack):
    w, Om, N = pack[+1]; w2, Om2, N2 = pack[-1]
    return 0.5*(np.einsum('t,tz,tw->zw', w, N.conj(), N) + np.einsum('t,tz,tw->zw', w2, N2, N2.conj()))
def logw(x):                                                       # trapezoid weights for int f dx on a log grid
    u = np.log(x); du = np.diff(u); wu = np.zeros_like(x); wu[:-1] += du/2; wu[1:] += du/2; return wu*x
def herm(X): return 0.5*(X + np.conj(np.swapaxes(X, -1, -2)))

# ---------------- 3. fluctuations = response ----------------
k3 = 0.02
pk3, sea3 = lines(states, Lz, k3, with_sea=True)
C_cols = corr_of(pk3)
xg = np.exp(np.linspace(np.log(1e-12), np.log(1e8), 500)); wx = logw(xg)
chi3 = chi_of(pk3, xg)
C_resp = -(1/np.pi)*np.tensordot(wx, herm(chi3), axes=(0, 0))
fdt_err = np.abs(C_resp - C_cols).max()/np.abs(C_cols).max()
# the lines the filled sea blocks: from the occupied conduction columns down to the sea's columns, read at a slab site
zs = zc - 1
up = sum((pk3[s][0]*np.abs(pk3[s][2][:, zs])**2).sum() for s in (+1, -1))
down = sum((w*np.abs(Ns[:, zs])**2).sum() for s in (+1, -1) for w, Om_, Ns in sea3[s])
down_rate = np.concatenate([Om_ for s in (+1, -1) for _, Om_, _ in sea3[s]])
nlines = sum(len(pk3[s][1]) for s in (+1, -1)); minOm = min(pk3[s][1].min() for s in (+1, -1))
ok("3. the plate's fluctuations are its response: at k = %.2f the equal-time charge correlation from the columns equals "
   "-(1/pi) int_0^inf dxi Herm chi(i xi) to %.1e; all %d lines are upward (smallest rate %.2e). The filled sea blocks the "
   "downward lines: at a slab site they would carry %.1f%% of a conduction walker's line weight, at rates %.3f to %.3f "
   "(a lone conduction walker is a lamp; in the lowest state it is not)"
   % (k3, fdt_err, nlines, minOm, 100*down/(up + down), -np.abs(down_rate).max(), -np.abs(down_rate).min()),
   fdt_err < 1e-6 and minOm > 0 and (down_rate < 0).all())

# ---------------- the main run: per k, the plate's reflection, and the assembled pair ----------------
NK, NX = 96, 128
ks = np.exp(np.linspace(np.log(1e-6), np.log(0.6), NK)); wk = logw(ks)
xis = np.exp(np.linspace(np.log(1e-9), np.log(20.0), NX)); wxi = logw(xis)
dlist = [24, 32, 48, 64, 96, 128, 172]
win = np.abs(z - zcen) <= 11.5; nw = win.sum(); zw = z[win]
rR = np.zeros((NK, NX), complex); rL = np.zeros_like(rR); rRw = np.zeros_like(rR); rLw = np.zeros_like(rR)
pair_ratio = {d: np.zeros((NK, NX), complex) for d in dlist}
FHF_k = {d: np.zeros(NK) for d in dlist}; C12_k = {d: np.zeros(NK) for d in dlist}
chk4 = {}
t1 = time.time()
for ik, k in enumerate(ks):
    pk = lines(states, Lz, k)
    chi = chi_of(pk, xis)
    ka = kappa_of(k); cc = q2/(2*np.sinh(ka))
    V11 = q2*g_of(k, z[:, None] - z[None, :])
    ep = np.exp(ka*(z - zcen)); em = np.exp(-ka*(z - zcen))
    R = np.linalg.solve(np.eye(Lz)[None] - chi @ V11[None], chi)
    rR[ik] = cc*np.einsum('z,xzw,w->x', ep, R, ep); rL[ik] = cc*np.einsum('z,xzw,w->x', em, R, em)
    cw = chi[:, win][:, :, win]; Vw = V11[np.ix_(win, win)]
    Rw = np.linalg.solve(np.eye(nw)[None] - cw @ Vw[None], cw)
    rRw[ik] = cc*np.einsum('z,xzw,w->x', ep[win], Rw, ep[win]); rLw[ik] = cc*np.einsum('z,xzw,w->x', em[win], Rw, em[win])
    s1, l1 = np.linalg.slogdet(np.eye(nw)[None] - Vw[None] @ cw)
    for d in dlist:
        Zall = np.concatenate([zw, zw + d]); Vf = q2*g_of(k, Zall[:, None] - Zall[None, :])
        cf = np.zeros((NX, 2*nw, 2*nw), complex); cf[:, :nw, :nw] = cw; cf[:, nw:, nw:] = cw
        sf, lf = np.linalg.slogdet(np.eye(2*nw)[None] - Vf[None] @ cf)
        pair_ratio[d][ik] = (sf/s1**2)*np.exp(lf - 2*l1)
        Rf = np.linalg.solve(np.eye(2*nw)[None] - cf @ Vf[None], cf)
        C12 = -(1/np.pi)*np.tensordot(wxi, herm(Rf)[:, :nw, nw:], axes=(0, 0))
        FHF_k[d][ik] = np.real((ka*Vf[:nw, nw:]*C12).sum()); C12_k[d][ik] = np.real(C12.sum())
    if ik in (np.argmin(abs(ks - 0.002)),):
        chk4['k'] = k; chk4['pk'] = pk; chk4['chi'] = chi; chk4['ik'] = ik
print("      main run: %d transverse momenta x %d imaginary frequencies, %.0f s" % (NK, NX, time.time() - t1)); sys.stdout.flush()

# ---------------- 4. the plate is a metal ----------------
k4 = chk4['k']; ik4 = chk4['ik']; pk4 = chk4['pk']
x4 = np.array([0.3, 1.0, 3.0])*k4*vF
chi4 = chi_of(pk4, x4); sheet = chi4.sum(axis=(1, 2)); gas = -D*k4**2/(x4**2 + vF**2*k4**2)
dev_sheet = np.abs(sheet/gas - 1).max()
xr = np.array([0.3, 1.0, 3.0])*np.sqrt(a_pl*k4)
chir = chi_of(pk4, xr); ka4 = kappa_of(k4); cc4 = q2/(2*np.sinh(ka4)); V4 = q2*g_of(k4, z[:, None] - z[None, :])
R4 = np.linalg.solve(np.eye(Lz)[None] - chir @ V4[None], chir); ep4 = np.exp(ka4*(z - zcen)); em4 = np.exp(-ka4*(z - zcen))
r4R = cc4*np.einsum('z,xzw,w->x', ep4, R4, ep4); r4L = cc4*np.einsum('z,xzw,w->x', em4, R4, em4)
r_gas4 = -a_pl*k4/(xr**2 + vF**2*k4**2 + a_pl*k4); r_pl4 = -a_pl*k4/(xr**2 + a_pl*k4)
chim = chi_of(pk4, xr, flip=True); Pz = np.eye(Lz)[::-1]
mir_err = np.abs(chim - Pz @ chir @ Pz).max()/np.abs(chir).max()
pk_pi = lines(states, Lz, k4, npockets=1, center=np.pi); pk_0 = lines(states, Lz, k4, npockets=1)
rep_err = np.abs(chi_of(pk_pi, xr) - chi_of(pk_0, xr)).max()/np.abs(chi_of(pk_0, xr)).max()
KA = kappa_of(ks)[:, None]
def force(rr, rl, d, first=False):
    x = rr*rl*np.exp(-2*KA*d); f = x if first else x/(1 - x)
    return -(1/np.pi)*wk @ ((1/(2*np.pi))*((2*KA*f).real @ wxi))
lowk = ks < 0.15
trunc_r = max(np.abs(rRw - rR)[lowk].max(), np.abs(rLw - rL)[lowk].max())
trunc_F = max(abs(force(rRw, rLw, d)/force(rR, rL, d) - 1) for d in dlist)
ok("4. the plate is a metal: at k = %.4f its charge response summed over the slab equals the 1D electron gas's -D k^2/(xi^2 + vF^2 k^2) "
   "to %.1e (xi = 0.3, 1, 3 k vF; D from the walkers' band). Its reflection at xi = 0.3, 1, 3 sqrt(a k), a = q^2 D/2 = %.3e: "
   "r_R = %s, r_L = %s; 1D-gas sheet %s; plasma sheet -a k/(xi^2 + a k) %s. The pocket at p ~ pi repeats p ~ 0 to %.0e; "
   "the response at -k is the mirror image of that at k to %.0e (so r_R(-k) = r_L(k) and the two signs of k give the same round trip); "
   "keeping 24 sites of each plate changes r by %.1e at k < 0.15 and the force by %.0e at every d"
   % (k4, dev_sheet, a_pl, np.round(r4R.real, 4).tolist(), np.round(r4L.real, 4).tolist(), np.round(r_gas4, 4).tolist(),
      np.round(r_pl4, 4).tolist(), rep_err, mir_err, trunc_r, trunc_F),
   dev_sheet < 1e-3 and np.abs(r4R.real - r_gas4).max() < 0.01 and np.abs(r4L.real - r_gas4).max() < 0.01 and rep_err < 1e-10 and mir_err < 1e-10
   and trunc_F < 1e-5)

# ---------------- 5. only round trips cross the gap ----------------
dev5 = 0.0
for d in dlist:
    x = rRw*rLw*np.exp(-2*kappa_of(ks)[:, None]*d)
    dev5 = max(dev5, np.abs(pair_ratio[d] - (1 - x)).max())
# one box holding both slabs' walkers
dB = 40; LB = Lz + dB; zB = np.arange(LB)
slabB = ((zB >= zc - W//2) & (zB < zc + W//2)) | ((zB >= zc - W//2 + dB) & (zB < zc + W//2 + dB))
UB, statesB = make_walk(LB, np.where(slabB, m_in, m_out))
winB = np.concatenate([np.where(win)[0], np.where(win)[0] + dB])
rows5 = []; blk = 0.0
for k in (0.01, 0.03):
    pkB = lines(statesB, LB, k); xB = np.array([1e-3, 3e-3, 1e-2])
    chB = chi_of(pkB, xB)[:, winB][:, :, winB]
    blk = max(blk, np.abs(chB[:, :nw, nw:]).max()/np.abs(chB).max())
    pk1 = lines(states, Lz, k); ch1 = chi_of(pk1, xB)[:, win][:, :, win]
    ka = kappa_of(k); cc = q2/(2*np.sinh(ka)); Zall = np.concatenate([zw, zw + dB]); Vf = q2*g_of(k, Zall[:, None] - Zall[None, :]); Vw = Vf[:nw, :nw]
    Rw1 = np.linalg.solve(np.eye(nw)[None] - ch1 @ Vw[None], ch1)
    rr = cc*np.einsum('z,xzw,w->x', np.exp(ka*(zw - zcen)), Rw1, np.exp(ka*(zw - zcen))); rl = cc*np.einsum('z,xzw,w->x', np.exp(-ka*(zw - zcen)), Rw1, np.exp(-ka*(zw - zcen)))
    sB, lB = np.linalg.slogdet(np.eye(2*nw)[None] - Vf[None] @ chB); s1, l1 = np.linalg.slogdet(np.eye(nw)[None] - Vw[None] @ ch1)
    ratio = (sB/s1**2)*np.exp(lB - 2*l1); pred = 1 - rr*rl*np.exp(-2*ka*dB)
    rows5.append((k, np.abs(ratio - pred).max()/np.abs(1 - pred).max(), np.abs(1 - pred).max()))
ok("5. only round trips cross the gap: assembled plates, every k, xi and d = %s: det(1 - v chi)_both / (det_1 det_2) = 1 - r_R r_L e^{-2 kappa d} "
   "to %.1e. One box holding both slabs' walkers (d = %d, %d sites): the walkers' response has no cross-gap part (%.0e of its size), and "
   "the same ratio holds to %s of the round-trip term (%s)"
   % (dlist, dev5, dB, LB, blk, ", ".join("%.0e" % e for _, e, _ in rows5), ", ".join("k = %.2f: 1 - ratio ~ %.1e" % (k, m) for k, _, m in rows5)),
   dev5 < 1e-10 and blk < 1e-6 and all(e < 1e-3 for _, e, _ in rows5))

# ---------------- 6. the force two ways ----------------
rows6 = []
for d in dlist:
    Fe = force(rRw, rLw, d); Fh = (1/np.pi)*wk @ FHF_k[d]; F1 = force(rR, rL, d, first=True); Ff = force(rR, rL, d)
    rows6.append((d, Fe, Fh, F1/Ff))
dev6 = max(abs(Fh/Fe - 1) for _, Fe, Fh, _ in rows6)
ikc = np.argmin(abs(ks - 0.05)); kc = ks[ikc]
rate = -np.log(C12_k[96][ikc]/C12_k[48][ikc])/48
sign_ok = all(C12_k[d][ikc] < 0 for d in dlist) and all(Fe < 0 for _, Fe, _, _ in rows6)
ok("6. the force two ways, d = %s: the push from the gap-induced correlation, -<dV_12/dd>, equals -dE/dd to %.1e (F = %s). "
   "The correlation is negative (opposite charges face each other) and at k = %.4f falls at %.5f per site from d = 48 to 96, against kappa = %.5f "
   "(one crossing); the force per k falls as e^{-2 kappa d} (a round trip). First round trip's share of the force: %s"
   % (dlist, dev6, ", ".join("%.3e" % Fe for _, Fe, _, _ in rows6), kc, rate, kappa_of(kc), ", ".join("%.3f" % s for *_, s in rows6)),
   dev6 < 1e-8 and sign_ok and abs(rate/kappa_of(kc) - 1) < 5e-3)

# ---------------- 7. against Lifshitz (near field) ----------------
def lifshitz_F(rfun, d):
    K = ks[:, None]; X = xis[None, :]; r = rfun(K, X); x = r*r*np.exp(-2*K*d)
    return -(1/np.pi)*wk @ ((1/(2*np.pi))*((2*K*x/(1 - x)) @ wxi))
r_gas = lambda K, X: -a_pl*K/(X**2 + vF**2*K**2 + a_pl*K)
r_pl = lambda K, X: -a_pl*K/(X**2 + a_pl*K)
# the same 1D gas, spread across the slab with the conduction walkers' own density profile n(z) at the two Fermi points
nprof = sum((np.abs(states(s*pF)[1][:, 0].reshape(Lz, 2))**2).sum(1) for s in (+1, -1)); nprof /= nprof.sum()
def r_thick(kk, xx):
    kak = kappa_of(kk); cck = q2/(2*np.sinh(kak))
    Fself = np.array([(nprof[:, None]*nprof[None, :]*np.exp(-ka*np.abs(z[:, None] - z[None, :]))).sum() for ka in kak])
    GR = np.array([(nprof*np.exp(ka*(z - zcen))).sum() for ka in kak]); GL = np.array([(nprof*np.exp(-ka*(z - zcen))).sum() for ka in kak])
    chis = -D*kk[:, None]**2/(xx[None, :]**2 + vF**2*kk[:, None]**2)
    den = 1 - cck[:, None]*Fself[:, None]*chis
    return cck[:, None]*GR[:, None]**2*chis/den, cck[:, None]*GL[:, None]**2*chis/den
def force_on(rr, rl, d, kk, xx):
    K_ = kappa_of(kk)[:, None]; x = rr*rl*np.exp(-2*K_*d)
    return -(1/np.pi)*logw(kk) @ ((1/(2*np.pi))*((2*K_*x/(1 - x)).real @ logw(xx)))
rth_R, rth_L = r_thick(ks, xis)
ks2 = np.exp(np.linspace(np.log(1e-6), np.log(0.6), 2*NK - 1)); xis2 = np.exp(np.linspace(np.log(1e-9), np.log(20.0), 2*NX - 1))
r2R, r2L = r_thick(ks2, xis2)
grid_dev = max(abs(force_on(rth_R, rth_L, d, ks, xis)/force_on(r2R, r2L, d, ks2, xis2) - 1) for d in dlist)
ell = (nprof[:, None]*nprof[None, :]*np.abs(z[:, None] - z[None, :])).sum()
# plasma sheet in closed form: E = -c_pl sqrt(a) d^{-3/2}, c_pl = -(1/2pi^2) int ds int dt ln(1 - s^2 e^{-2s}/(t^2+s)^2)
sg = np.exp(np.linspace(np.log(1e-6), np.log(40), 400)); tg = np.exp(np.linspace(np.log(1e-6), np.log(1e3), 400))
S_, T_ = np.meshgrid(sg, tg, indexing='ij'); c_pl = -(1/(2*np.pi**2))*logw(sg) @ (np.log(1 - S_**2*np.exp(-2*S_)/(T_**2 + S_)**2) @ logw(tg))
rows7 = []
for d in dlist:
    Fw = force(rR, rL, d); Ft = force(rth_R, rth_L, d); Fg = lifshitz_F(r_gas, d); Fp = lifshitz_F(r_pl, d); Fc = -1.5*c_pl*np.sqrt(a_pl)*d**-2.5
    rows7.append((d, a_pl*d, vF**2/(a_pl*d), Fw, Fw/Ft, Fw/Fg, Fw/Fp, Fp/Fc))
print("      d    a d (light)  vF^2/(a d)   F_walkers    F/F_thick   F/F_sheet   F/F_plasma   (plasma: grid/closed form)")
for d, ad, fm, Fw, rt, rg, rp, pc in rows7:
    print("    %4d    %7.4f     %7.4f    %+.4e    %.4f      %.4f      %.4f       %.6f" % (d, ad, fm, Fw, rt, rg, rp, pc))
th_dev = [abs(r[4] - 1) for r in rows7]; gas_dev = [r[5] - 1 for r in rows7]; pl = [r[6] for r in rows7]
ok("7. against Lifshitz's near-field formula, with the plates' material read from their walkers (n(z), vF, D; nothing fitted): "
   "(a) a 1D electron gas spread across the slab with the conduction walkers' own density profile: agreement to %.1f%% at d = %d, "
   "%.2f%% at d = %d, %.3f%% at d >= %d; (b) the same gas as a zero-thickness sheet: %.1f%% at d = %d falling to %.1f%% at d = %d, "
   "the slab's thickness (<|z - z'|> = %.2f sites over the walkers' profile) against d; (c) the plasma sheet: the ratio climbs "
   "%.3f -> %.3f as vF^2/(a d) falls %.2f -> %.3f (the walkers' own motion, which the plasma model leaves out). The plasma-sheet "
   "grid integral matches its closed form -(3/2) c sqrt(a) d^{-5/2} (c = %.5f) to %.0e; doubling both grids moves the force by %.0e"
   % (100*th_dev[0], dlist[0], 100*th_dev[2], dlist[2], 100*max(th_dev[4:]), dlist[4], 100*gas_dev[0], dlist[0], 100*gas_dev[-1], dlist[-1],
      ell, pl[0], pl[-1], rows7[0][2], rows7[-1][2], c_pl, max(abs(r[-1] - 1) for r in rows7), grid_dev),
   th_dev[0] < 0.02 and max(th_dev[4:]) < 1e-3 and all(gas_dev[i + 1] < gas_dev[i] for i in range(len(gas_dev) - 1)) and gas_dev[-1] > 0
   and all(pl[i + 1] > pl[i] for i in range(len(pl) - 1)) and max(abs(r[-1] - 1) for r in rows7) < 1e-3 and grid_dev < 1e-4)

# ---------------- NOTE: the vacuum's own polarization ----------------
def alpha_sea(m, qy=0.01, n=240):
    pp = (np.arange(n) + 0.5)*2*np.pi/n - np.pi; P, K = np.meshgrid(pp, pp, indexing='ij')
    c, s = np.cos(m), np.sin(m); Mm2 = np.array([[c, -s], [s, c]])
    def bloch(P):
        S0, S1 = np.exp(-1j*K), np.exp(1j*K); cp, sp = np.cos(P), np.sin(P)
        Uk = np.zeros(P.shape + (2, 2), complex)
        C = np.array([[cp, -1j*sp], [-1j*sp, cp]]).transpose(2, 3, 0, 1)
        Sd = np.stack([S0, S1], -1)[..., :, None]*C
        Uk = np.einsum('ij,...jk->...ik', Mm2, Sd)
        lam, V = np.linalg.eig(Uk); e = -np.angle(lam); o = np.argsort(e, axis=-1)
        e = np.take_along_axis(e, o, -1); V = np.take_along_axis(V, o[..., None, :], -1)
        return e, V/np.linalg.norm(V, axis=-2, keepdims=True)
    e0, V0 = bloch(P); tot = 0.0
    for sgn in (+1, -1):
        e1, V1 = bloch(P + sgn*qy)
        ov = np.abs(np.einsum('...c,...c->...', V1[..., :, 1].conj(), V0[..., :, 0]))**2
        tot += (ov/(e1[..., 1] - e0[..., 0])).mean()
    return tot/qy**2                                                 # alpha per unit area: chi = -alpha q^2
al_out, al_in = alpha_sea(m_out), alpha_sea(m_in)
note("the sea's own lines were left out: they are the vacuum's polarization. For a uniform region the sea's static response is "
     "-alpha k^2 with alpha = %.4f (m = pi/4, outside) and %.4f (m = 0.30, slab): a dielectric constant 1 + q^2 alpha = 1 + %.1e outside, "
     "1 + %.1e in the slab. The outside value rescales q^2 everywhere (a renormalization of the charge: q^2 here is the screened one); "
     "the difference is a faint dielectric background in the slab. Either moves the force by under 1%%."
     % (al_out, al_in, q2*al_out, q2*al_in))
note("the carrier was used static. The moving end-to-end rule (the Fokker kernel) adds the gap's light-crossing time and the "
     "current-current (magnetic) term; both enter at order a d = %.3f-%.3f for these plates (column 2 of the table). The far-field "
     "(Casimir) end needs that rule; the near field is complete here." % (rows7[0][1], rows7[-1][1]))
print("      total %.0f s" % (time.time() - T0))
print("ALL PASS" if ok_all else "SOME FAILED"); raise SystemExit(0 if ok_all else 1)
