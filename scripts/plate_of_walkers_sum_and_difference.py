"""
A plate of walkers with exclusion: do its fluctuations (the 'sum') equal its response (the 'difference')?

The walk: a 1+1 ring of L sites, the real tick U = R(theta_x) S of III.15a (stream + right, - left, then the coin at the
arrival site). The plate: a slab where the coin angle (the mass, V.1's knob) is small, theta_in, inside a region where it is
large, theta_out; low-energy walkers are confined to the slab. Walkers with exclusion: a set of occupied columns z_a
(eigencolumns of the tick, quasi-energies eps_a in (-pi, pi]); the negative-energy states are the sea.

A reading J (the slab's total current, a site's density, a site's current) at two ticks, with J(t) = U_t^T J U_t:
    C(t,s) = sum_a z_a^dag J(t) J(s) z_a  -  sum_{a,a'} (z_a^dag J(t) z_a')(z_a'^dag J(s) z_a)
             [each walker's own marked out-and-return pairs]   [exchange pairs: two walkers swapping between the marks, VI.1's E]
           = Tr[P J(t) (1-P) J(s)]                                                        (Wick, P = occupied projector)
Its lines: Omega = eps_b - eps_a, weight W = n_a (1 - n_b) |J_ba|^2. At Omega > 0:
    sum (lamp + negative lamp)   = (W(Omega) + W(-Omega)) / 2      <- Re C: the fluctuations
    difference (absorbs - emits) = (W(Omega) - W(-Omega)) / 2      <- Im C: the response (Kubo)
Zero-temperature fluctuation-dissipation: sum = difference at every Omega > 0, i.e. no line with Omega < 0 (nothing steps down).

Checks:
  1. the slab and the sea: low-energy walkers are confined to the slab; the sea is the conjugate half of the walkers (real tick):
     P_+ = (1 + iK)/2 with K real, K^2 = -1 - the lattice form of III.8c's quarter turn J = L Omega^-1
  2. from the columns: C(t,s) computed by ticking each walker's column out and back with two marks, minus the exchange pairs,
     equals Tr[P J(t)(1-P)J(s)] and the line sum; its real part is the closed pairs within each copy, its imaginary part the
     cross term between the copies (both less the exchange pairs)
  3. the plate with its lowest states filled (sea + N walkers): sum = difference at every frequency, for all three readings -
     and why: each walker's own pairs DO contain steps down (each alone would shine); the exchange pairs carry exactly the same
     downward weight and cancel it, because every lower level is occupied
  4. an excited plate (one walker lifted two levels): sum exceeds difference by exactly its emission lines - the lamp -
     the lifted walker falling into the holes; which gaps a reading sees follows the levels' parity
  5. a hot plate (Fermi occupations at temperature T): sum / difference = coth(Omega / 2T) = 1 + 2 nbar(Omega), line by line -
     the excess over the ground-state '1' is Planck's occupation
  6. a lone walker with the sea empty is a lamp: most of its current fluctuation is stepping down into the empty sea, at about
     twice the mass frequency (the beat of III.8c (14)); filling the sea switches that lamp off
"""
import numpy as np, scipy.linalg as sl
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all = ok_all and bool(cond); print(("PASS  " if cond else "FAIL  ") + msg)

L = 240; xs = np.arange(L); slab = (xs >= 90) & (xs < 150)
th_in, th_out = 0.15, 0.60
theta = np.where(slab, th_in, th_out)
dim = 2*L
U = np.zeros((dim, dim))
for x in range(L):
    for s in (0, 1):
        y = (x + 1) % L if s == 0 else (x - 1) % L
        c, sn = np.cos(theta[y]), np.sin(theta[y])
        R = np.array([[c, sn], [-sn, c]])
        U[2*y, 2*x + s] = R[0, s]; U[2*y + 1, 2*x + s] = R[1, s]
assert np.allclose(U.T @ U, np.eye(dim))
Tm, Z = sl.schur(U.astype(complex), output='complex')
lam = np.diag(Tm); nonnormal = np.abs(Tm - np.diag(lam)).max()
eps = -np.angle(lam); order = np.argsort(eps); eps = eps[order]; Z = Z[:, order]
pos = eps > 0; neg = eps < 0
cond = np.where(pos)[0]                                   # positive states, ascending
w_slab = (np.abs(Z.reshape(L, 2, dim))**2).sum(1)[slab].sum(0)

# ---------- 1. slab and sea ----------
nb = 0
while w_slab[cond[nb]] > 0.95: nb += 1
Pp = Z[:, pos] @ Z[:, pos].conj().T
K = -1j*(2*Pp - np.eye(dim))
okK = np.abs(K.imag).max() < 1e-10 and np.abs(K.real @ K.real + np.eye(dim)).max() < 1e-9 and np.abs(Pp.real - 0.5*np.eye(dim)).max() < 1e-10
ok("1. the slab: the lowest %d walker levels have > 95%% of their weight inside it (quasi-energies %.3f .. %.3f, between theta_in = %.2f and theta_out = %.2f); %d positive and %d negative states (the sea), mirror images (max |eps_k + eps_-k| = %.1e); P_+ = (1 + iK)/2 with K real (%.1e), K^2 = -1 (%.1e): the sea is the conjugate half, split by the lattice form of III.8c's quarter turn J"
   % (nb, eps[cond[0]], eps[cond[nb - 1]], th_in, th_out, pos.sum(), neg.sum(), np.abs(np.sort(eps) + np.sort(eps)[::-1]).max(), np.abs(K.imag).max(), np.abs(K.real @ K.real + np.eye(dim)).max()),
   nb >= 8 and okK and nonnormal < 1e-10)

# readings
Jcur = np.diag(np.repeat(slab.astype(float), 2)*np.tile([1.0, -1.0], L))        # the slab's total current (sum over its sites), q c = 1
Jrho = np.zeros((dim, dim)); Jloc = np.zeros((dim, dim))
for x0 in (119, 120):                                                          # the two central sites: the slab's mirror axis is at 119.5
    Jrho[2*x0, 2*x0] = Jrho[2*x0 + 1, 2*x0 + 1] = 1.0                          # density at the centre (even under the mirror)
    Jloc[2*x0, 2*x0] = 1.0; Jloc[2*x0 + 1, 2*x0 + 1] = -1.0                    # current at the centre (odd)
# the mirror: x -> L-1-x with the turn T on the senses (T commutes with the coin and swaps the senses); check it is a symmetry
Mir = np.zeros((dim, dim)); Tt = np.array([[0, 1], [-1, 0.]])
for x in range(L): Mir[2*(L - 1 - x):2*(L - 1 - x) + 2, 2*x:2*x + 2] = Tt
assert np.allclose(Mir @ U, U @ Mir)
readings = {"slab current": Jcur, "density at the centre": Jrho, "current at the centre": Jloc}
Jab = {k: Z.conj().T @ J @ Z for k, J in readings.items()}
Om = eps[:, None] - eps[None, :]                                               # Om[b, a] = eps_b - eps_a
def lines(n, M):                                                               # W[b, a] = n_a (1 - n_b) |J_ba|^2
    return np.outer(1 - n, n)*np.abs(M)**2
def split(W, tol=1e-12):
    return W[Om > tol].sum(), W[Om < -tol].sum()                               # up (absorb), down (emit)

N = 4
n_ground = neg.astype(float); n_ground[cond[:N]] = 1.0

# ---------- 2. from the columns ----------
occ = np.where(n_ground > 0.5)[0]; P = Z[:, occ] @ Z[:, occ].conj().T
W0 = lines(n_ground, Jab["slab current"])
worst = 0.0; Ut = np.eye(dim); parts = None
for tau in range(13):
    Jt = Ut.T @ Jcur @ Ut                                                       # out tau ticks, mark, back (U^T = U^-1)
    own = np.trace(P @ Jt @ Jcur)                                               # sum_a z_a^dag J(tau) J(0) z_a
    exch = np.trace(P @ Jt @ P @ Jcur)                                          # sum_{a,a'} (z_a^dag J(tau) z_a')(z_a'^dag J(0) z_a)
    C_cols = own - exch
    C_wick = np.trace(P @ Jt @ (np.eye(dim) - P) @ Jcur)
    C_lines = (W0*np.exp(-1j*Om*tau)).sum()
    worst = max(worst, abs(C_cols - C_wick), abs(C_cols - C_lines))
    if tau == 5:
        u = Z[:, occ].real; v = Z[:, occ].imag; Mt = Jt @ Jcur
        within = 0.5*(np.einsum('ia,ij,ja->', u, Mt + Mt.T, u) + np.einsum('ia,ij,ja->', v, Mt + Mt.T, v))
        between = np.einsum('ia,ij,ja->', u, Mt - Mt.T, v)
        parts = (within, between, exch, C_cols)
    Ut = U @ Ut
w_, b_, e_, Cc = parts
ok("2. from the columns: each walker ticked out and back with two marks, minus the exchange pairs, = Tr[P J(t)(1-P)J(s)] = the line sum, tau = 0..12: max deviation %.1e. At tau = 5: closed pairs within each copy %.4f, cross term between the copies %.4f, exchange %.4f%+.4fi -> C = %.4f%+.4fi (Re = within - Re exch, Im = between - Im exch)"
   % (worst, w_, b_, e_.real, e_.imag, Cc.real, Cc.imag), worst < 1e-9 and abs(Cc.real - (w_ - e_.real)) < 1e-9 and abs(Cc.imag - (b_ - e_.imag)) < 1e-9)

# ---------- 3. the plate with its lowest states filled ----------
rows = []; mech = []
for k, M in Jab.items():
    up, down = split(lines(n_ground, M))
    self_ = np.outer(np.ones(dim), n_ground)*np.abs(M)**2                     # each walker's own pairs (no exclusion between walkers)
    ex = np.outer(n_ground, n_ground)*np.abs(M)**2                            # exchange pairs
    mech.append((k, split(self_)[1], split(ex)[1]))
    rows.append((k, up, down))
okg = all(d == 0.0 for _, _, d in rows) and all(abs(ds - de) <= 1e-12*max(ds, 1e-300) for _, ds, de in mech)
ok("3. lowest states filled (sea + %d walkers): no line steps down, so sum = difference at every frequency - %s. Why: each walker's own pairs step down with total weight %s, and the exchange pairs carry exactly %s; they cancel, because every lower level is occupied"
   % (N, "; ".join("%s: up %.4f, down %.1e" % r for r in rows), ", ".join("%.4f" % ds for _, ds, _ in mech), ", ".join("%.4f" % de for _, _, de in mech)), okg)
W = lines(n_ground, Jab["slab current"]); up_mask = Om > 1e-12
top = np.argsort(W[up_mask])[::-1][:4]; Omu = Om[up_mask][top]; Wu = W[up_mask][top]
print("      slab current, strongest lines of the filled plate:  " + ";  ".join("Omega = %.4f: sum %.4f, difference %.4f" % (o, w/2, w/2) for o, w in zip(Omu, Wu)))

# ---------- 4. an excited plate ----------
n_exc = n_ground.copy(); n_exc[cond[N - 1]] = 0.0; n_exc[cond[N + 1]] = 1.0        # the top walker lifted two levels
lifted = cond[N + 1]; holes = [cond[N - 1], cond[N]]
pred = {round(eps[lifted] - eps[h], 10) for h in holes}
seen = set(); rep4 = []; ok4 = True; by_reading = {}
for k in ("slab current", "density at the centre"):
    M = Jab[k]; W = lines(n_exc, M); up, down = split(W)
    emis = sorted((round(-Om[b, a], 10), W[b, a]) for b in range(dim) for a in range(dim) if Om[b, a] < -1e-12 and W[b, a] > 1e-12)
    from_lifted = sum(abs(M[h, lifted])**2 for h in holes)
    ok4 = ok4 and down > 1e-6 and abs(down - from_lifted) < 1e-12 and all(o in pred for o, _ in emis)
    seen |= {o for o, _ in emis}; by_reading[k] = {o for o, _ in emis}
    rep4.append("%s: %s" % (k, ", ".join("Omega = %.4f (weight %.4f)" % e for e in emis)))
ok4 = ok4 and seen == pred
ok4 = ok4 and len(by_reading["slab current"]) == 1 and len(by_reading["density at the centre"]) == 1 and by_reading["slab current"] != by_reading["density at the centre"]
ok("4. one walker lifted two levels: sum exceeds difference by its emission lines - %s - all of it the lifted walker falling into the two holes, at the gaps %s; the slab's mirror (x -> L-1-x with T on the senses) is a symmetry of the tick, so the odd current sees only the gap between levels of opposite parity and the even density only the other: the lamp"
   % ("; ".join(rep4), ", ".join("%.4f" % o for o in sorted(pred))), ok4)

# ---------- 5. a hot plate ----------
mu = 0.5*(eps[cond[N - 1]] + eps[cond[N]]); Tk = 0.6*(eps[cond[N]] - eps[cond[N - 1]])
n_hot = 1.0/(np.exp(np.clip((eps - mu)/Tk, -700, 700)) + 1.0)
W = lines(n_hot, Jab["slab current"])
devs = []; tab = []
for b in range(dim):
    for a in range(dim):
        if Om[b, a] > 1e-9 and W[b, a] > 1e-8:
            s_ = W[b, a] + W[a, b]; d_ = W[b, a] - W[a, b]
            ratio = s_/d_; target = 1/np.tanh(Om[b, a]/(2*Tk))
            devs.append(abs(ratio/target - 1)); tab.append((Om[b, a], ratio, target))
tab.sort(key=lambda r: -r[1])
ok("5. a hot plate (T = %.4f): sum / difference = coth(Omega/2T) = 1 + 2 nbar(Omega) on all %d lines, max relative deviation %.1e; e.g. %s - the excess over the cold plate's 1 is Planck's occupation"
   % (Tk, len(devs), max(devs), ", ".join("Omega %.4f: %.4f vs %.4f" % r for r in tab[:3])), max(devs) < 1e-8)

# ---------- 6. a lone walker ----------
n_lone = np.zeros(dim); n_lone[cond[0]] = 1.0
W = lines(n_lone, Jab["slab current"]); up, down = split(W)
dmask = Om < -1e-12
mean_down = (np.abs(Om)*W)[dmask].sum()/W[dmask].sum()
n_sea = neg.astype(float); n_sea[cond[0]] = 1.0
up2, down2 = split(lines(n_sea, Jab["slab current"]))
ok("6. a lone walker, sea empty: %.1f%% of its current fluctuation is stepping down into the empty sea, at a mean |Omega| = %.3f (twice its own quasi-energy is %.3f): it shines at about twice the mass frequency, the beat of III.8c (14). Fill the sea: downward weight %.1e - the lamp is off"
   % (100*down/(up + down), mean_down, 2*eps[cond[0]], down2), 0.92 < down/(up + down) < 0.95 and abs(mean_down/(2*eps[cond[0]]) - 1) < 0.02 and down2 == 0.0)   # 10 Nov 2026: was '> 0.5'; the text quotes 93% at twice the mass frequency
print("ALL PASS" if ok_all else "SOME FAILED"); raise SystemExit(0 if ok_all else 1)
