"""VI.7 (master v39, 10 Nov 2026): leading loop, lagging loop, hole - why the lagging loops' energies add up, from the loop and its reversal.
Parts 1-4 are cited in VI.7 (A)-(F).  Part 5 (the pair made at a vertex) is NOT in the master yet: held for the section on interactions.

Vocabulary and matrices are III.15e's: r = (u+, u-, v+, v-), route u and twin v;  A (advance), S (exchange), T (turn), J (route -> twin);
G(k) = -ck A J + a T the t-reading on a mode;  P+-(k) = (1 +- J G(k)/w)/2 the leading / lagging class;  lambda = (1,0,0,-1), ell = (1,0,0,1).
Energy with its sign is the sense in which the route-twin pair turns: G r = -(E/hbar) J r, i.e. the energy matrix is hbar J G(k).

  1. Two reversals of a loop, and their product.
     K = diag(1,1,-1,-1), the copy swap v -> -v of III.15e (3) and V.2 (H): the twin's quarter turn reversed (v = Tu -> v = -Tu = T^-1 u);
     S with t -> -t, "the process run backward" (the forward and the backward counter exchanged).
     Both carry the leading class onto the lagging class, both reverse a field; K S is the reversal of the motion alone.
     Consequence: the energies come in pairs +-, on a mode and on a ring, and their sum is 0.
  2. The count of configurations (exclusion: a slot holds 0 or 1), built minor by minor on a ring of 4 sites: 256 configurations.
     Counted from the empty cut, a configuration's energy is the sum over the loops present, with no constant.  Because the sum over all
     slots is 0, the lowest configuration (every lagging loop present) has  -(1/2) sum |eps|  and the count is symmetric about 0.
     Relative to it a missing lagging loop and a present leading loop cost the same.
  3. Three dimensions (IV.4b's real matrices, IV.5a's eight counters): the same K also reverses the axis, because the axis, like the
     energy, is a sense of turning read against J.  In IV.4's chiral notation K reads "conjugate, then the turn matrix on the pair of
     senses and the turn matrix on the pair of axis components".  As algebra, a column and its reversal antisymmetrized are VI.3 (0)'s
     T column (3e) - but that is NOT the pair a vertex makes: see 5.
  4. In motion: the settled combination is the equilibrium of a flip process in every frame.  III.15d (12)'s two-rate process (flip rate
     a e^eta running forward, a e^-eta running backward) has rates 0 and -2w; its equilibrium has the imbalance of senses tanh(eta) = v/c
     and a/gamma flips per unit clock; the visible flip rate of a block of length tau is 2 pi+ pi- (1 - e^{-2 w tau})/tau (simulated).
  5. The pair made at a vertex, three dimensions, many walkers (27 sites, 108 slots, a sea of 54).  V.2's turn on one site for one tick, and
     V.1's change of the coin on one site: the pair amplitudes by determinants, then the pair's axis part.  Result: the part proportional
     to T (the column with no net axis) is ZERO for both vertices, exactly; the pair is in the symmetric columns sigma_k T, with the vector
     along the momentum the field supplies (field) or along the relative momentum of the two (coin).  The expectation "the pair made at a
     corner is the T column" is refuted.
"""
import itertools, sys
import numpy as np, sympy as sp, scipy.linalg as sl
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all = ok_all and bool(cond); print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()
mx = lambda M: float(np.abs(np.asarray(M, dtype=complex)).max())

# =====================================================================================================================================
# 1. the two reversals, exact (sympy), III.15e's four matrices
# =====================================================================================================================================
A = sp.diag(1, -1, 1, -1)
S = sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
T = sp.Matrix([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]])
J = sp.Matrix([[0, 0, -1, 0], [0, 0, 0, -1], [1, 0, 0, 0], [0, 1, 0, 0]])
K = sp.diag(1, 1, -1, -1)
I4 = sp.eye(4); Z4 = sp.zeros(4, 4)
k, a, c, phi = sp.symbols('k a c phi', real=True); w = sp.sqrt(a**2 + c**2*k**2)
G = lambda kk: -c*kk*A*J + a*T
Pp = lambda kk: (I4 + J*G(kk)/w)/2
Pm = lambda kk: (I4 - J*G(kk)/w)/2
lam = sp.Matrix([1, 0, 0, -1]); ell = sp.Matrix([1, 0, 0, 1])
zero = lambda M: sp.simplify(M) == sp.zeros(*M.shape)

ok("1.0 the setting: A S = T; J commutes with A, S, T; G(k)^2 = -w^2; the energy matrix J G(k) squares to w^2 and has trace 0; "
   "T lambda = -J lambda (leading), T ell = +J ell (lagging)",
   A*S == T and all(J*M == M*J for M in (A, S, T)) and zero(G(k)*G(k) + w**2*I4) and zero((J*G(k))**2 - w**2*I4) and sp.simplify((J*G(k)).trace()) == 0
   and T*lam == -J*lam and T*ell == J*ell)

ok("1a. K = diag(1,1,-1,-1), the twin's turn reversed: K commutes with A, S, T and anticommutes with J; K lambda = ell; "
   "it turns the launch rule v = Tu into v = -Tu = T^-1 u",
   all(K*M == M*K for M in (A, S, T)) and K*J == -J*K and K*lam == ell
   and all((K*sp.Matrix([u1, u2, u2, -u1]))[2:, 0] == -sp.Matrix([u2, -u1]) for u1, u2 in [sp.symbols('u1 u2')]))

ok("1b. K on a mode: K G(k) K = G(-k), so K P+(k) K = P-(-k): a leading loop at k is carried onto a lagging loop at -k. "
   "In x the generator -c A d_x + a T contains no J at all, so K commutes with the whole step",
   zero(K*G(k)*K - G(-k)) and zero(K*Pp(k)*K - Pm(-k)))

t = sp.symbols('t', real=True)
ok("1c. S with t -> -t, the process run backward: S G(k) S = -G(k) and S commutes with J, so if r(t) solves d_t r = G r then S r(-t) solves "
   "it too, and S P+(k) S = P-(k): the same loop run backward is the lagging loop at the same k. S lambda = %s is in the lagging class"
   % list(S*lam),
   zero(S*G(k)*S + G(k)) and S*J == J*S and zero(S*Pp(k)*S - Pm(k)) and T*(S*lam) == J*(S*lam))

ok("1d. the product K S keeps the class and reverses the motion: (KS) P+(k) (KS)^-1 = P+(-k), (KS) P-(k) (KS)^-1 = P-(-k). "
   "So 'twin reversed' and 'process run backward' are the same reversal up to turning the motion round",
   zero(K*S*Pp(k)*S*K - Pp(-k)) and zero(K*S*Pm(k)*S*K - Pm(-k)))

Gf = lambda kk, ph: G(kk) + ph*J                                             # V.2's field: a turn of the route-twin pair at the rate phi
ok("1e. with V.2's field (a turn phi J of the pair): K G_phi(k) K = G_-phi(-k) and -S G_phi(k) S = G_-phi(k): both reversals reverse the field. "
   "The reversed loop is the loop of the opposite charge",
   zero(K*Gf(k, phi)*K - Gf(-k, -phi)) and zero(-S*Gf(k, phi)*S - Gf(k, -phi)))

# the walk's own clock: M_k = -ck A + a S; (5): the growing combination is the lagging loop.  Same pairing there.
Mk = -c*k*A + a*S
ev_t = sorted([complex(e) for e, m in (J*G(k)).subs({k: 0.7, a: 1.0, c: 1.0}).evalf().eigenvals().items() for _ in range(m)], key=lambda z: z.real)
ev_c = sorted([complex(e) for e, m in Mk.subs({k: 0.7, a: 1.0, c: 1.0}).evalf().eigenvals().items() for _ in range(m)], key=lambda z: z.real)
ok("1f. the same pairing in the walk's own clock: Tr M_k = 0, and the rates of M_k are +-w, each twice, the same list as the energies of J G(k): %s and %s (ck/a = 0.7)"
   % (["%+.4f" % z.real for z in ev_c], ["%+.4f" % z.real for z in ev_t]),
   Mk.trace() == 0 and max(abs(x + y) for x, y in zip(ev_c, ev_t[::-1])) < 1e-12)

# ---- on a ring: III.15a's tick (the U of VI.4 (1)), real and orthogonal on one column; route and twin together are U on a complex column
def tick(L, th):
    dim = 2*L; U = np.zeros((dim, dim)); cc, ss = np.cos(th), np.sin(th)
    for y in range(L):
        yp, ym = (y - 1) % L, (y + 1) % L
        U[2*y, 2*yp] += cc; U[2*y, 2*ym + 1] += ss; U[2*y + 1, 2*yp] += -ss; U[2*y + 1, 2*ym + 1] += cc
    return U
def phases(U): return np.sort(-np.angle(np.linalg.eigvals(U)))
rows = []
for L in (4, 7, 12):
    U = tick(L, 0.6); eps = phases(U)
    rows.append((L, abs(eps + eps[::-1]).max(), abs(eps.sum()), np.abs(np.abs(eps) - np.pi).min(), abs(np.linalg.det(U) - 1)))
ok("1g. on a ring (4, 7, 12 sites, coin angle 0.6) the tick is real, so K - conjugation - maps a column of quasi-energy eps onto one of -eps: "
   "the 2L quasi-energies pair off (largest mismatch %.0e) and sum to 0 (%.0e); none sits at +-pi (nearest %.2f away); det U = 1 (%.0e)"
   % (max(r[1] for r in rows), max(r[2] for r in rows), min(r[3] for r in rows), max(r[4] for r in rows)),
   all(r[1] < 1e-12 and r[2] < 1e-12 and r[3] > 0.1 and r[4] < 1e-12 for r in rows))

L = 6; U = tick(L, 0.6); rng = np.random.default_rng(3); th = 0.3*rng.normal(size=L)
Dp = np.diag(np.repeat(np.exp(1j*th), 2)); Dm = Dp.conj()
e_p, e_m = phases(U@Dp), phases(U@Dm)
ok("1h. with a field on the ring (V.2's turn e^{i theta(x)} of the pair on each site, random theta): the quasi-energies in the field theta are "
   "minus those in the field -theta (%.0e); in one field alone they no longer pair (mismatch %.2f)"
   % (abs(e_p + e_m[::-1]).max(), abs(e_p + e_p[::-1]).max()), abs(e_p + e_m[::-1]).max() < 1e-12 and abs(e_p + e_p[::-1]).max() > 1e-3)

# the process run backward on the lattice: exchange the two senses on every site, with the coin moved to the other side of the shift
Sx = np.kron(np.eye(L), np.array([[0, 1.], [1, 0]])); R = np.kron(np.eye(L), np.array([[np.cos(0.6), np.sin(0.6)], [-np.sin(0.6), np.cos(0.6)]]))
X = R@Sx
ok("1i. the process run backward on the ring: with X = (coin)(exchange of the two senses on every site), X U X^-1 = U^-1 to %.0e"
   % mx(X@U@np.linalg.inv(X) - U.T), mx(X@U@np.linalg.inv(X) - U.T) < 1e-13)

# =====================================================================================================================================
# 2. the count of configurations, minor by minor (VI.5 (1)), 4 sites: 8 slots, 256 configurations
# =====================================================================================================================================
L = 4; U = tick(L, 0.6).astype(complex); D = 2*L
ev, V = np.linalg.eig(U); eps = -np.angle(ev); order = np.argsort(eps); eps = eps[order]
def compound(U, N):                                                          # the N-walker evolution: all N x N minors of U
    idx = list(itertools.combinations(range(D), N))
    if N == 0: return np.ones((1, 1), complex)
    return np.array([[np.linalg.det(U[np.ix_(r, cl)]) for cl in idx] for r in idx])
worst = 0.0; nconf = 0; allE = []
for N in range(D + 1):
    CN = compound(U, N); got = list(np.linalg.eigvals(CN))
    sums = np.array([eps[list(O)].sum() for O in itertools.combinations(range(D), N)]); allE += list(sums); nconf += len(sums)
    for z in np.exp(-1j*sums):                                               # match the two lists of factors e^{-iE} one to one
        j = int(np.argmin(np.abs(np.array(got) - z))); worst = max(worst, abs(got[j] - z)); got.pop(j)
    if N == 0: empty_factor = CN[0, 0]
    if N == D: full_factor = CN[0, 0]
allE = np.array(allE)
ok("2a. the N-walker evolution built from the N x N minors of the tick, N = 0..8 (%d configurations): its factors are e^{-iE} with E the sum of the "
   "one-walker quasi-energies over the slots occupied, matched one to one to %.0e. No walker: the factor is 1 (energy 0, no constant to choose). "
   "Every slot occupied: det U = %.12f, energy 0 again" % (nconf, worst, full_factor.real),
   nconf == 256 and worst < 1e-10 and empty_factor == 1 and abs(full_factor - 1) < 1e-12)

E_sea = eps[eps < 0].sum(); half = -0.5*np.abs(eps).sum()
symm = abs(np.sort(allE) + np.sort(allE)[::-1]).max()
n_low = int((allE < E_sea + 1e-9).sum())
ok("2b. the lowest of the 256: every lagging (eps < 0) slot occupied, E = %.6f; it equals -(1/2) sum|eps| = %.6f (difference %.0e) because the "
   "sum over all slots is 0; it is the only configuration at that energy (%d); the 256 energies are symmetric about 0 (%.0e): "
   "the complement of a configuration has minus its energy" % (E_sea, half, abs(E_sea - half), n_low, symm),
   abs(E_sea - half) < 1e-12 and n_low == 1 and symm < 1e-12 and allE.min() > E_sea - 1e-12)

# relative to the sea: cost = sum |eps| over (leading slots present) + (lagging slots absent)
worst_rel = 0.0
for N in range(D + 1):
    for O in itertools.combinations(range(D), N):
        Oset = set(O); cost = sum(abs(eps[i]) for i in range(D) if (eps[i] > 0 and i in Oset) or (eps[i] < 0 and i not in Oset))
        worst_rel = max(worst_rel, abs(eps[list(O)].sum() - E_sea - cost))
pairs = [(i, D - 1 - i) for i in range(D//2)]                                 # eps sorted: slot i and slot D-1-i are a reversal pair
pair_dev = max(abs(eps[i] + eps[j]) for i, j in pairs)
i0, j0 = pairs[0]
four = [0.0, eps[j0], eps[i0], eps[i0] + eps[j0]]                             # neither, leading only, lagging only, both: from the empty cut
ok("2c. relative to that configuration every one of the 256 costs sum|eps| over (leading loops present) + (lagging loops absent), to %.0e - never "
   "negative. For one reversal pair (eps = %+.4f, %+.4f) the four fillings 'neither, leading only, lagging only, both' are %s from the empty cut "
   "and %s from the sea: 'neither' (a lagging loop missing) and 'both' (a leading loop added) cost the same"
   % (worst_rel, eps[i0], eps[j0], ["%+.4f" % x for x in four], ["%+.4f" % (x - eps[i0]) for x in four]),
   worst_rel < 1e-12 and pair_dev < 1e-12 and abs((four[0] - eps[i0]) - (four[3] - eps[i0])) < 1e-12)

# the same four fillings along the walk's own clock (VI.5 (5)): weights 1, e^{-w tau}, e^{+w tau}, 1 = e^{-E tau} with the same four E
kk, aa, tau = 0.7, 1.0, 3.0; ww = np.hypot(aa, kk)
M2 = np.array([[-kk, aa], [aa, kk]]); Lk = sl.expm(tau*M2); g = np.sort(np.linalg.eigvals(Lk).real)   # one column's clock transfer, rates -+w
clock = [1.0, g[0], g[1], np.linalg.det(Lk)]
ok("2d. the same four fillings counted along the walk's own clock, L_k(tau) = e^{tau M_k} (ck/a = 0.7, a tau = 3): weights %s = e^{-E tau} with "
   "E = 0, +w, -w, 0 (to %.0e). 'Both' weighs det L_k = 1 because Tr M_k = 0: the pairing again. The filling the long clock selects is the lagging one"
   % (["%.4f" % x for x in clock], max(abs(clock[1] - np.exp(-ww*tau)), abs(clock[2] - np.exp(ww*tau)), abs(clock[3] - 1))),
   max(abs(clock[1] - np.exp(-ww*tau)), abs(clock[2] - np.exp(ww*tau))) < 1e-9 and abs(clock[3] - 1) < 1e-12)

# a fixed number N of walkers: the settled configuration is the N lowest slots; half filling is the lowest of all
byN = {}
for N in range(D + 1):
    byN[N] = min(eps[list(O)].sum() for O in itertools.combinations(range(D), N))
lowestN = [np.sort(eps)[:N].sum() for N in range(D + 1)]
ok("2e. with the number of walkers fixed at N = 0..8 the lowest configuration is the N lowest slots, energies %s (in units of pi: %s); "
   "N = 0 and N = 8 give 0, every N from 1 to 7 a negative settled energy, and half filling, N = 4, the lowest"
   % (["%.4f" % byN[N] for N in range(D + 1)], ["%.3f" % (byN[N]/np.pi) for N in range(D + 1)]),
   max(abs(byN[N] - lowestN[N]) for N in range(D + 1)) < 1e-12 and abs(byN[0]) < 1e-15 and abs(byN[D]) < 1e-12
   and all(byN[N] < -1e-6 for N in range(1, D)) and min(byN, key=byN.get) == D//2)

# the same total said two ways (VI.5 (1)-(2)): every filling of the list of slots, or one possible walker round a closed history
Ub = 0.5*U                                                                    # any step with eigenvalues inside the unit circle, so that the series converges
snap = sum(np.trace(compound(Ub, N)) for N in range(D + 1))                   # the snapshot: all 256 fillings, minor by minor
detsum = np.linalg.det(np.eye(D) + Ub)
Up = np.eye(D, dtype=complex); series = 0.0
for n in range(1, 120):
    Up = Up@Ub; series += (-1)**(n + 1)*np.trace(Up)/n                        # one walker going round n times, with the sign of the n-cycle
tau_ = 3.0
pair_id = abs((1 + np.exp(tau_))*(1 + np.exp(-tau_)) - (1 + np.exp(-tau_) + np.exp(tau_) + 1))
ok("2f. the same total said two ways: summed over all 256 fillings (principal minors, every number of walkers) = det(1 + U) (%.0e) = "
   "exp[sum_n (-1)^(n+1) Tr U^n / n], one possible walker round a closed history n times (%.0e; the step scaled by 1/2 so that the series converges); "
   "for one pair of slots (1 + e^{a tau})(1 + e^{-a tau}) is the sum of the four fillings' weights (%.0e)"
   % (abs(snap - detsum), abs(np.exp(series) - detsum), pair_id), abs(snap - detsum) < 1e-12 and abs(np.exp(series) - detsum) < 1e-12 and pair_id < 1e-12)

# =====================================================================================================================================
# 3. three dimensions: IV.4b's real matrices, IV.5a's eight counters
# =====================================================================================================================================
Ax = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0.]]); Ay = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0.]])
Az = np.diag([1., -1., -1., 1.]); T4 = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0.]]); E4 = np.eye(4)
big = lambda M: np.kron(np.eye(2), M)                                         # the same matrix on route and on twin
J8 = np.kron(np.array([[0, -1.], [1, 0]]), E4); K8 = np.kron(np.diag([1., -1.]), E4)
As = [Ax, Ay, Az]
def G8(kv, ph=0.0): return -sum(kv[j]*big(As[j]) for j in range(3))@J8 + 1.0*big(T4) + ph*J8      # c = a = 1
Sig = [-J8@big(Ay@Az), -J8@big(Az@Ax), -J8@big(Ax@Ay)]                        # the axis: -J A_j A_k, real and symmetric (with J as i: -i A_j A_k)
kv = np.array([0.3, -0.5, 0.4]); w3 = np.sqrt(1 + kv@kv)
ok("3a. eight real counters: K = diag(1_4, -1_4) anticommutes with J, K G(k) K = G(-k) (%.0e) and with a field K G_phi(k) K = G_-phi(-k) (%.0e); "
   "the energy matrix J G(k) squares to w^2 (%.0e) and has trace %.0e: four +w and four -w"
   % (mx(K8@G8(kv)@K8 - G8(-kv)), mx(K8@G8(kv, 0.37)@K8 - G8(-kv, -0.37)), mx((J8@G8(kv))@(J8@G8(kv)) - w3**2*np.eye(8)), abs(np.trace(J8@G8(kv)))),
   mx(K8@J8 + J8@K8) == 0 and mx(K8@G8(kv)@K8 - G8(-kv)) < 1e-14 and mx(K8@G8(kv, 0.37)@K8 - G8(-kv, -0.37)) < 1e-14
   and mx((J8@G8(kv))@(J8@G8(kv)) - w3**2*np.eye(8)) < 1e-13 and abs(np.trace(J8@G8(kv))) < 1e-13)

alg = max(mx(Sig[i]@Sig[i] - np.eye(8)) for i in range(3)) + max(mx(Sig[i]@Sig[j] - J8@Sig[l]) for i, j, l in ((0, 1, 2), (1, 2, 0), (2, 0, 1)))
ok("3b. the axis Sigma_j = -J A_k A_l is real and symmetric, Sigma_j^2 = 1, Sigma_x Sigma_y = J Sigma_z and cyclic (%.0e); it commutes with J and with "
   "the mass turn T; and K Sigma_j K = -Sigma_j for all three: the reversal that reverses the energy reverses the axis"
   % alg, alg < 1e-14 and all(mx(s - s.T) == 0 and mx(s@J8 - J8@s) == 0 and mx(s@big(T4) - big(T4)@s) == 0 and mx(K8@s@K8 + s) == 0 for s in Sig))

kz = np.array([0, 0, 0.8]); wz = np.sqrt(1 + 0.64)
def proj(kvec, Esign, ssign): return (np.eye(8) + Esign*J8@G8(kvec)/wz)/2 @ (np.eye(8) + ssign*Sig[2])/2
devs = [mx(K8@proj(kz, e, s)@K8 - proj(-kz, -e, -s)) for e in (1, -1) for s in (1, -1)]
rank = [int(round(np.trace(proj(kz, e, s)))) for e in (1, -1) for s in (1, -1)]
ok("3c. for motion along z the classes are labelled by (sign of energy, axis along z), each a real plane (traces %s); K carries "
   "(E, k, s) onto (-E, -k, -s) for all four (%.0e). So a missing lagging loop with (-E, -k, -s) is a present loop with (+E, +k, +s) and the opposite charge"
   % (rank, max(devs)), max(devs) < 1e-14 and rank == [2, 2, 2, 2] and mx(Sig[2]@G8(kz) - G8(kz)@Sig[2]) < 1e-14)

# K in IV.4's chiral notation: psi = V z, z = u + i v.  V intertwines IV.4b's real matrices with IV.4's chiral ones.
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1]).astype(complex); T2 = np.array([[0, 1.], [-1, 0]])
Ach = [np.kron(sz, s) for s in (sx, sy, sz)]; Tch = np.kron(T2, np.eye(2)).astype(complex)      # A_k = diag(sigma_k, -sigma_k), T = turn between the copies
rows_ = []
for Mr, Mc in zip(As + [T4], Ach + [Tch]):                                    # V M_real = M_chiral V   as a linear system on vec(V)
    rows_.append(np.kron(np.eye(4), Mc) - np.kron(Mr.T, np.eye(4)))
_, sv, Vh = np.linalg.svd(np.vstack(rows_)); null = Vh[-1].conj().reshape(4, 4).T
Vmap = null*2/np.linalg.norm(null)                                             # unitary up to a phase
inter = max(mx(Vmap@Mr - Mc@Vmap) for Mr, Mc in zip(As + [T4], Ach + [Tch]))
C = Vmap@Vmap.T                                                                # z -> conj(z)  reads  psi -> C conj(psi)
TT = np.kron(T2, T2); ph = np.vdot(TT, C)/np.vdot(TT, TT)
ok("3d. in IV.4's chiral notation (psi = V z; V unique up to a phase: null space of dimension %d, V unitary to %.0e, intertwining to %.0e) "
   "K reads psi -> C conj(psi) with C = (phase) x [turn matrix on the two senses] x [turn matrix on the two axis components] to %.0e: "
   "conjugate, and apply T twice - once between the copies, once inside each"
   % (int((sv < 1e-10).sum()), mx(Vmap@Vmap.conj().T - np.eye(4)), inter, mx(C - ph*TT)),
   int((sv < 1e-10).sum()) == 1 and mx(Vmap@Vmap.conj().T - np.eye(4)) < 1e-12 and inter < 1e-12 and abs(abs(ph) - 1) < 1e-12 and mx(C - ph*TT) < 1e-12)

sig = [sx, sy, sz]
flip = max(mx(T2@s.conj()@T2.T + s) for s in sig)
worst0 = 0.0; worstJ = 0.0
for _ in range(20):
    chi = rng.normal(size=2) + 1j*rng.normal(size=2); rev = T2@chi.conj()
    pair = np.outer(chi, rev) - np.outer(rev, chi)                            # the odd combination of the two, VI.3 (0)
    worst0 = max(worst0, mx(pair + np.vdot(chi, chi).real*T2))
    for s in sig: worstJ = max(worstJ, mx(s@pair + pair@s.T))                 # (Sigma on walker 1 + Sigma on walker 2) on the pair
    ax = np.array([np.vdot(chi, s@chi).real for s in sig]); axr = np.array([np.vdot(rev, s@rev).real for s in sig]); worstJ = max(worstJ, mx(ax + axr))
ok("3e. on the axis pair alone: T sigma_k* T^-1 = -sigma_k (%.0e), so the reversed column T conj(chi) has the opposite axis; and for 20 random chi "
   "the odd combination of chi and its reversal is -|chi|^2 T (to %.0e): VI.3 (0)'s T column with full weight, whatever the axis, total axis 0 (%.0e)"
   % (flip, worst0, worstJ), flip < 1e-15 and worst0 < 1e-14 and worstJ < 1e-14)

# =====================================================================================================================================
# 4. in motion: the equilibrium of the moving observer's flip process, III.15d (12)
# =====================================================================================================================================
rows4 = []
for eta in (0.0, 0.5, 1.3, 2.0):
    aa = 1.0; ck = aa*np.sinh(eta); om = aa*np.cosh(eta)
    Mk = np.array([[-ck, aa], [aa, ck]]); Bm = np.diag([np.exp(-eta/2), np.exp(eta/2)])
    gen = np.array([[-aa*np.exp(eta), aa*np.exp(-eta)], [aa*np.exp(eta), -aa*np.exp(-eta)]])
    d12 = mx(gen - (-om*np.eye(2) + Bm@Mk@np.linalg.inv(Bm)))
    evg, Vg = np.linalg.eig(gen); i0 = int(np.argmax(evg.real)); pi = Vg[:, i0].real; pi = pi/pi.sum()
    evm, Vm = np.linalg.eigh(Mk)                                              # ascending: -w, +w
    grow, decay = Vm[:, 1], Vm[:, 0]
    rows4.append((eta, d12, abs(np.sort(evg.real)[1]), abs(np.sort(evg.real)[0] + 2*om), abs(pi[0] - np.exp(-eta)/(2*np.cosh(eta))),
                  abs((pi[1] - pi[0]) - np.tanh(eta)), abs(pi[0]*aa*np.exp(eta) + pi[1]*aa*np.exp(-eta) - aa/np.cosh(eta)),
                  grow[0]*grow[1] > 0, decay[0]*decay[1] < 0))
ok("4a. the moving observer's two-rate flip process (rapidity 0, 0.5, 1.3, 2): its generator is -w + B M_k B^-1 (%.0e), columns sum to 0; rates 0 and -2w "
   "(%.0e, %.0e); equilibrium (e^-eta, e^eta)/(2 cosh eta) (%.0e): imbalance of the two senses = tanh(eta) = v/c (%.0e), flips per unit clock = a/cosh(eta) "
   "= a/gamma (%.0e). The growing column of M_k has both entries of one sign and the decaying one a change of sign, at every rapidity"
   % tuple(max(r[i] for r in rows4) for i in range(1, 7)),
   all(r[1] < 1e-13 and r[2] < 1e-13 and r[3] < 1e-13 and r[4] < 1e-13 and r[5] < 1e-13 and r[6] < 1e-13 and r[7] and r[8] for r in rows4))

def visible(eta, aa, Ttot, Ns, runs, rng):                                    # simulate the two-rate process started at equilibrium; count blocks that end reversed
    rate = (aa*np.exp(eta), aa*np.exp(-eta)); pi_p = np.exp(-eta)/(2*np.cosh(eta)); acc = {N: 0.0 for N in Ns}
    for _ in range(runs):
        sgn = 0 if rng.random() < pi_p else 1; t = 0.0; times = [0.0]; states = [sgn]
        while True:
            t += rng.exponential(1/rate[sgn])
            if t >= Ttot: break
            sgn = 1 - sgn; times.append(t); states.append(sgn)
        times = np.array(times); states = np.array(states)
        for N in Ns:
            edges = np.linspace(0, Ttot, N + 1); st = states[np.searchsorted(times, edges, side='right') - 1]
            acc[N] += np.sum(st[1:] != st[:-1])
    return {N: acc[N]/runs for N in Ns}
eta, aa, Ttot, Ns = 0.8, 1.0, 20.0, (400, 40, 10, 2)
sim = visible(eta, aa, Ttot, Ns, 6000, np.random.default_rng(7)); om = aa*np.cosh(eta); pp_ = np.exp(-eta)/(2*np.cosh(eta))
form = {N: N*2*pp_*(1 - pp_)*(1 - np.exp(-2*om*Ttot/N)) for N in Ns}
ok("4b. blocks of the two-rate process at rapidity 0.8 (6000 runs, an interval of 20/a cut into 400, 40, 10, 2 blocks): blocks that come out reversed, "
   "simulated %s against 2 pi+ pi- (1 - e^{-2 w tau}) per block %s; true flips in the interval a T/gamma = %.2f. The visible flip rate falls as the "
   "block grows, at the rate 2w" % (["%.2f" % sim[N] for N in Ns], ["%.2f" % form[N] for N in Ns], aa*Ttot/np.cosh(eta)),
   all(abs(sim[N] - form[N]) < 0.03*form[N] + 0.03 for N in Ns) and sim[400] > sim[40] > sim[10] > sim[2])

Ns0 = (2000, 100, 20, 10, 4, 1)
sim0 = visible(0.0, 1.0, 20.0, Ns0, 6000, np.random.default_rng(8)); form0 = {N: N*(1 - np.exp(-2*20.0/N))/2 for N in Ns0}
gen0 = np.array([[-1., 1.], [1., -1.]]); ev0 = np.sort(np.linalg.eigvals(gen0).real)
Pt = sl.expm(1.7*gen0)@np.array([1., 0.]); Pt_formula = 0.5*np.array([1., 1.]) + 0.5*np.exp(-2*1.7)*np.array([1., -1.])
ok("4c. at rest, III.2's lines with the loss in: generator a [[-1, 1], [1, -1]], rates %.0f and %.0f (in units of a); from the forward sense "
   "P(tau) = (1,1)/2 + e^{-2 a tau} (1,-1)/2 (%.0e). An interval with a T = 20 cut into %s blocks: blocks that come out reversed, simulated %s, "
   "formula N (1 - e^{-2 a T/N})/2 = %s" % (ev0[1], ev0[0], mx(Pt - Pt_formula), list(Ns0), ["%.2f" % sim0[N] for N in Ns0], ["%.2f" % form0[N] for N in Ns0]),
   abs(ev0[1]) < 1e-14 and abs(ev0[0] + 2) < 1e-14 and mx(Pt - Pt_formula) < 1e-14 and all(abs(sim0[N] - form0[N]) < 0.03*form0[N] + 0.03 for N in Ns0))

# =====================================================================================================================================
# 5. the pair made at a vertex, three dimensions
# =====================================================================================================================================
a3 = 1.0
Hc = lambda pv: sum(pv[j]*As[j] for j in range(3)) + 1j*a3*T4                 # one complex column z = u + i v: i d_t z = H z (IV.5a (5), J written as i)
wc = lambda pv: np.sqrt(a3*a3 + pv@pv)
PIp = lambda pv: (np.eye(4) + Hc(pv)/wc(pv))/2
PIm = lambda pv: (np.eye(4) - Hc(pv)/wc(pv))/2
Sg = [-1j*Ay@Az, -1j*Az@Ax, -1j*Ax@Ay]
u_, s_, vh_ = np.linalg.svd(PIp(np.zeros(3))@(np.eye(4) + Sg[2])/2); up = u_[:, 0]; Brest = np.stack([up, Sg[0]@up], axis=1)   # rest columns, axis +z and -z
pauli = max(mx(Brest.conj().T@Sg[j]@Brest - sig[j]) for j in range(3))
def ubasis(pv): return PIp(pv)@Brest*np.sqrt(2*wc(pv)/(wc(pv) + a3))          # the rest columns carried to momentum p
def pair_axis(pv, qv, Vtx):
    """sea column at q -> empty column at p through the vertex Vtx.  The gap at q is a lagging loop missing (a hole) = a leading loop conj(o) at -q (check 3c).
    Returns the 2x2 axis part c: pair = sum c_ss' (loop, axis s, at p) x (loop, axis s', at -q)."""
    return ubasis(pv).conj().T@Vtx@PIm(qv)@ubasis(-qv).conj()
def split(cm): return np.trace(T2.T@cm)/2, np.array([np.trace((sg@T2).conj().T@cm)/2 for sg in sig])
rng5 = np.random.default_rng(11); worst_single = 0.0; size = 0.0; errf = {}; errc = {}
for sc in (1.0, 0.3, 0.03, 0.003):
    errf[sc] = errc[sc] = 0.0
    for _ in range(10):
        pv, qv = sc*rng5.normal(size=3), sc*rng5.normal(size=3)
        al_f, be_f = split(pair_axis(pv, qv, np.eye(4))); al_c, be_c = split(pair_axis(pv, qv, 1j*T4))
        worst_single = max(worst_single, abs(al_f), abs(al_c)); size = max(size, np.abs(be_f).max())
        errf[sc] = max(errf[sc], mx(be_f - (pv - qv)/(2*a3))/sc); errc[sc] = max(errc[sc], mx(be_c + (pv + qv)/(2*a3))/sc)
worst_field, worst_coin = errf[0.003], errc[0.003]
ok("5a. the axis part of the pair, from the projectors (rest columns carry Sigma as the Pauli matrices to %.0e; 40 random pairs of momenta, |p| from 0.003 to 1 "
   "in units of mc): the part along T, the column with no net axis, is %.0e for V.2's field vertex and for V.1's coin vertex alike - zero, at every "
   "momentum. The pair is in the symmetric columns sigma_k T: for the field, vector (p - q)/2mc, the momentum the field supplies (relative error %.0e at "
   "|p| ~ 0.003 mc, %.0e at 0.03: the next order); for the coin, vector -(p + q)/2mc, the relative momentum of the leading loop and the hole (%.0e, %.0e)"
   % (pauli, worst_single, worst_field, errf[0.03], worst_coin, errc[0.03]),
   pauli < 1e-14 and worst_single < 1e-15 and size > 0.1 and worst_field < 1e-4 and worst_coin < 1e-4 and 30 < errf[0.03]/worst_field < 300 and 30 < errc[0.03]/worst_coin < 300)

# why the T part is zero: the mirror.  The pair's axis part changes sign when both momenta are reversed, and turns as a vector under rotations;
# a T part would be a number built from two vectors that changes sign in the mirror, and there is none.
odd = 0.0; rot = 0.0
for _ in range(10):
    pv, qv = 0.5*rng5.normal(size=3), 0.5*rng5.normal(size=3)
    for Vtx in (np.eye(4), 1j*T4):
        c1 = pair_axis(pv, qv, Vtx); odd = max(odd, mx(pair_axis(-pv, -qv, Vtx) + c1))
        ang = rng5.normal(size=3); th_ = np.linalg.norm(ang); nh = ang/th_
        Kx = np.array([[0, -nh[2], nh[1]], [nh[2], 0, -nh[0]], [-nh[1], nh[0], 0]]); Rm = np.eye(3) + np.sin(th_)*Kx + (1 - np.cos(th_))*Kx@Kx
        D2 = sl.expm(-0.5j*th_*sum(nh[j]*sig[j] for j in range(3)))
        rot = max(rot, mx(pair_axis(Rm@pv, Rm@qv, Vtx) - D2@c1@D2.T))
ok("5d. why: with both momenta reversed the pair's axis part changes sign (%.0e), and under a rotation of the momenta it turns as two axis columns do, "
   "c -> D c D^T (%.0e). A part along T would be a number made from two vectors that changes sign in the mirror; there is none" % (odd, rot),
   odd < 1e-14 and rot < 1e-13)

# the many-walker run: 3 x 3 x 3 sites, the sea of 54 columns, one tick with the vertex on one site; amplitudes by determinants
n3 = 3; a3 = 4.0; N3 = n3**3
S1 = np.roll(np.eye(n3), 1, axis=0); d1 = (S1.T - S1)/2; E3 = np.eye(n3)
dxs = [np.kron(np.kron(d1, E3), E3), np.kron(np.kron(E3, d1), E3), np.kron(np.kron(E3, E3), d1)]
Hl = sum(np.kron(-1j*dxs[j], As[j]) for j in range(3)) + 1j*a3*np.kron(np.eye(N3), T4)       # Hermitian, 108 x 108, ordering (site, counter)
el, Vl = np.linalg.eigh(Hl); Phi = Vl[:, el < 0]; Emp = Vl[:, el > 0]; e_sea, e_emp = el[el < 0], el[el > 0]
Ul = (Vl*np.exp(-1j*el))@Vl.conj().T
x0 = 13; Psite = np.zeros((4*N3, 4*N3)); Psite[4*x0:4*x0 + 4, 4*x0:4*x0 + 4] = np.eye(4)
theta = 1e-3
def run(Vtx4):
    Gv = np.kron(np.eye(N3), Vtx4)@Psite                                      # the vertex on site x0: theta * Gv is the angle turned in one tick
    Wl = Ul@sl.expm(1j*theta*Gv)
    den = Phi.conj().T@Wl@Phi; R = (Emp.conj().T@Wl@Phi)@np.linalg.inv(den)   # R[e, o] = det(sea with o replaced by e | W | sea) / det(sea | W | sea)
    pred = 1j*theta*np.exp(-1j*e_emp)[:, None]*(Emp.conj().T@Gv@Phi)*np.exp(1j*e_sea)[None, :]
    brute = 0.0
    for (ie, io) in ((0, 0), (5, 17), (53, 2), (20, 40)):                     # the same ratio from the two determinants, written out
        Phi2 = Phi.copy(); Phi2[:, io] = Emp[:, ie]
        brute = max(brute, abs(np.linalg.det(Phi2.conj().T@Wl@Phi)/np.linalg.det(den) - R[ie, io]))
    Psi = Emp@R@Phi.conj().T                                                  # the pair: sum R[e,o] e x conj(o), rows = leading loop, columns = hole
    PiP = Emp@Emp.conj().T; PiM = Phi@Phi.conj().T
    return mx(R - pred), mx(R), brute, mx(Psi - 1j*theta*Ul@PiP@Gv@PiM@Ul.conj().T), Psi
out_f = run(np.eye(4)); out_c = run(1j*T4)
ok("5b. many walkers, 27 sites, the sea of %d columns, one tick, vertex of angle %.0e on one site: the %d pair amplitudes (determinant ratios; four of them "
   "also from the two determinants written out, %.0e) equal i theta e^{-i eps_e} <e|vertex|o> e^{+i eps_o} to %.0e (field) and %.0e (coin), against a largest "
   "amplitude of %.0e and %.0e: the pair is i theta U Pi_+ (vertex) Pi_- U^-1 (%.0e, %.0e), second order in theta"
   % (Phi.shape[1], theta, Phi.shape[1]*Emp.shape[1], max(out_f[2], out_c[2]), out_f[0], out_c[0],
      out_f[1], out_c[1], out_f[3], out_c[3]),
   Phi.shape[1] == 54 and max(out_f[2], out_c[2]) < 1e-12 and out_f[0] < 20*theta*out_f[1] and out_c[0] < 20*theta*out_c[1] and out_f[3] < 20*theta*out_f[1]
   and out_c[3] < 20*theta*out_c[1])

# the lattice pair, momentum by momentum, is the projector formula of 5a with p -> sin p (central differences)
ks = 2*np.pi*np.arange(n3)/n3
def plane(kv):                                                                # 108 x 4: a plane wave on the 27 sites times the four counters
    ph = np.array([np.exp(1j*(kv[0]*ix + kv[1]*iy + kv[2]*iz)) for ix in range(n3) for iy in range(n3) for iz in range(n3)])/np.sqrt(N3)
    return np.kron(ph[:, None], np.eye(4))
a3_keep = a3
worst_blk = 0.0; worst_s = 0.0; big_t = 0.0
PiP = Emp@Emp.conj().T; PiM = Phi@Phi.conj().T
for Vtx4 in (np.eye(4), 1j*T4):
    Gv = np.kron(np.eye(N3), Vtx4)@Psite; core = PiP@Gv@PiM
    for ip in itertools.product(range(n3), repeat=3):
        for iq in itertools.product(range(n3), repeat=3):
            kp, kq = ks[list(ip)], ks[list(iq)]; sp_, sq_ = np.sin(kp), np.sin(kq)
            blk = plane(kp).conj().T@core@plane(kq)*N3                        # 4 x 4 block: empty column at kp <- sea column at kq
            x0v = np.array([x0//9, (x0//3) % 3, x0 % 3]); phase = np.exp(-1j*(kp - kq)@x0v)
            worst_blk = max(worst_blk, mx(blk - phase*PIp(sp_)@Vtx4@PIm(sq_)))
            al, be = split(pair_axis(sp_, sq_, Vtx4)); worst_s = max(worst_s, abs(al)); big_t = max(big_t, np.abs(be).max())
ok("5c. the lattice pair of 5b, block by block (27 x 27 pairs of lattice momenta, both vertices): the block from kq to kp is Pi_+(sin kp) (vertex) Pi_-(sin kq) "
   "times the phase of the site, to %.0e - 5a's formula with sin k for k. Its part along T is %.0e in every block; its symmetric part reaches %.2f"
   % (worst_blk, worst_s, big_t), worst_blk < 1e-12 and worst_s < 1e-15 and big_t > 0.05)

print("ALL PASS" if ok_all else "SOME FAILED"); raise SystemExit(0 if ok_all else 1)
