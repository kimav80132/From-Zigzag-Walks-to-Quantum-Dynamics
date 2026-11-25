"""Exploration, 11 Nov 2026 (not cited in the master yet): the count of fillings as a partition function; free energy; three entropies.

Vocabulary: slots (a leading and a lagging one per k), fillings (which slots hold a walker), the walk's own clock tau = x/c, L_k(tau) = e^{tau M_k}
(III.15d), the settled configuration (VI.5, VI.7).  Units hbar = c = 1.

  1. One pair of slots.  Z = sum over the four fillings = det(1 + L_k) = (1 + e^{w tau})(1 + e^{-w tau}).  The chance that a slot of energy E holds a
     walker is 1/(e^{E tau} + 1): the Fermi function, with the clock length in the place of 1/temperature.  Free energy F = -(1/tau) log Z -> -w = E_0.
     The entropy of the four fillings, -sum p log p, equals log Z + tau <E> and falls from 2 log 2 (short clock) to 0 (long clock).
  2. Three pairs (64 fillings, minor by minor): Z is the product, the fillings are independent slot by slot, the entropy adds up.
  3. Many pairs, massless: the entropy of the count per unit extent is pi/(3 tau), and F - E_0 per unit extent is -pi/(6 tau^2): the number c = 1.
  4. A different entropy: the sense of ONE walker (VI.7 (3)).  It rises from 0 to log 2 while the distance to equilibrium (relative entropy) only falls;
     in motion the equilibrium is lopsided and the settled sense has entropy h((1 - v/c)/2) < log 2.
  5. A third entropy: the settled configuration seen from inside a region of l sites (what the region shares with the rest), from the restriction of
     Pi_- to the region.  Continuum lines on a ring: with almost no mass it grows as (2/3) log l - two edges in k, the one at k = 0 and the sharp cut at
     the largest k; a mass closes the first and leaves (1/3) log l; at fixed large l the entropy is (1/3) log(1/a) + const.  The lattice tick has two
     gaps, at quasi-energy 0 and at pi, both closed by the coin angle: its entropy levels off at (2/3) log(1/theta).  One third per edge: c = 1 again.
"""
import itertools, sys
import numpy as np, scipy.linalg as sl
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all = ok_all and bool(cond); print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()
h2 = lambda p: 0.0 if p <= 0 or p >= 1 else float(-p*np.log(p) - (1 - p)*np.log(1 - p))          # entropy of a yes/no with probability p

# =====================================================================================================================================
# 1. one pair of slots
# =====================================================================================================================================
def Mk(k, a): return np.array([[-k, a], [a, k]])
rows = []
for (k, a, tau) in ((0.0, 1.0, 0.3), (0.7, 1.0, 1.0), (2.0, 1.0, 2.5), (0.4, 0.5, 6.0)):
    w = np.hypot(k, a); L = sl.expm(tau*Mk(k, a))
    fill = [1.0, np.trace(L), np.linalg.det(L)]                                # no walker, one (sum over the two slots), two: the principal minors
    Z = sum(fill); ev = np.sort(np.linalg.eigvalsh(L))                          # e^{-w tau} (leading slot), e^{+w tau} (lagging slot)
    wts = np.array([1.0, ev[0], ev[1], ev[0]*ev[1]])                           # neither, leading only, lagging only, both
    p = wts/wts.sum(); p_lag = p[2] + p[3]; p_lead = p[1] + p[3]
    E4 = np.array([0.0, w, -w, 0.0]); Emean = float(p@E4)
    Sfill = float(-(p*np.log(p)).sum())
    dlogZ = (np.log(np.linalg.det(np.eye(2) + sl.expm((tau + 1e-5)*Mk(k, a)))) - np.log(np.linalg.det(np.eye(2) + sl.expm((tau - 1e-5)*Mk(k, a)))))/2e-5
    rows.append(dict(Z=abs(Z - (1 + np.exp(w*tau))*(1 + np.exp(-w*tau)))/Z, det=abs(Z - np.linalg.det(np.eye(2) + L))/Z,
                     fermi=max(abs(p_lag - 1/(np.exp(-w*tau) + 1)), abs(p_lead - 1/(np.exp(w*tau) + 1))), indep=abs(p[3] - p_lag*p_lead),
                     E=abs(Emean + w*np.tanh(w*tau/2)), dZ=abs(Emean + dlogZ), S=abs(Sfill - (np.log(Z) + tau*Emean)), S2=abs(Sfill - h2(p_lag) - h2(p_lead)),
                     F=-np.log(Z)/tau, w=w, tau=tau, Sfill=Sfill, p_lag=p_lag))
ok("1a. one pair of slots, four choices of (k, a, tau): the sum over the four fillings (principal minors of L_k) = det(1 + L_k) = (1 + e^{w tau})(1 + e^{-w tau}) "
   "(relative %.0e, %.0e)" % (max(r['det'] for r in rows), max(r['Z'] for r in rows)), all(r['det'] < 1e-12 and r['Z'] < 1e-12 for r in rows))
ok("1b. the chance that a slot holds a walker is 1/(e^{E tau} + 1), E = -w for the lagging slot and +w for the leading one (%.0e), and the two slots are filled "
   "independently (%.0e): the Fermi function, the clock length in the place of 1/temperature" % (max(r['fermi'] for r in rows), max(r['indep'] for r in rows)),
   all(r['fermi'] < 1e-12 and r['indep'] < 1e-12 for r in rows))
ok("1c. mean energy of the fillings = -w tanh(w tau/2) (%.0e) = -d log Z/d tau (%.0e); F = -(1/tau) log Z = %s against -w = %s: F -> E_0 as the clock grows "
   "(F + w = -(2/tau) log(1 + e^{-w tau}))" % (max(r['E'] for r in rows), max(r['dZ'] for r in rows), ["%.4f" % r['F'] for r in rows], ["%.4f" % -r['w'] for r in rows]),
   all(r['E'] < 1e-12 and r['dZ'] < 1e-7 and abs(r['F'] + r['w'] + 2/r['tau']*np.log(1 + np.exp(-r['w']*r['tau']))) < 1e-12 for r in rows))
taus = np.array([1e-3, 0.3, 1.0, 3.0, 10.0, 30.0]); Sc = [2*h2(1/(1 + np.exp(-t))) for t in taus]
ok("1d. the entropy of the four fillings, -sum p log p, = log Z + tau <E> (%.0e) = the two slots' yes/no entropies added (%.0e). At w tau = %s it is %s: "
   "2 log 2 = %.4f when the clock is short (the four fillings share the count), 0 when it is long (one filling has it all)"
   % (max(r['S'] for r in rows), max(r['S2'] for r in rows), ["%g" % t for t in taus], ["%.4f" % s for s in Sc], 2*np.log(2)),
   all(r['S'] < 1e-12 and r['S2'] < 1e-12 for r in rows) and abs(Sc[0] - 2*np.log(2)) < 1e-6 and Sc[-1] < 1e-10 and all(x > y for x, y in zip(Sc[:-1], Sc[1:])))

# =====================================================================================================================================
# 2. three pairs of slots, 64 fillings, minor by minor
# =====================================================================================================================================
ks = (0.0, 0.7, 2.0); a = 1.0; tau = 1.2; ws = [np.hypot(k, a) for k in ks]
U = sl.block_diag(*[sl.expm(tau*Mk(k, a)) for k in ks]); D = 6
def minor(idx): return 1.0 if len(idx) == 0 else np.linalg.det(U[np.ix_(idx, idx)])
Zs = sum(minor(list(O)) for N in range(D + 1) for O in itertools.combinations(range(D), N))
# in the basis of the slots themselves (each L_k diagonal) a filling's weight is the product of its slots' factors
fac = np.concatenate([[np.exp(-w*tau), np.exp(w*tau)] for w in ws]); Es = np.concatenate([[w, -w] for w in ws])
probs = {};
for N in range(D + 1):
    for O in itertools.combinations(range(D), N): probs[O] = np.prod(fac[list(O)])
Zd = sum(probs.values()); probs = {O: v/Zd for O, v in probs.items()}
S64 = -sum(v*np.log(v) for v in probs.values()); Sadd = sum(h2(1/(np.exp(E*tau) + 1)) for E in Es)
best = max(probs, key=probs.get); sea = tuple(i for i in range(D) if Es[i] < 0)
ok("2. three pairs of slots (k = 0, 0.7, 2; a tau = 1.2), 64 fillings: sum of the principal minors = det(1 + U) = product over the pairs (%.0e, %.0e); "
   "the entropy of the 64 = the six slots' yes/no entropies added, %.4f (%.0e); the likeliest filling is 'every lagging slot, no leading one' with "
   "probability %.3f = product of 1/(1 + e^{-w tau})^2 (%.0e)"
   % (abs(Zs - np.linalg.det(np.eye(D) + U))/Zs, abs(Zs - Zd)/Zs, S64, abs(S64 - Sadd), probs[best], abs(probs[best] - np.prod([1/(1 + np.exp(-w*tau))**2 for w in ws]))),
   abs(Zs - np.linalg.det(np.eye(D) + U))/Zs < 1e-12 and abs(Zs - Zd)/Zs < 1e-12 and abs(S64 - Sadd) < 1e-12 and best == sea
   and abs(probs[best] - np.prod([1/(1 + np.exp(-w*tau))**2 for w in ws])) < 1e-12)

# =====================================================================================================================================
# 3. many pairs, no mass: per unit extent
# =====================================================================================================================================
res3 = []
for tau in (0.5, 1.0, 3.0):
    ext = 4000.0*tau; n = np.arange(-200000, 200001); kk = 2*np.pi*n/ext       # the k's of an extent 'ext' (period in t); spacing 2 pi/ext
    x = np.abs(kk)*tau; x = x[x < 60]
    S = float((2*(np.log1p(np.exp(-x)) + x/(np.exp(x) + 1))).sum())           # two slots per k; yes/no entropy of f = 1/(1 + e^-x) is log(1 + e^-x) + x/(e^x + 1)
    dF = float((-2/tau*np.log1p(np.exp(-x))).sum())                             # F - E_0
    res3.append((tau, S/ext*tau/(np.pi/3), dF/ext*tau**2/(-np.pi/6)))
ok("3. no mass, many k (two slots each): entropy of the count per unit extent = pi/(3 tau) (ratio %s at tau = 0.5, 1, 3) and F - E_0 per unit extent = "
   "-pi/(6 tau^2) (ratio %s): the numbers of one species with c = 1" % (["%.5f" % r[1] for r in res3], ["%.5f" % r[2] for r in res3]),
   all(abs(r[1] - 1) < 2e-3 and abs(r[2] - 1) < 2e-3 for r in res3))

# =====================================================================================================================================
# 4. the sense of one walker
# =====================================================================================================================================
aa = 1.0; ts = np.linspace(0, 6, 121)
H = np.array([h2((1 + np.exp(-2*aa*t))/2) for t in ts]); Drel = np.log(2) - H
ok("4a. one walker at rest, started in the forward sense, VI.7 (3): the entropy of its sense rises from %.4f to %.4f (log 2 = %.4f), every step up; "
   "its distance from equilibrium, sum P log(P/pi), is log 2 minus that and only falls" % (H[0], H[-1], np.log(2)),
   H[0] == 0 and abs(H[-1] - np.log(2)) < 1e-9 and np.all(np.diff(H) > 0) and np.all(np.diff(Drel) < 0))
worst_mono = True; eq = []
for eta in (0.5, 1.3, 2.0):
    gen = np.array([[-aa*np.exp(eta), aa*np.exp(-eta)], [aa*np.exp(eta), -aa*np.exp(-eta)]]); pi = np.array([np.exp(-eta), np.exp(eta)])/(2*np.cosh(eta))
    for start in (np.array([1., 0.]), np.array([0., 1.])):
        Dk = []
        for t in np.linspace(0, 4, 81):
            P = sl.expm(t*gen)@start; Dk.append(sum(P[i]*np.log(P[i]/pi[i]) for i in range(2) if P[i] > 1e-300))
        worst_mono = worst_mono and bool(np.all(np.diff(Dk) < 1e-13)) and Dk[-1] < 1e-3
    eq.append((eta, h2(pi[0]), h2((1 - np.tanh(eta))/2)))
ok("4b. in motion (III.15d (12), rapidity 0.5, 1.3, 2): the distance from equilibrium only falls, from either sense; at equilibrium the sense is not "
   "forgotten entirely - its entropy is h((1 - v/c)/2) = %s, below log 2, and 0 at the speed of light" % ["%.4f" % e[1] for e in eq],
   worst_mono and all(abs(e[1] - e[2]) < 1e-12 and e[1] < np.log(2) for e in eq))

# =====================================================================================================================================
# 5. the settled configuration seen from inside a region
# =====================================================================================================================================
def region_entropy(C, l):
    nu = np.linalg.eigvalsh(C[:2*l, :2*l]); nu = np.clip(nu, 1e-15, 1 - 1e-15)
    return float(-(nu*np.log(nu) + (1 - nu)*np.log(1 - nu)).sum())
def assemble(P, kk, N):
    F = np.exp(1j*np.outer(np.arange(N), kk))/np.sqrt(N); C = np.zeros((N, 2, N, 2), complex)
    for a_ in range(2):
        for b_ in range(2): C[:, a_, :, b_] = (F*P[:, a_, b_][None, :])@F.conj().T       # C[x a, y b] = sum_k F[x,k] P[k,a,b] conj(F[y,k])
    C = C.reshape(2*N, 2*N); return (C + C.conj().T)/2
def lines_C(N, a):
    """the continuum lines on a ring of N sites (N odd): for each k = 2 pi n/N, |n| <= (N-1)/2, the lower column of h(k) = k sigma_z - a sigma_y."""
    kk = 2*np.pi*np.arange(-(N - 1)//2, (N - 1)//2 + 1)/N; w = np.sqrt(kk*kk + a*a); P = np.zeros((N, 2, 2), complex)
    P[:, 0, 0] = (1 - kk/w)/2; P[:, 1, 1] = (1 + kk/w)/2; P[:, 0, 1] = -1j*a/(2*w); P[:, 1, 0] = 1j*a/(2*w)      # (1 - h/w)/2
    return assemble(P, kk, N)
def tick_C(L, th):
    """III.15a's tick on a ring: per k the 2 x 2 step; the columns with quasi-energy in (-pi, 0) are the positive part of (U - U^dagger)/2i."""
    kk = 2*np.pi*np.arange(L)/L; c, s = np.cos(th), np.sin(th); P = np.zeros((L, 2, 2), complex); gaps = []
    for i, k in enumerate(kk):
        Uk = np.array([[c*np.exp(-1j*k), s*np.exp(1j*k)], [-s*np.exp(-1j*k), c*np.exp(1j*k)]]); Hs = (Uk - Uk.conj().T)/2j
        e, v = np.linalg.eigh(Hs); P[i] = np.outer(v[:, 1], v[:, 1].conj()); gaps.append(-np.angle(np.linalg.eigvals(Uk)))
    return assemble(P, kk, L), np.array(gaps)

N = 601; C0 = lines_C(N, 1e-4)
pure = np.abs(C0@C0 - C0).max(); occ = np.abs(np.diag(C0).real - 0.5).max()
ls = np.array([4, 8, 16, 32, 64, 128, 200]); S0 = np.array([region_entropy(C0, l) for l in ls]); chord = np.log(N/np.pi*np.sin(np.pi*ls/N))
slope0 = np.polyfit(chord, S0, 1)[0]; comp = abs(region_entropy(C0, 200) - region_entropy(np.roll(np.roll(C0, -400, 0), -400, 1), 401))
ok("5a. continuum lines, 601 sites, almost no mass (a = 1e-4): Pi_- is a projector (%.0e) with every slot half filled (%.0e); the region and the rest "
   "have the same entropy (%.0e); for l = 4 ... 200 sites the entropy is %s, a straight line in log[(N/pi) sin(pi l/N)] of slope %.4f = 2/3: two edges "
   "in k, the one at k = 0 and the sharp cut at the largest k, one third each"
   % (pure, occ, comp, ["%.3f" % s for s in S0], slope0), pure < 1e-12 and occ < 1e-12 and comp < 1e-8 and abs(slope0 - 2/3) < 0.01)

Cm = lines_C(N, 0.05); lm = np.array([2, 4, 8, 64, 128]); Sm = np.array([region_entropy(Cm, l) for l in lm])
sl_short = (Sm[1] - Sm[0])/np.log(2); sl_long = (Sm[4] - Sm[3])/np.log(2)
am = np.array([0.02, 0.04, 0.08]); Sa = np.array([region_entropy(lines_C(N, x), 200) for x in am]); slope_a = np.polyfit(np.log(1/am), Sa, 1)[0]
ok("5b. with a mass (a = 0.05, Compton length 20 sites) the slope in log l is %.3f well inside the Compton length and %.3f well outside: the mass closes the "
   "edge at k = 0 and leaves the cut's one third. At l = 200 the entropy is %s for a = 0.02, 0.04, 0.08: slope %.3f in log(1/a) = 1/3. Through the edge "
   "at k = 0 the region shares (1/3) log(Compton length) with the rest, and no more however large it is" % (sl_short, sl_long, ["%.3f" % s for s in Sa], slope_a),
   sl_short > 0.6 and abs(sl_long - 1/3) < 0.04 and abs(slope_a - 1/3) < 0.02)

Lt = 1200; ths = (0.01, 0.02, 0.04); Ssat = []; flat = []; gap_ok = True
for th in ths:
    Ct, eps_k = tick_C(Lt, th)
    Ssat.append(region_entropy(Ct, 512)); flat.append(abs(region_entropy(Ct, 512) - region_entropy(Ct, 400)))
    near0 = np.abs(eps_k).min(); nearpi = (np.pi - np.abs(eps_k)).min(); gap_ok = gap_ok and abs(near0 - th) < 1e-9 and abs(nearpi - th) < 1e-9
slope_t = [(Ssat[i] - Ssat[i + 1])/np.log(2) for i in range(2)]
ok("5c. the lattice tick, 1200 sites, coin angle 0.01, 0.02, 0.04: its quasi-energies leave two gaps, around 0 and around pi, each of half-width the coin "
   "angle (to 1e-9) - the tick's sublattice of VI.4 (F); the region's entropy levels off (change from l = 400 to 512: %.0e) at %s, slope %.3f, %.3f in "
   "log(1/angle), towards 2/3: one third for each gap" % (max(flat), ["%.3f" % s for s in Ssat], slope_t[0], slope_t[1]),
   gap_ok and max(flat) < 1e-2 and all(abs(s - 2/3) < 0.03 for s in slope_t))

# =====================================================================================================================================
# 6. what the clock length is in the time reading, and the exchange of size and temperature
# =====================================================================================================================================
import sympy as sp, mpmath as mp
A4 = sp.diag(1, -1, 1, -1); S4 = sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]); T4 = A4*S4
J4 = sp.Matrix([[0, 0, -1, 0], [0, 0, 0, -1], [1, 0, 0, 0], [0, 1, 0, 0]]); D4 = sp.Matrix([[1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0]])
ksym, asym = sp.symbols('k a', real=True)
G4 = -ksym*A4*J4 + asym*T4; M4 = -ksym*A4 + asym*S4
exact = sp.simplify(D4.inv()*(J4*G4)*D4 + M4) == sp.zeros(4, 4) and D4*J4 == J4*D4 and sp.simplify(D4.inv()*G4*D4 - J4*M4) == sp.zeros(4, 4)
num = 0.0; zdev = 0.0
for (kv, av, beta) in ((0.0, 1.0, 0.7), (0.7, 1.0, 1.9), (2.0, 0.5, 3.0)):
    Gn = np.array(G4.subs({ksym: kv, asym: av}), float); Mn = np.array(M4.subs({ksym: kv, asym: av}), float); Jn = np.array(J4, float); Dn = np.array(D4, float)
    boltz = sl.expm(-beta*Jn@Gn); clock = Dn@sl.expm(beta*Mn)@np.linalg.inv(Dn)
    num = max(num, np.abs(boltz - clock).max()/np.abs(boltz).max())
    z2 = np.linalg.det(np.eye(2) + sl.expm(beta*Mk(kv, av)))
    zdev = max(zdev, abs(np.linalg.det(np.eye(4) + boltz) - z2**2)/z2**2)       # four real counters = the two slots, each written twice (route and twin)
ok("6a. the clock length IS the inverse temperature of the time reading: with III.15e's D (which commutes with J), D^-1 (J G(k)) D = -M_k exactly "
   "(III.15e (4) times J), so the Boltzmann factor e^{-beta J G(k)} of the mode k is D L_k(beta) D^-1, the walk's own count over a clock of length beta "
   "(three cases, relative %.0e), and det(1 + e^{-beta J G}) = det(1 + L_k(beta))^2 on the four real counters (%.0e). No continuation in i: the same k is "
   "the tilt in the clock and the wavenumber in t" % (num, zdev), exact and num < 1e-12 and zdev < 1e-12)

mp.mp.dps = 30
def chan(beta, L):                                                            # no mass: -beta E_0(L) + the rest, E_0(L) = -pi/(6L) (half-integer modes, +1/24)
    return beta*mp.pi/(6*L) + 4*mp.nsum(lambda m: mp.log(1 + mp.e**(-beta*(2*m + 1)*mp.pi/L)), [0, mp.inf])
pairs6 = [(mp.mpf(1), mp.mpf(2)), (mp.mpf('0.3'), mp.mpf('1.7')), (mp.mpf(5), mp.mpf('0.8'))]
dev6 = max(abs(chan(b, L) - chan(L, b))/abs(chan(b, L)) for b, L in pairs6)
half = mp.nsum(lambda n: (n + mp.mpf(1)/2)*mp.e**(-mp.mpf('0.001')*(n + mp.mpf(1)/2)), [0, mp.inf]) - 1/mp.mpf('0.001')**2
ok("6b. size and temperature exchange, no mass: a ring of length L at inverse temperature beta, and a clock of length L with the tilts k_n = (2n+1) pi/beta, "
   "give the same log Z: beta pi/(6L) + 4 sum_m log(1 + e^{-beta (2m+1) pi/L}) is unchanged when beta and L are exchanged (three pairs, relative %.0e). "
   "The first term is -beta E_0(L), E_0 = -pi/(6L), from the finite part +1/24 of the half-integer sum (here %s); the ring's settled energy in one reading is "
   "the other's thermal free energy -pi L/(6 beta^2)" % (float(dev6), mp.nstr(half, 8)), dev6 < mp.mpf(10)**-25 and abs(half - mp.mpf(1)/24) < 1e-6)

def ddiff(beta, L, a1, a2, Mmax=2_000_000):                                   # the difference between two masses converges (terms ~ 1/m^3): plain sum + tail
    pm = (2*np.arange(Mmax) + 1)*np.pi/L
    f = lambda a_: beta/np.sqrt(a_**2 + pm**2)*np.tanh(beta*np.sqrt(a_**2 + pm**2)/2)
    tail = -beta*(a1**2 - a2**2)/2*(L/np.pi)**3/(16*Mmax**2)                  # sum_{m >= M} 1/(2m+1)^3 = 1/(16 M^2) to leading order
    return float((f(a1) - f(a2)).sum() + tail)
pairs6f = [(float(b), float(L)) for b, L in pairs6]
dev6c = max(abs(ddiff(b, L, 0.5, 1.5) - ddiff(L, b, 0.5, 1.5))/abs(ddiff(b, L, 0.5, 1.5)) for b, L in pairs6f)
ok("6c. with mass: d log Z/d(a^2), taken between a = 0.5 and a = 1.5 so that it converges, is the same from sum_m (beta/2 eps_m) tanh(beta eps_m/2) over the ring's "
   "wavenumbers (2m+1) pi/L and from sum_n (L/2 w_n) tanh(L w_n/2) over the tilts (2n+1) pi/beta (three pairs, relative %.0e): one double sum over "
   "1/(a^2 + p_m^2 + k_n^2), done in either order" % float(dev6c), dev6c < 1e-9)

xs = [mp.mpf('0.4'), mp.mpf(2), mp.mpf(7)]
dev6d = max(abs(mp.log((1 + mp.e**x)*(1 + mp.e**(-x))) - 2*mp.log(2) - 2*mp.nsum(lambda n: mp.log(1 + x**2/((2*n + 1)*mp.pi)**2), [0, mp.inf])) for x in xs)
ok("6d. why the tilts are half-integers: (1 + e^x)(1 + e^-x) = 4 prod_n [1 + x^2/((2n+1) pi)^2] (to %.0e at x = 0.4, 2, 7). The + in 1 + e^x is VI.5 (A)'s "
   "exchange sign on a glued history; it puts the zeros at odd multiples of i pi, and those are the tilts of the partner reading" % float(dev6d), dev6d < mp.mpf(10)**-12)

print("ALL PASS" if ok_all else "SOME FAILED"); raise SystemExit(0 if ok_all else 1)
