"""
Creation and annihilation read off the walk (1+1 lattice, the t-reading of III.15a: shift, then the coin; U real orthogonal).

Nothing is added under the walk. A state at a time slice is a list of which columns hold a walker (exclusion, VI.1); the
amplitude between two such lists is VI.1's D - E generalized, a determinant of one-walker readings U_{XA}. The checks
below read the operators, their algebra, Feynman's rule, the vacuum energy and pair creation off that.

Checks:
  1. the determinant is what exclusion makes of the sum over zigzags: the antisymmetric two- and three-walker columns of
     VI.1, evolved tick by tick with one U per index, equal at every configuration the 2x2 / 3x3 determinant of U^t readings
  2. creation = adding a row to the determinant (the Laplace expansion), annihilation = contracting one index; from that
     alone  c_w c+_v + c+_v c_w = <w,v>  and  c+_v c+_v = 0 (exclusion), on random columns and random many-walker states
  3. the vacuum is the filled sea; Feynman's rule falls out of it:
        <sea + walker at X | U^t | sea + walker at A> / <sea|U^t|sea>  =  (U^t Pi_+)_{XA}        (the positive half, forward)
        <sea + walker in a sea column | U^t | sea + walker at A> = 0                             (blocked: exclusion)
        <hole at X | U^t | hole at A> / <sea|U^t|sea>  =  conj[(U^t Pi_-)_{XA}] = (Pi_- U^{-t})_{AX}  (a negative column from X back to A)
     the same ratios by brute force on 4 sites (the rank-5 antisymmetric tensor evolved index by index: the list of all
     configuration amplitudes, no determinant used)
  4. the vacuum-to-vacuum amplitude is e^{-i E_0 t}, E_0 = sum of the sea's quasi-energies = -(1/2) Tr|H|, and -(1/2)Tr|H| is a
     sum over closed loops: |eps| = pi/2 - (4/pi) sum_{k odd} cos(k eps)/k^2, so E_0 = -dim pi/4 + (2/pi) sum_{k odd} Tr U^k / k^2,
     Tr U^k the signed count of closed k-tick loops (the Casimir doc's "coin's version", line 2, now derived from the sea)
  5. pair creation is V.2's turn applied between configurations: a phase e^{i theta} on one site for one tick takes the sea to
     (sea with sea column o replaced by empty column e) with amplitude i theta N_eo e^{-i(eps_e-eps_o)}, N_eo = sum_c conj(e_c(x0)) o_c(x0): the
     "lines" of the plates runs; and sum over pairs |amp|^2 = 1 - |vacuum persistence|^2 to O(theta^3)
"""
import numpy as np, itertools, sys, scipy.linalg as sl
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all = ok_all and bool(cond); print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()
rng = np.random.default_rng(5)

def tick_matrix(L, th):
    dim = 2*L; U = np.zeros((dim, dim))
    R = np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]])
    for x in range(L):
        for s, y in ((0, (x + 1) % L), (1, (x - 1) % L)):
            U[2*y, 2*x + s] = R[0, s]; U[2*y + 1, 2*x + s] = R[1, s]
    assert np.allclose(U.T @ U, np.eye(dim)); return U
def wedge(cols):                                     # sum_sigma sgn(sigma) a_{sigma 1} x ... x a_{sigma n}
    n = len(cols); dim = len(cols[0]); T = np.zeros((dim,)*n, complex)
    for perm in itertools.permutations(range(n)):
        sgn = np.linalg.det(np.eye(n)[list(perm)])
        t = cols[perm[0]]
        for j in perm[1:]: t = np.multiply.outer(t, cols[j])
        T += sgn*t
    return T
def evolve(T, U, t):
    for _ in range(t):
        for ax in range(T.ndim): T = np.tensordot(U, T, axes=([1], [ax])); T = np.moveaxis(T, 0, ax)
    return T

# ---------------- 1. two and three walkers ----------------
L = 6; th = 0.3; U = tick_matrix(L, th); dim = 2*L; t = 7
f, g, h = (rng.normal(size=dim) for _ in range(3))
C2 = evolve(wedge([f, g]), U, t); C3 = evolve(wedge([f, g, h]), U, t)
Ut = np.linalg.matrix_power(U, t); F, G, H = Ut @ f, Ut @ g, Ut @ h
d2 = max(abs(C2[X, Y] - np.linalg.det(np.array([[F[X], G[X]], [F[Y], G[Y]]]))) for X in range(dim) for Y in range(dim))
d3 = max(abs(C3[X, Y, Z] - np.linalg.det(np.array([[F[X], G[X], H[X]], [F[Y], G[Y], H[Y]], [F[Z], G[Z], H[Z]]])))
         for X in range(dim) for Y in range(dim) for Z in range(dim))
ok("1. two and three walkers (12 columns, %d ticks): the antisymmetric columns evolved index by index equal, at every configuration, "
   "the determinant of one-walker readings, to %.1e and %.1e" % (t, d2, d3), d2 < 1e-12 and d3 < 1e-12)

# ---------------- 2. creation and annihilation from the determinant ----------------
def create(v, T):                                    # v wedge T: the Laplace expansion along the new row
    n = T.ndim; out = np.multiply.outer(v, T)
    for j in range(1, n + 1):                        # move the new index to position j with the sign of the transposition chain
        out = out + (-1)**j*np.moveaxis(np.multiply.outer(v, T), 0, j)
    return out
def annihilate(w, T):                                # contraction of the first index with <w|
    return np.tensordot(np.conj(w), T, axes=([0], [0]))
T2 = wedge([f, g]); v, w = rng.normal(size=dim) + 1j*rng.normal(size=dim), rng.normal(size=dim) + 1j*rng.normal(size=dim)
lap = max(abs(create(v, T2) - wedge([v, f, g])).max(), abs(create(v, create(w, T2)) + create(w, create(v, T2))).max())
car = abs(annihilate(w, create(v, T2)) + create(v, annihilate(w, T2)) - (np.conj(w) @ v)*T2).max()
exc = abs(create(v, create(v, T2))).max()
ok("2. creation = the Laplace expansion (v wedge T) to %.0e, anticommuting to %.0e; from it {c_w, c+_v} = <w,v> to %.1e and c+_v c+_v = 0 to %.0e: "
   "the operator algebra is the determinant's antisymmetry, i.e. exclusion in every state" % (lap, lap, car, exc), lap < 1e-12 and car < 1e-12 and exc < 1e-12)

# ---------------- 3. the sea, and Feynman's rule ----------------
def sea_of(U):
    Tm, Z = sl.schur(U.astype(complex), output='complex'); eps = -np.angle(np.diag(Tm))
    neg = eps < 0; pos = eps > 0
    return Z[:, neg], eps[neg], Z[:, pos] @ Z[:, pos].conj().T, Z[:, neg] @ Z[:, neg].conj().T
L = 40; U = tick_matrix(L, th); dim = 2*L; Phi, eps_neg, Pp, Pm = sea_of(U); N = Phi.shape[1]
t = 9; Ut = np.linalg.matrix_power(U, t)
def slater_amp(B, A): return np.linalg.det(B.conj().T @ Ut @ A)      # <b1^...^bn| U^t |a1^...^an>
den = slater_amp(Phi, Phi)
e = np.eye(dim); worst_p = 0.0; worst_b = 0.0
for A in rng.integers(0, dim, 6):
    for X in rng.integers(0, dim, 6):
        r = slater_amp(np.column_stack([Phi, e[:, X]]), np.column_stack([Phi, e[:, A]]))/den
        worst_p = max(worst_p, abs(r - (Ut @ Pp)[X, A]))
    worst_b = max(worst_b, abs(slater_amp(np.column_stack([Phi, Phi[:, 3]]), np.column_stack([Phi, e[:, A]]))))
# holes: c_{e_A} |sea> = sum_j (-1)^{j-1} <e_A, phi_j> |sea without j>
M = Phi.conj().T @ Ut @ Phi
def minor(M, i, j): return np.linalg.det(np.delete(np.delete(M, i, 0), j, 1))
minors = np.array([[minor(M, i, j) for j in range(N)] for i in range(N)])
worst_h = 0.0
for A in rng.integers(0, dim, 4):
    for X in rng.integers(0, dim, 4):
        amp = sum((-1)**(i + j)*np.conj(Phi[X, i])*Phi[A, j]*minors[i, j] for i in range(N) for j in range(N))/den
        worst_h = max(worst_h, abs(amp - np.conj((Ut @ Pm)[X, A])), abs(amp - (Pm @ np.linalg.matrix_power(U, -t))[A, X]))
# brute force on 4 sites: the rank-5 tensor is the list of all configuration amplitudes
Ls = 4; Us = tick_matrix(Ls, th); dims = 2*Ls; Phis, _, Pps, _ = sea_of(Us); Ns = Phis.shape[1]; ts = 5
Uts = np.linalg.matrix_power(Us, ts); es = np.eye(dims); worst_bf = 0.0
Tsea = evolve(wedge([Phis[:, j] for j in range(Ns)]), Us, ts)
den_bf = np.vdot(wedge([Phis[:, j] for j in range(Ns)]), Tsea)
for A in range(dims):
    T5 = evolve(wedge([Phis[:, j] for j in range(Ns)] + [es[:, A]]), Us, ts)
    for X in range(dims):
        num_bf = np.vdot(wedge([Phis[:, j] for j in range(Ns)] + [es[:, X]]), T5)
        worst_bf = max(worst_bf, abs(num_bf/den_bf/(Ns + 1) - (Uts @ Pps)[X, A]))     # (n+1)!/n! from the unnormalized wedges
ok("3. the sea (%d columns of 40 sites, %d ticks): <sea + walker at X|U^t|sea + walker at A>/<sea|U^t|sea> = (U^t Pi_+)_XA to %.1e; "
   "onto a sea column the amplitude is %.0e (blocked); <hole at X|U^t|hole at A>/<sea|U^t|sea> = conj[(U^t Pi_-)_XA] = (Pi_- U^-t)_AX to %.1e: "
   "the negative column runs from X back to A. Brute force on 4 sites (the rank-5 tensor evolved index by index, %d configurations, no determinant): "
   "the same ratio to %.1e. Feynman's rule is exclusion plus the filled sea"
   % (N, t, worst_p, worst_b, worst_h, int(np.math.comb(dims, Ns + 1)) if hasattr(np, 'math') else 56, worst_bf),
   worst_p < 1e-11 and worst_b < 1e-12 and worst_h < 1e-10 and worst_bf < 1e-10)

# ---------------- 4. the vacuum energy as closed loops ----------------
E0 = eps_neg.sum()
phase_dev = max(abs(np.linalg.det(Phi.conj().T @ np.linalg.matrix_power(U, tt) @ Phi) - np.exp(-1j*E0*tt)) for tt in (1, 5, 13))
Tm, Z = sl.schur(U.astype(complex), output='complex'); eps_all = -np.angle(np.diag(Tm))
half_trH = -0.5*np.abs(eps_all).sum()
K = 4001; trU = []; P = np.eye(dim)
for k in range(1, K + 1):
    P = P @ U; trU.append(np.trace(P))
trU = np.array(trU); ks = np.arange(1, K + 1); odd = ks % 2 == 1
loops = -dim*np.pi/4 + (2/np.pi)*(trU[odd]/ks[odd]**2).sum()
tail = (2/np.pi)*np.abs(trU).max()*(1/(2*K))                           # the omitted tail is below this
ok("4. <sea|U^t|sea> = e^{-i E_0 t} with E_0 = sum of the sea's quasi-energies = %.6f, to %.0e; E_0 = -(1/2)Tr|H| = %.6f; "
   "and -(1/2)Tr|H| as a sum over closed loops, -dim pi/4 + (2/pi) sum_{k odd} Tr U^k/k^2 with Tr U^k the signed count of closed k-tick loops, "
   "k <= %d: %.6f (the omitted tail is below %.1e). The vacuum's energy is the loop sum with weight 1/k^2, from the sea"
   % (E0, phase_dev, half_trH, K, loops, tail), phase_dev < 1e-10 and abs(E0 - half_trH) < 1e-12 and abs(loops - E0) < tail + 1e-6)
odd_max = np.abs(trU[odd]).max()
ok("4'. what check 4 tests of the loop sum: on this ring (an even number of sites) every Tr U^k with k odd is 0 (largest %.0e), so the loop sum is empty and "
   "E_0 = -dim pi/4 to %.0e. The loop term itself is tested on odd rings in review_v33_loop_sum_on_odd_ring.py" % (odd_max, abs(E0 + dim*np.pi/4)),
   odd_max < 1e-9 and abs(E0 + dim*np.pi/4) < 1e-10)

# ---------------- 5. pair creation by the turn ----------------
x0 = 7; theta = 0.02
V = np.eye(dim, dtype=complex); V[2*x0, 2*x0] = V[2*x0 + 1, 2*x0 + 1] = np.exp(1j*theta)      # V.2's turn at one site, one tick, on the complex column
UV = U @ V
Tm, Z = sl.schur(U.astype(complex), output='complex'); eps_all = -np.angle(np.diag(Tm)); negi = np.where(eps_all < 0)[0]; posi = np.where(eps_all > 0)[0]
Phi = Z[:, negi]; E = Z[:, posi]
den = np.linalg.det(Phi.conj().T @ UV @ Phi)
worst5 = 0.0; tot = 0.0; big5 = 0.0
for o in range(len(negi)):
    for ie in range(len(posi)):
        B = Phi.copy(); B[:, o] = E[:, ie]
        amp = np.linalg.det(B.conj().T @ UV @ Phi)/den
        N_eo = np.conj(E[2*x0:2*x0 + 2, ie]) @ Phi[2*x0:2*x0 + 2, o]
        pred = 1j*theta*N_eo*np.exp(-1j*(eps_all[posi[ie]] - eps_all[negi[o]]))      # first order; the tick is turn then step, so the pair carries one tick of its rate
        worst5 = max(worst5, abs(amp - pred)); tot += abs(amp)**2; big5 = max(big5, abs(amp))
persist = abs(np.linalg.det(Phi.conj().T @ UV @ Phi)/np.linalg.det(Phi.conj().T @ U @ Phi))**2
ok("5. V.2's turn, e^{i theta} on one site for one tick, theta = %.2f: the amplitude from the sea to (sea with column o replaced by empty column e) "
   "(relative to the vacuum's own amplitude) equals i theta N_eo e^{-i(eps_e - eps_o)}, N_eo = sum_c conj(e_c(x0)) o_c(x0), to %.1e against a largest amplitude of %.1e "
   "(relative %.0e, second order in theta): pair creation is the vertex between configurations, "
   "and these are the lines of the plates runs; sum over pairs |amp|^2 = %.6e against 1 - |vacuum persistence|^2 = %.6e (relative difference %.0e)"
   % (theta, worst5, big5, worst5/big5, tot, 1 - persist, abs(tot/(1 - persist) - 1)),
   worst5 < 3*theta*big5 and abs(tot/(1 - persist) - 1) < 3*theta**2)       # 10 Nov 2026: both tolerances were absolute and larger than the amplitudes themselves
print("ALL PASS" if ok_all else "SOME FAILED"); raise SystemExit(0 if ok_all else 1)
