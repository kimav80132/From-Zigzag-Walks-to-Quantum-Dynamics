"""Step 1 of the correction pass, 11 Nov 2026 (exploration; not cited in the master): ONE WALKER, ONE FRAME - the gaps closed by machine.

The walker's state in 3+1 is eight real numbers.  Two namings of the same eight:
   lists:  r = (u, v), two lists of IV.4b's four counters, obeying IV.5a (2):  d_t r = -c AA_k d_k r + a TT r      (AA_k, TT = A_k, T on each list);
   frame:  M = m0 1 + m_yz R_yz + m_zx R_zx + m_xy R_xy + m_x A_x + m_y A_y + m_z A_z + m_w w,   r = (M u0, M v0),
           R_jk = A_j A_k the turn of the space plane (j, k); A_k the turn of the plane (t, k) (a run); w = A_x A_y A_z;
           reference r0 = (u0, v0) = (e1, T e1): the walker at rest launched in the forward counter.
B = R_yz is the turn of the reference's own space plane across its spin axis.  rev(M) = M with the six plane turns reversed; M~ = T M T^-1.

  0. The dictionary m -> r is one-to-one, and Phi^T Phi = 2: the record r.r is twice the sum of the squares of the eight frame numbers.
  1. The matrices of the frame form, each 8 x 8 with entries 0, +-1:  LL_k (a run from the left), BB (the turn B from the right), PP (the mirror:
     +1 on [1, rotations], -1 on [runs, w]).  LL_k symmetric, square 1, anticommuting; BB antisymmetric, square -1, commuting with the LL_k;
     PP anticommutes with the LL_k and commutes with BB.  The mass flip on the frame is NN = -BB PP.
  2. The dictionary carries one form into the other:  Phi LL_k = AA_k Phi,  Phi NN = TT Phi,  Phi BB = JJ Phi.
  3. THE WHOLE EQUATION:   d_t m = -c LL_k d_k m - a BB PP m   (i.e. d_t M = -c A_k d_k M - a M~ B)   against the eight lines: mode by mode
     (random wavevectors) and on a periodic 3-d lattice from random real data.
  4. The record: density r.r = 2 m.m = (1/2) Tr M^T M; flow r.AA_k r = 2 m.LL_k m.
  5. THE FIELD: V.2's turn of the pair is the frame turned in its own spin plane, from the right:
         d_t m = -c LL_k (d_k - (e/hbar) Apot_k BB) m - a BB PP m - (e Phi/hbar) BB m ;
     run on a lattice with potentials against the pair form; a uniform change of clock-zero, m -> e^{chi BB} m, shifts Phi by -(hbar/e) d_t chi;
     BB commutes with every term, so VI.2 (B)'s argument holds as it stands; +-BB are the only weight-keeping turns the step does not read.
  6. OPPOSITE CHARGE: K (the second list reversed) is M -> M A_z; it reverses BB and commutes with the free step, so it reverses the field.
  7. THE TWO CLASSES: at rest the leading class is the mirror-even frames [1, rotations] and the lagging class the mirror-odd ones [runs, w];
     the first turn as e^{-a t BB}, the second as e^{+a t BB}.  In motion: the leading class at wavevector k is the rest class run to that speed.
  8. THE EIGHTH NUMBER: M rev(M) = alpha + gamma w = sigma e^{w beta}; no plane turn (from the left) and no turn in the spin plane (from the
     right) changes it; (alpha, gamma) are the two forms of `walker_at_rest_needs_the_pair.py` 4: s = -2 alpha, p = +2 gamma.  Every frame with a
     rest frame is  sqrt(sigma) e^{w beta/2} R  with R a product of plane turns (six numbers).  beta = 0: leading; beta = pi: lagging; the mirror
     sends beta to -beta, K sends it to beta + pi.
  9. Nothing depends on the reference: for any other walker at rest the turn B is the turn of the plane across ITS spin axis.
 10. The 1+1 stand-in: the frames [1, R_yz, A_y, A_z] are the pairs living on counters 1 and 3, the eight lines close on them when nothing depends
     on x and y, and there they are III.15e (2)'s four lines with (u+, u-, v+, v-) = (u1, u3, v1, v3).
 11. The record is unique: of the quadratic densities on the eight numbers, three are conserved with a flux and one of them is unchanged by the turn
     in the spin plane - the sum of the squares.
"""
import sys, itertools
import numpy as np, scipy.linalg as sl
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all = ok_all and bool(cond); print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()
rng = np.random.default_rng(21)
Ax = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]], float)
Ay = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]], float)
Az = np.diag([1., -1, -1, 1]); T = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]], float); I4 = np.eye(4); Ti = np.linalg.inv(T)
Ryz, Rzx, Rxy = Ay@Az, Az@Ax, Ax@Ay; w = Ax@Ay@Az
names = ['1', 'R_yz', 'R_zx', 'R_xy', 'A_x', 'A_y', 'A_z', 'w']
basis = [I4, Ryz, Rzx, Rxy, Ax, Ay, Az, w]; Emat = np.array([m.ravel() for m in basis]).T                # 16 x 8
coef = lambda X: np.linalg.lstsq(Emat, X.ravel(), rcond=None)[0]
frame = lambda m: sum(c_*E_ for c_, E_ in zip(m, basis))
inspan = lambda X: np.abs(Emat@coef(X) - X.ravel()).max() < 1e-12
u0 = np.array([1., 0, 0, 0]); v0 = T@u0; r0 = np.concatenate([u0, v0]); B = Ryz
d = lambda X: np.kron(np.eye(2), X)                                            # the same 4 x 4 matrix on each list
AA = [d(Ax), d(Ay), d(Az)]; TT = d(T); JJ = np.block([[np.zeros((4, 4)), -I4], [I4, np.zeros((4, 4))]]); KK = np.diag([1., 1, 1, 1, -1, -1, -1, -1])
Phi = np.array([np.concatenate([E_@u0, E_@v0]) for E_ in basis]).T             # r = Phi m
opmat = lambda f: np.array([coef(f(E_)) for E_ in basis]).T                   # the 8 x 8 matrix of a linear map on frames
LL = [opmat(lambda M, A=A: A@M) for A in (Ax, Ay, Az)]; BB = opmat(lambda M: M@B); PP = opmat(lambda M: T@M@Ti); NN = opmat(lambda M: -(T@M@Ti)@B)
rev = lambda M: frame(coef(M)*np.array([1, -1, -1, -1, -1, -1, -1, 1]))
show = lambda X: "[" + "; ".join(" ".join("%2d" % int(round(x)) for x in row) for row in X) + "]"
is_int = lambda X: np.abs(X - np.round(X)).max() < 1e-12 and np.abs(X).max() < 1 + 1e-12

# ---- 0. the dictionary
closed = all(inspan(P_@Q_) for P_ in basis for Q_ in basis)
ok("0. the eight frame matrices [%s] are closed under products; the dictionary r = Phi m is one-to-one (det Phi = %.0f) and Phi^T Phi = 2 x 1 (%.0e): "
   "the record r.r is twice the sum of the squares of the eight frame numbers"
   % (", ".join(names), np.linalg.det(Phi), np.abs(Phi.T@Phi - 2*np.eye(8)).max()), closed and abs(np.linalg.det(Phi)) > 1 and np.abs(Phi.T@Phi - 2*np.eye(8)).max() < 1e-12)

# ---- 1. the matrices of the frame form
acom = lambda X, Y: X@Y + Y@X; com = lambda X, Y: X@Y - Y@X; Z8 = np.zeros((8, 8)); I8 = np.eye(8)
props = all(is_int(X) for X in LL + [BB, PP, NN]) \
    and all(np.allclose(L_, L_.T) and np.allclose(L_@L_, I8) for L_ in LL) and all(np.allclose(acom(LL[i], LL[j]), Z8) for i in range(3) for j in range(i)) \
    and np.allclose(BB, -BB.T) and np.allclose(BB@BB, -I8) and all(np.allclose(com(BB, L_), Z8) for L_ in LL) \
    and np.allclose(PP, np.diag([1, 1, 1, 1, -1, -1, -1, -1])) and all(np.allclose(acom(PP, L_), Z8) for L_ in LL) and np.allclose(com(PP, BB), Z8) \
    and np.allclose(NN, -BB@PP) and np.allclose(NN, -NN.T) and np.allclose(NN@NN, -I8) and all(np.allclose(acom(NN, L_), Z8) for L_ in LL)
print("      order of the eight numbers: (m0, m_yz, m_zx, m_xy, m_x, m_y, m_z, m_w)")
for nm, X in (("LL_x", LL[0]), ("LL_y", LL[1]), ("LL_z", LL[2]), ("BB  ", BB), ("PP  ", PP), ("NN = -BB PP", NN)): print("      %s = %s" % (nm, show(X)))
pairs = sorted(set(tuple(sorted((j, int(np.argmax(np.abs(BB[:, j])))))) for j in range(8)))
ok("1. the matrices of the frame form (printed above), all with entries 0, +-1: LL_k symmetric, square 1, anticommuting; BB antisymmetric, square -1, "
   "commuting with the LL_k; PP = diag(1,1,1,1,-1,-1,-1,-1) anticommuting with the LL_k and commuting with BB; the mass flip NN = -BB PP is "
   "antisymmetric, squares to -1 and anticommutes with the LL_k. BB turns four pairs of frame numbers into each other: %s"
   % ", ".join("(%s, %s)" % (names[i], names[j]) for i, j in pairs), props and len(pairs) == 4)

# ---- 2. the dictionary carries one form into the other
inter = max(max(np.abs(Phi@L_ - A_@Phi).max() for L_, A_ in zip(LL, AA)), np.abs(Phi@NN - TT@Phi).max(), np.abs(Phi@BB - JJ@Phi).max())
ok("2. Phi LL_k = AA_k Phi, Phi NN = TT Phi, Phi BB = JJ Phi (%.0e): motion is a run from the left, the mass flip is the mirror frame turned in its "
   "spin plane, the turn between the two lists is the frame turned in its spin plane" % inter, inter < 1e-12)

# ---- 3. the whole equation
c_, a_ = 1.0, 0.8
gen_m = lambda kv: -1j*c_*sum(k_*L_ for k_, L_ in zip(kv, LL)) + a_*NN
gen_r = lambda kv: -1j*c_*sum(k_*A_ for k_, A_ in zip(kv, AA)) + a_*TT
dev_modes = max(np.abs(Phi@sl.expm(1.3*gen_m(kv)) - sl.expm(1.3*gen_r(kv))@Phi).max() for kv in rng.normal(size=(40, 3)))
n = 5; kgrid = 2*np.pi*np.fft.fftfreq(n, d=1.0)
m_start = rng.normal(size=(8, n, n, n)); r_start = np.einsum('ij,jxyz->ixyz', Phi, m_start)
def evolve(field, gen, tt):
    fh = np.fft.fftn(field, axes=(1, 2, 3)); out = np.zeros_like(fh)
    for ix, iy, iz in itertools.product(range(n), repeat=3):
        out[:, ix, iy, iz] = sl.expm(tt*gen((kgrid[ix], kgrid[iy], kgrid[iz])))@fh[:, ix, iy, iz]
    res = np.fft.ifftn(out, axes=(1, 2, 3)); return res.real, np.abs(res.imag).max()
m_t, im1 = evolve(m_start, gen_m, 1.3); r_t, im2 = evolve(r_start, gen_r, 1.3)
dev_lat = np.abs(np.einsum('ij,jxyz->ixyz', Phi, m_t) - r_t).max(); wt = abs((r_t**2).sum() - (r_start**2).sum())/(r_start**2).sum()
ok("3. THE WHOLE EQUATION in frame form, d_t m = -c LL_k d_k m - a BB PP m (that is, d_t M = -c A_k d_k M - a M~ B), against the eight lines of "
   "IV.5a (2): 40 random wavevectors, t = 1.3: %.0e; a periodic %d x %d x %d lattice from random real data: the two evolutions agree to %.0e at every "
   "site and counter, both stay real (%.0e), the weight is kept (%.0e)" % (dev_modes, n, n, n, dev_lat, max(im1, im2), wt),
   dev_modes < 1e-12 and dev_lat < 1e-12 and max(im1, im2) < 1e-12 and wt < 1e-12)

# ---- 4. the record
dev_rec = 0.0
for _ in range(100):
    m = rng.normal(size=8); r = Phi@m; M = frame(m)
    dev_rec = max(dev_rec, abs(r@r - 2*m@m), abs(r@r - 0.5*np.trace(M.T@M)), max(abs(r@A_@r - 2*m@L_@m) for A_, L_ in zip(AA, LL)))
ok("4. the record on the frame: density r.r = 2 m.m = (1/2) Tr M^T M, flow r.AA_k r = 2 m.LL_k m (%.0e)" % dev_rec, dev_rec < 1e-12)

# ---- 5. the field
ns = 24; Dm = (np.roll(np.eye(ns), 1, axis=1) - np.roll(np.eye(ns), -1, axis=1))/2.0                     # central difference on a ring (antisymmetric)
zz = np.arange(ns); pot = 0.6*np.cos(2*np.pi*zz/ns) + 0.2*np.sin(4*np.pi*zz/ns); apot = 0.4*np.sin(2*np.pi*zz/ns); eh = 1.0
def big(Lz, Nm, Bm, phi, av, shift=0.0):
    return (-c_*np.kron(Lz, Dm) + a_*np.kron(Nm, np.eye(ns)) - eh*np.kron(Bm, np.diag(phi + shift)) + eh*c_*np.kron(Lz@Bm, np.diag(av)))
Gm = big(LL[2], NN, BB, pot, apot); Gr = big(AA[2], TT, JJ, pot, apot); PhiBig = np.kron(Phi, np.eye(ns))
dev_field = np.abs(PhiBig@sl.expm(0.9*Gm) - sl.expm(0.9*Gr)@PhiBig).max(); antis = max(np.abs(Gm + Gm.T).max(), np.abs(Gr + Gr.T).max())
m0 = rng.normal(size=8*ns); nu = 0.37; tt = 0.9
gauge = {sg: np.abs(np.kron(sl.expm(sg*nu*tt*BB), np.eye(ns))@(sl.expm(tt*Gm)@m0) - sl.expm(tt*big(LL[2], NN, BB, pot, apot, shift=nu))@m0).max() for sg in (1, -1)}
sgn = min(gauge, key=gauge.get)
rows = []
for i in range(8):
    for j in range(8):
        E_ = np.zeros((8, 8)); E_[i, j] = 1.0; rows.append(np.concatenate([com(E_, G_).ravel() for G_ in LL + [NN]]))
comm_dim = 64 - np.linalg.matrix_rank(np.array(rows).T, tol=1e-9)
KKm = opmat(lambda M: M@Az)
in_comm = all(np.allclose(com(X, G_), Z8) for X in (BB, KKm, BB@KKm) for G_ in LL + [NN])
ok("5. THE FIELD on the frame: d_t m = -c LL_k (d_k - (e/hbar) Apot_k BB) m - a BB PP m - (e Phi/hbar) BB m. On a ring of %d sites with a scalar and a "
   "vector potential the frame form and the pair form (V.2 (E)-(F) with the turn JJ) agree to %.0e; both generators are antisymmetric (%.0e). Turning "
   "every frame in its spin plane by an angle that grows as nu t, m -> e^{%+d nu t BB} m, is the same as adding hbar nu/e to Phi (%.0e): the clock-zero "
   "of VI.2 (B) is the direction chosen in the spin plane. The matrices commuting with the whole free step on m are %d: 1, BB, K, BB K; so +-BB are "
   "the only weight-keeping turns the step does not read" % (ns, dev_field, antis, sgn, gauge[sgn], comm_dim),
   dev_field < 1e-11 and antis < 1e-12 and gauge[sgn] < 1e-11 and gauge[-sgn] > 1e-3 and comm_dim == 4 and in_comm)

# ---- 6. opposite charge
c6 = np.abs(Phi@KKm - KK@Phi).max(); rev_field = np.abs(np.kron(KKm, np.eye(ns))@Gm@np.kron(np.linalg.inv(KKm), np.eye(ns)) - big(LL[2], NN, BB, -pot, -apot)).max()
ok("6. OPPOSITE CHARGE: reversing the second list (VI.7's K) is, on the frame, M -> M A_z (%.0e); it reverses the spin-plane turn, K BB = -BB K (%s), "
   "and commutes with the free step, so the frame M A_z in the field (Phi, Apot) obeys the equation of M in the field (-Phi, -Apot) (%.0e)"
   % (c6, np.allclose(acom(KKm, BB), Z8), rev_field), c6 < 1e-12 and np.allclose(acom(KKm, BB), Z8) and rev_field < 1e-12)

# ---- 7. the two classes
Plead = (I8 + JJ@TT)/2                                                        # TT r = -JJ r : the second list a quarter turn ahead (III.15e (3))
c7a = np.abs(Phi@((I8 + PP)/2)@np.linalg.inv(Phi) - Plead).max()
turn = max(np.abs(sl.expm(a_*t_*NN)@((I8 + s_*PP)/2) - sl.expm(-s_*a_*t_*BB)@((I8 + s_*PP)/2)).max() for t_ in (0.4, 1.7) for s_ in (1, -1))
dev7 = 0.0
for _ in range(40):
    nv = rng.normal(size=3); nv /= np.linalg.norm(nv); eta = rng.uniform(0.1, 1.6); kk = a_*np.sinh(eta); ww = a_*np.cosh(eta)
    Q = frame(np.concatenate([rng.normal(size=4), np.zeros(4)]))                # a rest frame of the leading class: a size and rotations
    r = np.concatenate([sl.expm(eta*sum(x_*A_ for x_, A_ in zip(nv, (Ax, Ay, Az)))/2)@Q@u0, sl.expm(eta*sum(x_*A_ for x_, A_ in zip(nv, (Ax, Ay, Az)))/2)@Q@v0])
    G = -kk*sum(x_*A_ for x_, A_ in zip(nv, AA))@JJ + a_*TT; dev7 = max(dev7, np.abs(((I8 + JJ@G/ww)/2)@r - r).max())
ok("7. THE TWO CLASSES: at rest the leading class (second list a quarter turn ahead) is exactly the mirror-even frames, a size and the three rotations "
   "(%.0e), and the lagging class the mirror-odd ones, the three runs and w; the first turn as e^{-a t BB}, the second as e^{+a t BB} (%.0e): the sign of "
   "the energy is the sense in which the frame turns in its spin plane. In motion: the rest class run to rapidity eta in any direction is III.15e's "
   "leading class at ck = a sinh(eta) along that direction (40 random, %.0e)" % (c7a, turn, dev7), c7a < 1e-12 and turn < 1e-12 and dev7 < 1e-12)

# ---- 8. the eighth number
def sp_forms(m):
    r = Phi@m; u, v = r[:4], r[4:]; return 2*u@T@v, 2*u@T@w@v
dev8 = 0.0; inv_left = 0.0; inv_right = 0.0; canon = 0.0; mir = 0.0; kk8 = 0.0
gens6 = [Ryz, Rzx, Rxy, Ax, Ay, Az]
for _ in range(100):
    m = rng.normal(size=8); M = frame(m); N = M@rev(M); cN = coef(N); al, ga = cN[0], cN[7]
    dev8 = max(dev8, np.abs(np.delete(cN, [0, 7])).max(), abs(sp_forms(m)[0] + 2*al), abs(sp_forms(m)[1] - 2*ga))
    Lrot = I4
    for _ in range(3): Lrot = Lrot@sl.expm(sum(x_*g_ for x_, g_ in zip(rng.normal(size=6)*0.6, gens6))/2)
    cL = coef(Lrot@M@rev(Lrot@M)); inv_left = max(inv_left, abs(cL[0] - al), abs(cL[7] - ga))
    cR = coef(M@sl.expm(0.7*B)@rev(M@sl.expm(0.7*B))); inv_right = max(inv_right, abs(cR[0] - al), abs(cR[7] - ga))
    sig = np.hypot(al, ga); beta = np.arctan2(ga, al)
    R = (np.cos(beta/2)*I4 - np.sin(beta/2)*w)@M/np.sqrt(sig)                  # take out sqrt(sigma) e^{w beta/2}
    canon = max(canon, np.abs(rev(R)@R - I4).max(), np.abs(np.sqrt(sig)*(np.cos(beta/2)*I4 + np.sin(beta/2)*w)@R - M).max())
    cM = coef(T@M@Ti@rev(T@M@Ti)); mir = max(mir, abs(cM[0] - al), abs(cM[7] + ga))
    cK = coef(M@Az@rev(M@Az)); kk8 = max(kk8, abs(cK[0] + al), abs(cK[7] + ga))
lead = coef(I4@rev(I4)); lag = coef(Az@rev(Az))
ok("8. THE EIGHTH NUMBER, 100 random frames: M rev(M) = alpha + gamma w and nothing else (%.0e), with alpha, gamma the two forms that decide whether the "
   "walker has a rest frame, 2 u.Tv = -2 alpha, 2 u.(Tw)v = +2 gamma. No plane turn applied to the frame changes them (%.0e), nor the turn in its spin "
   "plane (%.0e). Writing alpha + gamma w = sigma e^{w beta}: M = sqrt(sigma) e^{w beta/2} R with rev(R) R = 1, a product of plane turns (%.0e): the eight "
   "numbers are six plane turns, a size sigma and an angle beta. The walker at rest: leading beta = 0, lagging beta = pi (alpha = %+.0f, %+.0f); the "
   "mirror sends beta to -beta (%.0e); K sends it to beta + pi (%.0e)" % (dev8, inv_left, inv_right, canon, lead[0], lag[0], mir, kk8),
   max(dev8, inv_left, inv_right, canon, mir, kk8) < 1e-10 and abs(lead[0] - 1) < 1e-12 and abs(lag[0] + 1) < 1e-12)

# ---- 9. nothing depends on the reference
dev9 = 0.0; par = 1.0
Sig = [-JJ@d(R_) for R_ in (Ryz, Rzx, Rxy)]                                    # Sigma_x, Sigma_y, Sigma_z of VI.7 (8)
for _ in range(60):
    qv = rng.normal(size=4); qv /= np.linalg.norm(qv); Q = frame(np.concatenate([qv, np.zeros(4)]))      # another walker at rest (leading class)
    rq = np.concatenate([Q@u0, Q@v0]); Bq = Q@B@np.linalg.inv(Q)
    dev9 = max(dev9, np.abs(JJ@rq - np.concatenate([Bq@Q@u0, Bq@Q@v0])).max(), np.abs(np.concatenate([Q@u0, Q@v0])@AA[0]@rq) + abs(rq@AA[1]@rq) + abs(rq@AA[2]@rq))
    bq = coef(Bq)[1:4]; spin = np.array([rq@S_@rq for S_ in Sig]); par = min(par, abs(bq@spin)/(np.linalg.norm(bq)*np.linalg.norm(spin)))
    dev9 = max(dev9, np.abs(np.delete(coef(Bq), [1, 2, 3])).max())
ok("9. nothing depends on the reference: for 60 other walkers at rest (the reference turned by a random rotation) the turn the step does not read is "
   "again one space-plane turn applied to that walker, its flow is zero (%.0e), and the plane is the one across that walker's own spin axis (cosine "
   "between the plane's axis and the spin: %.12f)" % (dev9, par), dev9 < 1e-12 and par > 1 - 1e-12)

# ---- 10. the 1+1 stand-in
sub = [0, 1, 5, 6]                                                             # [1, R_yz, A_y, A_z]
cnt = [0, 2, 4, 6]                                                             # counters u1, u3, v1, v3
others = [i for i in range(8) if i not in sub]; ocnt = [i for i in range(8) if i not in cnt]
on13 = np.abs(Phi[np.ix_(ocnt, sub)]).max(); Ph4 = Phi[np.ix_(cnt, sub)]
closes = max(np.abs(X[np.ix_(others, sub)]).max() for X in (LL[2], NN, BB))
A11 = np.diag([1., -1, 1, -1]); S11 = np.kron(np.eye(2), np.array([[0, 1.], [1, 0]])); T11 = A11@S11
J11 = np.block([[np.zeros((2, 2)), -np.eye(2)], [np.eye(2), np.zeros((2, 2))]])
same = max(np.abs(Ph4@LL[2][np.ix_(sub, sub)] - A11@Ph4).max(), np.abs(Ph4@NN[np.ix_(sub, sub)] - T11@Ph4).max(), np.abs(Ph4@BB[np.ix_(sub, sub)] - J11@Ph4).max())
ok("10. the 1+1 stand-in: the frames [1, R_yz, A_y, A_z] are exactly the pairs that live on counters 1 and 3 (%.0e); with nothing depending on x and y "
   "the run along z, the mass flip and the spin-plane turn keep them (%.0e); and on (u+, u-, v+, v-) = (u1, u3, v1, v3) they are III.15e's "
   "direction matrix, turn matrix and J (%.0e): the 1+1 two-list walk is the 3+1 walker moving along one axis, its phase the turn in the plane (y, z)"
   % (on13, closes, same), on13 < 1e-12 and closes < 1e-12 and same < 1e-12 and abs(np.linalg.det(Ph4)) > 0.5)

# ---- 11. the record on the eight numbers is the only one
import sympy as sp
qs = sp.symbols('q0:36'); idx = [(i, j) for i in range(8) for j in range(i, 8)]; Q = sp.zeros(8, 8)
for s_, (i, j) in zip(qs, idx): Q[i, j] = s_; Q[j, i] = s_
toS = lambda X: sp.Matrix(np.round(X).astype(int).tolist())
def nfree(conds):
    sol_ = sp.solve([e_ for c__ in conds for e_ in c__], qs, dict=True)[0]; Qs = Q.subs(sol_); return Qs, len(Qs.free_symbols)
cons = [list(Q*toS(L_) - toS(L_)*Q) for L_ in LL] + [list(Q*toS(NN) - toS(NN)*Q)]            # conserved with a flux: Q commutes with LL_k and with NN
turn = [list(toS(BB).T*Q*toS(BB) - Q)]                                                      # unchanged by the turn in the spin plane
Qc, nc = nfree(cons); Qb, nb = nfree(cons + turn)
ok("11. the record on the eight numbers: the quadratic densities m.Q m that are conserved with a flux are %d (u.u, v.v, u.v in the list naming); those "
   "also unchanged by the turn in the spin plane: %d, Q = the unit matrix - the sum of the squares of the eight numbers and nothing else"
   % (nc, nb), nc == 3 and nb == 1 and sp.simplify(Qb/Qb[0, 0]) == sp.eye(8))

print("ALL PASS" if ok_all else "SOME FAILED"); raise SystemExit(0 if ok_all else 1)
