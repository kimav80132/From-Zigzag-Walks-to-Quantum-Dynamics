"""
Rebuilt 9 Nov 2026 against master v34 (check 4: the widths of III.8c (E)'s table, in units of c/a, on a longer lattice).
Ord's entwined paths, implemented in this document's own walk (1+1): the dictionary, the lock, and what it draws.

Ord & Mann, 'Entwined paths, difference equations and the Dirac equation' (quant-ph/0208004), Eq. (5): four real densities
phi^1..phi^4 (two 'envelopes' = the forward strand and the return strand of one entwined path, each split by light-cone
direction), obeying
    d_t phi^1 =  c d_z phi^1 - a phi^1 - a phi^2       d_t phi^3 =  c d_z phi^3 - a phi^3 + a phi^4
    d_t phi^2 = -c d_z phi^2 - a phi^2 + a phi^1       d_t phi^4 = -c d_z phi^4 - a phi^4 - a phi^3
with the loss removed by u = e^{at} phi ('only cosmetic'); the colouring rule 'the colour changes at every second corner';
the two envelopes 'start at the origin moving in opposite directions' with 'opposite colouring'.

Checks (hbar = c = 1):
  1. dictionary: Ord's left envelope IS III.5a (H) with the loss put back, with (Q+, Q-) = (phi^2, phi^1); the right envelope
     is the same lines with (Q'+, Q'-) = (phi^4, -phi^3). Two copies of the walk's column, nothing else. (matrix identity)
  2. colouring: the running product of III.15a (A)'s flip signs (+sin into the forward sense, -sin into the backward)
     changes sign at every second flip: Ord's colour rule is III.5a (J)'s sign, read along a path. And the enumeration of
     all paths with that sign reproduces the lattice column (III.15a (C)1) - the path level of the dictionary.
  3. the lock: the two envelopes start at the origin in opposite senses with opposite charge: u0 = e+ delta, v0 = -e- delta
     = T u0. So the entwined source is z0 = u0 + i v0 = (1, -i) delta: the REST-FRAME positive-energy column at a point
     (IV.8 (C)'s e+), i.e. the second copy is the first one flip on, LOCALLY - the rest-frame rule J = T applied at k = 0.
  4. the lock against the true partner J = L Omega^{-1} (III.8c): |(T - J) u0| / |u0| for sources of width sigma;
     the negative-energy weight of z = U(t)(1,-i) f for the same sources: O((c k / a)^2) - zero at rest, large for a point.
  5. what the entwined pair draws: U(t)(1,-i) delta = the real evolution kernel of III.15 (first column) and its partner
     (second column), NOT the positive-energy (Feynman) propagator: compare with P+ applied to the same source.
"""
import sys
import numpy as np, sympy as sp
np.set_printoptions(linewidth=160)
fails = 0
def ok(msg, cond):
    global fails
    print(("PASS  " if cond else "FAIL  ") + msg); fails += (not cond)

# ---------- 1. the dictionary, as matrices ----------
c_, a_, kx = sp.symbols('c a k', real=True)
# Ord Eq (5) in Fourier (d_z -> i k), on the column (phi1, phi2, phi3, phi4): d_t phi = M_ord phi
M_ord = sp.Matrix([[ sp.I*c_*kx - a_, -a_, 0, 0],
                   [ a_, -sp.I*c_*kx - a_, 0, 0],
                   [ 0, 0,  sp.I*c_*kx - a_,  a_],
                   [ 0, 0, -a_, -sp.I*c_*kx - a_]])
# III.5a (H) with the loss put back, t-reading, on (Q+, Q-): d_t Q = [-c sz d_x + a T - a] Q
sz = sp.diag(1, -1); T = sp.Matrix([[0, 1], [-1, 0]])
M_doc = -c_*sz*(sp.I*kx) + a_*T - a_*sp.eye(2)
# relabelling: (Q+, Q-) = (phi2, phi1); (Q'+, Q'-) = (phi4, -phi3)
P = sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]])   # Q = P phi
M_two = sp.diag(M_doc, M_doc)
ok("1. Ord's Eq (5) = two copies of III.5a (H) with the loss, under (Q+,Q-)=(phi2,phi1), (Q'+,Q'-)=(phi4,-phi3): P M_ord P^-1 - diag(M,M) = %s"
   % sp.simplify(P*M_ord*P.inv() - M_two).tolist(), sp.simplify(P*M_ord*P.inv() - M_two) == sp.zeros(4, 4))

# ---------- the lattice tick of the walk, t-reading (as in second_copy_is_the_hilbert_transform.py) ----------
a = 0.5; dt = 0.3; Lx = 1024
R = np.array([[np.cos(a*dt), np.sin(a*dt)], [-np.sin(a*dt), np.cos(a*dt)]])   # cos 1 + sin T
def tick(col):                                                 # col[x, s]: + moves right, - moves left, then the coin
    out = np.empty_like(col); out[:, 0] = np.roll(col[:, 0], 1); out[:, 1] = np.roll(col[:, 1], -1); return out @ R.T
def evolve(col, n):
    for _ in range(n): col = tick(col)
    return col

# ---------- 2. the colouring rule is the running flip sign; enumeration = the column ----------
# III.15a (A): per tick, keep -> cos; flip into + -> +sin; flip into - -> -sin. A path's sign is the product of its flip signs.
signs = []
for nflips in range(1, 9):
    s = 1; seq = []
    sense = +1
    for f in range(nflips):
        sense = -sense; s *= (+1 if sense > 0 else -1); seq.append(s)
    signs = seq
ok("2a. the running product of flip signs along a path starting in the + sense: %s - it changes at every SECOND flip (Ord's colour rule)" % signs,
   signs == [-1, -1, 1, 1, -1, -1, 1, 1])
# enumeration on a small lattice: sum over all 2^n tick-choices with the signed weights = lattice column from e+ at the origin
n = 7; cs, sn = np.cos(a*dt), np.sin(a*dt)
col_enum = {}
for bits in range(2**n):
    x = 0; sense = +1; w = 1.0
    for i in range(n):
        x += sense                                             # advance: + moves right, - moves left (one site per tick)
        if (bits >> i) & 1:                                    # flip
            sense = -sense; w *= (sn if sense > 0 else -sn)    # +sin into forward, -sin into backward  (Ord: colour)
        else:
            w *= cs
    col_enum[(x, sense)] = col_enum.get((x, sense), 0.0) + w
col = np.zeros((Lx, 2)); col[0, 0] = 1.0; col = evolve(col, n)
dev = max(abs(col_enum.get((x, s), 0.0) - col[x % Lx, 0 if s > 0 else 1]) for x in range(-n, n+1) for s in (+1, -1))
ok("2b. enumeration of all %d signed paths (Ord's coloured counting) reproduces the lattice column tick by tick: max dev %.1e" % (2**n, dev), dev < 1e-12)

# ---------- 3. the lock: opposite senses, opposite charge, at the origin ----------
e_plus = np.array([1.0, 0.0]); e_minus = np.array([0.0, 1.0]); Tm = np.array([[0, 1], [-1, 0.]])
u0 = e_plus; v0 = -e_minus                                       # envelope 1 starts +, blue; envelope 2 starts -, red (-1)
ok("3. the entwined source: u0 = e+, v0 = -e- = T u0 (%s): the second envelope is the first one flip on, at the point; z0 = u0 + i v0 = (1, -i) = the rest-frame positive-energy column e+ of IV.8 (C)"
   % np.allclose(v0, Tm @ u0), np.allclose(v0, Tm @ u0))

# ---------- 4. the lock against J = L Omega^{-1}, by source width ----------
xs = np.arange(Lx) - Lx//2; k = 2*np.pi*np.fft.fftfreq(Lx, 1.0)
# generator per mode of the continuum lines L = -sz d_x + a T  (the lattice tick's own generator differs at O(dt); use the continuum J)
def J_apply(col):                                                 # J = L / Omega mode by mode, on a real column
    ck = np.fft.fft(col, axis=0); out = np.empty_like(ck)
    for i, kk in enumerate(k):
        Lk = np.array([[-1j*kk, a], [-a, 1j*kk]]); Om = np.sqrt(a*a + kk*kk); out[i] = (Lk/Om) @ ck[i]
    return np.fft.ifft(out, axis=0).real
def Pneg_weight(z):                                               # weight of the negative-energy part of a complex column
    ck = np.fft.fft(z, axis=0); w = 0.0; tot = 0.0
    for i, kk in enumerate(k):
        H = np.array([[kk, 1j*a], [-1j*a, -kk]]); ev, V = np.linalg.eigh(H); vm = V[:, 0]
        w += abs(np.vdot(vm, ck[i]))**2; tot += np.vdot(ck[i], ck[i]).real
    return w/tot
rows = []
lam = 1/a                                                        # the Compton length c/a, in sites
for sig_c in (0.0, 0.5, 1.5, 5.0, 15.0):                          # widths in units of c/a, as III.8c (E)'s table quotes them
    sig = sig_c*lam
    f = np.exp(-xs**2/(2*sig**2)) if sig > 0 else (xs == 0).astype(float)
    u0c = np.outer(f, e_plus); v_lock = np.outer(f, Tm @ e_plus); v_true = J_apply(u0c)
    mis = np.linalg.norm(v_lock - v_true)/np.linalg.norm(u0c)
    wneg = Pneg_weight(u0c + 1j*v_lock)
    rows.append((sig_c, mis, wneg))
print("      source width sigma, in units of c/a (= %.1f sites):  |T u0 - J u0|/|u0|   negative-energy weight of (1,-i) f" % lam)
for sig_c, mis, wneg in rows: print("      sigma = %5.1f                      %.3f                  %.5f" % (sig_c, mis, wneg))
r = {s_: (m_, w_) for s_, m_, w_ in rows}
ok("4. the lock v0 = T u0 against the true partner J u0, III.8c (E)'s table: point %.2f, %.2f; width 5: %.2f, %.3f; width 15: %.3f, %.4f "
   "(the text: 1.09, 0.30; 0.14, 0.005; 0.047, 0.0006). First order and second order in c/(a sigma): mismatch x sigma = %.3f, %.3f; "
   "weight x sigma^2 = %.3f, %.3f at widths 5 and 15 (1/sqrt2 = %.3f and 1/8 = %.3f for a Gaussian)"
   % (r[0.0][0], r[0.0][1], r[5.0][0], r[5.0][1], r[15.0][0], r[15.0][1], r[5.0][0]*5, r[15.0][0]*15, r[5.0][1]*25, r[15.0][1]*225, 2**-0.5, 0.125),
   abs(r[0.0][0] - 1.09) < 0.005 and abs(r[0.0][1] - 0.30) < 0.005 and abs(r[5.0][0] - 0.14) < 0.005 and abs(r[5.0][1] - 0.005) < 5e-4
   and abs(r[15.0][0] - 0.047) < 5e-4 and abs(r[15.0][1] - 0.0006) < 5e-5 and abs(r[15.0][0]*15/(r[5.0][0]*5) - 1) < 0.05
   and abs(r[15.0][1]*225/(r[5.0][1]*25) - 1) < 0.1)

# ---------- 5. what the entwined pair draws ----------
nt = 40
u_t = evolve(np.outer((xs == 0).astype(float), e_plus), nt); v_t = evolve(np.outer((xs == 0).astype(float), Tm @ e_plus), nt)
z_ord = u_t + 1j*v_t
# the positive-energy evolution of the same source (Feynman's rule: keep the positive-energy half, IV.10 bottom row)
z0 = np.outer((xs == 0).astype(float), np.array([1.0, -1j]))
ck = np.fft.fft(z0, axis=0); out = np.empty_like(ck)
for i, kk in enumerate(k):
    H = np.array([[kk, 1j*a], [-1j*a, -kk]]); ev, V = np.linalg.eigh(H); Pp = np.outer(V[:, 1], V[:, 1].conj()); out[i] = Pp @ ck[i]
zF0 = np.fft.ifft(out, axis=0)
# evolve zF0 with the same lattice tick (complex column = two real columns): the Feynman column at nt
zF_t = evolve(zF0.real, nt) + 1j*evolve(zF0.imag, nt)
rel = np.linalg.norm(z_ord - zF_t)/np.linalg.norm(z_ord)
# and: u_t is III.15's real evolution kernel (first column of U), v_t its partner column: check u_t == U e+ trivially; show the kernel's two legs
ok("5. the entwined pair from a point draws U(t)(1,-i)delta: the real evolution kernel of III.15 and its partner (both real kernels), not the positive-energy propagator: relative distance to Feynman's column %.2f (it is the negative-energy admixture of check 4 carried along)" % rel, rel > 0.2)
print("NOTE  in the walk's words: Ord's two envelopes are the out-leg and a second strand launched at the origin in the opposite sense with the opposite charge;")
print("      at the ensemble level that is two copies of the walk's column (IV.5a) with v0 = T u0 - the rest-frame rule J = T applied at a point.")
print("      It is the one-sign pair (u, J u) exactly at rest and to O((ck/a)^2) for a source wider than the Compton length; for a point source it is not.")
print("ALL PASS" if not fails else "%d FAILED" % fails); sys.exit(1 if fails else 0)
