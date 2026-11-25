"""
The ribbon's Z and the closed loop: exchange partners under exclusion (1+1 lattice walk, t-reading, 40 sites).

One walker added to the filled sea, a static field Phi(x) on (V.2's time-link angle e^{i q Phi(x)} per tick), t ticks.
The exact amplitude, normalized by the vacuum's, is a Schur complement:
    <sea + X | U_Phi^t | sea + A> / <sea | U_Phi^t | sea>
      =  (U_Phi^t)_{XA}                                                   the ribbon alone: one-walker Dirac in the field, Z's included
       - e_X^+ U_Phi^t Phi (Phi^+ U_Phi^t Phi)^{-1} Phi^+ U_Phi^t e_A      the exchange with the sea: the closed loop
Checks:
  1. the identity to roundoff; with the field off the exchange term is (U^t Pi_-)_{XA}: Feynman's blocked negative half
  2. with the field on, the ribbon alone (Z's included) lands on columns the sea occupies with weight w > 0 - exclusion violated
     in the final state; after the exchange term the amplitude onto every occupied column is 0 exactly. The loop is what
     exclusion demands of the Z, and its field-dependent part has the Z's size at every q.
  3. two different ribbons stay two: <sea + X + Y | U_Phi^t | sea + A + B> / <sea|U_Phi^t|sea> = det [[G_XA, G_XB], [G_YA, G_YB]]
     with G the sea-dressed one-ribbon amplitude of check 1. Their mutual exchange is VI.1's E; each one's exchange with the
     sea is its own loop; the two do not mix.
NOTE: what decides whether the loop is the Z's partner: whether the Z's forward branch is a walker on the slice.
"""
import numpy as np, scipy.linalg as sl, sys
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all = ok_all and bool(cond); print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()
def note(msg): print("NOTE  " + msg); sys.stdout.flush()
rng = np.random.default_rng(2)
L = 40; a = 0.3; dim = 2*L; t = 9
U = np.zeros((dim, dim)); R = np.array([[np.cos(a), np.sin(a)], [-np.sin(a), np.cos(a)]])
for x in range(L):
    for s, y in ((0, (x + 1) % L), (1, (x - 1) % L)):
        U[2*y, 2*x + s] = R[0, s]; U[2*y + 1, 2*x + s] = R[1, s]
Tm, Z = sl.schur(U.astype(complex), output='complex'); eps = -np.angle(np.diag(Tm))
Phi = Z[:, eps < 0]; Pm = Phi @ Phi.conj().T; Pp = np.eye(dim) - Pm; N = Phi.shape[1]
xs = np.arange(L); bump = np.exp(-((xs - 20)**2)/(2*4.0**2))
def UPhi_t(q):
    V = np.diag(np.repeat(np.exp(1j*q*bump), 2)); return np.linalg.matrix_power(U @ V, t)
def dressed(Ut):                                   # G_XA for all X, A: direct minus exchange with the sea
    M = Phi.conj().T @ Ut @ Phi
    return Ut - (Ut @ Phi) @ np.linalg.solve(M, Phi.conj().T @ Ut), Ut, (Ut @ Phi) @ np.linalg.solve(M, Phi.conj().T @ Ut)
e = np.eye(dim)
def slater(B, A, Ut): return np.linalg.det(B.conj().T @ Ut @ A)

# 1
worst1 = 0.0; off = 0.0
for q in (0.0, 0.4):
    Ut = UPhi_t(q); G, direct, exch = dressed(Ut); den = slater(Phi, Phi, Ut)
    for A in rng.integers(0, dim, 5):
        for X in rng.integers(0, dim, 5):
            r = slater(np.column_stack([Phi, e[:, X]]), np.column_stack([Phi, e[:, A]]), Ut)/den
            worst1 = max(worst1, abs(r - G[X, A]))
    if q == 0.0: off = np.abs(exch - np.linalg.matrix_power(U, t) @ Pm).max()
ok("1. the exact one-walker amplitude in the sea = (ribbon alone, Z's included) - (exchange with the sea), to %.1e with and without the field; "
   "field off, the exchange term is (U^t Pi_-)_XA to %.0e: Feynman's blocked negative half is the loop at zero field" % (worst1, off), worst1 < 1e-11 and off < 1e-12)

# 2
rows = []
for q in (0.1, 0.2, 0.4):
    Ut = UPhi_t(q); G, direct, exch = dressed(Ut)
    A = 2*20                                           # launch at the bump, forward sense; positive-energy part of it
    psi0 = Pp @ e[:, A]; psi0 /= np.linalg.norm(psi0)
    w_direct = np.linalg.norm(Phi.conj().T @ (direct @ psi0))**2     # weight the ribbon alone puts on occupied columns
    w_dressed = np.linalg.norm(Phi.conj().T @ (G @ psi0))            # after the exchange term
    # the Z's size: the part of the ribbon alone that passed through a negative column (first order in q): direct - its Pi_+ projected evolution
    z_part = np.linalg.norm(Pm @ direct @ psi0)
    loop_part = np.linalg.norm(exch @ psi0 - np.linalg.matrix_power(U, t) @ Pm @ psi0)   # the field-dependent part of the exchange
    rows.append((q, w_direct, w_dressed, z_part, loop_part))
ok("2. a positive walker launched into the field: the ribbon alone lands on occupied columns with weight %s for q = %s (its Z's, first order in q: %s); "
   "after the exchange term that weight is %s - exclusion holds in the final state exactly. The field-dependent part of the exchange term (the loop) is %s: "
   "the same size as the Z at every q (ratio %s) - the loop is the Z's partner, term for term"
   % (["%.2e" % r[1] for r in rows], [r[0] for r in rows], ["%.3f" % r[3] for r in rows], ["%.0e" % r[2] for r in rows], ["%.3f" % r[4] for r in rows],
      ["%.3f" % (r[4]/r[3]) for r in rows]),
   all(r[1] > 1e-5 and r[2] < 1e-12 and abs(r[4]/r[3] - 1) < 0.05 for r in rows))

# 3
Ut = UPhi_t(0.4); G, direct, exch = dressed(Ut); den = slater(Phi, Phi, Ut); worst3 = 0.0; mix = 0.0
for _ in range(6):
    A, B, X, Y = rng.integers(0, dim, 4)
    if A == B or X == Y: continue
    r2 = slater(np.column_stack([Phi, e[:, X], e[:, Y]]), np.column_stack([Phi, e[:, A], e[:, B]]), Ut)/den
    d2 = G[X, A]*G[Y, B] - G[X, B]*G[Y, A]
    worst3 = max(worst3, abs(r2 - d2))
    naive = direct[X, A]*direct[Y, B] - direct[X, B]*direct[Y, A]
    mix = max(mix, abs(r2 - naive))
ok("3. two different ribbons: the exact two-walker amplitude in the sea equals det[[G_XA, G_XB],[G_YA, G_YB]] of the sea-dressed one-ribbon amplitudes to %.1e. "
   "Their mutual exchange is VI.1's E (the off-diagonal product); each ribbon's loop is inside its own G. The undressed determinant (ribbons alone, Z's but no loops) "
   "differs by up to %.3f: the loops are not the other ribbon" % (worst3, mix), worst3 < 1e-10 and mix > 1e-3)
note("what the identity says: the loop is the exchange partner of the ribbon's Z - it carries the exchange sign and it is exactly what removes the Z's weight "
     "from occupied columns. It is the Z's partner only if the Z's forward branch is a walker on the slice that exclusion applies to (the sea reading); "
     "in a one-ribbon reading the Z is present and the loop is absent, and then the final state can sit in an occupied column. Two different ribbons "
     "remain two: their exchange is E, and neither is the other's loop.")
print("ALL PASS" if ok_all else "SOME FAILED"); raise SystemExit(0 if ok_all else 1)
