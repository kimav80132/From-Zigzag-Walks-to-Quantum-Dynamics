"""IV.4c (G): the identities every solution of the two runners' conditions obeys, derived.   Recreated 9 Nov 2026 against master v34.

Conditions (IV.4c (C), (F)), with A_z the advance (A_z^2 = 1), T the mass turn (T^2 = -1, {A_z, T} = 0) and A_z T the exchange:
   (C1)  W_x = [K_x, A_z]     (C2)  [K_x, W_x] = A_z     (C3)  [K_x, A_z T] = 0     (C4)  [K_x, W_y] = 0       and the same with x <-> y.
1. The operator identities, symbolically, for non-commuting symbols (no matrices assumed):
     a. {A, [K, A]} = [K, A^2]                         so (C1) and A_z^2 = 1 give {A_z, W_x} = 0;
     b. [K, {A, W}] = {[K, A], W} + {A, [K, W]}        so commuting a. with K_x and using (C1), (C2):  0 = 2 W_x^2 + 2 A_z^2,  W_x^2 = -1;
        the same line with W -> W_y and (C4):           0 = {W_x, W_y};
     c. with A_x = -A_z W_x, A_y = -A_z W_y:  A_x^2 = 1, {A_x, A_z} = 0, {A_x, A_y} = 0, W_y W_x = A_x A_y  (reduced with a.-b.).
2. On solutions found from random starts (A_z and T fixed to IV.4b's): everything in 1, and the rest of (G)'s argument:
     d. K_x - A_x/2 commutes with A_x, A_y, A_z, so it is alpha 1 + beta omega, omega = A_z A_x A_y (the commutant has dimension 2);
     e. (C3) gives BOTH {A_x, T} = 0 and beta = 0;
     f. 2 [K_x, K_y] = A_x A_y = W_y W_x = -W_x W_y.
3. What v34's text says at this point and what is the case: the text has "(C3) ... restricts nothing else, since on these counters
   A_z A_x A_y = -+ A_z T".  On IV.4b's matrices A_z A_x A_y and A_z T are different matrices and ANTIcommute; A_z T is symmetric with
   square +1, omega is antisymmetric with square -1.  (C3) does restrict: it removes the omega part.  The four-line proof:
     T has no part commuting with omega (such a part would lie in the algebra of A_x, A_y, A_z, anticommute with A_z, and by the
     even half of (C3) for x and for y anticommute with A_x and A_y: nothing does), so T omega = -omega T; the odd half of (C3) is
     (1/2){A_x, T} = 2 beta omega T; multiply by T^{-1}: (1/2)(A_x + T A_x T^{-1}) = 2 beta omega; conjugate by T: the left side is
     unchanged, the right side changes sign; so beta = 0 and {A_x, T} = 0.
   (G)'s conclusions do not depend on this: f. holds for any beta."""
import sys
import numpy as np, sympy as sp
from scipy.optimize import least_squares
fails = 0
def ok(name, cond):
    global fails
    print(("PASS  " if cond else "FAIL  ") + name); fails += (not cond)
comm = lambda X, Y: X*Y - Y*X
acom = lambda X, Y: X*Y + Y*X

# ---- 1. identities in the free algebra
K, A, W, Wy = sp.symbols('K A W W_y', commutative=False)
ok("1a. {A, [K, A]} = [K, A^2] for non-commuting symbols", sp.expand(acom(A, comm(K, A)) - comm(K, A*A)) == 0)
ok("1b. [K, {A, W}] = {[K, A], W} + {A, [K, W]} (the commutator is a derivation of the anticommutator)",
   sp.expand(comm(K, acom(A, W)) - acom(comm(K, A), W) - acom(A, comm(K, W))) == 0)
# the consequences, with the relations put in by hand as substitutions on words: A^2 = 1, W = [K, A], [K, W] = A, [K, Wy] = 0
Wx = comm(K, A)
step_b = sp.expand(acom(Wx, Wx) + acom(A, A))                       # 1b with W = W_x:  {[K,A], W_x} + {A, [K,W_x]} -> 2 W_x^2 + 2 A^2
ok("1b'. with (C1) W_x = [K, A_z] and (C2) [K, W_x] = A_z, line 1b reads 0 = 2 W_x^2 + 2 A_z^2: W_x^2 = -1 once A_z^2 = 1",
   sp.expand(step_b - 2*Wx*Wx - 2*A*A) == 0)
# 1c on generic matrices obeying the derived relations is checked in 2 (it needs the relations, not free symbols)

# ---- 2. solutions from random starts
Az = np.diag([1., -1., -1., 1.]); T = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0.]]); Ex = Az @ T; I4 = np.eye(4)
def conds(v):
    Kx, Ky = v[:16].reshape(4, 4), v[16:].reshape(4, 4); Wx_, Wy_ = Kx@Az - Az@Kx, Ky@Az - Az@Ky
    r = [Kx@Wx_ - Wx_@Kx - Az, Kx@Ex - Ex@Kx, Kx@Wy_ - Wy_@Kx, Ky@Wy_ - Wy_@Ky - Az, Ky@Ex - Ex@Ky, Ky@Wx_ - Wx_@Ky]
    return np.concatenate([m.ravel() for m in r])
rng = np.random.default_rng(5); sols = []
for _ in range(240):
    s = least_squares(conds, rng.normal(size=32), xtol=1e-14, ftol=1e-14, gtol=1e-14)
    if np.linalg.norm(conds(s.x)) < 1e-9: sols.append(s.x)
mx = lambda M: np.abs(M).max()
w = dict(a=0, b=0, b2=0, c1=0, c2=0, c3=0, c4=0, d1=0, d2=0, e1=0, e2=0, f=0, asym=0, om=0)
for v in sols:
    Kx, Ky = v[:16].reshape(4, 4), v[16:].reshape(4, 4); Wx_, Wy_ = Kx@Az - Az@Kx, Ky@Az - Az@Ky
    Ax, Ay = -Az@Wx_, -Az@Wy_; om = Az@Ax@Ay
    w['a'] = max(w['a'], mx(Az@Wx_ + Wx_@Az), mx(Az@Wy_ + Wy_@Az))
    w['b'] = max(w['b'], mx(Wx_@Wx_ + I4), mx(Wy_@Wy_ + I4)); w['b2'] = max(w['b2'], mx(Wx_@Wy_ + Wy_@Wx_))
    w['c1'] = max(w['c1'], mx(Ax@Ax - I4), mx(Ay@Ay - I4)); w['c2'] = max(w['c2'], mx(Ax@Az + Az@Ax), mx(Ay@Az + Az@Ay))
    w['c3'] = max(w['c3'], mx(Ax@Ay + Ay@Ax)); w['c4'] = max(w['c4'], mx(Wy_@Wx_ - Ax@Ay))
    D = Kx - Ax/2
    w['d1'] = max(w['d1'], mx(D@Ax - Ax@D), mx(D@Ay - Ay@D), mx(D@Az - Az@D))
    al = np.trace(D)/4; be = np.trace(D@om.T)/np.trace(om@om.T)
    w['d2'] = max(w['d2'], mx(D - al*I4 - be*om)); w['e1'] = max(w['e1'], mx(Ax@T + T@Ax), mx(Ay@T + T@Ay)); w['e2'] = max(w['e2'], abs(be))
    w['f'] = max(w['f'], mx(2*(Kx@Ky - Ky@Kx) - Ax@Ay), mx(Wy_@Wx_ - 2*(Kx@Ky - Ky@Kx)), mx(Wx_@Wy_ + 2*(Kx@Ky - Ky@Kx)))
    w['asym'] = max(w['asym'], mx(Ax - Ax.T)); w['om'] = max(w['om'], mx(om@om + I4), mx(om@Ex + Ex@om), mx(om@T + T@om))
n = len(sols)
ok("2. %d solutions of (C1)-(C4), (C1y)-(C4y) from 240 random starts" % n, n > 50)
ok("2a-b. on all of them: {A_z, W} = 0 (%.0e), W_x^2 = W_y^2 = -1 (%.0e), {W_x, W_y} = 0 (%.0e)" % (w['a'], w['b'], w['b2']),
   max(w['a'], w['b'], w['b2']) < 1e-8)
ok("2c. A_x = -A_z W_x, A_y = -A_z W_y: squares 1 (%.0e), anticommute with A_z (%.0e) and with each other (%.0e), W_y W_x = A_x A_y (%.0e). "
   "They need not be symmetric (largest |A_x - A_x^T| found: %.2f): a solution is IV.4b's up to a change of basis that keeps A_z and A_z T"
   % (w['c1'], w['c2'], w['c3'], w['c4'], w['asym']), max(w['c1'], w['c2'], w['c3'], w['c4']) < 1e-8)
ok("2d. K_x - A_x/2 commutes with A_x, A_y, A_z (%.0e) and equals alpha 1 + beta omega, omega = A_z A_x A_y (remainder %.0e)" % (w['d1'], w['d2']),
   max(w['d1'], w['d2']) < 1e-8)
ok("2e. (C3) gives both: {A_x, T} = {A_y, T} = 0 (%.0e), and beta = 0 on every solution (largest |beta| %.0e): K_x = A_x/2 + alpha 1" % (w['e1'], w['e2']),
   max(w['e1'], w['e2']) < 1e-8)
ok("2f. 2 [K_x, K_y] = A_x A_y = W_y W_x = -W_x W_y (%.0e)" % w['f'], w['f'] < 1e-8)

# ---- 3. the commutant, and the text's aside
Axd = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0.]]); Ayd = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0.]])
def commutant_dim(mats):
    Amat = np.vstack([np.kron(I4, M) - np.kron(M.T, I4) for M in mats]); return 16 - np.linalg.matrix_rank(Amat, tol=1e-10)
omd = Az@Axd@Ayd
ok("3a. the commutant of IV.4b's A_x, A_y, A_z has dimension %d: the multiples of 1 and of omega = A_z A_x A_y, omega^2 = -1, omega antisymmetric"
   % commutant_dim([Axd, Ayd, Az]), commutant_dim([Axd, Ayd, Az]) == 2 and np.allclose(omd@omd, -I4) and np.allclose(omd, -omd.T))
ok("3b. 'A_z A_x A_y = -+ A_z T' (the text before v35) does not hold: A_z T is symmetric with square +1, omega antisymmetric with square -1; "
   "|omega - A_z T| = %.0f, |omega + A_z T| = %.0f; they anticommute (%.0e), on IV.4b's matrices and on every solution (%.0e)"
   % (mx(omd - Ex), mx(omd + Ex), mx(omd@Ex + Ex@omd), w['om']),
   mx(omd - Ex) > 0.5 and mx(omd + Ex) > 0.5 and mx(omd@Ex + Ex@omd) < 1e-12 and w['om'] < 1e-8 and np.allclose(Ex, Ex.T) and np.allclose(Ex@Ex, I4))
bad = []
for be in (0.0, 0.3, -1.1):
    Kx = Axd/2 + 0.2*I4 + be*omd; Wx_ = Kx@Az - Az@Kx
    bad.append((be, mx(Kx@Wx_ - Wx_@Kx - Az), mx(Kx@Ex - Ex@Kx)))
ok("3c. K_x = A_x/2 + 0.2 + beta omega on IV.4b's matrices: (C1), (C2) hold for every beta (%.0e); (C3) fails by %.1f, %.1f at beta = %.1f, %.1f "
   "and holds at beta = 0 (%.0e): (C3) is what removes the omega part" % (max(b[1] for b in bad), bad[1][2], bad[2][2], bad[1][0], bad[2][0], bad[0][2]),
   max(b[1] for b in bad) < 1e-12 and bad[0][2] < 1e-12 and bad[1][2] > 0.1 and bad[2][2] > 0.1)
# the proof's two steps on the document's matrices: T anticommutes with omega; conjugating (1/2)(A_x + T A_x T^-1) = 2 beta omega by T flips the right side only
Ti = np.linalg.inv(T); lhs = (Axd + T@Axd@Ti)/2
ok("3d. the proof's steps: T omega = -omega T (%.0e); T [ (A_x + T A_x T^-1)/2 ] T^-1 is the same matrix (%.0e) while T omega T^-1 = -omega (%.0e): "
   "so 2 beta omega = -2 beta omega, beta = 0" % (mx(T@omd + omd@T), mx(T@lhs@Ti - lhs), mx(T@omd@Ti + omd)),
   mx(T@omd + omd@T) < 1e-12 and mx(T@lhs@Ti - lhs) < 1e-12 and mx(T@omd@Ti + omd) < 1e-12)
print("ALL PASS" if not fails else "%d FAILED" % fails); sys.exit(1 if fails else 0)
