"""III.15b: the two readings of the walk evolve data off two different cuts, and only one of them is causal.
x-reading (III.14): data on the timelike cut x = 0 (a worldline), kernel supported at x^2 - c^2 t^2 >= 0: spacelike from the source, including t < 0.
t-reading (III.15): data on the spacelike cut t = 0, kernel supported at c^2 t^2 - x^2 >= 0, t > 0: the future cone.
1. the x-clock record at a place x_D integrates to 1 over t (every path crosses once), and a finite fraction of it sits at t < 0: before the emission.
2. the interval sign of every point of the x-reading's support is spacelike (or null); of the t-reading's, timelike (or null) and future.
3. the noise matrices of the two readings (V.3) are the same metric with the clock's sign flipped: diag(+1 on the clock direction's partner...): in the clock z the
   spreading direction is t (III.5: a heat equation in t), in the clock t the spreading directions are space."""
import numpy as np
from scipy.special import i0, i1
from scipy.integrate import quad
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
a, c = 1.0, 1.0
def bulk(t, x):              # III.14 without the delta on the light line
    if x <= c*abs(t): return 0.0
    zeta = a/c*np.sqrt(x**2 - c**2*t**2)
    Pp = np.exp(-a*x/c)*(a/2)*(x + c*t)/np.sqrt(x**2 - c**2*t**2)*i1(zeta)
    Pm = np.exp(-a*x/c)*(a/2)*i0(zeta)
    return Pp + Pm
for xD in (1.0, 3.0):
    tot = quad(bulk, -xD/c, xD/c, args=(xD,), limit=200)[0] + np.exp(-a*xD/c)      # + the delta on the light line
    before = quad(bulk, -xD/c, 0.0, args=(xD,), limit=200)[0]
    ok(f"1. x_D = {xD}: x-clock record over all t = {tot:.6f} (= 1, every path crosses once); fraction registered at t < 0, before the emission: {before:.4f}",
       abs(tot - 1) < 1e-6 and before > 0.05)
# 2. supports vs interval
pts_x = [(t, x) for x in (0.5, 2.0) for t in np.linspace(-x/c, x/c, 7)]
ok("2. every point of the x-reading's support has x^2 - c^2 t^2 >= 0 (spacelike or null), half of them at t < 0; the t-reading's support c t >= |x| is timelike-or-null and t > 0",
   all(x**2 - c**2*t**2 >= -1e-12 for t, x in pts_x) and sum(t < 0 for t, x in pts_x) == 6
   and all(c**2*t**2 - x**2 >= -1e-12 and t > 0 for t in (0.5, 2.0) for x in np.linspace(-c*t, c*t, 7)))
# 3. the noise matrices of V.3: gamma^mu gamma^nu symmetrized, in the clock z and in the clock t (IV.4b real matrices)
import sympy as sp
Az = sp.diag(1, -1, -1, 1); Ax = sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]])
Ay = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]); T = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]])
# clock t (IV.4b (H)): d_t Q = -c(A_x d_x + A_y d_y + A_z d_z)Q + aTQ; the matrices on the non-clock directions anticommute, so the square is diagonal:
sq = lambda Ms: [(M*M)[0, 0] for M in Ms] if all((Ms[i]*Ms[j] + Ms[j]*Ms[i]) == sp.zeros(4) for i in range(len(Ms)) for j in range(i)) else None
ok("3a. clock t: (A_x, A_y, A_z) anticommute and square to (+1, +1, +1): d_t^2 = c^2 (d_x^2 + d_y^2 + d_z^2) - a^2, every non-clock direction spreads with the same sign",
   sq([Ax, Ay, Az]) == [1, 1, 1] and (T*T) == -sp.eye(4))
# clock z (IV.4b (F)): c d_z Q = -A_z d_t Q + c A_xA_z d_x Q + c A_yA_z d_y Q + a A_zT Q: the matrices on (t, x, y) are (A_z, A_xA_z, A_yA_z)
ok("3b. clock z: (A_z, A_xA_z, A_yA_z) anticommute and square to (+1, -1, -1): c^2 d_z^2 = d_t^2 - c^2 d_x^2 - c^2 d_y^2 + a^2: "
   "t spreads with the sign a space direction had in the clock t, and x, y with the opposite sign. The same metric, the clock's sign moved: this is where the Lorentzian sign sits",
   sq([Az, Ax*Az, Ay*Az]) == [1, -1, -1] and ((Az*T)*(Az*T)) == sp.eye(4))
