"""Bondi's radar and the two legs of an event: what is identical, checked.
Radar from the source A = (0,0): a signal leaves at clock time t1, bounces at the event X = (t,x), returns at t2. Then t = (t1+t2)/2, x = c(t2-t1)/2 and
1. s^2 = c^2 t^2 - x^2 = c^2 t1 t2: the interval is the product of the two leg times. Spacelike (s^2 < 0) iff t1 t2 < 0: one leg in the past of the trigger.
2. the support of the walk's own-clock kernel (III.14), x >= c|t|, is exactly the radar bracket t1 <= 0 <= t2; the t-reading's support (III.15), ct >= |x|, is t1, t2 >= 0:
   III.15b's input "a record is data on a spacelike cut" is "both radar legs after the trigger".
3. the Bessel argument of both kernels is zeta = (a/c) s: the flip rate times the interval - the interval counted in flips. Across the cone it turns imaginary.
4. the content of IV.4d: rho^2 - j^2 = (2 Q+ Q-)^2, the product of the forward and backward weights, invariant under the run Q+- -> e^{+-eta/2} Q+-, as t1 t2 is
   under a boost (t1 -> k t1, t2 -> t2/k). The two senses are the two radar legs of the content; the mass term couples to their product (III.10).
5. with a link field the two-leg step stays orthogonal: U[theta]^T U[theta] = 1. The field changes what a loop reads, not the determinant: VI.1 (B) holds with a field."""
import sympy as sp, numpy as np
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
t1, t2, c, a, eta, k = sp.symbols('t1 t2 c a eta k', real=True)
t, x = (t1 + t2)/2, c*(t2 - t1)/2
ok("1. s^2 = c^2 t^2 - x^2 = c^2 t1 t2", sp.simplify(c**2*t**2 - x**2 - c**2*t1*t2) == 0)
rng = np.random.default_rng(0); pts = rng.normal(size=(20000, 2))
tt, xx = pts[:, 0], pts[:, 1]; T1, T2 = tt - xx, tt + xx          # c = 1
ok("2. x >= c|t| (III.14's support) <=> t1 <= 0 <= t2 (one radar leg in the past); ct >= |x| (III.15's support) <=> t1, t2 >= 0 (both after the trigger)",
   np.array_equal(np.abs(tt) <= xx, (T1 <= 0) & (T2 >= 0)) and np.array_equal(tt >= np.abs(xx), (T1 >= 0) & (T2 >= 0)))
s = sp.sqrt(c**2*t1*t2)
ok("3. zeta = (a/c) sqrt(x^2 - c^2 t^2) = (a/c) sqrt(-c^2 t1 t2): flip rate times interval; real outside the cone (t1 t2 < 0), imaginary inside",
   sp.simplify((a/c)*sp.sqrt(x**2 - c**2*t**2) - (a/c)*sp.sqrt(-c**2*t1*t2)) == 0)
Qp, Qm = sp.symbols('Q_p Q_m', real=True)
rho, j = Qp**2 + Qm**2, Qp**2 - Qm**2
ok("4. rho^2 - j^2 = (2 Q+ Q-)^2, and Q+ Q- is invariant under the run Q+ -> e^{eta/2} Q+, Q- -> e^{-eta/2} Q-, as t1 t2 is under t1 -> k t1, t2 -> t2/k",
   sp.simplify(rho**2 - j**2 - (2*Qp*Qm)**2) == 0 and sp.simplify(sp.exp(eta/2)*Qp*sp.exp(-eta/2)*Qm - Qp*Qm) == 0 and sp.simplify(k*t1*t2/k - t1*t2) == 0)
# 5. orthogonality with a link field (lattice of pair_angle_is_a_connection.py)
exec(open('pair_angle_is_a_connection.py').read().split('A = idx(0, +1, 0)')[0].replace('ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)', ''))
ok("5. U[links]^T U[links] = 1 for a random link field: the field keeps the determinant, it changes only what a loop reads", np.allclose(U_th.T @ U_th, np.eye(dim)))
