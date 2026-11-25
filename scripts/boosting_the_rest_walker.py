"""Exploration, 12 Nov 2026 (not cited in the master): the eight numbers and their equation MADE by running (boosting) the walker at rest.

No second list is used anywhere.  IV.4b's real matrices; runs A_k (planes of time and one space direction), rotations R_jk = A_j A_k, w = A_x A_y A_z.
The construction, step by step:
  1. The walker at rest is one list of four counters, and its clock is a turn in a space plane: for every list q there is one space plane, with
     unit normal n(q), such that the mass flip acts on q as a quarter turn of that plane, T q = -(n.R) q; so q(t) = (cos at - sin at n.R) q, and
     the list's flow goes round in that plane at c.  A list is a size, an axis in space and a clock angle about it: q = Q u0, Q = a size and
     the three rotations (four numbers).
  2. A list cannot tell a runner how the walker moves: its flow is always at c along some direction j^, and it is unchanged, up to its size, by
     any run along j^: e^{eta j^.A/2} q = e^{eta/2} q.  Two combinations of the six plane turns do nothing to it.
  3. So keep the run with the walker.  The walker at rest, Q, then turning with its clock, then run to rapidity eta along n (the runner reads
     the walker's clock in his own coordinates, tau = t cosh(eta) - n.x sinh(eta)/c, III.8c (12)):
           M(t, x) = (cosh(eta/2) + sinh(eta/2) n.A) Q (cos a tau - sin a tau R_yz).
     With Q = 1 the product has components on six of the eight matrices: 1, R_yz, A_x, A_y, A_z, w (formulas checked).
  4. Successive runs in different directions bring in the other rotations (IV.4c (G)), and nothing beyond the eight [1, 3 rotations, 3 runs, w].
  5. The equation, by eliminating the run.  From 3:  d_t M = -a cosh(eta) M R_yz,  c d_k M = +a sinh(eta) n_k M R_yz,  and the mirror image
     (every run reversed: the same walker for the runner going the other way) is  M~ = (cosh(eta) - sinh(eta) n.A) M.  Together:
           d_t M + c A_k d_k M = -a M~ R_yz,
     with no eta and no n left in it.  For any product of runs the same, with (cosh eta, sinh eta n) replaced by the walker's velocity read off
     its own density and flow.  Checked also by finite differences on M(t, x) itself.
  6. The state now tells the runner how the walker moves: flow/density = tanh(eta) n, and the four new numbers are tanh(eta/2) (n.A) times the
     four of the walker at rest.
  7. The frame and its mirror image are what a run stretches oppositely, as the forward and the backward counter are in 1+1 (III.15c); the mirror
     image obeys the mirrored line, d_t M~ - c A_k d_k M~ = -a M R_yz.  In 1+1 the same construction, s B(eta) E(theta), gives the four numbers.
  8. Sums: every one of the eight-number states is a sum of two walkers at rest run in opposite directions; the equation of 5 is linear, so it
     holds for every sum.  (A sum of WAYS of one walker, not two walkers.)
  9. No second walker: the two lists are one walker's counters for the two senses in which it can be started (the two columns of the one step
     matrix: III.15's "two counters, four answers"); turning the lists into each other = starting its clock at another angle.
"""
import sys
import numpy as np, scipy.linalg as sl, sympy as sp
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all = ok_all and bool(cond); print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()
rng = np.random.default_rng(3)
Ax = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]], float)
Ay = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]], float)
Az = np.diag([1., -1, -1, 1]); T = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]], float); I4 = np.eye(4); Ti = np.linalg.inv(T)
A3 = [Ax, Ay, Az]; Ryz, Rzx, Rxy = Ay@Az, Az@Ax, Ax@Ay; R3 = [Ryz, Rzx, Rxy]; w = Ax@Ay@Az
names = ['1', 'R_yz', 'R_zx', 'R_xy', 'A_x', 'A_y', 'A_z', 'w']; basis = [I4] + R3 + A3 + [w]; Emat = np.array([m.ravel() for m in basis]).T
coef = lambda X: np.linalg.lstsq(Emat, X.ravel(), rcond=None)[0]
inspan = lambda X: np.abs(Emat@coef(X) - X.ravel()).max() < 1e-11
frame = lambda m: sum(c_*E_ for c_, E_ in zip(m, basis))
tilde = lambda M: T@M@Ti
u0 = np.array([1., 0, 0, 0]); v0 = T@u0
dot3 = lambda nv, mats: sum(x_*m_ for x_, m_ in zip(nv, mats))
run = lambda eta, nv: np.cosh(eta/2)*I4 + np.sinh(eta/2)*dot3(nv, A3)
unit = lambda: (lambda v: v/np.linalg.norm(v))(rng.normal(size=3))
rec1 = lambda q: (q@q, np.array([q@A_@q for A_ in A3]))                        # density and flow of one list
def rec(M):                                                                    # density and flow of the walker with frame M
    a_, b_ = M@u0, M@v0; return a_@a_ + b_@b_, np.array([a_@A_@a_ + b_@A_@b_ for A_ in A3])

# ---- 1. the walker at rest: one list, a clock turning in a space plane
d1 = 0.0; d1b = 0.0; d1c = 0.0
for _ in range(100):
    q = rng.normal(size=4); nq = np.array([-(R_@q)@(T@q) for R_ in R3])/(q@q)
    d1 = max(d1, abs(nq@nq - 1), np.abs(T@q + dot3(nq, R3)@q).max())
    for t_ in (0.3, 1.9):
        qt = sl.expm(0.7*t_*T)@q; d1b = max(d1b, np.abs(qt - (np.cos(0.7*t_)*I4 - np.sin(0.7*t_)*dot3(nq, R3))@q).max())
        rho, j = rec1(qt); d1c = max(d1c, abs(j@nq)/rho, abs(j@j - rho*rho)/rho**2)
Mq = np.array([E_@u0 for E_ in basis[:4]]).T
ok("1. the walker at rest is one list, and its clock is a turn in a space plane: for 100 random lists q there is one unit normal n(q) with "
   "T q = -(n.R) q (%.0e), so q(t) = (cos at - sin at n.R) q (%.0e); the list's flow stays in that plane and has the size c at every moment (%.0e). "
   "Every list is q = Q u0 with Q a size and the three rotations, once (det %.0f): a size, an axis in space and a clock angle about it"
   % (d1, d1b, d1c, np.linalg.det(Mq)), max(d1, d1b, d1c) < 1e-12 and abs(np.linalg.det(Mq)) > 0.5)

# ---- 2. a list cannot tell a runner how the walker moves
d2 = 0.0; ranks = set()
six = A3 + R3
for _ in range(100):
    q = rng.normal(size=4); rho, j = rec1(q); jh = j/np.linalg.norm(j); eta = rng.uniform(-1.5, 1.5)
    d2 = max(d2, np.abs(dot3(jh, A3)@q - q).max(), np.abs(run(eta, jh)@q - np.exp(eta/2)*q).max())
    ranks.add(int(np.linalg.matrix_rank(np.array([G_@q for G_ in six]).T, tol=1e-9)))
ok("2. a list cannot tell a runner how the walker moves: with j^ the direction of its flow, (j^.A) q = q, so a run of ANY rapidity along j^ returns "
   "the same list times e^{eta/2} (%.0e); and the six plane turns move a list in only %s independent ways. A runner cannot read his own speed "
   "relative to the walker off four counters" % (d2, sorted(ranks)), d2 < 1e-12 and ranks == {4})

# ---- 3. one run of the walker at rest: six of the eight numbers
ch, sh, cs, sn, nx, ny, nz = sp.symbols('ch sh cs sn n_x n_y n_z', real=True)
S = lambda X: sp.Matrix(np.round(X).astype(int).tolist())
Msym = (ch*sp.eye(4) + sh*(nx*S(Ax) + ny*S(Ay) + nz*S(Az)))*(cs*sp.eye(4) - sn*S(Ryz))
Esym = sp.Matrix([list(S(E_)) for E_ in basis]).T
csym = ((Esym.T*Esym).inv()*Esym.T*sp.Matrix(list(Msym.applyfunc(sp.expand)))).applyfunc(sp.expand)
want = sp.Matrix([ch*cs, -ch*sn, 0, 0, sh*nx*cs, sh*(ny*cs + nz*sn), sh*(nz*cs - ny*sn), -sh*nx*sn])
resid = (Esym*csym - sp.Matrix(list(Msym.applyfunc(sp.expand)))).applyfunc(sp.expand)
ok("3. the walker at rest (reference, Q = 1), turning with its clock, then run: (ch + sh n.A)(cos - sin R_yz), ch = cosh(eta/2), sh = sinh(eta/2), "
   "cos, sin of a tau. Multiplied out it is  m0 = ch cos, m_yz = -ch sin, m_x = sh n_x cos, m_y = sh (n_y cos + n_z sin), m_z = sh (n_z cos - n_y sin), "
   "m_w = -sh n_x sin, m_zx = m_xy = 0: one run puts the walker on six of the eight - the size, its spin plane, the three planes of time and space, and w",
   (csym - want).applyfunc(sp.expand) == sp.zeros(8, 1) and resid == sp.zeros(16, 1))

# ---- 4. successive runs: the other planes, and nothing beyond the eight
two = run(0.8, [0, 1, 0])@run(0.5, [1, 0, 0]); c2 = coef(two)
samples = []
for _ in range(40):
    M = I4
    for _ in range(3): M = run(rng.uniform(-1, 1), unit())@M
    M = M@(np.cos(0.4)*I4 - np.sin(0.4)*Ryz); samples.append(M)
rank8 = np.linalg.matrix_rank(np.array([coef(M) for M in samples])); allin = all(inspan(M) for M in samples)
ok("4. successive runs: a run along x then a run along y has a part on the plane (x, y), coefficient %.4f = -sinh(0.25) sinh(0.4) (%.0e) - the turn two "
   "runs leave (IV.4c (G)); forty walkers at rest run three times in random directions lie on all eight matrices (rank %d) and on nothing else (%s)"
   % (c2[3], abs(c2[3] + np.sinh(0.25)*np.sinh(0.4)), rank8, allin), abs(c2[3] + np.sinh(0.25)*np.sinh(0.4)) < 1e-12 and rank8 == 8 and allin)

# ---- 5. the equation, by eliminating the run
a_ = 0.9
def walker(Rm, Q, uvel):                                                       # M(t, x) for the run R, the rest frame Q, velocity (u0, u) (c = 1)
    return lambda t_, x_: Rm@Q@(np.cos(a_*(uvel[0]*t_ - uvel[1:]@x_))*I4 - np.sin(a_*(uvel[0]*t_ - uvel[1:]@x_))*Ryz)
d5a = d5b = d5c = d5d = 0.0
for trial in range(40):
    Q = frame(np.concatenate([rng.normal(size=4), np.zeros(4)]))
    if trial < 20:
        eta = rng.uniform(0.1, 1.5); nv = unit(); Rm = run(eta, nv); uvel = np.concatenate([[np.cosh(eta)], np.sinh(eta)*nv])
        d5a = max(d5a, np.abs(tilde(Rm@Q) - (np.cosh(eta)*I4 - np.sinh(eta)*dot3(nv, A3))@Rm@Q).max())
    else:
        Rm = run(rng.uniform(0.1, 1), unit())@run(rng.uniform(0.1, 1), unit())@run(rng.uniform(0.1, 1), unit())
        rho, j = rec(Rm); nrm = np.sqrt(rho*rho - j@j); uvel = np.concatenate([[rho/nrm], j/nrm])       # the velocity, read off density and flow
    d5b = max(d5b, np.abs(tilde(Rm@Q) - (uvel[0]*I4 - dot3(uvel[1:], A3))@Rm@Q).max())
    Mf = walker(Rm, Q, uvel); t0 = rng.normal(); x0 = rng.normal(size=3); h = 1e-5; M0 = Mf(t0, x0)
    dt = (Mf(t0 + h, x0) - Mf(t0 - h, x0))/(2*h); dx = [(Mf(t0, x0 + h*e_) - Mf(t0, x0 - h*e_))/(2*h) for e_ in np.eye(3)]
    d5c = max(d5c, np.abs(dt + a_*uvel[0]*M0@Ryz).max(), max(np.abs(dx[k_] - a_*uvel[1 + k_]*M0@Ryz).max() for k_ in range(3)))
    d5d = max(d5d, np.abs(dt + sum(A_@d_ for A_, d_ in zip(A3, dx)) + a_*tilde(M0)@Ryz).max())
ok("5. THE EQUATION, by eliminating the run. For the walker of 3: d_t M = -a cosh(eta) M R_yz and c d_k M = a sinh(eta) n_k M R_yz (finite differences, "
   "%.0e), and its mirror image - every run reversed, the same walker for the runner going the other way - is M~ = (cosh(eta) - sinh(eta) n.A) M (%.0e). "
   "Put together: d_t M + c A_k d_k M = -a M~ R_yz, with no eta and no n in it (residual by finite differences on M(t, x) itself: %.0e). For walkers run "
   "three times in different directions the same holds, with the velocity read off the walker's own density and flow (%.0e)" % (d5c, d5a, d5d, d5b),
   d5a < 1e-12 and d5b < 1e-11 and d5c < 1e-7 and d5d < 1e-7)

# ---- 6. the state now tells the runner how the walker moves
d6 = 0.0; one_list = []
for _ in range(60):
    eta = rng.uniform(0.1, 1.5); nv = unit(); Q = frame(np.concatenate([rng.normal(size=4), np.zeros(4)])); M = run(eta, nv)@Q
    rho, j = rec(M); d6 = max(d6, np.abs(j/rho - np.tanh(eta)*nv).max())
    cM = coef(M); Mev = frame(np.concatenate([cM[:4], np.zeros(4)])); Mod = frame(np.concatenate([np.zeros(4), cM[4:]]))
    d6 = max(d6, np.abs(Mod - np.tanh(eta/2)*dot3(nv, A3)@Mev).max(), np.abs(Mev - np.cosh(eta/2)*Q).max())
    r1, j1 = rec1(M@u0); one_list.append(np.linalg.norm(j1)/r1)
ok("6. the state tells the runner how the walker moves, 60 random: flow/density of the run walker = tanh(eta) n (%.0e); its four numbers on [1, rotations] "
   "are cosh(eta/2) times those of the walker at rest, and its four on [runs, w] are tanh(eta/2) (n.A) times those: the new four numbers are the record "
   "of the run. The first list of the same walker, alone, has flow/density = %.6f whatever eta is" % (d6, float(np.mean(one_list))),
   d6 < 1e-12 and abs(np.mean(one_list) - 1) < 1e-12 and np.ptp(one_list) < 1e-12)

# ---- 7. the frame and its mirror image: what a run stretches oppositely
d7 = 0.0
for _ in range(30):
    eta = rng.uniform(0.1, 1.5); nv = unit(); L = run(rng.uniform(0.1, 1), unit()); M = frame(rng.normal(size=8))
    d7 = max(d7, np.abs(tilde(run(eta, nv)@M) - run(-eta, nv)@tilde(M)).max(), np.abs(tilde(tilde(M)) - M).max(), np.abs(tilde(M@Ryz) - tilde(M)@Ryz).max(),
             max(np.abs(tilde(A_@M) + A_@tilde(M)).max() for A_ in A3))
A2 = np.diag([1., -1]); T2 = np.array([[0, 1.], [-1, 0]]); I2 = np.eye(2); e_p = np.array([1., 0])
d7b = 0.0; W4 = []
for _ in range(30):
    s_, eta, th = rng.uniform(0.5, 2), rng.uniform(-1, 1), rng.uniform(0, 6)
    W = s_*np.diag([np.exp(eta/2), np.exp(-eta/2)])@(np.cos(th)*I2 + np.sin(th)*T2); W4.append(W.ravel())
    um = s_*np.array([np.exp(eta/2)*np.cos(th), -np.exp(-eta/2)*np.sin(th)]); vm = s_*np.array([-np.exp(eta/2)*np.sin(th), -np.exp(-eta/2)*np.cos(th)])
    d7b = max(d7b, np.abs(W@e_p - um).max(), np.abs(W@(T2@e_p) - vm).max(),
              np.abs(T2@W@np.linalg.inv(T2) - s_*np.diag([np.exp(-eta/2), np.exp(eta/2)])@(np.cos(th)*I2 + np.sin(th)*T2)).max())
ok("7. the frame and its mirror image are what a run stretches oppositely: (run M)~ = (the run the other way) M~, (A_k M)~ = -A_k M~, (M R_yz)~ = M~ R_yz "
   "(%.0e); so the mirrored line of 5 is d_t M~ - c A_k d_k M~ = -a M R_yz: each is fed by the other, one moving with +A, one with -A - the form of the "
   "1+1 lines for the forward and the backward counter. In 1+1 the same construction, s (run)(clock) = s diag(e^{eta/2}, e^{-eta/2})(cos th + sin th T), "
   "applied to the forward counter and to the flipped one, is the moving walker of III.8c (12) with its second column (%.0e), fills all four numbers "
   "(rank %d), and its mirror image is the run reversed" % (d7, d7b, np.linalg.matrix_rank(np.array(W4))), d7 < 1e-12 and d7b < 1e-12 and np.linalg.matrix_rank(np.array(W4)) == 4)

# ---- 8. sums
eta = 0.9; Rp, Rm_ = run(eta, [1, 0, 0]), run(eta, [-1, 0, 0])
cols = [coef(Rp@E_) for E_ in basis[:4]] + [coef(Rm_@E_) for E_ in basis[:4]]; Smat = np.array(cols).T
d8 = 0.0
for _ in range(20):
    m = rng.normal(size=8); sol = np.linalg.solve(Smat, m); Q1 = frame(np.concatenate([sol[:4], np.zeros(4)])); Q2 = frame(np.concatenate([sol[4:], np.zeros(4)]))
    d8 = max(d8, np.abs(Rp@Q1 + Rm_@Q2 - frame(m)).max())
ok("8. sums: every state of eight numbers is the sum of two walkers at rest (each a size, an axis and a clock angle) run with the same rapidity in "
   "opposite directions (det %.2f, 20 random states: %.0e). The equation of 5 is linear, so it holds for every such sum - for every state"
   % (np.linalg.det(Smat), d8), abs(np.linalg.det(Smat)) > 1e-3 and d8 < 1e-12)

# ---- 9. no second walker: the eight numbers are ONE walker's counters for the two ways its clock can be started
def tick(L, th):                                                               # III.15a's tick in 1+1: shift, then the coin
    D = 2*L; Um = np.zeros((D, D)); cc, ss = np.cos(th), np.sin(th)
    for y in range(L):
        yp, ym = (y - 1) % L, (y + 1) % L
        Um[2*y, 2*yp] += cc; Um[2*y, 2*ym + 1] += ss; Um[2*y + 1, 2*yp] += -ss; Um[2*y + 1, 2*ym + 1] += cc
    return Um
Lr = 40; UN = np.linalg.matrix_power(tick(Lr, 0.3), 12); o = Lr//2; Tbig = np.kron(np.eye(Lr), T2)
ef = np.zeros(2*Lr); ef[2*o] = 1.0; eb = np.zeros(2*Lr); eb[2*o + 1] = 1.0     # the walker started in the forward sense / in the backward sense
u11 = UN@ef; v11 = UN@(Tbig@ef)
kern = max(np.abs(v11 + UN@eb).max(), np.abs(u11 - UN[:, 2*o]).max(), np.abs(v11 + UN[:, 2*o + 1]).max())
d9 = 0.0
for th_ in (0.4, 1.3, 2.9):                                                    # turning the two lists into each other = starting the clock at another angle
    start = np.cos(th_)*ef + np.sin(th_)*eb                                    # = e^{-th T} (forward start)
    d9 = max(d9, np.abs((np.cos(th_)*u11 - np.sin(th_)*v11) - UN@start).max(), np.abs((np.sin(th_)*u11 + np.cos(th_)*v11) - UN@(Tbig@start)).max())
d9b = 0.0
for _ in range(20):                                                            # the same in 3+1, mode by mode
    kv = rng.normal(size=3); Um3 = sl.expm(1.1*(-1j*dot3(kv, A3) + 0.9*T)); th_ = rng.uniform(0, 6)
    u3, v3 = Um3@u0, Um3@(T@u0); start = (np.cos(th_)*I4 - np.sin(th_)*T)@u0
    d9b = max(d9b, np.abs((np.cos(th_)*u3 - np.sin(th_)*v3) - Um3@start).max(), np.abs((np.sin(th_)*u3 + np.cos(th_)*v3) - Um3@(T@start)).max())
ok("9. no second walker. In 1+1, 12 ticks on a ring: the first list is the walker started in the forward sense, the second is minus the SAME walker "
   "started in the backward sense - the two columns of the one step matrix, 'two counters, four answers' (III.15) (%.0e). Turning the two lists into "
   "each other by an angle is the same walker with its clock started at that angle, e^{-theta T} (forward start) (%.0e; in 3+1, 20 random modes: %.0e). "
   "The second list is how the walker was started, which the step never reads" % (kern, d9, d9b), kern < 1e-14 and d9 < 1e-13 and d9b < 1e-12)

print("ALL PASS" if ok_all else "SOME FAILED"); raise SystemExit(0 if ok_all else 1)
