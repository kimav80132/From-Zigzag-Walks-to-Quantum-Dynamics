"""III.15a (A): Feynman's checkerboard is the real walk read with t as the clock, with no imaginary unit.
Lattice: column (Q+, Q-), real; per tick dt each sense moves c dt its way, then the coin acts, Q -> (cos(a dt) 1 + sin(a dt) T) Q.
The reading at a site is the sum over every zigzag from (0, 0, +).  Compared with III.15's kernel in the real basis, inside the cone:
    K_-+ = -(a/2c) J0(zeta'),   K_++(bulk) = -(a/2c) (ct + x)/sqrt(c^2 t^2 - x^2) J1(zeta'),   zeta' = (a/c) sqrt(c^2 t^2 - x^2).
(Rebuilt 10 Nov 2026: the earlier version printed its sign check and its reality check without testing them, and tested one component at one N.)"""
import numpy as np
from scipy.special import j0, j1
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
a, c, t = 1.0, 1.0, 4.0
T = np.array([[0., 1.], [-1., 0.]])
def run(N):
    dt = t/N; dx = c*dt; R = np.cos(a*dt)*np.eye(2) + np.sin(a*dt)*T
    X = 2*N + 1; Q = np.zeros((2, X)); Q[0, N] = 1.0
    for _ in range(N):
        Q = R @ np.vstack([np.roll(Q[0], 1), np.roll(Q[1], -1)])          # + moves right, - moves left, then the coin
    n = np.arange(X) - N; xs = n*dx
    occ = (n % 2) == (N % 2)                                               # the sites the walker can be on after N ticks
    return Q, xs, occ, dx
errs = {}
for N in (400, 800, 1600):
    Q, xs, occ, dx = run(N)
    ins = occ & (np.abs(xs) < 0.9*c*t)
    Km_of = lambda sh: -(a/(2*c))*j0((a/c)*np.sqrt(c**2*t**2 - (xs[ins] + sh)**2))
    z = (a/c)*np.sqrt(c**2*t**2 - xs[ins]**2)
    Km = Q[1][ins]/(2*dx); Kp = Q[0][ins]/(2*dx)                          # per unit length: occupied sites are 2 dx apart
    Kp_doc = -(a/(2*c))*(c*t + xs[ins])/np.sqrt(c**2*t**2 - xs[ins]**2)*j1(z)
    top = np.max(np.abs(Km_of(0)))
    errs[N] = (np.max(np.abs(Km - Km_of(0)))/top, np.max(np.abs(Kp - Kp_doc))/np.max(np.abs(Kp_doc)), np.max(np.abs(Km - Km_of(-dx)))/top)
    print(f"      N = {N:5d}: max error / max value   K_-+ {errs[N][0]:.2e}   K_++ {errs[N][1]:.2e}   K_-+ read one site over {errs[N][2]:.2e}")
ok(f"1. the real signed path sum lands on -(a/2c) J0(zeta'): {100*errs[1600][0]:.1f}% of the largest value at N = 1600, t = 4/a, |x| < 0.9 ct, site for site", errs[1600][0] < 0.01)
ok(f"2. the + reading lands on -(a/2c)(ct + x) J1/sqrt(c^2t^2 - x^2) to {errs[1600][1]:.0e}", errs[1600][1] < 1e-4)
ok("3. the 0.7% is first order in the tick (halving the tick halves it) and it is a one-site offset of the - counter, which is written after the shift: "
   f"read one site over, the - reading agrees to {errs[1600][2]:.0e}, second order like the + reading (a quarter per halving)",
   all(0.45 < errs[2*N][0]/errs[N][0] < 0.55 and 0.2 < errs[2*N][1]/errs[N][1] < 0.3 and 0.2 < errs[2*N][2]/errs[N][2] < 0.3 for N in (400, 800)))
Q, xs, occ, dx = run(1600)
ins = occ & (np.abs(xs) < 0.9*c*t); z = (a/c)*np.sqrt(c**2*t**2 - xs[ins]**2)
ok("4. the sign: the path sum has the sign of -(a/2c) J0 at every site (away from J0's zeros); no i was used, the sign is the sign of T",
   np.all(np.sign(Q[1][ins])[np.abs(j0(z)) > 0.05] == np.sign(-j0(z))[np.abs(j0(z)) > 0.05]))
edge = occ & (xs < -0.995*c*t) & (np.abs(xs) < c*t - 2*dx)               # just inside the cone on the - side, zeta' -> 0
ok(f"5. just inside the cone the - reading is -(a/2c): {Q[1][edge][-1]/(2*dx):+.4f} against {-(a/(2*c)):+.4f} (the first flip, + into -, weighs -a)",
   abs(Q[1][edge][-1]/(2*dx) + a/(2*c)) < 0.02)
ok(f"6. the length is kept: sum Q^2 = {np.sum(Q**2):.12f}; the no-flip weight on the light line is cos^N(a dt) = {Q[0, -1]:.6f}",
   abs(np.sum(Q**2) - 1) < 1e-12 and abs(Q[0, -1] - np.cos(a*t/1600)**1600) < 1e-12)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
