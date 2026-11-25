"""IV.2a, "a test that was run": III.1 with the radius r as its clock, divided by 4 pi r^2, against the 3+1 source answer
G(s) = m K_1(ms)/(4 pi^2 s), s = sqrt(r^2 - c^2 t^2), m = a/c.  c = 1.  The text says it fails, and how.  Each statement is tested.
(Rebuilt 10 Nov 2026: the earlier version printed a table with no test, and its table had no two points with the same s.)"""
import numpy as np
from scipy.special import i0, i1, k1
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
rng = np.random.default_rng(7)
a, c = 1.0, 1.0; m = a/c
def mc_kernel(R, n=1_000_000, nb=20):
    """III.1 run to the clock value R: exponential waits between flips, sense-averaged start; histogram of t on arrival."""
    t = np.zeros(n); eps = rng.choice([-1.0, 1.0], size=n); x = np.zeros(n); alive = np.ones(n, bool); nflips = np.zeros(n, int)
    while alive.any():
        idx = np.nonzero(alive)[0]
        dx = rng.exponential(c/a, size=idx.size)
        reach = x[idx] + dx >= R
        r_ = idx[reach]; t[r_] += eps[r_]*(R - x[r_])/c; x[r_] = R; alive[r_] = False
        s_ = idx[~reach]; t[s_] += eps[s_]*dx[~reach]/c; x[s_] += dx[~reach]; eps[s_] = -eps[s_]; nflips[s_] += 1
    bins = np.linspace(-R/c, R/c, nb + 1)
    h, _ = np.histogram(t[nflips > 0], bins=bins)            # the no-flip paths sit on the light line: the delta term, kept apart
    return 0.5*(bins[1:] + bins[:-1]), h/n/(bins[1] - bins[0]), np.sqrt(np.maximum(h, 1))/n/(bins[1] - bins[0]), np.mean(nflips == 0)
def exact(t, R):                                             # III.14, bulk, sense-averaged start, with its loss factor
    s = np.sqrt(R*R - c*c*t*t); z = (a/c)*s
    return np.exp(-a*R/c)*(a/2)*(i0(z) + (R/s)*i1(z))
G = lambda s: m*k1(m*s)/(4*np.pi**2*s)
from scipy.integrate import quad
worst = 0.0
print("  r      t       s     MC P/(4 pi r^2)   III.14/(4 pi r^2)   G(s)")
for R in (2.0, 4.0):
    tc, P, dP, p0 = mc_kernel(R); w = tc[1] - tc[0]
    for k in range(1, 19):
        ex = quad(exact, tc[k] - w/2, tc[k] + w/2, args=(R,))[0]/w
        worst = max(worst, abs(P[k] - ex)/dP[k])
        if k in (2, 6, 9):
            s = np.sqrt(R*R - tc[k]**2)
            print(f" {R:.1f}  {tc[k]:+.3f}  {s:.3f}     {P[k]/(4*np.pi*R**2):.5f}          {ex/(4*np.pi*R**2):.5f}          {G(s):.5f}")
    worst0 = abs(p0 - np.exp(-a*R/c))
    ok(f"1. r = {R}: the Monte Carlo of III.1 in the clock r is III.14 (bins: largest deviation {worst:.1f} standard errors; "
       f"no-flip fraction {p0:.4f} against e^(-ar/c) = {np.exp(-a*R/c):.4f})", worst < 4 and worst0 < 2e-3)
# 2. a Green's function is a function of s alone; the radial kernel over 4 pi r^2 is not: same s, two radii
s0 = 1.5; vals = [exact(np.sqrt(R*R - s0*s0), R)/(4*np.pi*R**2) for R in (2.0, 4.0)]
ok(f"2. at the same s = {s0}: P/(4 pi r^2) = {vals[0]:.5f} at r = 2 and {vals[1]:.5f} at r = 4 (ratio {vals[0]/vals[1]:.1f}); G(s) = {G(s0):.5f} at both. "
   "The radial kernel over the Jacobian is not a function of the interval", vals[0]/vals[1] > 3)
# 3. along a fixed r it grows with s while the Green's function falls
ts = np.linspace(0.0, 1.9, 40)[::-1]; ss = np.sqrt(4.0 - ts**2); pr = exact(ts, 2.0)/(4*np.pi*4.0)
ok(f"3. along r = 2: as s goes from {ss[0]:.2f} to {ss[-1]:.2f}, P/(4 pi r^2) rises from {pr[0]:.5f} to {pr[-1]:.5f} while G(s) falls from {G(ss[0]):.5f} to {G(ss[-1]):.5f}",
   np.all(np.diff(pr) > 0) and np.all(np.diff(G(ss)) < 0))
# 4. the loss factor knows the frame: e^{-ar/c} is not a function of s
ok(f"4. the loss factor e^(-ar/c) at the same s = {s0}: {np.exp(-2.0):.4f} at r = 2, {np.exp(-4.0):.4f} at r = 4", abs(np.exp(-2.0) - np.exp(-4.0)) > 0.1)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
