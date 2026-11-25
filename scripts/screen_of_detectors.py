"""A screen: a detector bit at every site. The packet arrives; each bit is turned by g per tick of presence.
1. Weak coupling: every bit is shifted a little; the shifts form the smooth profile of the time-integrated |psi|^2. Continuous, no click.
2. Perfect absorbers (g = pi/2, branch removed): the joint state becomes a sum of branches 'absorbed at D' with lengths
   |psi(D)|^2-dwell; they add to 1. Which branch is the run: not in the dynamics.
3. The walk's own discreteness: one path of III.1 is absorbed at exactly one place; over many paths the places are distributed as P.
"""
import numpy as np
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
# (Rebuilt 10 Nov 2026: the earlier version printed three lines and tested nothing; in 3 it did not compare the hits with P+ + P-.)
N, T, p_ = 240, 120, 0.1
xs = np.arange(N); th = np.arcsin(p_)
g0 = np.exp(-((xs - 60)/10.0)**2)*np.cos(0.45*xs); g1 = np.exp(-((xs - 60)/10.0)**2)*np.sin(0.45*xs)
psi_p = (g0 + 1j*g1).astype(complex); psi_m = 0.3*(g1 - 1j*g0)
n = np.sqrt((np.abs(psi_p)**2 + np.abs(psi_m)**2).sum()); psi_p /= n; psi_m /= n
screen = slice(130, 200); g = 0.02
dwell = np.zeros(N); absorbed = np.zeros(N); Pp, Pm = psi_p.copy(), psi_m.copy()
for _ in range(T):
    Pp, Pm = np.roll(Pp, 1), np.roll(Pm, -1); Pp, Pm = np.cos(th)*Pp + np.sin(th)*Pm, -np.sin(th)*Pp + np.cos(th)*Pm
    rho = np.abs(Pp)**2 + np.abs(Pm)**2; dwell[screen] += rho[screen]
shift = g*dwell
sc = shift[screen]; rough = np.abs(np.diff(sc, 2)).max()/sc.max()
ok(f"1. weak screen: {int((sc > 1e-6).sum())} of {screen.stop - screen.start} detectors shifted, each by g x its own dwell; largest shift {shift.max():.4f}; "
   f"the profile is smooth (second differences at most {100*rough:.1f}% of the peak); sum of shifts {shift.sum():.4f} = g x {dwell.sum():.3f}",
   (sc > 1e-6).all() and rough < 0.05 and abs(dwell.sum() - 47.7) < 0.05)
Pp, Pm = psi_p.copy(), psi_m.copy()
for _ in range(T):
    Pp, Pm = np.roll(Pp, 1), np.roll(Pm, -1); Pp, Pm = np.cos(th)*Pp + np.sin(th)*Pm, -np.sin(th)*Pp + np.cos(th)*Pm
    taken = np.abs(Pp[screen])**2 + np.abs(Pm[screen])**2; absorbed[screen] += taken; Pp[screen] = 0; Pm[screen] = 0
left = (np.abs(Pp)**2 + np.abs(Pm)**2).sum()
ok(f"2. perfect screen: branches 'absorbed at D' with squared lengths summing to {absorbed.sum():.6f}, plus {left:.6f} that never reached it; total {absorbed.sum() + left:.12f}; "
   f"every branch weight is >= 0 and the dynamics keeps them all (largest {absorbed.max():.4f} at D = {absorbed.argmax()})",
   abs(absorbed.sum() + left - 1) < 1e-12 and absorbed.min() >= 0 and abs(absorbed.sum() - 0.9846) < 1e-4)
rng = np.random.default_rng(0); M = 200_000; L = 70; pw = 0.1
t = np.zeros(M, int); e = np.ones(M, int)
for _ in range(L):
    fl = rng.random(M) < pw; t = np.where(fl, t, t + e); e = np.where(fl, -e, e)
# the exact lattice distribution of the arrival time on the screen x = L: III.2's step, P_e(t; x+1) = (1 - p) P_e(t - e; x) + p P_-e(t; x)
Tt = 2*L + 3; o = Tt//2; F = np.zeros((2, Tt)); F[0, o] = 1.0
for _ in range(L):
    F = np.array([(1 - pw)*np.roll(F[0], 1) + pw*F[1], pw*F[0] + (1 - pw)*np.roll(F[1], -1)])
Pexact = F.sum(0); counts = np.bincount(t + o, minlength=Tt)[:Tt]
big = Pexact*M > 20                                                    # bins with enough expected hits for a chi^2
chi2 = np.sum((counts[big] - Pexact[big]*M)**2/(Pexact[big]*M)); dof = int(big.sum()) - 1
ok(f"3. III.1's walk, a screen at x = {L} with a detector at every t: each of the {M:,} paths hits exactly one detector (hits counted: {counts.sum():,}); "
   f"the hits are distributed as P+ + P-: chi^2 = {chi2:.1f} on {dof} bins; mean {t.mean():+.2f} (exact {np.sum(Pexact*(np.arange(Tt) - o)):+.2f}), spread {t.std():.2f} (exact "
   f"{np.sqrt(np.sum(Pexact*(np.arange(Tt) - o)**2) - np.sum(Pexact*(np.arange(Tt) - o))**2):.2f})",
   counts.sum() == M and abs(Pexact.sum() - 1) < 1e-12 and chi2 < dof + 4*np.sqrt(2*dof))
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
