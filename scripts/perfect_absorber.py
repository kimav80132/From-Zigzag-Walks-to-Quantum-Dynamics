"""Does a perfect absorber conserve probability?  Dirac 1+1, t-clock lattice, detector at the left edge.
A. 'Perfect' as a boundary condition on the particle alone: left-movers that reach the edge leave; nothing comes back.
   Then the particle's norm DEcreases; what leaves per tick is exactly |psi_-(edge)|^2 = the flux c(|psi_-|^2) into the wall.
   Norm(particle) + sum(absorbed) = 1 at every tick: conservation holds once the absorbed flux is counted.
B. 'Perfect' as a unitary coupling to a detector bit (gauge_detector.py, g = pi/2): the joint norm is 1 at every tick exactly;
   the flipped branch is the record. Taking the flipped branch out (A) is the same bookkeeping.
C. A partial absorber (reflecting fraction) vs perfect: the record = 1 - reflected - in flight in every case.
D. What is NOT conserved: the particle's own norm after a perfect absorption (it drops to the survival fraction). That is
   not a loss of probability; it is the probability moved into the detector's record.
E. Compare with III.1's clock-x walk: a perfect detector at a place x_D is passed exactly once by every path, so the record
   there is int (P+ + P-) dt = 1 with no residue; the t-clock residue (in flight) is the part of the packet that moves away."""
import numpy as np
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
N, T, p_ = 260, 230, 0.1
xs = np.arange(N); th = np.arcsin(p_)
def init():
    g = np.exp(-((xs - 150)/12.0)**2)
    pp = (g*np.exp(-0.6j*xs)).astype(complex); pm = 0.9*g*np.exp(-0.6j*xs + 0.2j)       # heading left
    n = np.sqrt((np.abs(pp)**2 + np.abs(pm)**2).sum()); return pp/n, pm/n
def run(refl=0.0):
    pp, pm = init(); absorbed, flux_hist, escaped = 0.0, [], 0.0
    for _ in range(T):
        out = pm[0]                                      # left-mover on the edge site steps off
        leave = np.sqrt(1 - refl)*out; back = np.sqrt(refl)*out
        absorbed += np.abs(leave)**2; flux_hist.append(np.abs(leave)**2)
        pm[0] = 0.0; escaped += np.abs(pp[-1])**2                 # right-movers leaving through the far (open) edge
        pp, pm = np.roll(pp, 1), np.roll(pm, -1)
        pp[0] = back                                     # reflected part re-enters as a right-mover
        pm[-1] = 0.0
        pp, pm = np.cos(th)*pp + np.sin(th)*pm, -np.sin(th)*pp + np.cos(th)*pm
        if not hasattr(run, 'tot'): run.tot = []
        run.tot.append(absorbed + escaped + (np.abs(pp)**2 + np.abs(pm)**2).sum())
    run.escaped = escaped
    return absorbed, (np.abs(pp)**2 + np.abs(pm)**2).sum(), np.array(flux_hist)
run.tot = []; A_abs, A_left, flux = run(0.0); totA = np.array(run.tot)
A_esc = run.escaped
ok(f"A. perfect absorber: absorbed {A_abs:.5f} + escaped through the far edge {A_esc:.5f} + still in the box {A_left:.5f} = {A_abs + A_esc + A_left:.12f}; max deviation from 1 over all ticks {np.abs(totA - 1).max():.1e}",
   np.abs(totA - 1).max() < 1e-12)
ok(f"   the particle's own norm fell to {A_left:.5f}: not lost, moved into the record (D)", A_left < 0.5)
run.tot = []; B_abs, B_left, _ = run(0.3); totB = np.array(run.tot)
B_esc = run.escaped
ok(f"C. partial absorber (30% reflected): absorbed {B_abs:.5f} + escaped {B_esc:.5f} + in the box {B_left:.5f} = {B_abs + B_esc + B_left:.12f}, conserved at every tick ({np.abs(totB - 1).max():.1e})",
   np.abs(totB - 1).max() < 1e-12)
# B: the unitary detector-bit version with g = pi/2 is the same bookkeeping: flipped branch length per tick = |pm(edge)|^2
pp, pm = init(); Pp = np.stack([pp, 0*pp], 1); Pm = np.stack([pm, 0*pm], 1); rec = 0.0; same = True
for k in range(T):
    a0 = Pm[0, 0].copy(); Pm[0, 1] += a0; Pm[0, 0] = 0.0                 # g = pi/2 turn of the bit at the edge (cos 0, sin 1)
    rec_tick = np.abs(Pm[0, 1])**2; Pm[0, 1] = 0.0; rec += rec_tick     # the flipped branch is removed = recorded
    same &= abs(rec_tick - flux[k]) < 1e-12
    Pp, Pm = np.roll(Pp, 1, 0), np.roll(Pm, -1, 0); Pp[0] = 0; Pm[-1] = 0   # (right-edge outflow dropped here too: same as A)
    Pp, Pm = np.cos(th)*Pp + np.sin(th)*Pm, -np.sin(th)*Pp + np.cos(th)*Pm
ok(f"B. unitary bit at g = pi/2: record per tick identical to the absorbed flux of A ({same}); total record {rec:.5f}", same and abs(rec - A_abs) < 1e-12)
# E: clock-x walk, detector at a place: every path passes once
L, pw = 60, 0.1; Tt = 2*L + 3; o = Tt//2
F = np.zeros((2, Tt)); F[0, o] = 1.0
for _ in range(L):
    F = np.array([(1 - pw)*np.roll(F[0], 1) + pw*F[1], pw*F[0] + (1 - pw)*np.roll(F[1], -1)])
ok(f"E. III.1's walk, detector at the place x = {L}: every path passes once, record = sum_t (P+ + P-) = {F.sum():.12f}, no residue", abs(F.sum() - 1) < 1e-12)
