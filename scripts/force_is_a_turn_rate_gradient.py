"""VI.2 (D), "How a curvature becomes a push, and the three forces are one shape".   Rebuilt 9 Nov 2026 against master v34.

1. (symbolic) The loop angle of the time-space rectangle.  Link angles (V.1, V.2): a step dx at time t turns the pair by
   theta_x(t) = +(e/hbar) A_x dx, a tick dt at place x by theta_t(x) = -(e/hbar) Phi dt.  The text's loop is
       theta_box = theta_x(t) + theta_t(x+dx) - theta_x(t+dt) - theta_t(x) = (e/hbar)(-d_t A_x - d_x Phi) dt dx = (e/hbar) E_x dt dx.
   Exact for potentials linear in t and x; in general the leading term in dt dx.
2. (lattice, 1+1, the two real legs of V.2, one site per tick, so c = 1 in lattice units)  A walker at rest, wide against its Compton
   length, in Phi = -E x.  Measured separately:
     a. the phase the pair (Q_+, Q'_+) at each site advances per tick, a_eff(x): it is a + (e/hbar) Phi(x), slope -eE;
     b. the force d<p>/dt, from the walker's own lattice momentum: it equals -hbar d a_eff/dx.
   The coin a is in the rate and not in its gradient.
The earlier draft of this script traversed the loop the other way (it got -E), ran the lattice with a packet narrower than the
Compton length, and read the local rate over a fraction of a beat period of an unprojected packet (check 2d shows what that does).  Finite matrices and one lattice run; no fit
other than two straight lines."""
import sys
import numpy as np, sympy as sp
fails = 0
def ok(name, cond):
    global fails
    print(("PASS  " if cond else "FAIL  ") + name); fails += (not cond)

# ---- 1. the time-space plaquette, the text's orientation
t, x, dt, dx, e, hb = sp.symbols('t x dt dx e hbar', real=True)
Phi = sp.Function('Phi')(t, x); A = sp.Function('A')(t, x)
th_x = lambda tt, xx: +(e/hb)*A.subs({t: tt, x: xx}, simultaneous=True)*dx      # a step at time tt
th_t = lambda tt, xx: -(e/hb)*Phi.subs({t: tt, x: xx}, simultaneous=True)*dt    # a tick at place xx
box = th_x(t, x) + th_t(t, x + dx) - th_x(t + dt, x) - th_t(t, x)
lead = sp.diff(box, dt, dx).subs({dt: 0, dx: 0}).doit()                          # the coefficient of dt dx
Ex = -sp.diff(Phi, x) - sp.diff(A, t)
ok("1. theta_box = theta_x(t) + theta_t(x+dx) - theta_x(t+dt) - theta_t(x): its dt dx coefficient is (e/hbar) E_x, "
   "E_x = -d_x Phi - d_t A_x; and no dt-only or dx-only term survives", sp.simplify(lead - (e/hb)*Ex) == 0
   and sp.simplify(box.subs(dx, 0)) == 0 and sp.simplify(box.subs(dt, 0)) == 0)
p0, p1, p2, a0, a1, a2 = sp.symbols('p0 p1 p2 a0 a1 a2', real=True)
lin = {Phi: p0 + p1*t + p2*x, A: a0 + a1*t + a2*x}
box_lin = (+(e/hb)*(a0 + a1*t + a2*x)*dx - (e/hb)*(p0 + p1*t + p2*(x + dx))*dt
           - (e/hb)*(a0 + a1*(t + dt) + a2*x)*dx + (e/hb)*(p0 + p1*t + p2*x)*dt)
ok("1'. for potentials linear in t and x the loop angle is exactly (e/hbar) E_x dt dx", sp.expand(box_lin - (e/hb)*(-p2 - a1)*dt*dx) == 0)

# ---- 2. the lattice: one site per tick; columns (Q+, Q'+, Q-, Q'-)
def run(L, w, a=0.05, eE=2.0e-6, nt=400, navg=126):          # a = coin angle per tick: Compton length 1/a = 20 sites; navg = two beat periods
    xs = np.arange(L) - L/2
    ca, sa = np.cos(a), np.sin(a)
    th = eE*xs                                             # copy turn per tick: theta = -(e/hbar) Phi dt, Phi = -E x  ->  + eE x
    ct, st = np.cos(th), np.sin(th)
    def tick(c):
        qp, qpp, qm, qmp = np.roll(c[0], 1), np.roll(c[1], 1), np.roll(c[2], -1), np.roll(c[3], -1)   # shift: + right, - left
        n0, n2 = ca*qp + sa*qm, -sa*qp + ca*qm             # the coin, cos + sin T, on the unprimed copy
        n1, n3 = ca*qpp + sa*qmp, -sa*qpp + ca*qmp         # and on the primed copy
        return np.array([ct*n0 - st*n1, st*n0 + ct*n1, ct*n2 - st*n3, st*n2 + ct*n3])   # the turn of each pair by theta(x)
    g = np.exp(-xs**2/(2*w*w)); g /= np.sqrt((g*g).sum())
    # the walker at rest: the packet g times the rest column (1, -i)/sqrt2 of III.9's basis, projected onto the positive
    # quasi-energies of the free tick itself (per wavenumber the tick is the 2x2 matrix coin . shift), as VI.2 (G) does.
    # Without the projection the packet carries a little of the other sign of energy (III.8c (E)); it beats at 2a and
    # tilts the local phase by the packet's own slope, the osmotic term of III.11 (see check 2c).
    def project(psi):
        k = 2*np.pi*np.fft.fftfreq(L); pk = np.array([np.fft.fft(psi[0]), np.fft.fft(psi[1])])
        eps = np.arccos(ca*np.cos(k))                      # the tick's quasi-energy, cos eps = cos a cos k
        # eigenvector of R diag(e^{-ik}, e^{ik}) with eigenvalue e^{-i eps}: first row gives (ca e^{-ik} - e^{-i eps}) v0 + sa e^{ik} v1 = 0
        v0 = sa*np.exp(1j*k); v1 = np.exp(-1j*eps) - ca*np.exp(-1j*k)
        nrm = np.sqrt(abs(v0)**2 + abs(v1)**2); v0, v1 = v0/nrm, v1/nrm
        amp = np.conj(v0)*pk[0] + np.conj(v1)*pk[1]
        out = np.array([np.fft.ifft(v0*amp), np.fft.ifft(v1*amp)])
        return out, 1 - (np.abs(out)**2).sum()
    def measure(psi, navg=navg):
        col = np.array([psi[0].real, psi[0].imag, psi[1].real, psi[1].imag]); col /= np.sqrt((col**2).sum())
        def momentum(c):                                   # <p> = sum Im(psi^* d psi), central difference, both senses
            return sum(np.sum(np.conj(c[s] + 1j*c[s+1])*(np.roll(c[s] + 1j*c[s+1], -1) - np.roll(c[s] + 1j*c[s+1], 1))/2).imag for s in (0, 2))
        rates, p = [], [momentum(col)]
        for n in range(nt):
            before = col[0] + 1j*col[1]; col = tick(col); p.append(momentum(col))
            if n < navg: rates.append(-np.angle((col[0] + 1j*col[1])*np.conj(before)))
        a_eff = np.mean(rates, axis=0); m = np.abs(xs) < w
        slope, icpt = np.polyfit(xs[m], a_eff[m], 1)
        return slope, icpt, np.polyfit(np.arange(nt + 1), p, 1)[0], abs((col**2).sum() - 1)
    raw = np.array([g/np.sqrt(2), -1j*g/np.sqrt(2)])
    proj, removed = project(raw)
    return dict(a=a, eE=eE, removed=removed, proj=measure(proj), raw=measure(raw), raw40=measure(raw, 40), proj40=measure(proj, 40))

r1 = run(16384, 1000.0); r0 = run(8192, 500.0)
a, eE = r1['a'], r1['eE']; slope, icpt, force, dn = r1['proj']
ok("2a. the phase the pair advances per tick, site by site, is a + (e/hbar) Phi(x) = a - eE x (packet 50 Compton lengths wide): "
   "at the centre %.6f against a = %.6f; slope %.6e against -eE = %.6e, relative difference %.1e"
   % (icpt, a, slope, -eE, abs(slope/(-eE) - 1)), abs(icpt - a) < 2e-5 and abs(slope/(-eE) - 1) < 3e-4)
ok("2b. the force: d<p>/dt = %.6e against eE = %.6e (relative %.0e); against -d a_eff/dx = %.6e: %.3f%%; total length kept to %.0e"
   % (force, eE, abs(force/eE - 1), -slope, 100*abs(force/(-slope) - 1), dn), abs(force/eE - 1) < 1e-5 and abs(force/(-slope) - 1) < 3e-4 and dn < 1e-12)
d1, d0 = abs(r1['proj'][0]/(-eE) - 1), abs(r0['proj'][0]/(-eE) - 1)
ok("2c. what the difference is: at half the width (25 Compton lengths) it is %.1e, %.2f times larger - it goes as (Compton length/width)^2, "
   "half of it to within %.0f%%; the force itself is eE at both widths (%.0e, %.0e)"
   % (d0, d0/d1, 100*abs(d1/(0.5*(20/1000.0)**2) - 1), abs(r0['proj'][2]/eE - 1), abs(r1['proj'][2]/eE - 1)),
   3.5 < d0/d1 < 4.5 and abs(r0['proj'][2]/eE - 1) < 1e-5)
sr, sr40, sp40 = r0['raw'][0], r0['raw40'][0], r0['proj40'][0]
ok("2d. the beat. The same packet NOT projected onto one sign of energy (weight of the other sign %.1e) beats at 2a, period pi/a = %.1f ticks, "
   "and the beat's phase slopes with x: read over 40 ticks the local slope is %.3e, %.0f%% off -eE; read over two whole beats it is %.3e, "
   "%.2f%% off. Projected, 40 ticks already give %.3e (%.2f%% off). The force is eE in every case (%.0e)"
   % (r0['removed'], np.pi/a, sr40, 100*abs(sr40/(-eE) - 1), sr, 100*abs(sr/(-eE) - 1), sp40, 100*abs(sp40/(-eE) - 1), abs(r0['raw'][2]/eE - 1)),
   abs(sr40/(-eE) - 1) > 0.2 and abs(sr/(-eE) - 1) < 2e-3 and abs(sp40/(-eE) - 1) < 2e-3 and abs(r0['raw'][2]/eE - 1) < 1e-5)
print("NOTE  3. (a statement, not a test) V.5 (C): the detector is pushed by minus the gradient of its own clock shift; the Casimir document: "
      "the walls by -dE_vac/dL; here F = -d(hbar a_eff)/dx. Three gradients of a total turn rate with respect to a position.")
print("ALL PASS" if not fails else "%d FAILED" % fails); sys.exit(1 if fails else 0)
