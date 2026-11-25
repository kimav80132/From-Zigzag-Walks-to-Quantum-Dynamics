"""
Checks for drafts/the_sea_in_coordinates_slowly_2026_11_08.md. No complex continuation anywhere: III.14's real kernel outside
the light cone, its Laplace transform in the non-clock variable t, and the configuration count built from it.

III.14 (the walk in its own clock x, loss out, tau = x/c): for |t| < tau,
    Q_{++}(t;x) = delta(t - tau) + (a/2) sqrt((tau+t)/(tau-t)) I_1(a sqrt(tau^2 - t^2))
    Q_{--}(t;x) = delta(t + tau) + (a/2) sqrt((tau-t)/(tau+t)) I_1(a sqrt(tau^2 - t^2))
    Q_{-+}(t;x) = Q_{+-}(t;x) = (a/2) I_0(a sqrt(tau^2 - t^2))
(the inverse transforms of III.14's cosh(kappa x/c) + i omega sinh(kappa x/c)/kappa and a sinh(kappa x/c)/kappa).
Checks:
  1. the kernel is a probability with the loss put back: e^{-a tau} int dt [Q_{++} + Q_{-+}] = 1
  2. the Laplace transform in t at rate c k reproduces III.14's hat Q at omega = i c k, i.e. e^{tau M_k} with M_k = [[-ck, a],[a, ck]]:
     int dt e^{-ckt} Q_{++} = cosh(omega_k tau) - (ck/omega_k) sinh(omega_k tau), int dt e^{-ckt} Q_{-+} = (a/omega_k) sinh(omega_k tau)
  3. the same-sense count int dt e^{-ckt} [Q_{++} + Q_{--}] = 2 cosh(omega_k tau); at k = 0 it is the plain count 2 cosh(a tau)
  4. the configuration count det(1 + L_k(tau)) = 2 + 2 cosh(omega_k tau); the share of the growing eigen-direction -> 1; -(1/tau) log -> -omega_k
  5. the growing eigenvector of L_k(tau) is, under III.9's D = diag(1, i), the negative-energy eigenvector of h(k)
"""
import sys
import numpy as np
from scipy import integrate, special
def ok(msg, cond): print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()
a = 0.7; c = 1.0
def Qpp_int(t, tau): r = np.sqrt(tau*tau - t*t); return 0.5*a*np.sqrt((tau + t)/(tau - t))*special.i1(a*r)
def Qmm_int(t, tau): r = np.sqrt(tau*tau - t*t); return 0.5*a*np.sqrt((tau - t)/(tau + t))*special.i1(a*r)
def Qmp(t, tau): return 0.5*a*special.i0(a*np.sqrt(tau*tau - t*t))
def lap(f, tau, k, delta_at=None):          # int_{-tau}^{tau} e^{-ckt} f(t) dt  (+ the light-line delta if any)
    val, err = integrate.quad(lambda t: np.exp(-c*k*t)*f(t, tau), -tau, tau, limit=400)
    if delta_at is not None: val += np.exp(-c*k*delta_at)
    return val
tau = 1.3
# 1
tot = lap(Qpp_int, tau, 0.0, delta_at=tau) + lap(Qmp, tau, 0.0)
ok("1. loss put back: e^{-a tau} int dt [Q_++ + Q_-+] = %.10f (a probability: every path of III.1 is somewhere in |t| <= x/c)" % (np.exp(-a*tau)*tot),
   abs(np.exp(-a*tau)*tot - 1) < 1e-9)
# 2, 3
worst = 0.0
for k in (0.0, 0.8, 1.3, -0.8):
    w = np.sqrt(a*a + c*c*k*k)
    Lpp = lap(Qpp_int, tau, k, delta_at=tau); Lmm = lap(Qmm_int, tau, k, delta_at=-tau); Lmp = lap(Qmp, tau, k)
    worst = max(worst, abs(Lpp - (np.cosh(w*tau) - c*k/w*np.sinh(w*tau))), abs(Lmm - (np.cosh(w*tau) + c*k/w*np.sinh(w*tau))),
                abs(Lmp - a/w*np.sinh(w*tau)), abs(Lpp + Lmm - 2*np.cosh(w*tau)))
ok("2/3. Laplace transform in t of III.14's real kernel at rate ck = e^{tau M_k} entry by entry, and the same-sense count = 2cosh(omega_k tau), "
   "to %.1e for k = 0, 0.8, 1.3, -0.8 (k = 0: no weight, the plain count 2cosh(a tau) = %.6f)" % (worst, 2*np.cosh(a*tau)), worst < 1e-8)
# 4, 5
D = np.diag([1.0, 1j]); rows = []
for k in (0.0, 0.8):
    w = np.sqrt(a*a + c*c*k*k)
    for tt in (2.0, 8.0, 20.0):
        Lk = np.array([[np.cosh(w*tt) - c*k/w*np.sinh(w*tt), a/w*np.sinh(w*tt)], [a/w*np.sinh(w*tt), np.cosh(w*tt) + c*k/w*np.sinh(w*tt)]])  # = e^{tt M_k}
        Z = 1 + np.trace(Lk) + 1.0                                     # det(1 + L) = 1 + Tr L + det L, det L = 1
        rows.append((k, tt, Z, np.exp(w*tt)/Z, -np.log(Z)/tt, -w))
    ev, V = np.linalg.eig(Lk); grow = V[:, np.argmax(ev.real)]; hk = np.array([[c*k, 1j*a], [-1j*a, -c*k]])
    m = D @ grow; lam = (np.vdot(m, hk @ m)/np.vdot(m, m)).real; res = abs(hk @ m - lam*m).max()
    ok("5. k = %.1f: the growing eigenvector of L_k maps under D = diag(1, i) (III.9) to an eigenvector of h(k) with eigenvalue %+.4f = -omega_k (residual %.0e)"
       % (k, lam, res), res < 1e-10 and abs(lam + w) < 1e-10)
print("      configuration count det(1 + L_k(tau)) and the share of the growing direction")
for k, tt, Z, share, E0, mw in rows: print("        k = %.1f  tau = %5.1f   Z = %14.4f   share = %.6f   -(1/tau) log Z = %+.6f   (-omega_k = %+.6f)" % (k, tt, Z, share, E0, mw))
ok("4. the configuration count is carried by the growing direction (share -> 1) and -(1/tau) log Z -> -omega_k at tau = 20 for k = 0 and 0.8",
   all(abs(r[4] - r[5]) < 1e-6 and r[3] > 0.9999 for r in rows if r[1] == 20.0))
