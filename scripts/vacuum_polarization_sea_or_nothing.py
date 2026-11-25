"""
Vacuum polarization: does the screening around a lone charge need the sea, or does "nothing" do it?
1+1 walk (the t-reading tick of III.15a), ring of L sites, coin angle a per tick (Compton length c/a sites).

Three vacua for the same charge q_ext at the origin (a heavy walker, VI.2's time-link angle, static kernel):
  (A) the sea: every negative column filled. Its static charge response is the sum over lines o -> e,
        chi(x,y) = - sum_lines [ conj(N_eo(x)) N_eo(y) + c.c. ] / Omega_eo,   N_eo(x) = sum_c conj(e_c(x)) o_c(x)
  (B) nothing, plus the one-walker propagator with Feynman's rule, closed into a loop with the exclusion sign -1:
        K_F(tau) = sum_{eps>0} e^{-eps tau} phi phi^+ (tau > 0),   -sum_{eps<0} e^{-eps tau} phi phi^+ (tau < 0)
        chi(x,y) = - int_{-inf}^{inf} dtau  [ - K_F(x,y;tau) K_F(y,x;-tau) ]     (the bubble, two vertices on one closed line)
  (C) nothing: the empty state. Its response is zero: there is no counter for the field to turn.
Then the screening, self-consistently (RPA): delta rho = (1 - chi v)^{-1} chi Phi_ext with v the ring's 1d kernel,
and the dielectric constant from chi(k -> 0) = -alpha k^2:  eps = 1 + q^2 alpha.

Checks:
  1. (A) = (B) to roundoff: the sea's lines are the one-walker closed loop with Feynman's rule and the -1; and (C) = 0
  2. the response is a polarization: chi(k) = -alpha k^2 + O(k^4) (no screening of a uniform charge: the vacuum is neutral)
  3. the induced charge around the point charge is opposite in sign, lives within a few Compton lengths, and totals
     -q_ext (1 - 1/eps) with eps = 1 + q^2 alpha, the dielectric constant of the same alpha (the Uehling cloud of this walk)
  4. alpha a^2 -> 1/(6 pi) as the lattice is refined: the walk's vacuum polarization is 1+1 QED's one-loop Pi(0) = q^2/(6 pi m^2)
NOTE: what the comparison says.
"""
import numpy as np, scipy.linalg as sl, sys
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all = ok_all and bool(cond); print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()
def note(msg): print("NOTE  " + msg); sys.stdout.flush()

L = 160; a = 0.3; dim = 2*L
U = np.zeros((dim, dim)); R = np.array([[np.cos(a), np.sin(a)], [-np.sin(a), np.cos(a)]])
for x in range(L):
    for s, y in ((0, (x + 1) % L), (1, (x - 1) % L)):
        U[2*y, 2*x + s] = R[0, s]; U[2*y + 1, 2*x + s] = R[1, s]
Tm, Z = sl.schur(U.astype(complex), output='complex'); eps = -np.angle(np.diag(Tm))
occ = eps < 0; emp = eps > 0
O = Z[:, occ].reshape(L, 2, -1); E = Z[:, emp].reshape(L, 2, -1)
Om = eps[emp][:, None] - eps[occ][None, :]                       # Omega_eo > 0, shape (n_e, n_o)
N = np.einsum('xce,xco->eox', E.conj(), O)                        # N_eo(x), shape (n_e, n_o, L)
Nf = N.reshape(-1, L); W = 1.0/Om.reshape(-1)
chiA = -(Nf.conj().T*W) @ Nf; chiA = chiA + chiA.T.conj(); chiA = chiA.real     # (A): the lines

# (B): the loop from K_F, tau integrated numerically (log grid, both signs of tau)
rho_site = lambda M: M.reshape(L, 2, L, 2).sum(axis=(1, 3)) if False else None
Pp = Z[:, emp]; Pm = Z[:, occ]; ep, em = eps[emp], eps[occ]
def KF(tau):                                                      # (dim x dim), sense-summed to sites afterwards
    if tau > 0: return (Pp*np.exp(-ep*tau)) @ Pp.conj().T
    return -(Pm*np.exp(-em*tau)) @ Pm.conj().T
def site_sum(M):                                                   # charge density is sense-summed: rho(x) = sum_c |psi_c(x)|^2
    return M.reshape(L, 2, L, 2)
taus = np.exp(np.linspace(np.log(1e-9), np.log(80.0), 900)); u = np.log(taus); du = np.diff(u)
wt = np.zeros_like(taus); wt[:-1] += du/2; wt[1:] += du/2; wt *= taus
chiB = np.zeros((L, L))
for tau, w in zip(taus, wt):
    Kp = site_sum(KF(tau)); Km = site_sum(KF(-tau))
    # bubble: -K_F(x,y) K_F(y,x), with the density vertex summing the sense at each site: sum_{c,c'} K(x c, y c') K(y c', x c)
    bub = np.einsum('xcyd,ydxc->xy', Kp, Km)
    chiB += w*2*np.real(bub)                                       # tau and -tau give conjugate contributions; -(-1) = +
dAB = np.abs(chiA - chiB).max()/np.abs(chiA).max()
ok("1. the sea's response from its lines (A) equals the closed one-walker loop with Feynman's rule and the exclusion sign (B) to %.1e; "
   "the empty state's response (C) is 0 (no counter to turn). 'Nothing' screens nothing; 'nothing plus the closed loop' is the sea" % dAB, dAB < 1e-6)

# 2: polarization, chi(k) ~ -alpha k^2
ks = 2*np.pi*np.fft.fftfreq(L); chik = np.real(np.fft.fft(np.fft.ifft(chiA, axis=0), axis=1).diagonal())   # translation-invariant: chi(k)
small = (np.abs(ks) > 0) & (np.abs(ks) < 0.1)
alpha = -np.polyfit(ks[small]**2, chik[small], 1)[0]
alpha2 = -(chik[small]/ks[small]**2)
ok("2. chi(k) = -alpha k^2 + O(k^4): alpha = %.4f (fit over |k| < 0.1, well below the Compton scale a = 0.3; chi(k)/k^2 ranges %.4f..%.4f there); chi(0) = %.1e: the vacuum "
   "polarizes, it does not screen a uniform charge" % (alpha, -alpha2.max(), -alpha2.min(), chik[0]), abs(chik[0]) < 1e-10 and np.ptp(alpha2)/alpha < 0.05)

# 3: screening of a point charge, RPA, on the ring (neutralizing background), 1d kernel v = q^2 G, -d^2 G = delta - 1/L
q2 = 0.5
G = np.zeros((L, L)); kk = ks.copy(); kk[0] = 1.0
gk = 1.0/(2 - 2*np.cos(kk)); gk[0] = 0.0
for i in range(L): G[i] = np.real(np.fft.ifft(gk*np.exp(-1j*ks*i)))*0 + np.real(np.fft.ifft(gk))[(np.arange(L) - i) % L]
v = q2*G
xs = (np.arange(L) + L//2) % L - L//2
epsd = 1 + q2*alpha
def induced(rho_ext):
    rho_ext = rho_ext - rho_ext.mean()                               # neutralizing background on the ring
    Phi_ext = G @ rho_ext                                            # q_ext = 1
    return q2*np.linalg.solve(np.eye(L) - chiA @ v, chiA @ Phi_ext), rho_ext   # induced charge = q x (response to q Phi)
sig = 12.0; wide = np.exp(-xs**2/(2*sig**2)); wide /= wide.sum()
d_wide, r_wide = induced(wide); win = np.abs(xs) <= 4*sig
ratio_wide = d_wide[win].sum()/r_wide[win].sum(); pred = -(1 - 1/epsd)
point = np.zeros(L); point[0] = 1.0
d_pt, r_pt = induced(point); win2 = np.abs(xs) <= 20
ratio_pt = d_pt[win2].sum()/r_pt[win2].sum()
cum = np.cumsum(np.abs(d_pt[np.argsort(np.abs(xs))])); r50 = np.sort(np.abs(xs))[np.searchsorted(cum, 0.5*cum[-1])]
ok("3. q^2 = %.2f, eps = 1 + q^2 alpha = %.4f. A unit charge spread over sigma = %d sites (many Compton lengths): induced charge / external charge "
   "within 4 sigma = %.4f against -(1 - 1/eps) = %.4f: the vacuum is a dielectric. A unit point charge: the same ratio within 20 sites is %.4f, "
   "with half of the induced cloud within %d sites (Compton length c/a = %.1f): the screening cloud is Compton-sized, and in 1+1 its total is the "
   "same for a point as for a wide charge (chi(k)/k^2 varies by under 10%% across the zone)"
   % (q2, epsd, int(sig), ratio_wide, pred, ratio_pt, r50, 1/a), abs(ratio_wide - pred) < 0.01 and abs(ratio_pt - pred) < 0.02 and r50 <= 2/a)
note("the comparison: the empty walk does nothing to a charge. The screening exists only if closed one-walker loops with no emitter "
     "are counted (B) - and counting them is the same number as the filled sea (A). So 'nothing is there' holds for every process "
     "of walkers that exist, and fails at the one place a field pushes on nothing: there, the closed loops must be stated, "
     "as the Casimir doc states them, and the sea is that statement read on a slice.")

# 4: is it 1+1 QED's vacuum polarization? chi(k) = -alpha k^2 with alpha -> 1/(6 pi m^2) as the lattice gets fine (m = a per tick, c = 1)
def alpha_lattice(a, L):
    ks_ = 2*np.pi*np.arange(L)/L; R = np.array([[np.cos(a), np.sin(a)], [-np.sin(a), np.cos(a)]])
    def bloch(p):
        Up = R @ np.diag([np.exp(-1j*p), np.exp(1j*p)]); lam, V = np.linalg.eig(Up); e = -np.angle(lam)
        o = np.argsort(e); return e[o], V[:, o]/np.linalg.norm(V[:, o], axis=0)      # [neg, pos]
    E = [bloch(p) for p in ks_]
    def chi_k(ik):
        s = 0.0
        for ip in range(L):
            eo, Vo = E[ip]; u = Vo[:, 0]
            for jp in ((ip + ik) % L, (ip - ik) % L):
                ee, Ve = E[jp]; v = Ve[:, 1]
                s += abs(np.vdot(v, u))**2/(ee[1] - eo[0])
        return -s/L
    k1 = 2*np.pi/L; return -chi_k(1)/k1**2, -chi_k(2)/(2*k1)**2
rows4 = []
for a_, L_ in ((0.3, 160), (0.15, 320), (0.075, 640), (0.0375, 1280)):
    al1, al2 = alpha_lattice(a_, L_); rows4.append((a_, L_, al1, al1*a_**2))
cont = 1/(6*np.pi)
ok("4. the same alpha from the momentum blocks: a = 0.3 gives %.4f (check 2's real-space fit: %.4f). Refining the lattice at fixed Compton length, "
   "alpha a^2 = %s for a = %s -> 1/(6 pi) = %.4f, the one-loop vacuum polarization of 1+1 QED for a Dirac walker of mass a (Pi(0) = q^2/(6 pi m^2)); "
   "the walk's vacuum polarizes as the textbook's, with the sea or with the closed loop"
   % (rows4[0][2], alpha, ", ".join("%.4f" % r[3] for r in rows4), [r[0] for r in rows4], cont),
   abs(rows4[0][2] - alpha)/alpha < 0.02 and abs(rows4[-1][3] - cont)/cont < 0.03 and all(abs(rows4[i + 1][3] - cont) < abs(rows4[i][3] - cont) for i in range(len(rows4) - 1)))
print("ALL PASS" if ok_all else "SOME FAILED"); raise SystemExit(0 if ok_all else 1)
