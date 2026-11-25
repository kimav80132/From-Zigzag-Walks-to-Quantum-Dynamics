"""Three steps. 1+1, c = hbar = 1.
1. P IS A COUNT. III.2: P_eps(t;x) dt is the number of crossings of the worldline x = const with sense eps between t and t + dt, per path. Monte Carlo of III.1
   against III.14's formula (bulk and the delta mass), by sense.
2. IN SPACETIME IT IS NOT A COUNT.
   (a) the rule of III.1 is frame-bound: a + segment of path length l has x'-length l/k in the frame of rapidity eta and a - segment k l, so the frame-1 ensemble
       flips at rate a k per unit x' on + segments and a/k on - segments - not III.1 run in frame 2 (rate a on both). Monte Carlo: crossings of x' = const by the
       frame-1 ensemble, histogrammed in t', against III.14's P(t'; x') of the rule run in frame 2: they differ (chi^2 far above the Monte Carlo noise).
   (b) with the loss out, the lines for Q are covariant (half_rapidity_from_the_lines.py) - but Q'+- = e^{-+eta/2} Q+-, and no count per unit of anything
       transforms by a half power of k. P counts but is frame-bound; Q is covariant but does not count.
3. A COUNTING EVERY FRAME AGREES ON, built on the tracks. Boost a solution of the Q-lines numerically (Q' = e^{-eta A/2} Q at the boosted point) and check:
   (a) Q' solves the lines in the primed coordinates (finite differences), the covariance of 2(b) on a real packet;
   (b) the number of + segments crossing a stretch of a - track (a line v = const, u in [u1,u2]) is int 2Q+^2 du, and it is the same number in both frames
       (2Q+^2 -> 2Q+^2/k, du -> k du);
   (c) the total through a spacelike cut, int j^0 dx = 1, is the same in both frames; the proper density 2Q+Q- is the same number at corresponding points."""
import numpy as np
from scipy import special, interpolate
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
rng = np.random.default_rng(11); a = 1.0
# ---------- 1. P counts crossings of x = const ----------
def paths(npaths, X, dx=0.005):
    t = np.zeros(npaths); e = np.ones(npaths); n = int(round(X/dx))
    for _ in range(n):
        t += e*dx
        f = rng.random(npaths) < a*dx; e[f] = -e[f]
    return t, e
X = 2.0; tt, ee = paths(400000, X)
def P_bulk(t, x, eps):            # III.14 bulk, c = 1
    z = a*np.sqrt(np.maximum(x**2 - t**2, 1e-300))
    if eps > 0: return np.exp(-a*x)*0.5*a*(x + t)/np.sqrt(np.maximum(x**2 - t**2, 1e-300))*special.i1(z)
    return np.exp(-a*x)*0.5*a*special.i0(z)
bins = np.linspace(-X, X, 41) + 0.0025; mid = 0.5*(bins[1:] + bins[:-1]); w = bins[1] - bins[0]
hp, _ = np.histogram(tt[(ee > 0) & (tt < X - 1e-9)], bins); hm, _ = np.histogram(tt[ee < 0], bins)
hp = hp/len(tt)/w; hm = hm/len(tt)/w
delta_mass = np.mean((ee > 0) & (np.abs(tt - X) < 1e-9))
from scipy import integrate as _int
binavg = lambda eps: np.array([_int.quad(lambda t_: P_bulk(t_, X, eps), bins[i], bins[i+1], limit=200)[0]/w for i in range(len(bins)-1)])
sel = np.abs(mid) < 0.95*X
Pp_b, Pm_b = binavg(+1), binavg(-1)
err_p = np.max(np.abs(hp[sel] - Pp_b[sel]))/np.max(Pp_b[sel]); err_m = np.max(np.abs(hm[sel] - Pm_b[sel]))/np.max(Pm_b[sel])
ok(f"1. crossings of x = {X} by sense, Monte Carlo vs III.14: max rel. err {err_p:.3f} (+), {err_m:.3f} (-); light-line mass {delta_mass:.4f} vs e^-ax = {np.exp(-a*X):.4f}. P counts crossings of the clock cut",
   err_p < 0.04 and err_m < 0.04 and abs(delta_mass - np.exp(-a*X)) < 0.003)
# ---------- 2a. the same ensemble, counted on a boosted worldline ----------
eta = 0.6; g, b = np.cosh(eta), np.tanh(eta); k = np.exp(eta)
# crossings of the boosted worldline x' = Xp, detected on the fly: x' = g (x - b t)
Xp = 1.5; npaths = 300000; dx_ = 0.005; nst = int(round(6.0/dx_))
t = np.zeros(npaths); e = np.ones(npaths); x = 0.0
xp_prev = g*(x - b*t) - Xp
tpc_list, sense_list = [], []
nflip, lenp = [0, 0], [0.0, 0.0]
for i in range(nst):
    t_new = t + e*dx_; x_new = x + dx_
    xp_new = g*(x_new - b*t_new) - Xp
    c_ = np.sign(xp_new)*np.sign(xp_prev) < 0
    if c_.any():
        frac = xp_prev[c_]/(xp_prev[c_] - xp_new[c_]); tc = t[c_] + frac*(t_new[c_] - t[c_]); xc = x + frac*dx_
        tpc_list.append(g*(tc - b*xc)); sense_list.append(e[c_].copy())
    t, x, xp_prev = t_new, x_new, xp_new
    f = rng.random(npaths) < a*dx_
    nflip[0] += np.sum(f & (e > 0)); nflip[1] += np.sum(f & (e < 0))        # flips made on + and on - segments ...
    lenp[0] += np.sum(e > 0)*dx_/k; lenp[1] += np.sum(e < 0)*dx_*k          # ... and the x'-length travelled on each: dx' = dx/k on +, k dx on -
    e[f] = -e[f]
tpc = np.concatenate(tpc_list); sense = np.concatenate(sense_list)
binsp = np.linspace(-Xp, Xp, 31) + 0.0025; midp = 0.5*(binsp[1:] + binsp[:-1]); wp = binsp[1] - binsp[0]
hpp, _ = np.histogram(tpc[sense > 0], binsp); hpm, _ = np.histogram(tpc[sense < 0], binsp)
hpp = hpp/npaths/wp; hpm = hpm/npaths/wp
selp = np.abs(midp) < 0.85*Xp
Xsave = X; X = Xp; Pp2, Pm2 = np.array([_int.quad(lambda t_: P_bulk(t_, Xp, +1), binsp[i], binsp[i+1], limit=200)[0]/wp for i in range(len(binsp)-1)]), np.array([_int.quad(lambda t_: P_bulk(t_, Xp, -1), binsp[i], binsp[i+1], limit=200)[0]/wp for i in range(len(binsp)-1)]); X = Xsave
dev_p = np.max(np.abs(hpp[selp] - Pp2[selp]))/np.max(Pp2[selp]); dev_m = np.max(np.abs(hpm[selp] - Pm2[selp]))/np.max(Pm2[selp])
ok(f"2a. the frame-1 ensemble counted on the worldline x' = {Xp} of the frame eta = {eta}: max rel. deviation from III.1 run in frame 2 is {dev_p:.2f} (+), {dev_m:.2f} (-): "
   "the crossing count is a different count in a different frame. P is frame-bound", dev_p > 0.15 or dev_m > 0.15)
# the reason, measured: flips per unit x' by sense in the frame-1 ensemble
rp, rm = nflip[0]/lenp[0], nflip[1]/lenp[1]
ok(f"2a'. measured flips per unit x' in the frame-1 ensemble: {rp:.3f} on + segments (a k = {a*k:.3f}), {rm:.3f} on - segments (a/k = {a/k:.3f}); "
   f"III.1 run in frame 2 has {a:.3f} on both", abs(rp/(a*k) - 1) < 0.01 and abs(rm/(a/k) - 1) < 0.01)
# ---------- 3. the counting built on the tracks ----------
L_, N_ = 160.0, 4096; xg = (np.arange(N_) - N_/2)*(L_/N_); kk = 2*np.pi*np.fft.fftfreq(N_, d=L_/N_)
A2 = np.diag([1., -1.]); T2 = np.array([[0, 1], [-1, 0.]])
Q0 = np.stack([np.exp(-(xg - 3)**2/(2*4.0**2)), 0.3*np.exp(-(xg + 2)**2/(2*5.0**2))], 1); Q0 /= np.sqrt(np.sum(Q0**2)*(L_/N_))
Qh0 = np.fft.fft(Q0, axis=0)
def Q_at(t):                                   # G(k) = -i k A + a T, G^2 = -(k^2 + a^2): exp(G t) = cos(w t) + G sin(w t)/w
    w_ = np.sqrt(kk**2 + a**2)
    G = -1j*kk[:, None, None]*A2[None] + a*T2[None]
    E_ = np.cos(w_*t)[:, None, None]*np.eye(2)[None] + G*(np.sin(w_*t)/w_)[:, None, None]
    out = np.einsum('nij,nj->ni', E_, Qh0)
    return np.real(np.fft.ifft(out, axis=0))
ts = np.linspace(-24, 24, 481); Qgrid = np.stack([Q_at(t_) for t_ in ts])          # (nt, nx, 2)
interp = [interpolate.RegularGridInterpolator((ts, xg), Qgrid[:, :, c], method='cubic') for c in (0, 1)]
Lam_half = np.diag([np.exp(-eta/2), np.exp(eta/2)])                               # Q' = e^{-eta A/2} Q
def Q_frame2(tp, xp_):                                                            # Q'(t',x'): unprimed point t = g(t' + b x'), x = g(x' + b t')
    t = g*(tp + b*xp_); x = g*(xp_ + b*tp); pts = np.stack([t, x], -1)
    Q = np.stack([interp[0](pts), interp[1](pts)], -1); return Q @ Lam_half.T
# 3a. Q' solves the lines in the primed frame: d_t' Q' + A d_x' Q' - a T Q' = 0 by central differences
h = 0.02; tp0, xp0 = np.meshgrid(np.linspace(-2, 2, 9), np.linspace(-6, 8, 15), indexing='ij')
dtQ = (Q_frame2(tp0 + h, xp0) - Q_frame2(tp0 - h, xp0))/(2*h); dxQ = (Q_frame2(tp0, xp0 + h) - Q_frame2(tp0, xp0 - h))/(2*h)
res = dtQ + dxQ @ A2.T - a*(Q_frame2(tp0, xp0) @ T2.T)
ok(f"3a. Q' = e^(-eta A/2) Q(boosted point) solves the lines in the primed coordinates: max residual {np.max(np.abs(res)):.1e} against max |d_t Q'| {np.max(np.abs(dtQ)):.2f}",
   np.max(np.abs(res)) < 2e-3*np.max(np.abs(dtQ)))
# 3b. + crossings of a stretch of a - track: v = v0, u in [u1, u2]; density per unit u is 2Q+^2
v0, u1, u2 = 1.0, -9.0, 1.0; us = np.linspace(u1, u2, 801)
def Qpt_frame1(u, v): t = (u + v)/2; x = (v - u)/2; pts = np.stack([t, x], -1); return np.stack([interp[0](pts), interp[1](pts)], -1)
flux1 = np.trapezoid(2*Qpt_frame1(us, v0*np.ones_like(us))[:, 0]**2, us)
# the same stretch in frame 2: u' = k u, v' = v/k
usp = k*us; vp0 = v0/k
def Qpt_frame2(up, vp): tp = (up + vp)/2; xp_ = (vp - up)/2; return Q_frame2(tp, xp_)
flux2 = np.trapezoid(2*Qpt_frame2(usp, vp0*np.ones_like(usp))[:, 0]**2, usp)
ok(f"3b. + crossings of the - track v = {v0}, u in [{u1},{u2}]: frame 1 {flux1:.6f}, frame 2 {flux2:.6f} (2Q+^2 -> 2Q+^2/k, du -> k du): the same number",
   abs(flux1/flux2 - 1) < 2e-3)
# 3c. total through a spacelike cut, and the proper density at corresponding points
tot1 = np.sum(Qgrid[240]**2)*(L_/N_)                                               # t = 0
xps = np.linspace(-24, 24, 4801); Q2 = Q_frame2(np.zeros_like(xps), xps); tot2 = np.trapezoid(np.sum(Q2**2, 1), xps)
pt = np.array([0.7, 2.0]); Q1pt = np.stack([interp[0](pt[None]), interp[1](pt[None])], -1)[0]
tp_, xp_ = g*(pt[0] - b*pt[1]), g*(pt[1] - b*pt[0]); Q2pt = Q_frame2(np.array([tp_]), np.array([xp_]))[0]
ok(f"3c. int j^0 dx = {tot1:.6f} (frame 1), int j'^0 dx' = {tot2:.6f} (frame 2); proper density 2Q+Q- at one event: {2*Q1pt[0]*Q1pt[1]:.6f} and {2*Q2pt[0]*Q2pt[1]:.6f}",
   abs(tot1 - 1) < 1e-6 and abs(tot2 - 1) < 2e-3 and abs(2*Q1pt[0]*Q1pt[1] - 2*Q2pt[0]*Q2pt[1]) < 1e-9)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
