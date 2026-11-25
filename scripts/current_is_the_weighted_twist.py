"""
Rebuilt 9 Nov 2026 against master v34 (check 2 and check 4 no longer unwrap phases or use np.gradient).
III.11a: the exact relation between the current and the phase slope (Gordon's decomposition in the walk's variables).
Units kept: hbar = 1, c = 1.3, m = 0.8 (so a = m c^2 / hbar = 1.352).
  1. complex form, standard basis of III.9 (i hbar d_t psi = c sz p psi + m c^2 sx psi):
       m j = Re[psi^dag sx p psi] - (hbar/2c) d_t (psi^dag sy psi),     j = c psi^dag sz psi
  2. Madelung form: m j = rho sin(th) cos(phi) d_x Sbar + (hbar/2) rho sin(phi) d_x th - (hbar/2c) d_t(rho sin(th) sin(phi))
  3. one real column, walk's basis (d_t u = -c sz d_x u + a T u):  m j = hbar (u+ d_x u- - u- d_x u+) - (hbar/c) d_t(u+ u-)
  4. slow limit: m c cos(th) -> d_x Sbar on a slow positive-energy packet
  5. plane wave: cos th = pc/E, sin th = mc^2/E, m j = rho p / gamma
"""
import sys
import numpy as np
fails = 0
def ok(msg, cond):
    global fails
    print(("PASS  " if cond else "FAIL  ") + msg); fails += (not cond)
hb, c, m = 1.0, 1.3, 0.8; a = m*c*c/hb
L = 80.0; N = 2048; x = np.linspace(-L/2, L/2, N, endpoint=False); dx = x[1] - x[0]; k = 2*np.pi*np.fft.fftfreq(N, dx)
sz = np.diag([1., -1.]); sx = np.array([[0, 1], [1, 0.]]); sy = np.array([[0, -1j], [1j, 0]]); T = np.array([[0, 1], [-1, 0.]])
Hk = lambda kk: c*hb*kk*sz + m*c*c*sx                               # standard basis
def evolve(f0, t, Hf):
    fk = np.fft.fft(f0, axis=0); out = np.empty_like(fk, dtype=complex)
    for i, kk in enumerate(k):
        w, V = np.linalg.eigh(Hf(kk)); out[i] = V @ (np.exp(-1j*w*t/hb)*(V.conj().T @ fk[i]))
    return np.fft.ifft(out, axis=0)
dxs = lambda f: np.fft.ifft((1j*k)[:, None]*np.fft.fft(f, axis=0), axis=0) if f.ndim == 2 else np.fft.ifft(1j*k*np.fft.fft(f))
bil = lambda A, M, B: np.einsum('xi,ij,xj->x', A.conj(), M, B)

env = np.exp(-x**2/(2*5.0**2))
f0 = np.stack([env*(np.cos(0.7*x) + 0.3j*np.sin(1.1*x)), env*(0.5*np.sin(0.4*x) - 0.2j)], 1)   # mixed signs of energy
t0, h = 1.1, 1e-4
p = evolve(f0, t0, Hk); pp = evolve(f0, t0 + h, Hk); pm = evolve(f0, t0 - h, Hk)
j = c*bil(p, sz, p).real; rho = bil(p, np.eye(2), p).real
conv = bil(p, sx, -1j*hb*dxs(p)).real
Yt = (bil(pp, sy, pp).real - bil(pm, sy, pm).real)/(2*h)
r1 = abs(m*j - (conv - hb/(2*c)*Yt)).max()
ok("1. complex form: m j = Re[psi^dag sx p psi] - (hbar/2c) d_t(psi^dag sy psi) on a mixed-energy packet, max residual %.1e (scale %.2f)" % (r1, abs(m*j).max()), r1 < 1e-6)

# 2. the same identity in III.11's variables.  Every Madelung quantity is taken from the column by its definition, with spectral
#    derivatives and no unwrapping of a phase (the earlier draft unwrapped S+ and S- and used np.gradient, and failed at the points
#    where a component passes near zero and its phase jumps):
#      d_x S+- = hbar Im(psi+-^* d_x psi+-)/rho+-,  rho sin(th) cos(phi) = 2 Re(psi+^* psi-),  rho sin(th) sin(phi) = 2 Im(psi+^* psi-),
#      cos(th) = (rho+ - rho-)/rho,  d_x th = -d_x cos(th)/sin(th),  d_x rho+- = 2 Re(psi+-^* d_x psi+-).
rp = abs(p[:, 0])**2; rm = abs(p[:, 1])**2; rr = rp + rm
dp = dxs(p)
dSp = hb*(np.conj(p[:, 0])*dp[:, 0]).imag/rp; dSm = hb*(np.conj(p[:, 1])*dp[:, 1]).imag/rm; dSb = (dSp + dSm)/2
cross = np.conj(p[:, 0])*p[:, 1]                                    # sqrt(rho+ rho-) e^{i phi}
sth = 2*np.sqrt(rp*rm)/rr; cth = (rp - rm)/rr
drp = 2*(np.conj(p[:, 0])*dp[:, 0]).real; drm = 2*(np.conj(p[:, 1])*dp[:, 1]).real      # d_x rho+-, from the column's own derivative
dcth = ((drp - drm)*rr - (rp - rm)*(drp + drm))/rr**2; dth = -dcth/sth
rhs = 2*cross.real*dSb + hb/2*(2*cross.imag/sth)*dth - hb/(2*c)*Yt   # rho sin(phi) = 2 Im(cross)/sin(th)
mask = (rp > 1e-6*rr.max()) & (rm > 1e-6*rr.max())
r2 = abs((m*j - rhs)[mask]).max()
# and with the angles themselves, on the points where neither component is small (phases by atan2, no unwrapping needed for cos, sin)
th_ = np.arccos(np.clip(cth, -1, 1)); ph_ = np.angle(cross)
rhs_b = rr*np.sin(th_)*np.cos(ph_)*dSb + hb/2*rr*np.sin(ph_)*dth - hb/(2*c)*Yt
r2b = abs((m*j - rhs_b)[mask]).max()
ok("2. Madelung form (1) of III.11a: m j = rho sin(th) cos(phi) d_x Sbar + (hbar/2) rho sin(phi) d_x th - (hbar/2c) d_t(rho sin(th) sin(phi)), "
   "on the same mixed-energy packet, at the %d of %d points where both components exceed 1e-6 of the peak: max residual %.1e, and %.1e with th, phi "
   "as angles (scale %.2f)" % (mask.sum(), N, r2, r2b, abs(m*j).max()), r2 < 1e-6 and r2b < 1e-6 and mask.sum() > N//3)

# 3. one real column in the walk's basis
Gw = lambda kk: -c*1j*kk*sz + a*T
def evolve_real(u0, t):
    fk = np.fft.fft(u0, axis=0); out = np.empty_like(fk)
    for i, kk in enumerate(k):
        w, V = np.linalg.eig(Gw(kk)); out[i] = V @ (np.exp(w*t)*np.linalg.solve(V, fk[i]))
    return np.fft.ifft(out, axis=0).real
u0 = np.stack([env*np.cos(0.6*x), 0.4*env*np.sin(0.9*x)], 1)
u = evolve_real(u0, t0); up = evolve_real(u0, t0 + h); um = evolve_real(u0, t0 - h)
jc = c*(u[:, 0]**2 - u[:, 1]**2)
tw = hb*(u[:, 0]*dxs(u[:, 1]).real - u[:, 1]*dxs(u[:, 0]).real)
dprod = (up[:, 0]*up[:, 1] - um[:, 0]*um[:, 1])/(2*h)
r3 = abs(m*jc - (tw - hb/c*dprod)).max()
ok("3. one real column: m j = hbar (u+ d_x u- - u- d_x u+) - (hbar/c) d_t(u+ u-), max residual %.1e (scale %.2f): the current of one column is the twist between its senses plus the beat" % (r3, abs(m*jc).max()), r3 < 1e-6)

# 4. slow positive-energy packet
f1 = np.stack([env*np.exp(1j*0.1*x), 0*env], 1).astype(complex); fk = np.fft.fft(f1, axis=0)
for i, kk in enumerate(k):
    w, V = np.linalg.eigh(Hk(kk)); fk[i] = np.outer(V[:, 1], V[:, 1].conj()) @ fk[i]
q = np.fft.ifft(fk, axis=0); rp = abs(q[:, 0])**2; rm = abs(q[:, 1])**2
dq = dxs(q); dSb4 = hb*((np.conj(q[:, 0])*dq[:, 0]).imag/rp + (np.conj(q[:, 1])*dq[:, 1]).imag/rm)/2
cth4 = (rp - rm)/(rp + rm)
mask = (rp + rm) > 1e-2*(rp + rm).max()
r4 = abs(m*c*cth4 - dSb4)[mask].max()/abs(dSb4)[mask].max()
ok("4. slow limit (k = 0.1, kc/a = %.2f): max|m c cos th - d_x Sbar| / max|d_x Sbar| = %.3f, of order (v/c)^2 plus the packet's spread: III.11's box" % (0.1*c/a, r4), r4 < 0.1)

# 5. plane wave
kk = 0.9; E = np.sqrt((c*hb*kk)**2 + (m*c*c)**2); w, V = np.linalg.eigh(Hk(kk)); s = V[:, 1]
cth = (abs(s[0])**2 - abs(s[1])**2); sth = 2*abs(s[0]*s[1]); gam = E/(m*c*c)
rho_pw = 1.0; jpw = c*rho_pw*cth
ok("5. plane wave k = %.1f: cos th = %.6f vs pc/E = %.6f, sin th = %.6f vs mc^2/E = %.6f, m j = %.6f vs rho p / gamma = %.6f" % (
   kk, cth, c*hb*kk/E, sth, m*c*c/E, m*jpw, rho_pw*hb*kk/gam), abs(cth - c*hb*kk/E) < 1e-12 and abs(sth - m*c*c/E) < 1e-12 and abs(m*jpw - rho_pw*hb*kk/gam) < 1e-12)

print("ALL PASS" if not fails else "%d FAILED" % fails); sys.exit(1 if fails else 0)
