"""IV.5a (E): the phase rotates each pair of counters.   Recreated 9 Nov 2026 against master v34.

z = u + i v, two real columns of the walk.  The "phase" is not one angle between the two real vectors u and v:
  1. counterexample: z = e^{i theta} (1, i)/sqrt2 has u.v = 0 for every theta; what turns is each pair (u_j, v_j), by the same angle:
     atan2(v_j, u_j) = theta + arg chi_j.
  2. a common phase is the real rotation e^{theta J}, J the quarter-turn between the columns; with z = sqrt(rho) e^{i theta} chi,
     chi = r + i s:  u = sqrt(rho)(r cos theta - s sin theta),  v = sqrt(rho)(r sin theta + s cos theta).
  3. the momentum density in real fields, P_x = hbar sum_j (u_j d_x v_j - v_j d_x u_j) = Re[z^dag (-i hbar d_x) z], site by site;
     = rho d_x S + hbar rho Im(chi^dag d_x chi) for z = sqrt(rho) e^{iS/hbar} chi; integrated it is <p>.  1+1 and a 3+1 column with a
     varying axis.
  4. a plane wave of one sign of energy: v_j(x) = u_j(x - lambda/4) for every component (a quarter wavelength in x; a quarter period
     in t); moving the other way the winding and P_x change sign; each pair's angle turns in t at E/hbar.
Exact arithmetic on finite grids; derivatives are spectral."""
import sys
import numpy as np
fails = 0
def ok(name, cond):
    global fails
    print(("PASS  " if cond else "FAIL  ") + name); fails += (not cond)
mx = lambda M: np.abs(M).max()

# ---- 1
worst = 0.0; ang = 0.0
chi = np.array([1, 1j])/np.sqrt(2)
for th in np.linspace(0, 2*np.pi, 25):
    z = np.exp(1j*th)*chi; u, v = z.real, z.imag
    worst = max(worst, abs(u@v))
    ang = max(ang, mx(np.angle(np.exp(1j*(np.arctan2(v, u) - th - np.angle(chi))))))
ok("1. z = e^{i theta}(1, i)/sqrt2: u.v = 0 at every theta (%.0e) - the two real vectors stay at 90 degrees while the phase turns; "
   "each pair's angle atan2(v_j, u_j) = theta + arg chi_j (%.0e)" % (worst, ang), worst < 1e-15 and ang < 1e-14)

# ---- 2
rng = np.random.default_rng(1); n = 4
r, s = rng.normal(size=n), rng.normal(size=n); nrm = np.sqrt(r@r + s@s); r, s = r/nrm, s/nrm
rho, th = 1.7, 0.83
z = np.sqrt(rho)*np.exp(1j*th)*(r + 1j*s)
Jm = np.block([[np.zeros((n, n)), -np.eye(n)], [np.eye(n), np.zeros((n, n))]])
import scipy.linalg as sl
X0 = np.sqrt(rho)*np.concatenate([r, s]); Xt = sl.expm(th*Jm)@X0
rot = np.block([[np.cos(th)*np.eye(n), -np.sin(th)*np.eye(n)], [np.sin(th)*np.eye(n), np.cos(th)*np.eye(n)]])
ok("2. e^{theta J} = [[cos, -sin],[sin, cos]] in blocks (%.0e); applied to (r, s) it gives u = sqrt(rho)(r cos - s sin), v = sqrt(rho)(r sin + s cos) "
   "= the real and imaginary parts of sqrt(rho) e^{i theta}(r + i s) (%.0e)"
   % (mx(sl.expm(th*Jm) - rot), mx(Xt - np.concatenate([z.real, z.imag]))),
   mx(sl.expm(th*Jm) - rot) < 1e-14 and mx(Xt - np.concatenate([z.real, z.imag])) < 1e-14
   and mx(z.real - np.sqrt(rho)*(r*np.cos(th) - s*np.sin(th))) < 1e-14 and mx(z.imag - np.sqrt(rho)*(r*np.sin(th) + s*np.cos(th))) < 1e-14)

# ---- 3
hb = 1.0; N = 4096; Lx = 240.0; x = np.linspace(-Lx/2, Lx/2, N, endpoint=False); h = x[1] - x[0]; k = 2*np.pi*np.fft.fftfreq(N, h)
d = lambda f: np.fft.ifft(1j*k*np.fft.fft(f, axis=0).T).T if f.ndim == 1 else np.fft.ifft((1j*k)[:, None]*np.fft.fft(f, axis=0), axis=0)
# (i) a 1+1 packet of mean momentum 0.3, two components
env = np.exp(-x**2/(2*6.0**2)); z2 = np.stack([env*np.exp(1j*0.3*x), 0.6*env*np.exp(1j*0.3*x)*np.exp(0.2j)], axis=1); z2 /= np.sqrt((abs(z2)**2).sum()*h)
u, v = z2.real, z2.imag; du, dv = d(u.astype(complex)).real, d(v.astype(complex)).real
P_real = hb*(u*dv - v*du).sum(axis=1); P_cplx = (np.conj(z2)*(-1j*hb*d(z2))).sum(axis=1).real
p_mean = P_real.sum()*h
ok("3a. 1+1 packet: P_x = hbar sum_j (u_j d_x v_j - v_j d_x u_j) equals Re[z^dag(-i hbar d_x)z] site by site (%.0e); integrated: %.5f against "
   "the packet's <p> = %.5f" % (mx(P_real - P_cplx), p_mean, 0.3), mx(P_real - P_cplx) < 1e-12 and abs(p_mean - 0.3) < 1e-9)
# (ii) a column with a varying axis: z = sqrt(rho) e^{iS/hbar} chi(x), chi = (cos(th/2), e^{i ph} sin(th/2)) repeated on four components (two copies)
rho_x = np.exp(-x**2/(2*7.0**2)); S_x = 0.4*x + 0.3*np.sin(0.2*x); th_x = 1.0 + 0.5*np.tanh(x/5); ph_x = 0.7*np.sin(0.15*x)
chi_x = np.stack([np.cos(th_x/2), np.exp(1j*ph_x)*np.sin(th_x/2), 0*x, 0*x], axis=1)
z4 = np.sqrt(rho_x)[:, None]*np.exp(1j*S_x/hb)[:, None]*chi_x
u, v = z4.real, z4.imag; du, dv = d(u.astype(complex)).real, d(v.astype(complex)).real
P_real = hb*(u*dv - v*du).sum(axis=1); P_cplx = (np.conj(z4)*(-1j*hb*d(z4))).sum(axis=1).real
dS = 0.4 + 0.06*np.cos(0.2*x); dph = 0.105*np.cos(0.15*x)
P_fluid = rho_x*dS + hb*rho_x*np.sin(th_x/2)**2*dph                         # Im(chi^dag d_x chi) = sin^2(th/2) d_x phi (IV.8 (A))
ok("3b. a column with a varying axis: the same equality site by site (%.0e), and P_x = rho d_x S + hbar rho Im(chi^dag d_x chi) with "
   "Im(chi^dag d_x chi) = sin^2(theta/2) d_x phi (%.0e): the axis's own twist cannot be dropped (it is %.0f%% of P_x at the centre)"
   % (mx(P_real - P_cplx), mx(P_real - P_fluid), 100*abs(hb*np.sin(th_x[N//2]/2)**2*dph[N//2]/dS[N//2])),
   mx(P_real - P_cplx) < 1e-10 and mx(P_real - P_fluid) < 1e-8)

# ---- 4. plane wave of one sign of energy, 1+1: H = c sz p + m c^2 sx (standard basis)
c, m = 1.0, 0.8; a = m*c*c/hb
res = {}
for kk in (+0.6, -0.6):
    Hm = np.array([[c*hb*kk, m*c*c], [m*c*c, -c*hb*kk]]); wv, vv = np.linalg.eigh(Hm); col = vv[:, 1]; E = wv[1]
    lam = 2*np.pi/abs(kk); xs = np.linspace(0, 2*lam, 400, endpoint=False)
    zf = lambda xx, tt: col[None, :]*np.exp(1j*(kk*xx - E*tt/hb))[:, None]
    z0 = zf(xs, 0.0); shift = zf(xs - np.sign(kk)*lam/4, 0.0)               # displaced a quarter wavelength along the motion
    quarter_x = mx(z0.imag - shift.real)
    quarter_t = mx(z0.imag - zf(xs, +(2*np.pi*hb/E)/4).real)                # v(t) = u(t + quarter period): e^{-iEt}: Im at t = Re at t + T/4
    wind = np.gradient(np.unwrap(np.angle(z0[:, 0])), xs).mean()
    Px = hb*(z0.real*(kk*z0.real) + z0.imag*(kk*z0.imag)).sum(axis=1)       # u d_x v - v d_x u, with d_x v = k u, d_x u = -k v for e^{ikx}
    rate = -np.angle(zf(xs, 1e-3)[:, 0]*np.conj(z0[:, 0])).mean()/1e-3
    res[kk] = (quarter_x, quarter_t, wind, Px.mean(), rate, E)
ok("4. plane wave of one sign of energy: v_j(x) = u_j(x -+ lambda/4) along the motion (%.0e, %.0e) and v_j(t) = u_j(t + quarter period) (%.0e); "
   "the pair's angle winds along x at +k = %+.3f moving right and %+.3f moving left; P_x = %+.3f and %+.3f (it changes sign); each pair turns in t "
   "at %.6f against E/hbar = %.6f"
   % (res[0.6][0], res[-0.6][0], res[0.6][1], res[0.6][2], res[-0.6][2], res[0.6][3], res[-0.6][3], res[0.6][4], res[0.6][5]/hb),
   max(res[0.6][0], res[-0.6][0], res[0.6][1], res[-0.6][1]) < 1e-12 and abs(res[0.6][2] - 0.6) < 1e-6 and abs(res[-0.6][2] + 0.6) < 1e-6
   and abs(res[0.6][3] - 0.6) < 1e-12 and abs(res[-0.6][3] + 0.6) < 1e-12 and abs(res[0.6][4] - res[0.6][5]/hb) < 1e-6)
print("ALL PASS" if not fails else "%d FAILED" % fails); sys.exit(1 if fails else 0)
