"""III.8c: one pool of routes, read two ways.   Recreated 9 Nov 2026 against master v34 (check labels as the text cites them).

 (e1)-(e3)  III.8c (1)-(3), on the lattice tick U_tick = coin . shift (a ring of 24 sites, 9 ticks, a random real source f):
            alpha_b = e_b^T U f = f^T U^T e_b;  W_b = alpha_b^2,  sum_b W_b = |f|^2;
            r_b(t)^T q(t) = alpha_b on every time slice;  sum_b alpha_b r_b(t) = q(t).
 5          (11): at rest u(t) = (cos at, -sin at) = e^{aTt}(1,0); e^{-at} e^{aTt} is the Poisson sum over the number of flips; the angle of
            the summed rotation is a t (the mean number of flips), not (pi/2) a t.
 6, 5d      (12)-(13): the moving clock u = diag(e^{eta/2}, e^{-eta/2})(cos th, -sin th), th = a gamma (t - V x/c^2), solves the lines;
            th is read back from the two counters by undoing the rescaling (5d); S = -hbar th gives p = gamma m V, E = gamma m c^2, mass shell.
 8          (14): rho_u = cosh eta + sinh eta cos 2th, j_u = c (sinh eta + cosh eta cos 2th): one column beats at twice the clock.
 9          (15)-(16): L^2 = -Omega^2; J(k) = G(k)/omega_k is a real operator in x, commutes with the step, J^2 = -1; u + iJu = 2u_(+),
            u - iJu = 2u_(-) for every real column; at a = 0, J is the Hilbert transform (-i sgn k on the + sense, +i sgn k on the -).
 (c), (c')  at rest the twin Tu, the one-sign partner Ju and the retraced walk S u(-t) coincide up to the stated sign, S u(-t) = -Tu.
            Moving: the twin and the partner differ, and the twin is no longer the first column a quarter clock later; but the twin
            and the retraced walk stay the same column, S u(-t) = -(twin), at any speed (v34's text says the three differ: two do).
Units hbar = 1."""
import sys
import numpy as np, sympy as sp, scipy.linalg as sl
fails = 0
def ok(name, cond):
    global fails
    print(("PASS  " if cond else "FAIL  ") + name); fails += (not cond)
mx = lambda M: np.abs(M).max()
Az = np.diag([1., -1.]); T = np.array([[0, 1.], [-1, 0]]); S = np.array([[0, 1.], [1, 0]])

# ---- (e1)-(e3): the lattice
Ls, th, nt = 24, 0.31, 9
ct, st = np.cos(th), np.sin(th)
U1 = np.zeros((2*Ls, 2*Ls))
for y in range(Ls):
    U1[2*y, 2*((y - 1) % Ls)] += ct; U1[2*y, 2*((y + 1) % Ls) + 1] += st
    U1[2*y + 1, 2*((y - 1) % Ls)] += -st; U1[2*y + 1, 2*((y + 1) % Ls) + 1] += ct
rng = np.random.default_rng(2); f = rng.normal(size=2*Ls)
Upow = [np.linalg.matrix_power(U1, n) for n in range(nt + 1)]; U = Upow[nt]
alpha = U@f; alpha_back = np.array([f@(U.T@np.eye(2*Ls)[b]) for b in range(2*Ls)])
ok("(e1) alpha_b = e_b^T U f = f^T U^T e_b for every endpoint b (%.0e); W_b = alpha_b^2 and sum_b W_b = |f|^2 (%.0e); U^T U = 1 (%.0e)"
   % (mx(alpha - alpha_back), abs((alpha**2).sum() - f@f), mx(U.T@U - np.eye(2*Ls))),
   mx(alpha - alpha_back) < 1e-13 and abs((alpha**2).sum()/(f@f) - 1) < 1e-13 and mx(U.T@U - np.eye(2*Ls)) < 1e-13)
w2 = w3 = 0.0
for n in range(nt + 1):
    q = Upow[n]@f; Rn = Upow[nt - n].T                      # column b of Rn is r_b(t_n) = U(t_f, t_n)^T e_b
    w2 = max(w2, mx(Rn.T@q - alpha)); w3 = max(w3, mx(Rn@alpha - q))
ok("(e2) the overlap of the two legs, r_b(t)^T q(t), equals alpha_b on every one of the %d time slices (%.0e)" % (nt + 1, w2), w2 < 1e-13)
ok("(e3) all the returns, weighted by their outward amplitudes, rebuild the outward column: sum_b alpha_b r_b(t) = q(t) on every slice (%.0e)" % w3, w3 < 1e-13)

# ---- 5: the clock at rest
a, t = 0.7, 2.3
u_t = sl.expm(a*T*t)@np.array([1., 0.])
pois = sum(np.exp(-a*t)*(a*t)**n/sp.factorial(n)*np.linalg.matrix_power(T, n) for n in range(60)).astype(float)
ang = np.arctan2(-u_t[1], u_t[0])
ok("5. at rest u(t) = e^{aTt}(1,0) = (cos at, -sin at) (%.0e); e^{-at} e^{aTt} = sum_n e^{-at}(at)^n/n! T^n, the Poisson sum over n flips (%.0e); "
   "the angle of the summed rotation is %.6f = a t = %.6f, the mean number of flips, not (pi/2) a t = %.3f"
   % (mx(u_t - np.array([np.cos(a*t), -np.sin(a*t)])), mx(pois - np.exp(-a*t)*sl.expm(a*T*t)), ang, a*t, np.pi/2*a*t),
   mx(u_t - np.array([np.cos(a*t), -np.sin(a*t)])) < 1e-14 and mx(pois - np.exp(-a*t)*sl.expm(a*T*t)) < 1e-14 and abs(ang - a*t) < 1e-14)

# ---- 6, 5d, 8: the moving clock, symbolically
ts, xs, a_, c_, eta, m_ = sp.symbols('t x a c eta m', positive=True)
V = c_*sp.tanh(eta); gam = sp.cosh(eta); thet = a_*gam*(ts - V*xs/c_**2)
up, um = sp.exp(eta/2)*sp.cos(thet), -sp.exp(-eta/2)*sp.sin(thet)
line_p = sp.diff(up, ts) + c_*sp.diff(up, xs) - a_*um          # d_t u+ = -c d_x u+ + a u-
line_m = sp.diff(um, ts) - c_*sp.diff(um, xs) + a_*up          # d_t u- = +c d_x u- - a u+
ok("6. (12): u = diag(e^{eta/2}, e^{-eta/2})(cos th, -sin th), th = a gamma (t - V x/c^2), tanh eta = V/c, solves both lines of the t-reading exactly",
   sp.simplify(line_p.rewrite(sp.exp)) == 0 and sp.simplify(line_m.rewrite(sp.exp)) == 0)
back = sp.atan2(-sp.exp(eta/2)*um, sp.exp(-eta/2)*up)
num_ok = all(abs(float((sp.atan2(-sp.exp(eta/2)*um, sp.exp(-eta/2)*up) - thet).subs({ts: tv, xs: xv, a_: 0.7, c_: 1.3, eta: ev}))) < 1e-12
             for tv, xv, ev in ((0.3, 0.2, 0.4), (1.1, -0.5, 1.2), (0.05, 0.9, 0.1)))
ok("5d. the angle is read back from the two counters by undoing the rescaling: th = atan2(-e^{eta/2} u_-, e^{-eta/2} u_+) (three points, |th| < pi)", num_ok)
Sx = -thet                                                      # S = -hbar th, hbar = 1; a = m c^2/hbar
p = sp.diff(Sx, xs); E = -sp.diff(Sx, ts)
ok("6'. (13): S = -hbar th: p = d_x S = gamma m V, E = -d_t S = gamma m c^2, and -E^2/c^2 + p^2 = -m^2 c^2 (with a = m c^2/hbar)",
   sp.simplify((p - gam*m_*V).subs(a_, m_*c_**2)) == 0 and sp.simplify((E - gam*m_*c_**2).subs(a_, m_*c_**2)) == 0
   and sp.simplify((-E**2/c_**2 + p**2 + m_**2*c_**2).subs(a_, m_*c_**2).rewrite(sp.exp)) == 0)
rho_u = up**2 + um**2; j_u = c_*(up**2 - um**2)
ok("8. (14): rho_u = cosh eta + sinh eta cos 2th, j_u = c (sinh eta + cosh eta cos 2th): density and current of one column beat at twice the clock; "
   "at rest j_u = c cos 2at",
   sp.simplify((rho_u - sp.cosh(eta) - sp.sinh(eta)*sp.cos(2*thet)).rewrite(sp.exp)) == 0
   and sp.simplify((j_u - c_*(sp.sinh(eta) + sp.cosh(eta)*sp.cos(2*thet))).rewrite(sp.exp)) == 0
   and sp.simplify(j_u.subs(eta, 0) - c_*sp.cos(2*a_*ts)) == 0)

# ---- 9: the one-sign partner
k = sp.symbols('k', real=True); om = sp.sqrt(a_**2 + c_**2*k**2)
G = sp.Matrix([[-sp.I*c_*k, a_], [-a_, sp.I*c_*k]]); Jk = G/om
ok("9a. on a mode: G(k)^2 = -omega_k^2 (L^2 = -Omega^2, Klein-Gordon), J(k) = G(k)/omega_k, J^2 = -1, J commutes with G; at k = 0, J = T",
   sp.simplify(G*G + om**2*sp.eye(2)) == sp.zeros(2, 2) and sp.simplify(Jk*Jk + sp.eye(2)) == sp.zeros(2, 2)
   and sp.simplify(Jk*G - G*Jk) == sp.zeros(2, 2) and Jk.subs(k, 0) == sp.Matrix([[0, 1], [-1, 0]]))
N = 512; Lx = 60.0; xg = np.linspace(-Lx/2, Lx/2, N, endpoint=False); kg = 2*np.pi*np.fft.fftfreq(N, xg[1] - xg[0]); an, cn = 0.7, 1.0
def apply_k(col, fun):
    ck = np.fft.fft(col, axis=0); out = np.empty_like(ck, dtype=complex)
    for i, kk in enumerate(kg): out[i] = fun(kk)@ck[i]
    return np.fft.ifft(out, axis=0)
Gk = lambda kk: np.array([[-1j*cn*kk, an], [-an, 1j*cn*kk]]); wk = lambda kk: np.sqrt(an*an + cn*cn*kk*kk)
Jn = lambda kk: Gk(kk)/wk(kk)
def proj(sign):
    def f(kk):
        Hm = 1j*Gk(kk); ev, Vv = np.linalg.eigh(Hm); v_ = Vv[:, 1] if sign > 0 else Vv[:, 0]; return np.outer(v_, v_.conj())
    return f
u0 = np.stack([np.exp(-(xg - 1.0)**2/3.0)*np.cos(1.3*xg), 0.5*np.exp(-(xg + 2.0)**2/5.0)], axis=1)
Ju = apply_k(u0, Jn)
upl, umi = apply_k(u0, proj(+1)), apply_k(u0, proj(-1))
ok("9b. J is a real operator in x: J u is real for a real column (largest imaginary part %.0e); (16): u + iJu = 2u_(+) (%.0e), u - iJu = 2u_(-) (%.0e), "
   "u_(+-) the positive- and negative-energy parts of u; a real column is half of each (weights %.3f, %.3f)"
   % (mx(Ju.imag), mx(u0 + 1j*Ju - 2*upl), mx(u0 - 1j*Ju - 2*umi), (abs(upl)**2).sum()/(u0**2).sum(), (abs(umi)**2).sum()/(u0**2).sum()),
   mx(Ju.imag) < 1e-12 and mx(u0 + 1j*Ju - 2*upl) < 1e-12 and mx(u0 - 1j*Ju - 2*umi) < 1e-12
   and abs((abs(upl)**2).sum()/(u0**2).sum() - 0.5) < 1e-12)
g = np.exp(-xg**2/4.0)*np.cos(2.0*xg); hil = np.fft.ifft(-1j*np.sign(kg)*np.fft.fft(g)).real      # the Hilbert transform of g
J0 = apply_k(np.stack([g, g], axis=1), lambda kk: np.array([[-1j*np.sign(kk), 0], [0, 1j*np.sign(kk)]]) if kk != 0 else np.zeros((2, 2)))
ok("9c. at a = 0, J(k) = diag(-i sgn k, +i sgn k): on the + sense it is the Hilbert transform, on the - sense minus it (%.0e): "
   "the second copy of a massless column is its Hilbert transform" % max(mx(J0[:, 0].real - hil), mx(J0[:, 1].real + hil)),
   all(abs(complex(Jk.subs({a_: 0, c_: 1.3, k: kv})[i, i]) - sg*(-1j if i == 0 else 1j)) < 1e-14 for kv, sg in ((0.7, 1), (-2.1, -1)) for i in (0, 1))
   and max(mx(J0[:, 0].real - hil), mx(J0[:, 1].real + hil)) < 1e-12)

# ---- (c), (c'): the three columns
tt = sp.symbols('tt', real=True)
ur = sp.Matrix([sp.cos(a_*tt), -sp.sin(a_*tt)]); Tm = sp.Matrix([[0, 1], [-1, 0]]); Sm = sp.Matrix([[0, 1], [1, 0]])
ok("(c) at rest: the twin T u(t), the one-sign partner J u = T u (J = T at k = 0), and the walk retraced with the senses swapped, "
   "S u(-t) = -T u(t): the three coincide up to that sign", sp.simplify(Sm*ur.subs(tt, -tt) + Tm*ur) == sp.zeros(2, 1))
evo = lambda col, t_: apply_k(col, lambda kk: sl.expm(t_*Gk(kk)))
umov = np.stack([np.exp(-xg**2/8.0)*np.cos(0.9*xg), 0*xg], axis=1); t1 = 1.7
u_t1 = evo(umov, t1).real
twin = evo(umov@T.T, t1).real                                    # launched with v(0) = T u(0), then the same lines
partner = apply_k(u_t1, Jn).real                                 # J u(t)
retr = evo(umov, -t1).real@S.T                                   # S u(-t)
quart = evo(umov, t1 + np.pi/(2*an)).real                        # the first column a quarter of the Compton clock later
nrm = np.linalg.norm(u_t1)
d1, d2, d3, d4 = (np.linalg.norm(twin - partner)/nrm, np.linalg.norm(twin + retr)/nrm, np.linalg.norm(partner + retr)/nrm,
                  np.linalg.norm(twin - quart)/nrm)
# the same with a start that is not purely in the forward sense
ugen = np.stack([np.exp(-xg**2/8.0)*np.cos(0.9*xg), 0.6*np.exp(-(xg - 1)**2/6.0)], axis=1)
d5 = np.linalg.norm(evo(ugen@T.T, t1).real + evo(ugen, -t1).real@S.T)/np.linalg.norm(evo(ugen, t1).real)
ok("(c') moving (a packet launched in the forward sense, wavenumbers near 0.9, c k/a = %.1f; relative distances): twin against the one-sign partner "
   "J u: %.2f; twin against the first column a quarter of the Compton clock later: %.2f - different columns. But twin against minus the retraced "
   "walk, S u(-t): %.0e - the SAME column at any speed: S G S = -G, so S u(-t) = U(t) S u(0), and S u(0) = -T u(0) whenever the launch is purely "
   "in the forward sense. (With a start that has a backward part the two differ: %.2f.)" % (0.9*cn/an, d1, d4, d2, d5),
   d1 > 0.05 and d4 > 0.05 and d3 > 0.05 and d2 < 1e-12 and d5 > 0.05)
Gs = sp.Matrix([[-sp.I*c_*k, a_], [-a_, sp.I*c_*k]])
ok("(c'') the reason, as matrices: S G(k) S = -G(k) on every mode with x kept (S sz S = -sz, S T S = -T), and (S + T)(1, 0)^T = 0",
   sp.simplify(Sm*Gs*Sm + Gs) == sp.zeros(2, 2) and (Sm + Tm)*sp.Matrix([1, 0]) == sp.zeros(2, 1))
print("ALL PASS" if not fails else "%d FAILED" % fails); sys.exit(1 if fails else 0)
