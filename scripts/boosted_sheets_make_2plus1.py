"""Can the 2+1 column be built from POSITIVE 1+1 walkers plus special relativity?  (c = 1, z the walk's clock, loss out.)
1. Boost map (sympy): if Q0(tau, z) solves the 1+1 z-clock lines, then Q(t,x,z) = e^{eta sigma_x/2} Q0(gamma(t - beta x), z) solves the
   2+1 z-clock lines of IV.4c (E) -- a runner along x, K_x = A_x/2 on the counters, the Lorentz change on (t, x).  Exact.
2. The 1+1 piece is positive: Q0 = expected counts of the branching walk (keep; spawn the other sense at rate a). Lattice expectation vs
   a Monte Carlo of that rule.
3. One transverse boost already puts a sign on a reading: e^{eta sigma_x/2} = cosh(eta/2) + sinh(eta/2) sigma_x has a negative entry for
   eta < 0, and the runner's Q_- is negative wherever Q0+ > coth(|eta|/2) Q0-.  The sheet's paths are all positive; the reading is not a count.
4. The 2+1 z-clock walk (IV.4c (E), coins +-1/2) run on transverse data blows up as the spacing shrinks: the transverse term is a turn in the
   clock, with real growth rates +-k (IV.4b "what the turns do").  No positive z-clock process in 2+1, numerically.
5. The same functions solve the t-reading (2+1 Dirac).  The unitary t-clock automaton, started from a boosted sheet or a signed sum of them,
   stays on the analytic solution to O(h).  So: the 2+1 column = boosted positive 1+1 walks + the signs SR puts on the readings."""
import numpy as np, sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
# ---------- 1. boost map, on the full basis of 1+1 z-clock plane waves
t, x, z, eta, a, w = sp.symbols('t x z eta a omega', real=True)
g, b = sp.cosh(eta), sp.tanh(eta); tau = g*(t - b*x)
sx = sp.Matrix([[0, 1], [1, 0]]); E = sp.cosh(eta/2)*sp.eye(2) + sp.sinh(eta/2)*sx
A_, B_, k = sp.symbols('A B k')
sol = sp.solve([sp.I*k*A_ - sp.I*w*A_ - a*B_, sp.I*k*B_ + sp.I*w*B_ - a*A_], [A_, k], dict=True)   # d_z f = -d_tau f + a h ; d_z h = d_tau h + a f
res_ok = True
for s_ in sol:
    fe = s_[A_]*sp.exp(sp.I*(s_[k]*z - w*tau)); he = B_*sp.exp(sp.I*(s_[k]*z - w*tau)); Qe = E*sp.Matrix([fe, he])
    r1 = sp.diff(Qe[0], z) + sp.diff(Qe[0], t) - a*Qe[1] + sp.diff(Qe[1], x)
    r2 = sp.diff(Qe[1], z) - sp.diff(Qe[1], t) - a*Qe[0] - sp.diff(Qe[0], x)
    res_ok &= (sp.simplify(r1.rewrite(sp.exp)) == 0 and sp.simplify(r2.rewrite(sp.exp)) == 0)
ok("1. boost map: Q = e^{eta sigma_x/2} Q0(gamma(t - beta x), z) solves the 2+1 z-clock lines for every 1+1 plane wave (both branches)", res_ok)
# ---------- 2. the positive piece
rng = np.random.default_rng(0); av = 1.0
def branching(h, Z, T=16.0):
    nz = int(round(Z/h)); Lt = int(round(T/h)) + 1; tt = (np.arange(Lt) - Lt//2)*h
    Qp = np.exp(-(tt/0.4)**2); Qm = np.zeros(Lt); Q0 = np.zeros((nz + 1, Lt, 2)); Q0[0, :, 0] = Qp
    for n in range(nz):
        Qp, Qm = np.roll(Qp, 1) + av*h*Qm, np.roll(Qm, -1) + av*h*Qp; Q0[n + 1, :, 0] = Qp; Q0[n + 1, :, 1] = Qm
    return Q0, tt
Q0c, ttc = branching(0.05, 2.0, 12.0)
ok("2a. the 1+1 piece is positive on every z-cut (branching counts)", (Q0c >= 0).all())
R = 20000; Lt = len(ttc); nz = Q0c.shape[0] - 1; starts = rng.choice(Lt, size=R, p=Q0c[0, :, 0]/Q0c[0, :, 0].sum()); cp = np.zeros(Lt); cm = np.zeros(Lt)
for r in range(R):
    ti = np.array([starts[r]]); e = np.array([1])
    for n in range(nz):
        ti = ti + e; spn = rng.random(len(ti)) < av*0.05
        ti = np.concatenate([ti, ti[spn]]); e = np.concatenate([e, -e[spn]])
    np.add.at(cp, ti[e > 0] % Lt, 1); np.add.at(cm, ti[e < 0] % Lt, 1)
cp *= Q0c[0, :, 0].sum()/R; cm *= Q0c[0, :, 0].sum()/R
e1, e2 = np.abs(cp - Q0c[-1, :, 0]).sum()/Q0c[-1, :, 0].sum(), np.abs(cm - Q0c[-1, :, 1]).sum()/Q0c[-1, :, 1].sum()
ok("2b. Monte Carlo of the branching rule vs the lattice expectation at z = 2: rel. L1 %.3f (+), %.3f (-)" % (e1, e2), e1 < 0.05 and e2 < 0.05)
# ---------- helpers: the boosted sheet on the (x, z) grid at time t
def sheet(Q0, tt, etav, tv, xs):
    gv, bv = np.cosh(etav), np.tanh(etav); nz1 = Q0.shape[0]; out = np.zeros((len(xs), nz1, 2))
    for j, xv in enumerate(xs):
        tau_ = gv*(tv - bv*xv)
        fp = np.array([np.interp(tau_, tt, Q0[kk, :, 0], left=0, right=0) for kk in range(nz1)])
        fm = np.array([np.interp(tau_, tt, Q0[kk, :, 1], left=0, right=0) for kk in range(nz1)])
        out[j, :, 0] = np.cosh(etav/2)*fp + np.sinh(etav/2)*fm; out[j, :, 1] = np.sinh(etav/2)*fp + np.cosh(etav/2)*fm
    return out
# ---------- 3. one boost, one sign
h = 0.025; Q0, tt = branching(h, 3.0); xs = (np.arange(int(round(6/h)) + 1) - int(round(3/h)))*h
N = sheet(Q0, tt, -1.0, 0.0, xs); P = sheet(Q0, tt, +1.0, 0.0, xs)
fracN = np.mean(N[..., 1][np.abs(N[..., 1]) > 1e-3] < 0); fracP = np.mean(P[..., 1][np.abs(P[..., 1]) > 1e-3] < 0)
ok("3. runner at eta = -1 reads Q_- < 0 on %.0f%% of the sites where it is nonzero (eta = +1: %.0f%%); every path of the sheet is positive" % (100*fracN, 100*fracP),
   fracN > 0.1 and fracP < 0.01)
# ---------- 4. the z-clock 2+1 walk on transverse data: growth
def zwalk_err(h_):
    Q0_, tt_ = branching(h_, 2.0, 12.0); nz_ = Q0_.shape[0] - 1; xs_ = (np.arange(int(round(4/h_)) + 1) - int(round(2/h_)))*h_
    def step21(Qp, Qm):
        Tp = np.roll(Qp, 1, axis=0) + av*h_*Qm; Tm = np.roll(Qm, -1, axis=0) + av*h_*Qp
        Lp, Lm = np.roll(Tp, 1, axis=1), np.roll(Tm, 1, axis=1); Rp, Rm = np.roll(Tp, -1, axis=1), np.roll(Tm, -1, axis=1)
        return 0.5*Lp + 0.5*Lm + 0.5*Rp - 0.5*Rm, -0.5*Lp + 0.5*Lm + 0.5*Rp + 0.5*Rm
    z0 = nz_//2
    S0 = sheet(Q0_, tt_, 1.0, 0.0, xs_)[:, z0, :] + sheet(Q0_, tt_, -1.0, 0.0, xs_)[:, z0, :]     # data on the cut z = z0 ... as (x, t)? no: sheet gives (x, z); we need (t, x) at fixed z
    return None
# (t, x) data at fixed z: Q_eta(t, x, z0) = E Q0(gamma(t - beta x), z0)
def zwalk_err(h_):
    Q0_, tt_ = branching(h_, 2.0, 12.0); nz_ = Q0_.shape[0] - 1; xs_ = (np.arange(int(round(4/h_)) + 1) - int(round(2/h_)))*h_
    def field(etav, zi):
        gv, bv = np.cosh(etav), np.tanh(etav); op = np.zeros((len(tt_), len(xs_))); om = np.zeros_like(op)
        for j, xv in enumerate(xs_):
            tau_ = gv*(tt_ - bv*xv); fp = np.interp(tau_, tt_, Q0_[zi, :, 0], left=0, right=0); fm = np.interp(tau_, tt_, Q0_[zi, :, 1], left=0, right=0)
            op[:, j] = np.cosh(etav/2)*fp + np.sinh(etav/2)*fm; om[:, j] = np.sinh(etav/2)*fp + np.cosh(etav/2)*fm
        return op, om
    def step21(Qp, Qm):
        Tp = np.roll(Qp, 1, axis=0) + av*h_*Qm; Tm = np.roll(Qm, -1, axis=0) + av*h_*Qp
        Lp, Lm = np.roll(Tp, 1, axis=1), np.roll(Tm, 1, axis=1); Rp, Rm = np.roll(Tp, -1, axis=1), np.roll(Tm, -1, axis=1)
        return 0.5*Lp + 0.5*Lm + 0.5*Rp - 0.5*Rm, -0.5*Lp + 0.5*Lm + 0.5*Rp + 0.5*Rm
    z0 = nz_//2; A1, A2 = field(1.0, z0), field(-1.0, z0); Wp, Wm = A1[0] + A2[0], A1[1] + A2[1]
    for n in range(nz_ - z0): Wp, Wm = step21(Wp, Wm)
    B1, B2 = field(1.0, nz_), field(-1.0, nz_); Sp, Sm = B1[0] + B2[0], B1[1] + B2[1]
    it = slice(int(2/h_), len(tt_) - int(2/h_)); ix = slice(int(0.5/h_), len(xs_) - int(0.5/h_))
    return (np.abs(Wp[it, ix] - Sp[it, ix]).sum() + np.abs(Wm[it, ix] - Sm[it, ix]).sum())/(np.abs(Sp[it, ix]).sum() + np.abs(Sm[it, ix]).sum())
errs = [zwalk_err(h_) for h_ in (0.05, 0.025, 0.0125)]
ok("4. the 2+1 z-clock walk on transverse data: rel. error %.2f -> %.0f -> %.0e as the spacing halves twice: the transverse turn grows as e^{kz}, no z-clock process" % tuple(errs),
   errs[1] > 10*errs[0] and errs[2] > 10*errs[1])
# ---------- 5. the t-clock automaton (2+1 Dirac) on the boosted pieces
def step_t(Q):
    Qn = Q.copy(); Qn[:, :, 0] = np.roll(Q[:, :, 0], 1, axis=1); Qn[:, :, 1] = np.roll(Q[:, :, 1], -1, axis=1)   # Q+ moves +z, Q- moves -z
    u = (Qn[:, :, 0] + Qn[:, :, 1])/np.sqrt(2); v = (Qn[:, :, 0] - Qn[:, :, 1])/np.sqrt(2)                       # the sigma_x streams
    u = np.roll(u, 1, axis=0); v = np.roll(v, -1, axis=0)                                                          # +x and -x
    Qp_ = (u + v)/np.sqrt(2); Qm_ = (u - v)/np.sqrt(2); c, s = np.cos(av*h), np.sin(av*h)
    return np.stack([c*Qp_ + s*Qm_, -s*Qp_ + c*Qm_], axis=2)                                                       # the mass turn T
dT = 0.5; nT = int(round(dT/h)); nz1 = Q0.shape[0]
ix = slice(int(1/h), len(xs) - int(1/h)); iz = slice(int(0.8/h), nz1 - int(0.8/h))
out = []
for name, etas, coef in [("one sheet, eta = -1", [-1.0], [1.0]), ("sum of the +1 and -1 sheets", [1.0, -1.0], [1.0, 1.0]), ("their difference", [1.0, -1.0], [1.0, -1.0])]:
    Qa = sum(c_*sheet(Q0, tt, e_, 0.0, xs) for e_, c_ in zip(etas, coef)); Qb = sum(c_*sheet(Q0, tt, e_, dT, xs) for e_, c_ in zip(etas, coef))
    W = Qa.copy()
    for n in range(nT): W = step_t(W)
    out.append((name, np.abs(W[ix, iz] - Qb[ix, iz]).sum()/np.abs(Qb[ix, iz]).sum(), (W**2).sum()/(Qa**2).sum()))
ok("5. the t-clock automaton (unitary, 2+1 Dirac) from the boosted data, t = 0 -> 0.5, h = %.3f: rel. L1 " % h + ", ".join("%s %.3f" % (n_, e_) for n_, e_, _ in out) + "; length ratio %.6f" % out[1][2],
   all(e_ < 0.05 for _, e_, _ in out) and abs(out[1][2] - 1) < 1e-9)
print("    pieces: positive branching walks in their own frame.  Signs: from e^{eta sigma_x/2} on one runner's reading, and from the coefficients of the sum.")
print("    The signed sum is what the unitary 2+1 step evolves; the z-clock reading of it has no process (4).")
