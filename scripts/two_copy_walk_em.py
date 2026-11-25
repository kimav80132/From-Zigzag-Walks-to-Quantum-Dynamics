"""The charged walk, built step by step (1+1, clock x).  Two copies of the time-flip walk: counters P+^1, P-^1 (copy 1) and
P+^2, P-^2 (copy 2).  Column (P+^1, P+^2, P-^1, P-^2).  One clock step dx:
  advance : sense eps moves t by eps dx/c (carried by the counter)
  flip    : with probability a dx/c the sense reverses (copy kept)
  loss    : the complement
  turn    : the two copies turn into each other by the angle  theta_eps = -(e/hbar) [ eps Phi dx/c - A_x dx ]
            (the phase a step of the path picks up: -e Phi dt + e A_x dx with dt = eps dx/c, divided by hbar)
One step -> Taylor in dx -> the lines.  Then: loss out, read in t, pair the copies as one complex column -> Dirac with the field.
Checks: (1) the step's exact lines; (2) copy-length conserved by the turn; (3) the t-reading equals Dirac with Phi, A_x;
        (4) gauge change: Phi -> Phi - d_t chi, A_x -> A_x + d_x chi is the same walk with the copies turned by e chi/hbar;
        (5) a lattice run of the two-copy walk against the exact Dirac evolution with a potential step."""
import sympy as sp, numpy as np
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
t, x, dx = sp.symbols('t x dx', real=True); c, a, e, hb = sp.symbols('c a e hbar', positive=True)
Phi, Ax, chi = [sp.Function(n, real=True)(t, x) for n in ('Phi', 'Ax', 'chi')]
P = {(s, k): sp.Function(f'P{"p" if s > 0 else "m"}{k}')(t, x) for s in (1, -1) for k in (1, 2)}
R = lambda th: sp.Matrix([[sp.cos(th), -sp.sin(th)], [sp.sin(th), sp.cos(th)]])     # turn of (copy1, copy2)
def step(s):
    th = -(e/hb)*(s*Phi*dx/c - Ax*dx)
    stay = sp.Matrix([P[(s, 1)].subs(t, t - s*dx/c), P[(s, 2)].subs(t, t - s*dx/c)])      # advance, no flip
    flip = sp.Matrix([P[(-s, 1)], P[(-s, 2)]])                                              # flip, no advance
    return R(th)*((1 - a*dx/c)*stay + a*dx/c*flip)
first = lambda ex: sp.expand(sp.series(ex, dx, 0, 2).removeO())
lines = {}
for s in (1, -1):
    new = step(s)
    for k in (1, 2):
        lines[(s, k)] = sp.simplify((first(new[k - 1]) - P[(s, k)])*c/dx)                 # c d_x P
expect = lambda s, k: (-s*P[(s, k)].diff(t) - a*P[(s, k)] + a*P[(-s, k)]
                       + (-1)**(k + 1)*(e/hb)*(s*Phi - c*Ax)*P[(s, 3 - k)])               # turn by theta = -(e/hbar)(eps Phi - c Ax) dx/c: copy1 gains -theta copy2
ok("1. one step -> lines:  c d_x P_eps^1 = -eps d_t P^1 - a P^1 + a P_-eps^1 + (e/hbar)(eps Phi - c A_x) P_eps^2,  and copy 2 with - sign",
   all(sp.simplify(lines[(s, k)] - expect(s, k)) == 0 for s in (1, -1) for k in (1, 2)))
ok("2. the turn keeps (P^1)^2 + (P^2)^2 of each sense: R^T R = 1", sp.simplify(R(sp.Symbol('th')).T*R(sp.Symbol('th')) - sp.eye(2)) == sp.zeros(2))
# loss out, read in t
Q = {key: sp.Function('Q' + str(key[0]) + str(key[1]))(t, x) for key in P}
sub = {P[key]: sp.exp(-a*x/c)*Q[key] for key in P}
tread = {}
for s in (1, -1):
    for k in (1, 2):
        ln = sp.expand(sp.simplify((c*P[(s, k)].diff(x) - expect(s, k)).subs(sub).doit()*sp.exp(a*x/c)))
        tread[(s, k)] = sp.solve(ln, Q[(s, k)].diff(t))[0]
psi = {s: Q[(s, 1)] + sp.I*Q[(s, 2)] for s in (1, -1)}
dpsi = {s: tread[(s, 1)] + sp.I*tread[(s, 2)] for s in (1, -1)}
dirac = {s: -s*c*(psi[s].diff(x) - sp.I*e*Ax/hb*psi[s]) + s*a*psi[-s] - sp.I*e*Phi/hb*psi[s] for s in (1, -1)}
ok("3. copies paired as psi = Q^1 + i Q^2: d_t psi_eps = -eps c (d_x - i e A_x/hbar) psi_eps + eps a psi_-eps - i (e Phi/hbar) psi_eps  (Dirac with the field)",
   all(sp.simplify(sp.expand(dpsi[s] - dirac[s])) == 0 for s in (1, -1)))
# 4 gauge
g = sp.exp(sp.I*e*chi/hb)
psi_g = {s: g*psi[s] for s in (1, -1)}
dirac_g = {s: -s*c*(psi_g[s].diff(x) - sp.I*e*(Ax + chi.diff(x))/hb*psi_g[s]) + s*a*psi_g[-s] - sp.I*e*(Phi - chi.diff(t))/hb*psi_g[s] for s in (1, -1)}
ok("4. Phi -> Phi - d_t chi, A_x -> A_x + d_x chi with the copies turned by e chi/hbar: the same equation",
   all(sp.simplify(sp.expand((g*psi[s]).diff(t).subs({Q[(s, 1)].diff(t): tread[(s, 1)], Q[(s, 2)].diff(t): tread[(s, 2)]}) - dirac_g[s])) == 0 for s in (1, -1)))
# 4' (added 10 Nov 2026) V.2 (H): the real lines read in t, and the copies swapped
ok("4a. (E): d_t Q_eps = -eps c d_x Q_eps + eps a Q_-eps + (e/hbar)(Phi - eps c A_x) Q'_eps, and the primed line with the opposite sign on the field",
   all(sp.simplify(tread[(s_, 1)] - (-s_*c*Q[(s_, 1)].diff(x) + s_*a*Q[(-s_, 1)] + e/hb*(Phi - s_*c*Ax)*Q[(s_, 2)])) == 0
       and sp.simplify(tread[(s_, 2)] - (-s_*c*Q[(s_, 2)].diff(x) + s_*a*Q[(-s_, 2)] - e/hb*(Phi - s_*c*Ax)*Q[(s_, 1)])) == 0 for s_ in (1, -1)))
swap = {Q[(s_, k_)]: Q[(s_, 3 - k_)] for s_ in (1, -1) for k_ in (1, 2)}
ok("4b. (H): the copies swapped is the charge reversed: the lines for (Q', Q) are the lines for (Q, Q') with e -> -e; the free part does not see the swap",
   all(sp.simplify(tread[(s_, k_)].subs(swap, simultaneous=True) - tread[(s_, 3 - k_)].subs(e, -e)) == 0 for s_ in (1, -1) for k_ in (1, 2)))
# 5 lattice run (c = dx = dt = 1, a -> p per step), potential step Phi = V0 for x > 0, A_x = 0, read in t
p_, V0, N, T_, hbv = 0.08, 0.35, 300, 140, 1.0
xs = np.arange(-N, N + 1); Phiv = np.where(xs > 0, V0, 0.0)
def walk_step(Qp1, Qp2, Qm1, Qm2):                    # t-clock version of the same step: move, flip-turn, phase turn
    Qp1, Qp2 = np.roll(Qp1, 1), np.roll(Qp2, 1); Qm1, Qm2 = np.roll(Qm1, -1), np.roll(Qm2, -1)
    th = np.arcsin(p_)                                 # the exchange turn (loss out): (Q+, Q-) -> rotation by th
    Qp1, Qm1 = np.cos(th)*Qp1 + np.sin(th)*Qm1, -np.sin(th)*Qp1 + np.cos(th)*Qm1
    Qp2, Qm2 = np.cos(th)*Qp2 + np.sin(th)*Qm2, -np.sin(th)*Qp2 + np.cos(th)*Qm2
    ph = -Phiv/hbv                                     # turn between copies: angle -e Phi dt/hbar, same for both senses in the t-clock
    Qp1, Qp2 = np.cos(ph)*Qp1 - np.sin(ph)*Qp2, np.sin(ph)*Qp1 + np.cos(ph)*Qp2
    Qm1, Qm2 = np.cos(ph)*Qm1 - np.sin(ph)*Qm2, np.sin(ph)*Qm1 + np.cos(ph)*Qm2
    return Qp1, Qp2, Qm1, Qm2
g0 = np.exp(-((xs + 60)/12.0)**2)*np.cos(0.6*xs); g1 = np.exp(-((xs + 60)/12.0)**2)*np.sin(0.6*xs)
Qp1, Qp2, Qm1, Qm2 = g0.copy(), g1.copy(), 0.3*g1, -0.3*g0
psi_p, psi_m = Qp1 + 1j*Qp2, Qm1 + 1j*Qm2
norm0 = (np.abs(psi_p)**2 + np.abs(psi_m)**2).sum()
for _ in range(T_):
    Qp1, Qp2, Qm1, Qm2 = walk_step(Qp1, Qp2, Qm1, Qm2)
    psi_p, psi_m = np.roll(psi_p, 1), np.roll(psi_m, -1)                     # exact complex evolution, same splitting
    th = np.arcsin(p_); psi_p, psi_m = np.cos(th)*psi_p + np.sin(th)*psi_m, -np.sin(th)*psi_p + np.cos(th)*psi_m
    psi_p, psi_m = psi_p*np.exp(-1j*Phiv/hbv), psi_m*np.exp(-1j*Phiv/hbv)
ok(f"5a. lattice: two real copies vs one complex column, after {T_} ticks: max |difference| = {max(np.abs(Qp1 + 1j*Qp2 - psi_p).max(), np.abs(Qm1 + 1j*Qm2 - psi_m).max()):.2e}",
   np.abs(Qp1 + 1j*Qp2 - psi_p).max() < 1e-10)
ok(f"5b. total length over both copies conserved: {(Qp1**2 + Qp2**2 + Qm1**2 + Qm2**2).sum()/norm0:.12f}", abs((Qp1**2 + Qp2**2 + Qm1**2 + Qm2**2).sum()/norm0 - 1) < 1e-10)
refl = (np.abs(psi_p[xs <= 0])**2 + np.abs(psi_m[xs <= 0])**2).sum()/norm0
print(f"     potential step V0 = {V0}: fraction left of the step after the run = {refl:.3f} (partly reflected; a free packet would have crossed)")
ok(f"5c. the potential step reflects part of the packet: fraction left of the step {refl:.3f}, strictly between 0 and 1", 0.01 < refl < 0.99)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
