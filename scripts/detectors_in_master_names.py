"""Checks for the section drafts of 12 Nov 2026 (step 4): the detector reading written in the master's own names, u_1..u_4 (the original four
detectors, IV.4b's counters) and v_1..v_4 (the added four), order r = (u_1, u_2, u_3, u_4, v_1, v_2, v_3, v_4); and in one dimension
r = (u_+, u_-, v_+, v_-) (III.15e).  No relabelling of the master's equations is made; new variables only for the sums and differences.

  1. the table: which detector (time sense, orientation) each u_j, v_j is; v_j has u_j's time sense and u_j's orientation turned in the walker's own plane.
  2. the exchange of the forward and the backward detector of one orientation is J T (IV.5a's product), with square +1; at rest v = +-T u <=> J T r = +-r.
  3. sums and differences from u, v (the matrix), the record, the flow, alpha = -u.Tv, gamma = u.(T w)v in them.
  4. the copy swap K = diag(1,-1) carries sums onto differences; the commutant of the free step is 1, J, K, JK; turns that keep the record: +-J.
  5. one dimension: detectors, sums and differences, the two loops, W = -alpha, K_ = gamma-like; the clock hand: rho_eps D_eps theta_eps = a W - b_eps rho_eps.
  6. the effective flip rate at rest on (u, v): sums -(a + b) J, differences +(a - b) J.
  7. the clock hand in motion: W = -1, rho_+- = e^{+-eta} for one clock run to rapidity eta; the other clock has W = +1.
  8. V.4 as drafted: the 4x4 matrices of one dimension, the split of the potential's term, and the frame form.
"""
import sys, itertools
import numpy as np, sympy as sp
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all = ok_all and bool(cond); print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()
rng = np.random.default_rng(12)
Ax = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]], float); Ay = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]], float)
Az = np.diag([1., -1, -1, 1]); T = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]], float); I4 = np.eye(4); Z4 = np.zeros((4, 4))
Ryz, Rzx, Rxy, w = Ay@Az, Az@Ax, Ax@Ay, Ax@Ay@Az
u0 = np.array([1., 0, 0, 0]); v0 = T@u0
AA = [np.kron(np.eye(2), A_) for A_ in (Ax, Ay, Az)]; TT = np.kron(np.eye(2), T); J = np.block([[Z4, -I4], [I4, Z4]]); K = np.diag([1., 1, 1, 1, -1, -1, -1, -1])
names = ['u1', 'u2', 'u3', 'u4', 'v1', 'v2', 'v3', 'v4']
# ---- 1. the table
orient = {'as launched': I4, 'turned in (y,z)': Ryz, 'turned in (z,x)': Rzx, 'turned in (x,y)': Rxy}
table = {}
for nm, O in orient.items():
    for src, lab in ((u0, 'u'), (v0, 'v')):
        col = O@src; i = int(np.argmax(np.abs(col))); assert abs(abs(col[i]) - 1) < 1e-12
        table[lab + str(i + 1)] = (('forward' if i in (0, 3) else 'backward'), nm, int(np.sign(col[i])), lab)
# orientation of v_j: the walker with orientation O turned in its own plane is O R_yz; applied to u0 it lands on a v-slot: O R_yz u0 = -O v0
own = all(abs(abs((O@Ryz@u0)@(O@v0)) - 1) < 1e-12 for O in orient.values())
same_sense = all(table['u%d' % j][0] == table['v%d' % j][0] for j in (1, 2, 3, 4))
quat = lambda O: [k_ for k_, P in orient.items() if abs(abs(np.trace(P.T@O))/4 - 1) < 1e-12][0]
turned = all(quat(orient[table['u%d' % j][1]]@Ryz) == table['v%d' % j][1] for j in (1, 2, 3, 4))
print("      " + "; ".join("%s: %s, %s (%+d)" % (k_, *table[k_][:3]) for k_ in names))
ok("1. the eight counters as detectors: each u_j, v_j is one time sense and one orientation (printed, with the sign of the counter's convention). "
   "Counters 1 and 4 are forward and 2 and 3 backward in both sets (%s); v_j has the orientation of u_j turned in the walker's own plane, i.e. "
   "multiplied by R_yz on the right (%s, %s)" % (same_sense, turned, own), same_sense and turned and own and len(table) == 8)
# ---- 2. the exchange of forward and backward of one orientation
X = J@TT
pairs = {}
for k_, (sense, nm, sg, lab) in table.items(): pairs.setdefault(nm, {})[sense] = (k_, sg)
perm_ok = True
for nm, d_ in pairs.items():
    (kf, sf), (kb, sb) = d_['forward'], d_['backward']; ef = np.eye(8)[names.index(kf)]; eb = np.eye(8)[names.index(kb)]
    perm_ok = perm_ok and np.abs(X@(sf*ef) - sb*eb).max() < 1e-12 and np.abs(X@(sb*eb) - sf*ef).max() < 1e-12
rest = 0.0
for _ in range(20):
    u = rng.normal(size=4)
    for sg in (1, -1):
        r = np.concatenate([u, sg*T@u]); rest = max(rest, np.abs(X@r - sg*r).max())
ok("2. the exchange of the forward and the backward detector of one orientation is J T = [[0, -T], [T, 0]]: it carries each forward detector onto the "
   "backward detector of the same orientation and back, with the counters' signs (%s); (J T)^2 = +1 (%.0e); J T commutes with J and with T, and "
   "anticommutes with the three direction matrices (%.0e); v = +T u is J T r = +r and v = -T u is J T r = -r (%.0e). Pairs of one orientation: %s"
   % (perm_ok, np.abs(X@X - np.eye(8)).max(), max(np.abs(X@J - J@X).max(), np.abs(X@TT - TT@X).max(), max(np.abs(X@A_ + A_@X).max() for A_ in AA)), rest,
      "; ".join("%s: (%s, %s)" % (nm, d_['forward'][0], d_['backward'][0]) for nm, d_ in pairs.items())),
   perm_ok and np.abs(X@X - np.eye(8)).max() < 1e-12 and rest < 1e-12 and max(np.abs(X@A_ + A_@X).max() for A_ in AA) < 1e-12)
# ---- 3. sums and differences from (u, v)
basis = [I4, Ryz, Rzx, Rxy, w, Ax, Ay, Az]                                    # s_t, s_x, s_y, s_z | d_t, d_x, d_y, d_z
Phi = np.array([np.concatenate([E_@u0, E_@v0]) for E_ in basis]).T             # r = Phi n
Sm = np.linalg.inv(Phi)                                                        # n = Sm r
half = np.abs(np.abs(Sm[np.abs(Sm) > 1e-12]) - 0.5).max() < 1e-12 and np.abs(Phi.T@Phi - 2*np.eye(8)).max() < 1e-12
nl = ['s_t', 's_x', 's_y', 's_z', 'd_t', 'd_x', 'd_y', 'd_z']
for i in range(8): print("      %s = " % nl[i] + " ".join("%s1/2 %s" % ('+' if Sm[i, j] > 0 else '-', names[j]) for j in range(8) if abs(Sm[i, j]) > 1e-12))
Xn = Sm@X@Phi; d3 = 0.0
for _ in range(40):
    r = rng.normal(size=8); u, v = r[:4], r[4:]; n = Sm@r; st, s, dt, d = n[0], n[1:4], n[4], n[5:8]
    rho = r@r; j = np.array([r@A_@r for A_ in AA]); al = -(u@T@v); ga = u@T@w@v
    d3 = max(d3, abs(rho - 2*n@n), np.abs(j - 4*(st*d + dt*s - np.cross(s, d))).max(), abs(al - (st**2 + s@s - dt**2 - d@d)), abs(ga - 2*(st*dt - s@d)),
             abs(rho**2 - j@j - 4*(al**2 + ga**2)))
ok("3. sums and differences from the counters (printed): each is half the sum or half the difference of the two counters of one orientation, with the "
   "counters' signs (%s); J T keeps the four sums and reverses the four differences (%.0e); record = 2 (sum of the eight squares), flow = "
   "4c (s_t d + d_t s - s x d), alpha = -u.Tv = (sums squared) - (differences squared), gamma = u.(T w)v = 2 (s_t d_t - s.d), record^2 - flow^2/c^2 = "
   "4 (alpha^2 + gamma^2) (40 random, %.0e)" % (half, np.abs(Xn - np.diag([1., 1, 1, 1, -1, -1, -1, -1])).max(), d3),
   half and np.abs(Xn - np.diag([1., 1, 1, 1, -1, -1, -1, -1])).max() < 1e-12 and d3 < 1e-9)
# ---- 4. the copy swap, and the commutant
Kn = np.round(Sm@K@Phi); kswap = np.abs(Kn[:4, :4]).max() + np.abs(Kn[4:, 4:]).max()
gens = AA + [TT]; Msym = sp.Matrix(8, 8, sp.symbols('c0:64')); eqs = []
for G in gens:
    Gs = sp.Matrix(G.astype(int)); eqs += list(Msym*Gs - Gs*Msym)
sol = sp.solve(eqs, list(Msym), dict=True)[0]; comm_dim = len(Msym.subs(sol).free_symbols)
span = np.array([M_.ravel() for M_ in (np.eye(8), J, K, J@K)]); in_span = all(np.abs(M_@G - G@M_).max() < 1e-12 for M_ in (np.eye(8), J, K, J@K) for G in gens)
ok("4. counting the added four with the other sign, K = diag(1,1,1,1,-1,-1,-1,-1), carries every sum onto a difference and every difference onto a sum "
   "(%.0e); K J = -J K (%.0e). The matrices commuting with A_x, A_y, A_z and T on both sets are %d: 1, J, K, J K (%s)"
   % (kswap, np.abs(K@J + J@K).max(), comm_dim, in_span), kswap < 1e-12 and comm_dim == 4 and in_span)
# ---- 5. one dimension
a, b = sp.symbols('a b', real=True); up, um, vp, vm = [sp.Function(n_)(sp.Symbol('t'), sp.Symbol('x')) for n_ in ('u_p', 'u_m', 'v_p', 'v_m')]
t, x = sp.Symbol('t'), sp.Symbol('x'); c = 1
# III.8c (10a): D_+ r_+ = a r_- - b_+ J r_+, D_- r_- = -a r_+ - b_- J r_-, with J r = (-v, u); take b_+ = b_- = b (scalar potential)
rules = {sp.Derivative(up, t): -sp.Derivative(up, x) + a*um + b*vp, sp.Derivative(vp, t): -sp.Derivative(vp, x) + a*vm - b*up,
         sp.Derivative(um, t): sp.Derivative(um, x) - a*up + b*vm, sp.Derivative(vm, t): sp.Derivative(vm, x) - a*vp - b*um}
Dp = lambda f: (sp.diff(f, t) + sp.diff(f, x)).subs(rules); Dm = lambda f: (sp.diff(f, t) - sp.diff(f, x)).subs(rules)
W = up*vm - vp*um; rp, rm = up**2 + vp**2, um**2 + vm**2
hand_p = sp.simplify(up*Dp(vp) - vp*Dp(up) - (a*W - b*rp)); hand_m = sp.simplify(um*Dm(vm) - vm*Dm(um) - (a*W - b*rm))
# detectors in one dimension: F0 = u+, Bq = u-, Fq = v+, B0 = -v-;  alpha = F0 B0 + Fq Bq
alpha1 = up*(-vm) + vp*um
lam = np.array([1., 0, 0, -1]); ell = np.array([1., 0, 0, 1])                  # III.15e's leading and lagging loops
S1 = 0.5*np.array([[1, 0, 0, -1], [0, 1, 1, 0], [1, 0, 0, 1], [0, -1, 1, 0]])  # (sigma_0, sigma_q, delta_0, delta_q) from (u+, u-, v+, v-)
ok("5. one dimension, r = (u_+, u_-, v_+, v_-): u_+ forward with the hand at the start, u_- backward with the hand a quarter on, v_+ forward a "
   "quarter on, -v_- backward at the start. Sums (u_+ - v_-)/2, (v_+ + u_-)/2 and differences (u_+ + v_-)/2, (v_+ - u_-)/2: the leading loop "
   "(1,0,0,-1) is sums only %s, the lagging loop (1,0,0,1) differences only %s. W = u_+ v_- - v_+ u_- is minus the overlap of the forward and the "
   "backward detectors of the same hand position (%s). The clock hand: rho_eps D_eps theta_eps = a W - b rho_eps for both senses, from III.8c (10a) "
   "(symbolically: %s, %s)" % ((S1@lam).tolist(), (S1@ell).tolist(), sp.simplify(W + alpha1) == 0, hand_p == 0, hand_m == 0),
   np.abs(S1@lam - np.array([1., 0, 0, 0])).max() < 1e-12 and np.abs(S1@ell - np.array([0, 0, 1., 0])).max() < 1e-12 and sp.simplify(W + alpha1) == 0
   and hand_p == 0 and hand_m == 0)
# ---- 6. the effective flip rate at rest, on (u, v)
aa, bb = 0.8, 0.37; G = aa*TT - bb*J                                           # V.2 (E) at rest: the field's term is b (Q', -Q) = -b J (Q, Q'), b = e Phi/hbar
Ps, Pd = (np.eye(8) + X)/2, (np.eye(8) - X)/2
e6 = max(np.abs(G@Ps + (aa + bb)*J@Ps).max(), np.abs(G@Pd - (aa - bb)*J@Pd).max(), np.abs(TT + J@X).max())
Gs = (aa + 0.21)*TT                                                            # a change of the coin, V.1 first row
e6 = max(e6, np.abs(Gs@Ps + (aa + 0.21)*J@Ps).max(), np.abs(Gs@Pd - (aa + 0.21)*J@Pd).max())
plus_sign = np.abs((aa*TT + bb*J)@Ps + (aa - bb)*J@Ps).max()                   # V.4 as printed, "a T + (e Phi/hbar) J": gives a - b on the sums
ok("6. at rest with a scalar potential the generator is a T - b J (V.2 (E): the field's term is b (Q', -Q) = -b J (Q, Q'), b = e Phi/hbar). With "
   "T = -J (J T) it is -(a + b) J on the sums and +(a - b) J on the differences; a change of the coin, (a + V) T, is -(a + V) J on the sums and "
   "+(a + V) J on the differences (%.0e). V.4 prints the generator as 'a T + (e Phi/hbar) J': with IV.5a's J that would give a - e Phi/hbar on the "
   "sums (%.0e) against V.4's own a_eff = a + e Phi/hbar; the sign of J there is a slip, to correct with the section" % (e6, plus_sign),
   e6 < 1e-12 and plus_sign < 1e-12)
# ---- 7. the clock hand of each time sense: at rest, the other clock, and in motion
eta = sp.Symbol('eta', real=True); th = a*(sp.cosh(eta)*t - sp.sinh(eta)*x)                       # the walker's own time, c = 1
def hand(sgn):
    U = sp.Matrix([sp.exp(eta/2)*sp.cos(th), -sp.exp(-eta/2)*sp.sin(th)]) if sgn == 1 else sp.Matrix([sp.exp(eta/2)*sp.cos(th), sp.exp(-eta/2)*sp.sin(th)])
    Tm = sp.Matrix([[0, 1], [-1, 0]]); B_ = sp.diag(sp.exp(eta/2), sp.exp(-eta/2))
    rest_u = sp.Matrix([sp.cos(th), -sp.sin(th)]) if sgn == 1 else sp.Matrix([sp.cos(th), sp.sin(th)])
    U = B_*rest_u; V = B_*(sgn*Tm*rest_u)
    L = lambda q: sp.Matrix([-sp.diff(q[0], x) + a*q[1], sp.diff(q[1], x) - a*q[0]])
    if sgn == -1:                                                                                  # the other clock runs the other way in its own time
        U = U.subs(th, -th); V = V.subs(th, -th)
    solves = sp.simplify(sp.diff(U, t) - L(U)) == sp.zeros(2, 1) and sp.simplify(sp.diff(V, t) - L(V)) == sp.zeros(2, 1)
    Wv = sp.simplify(U[0]*V[1] - V[0]*U[1]); rp = sp.simplify(U[0]**2 + V[0]**2); rm = sp.simplify(U[1]**2 + V[1]**2)
    return solves, Wv, rp, rm
s1, W1, rp1, rm1 = hand(1); s2, W2, rp2, rm2 = hand(-1)
ok("7. the clock hand. One clock run to rapidity eta, u = diag(e^{eta/2}, e^{-eta/2}) (cos th, -sin th), v = the same on T (cos th, -sin th), th = a (t cosh eta "
   "- x sinh eta), solves the lines (%s) and has W = %s, rho_+ = %s, rho_- = %s: by III.8c (8) the forward hand turns at a e^{-eta} and the backward one at "
   "a e^{+eta} along their light lines, the flip rate times sqrt(rho_-/rho_+) and sqrt(rho_+/rho_-). The other clock (v = -T u, th -> -th) solves the lines "
   "(%s) with W = %s: the hands turn the other way" % (s1, W1, rp1, rm1, s2, W2),
   s1 and s2 and sp.simplify(W1 + 1) == 0 and sp.simplify(W2 - 1) == 0 and sp.simplify(rp1 - sp.exp(eta)) == 0 and sp.simplify(rm1 - sp.exp(-eta)) == 0)

# ---- 8. V.4 as drafted: the 4x4 matrices of one dimension written out, the split (3), and the frame form (4), (5)
T2 = np.array([[0, 1.], [-1, 0]]); TT1 = np.kron(np.eye(2), T2); J1 = np.array([[0, 0, -1, 0], [0, 0, 0, -1], [1, 0, 0, 0], [0, 1, 0, 0.]])
TT1_printed = np.array([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0.]]); X1 = np.array([[0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0.]])
G1 = aa*TT1 - bb*J1; Ps1, Pd1 = (np.eye(4) + X1)/2, (np.eye(4) - X1)/2
# V.2 (E) at rest, line by line: dQ+ = a Q- + b Q'+, dQ- = -a Q+ + b Q'-, dQ'+ = a Q'- - b Q+, dQ'- = -a Q'+ - b Q-
G_lines = np.array([[0, aa, bb, 0], [-aa, 0, 0, bb], [-bb, 0, 0, aa], [0, -bb, -aa, 0]])
e8 = max(np.abs(TT1 - TT1_printed).max(), np.abs(J1@TT1 - X1).max(), np.abs(X1@X1 - np.eye(4)).max(), np.abs(G1 - G_lines).max(),
         np.abs(TT1 + J1@X1).max(), np.abs(G1@Ps1 + (aa + bb)*J1@Ps1).max(), np.abs(G1@Pd1 - (aa - bb)*J1@Pd1).max(),
         np.abs(bb*J1 - (bb*J1@X1 + bb*J1@(np.eye(4) - X1))).max(), np.abs((np.eye(4) - X1)@Ps1).max(), np.abs((np.eye(4) - X1)@Pd1 - 2*Pd1).max(),
         np.abs(X1@lam - lam).max(), np.abs(X1@ell + ell).max())
basis = [I4, Ryz, Rzx, Rxy, w, Ax, Ay, Az]                                    # order of n = (s_t, s_x, s_y, s_z | d_t, d_x, d_y, d_z)
def frame(n_): return sum(c*B_ for c, B_ in zip(n_, basis))
e8f = 0.0
for _ in range(20):
    r_ = rng.normal(size=8); n_ = Sm@r_; M_ = frame(n_); Mt = T@M_@np.linalg.inv(T)
    e8f = max(e8f, np.abs(frame(Sm@(J@r_)) - M_@Ryz).max(), np.abs(frame(Sm@(TT@r_)) + Mt@Ryz).max(), np.abs(frame(Sm@(X@r_)) - Mt).max(),
              np.abs((M_ - Mt) - 2*(n_[5]*Ax + n_[6]*Ay + n_[7]*Az + n_[4]*w)).max(),
              np.abs(frame(Sm@((aa*TT - bb*J)@r_)) - (-aa*Mt@Ryz - bb*M_@Ryz)).max())
ok("8. V.4 as drafted. One dimension, r = (Q+, Q-, Q'+, Q'-): the printed 4x4 matrices of the time flip and the walker's own turn, their product "
   "(the exchange Q+ <-> -Q'-, Q'+ <-> Q-, square +1), the generator a T - b J = V.2 (E)'s four lines at rest, -(a + b) J on the part the exchange "
   "keeps and +(a - b) J on the part it reverses, the split b J = b J (JT) + b J (1 - JT) with the second term on the reversed part only; the leading "
   "loop is kept and the lagging loop reversed (%.0e). On the frame (20 random walkers): the potential's term is -b M R_yz, the mass term -a M~ R_yz, "
   "and M - M~ = 2 (d_x A_x + d_y A_y + d_z A_z + d_t w) (%.0e)" % (e8, e8f), e8 < 1e-12 and e8f < 1e-12)

print("ALL PASS" if ok_all else "SOME FAILED"); raise SystemExit(0 if ok_all else 1)
