"""III.15c (B), (C), (G): the counter row, the symmetric part of the generator, and the coincidence meaning of the square.
(B) a conserved pairing g.Q needs g(n) = U^T g(n+1); the stochastic coin has the fixed row 1, the rotation has none
    (left eigenvectors (1, +-i), eigenvalues e^{-+ i theta}); U^T = U^{-1} makes the moving row the walk read backward.
(C) d/dt sum Q^2 = Q.(L + L^T) Q : only the symmetric part of the generator moves a sum of squares.
    P generator: symmetric part a sigma_x - a, eigenvalues 0, -2a  -> sum P^2 decreases.
    loss out, clock x: symmetric part a sigma_x -> sum Q^2 not conserved; int (Q+^2 - Q-^2) dt conserved.
    same lines, clock t: symmetric part 0 -> sum Q^2 conserved.
(G) sum_x P(x,N)^2 = Pr(two independent copies coincide at N) = agreement of the walk with its come-back from where it arrived;
    decays for the coin (walk forgets), = 1 for the rotation (walk retraces)."""
import numpy as np
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
rng = np.random.default_rng(0)
L = 161; idx = lambda x, s: 2*(x % L) + s
S = np.zeros((2*L, 2*L))
for x in range(L): S[idx(x+1, 0), idx(x, 0)] = 1; S[idx(x-1, 1), idx(x, 1)] = 1
p, th = 0.2, 0.45
Mcoin = np.kron(np.eye(L), np.array([[1-p, p], [p, 1-p]])) @ S
R2 = np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]]); U = np.kron(np.eye(L), R2) @ S
one = np.ones(2*L)
# (B)
ok("B1. stochastic coin: 1^T M = 1^T (fixed counter row)", np.allclose(one @ Mcoin, one))
w, v = np.linalg.eig(R2.T)
ok("B2. rotation coin: left eigenvectors (1,+-i), eigenvalues e^{-+ i theta}; no real fixed row",
   np.allclose(sorted(w, key=np.angle), sorted([np.exp(-1j*th), np.exp(1j*th)], key=np.angle)) and not np.allclose(one[:2] @ R2, one[:2]))
e = np.zeros(2*L); e[idx(L//2, 0)] = 1; N = 40
Q = e.copy(); Qs = [Q.copy()]
for n in range(N): Q = U @ Q; Qs.append(Q.copy())
g = rng.normal(size=2*L); gs = [None]*(N+1); gs[N] = g
for n in range(N-1, -1, -1): gs[n] = U.T @ gs[n+1]
pair = [gs[n] @ Qs[n] for n in range(N+1)]
ok("B3. g(n) = U^T g(n+1): the pairing g.Q is the same on every cut", np.ptp(pair) < 1e-12)
ok("B4. U^T = U^{-1}: the moving row from the endpoint is the walk read backward; pairing at both ends = U_XA",
   np.allclose(U.T, np.linalg.inv(U)))
# (C) generators on the two-sense column (k-space, so the shift is a number): L_P, L_x (loss out, clock x), L_t (clock t)
k = 0.7; a = 1.0
A = np.diag([1., -1.]); sx = np.array([[0, 1], [1, 0.]]); T = np.array([[0, 1], [-1, 0.]])
L_P = -1j*k*A + a*sx - a*np.eye(2)          # d_t P = L_P P  (telegraph coin, t clock)
L_x = -1j*k*A + a*sx                        # loss out, read in the clock x (d_t -> i k on the transverse variable)
L_t = -1j*k*A + a*T                         # the same lines read in t (III.5a (J))
sym = lambda M: (M + M.conj().T)/2
ok("C1. P generator: symmetric part a sigma_x - a, eigenvalues 0 and -2a (sum P^2 can only decrease)",
   np.allclose(sorted(np.linalg.eigvalsh(sym(L_P))), [-2*a, 0]))
ok("C2. loss out, clock x: symmetric part a sigma_x (sum Q^2 not conserved); G = diag(1,-1): G L_x antisymmetric (int Q+^2 - Q-^2 conserved)",
   np.allclose(sym(L_x), a*sx) and np.allclose(sym(np.diag([1, -1]) @ L_x), 0))
ok("C3. same lines, clock t: symmetric part 0 (sum Q^2 conserved); the exchange becomes the turn by one multiplication, A sigma_x = T (III.5a (J))",
   np.allclose(sym(L_t), 0) and np.allclose(A @ sx, T))
# (G) coincidence
def K(M):
    P = e.copy(); out = []
    for n in range(61): out.append((P**2).sum()); P = M @ P
    return np.array(out)
Kc, Kr = K(Mcoin), K(U)
Rn = 100000
def mc(Nn):
    x = np.zeros((2, Rn), int); s = np.ones((2, Rn), int)
    for n in range(Nn):
        x += s; f = rng.random((2, Rn)) < p; s = np.where(f, -s, s)
    return np.mean((x[0] == x[1]) & (s[0] == s[1]))
ok("G1. coin: sum P(N)^2 = Pr(two independent copies coincide at N), MC at N=10,20: %.4f/%.4f vs %.4f/%.4f" % (mc(10), mc(20), Kc[10], Kc[20]),
   abs(mc(10) - Kc[10]) < 4e-3 and abs(mc(20) - Kc[20]) < 4e-3)
ok("G2. coin: it decays (1 -> %.3f at N=60); rotation: sum Q(N)^2 = 1 at every N" % Kc[60], Kc[60] < 0.03 and np.allclose(Kr, 1))
Ps = [e.copy()]
for n in range(N): Ps.append(Mcoin @ Ps[-1])
Pt = [None]*(N+1); Pt[N] = Ps[N].copy()
for n in range(N-1, -1, -1): Pt[n] = Mcoin.T @ Pt[n+1]
ag = [Ps[n] @ Pt[n] for n in range(N+1)]
ok("G3. come-back launched from the arrived distribution: agreement on every cut = sum P(N)^2", np.ptp(ag) < 1e-12 and abs(ag[0] - (Ps[N]**2).sum()) < 1e-12)
