"""The line-by-line content of III.15a (C), checked as matrix identities on the lattice.
U = the full N-tick evolution matrix on the space of (site, sense), real.
 1. U is orthogonal: U^T U = 1 (each tick is a shift composed with a rotation).
 2. outbound sum from A to X = U[X, A] (the sum over paths of the product of the selected matrix elements).
 3. return sum from X to A along the backward walk = (U^-1)[A, X] = (U^T)[A, X] = U[X, A]: the SAME number, by transposition. No 'same equation' argument needed.
 4. the backward step read with the senses swapped is the forward step: Ax (cos 1 - sin T) Ax = cos 1 + sin T.
 5. the double sum factorizes by distributivity alone: sum_ij w_i w_j = (sum_i w_i)(sum_j w_j); checked on the enumerated weights.
 6. event weight = U[X,A]^2; summed over X: sum_X U[X,A]^2 = (U^T U)[A,A] = 1.
 7. (C)4's last sentence: the same with a different real orthogonal matrix at every tick (coin angle varying from site to site and tick to tick)
    and any real start column f: W_X = (e_X^T U f)(f^T U^T e_X) = q_X^2, sum_X W_X = |f|^2.   (Check 7 added 10 Nov 2026.)"""
import numpy as np, itertools
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
N, a, dt = 6, 1.0, 0.35
c_, s_ = np.cos(a*dt), np.sin(a*dt)
T = np.array([[0, 1], [-1, 0]]); R = c_*np.eye(2) + s_*T; Ax = np.array([[0, 1], [1, 0]])
L = 2*N + 1; D = 2*L                       # state index: (sense, site) -> sense*L + site
def tick_matrix():
    M = np.zeros((D, D))
    for x in range(L):
        for e_i, e in enumerate((1, -1)):          # e_i = 0 for +, 1 for -
            xn = (x + e) % L                       # periodic sites: the step is exactly orthogonal; N < L/2 so no path wraps
            for f_i in range(2):                   # new sense index after the coin
                M[f_i*L + xn, e_i*L + x] += R[f_i, e_i]
    return M
Ut = tick_matrix(); U = np.linalg.matrix_power(Ut, N)
inner = U.T @ U
# orthogonality holds away from the edges (no paths reach the edge in N ticks from the center)
A = 0*L + N                                        # source: sense +, site N (center)
ok("1. U^T U = 1 (orthogonal step)", np.allclose(inner, np.eye(D)))
# 2. outbound path sum = U[X, A], by enumeration
def paths(nsteps):
    out = {}
    for flips in itertools.product([0, 1], repeat=nsteps):
        x, e, w = N, 1, 1.0
        for f in flips:
            x += e
            if f: w *= (R[1, 0] if e == 1 else R[0, 1]); e = -e
            else: w *= R[0, 0]
        out.setdefault((x, e), []).append(w)
    return out
out = paths(N)
ok("2. sum over outbound paths of the product of the selected matrix elements = U[X, A]",
   all(abs(sum(ws) - U[(0 if e == 1 else 1)*L + x, A]) < 1e-12 for (x, e), ws in out.items()))
# 3. return sum = (U^-1)[A, X] = (U^T)[A, X] = U[X, A]
Uinv = np.linalg.inv(U)
ok("3. (U^-1)[A, X] = (U^T)[A, X] = U[X, A] for every X: the return sum is the outbound sum, by transposition",
   np.allclose(Uinv[A, :], U.T[A, :]) and np.allclose(U.T[A, :], U[:, A]))
ok("4. Ax (cos 1 - sin T) Ax = cos 1 + sin T: the backward coin with the senses swapped is the forward coin", np.allclose(Ax @ R.T @ Ax, R))
# 5. distributivity on the enumerated weights
X0 = max(out, key=lambda k: len(out[k])); ws = out[X0]; back = ws[:]   # return weights: same multiset by 3 (checked in loop_closure.py)
ok("5. sum_ij w_i w_j = (sum_i w_i)(sum_j w_j) on the enumerated weights (distributivity, nothing else)",
   abs(sum(wi*wj for wi in ws for wj in back) - sum(ws)*sum(back)) < 1e-12)
ok(f"6. sum_X U[X,A]^2 = (U^T U)[A,A] = {np.sum(U[:, A]**2):.12f}", abs(np.sum(U[:, A]**2) - 1) < 1e-12)
# 7. tick matrices that vary from tick to tick and from site to site, any real start column
rng = np.random.default_rng(3)
def tick_var():
    M = np.zeros((D, D)); th = rng.uniform(-1.2, 1.2, size=L)             # a coin angle per site, drawn anew for this tick
    for x in range(L):
        Rx = np.array([[np.cos(th[(x + 1) % L]), np.sin(th[(x + 1) % L])], [-np.sin(th[(x + 1) % L]), np.cos(th[(x + 1) % L])]])
        Rxm = np.array([[np.cos(th[(x - 1) % L]), np.sin(th[(x - 1) % L])], [-np.sin(th[(x - 1) % L]), np.cos(th[(x - 1) % L])]])
        for f_i in range(2):
            M[f_i*L + (x + 1) % L, 0*L + x] += Rx[f_i, 0]                  # + moves right, then the coin of the site it lands on
            M[f_i*L + (x - 1) % L, 1*L + x] += Rxm[f_i, 1]                 # - moves left, then the coin of the site it lands on
    return M
ticks = [tick_var() for _ in range(9)]
Uv = np.eye(D)
for M in ticks: Uv = M @ Uv
f = rng.normal(size=D); q = Uv @ f
W = np.array([(np.eye(D)[X] @ Uv @ f)*(f @ Uv.T @ np.eye(D)[X]) for X in range(D)])
ok(f"7. nine different orthogonal ticks, random real start column: each tick has M^T M = 1, W_X = (e_X^T U f)(f^T U^T e_X) = q_X^2 at every X "
   f"(max diff {np.max(np.abs(W - q**2)):.1e}), sum_X W_X = {W.sum():.12f} = |f|^2 = {f @ f:.12f}",
   all(np.allclose(M.T @ M, np.eye(D)) for M in ticks) and np.max(np.abs(W - q**2)) < 1e-12 and abs(W.sum() - f @ f) < 1e-10 and np.all(W >= 0))
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
