"""Loop closure at the path level (1+1 real walk, t as clock, lattice).
Step: each sense moves its way by one site, then the coin R = cos(a dt) 1 + sin(a dt) T acts (T = [[0,1],[-1,0]]), real, orthogonal.
A path = a sequence of N choices (keep sense / flip). Signed weight: cos per keep, +sin for a flip - -> +, -sin for a flip + -> -.
Out leg: source (x=0, sense +) at t=0  ->  (x, eps) at t=N.   Back leg: the walk run backward: (x, eps) at t=N -> source at t=0,
with the backward step R^T = cos 1 - sin T (the inverse of an orthogonal step: same walk, flip signs reversed, shifts reversed).
A loop = an out leg followed by a back leg that lands on the source with the source's sense.  Closure weight of a loop = product of the two leg weights.
Claims:
 1. sum over ALL loops through (x, eps) of the closure weight = Q_eps(x)^2  (Q from the matrix evolution)
 2. sum over eps = Q+^2 + Q-^2 = rho(x): the count of IV.4d (B) is the total closure weight through x
 3. the loops that come back along the same path (q = p) give sum_p w_p^2, which is NOT Q^2: the loops that explore another way back are the interference
 4. sum over x of the closure weight = 1: every loop closes somewhere (orthogonality)
 5. none of this needs a second copy: one real column, one traversal out, one back."""
import numpy as np, itertools
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
N, a, dt = 8, 1.0, 0.35
c_, s_ = np.cos(a*dt), np.sin(a*dt)
T = np.array([[0, 1], [-1, 0]]); R = c_*np.eye(2) + s_*T; Rb = R.T
# ---- matrix evolution: Q_eps(x) after N steps from (x=0, +)
L = 2*N + 1; Q = np.zeros((2, L)); Q[0, N] = 1.0
for _ in range(N):
    Q = R @ np.vstack([np.roll(Q[0], 1), np.roll(Q[1], -1)])
# ---- path enumeration for the out leg: dict (x, eps) -> list of signed weights
def legs(step_matrix, start_x, start_eps, nsteps, direction):
    """direction=+1: forward in t (move then coin). direction=-1: backward in t (undo coin, then undo move)."""
    out = {}
    for flips in itertools.product([0, 1], repeat=nsteps):
        x, e, w = start_x, start_eps, 1.0
        for f in flips:
            if direction == +1:
                x += e                                   # + moves right, - moves left
                if f: w *= (step_matrix[1, 0] if e == 1 else step_matrix[0, 1]); e = -e     # + -> - : R[1,0] = -sin ; - -> + : R[0,1] = +sin
                else: w *= step_matrix[0, 0]
            else:
                if f: w *= (step_matrix[1, 0] if e == 1 else step_matrix[0, 1]); e = -e
                else: w *= step_matrix[0, 0]
                x -= e                                   # undo the move of the sense we are now in
        out.setdefault((x, e), []).append(w)
    return out
out = legs(R, 0, 1, N, +1)
# check the enumeration reproduces the matrix column
ok("0. path sum of the out leg = matrix column Q_eps(x) (max err %.1e)" % max(abs(sum(out.get((x - N, e), [0])) - Q[0 if e == 1 else 1, x]) for x in range(L) for e in (1, -1)),
   max(abs(sum(out.get((x - N, e), [0])) - Q[0 if e == 1 else 1, x]) for x in range(L) for e in (1, -1)) < 1e-12)
# ---- back legs from each (x, eps): run the backward walk, keep only those landing on the source (0, +)
closure, diag, total = {}, {}, 0.0
for (x, e), ws in out.items():
    back = legs(Rb, x, e, N, -1).get((0, 1), [])
    closure[(x, e)] = sum(ws)*sum(back)                 # sum over all (out, back) pairs of w_p * w_q
    diag[(x, e)] = sum(w*w for w in ws)                 # loops that retrace the same path (q = p, weight w_p^2)
    total += closure[(x, e)]
err1 = max(abs(closure[(x, e)] - Q[0 if e == 1 else 1, x + N]**2) for (x, e) in closure)
ok(f"1. closure weight through (x, eps) = Q_eps(x)^2 for every point and sense (max err {err1:.1e})", err1 < 1e-12)
rho = {x: closure.get((x, 1), 0) + closure.get((x, -1), 0) for x in range(-N, N + 1)}
ok("2. summed over the sense: closure through x = Q+^2 + Q-^2 = rho(x)", all(abs(rho[x] - (Q[0, x + N]**2 + Q[1, x + N]**2)) < 1e-12 for x in rho))
xm = max(closure, key=lambda k: abs(closure[k] - diag[k]))
ok(f"3. same-path loops alone give sum_p w_p^2, not Q^2: at {xm} closure {closure[xm]:+.5f} vs same-path {diag[xm]:+.5f}; the rest is loops that come back another way",
   abs(closure[xm] - diag[xm]) > 1e-3)
ok(f"4. every loop closes somewhere: total closure weight over x = {total:.12f}", abs(total - 1) < 1e-12)
print(f"      (remark, not a check: the enumeration used one real column, {len(list(itertools.product([0, 1], repeat=N)))} paths out and as many back; no second copy and no complex number appear above)")
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
