"""Where do the signs sit in the loop picture?  Same enumeration as loop_closure.py.
(a) backward-in-time crossings as negative probability?  -> the closure through (x, -) is Q_-^2 >= 0: a backward crossing closes with POSITIVE weight.
(b) the back leg as negative probability (linear)?       -> out-leg sum + (-back-leg sum) = Q - Q = 0 at every point: a detector that responds
    to each leg separately feels a push and an equal pull-back; nothing net. The closure (product) is Q^2. That is 'pushed then pulled back' of detector_as_force.py 3b.
(c) where the negative weights actually are: among loops (p, q) with q != p. Count them, sum them: negatives + positives = Q^2, and
    the negatives are exactly the loops whose two legs disagree in the parity of their U-turn signs."""
import numpy as np, itertools
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
N, a, dt = 8, 1.0, 0.35
c_, s_ = np.cos(a*dt), np.sin(a*dt)
T = np.array([[0, 1], [-1, 0]]); R = c_*np.eye(2) + s_*T; Rb = R.T
def legs(M, x0, e0, n, direction):
    """Returns {(x, e): [(weight, route)]}; a route is the tuple of states (x, e) visited, written in the order of increasing t for both directions,
    so an out leg and the back leg that retraces it have the same route."""
    out = {}
    for flips in itertools.product([0, 1], repeat=n):
        x, e, w = x0, e0, 1.0; route = [(x, e)]
        for f in flips:
            if direction == +1: x += e
            if f:
                w *= (M[1, 0] if e == 1 else M[0, 1]); e = -e
            else: w *= M[0, 0]
            if direction == -1: x -= e
            route.append((x, e))
        if direction == -1: route = route[::-1]
        out.setdefault((x, e), []).append((w, tuple(route)))
    return out
out = legs(R, 0, 1, N, +1)
L = 2*N + 1; Q = np.zeros((2, L)); Q[0, N] = 1.0
for _ in range(N): Q = R @ np.vstack([np.roll(Q[0], 1), np.roll(Q[1], -1)])
neg_total, pos_total, closure_minus_ok, linear_ok, nneg, nloops = 0.0, 0.0, True, True, 0, 0
n_retrace, retrace_pos, neg_all_differ, every_out_has_its_retrace = 0, True, True, True
for (x, e), legs_out in out.items():
    back = legs(Rb, x, e, N, -1).get((0, 1), [])
    ws = [w for w, _ in legs_out]
    S_out, S_back = sum(ws), sum(w for w, _ in back)
    if e == -1 and S_out*S_back < -1e-14: closure_minus_ok = False
    if abs(S_out - S_back) > 1e-12: linear_ok = False
    back_routes = {r: w for w, r in back}
    for wp, rp in legs_out:
        if rp not in back_routes or abs(back_routes[rp] - wp) > 1e-15: every_out_has_its_retrace = False
        for wq, rq in back:
            nloops += 1
            if rq == rp:
                n_retrace += 1
                if not wp*wq > 0: retrace_pos = False
            if wp*wq < 0:
                neg_total += wp*wq; nneg += 1
                if rq == rp: neg_all_differ = False
            else: pos_total += wp*wq
ok("(a) closure through (x, -), the backward crossings, is Q_-^2 >= 0 everywhere: a backward crossing is not a negative probability", closure_minus_ok)
ok("(b) linear: out-leg sum = back-leg sum = Q at every point, so 'forward +, backward -' cancels exactly (push, then equal pull-back); only the product survives", linear_ok)
tot = pos_total + neg_total; Qsq = (Q**2).sum()
ok(f"(c) of {nloops} loops, {nneg} have negative weight; positives {pos_total:+.4f} + negatives {neg_total:+.4f} = {tot:.6f} = sum of Q^2 = {Qsq:.6f}", abs(tot - Qsq) < 1e-10)
ok(f"(c') every out route has exactly its own retraced return among the back legs, with the same weight ({n_retrace} such loops, each w_p^2 > 0); "
   "every one of the negative loops comes back by a different route", every_out_has_its_retrace and n_retrace == 2**N and retrace_pos and neg_all_differ)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
