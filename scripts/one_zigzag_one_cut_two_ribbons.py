"""
Exact count of the paths of III.1 (clock x, Delta x = c Delta t = 1) for the hand example of
drafts/one_zigzag_one_cut_two_ribbons_by_hand.md. Nothing is solved or approximated: each path's weight is the product of its
coin factors with the loss out (III.5a (F): keep 1, flip beta), crossings are read on the cut t = 3/2, and the two-walker
column through one coin is formed as VI.1 (D) writes it, C -> M C M^T, with C antisymmetric.
"""
import itertools, sys
import numpy as np
def ok(msg, cond): print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()
beta = 0.3
def run(senses, x0, t0):                     # senses: the sense carried on each step; coin after each step sets the next
    x, t, pts = x0, t0, [(x0, t0)]
    for e in senses:
        x += 1; t += e; pts.append((x, t))
    flips = sum(1 for i in range(len(senses) - 1) if senses[i] != senses[i + 1]) + 0   # coins between steps
    return pts, flips
def weight(senses, arrival):                 # arrival: the sense read at the endpoint (one more coin after the last step)
    s = list(senses) + [arrival]
    return beta**sum(1 for i in range(len(s) - 1) if s[i] != s[i + 1])
def crossings(pts, tc):                      # (x at the cut, charge) for each step that crosses t = tc
    out = []
    for (x0, t0), (x1, t1) in zip(pts, pts[1:]):
        if (t0 - tc)*(t1 - tc) < 0:
            out.append((x0 + (tc - t0)/(t1 - t0), +1 if t1 > t0 else -1))
    return out

# 1. path 1
s1 = (+1, +1, -1, +1, +1); pts1, _ = run(s1, 0, 0); w1 = weight(s1, +1); cr1 = crossings(pts1, 1.5)
ok("path 1 %s ends at %s, weight %.4f = beta^2 = %.4f; crossings on t = 3/2 at x = %s with charges %s"
   % (s1, pts1[-1], w1, beta**2, [c[0] for c in cr1], [c[1] for c in cr1]),
   pts1[-1] == (5, 3) and abs(w1 - beta**2) < 1e-15 and [c for c in cr1] == [(1.5, 1), (2.5, -1), (3.5, 1)])

# 2. III.3 on every path from A: charge +1 on the cut, odd number of crossings; total weight with the loss out
tot = 0.0; good = True; ncross = 0
for rest in itertools.product((+1, -1), repeat=4):
    s = (+1,) + rest                          # A carries the forward sense on its first step
    pts, _ = run(s, 0, 0); cr = crossings(pts, 1.5)
    above = pts[-1][1] > 1.5                  # the path ends above the cut (its last crossing was upward)
    if above: ncross += 1; good &= (sum(c[1] for c in cr) == 1) and (len(cr) % 2 == 1)
    else: good &= (sum(c[1] for c in cr) == 0) and (len(cr) % 2 == 0)    # ends below: as many down as up crossings
    tot += weight(s, +1) + weight(s, -1)
ok("all 16 paths of III.1 from A (first sense +): the %d that end above the cut carry charge +1 on it with an odd number of crossings, "
   "the rest charge 0 with an even number (III.3: the charge on a cut is the charge that has passed it); sum of weights over paths and "
   "arrival senses = %.6f = (1 + beta)^5, five coins each (1 + beta) (the loss-out count grows with x, III.5a (J))" % (ncross, tot),
   good and abs(tot - (1 + beta)**5) < 1e-12)

# 3. the two-walker column through the coin at x = 3
w = beta*1.0                                 # path 1's first flip times ribbon 2's keep, before the coin at x = 3
I, J = 1, 0                                  # states at t = 1: index 0 = (1,+), 1 = (1,-)
C = np.zeros((2, 2)); C[I, J] = +w; C[J, I] = -w
M = np.array([[1.0, beta], [beta, 1.0]])     # rows new sense (+,-), columns old sense (+,-)
Cp = M @ C @ M.T
D = M[0, I]*C[I, J]*M[0, J]; E = M[0, J]*C[J, I]*M[0, I]
ok("two-walker column at x = 3: same-state entry C'_{(1,+),(1,+)} = %.1e, from D = %+.4f (path 1 flips into +, ribbon 2 keeps +) "
   "and E = %+.4f (the slots exchanged); surviving entries C'_{(1,-),(1,+)} = %+.4f = (1 - beta^2) w = %+.4f, "
   "C'_{(1,+),(1,-)} = %+.4f; C' antisymmetric to %.0e"
   % (Cp[0, 0], D, E, Cp[1, 0], (1 - beta**2)*w, Cp[0, 1], abs(Cp + Cp.T).max()),
   abs(Cp[0, 0]) < 1e-15 and abs(D - beta*w) < 1e-15 and abs(E + beta*w) < 1e-15 and abs(Cp[1, 0] - (1 - beta**2)*w) < 1e-15
   and abs(Cp + Cp.T).max() < 1e-15)
