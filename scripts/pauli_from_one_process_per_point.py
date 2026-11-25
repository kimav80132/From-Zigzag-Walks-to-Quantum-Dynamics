"""The rule: at one spacetime point there is at most one incoming and one outgoing process. Two paths never share a point.
Walk: 1+1 t-reading lattice (shift then coin), path weight = product of cos(a dt) per keep, +sin per flip into +, -sin per flip into - (III.15a (A)).
Two walkers from A1 = (0,+), A2 = (4,+) (same lattice parity), N ticks. Enumerate all path pairs; keep the vertex-disjoint ones (no common (t, x)).
Claims:
1. single-path sums reproduce U_{XA} (the matrix of pauli_test.py).
2. the vertex-disjoint sum, by endpoints (X1 for walker 1, X2 for walker 2), equals det [[U_{X1A1}, U_{X2A1}],[U_{X1A2}, U_{X2A2}]] = D - E when walker 1 ends left of walker 2,
   and is 0 otherwise (two paths cannot swap order without meeting): Lindstrom-Gessel-Viennot, with the walk's signed weights.
3. so the out-sum of the two-walker event is D - E, the event weight (D - E)^2, and at X1 = X2 it is 0: Pauli exclusion, from the rule and nothing else.
4. without the rule (all pairs) the out-sum is D + E by pairing: Bose. The rule is the whole difference."""
import numpy as np, itertools
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
L, N, a, dt = 20, 6, 1.0, 0.3
c_, s_ = np.cos(a*dt), np.sin(a*dt)
def paths(x0, e0):
    """all 2^N paths: list of (points [(t,x)], endpoint (x,e), weight)"""
    out = []
    for flips in itertools.product([0, 1], repeat=N):
        x, e, w, pts = x0, e0, 1.0, [(0, x0)]
        for t, f in enumerate(flips):
            x = (x + e) % L                      # shift: + right, - left
            if f: e = -e; w *= (s_ if e == +1 else -s_)   # flip into +: +sin, into -: -sin
            else: w *= c_
            pts.append((t+1, x))
        out.append((pts, (x, e), w))
    return out
P1, P2 = paths(0, +1), paths(4, +1)
# 1. single-path sums vs the matrix U
T2 = np.array([[0, 1], [-1, 0.]]); R = c_*np.eye(2) + s_*T2
dim = 2*L; idx = lambda x, e: 2*(x % L) + (0 if e == +1 else 1)
S = np.zeros((dim, dim))
for x in range(L): S[idx(x+1, +1), idx(x, +1)] = 1; S[idx(x-1, -1), idx(x, -1)] = 1
U = np.linalg.matrix_power(np.kron(np.eye(L), R) @ S, N)
def pathsum(P):
    acc = {}
    for pts, X, w in P: acc[X] = acc.get(X, 0) + w
    return acc
ps1, ps2 = pathsum(P1), pathsum(P2)
ok("1. the signed path sums equal the matrix elements U_{XA} for both walkers",
   all(abs(ps1.get((x, e), 0) - U[idx(x, e), idx(0, +1)]) < 1e-12 and abs(ps2.get((x, e), 0) - U[idx(x, e), idx(4, +1)]) < 1e-12 for x in range(L) for e in (+1, -1)))
# 2. vertex-disjoint pairs
disjoint, allpairs = {}, {}
for pts1, X1, w1 in P1:
    s1 = set(pts1)
    for pts2, X2, w2 in P2:
        allpairs[(X1, X2)] = allpairs.get((X1, X2), 0) + w1*w2
        if s1.isdisjoint(pts2): disjoint[(X1, X2)] = disjoint.get((X1, X2), 0) + w1*w2
Uel = lambda X, A: U[idx(*X), idx(*A)]
A1, A2 = (0, +1), (4, +1)
unwrap = lambda x: x if x <= L//2 else x - L            # positions near the origin, no wrap in N = 6 ticks
def det(X1, X2): return Uel(X1, A1)*Uel(X2, A2) - Uel(X2, A1)*Uel(X1, A2)
mism = 0; tot_pairs = 0
for X1 in [(x, e) for x in range(L) for e in (+1, -1)]:
    for X2 in [(x, e) for x in range(L) for e in (+1, -1)]:
        d = disjoint.get((X1, X2), 0.0)
        expect = det(X1, X2) if unwrap(X1[0]) < unwrap(X2[0]) else 0.0
        if abs(d) > 1e-14 or abs(expect) > 1e-14: tot_pairs += 1
        if abs(d - expect) > 1e-12: mism += 1
ok(f"2. the vertex-disjoint sum equals the determinant D - E at {tot_pairs - mism} of {tot_pairs} endpoint pairs and differs at {mism}: "
   "LGV needs the tail swap at the first meeting to keep the weight, and the walk's coin factor depends on the sense, so it does not always",
   mism > 0 and tot_pairs - mism > 3*mism)
# the involution test: every first meeting is head-on (a right-mover and a left-mover), and the tail swap changes the coin factor when the two tails leave with different senses
def wt(pts):
    sen = [(+1 if ((pts[i+1][1] - pts[i][1]) % L) == 1 else -1) for i in range(len(pts)-1)]
    w = 1.0
    for i in range(1, len(sen)): w *= c_ if sen[i] == sen[i-1] else (s_ if sen[i] == +1 else -s_)
    return w
head_on = keep = change = 0
for pts1, X1, w1 in P1:
    for pts2, X2, w2 in P2:
        common = [t for t in range(1, N+1) if pts1[t] == pts2[t]]
        if not common: continue
        t0 = common[0]; head_on += (pts1[t0][1] - pts1[t0-1][1]) % L != (pts2[t0][1] - pts2[t0-1][1]) % L
        a = wt(pts1[:t0+1] + pts2[t0+1:])*wt(pts2[:t0+1] + pts1[t0+1:]); b = wt(pts1)*wt(pts2)
        if abs(a - b) < 1e-12: keep += 1
        else: change += 1
ok(f"2'. every first meeting of two paths is head-on ({head_on} = {keep + change}); swapping the tails there keeps the product weight in {keep} cases and changes it in {change}: "
   "the coin at the meeting point is what breaks the determinant", head_on == keep + change and change > 0)
same = [(X, X) for X in [(x, e) for x in range(L) for e in (+1, -1)] if abs(allpairs.get((X, X), 0)) > 1e-12]
ok(f"3. same endpoint: the vertex-disjoint sum is 0 at all {len(same)} sites where the unconstrained sum is not: exclusion, from the rule alone",
   all(abs(disjoint.get((X, X), 0)) < 1e-14 for X, _ in same) and len(same) > 0)
D = lambda X1, X2: Uel(X1, A1)*Uel(X2, A2); E = lambda X1, X2: Uel(X2, A1)*Uel(X1, A2)
ok("4. without the rule, all pairs by endpoints give D, and the identical-walker event sums D + E: Bose",
   all(abs(allpairs.get((X1, X2), 0) - D(X1, X2)) < 1e-12 for X1 in [(x, e) for x in range(L) for e in (+1, -1)] for X2 in [(x, e) for x in range(L) for e in (+1, -1)]))
tot = sum(disjoint[k]**2 for k in disjoint)
totdet = sum(det(X1, X2)**2 for X1 in [(x, e) for x in range(L) for e in (+1, -1)] for X2 in [(x, e) for x in range(L) for e in (+1, -1)])/2
ok(f"5. event weights summed over endpoints: path rule {tot:.4f} (leaks), determinant rule {totdet:.4f} (= 1): the antisymmetric amplitude conserves the weight, the path constraint does not",
   tot < 0.9 and abs(totdet - 1) < 1e-9)
