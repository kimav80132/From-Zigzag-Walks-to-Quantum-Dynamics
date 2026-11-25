"""Exploration, 12 Nov 2026 (not cited in the master): the original argument (IV.4b: clock z, counters = detectors labelled by time sense) with the
set of detectors COMPLETED: for each time sense, one detector for each orientation of the walker - as the reference, or turned over in the plane
(y,z), (z,x), (x,y).  2 x 4 = 8 detectors, one walker.

Notation.  Orientation = a combination  s + p_yz R_yz + p_zx R_zx + p_xy R_xy  (R_jk = A_j A_k, a half-turn of the plane (j,k) as a turn of space).
Forward part of the walker: F = (F_s, F_yz, F_zx, F_xy); backward part: B = (B_s, B_yz, B_zx, B_xy).  The walker: M = P+ F + P- B, P+- = (1 +- A_z)/2.

  1. The eight detectors are the eight list counters of IV.5a, one for one: F = (u1, v1, v4, -u4), B = (-v3, u3, -u2, -v2); and the frame numbers
     of `one_walker_one_frame.py` are their half-sums and half-differences.
  2. The three rules of IV.4b, on the eight detectors (matrices with entries 0, +-1):
       advance: every forward detector moves with +z, every backward one with -z;
       time flip: forward <-> backward, with the orientation turned a quarter turn in the reference's spin plane (multiplied by R_yz on the right);
       sideways: a move along x takes forward <-> backward with the orientation turned in the plane (z,x) (R_zx on the left); along y, the plane (y,z).
     The lines:   d_t F = -c d_z F - c R_zx d_x B + c R_yz d_y B - a B R_yz ,     d_t B = +c d_z B + c R_zx d_x F - c R_yz d_y F - a F R_yz .
     They are the master's eight lines (IV.5a (2)), mode by mode.
  3. IV.4b's four counters are a closed half of the eight: {forward as-reference, forward xy, backward yz, backward zx}, and on it the rules are
     IV.4b's own.  The other half {forward yz, forward zx, backward as-reference, backward xy} is closed too.
  4. What the four lack: there, "backward" exists only for the orientations turned over in a plane containing z.  A list moves at c along the
     direction its frame points.  With the eight, every orientation has a forward and a backward detector, as the one orientation of 1+1 has; and
     the walker at rest is equal amounts in the forward and the backward detector of the SAME orientation (III.9, remark 5, in three dimensions).
  5. Along one axis: the detectors {as-reference, yz} x {forward, backward} are closed, and split into the two-counter walk of III.5a (H) twice:
     (forward as-reference, backward yz) and (forward yz, backward as-reference).
  6. Direction detectors and twist detectors: along any one axis the eight sort into four + and four - traversal detectors; the three sortings do
     not go together; the four orientations are products of two sideways turns (xy = x-turn times y-turn); the other four are reached only by the
     clock-plane quarter turn without a change of time sense (the field's turn).
  7. The same lines with x as the direction axis (sums and differences of detectors).
  8. The eight sorted by axis (t, x, y, z) and by sum / difference of the two time senses; the six planes as operations (signed permutations);
     the lines in that sorting.
  9. At rest B = F (or B = -F); run along the clock axis: forward x e^{eta/2}, backward x e^{-eta/2}.
 10. The field's turn, the opposite charge and the record on the detectors; the products of the orientations.
 11. Record, flow, alpha and gamma on the detectors in both sortings; the mirror image is forward <-> backward; the clock form of the eight lines;
     along the clock the eight are four two-counter walks.
 12. The lines in the sorting by axis, and the walker in motion: a clock at rest seen by detectors that run.
 13. The same lines with divergence, gradient and curl; the field's turn commutes with the rules, the opposite charge reverses it.
"""
import sys
import numpy as np, scipy.linalg as sl
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all = ok_all and bool(cond); print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()
rng = np.random.default_rng(8)
Ax = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]], float)
Ay = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]], float)
Az = np.diag([1., -1, -1, 1]); T = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]], float); I4 = np.eye(4)
Ryz, Rzx, Rxy = Ay@Az, Az@Ax, Ax@Ay; quat = [I4, Ryz, Rzx, Rxy]; Pp, Pm = (I4 + Az)/2, (I4 - Az)/2
u0 = np.array([1., 0, 0, 0]); v0 = T@u0
qmat = lambda p: sum(c_*E_ for c_, E_ in zip(p, quat))
def lists(F, B):                                                              # the walker with forward part F and backward part B, as the two lists
    M = Pp@qmat(F) + Pm@qmat(B); return np.concatenate([M@u0, M@v0])
D = np.array([lists(*(np.eye(8)[i][:4], np.eye(8)[i][4:])) for i in range(8)]).T        # r = D (F, B)
lab = ['F_s', 'F_yz', 'F_zx', 'F_xy', 'B_s', 'B_yz', 'B_zx', 'B_xy']; cnt = ['u1', 'u2', 'u3', 'u4', 'v1', 'v2', 'v3', 'v4']
one_to_one = all(np.sum(np.abs(D[:, j]) > 0.5) == 1 for j in range(8)) and np.abs(np.abs(D).sum(axis=0) - 1).max() < 1e-12 and abs(abs(np.linalg.det(D)) - 1) < 1e-12
pairs = ", ".join("%s = %s%s" % (lab[j], '-' if D[:, j].sum() < 0 else '', cnt[int(np.argmax(np.abs(D[:, j])))]) for j in range(8))
ok("1. the eight detectors are the eight list counters, one for one: %s. Forward detectors are the counters 1 and 4 of the two lists, backward ones "
   "the counters 2 and 3" % pairs, one_to_one)

# ---- 2. the three rules on the eight detectors
def Lq(X):                                                                    # X q as a 4 x 4 matrix on (s, yz, zx, xy)
    return np.array([np.linalg.lstsq(np.array([E_.ravel() for E_ in quat]).T, (X@E_).ravel(), rcond=None)[0] for E_ in quat]).T
def Rq(X):
    return np.array([np.linalg.lstsq(np.array([E_.ravel() for E_ in quat]).T, (E_@X).ravel(), rcond=None)[0] for E_ in quat]).T
Z = np.zeros((4, 4)); E4 = np.eye(4)
ADV = np.block([[E4, Z], [Z, -E4]])                                           # advance along z
FLIP = np.block([[Z, -Rq(Ryz)], [-Rq(Ryz), Z]])                               # d_t F = -a B R_yz, d_t B = -a F R_yz
SX = np.block([[Z, Lq(Rzx)], [-Lq(Rzx), Z]])                                  # A_x on (F, B):  -c SX d_x
SY = np.block([[Z, -Lq(Ryz)], [Lq(Ryz), Z]])                                  # A_y on (F, B):  -c SY d_y
AA = [np.kron(np.eye(2), A_) for A_ in (Ax, Ay, Az)]; TT = np.kron(np.eye(2), T)
inter = max(np.abs(D@SX - AA[0]@D).max(), np.abs(D@SY - AA[1]@D).max(), np.abs(D@ADV - AA[2]@D).max(), np.abs(D@FLIP - TT@D).max())
ints = all(np.abs(X - np.round(X)).max() < 1e-12 and np.abs(X).max() < 1 + 1e-12 for X in (ADV, FLIP, SX, SY))
dev_modes = max(np.abs(D@sl.expm(1.2*(-1j*(k_[0]*SX + k_[1]*SY + k_[2]*ADV) + 0.8*FLIP)) - sl.expm(1.2*(-1j*sum(x_*A_ for x_, A_ in zip(k_, AA)) + 0.8*TT))@D).max()
                for k_ in rng.normal(size=(30, 3)))
show = lambda X: "[" + "; ".join(" ".join("%2d" % int(round(x)) for x in row) for row in X) + "]"
print("      order (F_s, F_yz, F_zx, F_xy, B_s, B_yz, B_zx, B_xy)")
for nm, X in (("advance", ADV), ("flip   ", FLIP), ("along x", SX), ("along y", SY)): print("      %s = %s" % (nm, show(X)))
ok("2. the three rules on the eight detectors, matrices with entries 0, +-1 (printed): advance (forward with +z, backward with -z); time flip (forward <-> "
   "backward with the orientation multiplied by R_yz on the right); sideways (forward <-> backward with the orientation turned in the plane (z,x) or (y,z)). "
   "The lines d_t F = -c d_z F - c R_zx d_x B + c R_yz d_y B - a B R_yz, d_t B = +c d_z B + c R_zx d_x F - c R_yz d_y F - a F R_yz are the master's eight "
   "lines: matrix by matrix (%.0e), and 30 random modes evolved both ways (%.0e)" % (inter, dev_modes), ints and inter < 1e-12 and dev_modes < 1e-12)

# ---- 3. IV.4b's four counters are a closed half
half1 = [0, 3, 5, 6]; half2 = [1, 2, 4, 7]                                     # {F_s, F_xy, B_yz, B_zx} and {F_yz, F_zx, B_s, B_xy}
closed = max(np.abs(X[np.ix_(half2, half1)]).max() + np.abs(X[np.ix_(half1, half2)]).max() for X in (ADV, FLIP, SX, SY))
is_u = np.abs(D[np.ix_([4, 5, 6, 7], half1)]).max() < 1e-12 and np.abs(D[np.ix_([0, 1, 2, 3], half2)]).max() < 1e-12
Du = D[np.ix_([0, 1, 2, 3], half1)]
same4 = max(np.abs(Du@X[np.ix_(half1, half1)] - Y@Du).max() for X, Y in ((SX, Ax), (SY, Ay), (ADV, Az), (FLIP, T)))
ok("3. IV.4b's four counters are a closed half of the eight: {forward as-reference, forward xy, backward yz, backward zx} is kept by all three rules, and "
   "so is the other half {forward yz, forward zx, backward as-reference, backward xy} (%.0e); the first half is the list u, the second the list v (%s); on "
   "the first half the rules are IV.4b's A_x, A_y, A_z, T exactly (%.0e)" % (closed, is_u, same4), closed < 1e-12 and is_u and same4 < 1e-12)

# ---- 4. what the four lack
def rot3(Q):                                                                  # the turn of space that Q does to the three directions: Q^-1 A_k Q = sum_j Rm[k, j] A_j
    Qi = np.linalg.inv(Q); Am = np.array([A_.ravel() for A_ in (Ax, Ay, Az)]).T
    return np.array([np.linalg.lstsq(Am, (Qi@A_@Q).ravel(), rcond=None)[0] for A_ in (Ax, Ay, Az)])
d4 = 0.0
for _ in range(100):
    p = rng.normal(size=4); Q = qmat(p); q = Q@u0; rho = q@q; j = np.array([q@A_@q for A_ in (Ax, Ay, Az)])
    d4 = max(d4, np.abs(j/rho - rot3(Q)[:, 2]).max())                         # the image of the z direction
turned = [np.round(rot3(E_)[2, 2]) for E_ in quat]                             # does the orientation keep the z direction (+1) or reverse it (-1)?
rest = lists(np.array([1., 0, 0, 0]), np.array([1., 0, 0, 0])); fl = [rest@A_@rest for A_ in AA]
steady = max(np.abs([ (sl.expm(t_*0.8*TT)@rest)@A_@(sl.expm(t_*0.8*TT)@rest) for A_ in AA]).max() for t_ in (0.3, 1.1, 2.6))
ok("4. what the four lack. A single list q = Q u0 moves at c along the direction its orientation Q gives to the z axis (100 random: %.0e). Of the four "
   "orientations, as-reference and xy keep that direction (%+d, %+d) and yz, zx reverse it (%+d, %+d): so in the four-counter walk a detector is forward "
   "exactly when its orientation keeps the direction, and there is NO backward detector for an orientation that keeps it. With the eight, the walker at "
   "rest is equal amounts in the forward and the backward detector of the same orientation, F_s = B_s: flow %s, and it stays zero (%.0e)"
   % (d4, turned[0], turned[3], turned[1], turned[2], np.round(fl, 12).tolist(), steady),
   d4 < 1e-12 and turned == [1, -1, -1, 1] and np.abs(fl).max() < 1e-12 and steady < 1e-12)

# ---- 5. along one axis
sub = [0, 1, 4, 5]; oth = [2, 3, 6, 7]
c5 = max(np.abs(X[np.ix_(oth, sub)]).max() for X in (ADV, FLIP))
g = lambda X, idx: X[np.ix_(idx, idx)]
w1 = [0, 5]; w2 = [1, 4]                                                       # (F_s, B_yz) and (F_yz, B_s)
A2 = np.diag([1., -1]); T2 = np.array([[0, 1.], [-1, 0]])
split = np.abs(FLIP[np.ix_(w1, w2)]).max() + np.abs(FLIP[np.ix_(w2, w1)]).max()
s1 = np.diag([np.sign(D[:, 0].sum()), np.sign(D[:, 5].sum())]); s2 = np.diag([np.sign(D[:, 1].sum()), np.sign(D[:, 4].sum())])
same2 = max(np.abs(s1@g(FLIP, w1)@s1 - T2).max(), np.abs(g(ADV, w1) - A2).max(), np.abs(s2@g(FLIP, w2)@s2 - T2).max(), np.abs(g(ADV, w2) - A2).max())
ok("5. along one axis (nothing depends on x, y): the detectors {as-reference, yz} x {forward, backward} are kept by advance and flip (%.0e) and split into "
   "two walks that do not mix (%.0e): (forward as-reference, backward yz) and (forward yz, backward as-reference); each is III.5a (H)'s two-counter walk, "
   "direction matrix diag(1, -1) and turn matrix [[0, 1], [-1, 0]] (%.0e). The original backward counter is 'backward, and a quarter turn on in the "
   "clock plane'" % (c5, split, same2), c5 < 1e-12 and split < 1e-12 and same2 < 1e-12)

# ---- 6. direction detectors, and twist detectors
E8 = np.eye(8); idx = {n_: i for i, n_ in enumerate(lab)}
ev = lambda X: sorted(np.round(np.linalg.eigvalsh(X)).astype(int).tolist())
each_axis = ev(SX) == [-1]*4 + [1]*4 and ev(SY) == [-1]*4 + [1]*4 and ev(ADV) == [-1]*4 + [1]*4
not_together = np.abs(SX@ADV - ADV@SX).max() > 1 and np.abs(SX@SY - SY@SX).max() > 1 and np.abs(SY@ADV - ADV@SY).max() > 1
Fpart, Bpart = rng.normal(size=4), None
plus_x = np.concatenate([Lq(Rzx)@(b_ := rng.normal(size=4)), b_])              # forward part = the backward part turned in (z,x)
is_plus = np.abs(SX@plus_x - plus_x).max() < 1e-12
def where(vec):                                                               # which single detector (with sign) a vector sits on
    i = int(np.argmax(np.abs(vec))); assert abs(abs(vec[i]) - 1) < 1e-12 and np.abs(np.delete(vec, i)).max() < 1e-12; return ('-' if vec[i] < 0 else '+') + lab[i]
start = E8[idx['F_s']]
chain = (where(SX@start), where(SY@start), where(SY@SX@start), where(SX@SY@start), where(FLIP@start))
JJ = np.block([[np.zeros((4, 4)), -np.eye(4)], [np.eye(4), np.zeros((4, 4))]]); Jd = np.linalg.inv(D)@JJ@D          # the turn between the lists, on the detectors
field = (where(Jd@start), where(Jd@E8[idx['B_yz']]))
xy_is_product = np.abs(SX@SY - np.block([[Lq(Rxy), Z], [Z, Lq(Rxy)]])).max()
ok("6. direction detectors and twist detectors. Sorted along any ONE axis the eight are four detectors of + traversals and four of - traversals (x, y "
   "and z alike: %s), and a +x traversal is the backward part together with the same thing turned in (z,x) in the forward part (%s); but the three "
   "sortings do not go together (%s): direction gives the 2, not the 4. The 4 are products of two sideways turns: from forward-as-reference, a turn "
   "along x gives %s, along y gives %s, x then y gives %s and y then x gives %s - the xy orientation is the x-turn times the y-turn (%.0e). A time flip "
   "from the same start gives %s. The other four detectors are not reached this way: the turn between the two lists (what a field does) takes "
   "F_s to %s and B_yz to %s - the same quarter turn of the clock plane, WITHOUT the change of time sense"
   % (each_axis, is_plus, not_together, *chain[:4], xy_is_product, chain[4], *field),
   each_axis and is_plus and not_together and chain[0][1:] == 'B_zx' and chain[1][1:] == 'B_yz' and chain[2][1:] == 'F_xy' and chain[3][1:] == 'F_xy'
   and chain[2][0] != chain[3][0] and chain[4][1:] == 'B_yz' and xy_is_product < 1e-12 and field[0][1:] == 'F_yz' and field[1][1:] == 'B_s')

# ---- 7. the same lines with another axis as the direction axis
def left8(X):                                                                 # the frame multiplied by X on the left, as a matrix on the eight detectors (F, B)
    cols = []
    for i in range(8):
        f_, b_ = np.eye(8)[i][:4], np.eye(8)[i][4:]; M = X@(Pp@qmat(f_) + Pm@qmat(b_)); r = np.concatenate([M@u0, M@v0]); cols.append(np.linalg.solve(D, r))
    return np.array(cols).T
Ux = left8((I4 + Rzx)/np.sqrt(2))                                             # a quarter turn of space in the plane (z,x): sums and differences of detectors
orth = np.abs(Ux@Ux.T - np.eye(8)).max(); halves = sorted(set(np.round(np.abs(Ux[np.abs(Ux) > 1e-9]), 6).tolist()))
nSX, nSY, nADV, nFL = (Ux@X@Ux.T for X in (SX, SY, ADV, FLIP))
diag_x = np.abs(nSX - np.diag(np.diag(nSX))).max() < 1e-12 and sorted(np.round(np.diag(nSX)).astype(int).tolist()) == [-1]*4 + [1]*4
same_form = min(np.abs(nSX - sg*ADV).max() for sg in (1, -1)) < 1e-12 and min(np.abs(nADV - sg*SX).max() for sg in (1, -1)) < 1e-12 \
    and np.abs(nSY - SY).max() < 1e-12 and np.abs(nFL - FLIP).max() < 1e-12
ok("7. both sortings give the SAME eight lines. The x-direction detectors are sums and differences of the z-direction ones (a quarter turn of space in "
   "the plane (z,x): every entry of the change is 0 or +-%.4f, and it keeps the sum of squares, %.0e). In them the x-advance is diagonal (+1 four "
   "times, -1 four times: %s) and the four matrices are the old ones with x and z exchanged (up to a sign): advance along x <-> advance along z, the "
   "y-turn and the time flip unchanged (%s). One direction axis and the four orientations give the equation; which axis is a free choice"
   % (halves[-1], orth, diag_x, same_form), orth < 1e-12 and diag_x and same_form and abs(halves[-1] - 1/np.sqrt(2)) < 1e-6)

# ---- 8. sums and differences: the eight sorted by axis; the six planes as OPERATIONS; the lines in that sorting
w = Ax@Ay@Az; nbasis = [I4, Ryz, Rzx, Rxy, w, Ax, Ay, Az]; nlab = ['t,sum', 'x,sum', 'y,sum', 'z,sum', 't,dif', 'x,dif', 'y,dif', 'z,dif']
N8 = np.array([np.linalg.solve(D, np.concatenate([E_@u0, E_@v0])) for E_ in nbasis]).T      # (F, B) = N8 n
want = np.zeros((8, 8))                                                       # n in terms of (F, B): sums and differences of one orientation
for row, (i, j, sg) in enumerate(((0, 4, +1), (1, 5, +1), (2, 6, +1), (3, 7, +1), (3, 7, -1), (2, 6, -1), (5, 1, -1), (0, 4, -1))): want[row, i] = 0.5; want[row, j] = 0.5*sg
tab = np.abs(np.linalg.inv(N8) - want).max()
def perm(Xop):                                                                # the operation "multiply the walker on the left by Xop", on the eight n
    Mx = np.array([np.linalg.lstsq(np.array([b_.ravel() for b_ in nbasis]).T, (Xop@E_).ravel(), rcond=None)[0] for E_ in nbasis]).T
    assert np.abs(np.abs(Mx).sum(axis=0) - 1).max() < 1e-12 and np.abs(np.abs(Mx) - np.round(np.abs(Mx))).max() < 1e-12; return np.round(Mx)
ops = {'rot xy': Rxy, 'rot yz': Ryz, 'rot zx': Rzx, 'boost x': Ax, 'boost y': Ay, 'boost z': Az, '0123': w}
keeps = all(np.abs(perm(ops[k_])[:4, 4:]).max() == 0 for k_ in ('rot xy', 'rot yz', 'rot zx'))
swaps = all(np.abs(perm(ops[k_])[:4, :4]).max() == 0 for k_ in ('boost x', 'boost y', 'boost z', '0123'))
ax = lambda i: i % 4
def pairs_of(Mx): return sorted(set(tuple(sorted((ax(j), ax(int(np.argmax(np.abs(Mx[:, j])))))) ) for j in range(8)))
T_, X_, Y_, Z_ = 0, 1, 2, 3
pair_ok = pairs_of(perm(Rxy)) == [(T_, Z_), (X_, Y_)] and pairs_of(perm(Ryz)) == [(T_, X_), (Y_, Z_)] and pairs_of(perm(Rzx)) == [(T_, Y_), (X_, Z_)] \
    and pairs_of(perm(Ax)) == [(T_, X_), (Y_, Z_)] and pairs_of(perm(Ay)) == [(T_, Y_), (X_, Z_)] and pairs_of(perm(Az)) == [(T_, Z_), (X_, Y_)] \
    and pairs_of(perm(w)) == [(T_, T_), (X_, X_), (Y_, Y_), (Z_, Z_)]
# the lines in terms of Sigma = (F + B)/2, Delta = (F - B)/2 (orientation order 0, yz, zx, xy)
Hh = np.block([[E4, E4], [E4, -E4]])/2; Hi = np.linalg.inv(Hh)                 # (Sigma, Delta) = Hh (F, B)
gSD = lambda kv: Hh@(-1j*(kv[0]*SX + kv[1]*SY + kv[2]*ADV) + 0.8*FLIP)@Hi
def gSD_formula(kv):
    Dm = 1j*(kv[2]*E4 - kv[0]*Lq(Rzx) + kv[1]*Lq(Ryz)); Dp = 1j*(kv[2]*E4 + kv[0]*Lq(Rzx) - kv[1]*Lq(Ryz))
    return np.block([[-0.8*Rq(Ryz), -Dm], [-Dp, +0.8*Rq(Ryz)]])
sd_lines = max(np.abs(gSD(k_) - gSD_formula(k_)).max() for k_ in rng.normal(size=(20, 3)))
ok("8. sums and differences. The eight sorted by axis, n = (t, x, y, z sums | t, x, y, z differences), are half-sums and half-differences of the forward "
   "and the backward detector of ONE orientation each (%.0e). The six planes are operations on them, each moving every detector onto one other with a "
   "sign: the three rotations keep sums as sums and differences as differences (%s), the three boosts and 0123 exchange them (%s); a rotation in xy "
   "turns x<->y and t<->z, yz: y<->z and t<->x, zx: z<->x and t<->y; a boost along x exchanges t<->x and y<->z, along y: t<->y and z<->x, along z: t<->z and "
   "x<->y; 0123 keeps the axis (%s). The lines in this sorting: d_t Sigma = -c (d_z - R_zx d_x + R_yz d_y) Delta - a Sigma R_yz, "
   "d_t Delta = -c (d_z + R_zx d_x - R_yz d_y) Sigma + a Delta R_yz (20 random modes, %.0e)" % (tab, keeps, swaps, pair_ok, sd_lines),
   tab < 1e-12 and keeps and swaps and pair_ok and sd_lines < 1e-12)

# ---- 9. at rest, and run
Plead = (np.eye(8) + JJ@TT)/2                                                 # the leading class at rest (the clock turning one way)
Psame = np.block([[E4, E4], [E4, E4]])/2                                      # B = F
c9 = np.abs(np.linalg.inv(D)@Plead@D - Psame).max()
d9 = 0.0
for _ in range(40):
    Fq = rng.normal(size=4)
    for sg in (1, -1):
        r_ = lists(Fq, sg*Fq)
        d9 = max(d9, max(np.abs([(sl.expm(t_*0.8*TT)@r_)@A_@(sl.expm(t_*0.8*TT)@r_) for A_ in AA]).max() for t_ in (0.0, 0.7, 2.2)))
    eta = rng.uniform(0.1, 1.5); Mrun = sl.expm(eta*Az/2)@qmat(Fq); r_run = np.concatenate([Mrun@u0, Mrun@v0]); fb = np.linalg.solve(D, r_run)
    d9 = max(d9, np.abs(fb[:4] - np.exp(eta/2)*Fq).max(), np.abs(fb[4:] - np.exp(-eta/2)*Fq).max(), abs(r_run@AA[2]@r_run/(r_run@r_run) - np.tanh(eta)))
ok("9. at rest and run. The walker whose clock turns one way at rest is exactly B = F: the same amount in the forward and the backward detector of every "
   "orientation (%.0e); B = F and B = -F both have zero flow at all times. Run to rapidity eta along the clock axis, the forward detectors are multiplied "
   "by e^{eta/2} and the backward ones by e^{-eta/2}, so the differences are tanh(eta/2) times the sums, and flow/density = tanh(eta) (40 random, %.0e)"
   % (c9, d9), c9 < 1e-12 and d9 < 1e-12)

# ---- 10. the field's turn, the opposite charge, the record
turnR = np.block([[Rq(Ryz), Z], [Z, Rq(Ryz)]])                                 # orientation multiplied by R_yz on the right, time sense untouched
c10a = np.abs(Jd - turnR).max()
Kd = np.linalg.inv(D)@np.diag([1., 1, 1, 1, -1, -1, -1, -1])@D; added = np.ones(8); added[half2] = -1
c10b = np.abs(Kd - np.diag(added)).max()
mixes = np.abs(Jd[np.ix_(half1, half1)]).max() + np.abs(Jd[np.ix_(half2, half2)]).max()
rec_ok = max(abs(lists(f_[:4], f_[4:])@lists(f_[:4], f_[4:]) - f_@f_) for f_ in rng.normal(size=(30, 8)))
qt = {(a_, b_): where4 for a_ in range(4) for b_ in range(4) for where4 in [None]}
qn = ['1', 'R_yz', 'R_zx', 'R_xy']
def qprod(i, j):
    cq = np.linalg.lstsq(np.array([E_.ravel() for E_ in quat]).T, (quat[i]@quat[j]).ravel(), rcond=None)[0]; k_ = int(np.argmax(np.abs(cq)))
    return ('-' if cq[k_] < 0 else '') + qn[k_]
table = "; ".join("%s %s = %s" % (qn[i], qn[j], qprod(i, j)) for i in range(1, 4) for j in range(1, 4))
ok("10. the field's turn (V.2) on the detectors is the orientation multiplied by R_yz on the right with the time sense untouched (%.0e): it carries every "
   "one of the original four detectors onto one of the added four and back (%.0e). The opposite charge (VI.7's K) is the four added detectors counted "
   "with the other sign (%.0e). The record is the sum of the squares of the eight detectors (%.0e). Products of the orientations: %s"
   % (c10a, mixes, c10b, rec_ok, table), c10a < 1e-12 and mixes < 1e-12 and c10b < 1e-12 and rec_ok < 1e-12)

# ---- 11. record, flow and the two numbers that do not change under the operations, on the detectors; the clock form; four walks along the clock
Ni = np.linalg.inv(N8); d11 = 0.0
for _ in range(60):
    f_ = rng.normal(size=8); F_, B_ = f_[:4], f_[4:]; r_ = lists(F_, B_); n_ = Ni@f_; tp, xp, yp, zp, tm, xm, ym, zm = n_
    rho_ = r_@r_; jx, jy, jz = [r_@A_@r_ for A_ in AA]
    al_ = -(r_[:4]@T@r_[4:]); ga_ = r_[:4]@T@w@r_[4:]
    npl, nmi = np.array([xp, yp, zp]), np.array([xm, ym, zm]); cr = np.cross(npl, nmi)
    d11 = max(d11, abs(rho_ - f_@f_), abs(rho_ - 2*n_@n_),
              abs(jx - 2*F_@(Lq(Rzx)@B_)), abs(jy + 2*F_@(Lq(Ryz)@B_)), abs(jz - (F_@F_ - B_@B_)),
              np.abs(np.array([jx, jy, jz]) - 4*(tp*nmi + tm*npl - cr)).max(),
              abs(al_ - F_@B_), abs(al_ - (n_[:4]@n_[:4] - n_[4:]@n_[4:])), abs(ga_ - 2*(tp*tm - npl@nmi)),
              abs(ga_ - (B_[0]*F_[3] - B_[3]*F_[0] - B_[1]*F_[2] + B_[2]*F_[1])),
              abs(rho_**2 - jx**2 - jy**2 - jz**2 - 4*(al_**2 + ga_**2)))
half_light = max(abs(lists(*np.split(np.where(np.isin(np.arange(8), hf), f_, 0.0), 2))[:4]@T@lists(*np.split(np.where(np.isin(np.arange(8), hf), f_, 0.0), 2))[4:])
                 for hf in (half1, half2) for f_ in rng.normal(size=(10, 8)))
SWAP = np.block([[Z, E4], [E4, Z]])                                            # forward <-> backward, orientation kept
mir2 = 0.0
for _ in range(20):
    f_ = rng.normal(size=8); Mm = Pp@qmat(f_[:4]) + Pm@qmat(f_[4:]); Mt = T@Mm@np.linalg.inv(T); Ms = Pp@qmat(f_[4:]) + Pm@qmat(f_[:4]); mir2 = max(mir2, np.abs(Mt - Ms).max())
# clock form:  c d_z f = [ -ADV d_t  +  a ADV FLIP  -  c ADV SX d_x  -  c ADV SY d_y ] f
EXC, XT, YT = ADV@FLIP, -ADV@SX, -ADV@SY
shapes = np.abs(EXC - EXC.T).max() + np.abs(XT + XT.T).max() + np.abs(YT + YT.T).max() + np.abs(EXC@EXC - np.eye(8)).max() + np.abs(XT@XT + np.eye(8)).max() + np.abs(YT@YT + np.eye(8)).max()
clock4 = max(np.abs(Du@X[np.ix_(half1, half1)] - Y@Du).max() for X, Y in
             ((EXC, np.array([[0, 0, 1, 0], [0, 0, 0, -1], [1, 0, 0, 0], [0, -1, 0, 0.]])), (XT, np.array([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0.]])),
              (YT, np.array([[0, 0, -1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, -1, 0, 0.]]))))
walks = [(0, 5), (1, 4), (2, 7), (3, 6)]                                       # (F_s, B_yz), (F_yz, B_s), (F_zx, B_xy), (F_xy, B_zx)
four_walks = sum(np.abs(FLIP[np.ix_(list(wa), [k_ for k_ in range(8) if k_ not in wa])]).max() for wa in walks)
tw = max(min(np.abs(sg*FLIP[np.ix_(list(wa), list(wa))] - T2).max() for sg in (1, -1)) for wa in walks)
ok("11. on the detectors. Record = F.F + B.B = 2 (sum of the eight sorted numbers squared); flow along z = c (F.F - B.B), along x = 2c F.(R_zx B), along y = "
   "-2c F.(R_yz B); in the sorting by axis, flow_k = 4c [t+ k- + t- k+ - (n+ x n-)_k]; alpha = F.B = (sums squared) - (differences squared); gamma = "
   "2 (t+ t- - x+ x- - y+ y- - z+ z-); record^2 - flow^2/c^2 = 4 (alpha^2 + gamma^2) (60 random, %.0e). On either half of the detectors alone alpha = 0 (%.0e). "
   "The mirror image (runs reversed) is the exchange of the forward and the backward part, orientation kept (%.0e). Clock form: exchange symmetric with "
   "square +1, x-turn and y-turn antisymmetric with square -1, and on IV.4b's four they are IV.4b (A)'s own matrices (%.0e, %.0e). With nothing depending on "
   "x, y the eight are four two-counter walks of III.5a that do not mix: (F_s, B_yz), (F_yz, B_s), (F_zx, B_xy), (F_xy, B_zx) (%.0e, %.0e)"
   % (d11, half_light, mir2, shapes, clock4, four_walks, tw), d11 < 1e-9 and half_light < 1e-12 and mir2 < 1e-12 and shapes < 1e-12 and clock4 < 1e-12
   and four_walks < 1e-12 and tw < 1e-12)

# ---- 12. the walker in motion, in the sorting by axis: the lines are what a clock at rest looks like to detectors that run
Aop = [Ni@X@N8 for X in (SX, SY, ADV)]; Bop = Ni@turnR@N8; Eop = np.diag([1., 1, 1, 1, -1, -1, -1, -1]); Fop = Ni@FLIP@N8
alg = max(np.abs(Fop + Bop@Eop).max(), max(np.abs(Bop@X - X@Bop).max() for X in Aop), max(np.abs(Eop@X + X@Eop).max() for X in Aop), np.abs(Bop@Eop - Eop@Bop).max(),
          max(np.abs(Aop[i]@Aop[j] + Aop[j]@Aop[i] - 2*(i == j)*np.eye(8)).max() for i in range(3) for j in range(3)),
          max(np.abs(Aop[k_] - perm(A_)).max() for k_, A_ in enumerate((Ax, Ay, Az))))
aa = 0.8; d12 = 0.0; d12fd = 0.0; inv_ag = 0.0
d12b = 0.0
for _ in range(30):
    n0 = np.concatenate([rng.normal(size=4), np.zeros(4)]); eta = rng.uniform(0.1, 1.6); nh = rng.normal(size=3); nh /= np.linalg.norm(nh); Nn = sum(c_*X for c_, X in zip(nh, Aop))
    def nfield(t_, x_):
        tau = t_*np.cosh(eta) - (nh@x_)*np.sinh(eta)
        return sl.expm(eta/2*Nn)@(np.cos(aa*tau)*np.eye(8) - np.sin(aa*tau)*Bop)@n0
    t_, x_ = rng.normal(), rng.normal(size=3); n_ = nfield(t_, x_)
    dt_ = -aa*np.cosh(eta)*Bop@n_; dk_ = [aa*np.sinh(eta)*nh[k_]*Bop@n_ for k_ in range(3)]
    m0 = np.concatenate([np.zeros(4), rng.normal(size=4)]); tau_ = t_*np.cosh(eta) - (nh@x_)*np.sinh(eta)           # differences only at rest: the clock turns the other way
    m_ = sl.expm(eta/2*Nn)@(np.cos(aa*tau_)*np.eye(8) + np.sin(aa*tau_)*Bop)@m0
    d12b = max(d12b, np.abs(aa*np.cosh(eta)*Bop@m_ - sum(Aop[k_]@(aa*np.sinh(eta)*nh[k_]*Bop@m_) for k_ in range(3)) + aa*Bop@Eop@m_).max(),
               abs((m_[:4]@m_[:4] - m_[4:]@m_[4:]) + m0@m0), abs((n_[:4]@n_[:4] - n_[4:]@n_[4:]) - n0@n0))
    d12 = max(d12, np.abs(dt_ + sum(Aop[k_]@dk_[k_] for k_ in range(3)) + aa*Bop@Eop@n_).max(), np.abs(sl.expm(-eta*Nn) - (np.cosh(eta)*np.eye(8) - np.sinh(eta)*Nn)).max(),
              np.abs(n_[4:] - np.tanh(eta/2)*(Nn@np.concatenate([n_[:4], np.zeros(4)]))[4:]).max())
    h_ = 1e-5; fd_t = (nfield(t_ + h_, x_) - nfield(t_ - h_, x_))/(2*h_); fd_k = [(nfield(t_, x_ + h_*np.eye(3)[k_]) - nfield(t_, x_ - h_*np.eye(3)[k_]))/(2*h_) for k_ in range(3)]
    d12fd = max(d12fd, np.abs(fd_t - dt_).max(), max(np.abs(fd_k[k_] - dk_[k_]).max() for k_ in range(3)))
    f_ = N8@n_; r_ = D@f_; vel = np.array([r_@A_@r_ for A_ in AA])/(r_@r_); d12 = max(d12, np.abs(vel - np.tanh(eta)*nh).max())
    g_ = rng.normal(size=8); th = rng.normal(); op = sl.expm(eta/2*Nn)@sl.expm(th/2*Aop[0]@Aop[1])@sl.expm(-0.4*Aop[1]@Aop[2]); g2 = op@g_
    alga = lambda v_: (v_[:4]@v_[:4] - v_[4:]@v_[4:], 2*(v_[0]*v_[4] - v_[1:4]@v_[5:8]))
    inv_ag = max(inv_ag, abs(alga(g_)[0] - alga(g2)[0]), abs(alga(g_)[1] - alga(g2)[1]))
ok("12. the walker in motion (units c = 1). In the sorting by axis the lines are d_t n = -(AA_x d_x + AA_y d_y + AA_z d_z) n - a BB EE n, where AA_k is "
   "the boost operation along k of check 8, BB the quarter turn in the walker's own plane and EE = diag(1,1,1,1,-1,-1,-1,-1); BB commutes with the "
   "three AA_k, EE anticommutes with them, the AA_k square to 1 and anticommute (%.0e). A walker at rest with sums only, its clock turning, seen by detectors "
   "that run with rapidity eta along a direction: n = exp(eta/2 n.AA) (cos a tau - sin a tau BB) n0, tau = t cosh eta - n.x sinh eta. Then d_t n = "
   "-a cosh eta BB n, d_k n = a sinh eta n_k BB n, and the lines hold; the differences are tanh(eta/2) times the sums moved by the run; flow/record = "
   "tanh(eta) along the direction (30 random, %.0e; derivatives against finite differences %.0e). alpha and gamma of check 11 are unchanged by runs and "
   "rotations of any size (%.0e). The walker with differences only at rest, its clock turning the other way, satisfies the same lines when run; "
   "alpha is + (sums squared) for the first and - (differences squared) for the second at every rapidity (%.0e)" % (alg, d12, d12fd, inv_ag, d12b),
   alg < 1e-12 and d12 < 1e-9 and d12fd < 1e-6 and inv_ag < 1e-9 and d12b < 1e-9)

# ---- 13. the lines in the sorting by axis, written with divergence, gradient and curl; and the field's turn commutes with the three rules
d13 = 0.0
for _ in range(40):
    kv = rng.normal(size=3); n_ = rng.normal(size=8); tp, npl, tm, nmi = n_[0], n_[1:4], n_[4], n_[5:8]
    lhs = -sum(kv[k_]*Aop[k_] for k_ in range(3))@n_                             # the term -c AA_k d_k n, with d_k -> k_k, c = 1
    rhs = np.concatenate([[-kv@nmi], -kv*tm - np.cross(kv, nmi), [-kv@npl], -kv*tp + np.cross(kv, npl)])
    flip_ = np.array([npl[0], -tp, npl[2], -npl[1], nmi[0], -tm, -nmi[2], nmi[1]])
    d13 = max(d13, np.abs(lhs - rhs).max(), np.abs(Fop@n_ - flip_).max())
comm = max(np.abs(turnR@X - X@turnR).max() for X in (ADV, SX, SY, FLIP)); Cd = np.diag(added)
ch = max(np.abs(Cd@X - X@Cd).max() for X in (ADV, SX, SY, FLIP)) + np.abs(Cd@turnR + turnR@Cd).max()
Cn = np.round(Ni@Cd@N8); c_swaps = np.abs(Cn[:4, :4]).max() + np.abs(Cn[4:, 4:]).max()      # the opposite charge, in the sorted order: sums <-> differences
ok("13. in the sorting by axis, with n+ = (x+, y+, z+) and n- = (x-, y-, z-), the lines read  d_t t+ = -c div n- + a x+,  d_t n+ = -c grad t- - c curl n- "
   "+ a (-t+, z+, -y+),  d_t t- = -c div n+ + a x-,  d_t n- = -c grad t+ + c curl n+ + a (-t-, -z-, y-)  (40 random, %.0e). The quarter turn in the "
   "walker's own plane with the time sense untouched commutes with advance, both sideways rules and the time flip (%.0e); counting the four added "
   "detectors with the other sign commutes with the rules and reverses that turn (%.0e); in the sorted order it carries every sum onto a difference "
   "and every difference onto a sum (%.0e): the walker whose clock turns the other way is the walker of the opposite charge"
   % (d13, comm, ch, c_swaps), d13 < 1e-12 and comm < 1e-12 and ch < 1e-12 and c_swaps < 1e-12)
print("      opposite charge, sorted order (s_t, s_x, s_y, s_z | d_t, d_x, d_y, d_z): " + show(Cn))

# ---- 14. the eight as products of the three direction detectors (the frame reading is this change of variable)
prod_names = ['1', 'y.z', 'z.x', 'x.y', 'x.y.z', 'x', 'y', 'z']
reach = [np.eye(8), Aop[1]@Aop[2], Aop[2]@Aop[0], Aop[0]@Aop[1], Aop[0]@Aop[1]@Aop[2], Aop[0], Aop[1], Aop[2]]
e_st = np.eye(8)[0]; lands = [int(np.argmax(np.abs(R_@e_st))) for R_ in reach]; sgn = [float((R_@e_st)[i]) for R_, i in zip(reach, lands)]
each_once = lands == list(range(8)) and all(abs(x - 1) < 1e-12 for x in sgn)
d14 = 0.0
for _ in range(40):
    e1, e2, e3 = rng.uniform(-1.2, 1.2, size=3); t1, t2, t3 = np.tanh(e1/2), np.tanh(e2/2), np.tanh(e3/2)
    n_ = sl.expm(e1/2*Aop[0])@sl.expm(e2/2*Aop[1])@sl.expm(e3/2*Aop[2])@e_st           # the reference, run along z, then y, then x
    st, sx, sy, sz, dt, dx, dy, dz = n_/n_[0]
    d14 = max(d14, abs(dx - t1), abs(dy - t2), abs(dz - t3), abs(sz - t1*t2), abs(sx - t2*t3), abs(sy + t1*t3), abs(dt - t1*t2*t3))
    f_ = N8@n_; Mm = Pp@qmat(f_[:4]) + Pm@qmat(f_[4:]); d14 = max(d14, np.abs(Mm - sum(c_*E_ for c_, E_ in zip(n_, nbasis))).max())
Mof = lambda v_: sum(c_*E_ for c_, E_ in zip(v_, nbasis)); Ti = np.linalg.inv(T); d14f = 0.0
for _ in range(40):
    v_ = rng.normal(size=8); Mv = Mof(v_); f_ = N8@v_; r_ = D@f_
    rev = Mof(v_*np.array([1, -1, -1, -1, 1, -1, -1, -1])); al_ = v_[:4]@v_[:4] - v_[4:]@v_[4:]; ga_ = 2*(v_[0]*v_[4] - v_[1:4]@v_[5:8])
    d14f = max(d14f, max(np.abs(Mof(Aop[k_]@v_) - A_@Mv).max() for k_, A_ in enumerate((Ax, Ay, Az))), np.abs(Mof(Bop@v_) - Mv@Ryz).max(),
               np.abs(Mof(Eop@v_) - T@Mv@Ti).max(), np.abs(Mof(Fop@v_) + T@Mv@Ti@Ryz).max(), np.abs(Mof(Cn@v_) - Mv@Az).max(),
               abs(r_@r_ - 0.5*np.trace(Mv.T@Mv)), np.abs(Mv@rev - al_*I4 - ga_*w).max())
ok("14. products of the three direction detectors. From the detector s_t alone, the eight products 1, y.z, z.x, x.y, x.y.z, x, y, z of the three runs "
   "land on the eight sorted detectors s_t, s_x, s_y, s_z, d_t, d_x, d_y, d_z, each once, with sign + (%s). The walker written on these eight products "
   "is draft 1's frame, M = s_t + s_x R_yz + s_y R_zx + s_z R_xy + d_t w + d_x A_x + d_y A_y + d_z A_z. And for the reference run along z, then y, then "
   "x, the amounts are literally products: d_x/s_t = t1, d_y/s_t = t2, d_z/s_t = t3, s_z/s_t = t1 t2, s_x/s_t = t2 t3, s_y/s_t = -t1 t3, d_t/s_t = "
   "t1 t2 t3, with t_k = tanh(eta_k/2) (40 random, %.0e). The rules on the frame: a run is A_k M; the walker's own turn is M R_yz; the exchange of "
   "forward and backward is the mirror image T M T^-1 (runs and w reversed); the time flip is -(T M T^-1) R_yz; the opposite charge is M A_z; the "
   "record is (1/2) Tr M^T M; M times M with its six plane turns reversed is alpha + gamma w (40 random, %.0e)" % (each_once, d14, d14f),
   each_once and d14 < 1e-12 and d14f < 1e-12)

# ---- 15. the effective flip rate (a potential against a change of flip rate), and the null form of the lines
phi, Vv = 0.37, 0.21; Jn = Bop; Xn = Eop; Ps, Pd = (np.eye(8) + Xn)/2, (np.eye(8) - Xn)/2
G_em = -aa*Jn@Xn - phi*Jn; G_sc = -(aa + Vv)*Jn@Xn                                     # at rest: mass + scalar potential; mass with the coin changed
e15 = max(np.abs(G_em@Ps + (aa + phi)*Jn@Ps).max(), np.abs(G_em@Pd - (aa - phi)*Jn@Pd).max(),
          np.abs(G_sc@Ps + (aa + Vv)*Jn@Ps).max(), np.abs(G_sc@Pd - (aa + Vv)*Jn@Pd).max(),
          np.abs(phi*Jn - phi*Jn@Xn - 2*phi*Jn@Pd).max(), np.abs(Fop + Jn@Xn).max())
n15 = 0.0
for kv in rng.normal(size=(20, 3)):                                                    # (d_t + c d_z) F and (d_t - c d_z) B, on a mode exp(i k.x), c = 1
    Gk = -1j*(kv[0]*SX + kv[1]*SY + kv[2]*ADV) + aa*FLIP
    rhs_F = -(1j*kv[0]*Lq(Rzx)) + 1j*kv[1]*Lq(Ryz)                                        # -c L_zx d_x B + c L_yz d_y B
    n15 = max(n15, np.abs((Gk + 1j*kv[2]*np.eye(8))[:4, 4:] - (rhs_F - aa*Rq(Ryz))).max(), np.abs((Gk + 1j*kv[2]*np.eye(8))[:4, :4]).max(),
              np.abs((Gk - 1j*kv[2]*np.eye(8))[4:, :4] - (-rhs_F - aa*Rq(Ryz))).max(), np.abs((Gk - 1j*kv[2]*np.eye(8))[4:, 4:]).max())
ok("15. the effective flip rate. At rest with a scalar potential, d_t n = (-a JJ XX - (e Phi/hbar) JJ) n: on the sums it is -(a + e Phi/hbar) JJ, on the "
   "differences +(a - e Phi/hbar) JJ; a change of the coin, a -> a + V/hbar in the mass term alone, gives -(a + V/hbar) JJ on the sums and "
   "+(a + V/hbar) JJ on the differences; and (e Phi/hbar) JJ = (e Phi/hbar) JJ XX + 2 (e Phi/hbar) JJ (projection on the differences) (%.0e). The "
   "null form: (d_t + c d_z) F = -c L_zx d_x B + c L_yz d_y B - a J B and (d_t - c d_z) B = +c L_zx d_x F - c L_yz d_y F - a J F, with no F on the "
   "right of the first and no B on the right of the second (20 random modes, %.0e)" % (e15, n15), e15 < 1e-12 and n15 < 1e-12)

print("ALL PASS" if ok_all else "SOME FAILED"); raise SystemExit(0 if ok_all else 1)
