"""
VI.1 (D): why "one process per state" has to be stated for every column, not for sites and senses only.

Two-walker column C_{IJ}, free step C -> U C U^T (one orthogonal U per walker). A rule stated on the diagonal of one
basis, C_{II} = 0 for the site-sense columns e_I, is not kept by one tick:
    (U C U^T)_{II} = sum_{J,K} U_{IJ} C_{JK} U_{IK}
involves the off-diagonal entries C_{JK}, J != K, because the coin sends e_I to a combination of several e_J.
The rule stated for every column, v^T C v = 0 for all v  <=>  C^T = -C, is kept exactly.

Checks (1+1 lattice of III.15a: shift then coin, 16 sites, 20 ticks):
  1. symmetric C with C_{II} = 0 for every site-sense I (random otherwise): after one tick max |C_{II}| is O(1)
  2. antisymmetric C: stays antisymmetric and C_{II} = 0 after every tick, to roundoff
  3. the diagonal rule in the site-sense basis and the same rule in the tick-rotated basis U e_I differ for the same C:
     the diagonal rule is a statement about which basis, the antisymmetric rule about none
"""
import numpy as np, sys
def ok(msg, cond): print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()
L = 16; th = 0.3; dim = 2*L; U = np.zeros((dim, dim))
for x in range(L):
    for s, y in ((0, (x + 1) % L), (1, (x - 1) % L)):
        R = np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]])
        U[2*y, 2*x + s] = R[0, s]; U[2*y + 1, 2*x + s] = R[1, s]
assert np.allclose(U.T @ U, np.eye(dim))
rng = np.random.default_rng(1)
A = rng.normal(size=(dim, dim)); Cs = A + A.T; np.fill_diagonal(Cs, 0.0); Cs /= np.abs(Cs).max()
Ca = A - A.T; Ca /= np.abs(Ca).max()
C1 = U @ Cs @ U.T
ok("1. symmetric column with no weight on any same site-and-sense pair: after one tick the same-state entries are max |C_II| = %.3f (of entries of size 1): the diagonal rule is not kept" % np.abs(np.diag(C1)).max(), np.abs(np.diag(C1)).max() > 0.05)
C = Ca.copy(); worst = 0.0
for t in range(20):
    C = U @ C @ U.T; worst = max(worst, np.abs(C + C.T).max(), np.abs(np.diag(C)).max())
ok("2. antisymmetric column (v^T C v = 0 for every v): after 20 ticks max |C + C^T|, max |C_II| = %.1e: kept exactly" % worst, worst < 1e-13)
# 3: the diagonal rule read in the rotated basis
diag_rot = np.array([U[:, i] @ Cs @ U[:, i] for i in range(dim)])
ok("3. the same symmetric column read in the basis the tick makes, (U e_I)^T C (U e_I): max %.3f, so 'no two walkers on one site and sense' holds in one basis and fails in the next; only C^T = -C says it in every basis" % np.abs(diag_rot).max(), np.abs(diag_rot).max() > 0.05)
