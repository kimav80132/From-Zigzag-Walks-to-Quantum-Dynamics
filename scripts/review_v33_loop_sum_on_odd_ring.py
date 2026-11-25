"""
Review of master v33, VI.4 (D), eq. (8):

    E_0 = sum of the sea's quasi-energies = -(dim/4) pi + (2/pi) sum_{k odd} Tr U^k / k^2 .

Finding: on a ring with an EVEN number of sites every odd power of the tick moves a walker by an odd
number of sites, so Tr U^k = 0 for every odd k and the loop sum is exactly zero.  The number the text
quotes for (8), -62.831853 on 40 sites, is -(80/4) pi with no loop in it: that check does not test the
loop sum.  On a ring with an ODD number of sites the closed loops are the windings round the ring, the
loop sum is not zero, and (8) holds with it.  This script shows both.  Exact finite matrices only.
"""
import numpy as np


def tick(L, th):
    """III.15a (A): shift, then the coin.  Column order (Q+(1), Q-(1), Q+(2), Q-(2), ...)."""
    U = np.zeros((2*L, 2*L)); c, s = np.cos(th), np.sin(th)
    for y in range(L):
        U[2*y,     2*((y - 1) % L)] += c;  U[2*y,     2*((y + 1) % L) + 1] += s
        U[2*y + 1, 2*((y - 1) % L)] += -s; U[2*y + 1, 2*((y + 1) % L) + 1] += c
    return U


ok_all = True
print(" L   dim   E_0 (sea)      -(dim/4)pi     loop sum       -(dim/4)pi + loop sum   first odd k with Tr U^k != 0")
for L, th in ((40, 0.3), (12, 0.3), (7, 0.3), (9, 0.6), (41, 0.3)):
    U = tick(L, th); dim = 2*L
    assert np.allclose(U.T @ U, np.eye(dim))
    eps = -np.angle(np.linalg.eigvals(U))                       # U phi = e^{-i eps} phi
    E0 = eps[eps < 0].sum()
    ks = np.arange(1, 400001, 2)
    tr = np.cos(np.outer(ks, eps)).sum(axis=1)                  # Tr U^k = sum_j e^{-i k eps_j}, real
    loop = (2/np.pi)*(tr/ks**2).sum()
    tail = (2/np.pi)*dim/(2*ks[-1])                             # bound on the omitted tail
    nz = ks[np.abs(tr) > 1e-6]                                  # eigenvalue round-off at k ~ 1e5 is ~1e-9
    first = int(nz[0]) if len(nz) else None
    # direct count for small k: Tr U^k by matrix powers
    P = np.eye(dim); direct = []
    for k in range(1, 2*L + 2):
        P = P @ U
        if k % 2: direct.append(np.trace(P))
    assert np.allclose(direct, tr[:len(direct)], atol=1e-9)
    good = abs(-dim*np.pi/4 + loop - E0) < tail + 1e-9
    if L % 2 == 0:
        good = good and first is None and abs(loop) < 1e-9 and max(abs(np.array(direct))) < 1e-12   # even ring: no odd closed loop
    else:
        good = good and first == L and abs(loop) > 1e-9         # odd ring: the first closed odd loop winds once
    ok_all &= good
    print(f"{L:3d}  {dim:3d}  {E0:12.6f}  {-dim*np.pi/4:12.6f}  {loop:+12.3e}  {-dim*np.pi/4 + loop:12.6f}"
          f"            {first}      {'PASS' if good else 'FAIL'}")

print()
print("even ring: the loop sum is exactly 0 and E_0 = -(dim/4) pi  (the text's -62.831853 is this case)")
print("odd ring : the closed odd loops are the windings round the ring (first at k = L); (8) holds with them")
print("ALL PASS" if ok_all else "FAIL")
raise SystemExit(0 if ok_all else 1)
