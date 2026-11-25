"""IV.5a: the i is bookkeeping - eight real counters and the quarter-turn J between two columns.   Recreated 9 Nov 2026 (master v34).

  0. IV.4b's real matrices: A_k^2 = 1, T^2 = -1, all six pairs anticommute, A_k symmetric, T antisymmetric.
  1. (3)-(4): J = [[0, -1],[1, 0]] (blocks of 4): J^2 = -1, J X = (-v, u), and J commutes with diag(A_k, A_k) and diag(T, T).
  2. On a 7^3 periodic lattice (central differences), random real u, v: the 8-counter real evolution  d_t X = diag(L, L) X  and the
     complex evolution  i hbar d_t z = [c sum A_k (-i hbar d_k) + i hbar a T] z  of z = u + i v agree; and (5), hbar J d_t X = [...] X,
     is the same statement with J where the complex form has i.
  3a. sum(u^2 + v^2) = sum |z|^2 is conserved.   3b. the slice: v = 0 at t = 0 stays 0; real data stay real in the complex evolution.
  4. T is not a complex structure for the walker: it squares to -1 but anticommutes with the A_k; and no real 4x4 matrix other than a
     multiple of 1 commutes with all of A_x, A_y, A_z, T (the commutant has dimension 1), so no real J lives on one column.
Finite matrices; the evolutions are matrix exponentials of the lattice generator."""
import sys
import numpy as np, scipy.linalg as sl
fails = 0
def ok(name, cond):
    global fails
    print(("PASS  " if cond else "FAIL  ") + name); fails += (not cond)
Ax = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0.]]); Ay = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0.]])
Az = np.diag([1., -1., -1., 1.]); T = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0.]]); I4 = np.eye(4)
mx = lambda M: np.abs(M).max()
As = [Ax, Ay, Az]
ok("0. A_k^2 = 1, T^2 = -1, the six pairs anticommute, A_k symmetric, T antisymmetric",
   all(np.array_equal(A@A, I4) for A in As) and np.array_equal(T@T, -I4)
   and all(mx(P@Q + Q@P) == 0 for i, P in enumerate(As + [T]) for Q in (As + [T])[i+1:])
   and all(np.array_equal(A, A.T) for A in As) and np.array_equal(T, -T.T))
J = np.block([[0*I4, -I4], [I4, 0*I4]]); big = lambda M: np.block([[M, 0*M], [0*M, M]])
X = np.arange(1., 9.)
ok("1. J^2 = -1_8; J (u, v) = (-v, u); J commutes with diag(A_k, A_k) and diag(T, T); J T = T J = [[0, -T],[T, 0]]",
   np.array_equal(J@J, -np.eye(8)) and np.array_equal(J@X, np.concatenate([-X[4:], X[:4]]))
   and all(np.array_equal(J@big(M), big(M)@J) for M in As + [T]) and np.array_equal(J@big(T), np.block([[0*T, -T], [T, 0*T]])))

# ---- 2. the lattice
n, c, a, hb, t = 7, 1.0, 0.8, 1.0, 1.3
S1 = np.roll(np.eye(n), 1, axis=0); d1 = (S1.T - S1)/2                      # central difference on a ring of 7: (d f)(x) = (f(x+1) - f(x-1))/2
E = np.eye(n)
dx, dy, dz = np.kron(np.kron(d1, E), E), np.kron(np.kron(E, d1), E), np.kron(np.kron(E, E), d1)
N3 = n**3
Lgen = -c*(np.kron(dx, Ax) + np.kron(dy, Ay) + np.kron(dz, Az)) + a*np.kron(np.eye(N3), T)     # the generator on one real column, 1372 x 1372
ok("2a. the lattice generator L = -c sum A_k d_k + a T is real and antisymmetric (%d x %d)" % Lgen.shape, mx(Lgen + Lgen.T) < 1e-14 and np.isrealobj(Lgen))
rng = np.random.default_rng(0); u0 = rng.normal(size=4*N3); v0 = rng.normal(size=4*N3)
U = sl.expm(t*Lgen)                                                         # real, orthogonal
u1, v1 = U@u0, U@v0                                                         # the 8-counter evolution is block diagonal: one U per column
# the complex form: H = c sum A_k (-i hbar d_k) + i hbar a T, Hermitian; z(t) = exp(-i H t / hbar) z(0)
H = c*(np.kron(-1j*hb*dx, Ax) + np.kron(-1j*hb*dy, Ay) + np.kron(-1j*hb*dz, Az)) + 1j*hb*a*np.kron(np.eye(N3), T)
w, V = np.linalg.eigh(H); z1 = V@(np.exp(-1j*w*t/hb)*(V.conj().T@(u0 + 1j*v0)))
dev = mx(z1 - (u1 + 1j*v1))
ok("2b. H is Hermitian (%.0e); after t = %.1f the complex evolution of z = u + i v and the two real evolutions agree to %.1e"
   % (mx(H - H.conj().T), t, dev), mx(H - H.conj().T) < 1e-14 and dev < 1e-12)
# (5): hbar J d_t X = [-c hbar sum J A_k d_k + hbar a J T] X  is  J applied to (2); check as matrices on the 8-counter lattice column
G8 = np.kron(np.eye(2), Lgen)                                               # ordering (copy, site, counter)
J8 = np.kron(np.array([[0, -1.], [1, 0]]), np.eye(4*N3))
bigL = lambda M: np.kron(np.eye(2), M)
rhs5 = -c*hb*(J8@bigL(np.kron(dx, Ax)) + J8@bigL(np.kron(dy, Ay)) + J8@bigL(np.kron(dz, Az))) + hb*a*J8@bigL(np.kron(np.eye(N3), T))
ok("2c. (5): hbar J d_t X = [-c hbar sum J A_k d_k + hbar a J T] X is J times (2), every matrix real (%.0e); and J commutes with the lattice generator (%.0e)"
   % (mx(hb*J8@G8 - rhs5), mx(J8@G8 - G8@J8)), mx(hb*J8@G8 - rhs5) < 1e-13 and mx(J8@G8 - G8@J8) < 1e-13 and np.isrealobj(rhs5))
n0 = u0@u0 + v0@v0; n1 = u1@u1 + v1@v1; nz = np.vdot(z1, z1).real
ok("3a. sum(u^2 + v^2) = sum |z|^2 conserved: relative change %.1e (real form), %.1e (complex form)" % (abs(n1/n0 - 1), abs(nz/n0 - 1)),
   abs(n1/n0 - 1) < 1e-13 and abs(nz/n0 - 1) < 1e-13)
zr = V@(np.exp(-1j*w*t/hb)*(V.conj().T@u0.astype(complex)))
ok("3b. the slice: with v = 0 at t = 0 the second column stays 0 (exactly: the evolution is block diagonal), and the complex evolution of real "
   "data stays real: largest imaginary part %.1e" % mx(zr.imag), mx(U@(0*v0)) == 0 and mx(zr.imag) < 1e-12 and mx(zr.real - u1) < 1e-12)
def commutant_dim(mats):
    Am = np.vstack([np.kron(I4, M) - np.kron(M.T, I4) for M in mats]); return 16 - np.linalg.matrix_rank(Am, tol=1e-10)
ok("4. T^2 = -1 but T anticommutes with each A_k, so reading i = T makes the direction matrices antilinear; the real matrices commuting with all "
   "of A_x, A_y, A_z, T are the multiples of 1 only (commutant dimension %d; of A_x, A_y, A_z alone: %d): the i has to be a second column"
   % (commutant_dim(As + [T]), commutant_dim(As)), commutant_dim(As + [T]) == 1 and all(mx(T@A + A@T) == 0 for A in As))
print("ALL PASS" if not fails else "%d FAILED" % fails); sys.exit(1 if fails else 0)
