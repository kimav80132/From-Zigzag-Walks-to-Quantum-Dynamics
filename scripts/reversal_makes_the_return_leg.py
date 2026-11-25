"""Does a process that reverses its time sense give the event its return leg? Two walks, run backward.
S (III.1, read in t): tick = coin R = cos(a dt) 1 + sin(a dt) T after the shift. T is the quarter-turn, so the tick is a rotation.
T (telegraph, real rate): tick = (1 - a dt) 1 + a dt sigma_x after the shift: a stochastic matrix, not a rotation.
1. S: the inverse tick is the transpose, S^{-1} R^T, and with the senses swapped (P = sigma_x on the senses) it is the forward walk again,
   P U^{-N} P = (S R)^N: the same shift, the same coin, in the other order. A walker running backward along its evolution is a walker of the same
   kind with every segment's time sense reversed. So the return leg EXISTS, as a walk, and its sum is the out-sum: (U^{-1})_{AX} = U_{XA}.
2. T: the inverse tick is (1 + a dt) 1 - a dt sigma_x + O(dt^2): a coin with a NEGATIVE flip weight. Run backward the telegraph is not a walk, and
   (U_T^{-1})_{AX} != (U_T)_{XA}: the return sum is not the out sum, the event weight out x return is not a square, and some of them are negative.
3. the chain: the coin in t (forced by the sign of Klein-Gordon, III.6, i.e. by eta^{00} = -1) -> the exchange is a quarter-turn (III.8) -> the tick is a
   rotation -> the inverse is the transpose -> the return leg is the same walk with the senses swapped -> the event is out x return = Q^2 (III.15a)."""
import numpy as np
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
L, N, a, dt = 24, 8, 1.0, 0.25
idx = lambda x, e: 2*(x % L) + (0 if e == +1 else 1); dim = 2*L
S = np.zeros((dim, dim))
for x in range(L): S[idx(x+1, +1), idx(x, +1)] = 1; S[idx(x-1, -1), idx(x, -1)] = 1
P = np.kron(np.eye(L), np.array([[0, 1], [1, 0.]]))                       # swap the senses
T2 = np.array([[0, 1], [-1, 0.]]); sx = np.array([[0, 1], [1, 0.]])
R_S = np.kron(np.eye(L), np.cos(a*dt)*np.eye(2) + np.sin(a*dt)*T2)
R_T = np.kron(np.eye(L), (1 - a*dt)*np.eye(2) + a*dt*sx)
U_S, U_T = np.linalg.matrix_power(R_S @ S, N), np.linalg.matrix_power(R_T @ S, N)
A = idx(0, +1)
ok("1a. S: P U^{-N} P = (S R)^N exactly - the backward walk is the forward walk with the senses swapped (shift and coin in the other order)",
   np.allclose(P @ np.linalg.inv(U_S) @ P, np.linalg.matrix_power(S @ R_S, N)))
ok("1b. S: (U^{-1})_{AX} = U_{XA} for every X: the return sum is the out sum; out x return = U_{XA}^2 >= 0, summing to 1",
   np.allclose(np.linalg.inv(U_S)[A, :], U_S[:, A]) and abs(np.sum(U_S[:, A]*np.linalg.inv(U_S)[A, :]) - 1) < 1e-12)
Rinv_T = np.linalg.inv(R_T)
ok(f"2a. T: the inverse coin has a negative flip weight ({Rinv_T[0, 1]:.4f}, against +{R_T[0, 1]:.4f} forward): run backward the telegraph is not a walk",
   Rinv_T[0, 1] < 0)
w = U_T[:, A]*np.linalg.inv(U_T)[A, :]
ok(f"2b. T: (U^{-1})_{{AX}} != U_{{XA}}; the products out x return are not squares - {np.sum(w < -1e-12)} of {np.sum(np.abs(w) > 1e-12)} are negative (they still sum to 1, as any inverse does)",
   not np.allclose(np.linalg.inv(U_T)[A, :], U_T[:, A]) and np.sum(w < -1e-12) > 0 and abs(np.sum(w) - 1) < 1e-12)
ok("3. the one difference: R_S^T R_S = 1 (a rotation), R_T^T R_T != 1 (a stochastic matrix); the quarter-turn is what III.8 derived from conserving Q+^2 + Q-^2",
   np.allclose(R_S.T @ R_S, np.eye(dim)) and not np.allclose(R_T.T @ R_T, np.eye(dim)))
