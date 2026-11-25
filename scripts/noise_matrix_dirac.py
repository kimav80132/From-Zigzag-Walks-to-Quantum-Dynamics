"""The 'noise matrix' of the Dirac process: the two-step second moment of the step matrices, N_ij = (M_i M_j + M_j M_i)/2.
In Part II the kick had <xi^mu xi^nu> = eta^{mu nu} (time kick imaginary). Here each kind of move carries a matrix;
two moves in a row give products of matrices, and the symmetric part is the 'variance'.
Readings: (A) 1+1, clock x (III.5a (F)); (B) 1+1, clock t (III.5a (H)); (C) 3+1, clock z (IV.4b (F)); (D) 3+1, clock t (IV.4b (I))."""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
def N(Ms):
    k = len(Ms); I = sp.eye(Ms[0].shape[0]); out = sp.zeros(k)
    for i in range(k):
        for j in range(k):
            S = (Ms[i]*Ms[j] + Ms[j]*Ms[i])/2
            assert S == S[0, 0]*I, "not a multiple of 1"
            out[i, j] = S[0, 0]
    return out
A = sp.diag(1, -1); T = sp.Matrix([[0, 1], [-1, 0]])
nA = N([A, A*T]); ok(f"(A) 1+1 clock x: moves (advance in t, exchange): N = {nA.tolist()}", nA == sp.diag(1, 1))
nB = N([A, T]);   ok(f"(B) 1+1 clock t: moves (advance in x, turn):      N = {nB.tolist()}", nB == sp.diag(1, -1))
Az = sp.diag(1, -1, -1, 1)
Ax = sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]])
Ay = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
T4 = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]])
nC = N([Az, Ax*Az, Ay*Az, Az*T4]); ok(f"(C) 3+1 clock z: moves (advance in t, x-turn, y-turn, exchange): N = {nC.tolist()}", nC == sp.diag(1, -1, -1, 1))
nD = N([Ax, Ay, Az, T4]);           ok(f"(D) 3+1 clock t: moves (along x, y, z, turn):                   N = {nD.tolist()}", nD == sp.diag(1, 1, 1, -1))
# the squares of the operators: Klein-Gordon in each clock
t, x, y, z = sp.symbols('t x y z'); a, c = sp.symbols('a c', positive=True)
Q = sp.Matrix([sp.Function(f'Q{i}')(t, x, y, z) for i in range(1, 5)])
Lz = lambda q: -Az*q.diff(t) + c*Ax*Az*q.diff(x) + c*Ay*Az*q.diff(y) + a*Az*T4*q           # c d_z Q
Lt = lambda q: -c*Ax*q.diff(x) - c*Ay*q.diff(y) - c*Az*q.diff(z) + a*T4*q                  # d_t Q
ok("    clock z: (c d_z)^2 Q = d_t^2 Q - c^2 (d_x^2 + d_y^2) Q + a^2 Q", sp.expand(Lz(Lz(Q)) - (Q.diff(t, 2) - c**2*(Q.diff(x, 2) + Q.diff(y, 2)) + a**2*Q)) == sp.zeros(4, 1))
ok("    clock t: d_t^2 Q = c^2 (d_x^2 + d_y^2 + d_z^2) Q - a^2 Q", sp.expand(Lt(Lt(Q)) - (c**2*(Q.diff(x, 2) + Q.diff(y, 2) + Q.diff(z, 2)) - a**2*Q)) == sp.zeros(4, 1))
