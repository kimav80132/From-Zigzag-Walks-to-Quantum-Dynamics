"""Can the real two-counter walk feel an electric potential?  A potential V adds (V/hbar) J Q to d_t Q, with J a real matrix
that turns like i (J^2 = -1) and commutes with everything in the lines (A_z and the turn T), so that it acts as a phase.
1. two counters: no real J with J^2 = -1 commuting with A_z and T
2. two copies of the walk (four counters, the copies side by side): J = 1 (x) eps, eps = [[0,1],[-1,0]] between the copies, works
"""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
Az = sp.diag(1, -1); T = sp.Matrix([[0, 1], [-1, 0]])
j = sp.symbols('j1:5', real=True); J = sp.Matrix([[j[0], j[1]], [j[2], j[3]]])
sol = sp.solve(list(J*Az - Az*J) + list(J*T - T*J) + list(J*J + sp.eye(2)), j, dict=True)
ok(f"1. two counters: real J with J^2 = -1, [J, A_z] = 0, [J, T] = 0 -> solutions {sol}", sol == [])
eps = sp.Matrix([[0, 1], [-1, 0]]); K = lambda A, B: sp.kronecker_product(A, B)
Az4, T4, J4 = K(Az, sp.eye(2)), K(T, sp.eye(2)), K(sp.eye(2), eps)
ok("2. two copies: J = 1 (x) eps squares to -1 and commutes with A_z and T", J4*J4 == -sp.eye(4) and J4*Az4 == Az4*J4 and J4*T4 == T4*J4)
