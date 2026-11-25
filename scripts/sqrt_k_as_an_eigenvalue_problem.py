"""The square root of k as an eigenvalue problem. A = diag(1,-1), T = [[0,1],[-1,0]], AT = -TA.
1. The boost on the null lines is the matrix e^{eta A} = diag(k, 1/k): eigenvectors = the two senses, eigenvalues k and 1/k (the FULL factors).
   This is what the derivative operator picks up: d_t + A d_x = e^{-eta A} (d_t' + A d_x').
2. The column must be rescaled by a matrix S that commutes with A (so the streaming keeps its form): S = e^{s A}, eigenvalues e^{s}, e^{-s} on the two senses.
   The mass term then becomes a S^{-1} e^{eta A} T S. Because T swaps the senses (AT = -TA), S^{-1} e^{eta A} T S = e^{(eta - 2s) A} T:
   an off-diagonal matrix sees the RATIO of the eigenvalues of S, squared. Form invariance is eta - 2s = 0, i.e.
       S^2 = e^{eta A}:  the column's boost is the square root of the null lines' boost. Eigenvalues sqrt(k), 1/sqrt(k).
3. III.8a's "once or twice": a vector is represented by a matrix hit from both sides, S V S^{-1}; with V = [[0, v],[u, 0]] (gamma.x in III.13's basis),
   S V S^{-1} = [[0, e^{2s} v],[e^{-2s} u, 0]]: the vector's null parts get e^{+-2s} = k^{+-1} when the column gets e^{+-s} = k^{+-1/2}.
   Same S, applied once to a column and twice to a vector."""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
eta, s, u, v = sp.symbols('eta s u v', real=True)
A = sp.diag(1, -1); T = sp.Matrix([[0, 1], [-1, 0]])
expA = lambda z: sp.diag(sp.exp(z), sp.exp(-z))                 # e^{z A}
k = sp.exp(eta)
ok("1. e^{eta A} = diag(k, 1/k): eigenvectors (1,0), (0,1) = the two senses; eigenvalues k and 1/k",
   expA(eta) == sp.diag(k, 1/k) and set(expA(eta).eigenvals().keys()) == {k, 1/k})
S = expA(s)
conj = sp.simplify(S.inv()*expA(eta)*T*S)
ok("2a. S^{-1} e^{eta A} T S = e^{(eta - 2s) A} T: the off-diagonal T sees the ratio of S's eigenvalues", sp.simplify(conj - expA(eta - 2*s)*T) == sp.zeros(2, 2))
ok("2b. form invariance <=> eta - 2s = 0 <=> S^2 = e^{eta A}: eigenvalues of S are sqrt(k) and 1/sqrt(k)",
   sp.simplify((S*S - expA(eta)).subs(s, eta/2)) == sp.zeros(2, 2) and set(S.subs(s, eta/2).eigenvals().keys()) == {sp.sqrt(k), 1/sp.sqrt(k)})
V = sp.Matrix([[0, v], [u, 0]])
ok("3. S V S^{-1} = [[0, e^{2s} v],[e^{-2s} u, 0]]: a vector, hit from both sides, gets k and 1/k while the column, hit once, gets sqrt(k) and 1/sqrt(k)",
   sp.simplify(S*V*S.inv() - sp.Matrix([[0, sp.exp(2*s)*v], [sp.exp(-2*s)*u, 0]])) == sp.zeros(2, 2))
