"""Energy goes up by e Phi; does the coin? 1+1 Dirac with Phi, A_x: E = e Phi +- sqrt(c^2 (p - e A_x)^2 + m^2 c^4).
1. mass shell with the field: (E - e Phi)^2 - c^2 (p - e A_x)^2 = m^2 c^4: the kinetic 4-momentum has invariant length m whatever the field.
2. two-copy walk at rest: the generator is a T (sense turn, the coin) + (e Phi/hbar) J (copy turn); [T, J] = 0 so the two turns run
   independently: exp(t(aT + wJ)) = exp(a t T) exp(w t J). Total phase rate = a + e Phi/hbar = E/hbar; the coin's share stays a.
3. under a boost the copy-turn per step, -e A_mu dx^mu/hbar, is a scalar (Phi dt and A_x dx trade places but the sum is invariant),
   and the coin a per proper time is a scalar: neither can feed the other."""
import sympy as sp
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
p, m, c, e, hb, Phi, Ax, t, w, a = sp.symbols('p m c e hbar Phi A_x t omega a', positive=True)
sz, sx = sp.diag(1, -1), sp.Matrix([[0, 1], [1, 0]])
H = c*(p - e*Ax)*sz + m*c**2*sx + e*Phi*sp.eye(2)
Es = list(H.eigenvals().keys())
ok("1. (E - e Phi)^2 - c^2 (p - e A_x)^2 = m^2 c^4 for both branches", all(sp.simplify((E - e*Phi)**2 - c**2*(p - e*Ax)**2 - m**2*c**4) == 0 for E in Es))
K = lambda X_, Y_: sp.kronecker_product(X_, Y_); T = sp.Matrix([[0, 1], [-1, 0]]); eps = sp.Matrix([[0, -1], [1, 0]])
T4, J4 = K(T, sp.eye(2)), K(sp.eye(2), eps)
G = a*T4 + w*J4
import numpy as np
from scipy.linalg import expm
T4n, J4n = np.array(T4, float), np.array(J4, float); an, wn, tn = 0.7, 0.3, 1.9
Gn = an*T4n + wn*J4n
ok("2. [T, J] = 0 and exp(t G) = exp(a t T) exp(w t J): the sense turn (coin, rate a) and the copy turn (rate e Phi/hbar) run independently",
   T4*J4 == J4*T4 and np.abs(expm(tn*Gn) - expm(an*tn*T4n) @ expm(wn*tn*J4n)).max() < 1e-12)
eta = sp.symbols('eta', real=True); dt, dx = sp.symbols('dt dx', real=True)
# boost of the 4-vector (Phi/c, A_x) and of the step (c dt, dx); the contraction is invariant
L = sp.Matrix([[sp.cosh(eta), sp.sinh(eta)], [sp.sinh(eta), sp.cosh(eta)]]); g = sp.diag(-1, 1)
Avec = sp.Matrix([Phi/c, Ax]); dxv = sp.Matrix([c*dt, dx])
inv = (Avec.T*g*dxv)[0]; inv_b = ((L*Avec).T*g*(L*dxv))[0]
ok("3. the copy-turn per step -e A_mu dx^mu = e(Phi dt - A_x dx) is unchanged by a boost (Phi dt and A_x dx trade places)",
   sp.simplify(sp.expand((inv_b - inv).rewrite(sp.exp))) == 0 and sp.simplify(inv - (-Phi*dt + Ax*dx)) == 0)
