"""walk_lines_against_ord_mann_eq6.py

Compares the walk's lines with Ord and Mann, "Entwined Paths, Difference Equations and the Dirac Equation"
(quant-ph/0208004), equations (5), (6) and (7), as printed on pages 13-14 of that paper.

The walk (master III.2, III.5a): clock x, two POSITIVE probabilities P+ (time forward), P- (time backward):
    c d_x P+ + d_t P+ = a (P- - P+)
    c d_x P- - d_t P- = a (P+ - P-)
Loss out along the clock: P = exp(-a x / c) Q.

Ord and Mann (6), first pair (left envelope), u = exp(a t) phi:
    d_t u1 =  c d_z u1 - a u2
    d_t u2 = -c d_z u2 + a u1
"""
import sympy as sp
x, t, a, c = sp.symbols('x t a c', positive=True)
Qp, Qm = sp.Function('Qp')(x, t), sp.Function('Qm')(x, t)
ok = True
def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name); ok = ok and bool(cond)

# 1. the walk's lines for P, with P = exp(-a x/c) Q, give lines for Q with no loss term
E = sp.exp(-a*x/c)
Pp, Pm = E*Qp, E*Qm
L1 = sp.simplify((c*sp.diff(Pp, x) + sp.diff(Pp, t) - a*(Pm - Pp))/E)
L2 = sp.simplify((c*sp.diff(Pm, x) - sp.diff(Pm, t) - a*(Pp - Pm))/E)
check("1 loss out: c d_x Q+ + d_t Q+ = a Q-", sp.simplify(L1 - (c*sp.diff(Qp, x) + sp.diff(Qp, t) - a*Qm)) == 0)
check("2 loss out: c d_x Q- - d_t Q- = a Q+", sp.simplify(L2 - (c*sp.diff(Qm, x) - sp.diff(Qm, t) - a*Qp)) == 0)

# 2. solved for d_t (the reading with t as the clock)
dtQp = sp.solve(sp.Eq(L1, 0), sp.diff(Qp, t))[0]
dtQm = sp.solve(sp.Eq(L2, 0), sp.diff(Qm, t))[0]
check("3 read in t: d_t Q+ = -c d_x Q+ + a Q-", sp.simplify(dtQp - (-c*sp.diff(Qp, x) + a*Qm)) == 0)
check("4 read in t: d_t Q- = +c d_x Q- - a Q+", sp.simplify(dtQm - (c*sp.diff(Qm, x) - a*Qp)) == 0)

# 3. Ord and Mann (6), first pair, with u1 = Q-, u2 = Q+ and z = x
u1, u2 = Qm, Qp
om1 = c*sp.diff(u1, x) - a*u2          # d_t u1
om2 = -c*sp.diff(u2, x) + a*u1         # d_t u2
check("5 Ord-Mann (6) line 1 is the walk's line for Q- (u1 = Q-)", sp.simplify(om1 - dtQm) == 0)
check("6 Ord-Mann (6) line 2 is the walk's line for Q+ (u2 = Q+)", sp.simplify(om2 - dtQp) == 0)

# 4. the coupling matrix on (u1, u2) is antisymmetric and squares to minus the identity; Kac's would be symmetric
T = sp.Matrix([[0, -1], [1, 0]])
check("7 coupling on (u1,u2) is antisymmetric, square = -1", T.T == -T and T*T == -sp.eye(2))

# 5. Ord-Mann (6), second pair (right envelope): same lines with a -> -a
u3, u4 = sp.Function('u3')(x, t), sp.Function('u4')(x, t)
second = sp.Matrix([c*sp.diff(u3, x) + a*u4, -c*sp.diff(u4, x) - a*u3])
first_with_minus_a = sp.Matrix([c*sp.diff(u3, x) - (-a)*u4, -c*sp.diff(u4, x) + (-a)*u3])
check("8 Ord-Mann second pair = first pair with a -> -a", sp.simplify(second - first_with_minus_a) == sp.zeros(2, 1))

# 6. where the loss sits: Ord-Mann remove exp(-a t) (loss along t); the walk removes exp(-a x/c) (loss along x)
phi1, phi2 = sp.exp(-a*t)*u1, sp.exp(-a*t)*u2
om5 = sp.simplify(sp.diff(phi1, t) - (c*sp.diff(phi1, x) - a*phi1 - a*phi2))
check("9 Ord-Mann (5) line 1 follows from (6) line 1 with phi = exp(-a t) u", sp.simplify(om5.subs(sp.diff(Qm, t), dtQm)) == 0)
print("ALL PASS" if ok else "SOME FAILED")
