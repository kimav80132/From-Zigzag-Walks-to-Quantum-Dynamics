"""2+1, clock z (the path parameter). Counters: + (forward in t), - (backward in t). Loss taken out.
One step dz is made of three kinds of move, each with a weight and a displacement in (t, x):
  time move : counter + shifts t by +dz/c, counter - by -dz/c        (carried by the counter, weight 1)
  switch    : from the other counter, half from x-ell, half from x+ell, weights +-1/2   (ell = dz)
  flip      : from the other counter, weight a dz/c, no displacement
For each move: total weight (zeroth moment) and signed displacement (first moment).
Then the two-step moments: for every ordered pair of consecutive moves, product of their signed displacements.
Claim: summed over consecutive pairs, (t,t) gives +, (x,x) gives -, mixed pairs cancel between the two orders,
       two flips give +a^2: i.e. c^2 dz^2 = dt^2 - c^2 dx^2 + a^2, Klein-Gordon read along z, with the Lorentzian minus
       produced by the moves themselves."""
import sympy as sp
c, a, ell = sp.symbols('c a ell', positive=True)
dz = ell
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
# moves: (from, to, weight, dt, dx)
moves = []
for k, s in (('+', +1), ('-', -1)):
    moves.append(dict(kind='time', frm=k, to=k, w=1, dt=s*dz/c, dx=0))
# switches into '+' from '-': +1/2 from x-ell (displacement +ell), -1/2 from x+ell (displacement -ell)
# switches into '-' from '+': +1/2 from x+ell (displacement -ell), -1/2 from x-ell (displacement +ell)
for to, frm, sign in (('+', '-', +1), ('-', '+', -1)):
    moves.append(dict(kind='switch', frm=frm, to=to, w=sp.Rational(1, 2)*sign, dt=0, dx=+ell))
    moves.append(dict(kind='switch', frm=frm, to=to, w=-sp.Rational(1, 2)*sign, dt=0, dx=-ell))
for to, frm in (('+', '-'), ('-', '+')):
    moves.append(dict(kind='flip', frm=frm, to=to, w=a*dz/c * (1 if to == '+' else 1), dt=0, dx=0))
# first moments of a switch
def switch_moments(to):
    ms = [m for m in moves if m['kind'] == 'switch' and m['to'] == to]
    return sum(m['w'] for m in ms), sum(m['w']*m['dx'] for m in ms)
w0p, d1p = switch_moments('+'); w0m, d1m = switch_moments('-')
ok(f"a switch carries zero total weight ({w0p}, {w0m}) but a signed displacement: +ell entering '+', -ell entering '-' ({d1p}, {d1m})",
   w0p == 0 and w0m == 0 and sp.simplify(d1p - ell) == 0 and sp.simplify(d1m + ell) == 0)
# collapse each move kind to its first moment as an operator on (Q+, Q-): coefficient of d_t, d_x, and constant
def first_moment(kind):
    M = {}
    for m in moves:
        if m['kind'] != kind: continue
        key = (m['to'], m['frm'])
        M.setdefault(key, [0, 0, 0])
        M[key][0] += m['w']*m['dt']; M[key][1] += m['w']*m['dx']; M[key][2] += m['w'] if kind == 'flip' else 0
    return M
T, X, F = first_moment('time'), first_moment('switch'), first_moment('flip')
idx = {'+': 0, '-': 1}
def mat(M, comp):
    A = sp.zeros(2)
    for (to, frm), v in M.items(): A[idx[to], idx[frm]] += v[comp]
    return A
Tt = mat(T, 0) * c / dz          # time displacement per unit z, times c  -> diag(+1, -1)
Xx = mat(X, 1) / dz              # signed x displacement per unit z       -> entering '+': +1, entering '-': -1
Fa = mat(F, 2) * c / dz          # flips per unit z, times c              -> a
ok("time moves: diag(+1, -1) (counter + goes forward, counter - backward, each carried)", Tt == sp.diag(1, -1))
ok("switch moves: [[0, +1], [-1, 0]] (entering + displaces +, entering - displaces -)", Xx == sp.Matrix([[0, 1], [-1, 0]]))
ok("flip moves: a [[0,1],[1,0]]", sp.simplify(Fa - a*sp.Matrix([[0, 1], [1, 0]])) == sp.zeros(2))
# two consecutive moves: products of first moments, both orders
ok("two time moves in a row: +1 (same counter twice: (+)(+) or (-)(-))", Tt*Tt == sp.eye(2))
ok("two switches in a row: -1 (always enter one of each counter: (+ell)(-ell))", Xx*Xx == -sp.eye(2))
ok("two flips in a row: +a^2", sp.simplify(Fa*Fa - a**2*sp.eye(2)) == sp.zeros(2))
ok("time-then-switch cancels switch-then-time", Tt*Xx + Xx*Tt == sp.zeros(2))
ok("time-then-flip cancels flip-then-time", Tt*Fa + Fa*Tt == sp.zeros(2))
ok("switch-then-flip cancels flip-then-switch", Xx*Fa + Fa*Xx == sp.zeros(2))
# assemble: c dz Q = (-Tt dt + c Xx (-dx)... ) -- check against the 2+1 lines and the square
t, x, z = sp.symbols('t x z')
Q = sp.Matrix([sp.Function('Qp')(t, x, z), sp.Function('Qm')(t, x, z)])
L = lambda g: -Tt*g.diff(t) - c*Xx*g.diff(x) + Fa*g          # displacement d -> factor (1 - d.grad): first moment enters with a minus
lines = sp.Matrix([-Q[0].diff(t) + a*Q[1] - c*Q[1].diff(x), +Q[1].diff(t) + a*Q[0] + c*Q[0].diff(x)])
ok("these first moments reproduce the 2+1 lines (loss taken out)", sp.simplify(L(Q) - lines) == sp.zeros(2, 1))
ok("two steps (L applied twice) = dt^2 - c^2 dx^2 + a^2: Klein-Gordon along z, Lorentzian minus from the switches",
   sp.simplify(L(L(Q)) - (Q.diff(t, 2) - c**2*Q.diff(x, 2) + a**2*Q)) == sp.zeros(2, 1))
