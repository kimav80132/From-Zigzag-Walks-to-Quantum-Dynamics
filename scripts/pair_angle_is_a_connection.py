"""VI.2 (A)-(C): the angle of the pair (out-leg, return leg) is a connection, and only loops see it.
1+1 two-copy lattice of V.2 read in t: column per site (Q+, Q+', Q-, Q-'); one tick = shift (+ right, - left) then the coin
R = cos(a dt) 1 + sin(a dt) T on the senses, the same for both legs; the field = a turn of (leg, leg') by the angle theta on each link.
1. a GLOBAL turn J of the pair commutes with the free step: the coin and the shift act on both legs alike, so the absolute angle of the
   pair at an event is not a dynamical variable (the noise belongs to the event, VI.1 (B): it cannot tell the legs apart).
2. a LOCAL turn theta(x,t) does not commute with the free step: the mismatch is the difference of theta across one link. So a step must
   carry an angle per link, theta_link, to say how the pair at one event compares with the pair at the next: a connection.
3. gauge covariance, exactly on the lattice: U[theta_link + d theta] = R(theta_after) U[theta_link] R(theta_before)^{-1}. Readings turn
   locally; the event weight Q^2 + Q'^2 (III.15a, V.5 (B)) does not change at all.
4. only loops are seen: a pure-gauge link field (every plaquette sum zero) gives the free event weights to machine precision; a field with
   flux changes them. The event weight depends on the link angles only through plaquette sums = the flux F.
5. the return leg: the inverse step is the transpose, so coming back along the same links undoes the angles exactly (loop angle 0); out
   one way and back another encloses the flux: an event behind a flux line sees the flux with period 2 pi (Aharonov-Bohm on the lattice)."""
import numpy as np
ok = lambda n, b: print(("PASS  " if b else "FAIL  ") + n)
rng = np.random.default_rng(3)
L, N, a, dt = 20, 9, 1.0, 0.3
T2 = np.array([[0, 1], [-1, 0.]]); Rc = np.cos(a*dt)*np.eye(2) + np.sin(a*dt)*T2        # coin on the senses
J2 = np.array([[0, -1], [1, 0.]])                                                       # quarter-turn of (leg, leg')
Rl = lambda th: np.cos(th)*np.eye(2) + np.sin(th)*J2                                    # turn of the pair by th
dim = 4*L
idx = lambda x, e, l: 4*(x % L) + (0 if e == +1 else 2) + l                             # (site, sense, leg)
# free shift and coin
S = np.zeros((dim, dim))
for x in range(L):
    for l in (0, 1):
        S[idx(x+1, +1, l), idx(x, +1, l)] = 1; S[idx(x-1, -1, l), idx(x, -1, l)] = 1
C = np.kron(np.eye(L), np.kron(Rc, np.eye(2)))                                           # coin on senses, identity on legs
Jg = np.kron(np.eye(2*L), J2)                                                            # global quarter-turn of every pair
U0 = C @ S
ok("1. the global turn of the pair commutes with the free step (coin and shift act on both legs alike): the absolute angle is not dynamical",
   np.allclose(Jg @ U0, U0 @ Jg))
# local turn
def Rloc(theta):                                                                         # theta[x]: angle of the pair at site x (both senses)
    M = np.zeros((dim, dim))
    for x in range(L):
        for e in (+1, -1):
            b = idx(x, e, 0); M[b:b+2, b:b+2] = Rl(theta[x])
    return M
theta = rng.normal(size=L)
mis = Rloc(theta).T @ U0 @ Rloc(theta) - U0
ok("2. a local turn does not commute with the free step; the mismatch is first order in the difference of the angle across a link",
   np.linalg.norm(mis) > 1e-3 and abs(np.linalg.norm(Rloc(theta/100).T @ U0 @ Rloc(theta/100) - U0)/np.linalg.norm(Rloc(theta/200).T @ U0 @ Rloc(theta/200) - U0) - 2) < 0.02)
# step with link angles: a link (x -> x+e) carries an angle th[x, e]; the pair is turned by it as it crosses the link
def step_with_links(th):
    Sg = np.zeros((dim, dim))
    for x in range(L):
        for e in (+1, -1):
            b0 = idx(x, e, 0); b1 = idx(x+e, e, 0)
            Sg[b1:b1+2, b0:b0+2] = Rl(th[(x, e)])
    return C @ Sg
th = {(x, e): rng.normal() for x in range(L) for e in (+1, -1)}
th_g = {(x, e): th[(x, e)] + theta[(x+e) % L] - theta[x] for x in range(L) for e in (+1, -1)}   # gauge-shifted links
U_th, U_g = step_with_links(th), step_with_links(th_g)
ok("3. gauge covariance on the lattice: U[links + d theta] = R(theta) U[links] R(theta)^T exactly",
   np.allclose(U_g, Rloc(theta) @ U_th @ Rloc(theta).T))
A = idx(0, +1, 0); v = np.zeros(dim); v[A] = 1
weight = lambda U: np.array([np.sum((np.linalg.matrix_power(U, N) @ v)[4*x:4*x+4]**2) for x in range(L)])   # Q^2 + Q'^2 per site
ok("3'. the event weights Q^2 + Q'^2 at every absorber are unchanged by the local turn", np.allclose(weight(U_th), weight(U_g)))
# 4. pure gauge vs flux. On this lattice a plaquette is: + link x->x+1 at one tick, then the - link back x+1->x at the next tick, etc.
# a pure-gauge field is th = theta[x+e] - theta[x]; its loop sums vanish. Compare with the free step.
th_pure = {(x, e): theta[(x+e) % L] - theta[x] for x in range(L) for e in (+1, -1)}
th_flux = {(x, e): (0.7 if (x, e) == (5, +1) else 0.0) for x in range(L) for e in (+1, -1)}       # one link carries 0.7: flux through every loop crossing it once
ok("4. a pure-gauge link field gives the free event weights to machine precision; one link carrying an angle changes them",
   np.allclose(weight(step_with_links(th_pure)), weight(U0)) and not np.allclose(weight(step_with_links(th_flux)), weight(U0)))
# 5. the return leg undoes the angle: U^T carries R(th)^T = R(-th) on each link. Loop out-and-back along the same link: angle 0.
Ut = U_th.T
ok("5a. the inverse step is the transpose and carries the opposite angle on every link: out and back along the same path is the free loop, "
   "(U^T U)_{AA} = 1 with any link field", abs((Ut @ U_th)[A, A] - 1) < 1e-12)
# Aharonov-Bohm: a time-independent angle alpha on the link (5,+) only (and -alpha on (6,-), i.e. a pure flux line between sites 5 and 6
# in the sense that going right through it gives +alpha and going left back through it gives -alpha... that is a pure gauge).
# A flux is a link angle that is NOT undone by the reverse link: alpha on (5,+) and +alpha on (6,-) (same sense of rotation both ways).
def ab_weights(alpha):
    thA = {(x, e): 0.0 for x in range(L) for e in (+1, -1)}
    thA[(5, +1)] = alpha; thA[(6, -1)] = alpha
    return weight(step_with_links(thA))
w0, w1, w2 = ab_weights(0.0), ab_weights(1.3), ab_weights(1.3 + 2*np.pi)
ok("5b. an angle that both crossings of a line pick up with the same sign is a flux: the event weights change with it and are periodic "
   "in it with period 2 pi (Aharonov-Bohm), while a pure gauge (opposite signs) does nothing",
   not np.allclose(w0, w1) and np.allclose(w1, w2))
