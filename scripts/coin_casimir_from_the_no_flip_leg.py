"""Casimir document, "The displayed bridge: the coin's number from its own kernel".   Recreated 9 Nov 2026 (doc rev 68).  hbar = c = 1.

Walls at 0 and L: psi_-(0) = r_0 psi_+(0), psi_+(L) = r_L psi_-(L).  Unfolded, phi(s) = psi_+(s) on [0, L), r_L psi_-(2L - s) on [L, 2L):
one right-mover on a circle of length 2L with phi(s + 2L) = r_L r_0 phi(s); r_L r_0 = (wall matrix)^2 = -1 for the turn matrix, +1 for
the half-turn.
  1. the modes e^{2ikL} = r_L r_0: (n + 1/2) pi/L for -1, n pi/L for +1; the wall matrices' squares.
  2. the no-flip kernel G_0(u) = 1/(2 pi i u) with the walls by images: sum_m (r_L r_0)^m G_0(u + 2mL)
        = (pi/2L)/(2 pi i sin(pi u/2L))  for -1,   (pi/2L)/(2 pi i tan(pi u/2L))  for +1.
  3. the energy bilinear of a right-mover, -i psi^dag d_s psi, point-split on G: the m = 0 image is the bulk; every other image gives
     (r_L r_0)^m / (2 pi (2mL)^2).  Summed: <T> = -pi/(48 L^2) (turn matrix), +pi/(24 L^2) (half-turn); the same from the closed forms of 2.
  4. over the circle of length 2L: E = -pi/(24 L) = -0.130900 and +pi/(12 L) = +0.261799 at L = 1 - the wall table - and the same two numbers
     from the modes of 1 with a smooth cutoff (the filled lower half of the walker's spectrum, -sum of the positive modes)."""
import sys
import mpmath as mp, sympy as sp
mp.mp.dps = 30
fails = 0
def ok(name, cond):
    global fails
    print(("PASS  " if cond else "FAIL  ") + name); fails += (not cond)
T = sp.Matrix([[0, 1], [-1, 0]]); Sx = sp.Matrix([[0, 1], [1, 0]])
L = mp.mpf(1)
# 1
k_half = [(n + mp.mpf(1)/2)*mp.pi/L for n in range(-3, 4)]; k_int = [n*mp.pi/L for n in range(-3, 4)]
ok("1. (turn matrix)^2 = -1, (half-turn)^2 = +1; e^{2ikL} = -1 at k = (n + 1/2) pi/L (%.0e), = +1 at k = n pi/L (%.0e)"
   % (max(abs(mp.e**(2j*k*L) + 1) for k in k_half), max(abs(mp.e**(2j*k*L) - 1) for k in k_int)),
   T*T == -sp.eye(2) and Sx*Sx == sp.eye(2) and max(abs(mp.e**(2j*k*L) + 1) for k in k_half) < 1e-25 and max(abs(mp.e**(2j*k*L) - 1) for k in k_int) < 1e-25)
# 2: image sums, summed symmetrically (m and -m together), against the closed forms
G0 = lambda u: 1/(2j*mp.pi*u)
def images(u, sgn):
    return G0(u) + mp.nsum(lambda m: sgn**int(m)*(G0(u + 2*m*L) + G0(u - 2*m*L)), [1, mp.inf])
worst = mp.mpf(0)
for u in (mp.mpf('0.13'), mp.mpf('0.7'), mp.mpf('1.9')):
    worst = max(worst, abs(images(u, -1) - (mp.pi/(2*L))/(2j*mp.pi*mp.sin(mp.pi*u/(2*L)))),
                abs(images(u, +1) - (mp.pi/(2*L))/(2j*mp.pi*mp.tan(mp.pi*u/(2*L)))))
ok("2. sum_m (r_L r_0)^m G_0(u + 2mL) = (pi/2L)/(2 pi i sin(pi u/2L)) for -1 and (pi/2L)/(2 pi i tan(pi u/2L)) for +1, at three u (%.0e)" % worst, worst < 1e-20)
# 3: -i d_u of each image at u -> 0
u_ = sp.symbols('u'); m_, L_ = sp.symbols('m L', positive=True)
img = sp.simplify((-sp.I*sp.diff(1/(2*sp.pi*sp.I*(u_ + 2*m_*L_)), u_)).subs(u_, 0))
T_turn = mp.nsum(lambda m: 2*(-1)**int(m)/(2*mp.pi*(2*m*L)**2), [1, mp.inf]); T_half = mp.nsum(lambda m: 2/(2*mp.pi*(2*m*L)**2), [1, mp.inf])
# the same from the closed forms: -i d_u [G - G_0] at u = 0
uu, LL = sp.symbols('u L', positive=True)
c_turn = sp.limit(-sp.I*sp.diff((sp.pi/(2*LL))/(2*sp.pi*sp.I*sp.sin(sp.pi*uu/(2*LL))) - 1/(2*sp.pi*sp.I*uu), uu), uu, 0)
c_half = sp.limit(-sp.I*sp.diff((sp.pi/(2*LL))/(2*sp.pi*sp.I*sp.tan(sp.pi*uu/(2*LL))) - 1/(2*sp.pi*sp.I*uu), uu), uu, 0)
ok("3. -i d_u of the m-th image at u = 0 is %s (times (r_L r_0)^m); summed over m != 0: <T> = %s (turn matrix), %s (half-turn), against "
   "-pi/(48 L^2) = %s and +pi/(24 L^2) = %s; from the closed forms of 2: %s and %s"
   % (img, mp.nstr(T_turn, 12), mp.nstr(T_half, 12), mp.nstr(-mp.pi/48, 12), mp.nstr(mp.pi/24, 12), c_turn, c_half),
   sp.simplify(img - 1/(2*sp.pi*(2*m_*L_)**2)) == 0 and abs(T_turn + mp.pi/48) < 1e-20 and abs(T_half - mp.pi/24) < 1e-20
   and sp.simplify(c_turn + sp.pi/(48*LL**2)) == 0 and sp.simplify(c_half - sp.pi/(24*LL**2)) == 0)
# 4: energies, and the modes with a smooth cutoff
E_turn, E_half = 2*L*T_turn, 2*L*T_half
def mode_energy(shift, eps):                    # -(sum of the positive modes) with the cutoff e^{-eps k}, bulk 1/eps^2 term removed (it is proportional to L)
    s = mp.nsum(lambda n: (n + shift)*mp.pi/L*mp.e**(-eps*(n + shift)*mp.pi/L), [0, mp.inf])
    return -(s - L/(mp.pi*eps**2))
m_turn = [mode_energy(mp.mpf(1)/2, e) for e in (mp.mpf('0.02'), mp.mpf('0.01'), mp.mpf('0.005'))]
m_half = [mode_energy(mp.mpf(0), e) for e in (mp.mpf('0.02'), mp.mpf('0.01'), mp.mpf('0.005'))]
ok("4. E = 2L <T>: %s (turn matrix) and %s (half-turn) at L = 1, against -pi/24 = -0.130900 and +pi/12 = +0.261799; the modes of 1 with a "
   "smooth cutoff, bulk removed: %s and %s at cutoff 0.005 (converging as the cutoff squared)"
   % (mp.nstr(E_turn, 9), mp.nstr(E_half, 9), mp.nstr(m_turn[-1], 9), mp.nstr(m_half[-1], 9)),
   abs(E_turn + mp.pi/24) < 1e-20 and abs(E_half - mp.pi/12) < 1e-20 and abs(m_turn[-1] + mp.pi/24) < 1e-4 and abs(m_half[-1] - mp.pi/12) < 1e-4
   and abs(m_turn[0] + mp.pi/24) > abs(m_turn[-1] + mp.pi/24))
print("ALL PASS" if not fails else "%d FAILED" % fails); sys.exit(1 if fails else 0)
