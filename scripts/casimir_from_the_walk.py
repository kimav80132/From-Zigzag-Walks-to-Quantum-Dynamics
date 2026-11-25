"""casimir_from_the_walk.py  --  recreated 9 Nov 2026 (the original is not on this disk).
Check numbers are the ones the companion document cites: 1a, 1b, 2, 3, 4.  The original had 7 checks; this recreation
prints 9 lines (2 and 3 are each split into three), covering the same statements.

 1a  ln det(1 - zP) = - sum_n z^n Tr P^n / n for the walk's step matrix P on a ring.
 1b  Tr P^n / n is one count per closed loop divided by its symmetry: closed n-step walks on the ring, enumerated,
     grouped into loops (a loop = a closed walk up to its starting point).
 2   mode sums with a smooth cutoff: turn-matrix wall, omega = (n+1/2) pi/L  ->  -pi/24L ;
     half-turn wall, omega = n pi/L  ->  +pi/12L.   (hbar = c = 1; E = - sum of the mode frequencies, L-linear part removed.)
 3   the kernel route, massless scalar between walls at 0 and L, G0 = -(1/4pi) ln|(x-x')^2 - (t-t')^2|:
     on the direct images (d_t d_t' + d_x d_x') G0 = -1/(pi u^2), u = 2nL; on the reflected images exactly 0;
     <T00> = -pi/(24 L^2), E = -pi/(24 L).
 4   the force: <T_xx> at the wall = -pi/(24 L^2) = -dE/dL.
"""
import itertools, numpy as np, sympy as sp, mpmath as mp
ok=True
def report(name,cond,extra=""):
    global ok; ok&=bool(cond); print(f"{name}: {'PASS' if cond else 'FAIL'}  {extra}")

# ---- 1a, 1b : the 1/n of a loop is counting
M=7; P=np.zeros((M,M))
for i in range(M): P[i,(i+1)%M]=0.5; P[i,(i-1)%M]=0.5          # Part II's step on a ring of 7 sites
z=0.8
lhs=np.log(np.linalg.det(np.eye(M)-z*P))
rhs=-sum(z**n*np.trace(np.linalg.matrix_power(P,n))/n for n in range(1,400))
report("1a. ln det(1 - zP) = - sum_n z^n Tr P^n / n",abs(lhs-rhs)<1e-12,f"({lhs:.12f} both sides, z={z})")
good=True; shown=[]
for n in (4,6,7,8,9):
    walks=[]                                                    # closed n-step walks: (start, steps)
    for start in range(M):
        for steps in itertools.product((+1,-1),repeat=n):
            if sum(steps)%M==0: walks.append((start,steps))
    # a loop = the set of walks obtained from one another by moving the starting point along the walk
    seen=set(); weighted=0.0
    for start,steps in walks:
        orbit=set(); pos=start
        for k in range(n):
            orbit.add((pos,steps[k:]+steps[:k])); pos=(pos+steps[k])%M
        key=min(orbit)
        if key in seen: continue
        seen.add(key); weighted+=len(orbit)/n                   # 1/symmetry: a loop traversed r times has n/r starting points
    count=len(walks)
    good&=np.isclose(count*0.5**n,np.trace(np.linalg.matrix_power(P,n))) and np.isclose(weighted,count/n)
    shown.append(f"n={n}: {count} closed walks, {len(seen)} loops")
report("1b. Tr P^n / n = sum over closed loops of (weight / symmetry)",good,"("+"; ".join(shown)+")")

# ---- 2 : modes with a smooth cutoff
mp.mp.dps=30
def E_modes(L,eps,shift):
    # E = - sum_n omega_n exp(-eps omega_n), omega_n = (n+shift) pi / L  -- closed form of the geometric sum, no asymptotics
    s=eps*mp.pi/L; q=mp.e**(-s)
    tot=q**shift*(shift+(1-shift)*q)/(1-q)**2                    # sum_{n>=0} (n+shift) q^(n+shift)
    return -(mp.pi/L)*tot
def finite_part(L,eps,shift): return E_modes(L,eps,shift)+L/(mp.pi*eps**2)   # remove the part proportional to L (there without walls)
L=mp.mpf(1)
vals=[finite_part(L,mp.mpf(e),mp.mpf(1)/2) for e in ('0.1','0.01','0.001')]
report("2a. turn-matrix wall, omega=(n+1/2)pi/L: finite part -> -pi/24L",abs(vals[-1]+mp.pi/24)<1e-6,"("+", ".join(mp.nstr(v,8) for v in vals)+f" ; -pi/24 = {mp.nstr(-mp.pi/24,8)})")
vals=[finite_part(L,mp.mpf(e),mp.mpf(0)) for e in ('0.1','0.01','0.001')]
report("2b. half-turn wall, omega=n pi/L: finite part -> +pi/12L",abs(vals[-1]-mp.pi/12)<1e-6,"("+", ".join(mp.nstr(v,8) for v in vals)+f" ; +pi/12 = {mp.nstr(mp.pi/12,8)})")
# brute sum as a cross-check of the closed form
brute=-sum((n+0.5)*np.pi*np.exp(-0.05*(n+0.5)*np.pi) for n in range(4000))
report("2c. closed form of the cut-off sum against the term-by-term sum",abs(brute-float(E_modes(1,mp.mpf('0.05'),mp.mpf(1)/2)))<1e-9)

# ---- 3 : the kernel and the square rule
t,tp,x,xp,u,Ls=sp.symbols('t tp x xp u L',real=True)
def G0(dx,dt): return -sp.log(dx**2-dt**2)/(4*sp.pi)
op=lambda G: sp.diff(G,t,tp)+sp.diff(G,x,xp)
direct=sp.simplify(op(G0(x-xp-u,t-tp)).subs({tp:t,xp:x}))
reflected=sp.simplify(op(G0(x+xp-u,t-tp)))                       # before coincidence: identically zero
report("3a. direct image at distance u: (d_t d_t' + d_x d_x') G0 = -1/(pi u^2)",sp.simplify(direct+1/(sp.pi*u**2))==0,f"({direct})")
report("3b. reflected image: (d_t d_t' + d_x d_x') G0 = 0 identically, at every x, x', t, t'",reflected==0)
n=sp.symbols('n',integer=True,positive=True)
T00=sp.Rational(1,2)*2*sp.summation(-1/(sp.pi*(2*n*Ls)**2),(n,1,sp.oo))    # 1/2 (..)G at coincidence, images n = +-1, +-2, ...
E=sp.simplify(T00*Ls)
report("3c. <T00> = -pi/(24 L^2), uniform; E(L) = -pi/(24 L)",sp.simplify(T00+sp.pi/(24*Ls**2))==0 and sp.simplify(E+sp.pi/(24*Ls))==0,f"(T00 = {sp.simplify(T00)}, E = {E})")

# ---- 4 : the force on the wall
# massless in 1+1: T_xx = 1/2 (d_t d_t' + d_x d_x') G at coincidence, the same bilinear as T00; the reflected images add 0 at the wall too (3b)
Txx=T00
report("4.  <T_xx> at the wall = -pi/(24 L^2) = -dE/dL",sp.simplify(Txx-(-sp.diff(E,Ls)))==0 and sp.simplify(Txx-(-sp.pi/(24*Ls**2)))==0,f"(-dE/dL = {sp.simplify(-sp.diff(E,Ls))})")
print("ALL PASS" if ok else "SOME FAIL")
raise SystemExit(0 if ok else 1)
