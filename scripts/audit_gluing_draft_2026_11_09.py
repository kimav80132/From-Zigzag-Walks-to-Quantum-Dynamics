"""Recheck of drafts/gluing_the_route_to_itself_2026_11_07.md, 9 Nov 2026. Exact algebra / finite sums only."""
import numpy as np, sympy as sp
ok=lambda b: "PASS" if b else "FAIL"
a=1.0
print("=== 1. the winding series of section 1 in the Euclidean regime of section 3 ===")
beta=3.0; x=np.exp(beta*a)                      # the growing slot's eigenvalue of U = e^{-beta H}
part=[sum((-1)**(n+1)*x**n/n for n in range(1,N+1)) for N in (5,10,20)]
print("  sum_n (-1)^(n+1) x^n/n, x=e^{beta a}=%.2f, partial sums N=5,10,20: %s  -> diverges; log(1+x) = %.6f"%(x,np.array(part),np.log1p(x)))
conv=beta*a+sum((-1)**(n+1)*np.exp(-n*beta*a)/n for n in range(1,60))
print("  refilled form  beta*a + sum_n (-1)^(n+1) e^{-n beta a}/n = %.12f  vs log(1+e^{beta a}) = %.12f :"%(conv,np.log1p(x)),ok(abs(conv-np.log1p(x))<1e-12))
tot=beta*a+2*sum((-1)**(n+1)*np.exp(-n*beta*a)/n for n in range(1,60))
print("  whole rest mode: log Z = beta*a + 2 sum_n (-1)^(n+1) e^{-n beta a}/n :",ok(abs(tot-np.log(2+2*np.cosh(beta*a)))<1e-12),
      " (after the growing slot is filled, BOTH remaining series are windings that cost e^{-n beta a}: a particle and a hole)")
print()
print("=== 2. tr e^{s M_k} = 2 cosh(omega s) as a sum over turn placements: proof by Laplace transform in s ===")
p,k,lam,s=sp.symbols('p k lambda s',positive=True); m=sp.symbols('m',integer=True,nonnegative=True)
# 2m turns, start +: (m+1) stretches e^{-k u} and m stretches e^{+k u}; start -: the mirror. Laplace in s of an ordered convolution = product
term=lam**(2*m)*((p+k)**(-(m+1))*(p-k)**(-m)+(p-k)**(-(m+1))*(p+k)**(-m))
ratio=lam**2/(p**2-k**2); pref=2*p/(p**2-k**2)
is_geom=all(sp.simplify(term.subs(m,mm)-pref*ratio**mm)==0 for mm in range(0,6))      # each term = prefactor * ratio^m
target=2*p/(p**2-k**2-lam**2)                    # Laplace transform of 2 cosh(s sqrt(k^2+lam^2))
print("  each 2m-turn term = [2p/(p^2-k^2)] * [lambda^2/(p^2-k^2)]^m (m=0..5):",ok(is_geom))
print("  geometric sum = 2p/(p^2-k^2-lambda^2) = Laplace transform of 2cosh(s sqrt(k^2+lambda^2)):",
      ok(sp.simplify(pref/(1-ratio)-target)==0 and sp.simplify(sp.laplace_transform(2*sp.cosh(s*sp.sqrt(k**2+lam**2)),s,p,noconds=True)-target)==0))
two=sp.integrate((s-sp.Symbol('u'))*2*sp.cosh(k*(s-2*sp.Symbol('u'))),(sp.Symbol('u'),0,s))
print("  the two-turn term done directly, integral of (s-u) 2cosh(k(s-2u)) du = s sinh(ks)/k:",ok(sp.simplify(two-s*sp.sinh(k*s)/k)==0))
print()
print("=== 3. section 6b (ii): is prod_k det(1+e^{s M_k}) the count of non-crossing lines in the real (x,t) plane? ===")
# a line ensemble in the real plane has states (t,eps) on a clock cut; its transfer operator in t-frequency nu is M(nu)=[[i nu,a],[a,-i nu]]
nus=np.linspace(-6,6,240001); dnu=nus[1]-nus[0]
rate=np.where(np.abs(nus)<a,np.sqrt(np.clip(a*a-nus**2,0,None)),0.0)
print("  real-plane line count: growth per unit clock and per unit t-length = int_{|nu|<a} sqrt(a^2-nu^2) dnu/2pi = %.6f  (a^2/4 = %.6f) : bounded"%((rate.sum()*dnu)/(2*np.pi),a*a/4))
for K in (5,20,80):
    ks=np.linspace(-K,K,400001); print("  mode sum over |k|<%d: int sqrt(a^2+k^2) dk/2pi = %.3f  (grows with the cutoff)"%(K,np.sqrt(a*a+ks**2).sum()*(ks[1]-ks[0])/(2*np.pi)))
print("  -> two different objects. The draft's product over k is the trace over the t-side modes at imaginary time (equivalently the")
print("     x-transfer in the plane (x, imaginary t)); each FACTOR is a positive count along the walk's clock (section 6, 6c); the PRODUCT")
print("     is not shown to be a count of non-crossing lines in the real plane, and that count has a different, bounded rate.")
print()
print("=== 4. section 6c: the measure of the frame sum ===")
eta=sp.symbols('eta',real=True); aa,cc=sp.symbols('a c',positive=True)
dk_deta=sp.diff(aa*sp.sinh(eta)/cc,eta)
print("  dk/d eta =",dk_deta,"  -> the box sum over k weights frames by cosh(eta); III.16's frame sum is uniform in eta.")
print("  per unit length: -sum_k hbar omega_k = -(hbar a^2/2 pi c) * int cosh^2(eta) d eta : the invariant measure d eta times cosh^2(eta)")
