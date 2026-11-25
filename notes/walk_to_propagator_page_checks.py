"""Checks behind the page "From a random walk to the Dirac propagator, in 1+1" (25 Nov 2026).   Master: III.1, III.2, III.5a, III.14, III.15, IV.10.
 1. the step rule, expanded to first order in the step, gives the two lines  c d_x P_e + e d_t P_e = a (P_-e - P_e)
 2. the walk itself (probabilities, clock x), run on a lattice, against III.14's box: P_- = e^(-ax/c) (a/2) I_0(zeta)
 3. reading in t: the sum over zigzags, each corner weighing -a dt (forward into backward) or +a dt (backward into forward), on a lattice, against
        Q_- = -(a/2c) J_0(zeta'),    Q_+ bulk = -(a/2c) (ct+x)/sqrt(c^2t^2-x^2) J_1(zeta'),    zeta' = (a/c) sqrt(c^2 t^2 - x^2)
 4. corner by corner: the zigzags with 2n+1 corners from (0,0) to (t,x) number (1/2c) r^n l^n / (n!)^2 per unit x, r = (ct+x)/2c, l = (ct-x)/2c
    (time spent going right, going left); each weighs -a (-a^2)^n; the sum over n is the series of -(a/2c) J_0
 5. on one wavelength: the walk's generator G = -ick sigma_z + a T has G G = -(c^2k^2 + a^2); e^(tG) = cos(wt) + (G/w) sin(wt); the lattice kernel of 3, summed
    against e^(-ikx), gives the same column
 6. usual symbols: e^(tG) = e^(-iHt/hbar), H = c sigma_z p - m c^2 sigma_y; = e^(-iEt) Pi_+ + e^(+iEt) Pi_-, Pi_+- = (1 +- H/E)/2
 7. Feynman's rule: K_F = theta(t) e^(-iEt) Pi_+ - theta(-t) e^(+iEt) Pi_- = (i d_t + H) Delta_F, Delta_F = int dk/2pi e^(ikx) e^(-iw|t|)/(2w);
    Delta_F = (1/2 pi c) K_0((a/c) sqrt(x^2 - c^2t^2)) outside the cone, -(i/4c) H0^(2)((a/c) sqrt(c^2t^2 - x^2)) inside;
    the evolution kernel is (i d_t + H) D with D = (1/2c) J_0(zeta') inside: Feynman's rule keeps half of the J_0 and adds the Y_0
"""
import numpy as np, sympy as sp
from scipy.special import i0, j0, j1, k0, hankel2, y0
from scipy.integrate import quad
from math import factorial
fails=[]
def ok(m,c):
    print(("PASS  " if c else "FAIL  ")+m)
    if not c: fails.append(m)
# 1
t,x,dx,a,c=sp.symbols('t x dx a c',positive=True); res=True
for e in (1,-1):
    P=sp.Function('P'); Pm=sp.Function('Pm')
    rhs=(1-a*dx/c)*P(t-e*dx/c,x)+a*dx/c*Pm(t,x)
    lhs=P(t,x+dx)
    d=sp.series(lhs-rhs,dx,0,2).removeO().doit()
    line=sp.simplify(d*c/dx-(c*sp.diff(P(t,x),x)+e*sp.diff(P(t,x),t)-a*(Pm(t,x)-P(t,x))))
    res&=(line==0)
ok("1. step rule to first order in the step:  c d_x P_e + e d_t P_e = a (P_-e - P_e), both senses",res)
# 2 the walk itself: clock x, exact probabilities on a lattice
av,cv=1.3,2.0; X=1.2; N=1200; h=X/N; p=av*h/cv
Pp=np.zeros(2*N+3); Pm_=np.zeros(2*N+3); Pp[N+1]=1.0
for n in range(N):
    Pp=np.roll(Pp,1); Pm_=np.roll(Pm_,-1)          # advance: t -> t + e h/c
    Pp,Pm_=(1-p)*Pp+p*Pm_,(1-p)*Pm_+p*Pp           # flip
idx=np.arange(2*N+3)-(N+1); tt=idx*h/cv
errs=[]
for tq in (-0.3,0.0,0.25):
    j=np.argmin(np.abs(tt-tq)); j+= (j-(N+1)+N)%2   # reachable parity
    dens=Pm_[j]/(2*h/cv); zeta=av/cv*np.sqrt(X**2-cv**2*tt[j]**2)
    errs.append(abs(dens/(np.exp(-av*X/cv)*av/2*i0(zeta))-1))
ok("2. the walk on a lattice (a = 1.3, c = 2, x = 1.2): P_- against e^(-ax/c) (a/2) I_0(zeta) at t = -0.3, 0, 0.25: off by %s"%np.round(errs,4), max(errs)<0.01)
# 3 reading in t
T=2.0; N=2000; dt=T/N; h=cv*dt
Qp=np.zeros(2*N+3); Qm=np.zeros(2*N+3); Qp[N+1]=1.0
for n in range(N):
    Qp=np.roll(Qp,1); Qm=np.roll(Qm,-1)
    Qp,Qm=Qp+av*dt*Qm,Qm-av*dt*Qp
xx=(np.arange(2*N+3)-(N+1))*h; e3=[]
for xq in (-1.5,0.0,2.2):
    j=np.argmin(np.abs(xx-xq)); j+=(j-(N+1)+N)%2
    z=av/cv*np.sqrt(cv**2*T**2-xx[j]**2)
    e3.append(abs(Qm[j]/(2*h)-(-(av/(2*cv))*j0(z))))
    e3.append(abs(Qp[j]/(2*h)-(-(av/(2*cv))*(cv*T+xx[j])/np.sqrt(cv**2*T**2-xx[j]**2)*j1(z))))
ok("3. reading in t (t = 2): the lattice sum over zigzags against -(a/2c) J_0 and -(a/2c)(ct+x)/sqrt(c^2t^2-x^2) J_1 at x = -1.5, 0, 2.2: largest miss %.4f"%max(e3), max(e3)<0.004)
# 4 corner by corner: count lattice zigzags with R = 2n+1 corners ending backward
from math import comb
Nn=40000; dtn=T/Nn; kR=26000; kL=Nn-kR                     # kR ticks going right, kL going left; x = (kR-kL) c dt
r=kR*dtn; l=kL*dtn; e4=[]
for n in (0,1,2,3):
    # n+1 right segments (each >= 1 tick) and n+1 left segments (each >= 1 tick): compositions
    count=comb(kR-1,n)*comb(kL-1,n)
    lattice=count*dtn**(2*n+1)/(2*cv*dtn)                # x (weight a dt per corner stripped) per unit x: reachable sites are 2 c dt apart
    cont=(1/(2*cv))*r**n*l**n/factorial(n)**2
    e4.append(abs(lattice/cont-1))
series=sum(-(av/(2*cv))*(-av**2*r*l)**n/factorial(n)**2 for n in range(40)); zz=av/cv*np.sqrt(cv**2*T**2-((kR-kL)*cv*dtn)**2)
ok("4. zigzags with 2n+1 corners: the number is (1/2c) r^n l^n/(n!)^2 per unit x (lattice count off by %s for n = 0..3); weighted -a(-a^2)^n and summed: %.6f = -(a/2c) J_0 = %.6f"
   %(np.round(e4,5),series,-(av/(2*cv))*j0(zz)), max(e4)<0.001 and abs(series+(av/(2*cv))*j0(zz))<1e-12)
# 5 one wavelength
k=0.9; w=np.sqrt(cv**2*k**2+av**2); sz=np.diag([1.,-1]); Tm=np.array([[0,1],[-1,0.]]); G=-1j*cv*k*sz+av*Tm
from scipy.linalg import expm
col=np.array([np.sum(Qp*np.exp(-1j*k*xx)),np.sum(Qm*np.exp(-1j*k*xx))])
ok("5. G G = -(c^2k^2 + a^2); e^(tG) = cos(wt) + (G/w) sin(wt); the lattice kernel summed against e^(-ikx) gives e^(tG)(1,0) (miss %.4f)"%np.max(np.abs(col-expm(T*G)[:,0])),
   np.allclose(G@G,-w*w*np.eye(2)) and np.allclose(expm(T*G),np.cos(w*T)*np.eye(2)+G/w*np.sin(w*T)) and np.max(np.abs(col-expm(T*G)[:,0]))<0.004)
# 6 usual symbols (hbar = 1, m c^2 = a)
sy=np.array([[0,-1j],[1j,0]]); H=cv*k*sz-av*sy; Pip=(np.eye(2)+H/w)/2; Pim=(np.eye(2)-H/w)/2
ok("6. e^(tG) = e^(-iHt); = e^(-iEt) Pi_+ + e^(+iEt) Pi_-; Pi_+ Pi_+ = Pi_+, Pi_+ Pi_- = 0",
   np.allclose(expm(T*G),expm(-1j*H*T)) and np.allclose(expm(-1j*H*T),np.exp(-1j*w*T)*Pip+np.exp(1j*w*T)*Pim) and np.allclose(Pip@Pip,Pip) and np.allclose(Pip@Pim,0))
# 7 Feynman
f=lambda s: np.exp(-1j*w*abs(s))/(2*w); d=1e-6; r7=True
for s in (0.8,-0.8):
    lhs=(1j*(f(s+d)-f(s-d))/(2*d))*np.eye(2)+H*f(s)
    want=np.exp(-1j*w*s)*Pip if s>0 else -np.exp(1j*w*s)*Pim
    r7&=np.allclose(lhs,want,atol=1e-6)
xo=0.7
out=quad(lambda q: 1/(2*np.sqrt(cv**2*q**2+av**2))/np.pi,0,np.inf,weight='cos',wvar=xo)[0]
tin=1.1; eps=0.02; tz=tin-1j*eps
ins=quad(lambda e: np.real(np.exp(-1j*av*tz*np.cosh(e))),0,12,limit=2000)[0]+1j*quad(lambda e: np.imag(np.exp(-1j*av*tz*np.cosh(e))),0,12,limit=2000)[0]
ins=ins/(2*np.pi*cv)
ok("7. K_F = (i d_t + H) Delta_F on one wavelength, both signs of t; Delta_F at t = 0: %.6f = K_0(a x/c)/(2 pi c) = %.6f; inside (x = 0, t slightly lowered into the complex plane): matches -(i/4c) H0^(2)(a t) to %.1e;"
   " -(i/4c) H0^(2) = -(i/4c) J_0 - (1/4c) Y_0"%(out,k0(av*xo/cv)/(2*np.pi*cv),abs(ins-(-1j/(4*cv))*hankel2(0,av*tz))),
   r7 and abs(out-k0(av*xo/cv)/(2*np.pi*cv))<1e-6 and abs(ins-(-1j/(4*cv))*hankel2(0,av*tz))<1e-6 and np.isclose(hankel2(0,1.3),j0(1.3)-1j*y0(1.3)))
# 8 renaming the backward counter
Dm=np.diag([1,1j]); ren=Dm@(av*Tm)@np.linalg.inv(Dm)
ok("8. renaming the backward counter Q_- -> i Q_-: the corner matrix a T becomes -i a sigma_x: every corner weighs -i a dt", np.allclose(ren,-1j*av*np.array([[0,1],[1,0]])))
# 9 Delta_F at general points, time lowered slightly into the complex plane
def DF(tq,xq,eps=0.02):
    e=np.linspace(-9.5,9.5,800001); tz=tq-1j*eps
    f=np.exp(1j*av/cv*xq*np.sinh(e)-1j*av*tz*np.cosh(e)); return np.trapezoid(f,e)/(4*np.pi*cv), tz
r9=[]
for (tq,xq) in ((1.1,0.9),(0.9,-1.2)):          # inside: c t = 2.2 > 0.9 ; and c t = 1.8 > 1.2
    num,tz=DF(tq,xq); r9.append(abs(num-(-1j/(4*cv))*hankel2(0,av/cv*np.sqrt(cv**2*tz**2-xq**2))))
from scipy.special import kv
for (tq,xq) in ((0.3,1.7),(0.5,-2.4)):          # outside
    num,tz=DF(tq,xq); r9.append(abs(num-kv(0,av/cv*np.sqrt(xq**2-cv**2*tz**2))/(2*np.pi*cv)))
ok("9. Delta_F at two points inside the cone against -(i/4c) H0^(2)(zeta) and two outside against K_0(zeta')/(2 pi c): largest miss %.1e"%max(r9), max(r9)<1e-6)
print("\nALL PASS" if not fails else "\nFAILED %d"%len(fails))
