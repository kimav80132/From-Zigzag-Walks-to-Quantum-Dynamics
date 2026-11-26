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
# 10 step 11 in real numbers: Q_pm = A_pm cos kx - B_pm sin kx, column (A+, A-, B+, B-)
from scipy.linalg import expm as _expm
kk=0.9; K4=np.array([[0,0,-1,0],[0,0,0,1],[1,0,0,0],[0,-1,0,0]],float); T4=np.array([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]],float)
G4=-cv*kk*K4+av*T4; w4=np.sqrt(cv**2*kk**2+av**2); I4=np.eye(4)
# the matrices from the lines of step 6: d_t Q+ = -c d_x Q+ + a Q-, d_t Q- = +c d_x Q- - a Q+, with d_x (A,B) = k (-B,A)
xs_=sp.symbols('x_'); Ap,Am,Bp,Bm,ks,cs,as_=sp.symbols('Ap Am Bp Bm k c a')
Qp_=Ap*sp.cos(ks*xs_)-Bp*sp.sin(ks*xs_); Qm_=Am*sp.cos(ks*xs_)-Bm*sp.sin(ks*xs_)
r1=sp.expand(-cs*sp.diff(Qp_,xs_)+as_*Qm_); r2=sp.expand(cs*sp.diff(Qm_,xs_)-as_*Qp_)
rows=[r1.coeff(sp.cos(ks*xs_)),r2.coeff(sp.cos(ks*xs_)),-r1.coeff(sp.sin(ks*xs_)),-r2.coeff(sp.sin(ks*xs_))]
Gs=sp.Matrix([[sp.diff(r,v) for v in (Ap,Am,Bp,Bm)] for r in rows])
Kpage=sp.Matrix(K4.astype(int)); Tpage=sp.Matrix(T4.astype(int))
c10=(sp.simplify(Gs-(-cs*ks*Kpage+as_*Tpage))==sp.zeros(4,4)) and np.allclose(K4@K4,-I4) and np.allclose(T4@T4,-I4) and np.allclose(K4@T4+T4@K4,0) \
    and np.allclose(G4@G4,-w4**2*I4) and np.allclose(_expm(0.7*G4),np.cos(w4*0.7)*I4+G4/w4*np.sin(w4*0.7))
ok("10. one pair alone (the version in the usual symbols, step 11): the 4x4 matrix G on (A+,A-,B+,B-) is the one the lines of step 6 give; each of its two matrices squares to -1 and they anticommute; G G = -(c^2k^2+a^2); e^(tG) = cos(wt) + (G/w) sin(wt)",c10)
# 11 the same in the usual complex form: z_pm = A_pm + i B_pm, G -> -i c k sigma_z + a T = -i H
sz=np.diag([1.,-1.]); T2=np.array([[0,1],[-1,0]],float); G2=-1j*cv*kk*sz+av*T2
v0=np.array([0.3,-1.1,0.8,0.5]); v1=_expm(0.7*G4)@v0; z0=v0[:2]+1j*v0[2:]; z1=_expm(0.7*G2)@z0
Hm=cv*kk*sz-av*np.array([[0,-1j],[1j,0]])
ok("11. the same version, step 12: with z = A + iB the real 4x4 turn is the complex 2x2 one, and G = -iH (H = c sigma_z p - m c^2 sigma_y, hbar = 1)", np.allclose(v1[:2]+1j*v1[2:],z1) and np.allclose(G2,-1j*Hm))
# ===== steps 8 and 11 to 15 on four real counters r = (u+, u-, v+, v-): every line against the usual complex form (psi = u + i v, JJ = i)
from scipy.linalg import expm
from scipy.special import j0, y0, k0, hankel2
def check(n,cnd):
    num,rest=n.split(' ',1); ok("%d. %s"%(int(num)+11,rest),cnd)
I4=np.eye(4)
AA=np.diag([1.,-1.,1.,-1.])
TT=np.array([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]],float)
JJ=np.array([[0,0,-1,0],[0,0,0,-1],[1,0,0,0],[0,1,0,0]],float)
check("1 AA AA = 1, TT TT = -1, JJ JJ = -1; AA TT = -TT AA; JJ commutes with AA and with TT; JJ sends (u, v) to (-v, u)",
      np.allclose(AA@AA,I4) and np.allclose(TT@TT,-I4) and np.allclose(JJ@JJ,-I4) and np.allclose(AA@TT,-TT@AA)
      and np.allclose(JJ@AA,AA@JJ) and np.allclose(JJ@TT,TT@JJ) and np.allclose(JJ@np.array([1,2,3,4.]),[-3,-4,1,2]))
# 2 a wave: r(x) = (cos kx + JJ sin kx) r_k ; d_x r = k JJ r
x,k=sp.symbols('x k',real=True); Js=sp.Matrix(JJ.astype(int)); Wv=sp.cos(k*x)*sp.eye(4)+sp.sin(k*x)*Js
check("2 a wave r(x) = (cos kx + JJ sin kx) r_k has d_x r = k JJ r", sp.simplify(sp.diff(Wv,x)-k*Js*Wv)==sp.zeros(4,4))
a,c,kk=1.3,2.0,0.9; w=np.sqrt(c*c*kk*kk+a*a)
G=-c*kk*AA@JJ+a*TT
Gpage=np.array([[0,a,c*kk,0],[-a,0,0,-c*kk],[-c*kk,0,0,a],[0,c*kk,-a,0]])
check("3 G = -ck AA JJ + a TT is the matrix written on the page; it is antisymmetric; G G = -(c^2k^2 + a^2)",
      np.allclose(G,Gpage) and np.allclose(G.T,-G) and np.allclose(G@G,-w*w*I4))
t=0.7
check("4 e^(tG) = cos(wt) + (G/w) sin(wt), all real", np.allclose(expm(t*G),np.cos(w*t)*I4+G/w*np.sin(w*t)))
JG=JJ@G; JGpage=np.array([[c*kk,0,0,-a],[0,-c*kk,a,0],[0,a,c*kk,0],[-a,0,0,-c*kk]])
Pp,Pm=(I4+JG/w)/2,(I4-JG/w)/2
check("5 JJ G = ck AA + a JJ TT is the matrix on the page; it is symmetric and its square is w^2",
      np.allclose(JG,JGpage) and np.allclose(JG,c*kk*AA+a*JJ@TT) and np.allclose(JG.T,JG) and np.allclose(JG@JG,w*w*I4))
check("6 the two halves: Pi+ + Pi- = 1, each applied twice is itself, Pi+ Pi- = 0, both symmetric, each of rank 2",
      np.allclose(Pp+Pm,I4) and np.allclose(Pp@Pp,Pp) and np.allclose(Pm@Pm,Pm) and np.allclose(Pp@Pm,0) and np.allclose(Pp,Pp.T)
      and np.linalg.matrix_rank(Pp)==2 and np.linalg.matrix_rank(Pm)==2)
check("7 on the half Pi+ the step is G = -w JJ, on the half Pi- it is G = +w JJ", np.allclose(G@Pp,-w*JJ@Pp) and np.allclose(G@Pm,w*JJ@Pm))
Rm=lambda s: np.cos(w*s)*I4-JJ*np.sin(w*s); Rp=lambda s: np.cos(w*s)*I4+JJ*np.sin(w*s)
check("8 e^(tG) = (cos wt - JJ sin wt) Pi+ + (cos wt + JJ sin wt) Pi-", np.allclose(expm(t*G),Rm(t)@Pp+Rp(t)@Pm))
# dictionary with the usual symbols: a complex 2x2 matrix X + iY acts on psi = u + iv as the real 4x4 [[X,-Y],[Y,X]]
def real_of(C): X,Y=C.real,C.imag; return np.block([[X,-Y],[Y,X]])
sz=np.diag([1.,-1.]); sy=np.array([[0,-1j],[1j,0]]); H=c*kk*sz-a*sy; E=w; pip,pim=(np.eye(2)+H/E)/2,(np.eye(2)-H/E)/2
check("9 dictionary: JJ is i; AA is sigma_z; TT is i sigma_y; G is -iH; JJ G is H (H = c sigma_z p - m c^2 sigma_y, hbar = 1)",
      np.allclose(real_of(1j*np.eye(2)),JJ) and np.allclose(real_of(sz+0j),AA) and np.allclose(real_of(1j*sy),TT)
      and np.allclose(real_of(-1j*H),G) and np.allclose(real_of(H),JG))
check("10 dictionary: e^(tG) is e^(-iHt); Pi+ and Pi- are (1 +- H/E)/2", np.allclose(real_of(expm(-1j*H*t)),expm(t*G)) and np.allclose(real_of(pip),Pp) and np.allclose(real_of(pim),Pm))
KFc=lambda s: (s>0)*np.exp(-1j*E*s)*pip-(s<0)*np.exp(1j*E*s)*pim
KFr=lambda s: (s>0)*Rm(s)@Pp-(s<0)*Rp(s)@Pm
check("11 Feynman's rule on real counters, theta(t)(cos wt - JJ sin wt)Pi+ - theta(-t)(cos wt + JJ sin wt)Pi-, is the usual one, for t = 0.7 and t = -0.4",
      np.allclose(real_of(KFc(0.7)),KFr(0.7)) and np.allclose(real_of(KFc(-0.4)),KFr(-0.4)))
Hm=c*(-kk)*sz-a*sy
check("12 on ONE real pair the halves do not separate: the conjugate of Pi+(k) is Pi-(-k)", np.allclose(np.conj(pip),(np.eye(2)-Hm/E)/2))
# closed form in x: Delta_F = -(i/4c) H0^(2)(zeta) inside, K0(zeta')/(2 pi c) outside; i Delta_F in real form
for z in (0.4,1.7,5.2):
    iD=1j*(-(1j/(4*c))*hankel2(0,z)); F=(j0(z)*I4-JJ*y0(z))/(4*c)
    assert np.allclose(real_of(iD*np.eye(2)),F)
iDo=1j*k0(0.8)/(2*np.pi*c)
check("13 i Delta_F is (1/4c)[J_0 - JJ Y_0] inside the light cone and JJ K_0/(2 pi c) outside", np.allclose(real_of(iDo*np.eye(2)),JJ*k0(0.8)/(2*np.pi*c)))
# K_F = (i d_t + H) Delta_F = (d_t + G_x)(i Delta_F), G_x = -c AA d_x + a TT : the operator identity, by finite differences at two points
def DF(tq,xq):
    s2=c*c*tq*tq-xq*xq
    return -(1j/(4*c))*hankel2(0,a/c*np.sqrt(s2)) if s2>0 else k0(a/c*np.sqrt(-s2))/(2*np.pi*c)
def Freal(tq,xq):
    s2=c*c*tq*tq-xq*xq
    if s2>0: z=a/c*np.sqrt(s2); return (j0(z)*I4-JJ*y0(z))/(4*c)
    return JJ*k0(a/c*np.sqrt(-s2))/(2*np.pi*c)
h=1e-5; good=True
for (tq,xq) in ((1.1,0.9),(0.4,1.6)):
    dtD=(DF(tq+h,xq)-DF(tq-h,xq))/(2*h); dxD=(DF(tq,xq+h)-DF(tq,xq-h))/(2*h)
    Kc=1j*dtD*np.eye(2)+(-1j*c*sz*dxD-a*sy*DF(tq,xq))
    dtF=(Freal(tq+h,xq)-Freal(tq-h,xq))/(2*h); dxF=(Freal(tq,xq+h)-Freal(tq,xq-h))/(2*h)
    Kr=dtF-c*AA@dxF+a*TT@Freal(tq,xq)
    good&=np.allclose(real_of(Kc),Kr,atol=1e-6)
check("14 (i d_t + H) Delta_F in the usual symbols is (d_t + G_x) F on real counters, inside (t = 1.1, x = 0.9) and outside (t = 0.4, x = 1.6) the light cone", good)
print("\nALL PASS" if not fails else "\nFAILED %d"%len(fails))
