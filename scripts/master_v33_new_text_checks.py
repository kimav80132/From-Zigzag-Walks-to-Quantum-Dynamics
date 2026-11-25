"""master_v33_new_text_checks.py -- 9 Nov 2026.
Checks of the sentences that were written new when the pending drafts were put into the master (v33), i.e. the
statements of III.15d, III.15e, VI.5, VI.6 that are not already covered, in these words, by the drafts' own scripts.
Exact linear algebra on small matrices; nothing is solved by a method or approximated.

 1  III.15e (C): the leading and lagging loops written out in time: twin = route a quarter period on / back.
 2  III.15e (C): the two classes, T r = -/+ J r; P_+ is the leading class; a bare route is half of each.
 3  III.15e (C): the launch rule v(0)=Tu(0) leaves (omega-a)/2omega = (E-mc^2)/2E in the other class: 0, 0.090, 0.276.
 4  III.15d (10): the decaying eigenvector of M_k is the rest solution (1,-1)e^{-ax/c} seen from the frame -eta.
 5  III.15d (12): -omega + B M_k B^-1 = [[-a e^eta, a e^-eta],[a e^eta, -a e^-eta]], columns summing to zero.
 6  VI.5 (3): log(1+e^{a tau}) = a tau + sum (-1)^(n+1) e^{-n a tau}/n ; the shares 0.64, 0.993, 0.999998.
 7  VI.6 (7): e^{tG} = e^{-i w t} Pi_+ + e^{+i w t} Pi_- ; both K and K_F jump by 1 at t=0; the complex Pi_+- are P_+- with J -> i.
 8  VI.6 (F): the charge at a source: both legs leaving gives q(w_+ - w_-) != q; one leaving, one arriving gives q.
"""
import numpy as np
from scipy.linalg import expm
ok=True
def report(name,cond,extra=""):
    global ok; ok&=bool(cond); print(f"{name}: {'PASS' if cond else 'FAIL'}  {extra}")
a=1.0; c=1.0
A=np.diag([1.,-1,1,-1])
S=np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0.]])
T4=np.array([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0.]])
J4=np.array([[0,0,-1,0],[0,0,0,-1],[1,0,0,0],[0,1,0,0.]])
lam=np.array([1.,0,0,-1]); ell=np.array([1.,0,0,1])
I4=np.eye(4)
# 1
good=True
for t in np.linspace(0,7,29):
    rl=expm(a*T4*t)@lam; re=expm(a*T4*t)@ell
    u=lambda s: np.array([np.cos(a*s),-np.sin(a*s)])
    good&=np.allclose(rl[:2],u(t)) and np.allclose(rl[2:],u(t+np.pi/(2*a))) and np.allclose(re[:2],u(t)) and np.allclose(re[2:],u(t-np.pi/(2*a)))
report("1. leading loop: v(t) = u(t + pi/2a); lagging loop: v(t) = u(t - pi/2a); u(t) = (cos at, -sin at) for both",good)
# 2
G=lambda k: -c*k*A@J4+a*T4
P=lambda k,s: 0.5*(I4+s*J4@G(k)/np.sqrt(a*a+c*c*k*k))
good=np.allclose(T4@lam,-J4@lam) and np.allclose(T4@ell,J4@ell) and np.allclose(P(0,+1)@lam,lam) and np.allclose(P(0,-1)@ell,ell) \
     and np.allclose(np.array([1.,0,0,0]),0.5*lam+0.5*ell) and np.allclose(np.diag([1.,1,-1,-1])@lam,ell) \
     and np.allclose(expm(a*T4*0.37)@lam,expm(-a*J4*0.37)@lam) and np.allclose(expm(a*T4*0.37)@ell,expm(+a*J4*0.37)@ell)
for k in (0,0.7,2.0):
    for s in (+1,-1): good&=np.allclose(P(k,s)@P(k,s),P(k,s)) and np.allclose(P(k,+1)@P(k,-1),0)
report("2. T lam = -J lam, T ell = +J ell; lam turns as e^{-atJ}, ell as e^{+atJ}; P_+ is the leading class; (1,0,0,0) = lam/2 + ell/2; copy swap maps lam to ell",good)
# 3
ws=[]
for k in (0,0.7,2.0):
    w=np.sqrt(a*a+c*c*k*k); r=lam/np.linalg.norm(lam); ws.append(np.linalg.norm(P(k,-1)@r)**2)
    ok_k=np.isclose(ws[-1],(w-a)/(2*w))
    good&=ok_k
report("3. launch rule's weight in the other class = (omega-a)/2omega",good,"("+", ".join(f"{x:.3f}" for x in ws)+" at ck/a = 0, 0.7, 2)")
# 4, 5
good4=True; good5=True
for eta in (0.3,1.2,-0.7):
    ck=a*np.sinh(eta); w=a*np.cosh(eta); M=np.array([[-ck,a],[a,ck]])
    vdec=np.array([1,-np.exp(-eta)]); good4&=np.allclose(M@vdec,-w*vdec)
    Bm=np.diag([np.exp(+eta/2),np.exp(-eta/2)])          # boost by -eta on the readings
    col=Bm@np.array([1,-1.]); good4&=np.allclose(col/col[0],vdec)
    # exponent: rest solution e^{-a x/c}, boost by -eta: x = gamma(x' - v t') -> e^{-w x'} e^{+ck t'} : slope +ck, as M_k's solutions have
    g=np.cosh(eta); v=np.tanh(eta); good4&=np.isclose(-a*g,-w) and np.isclose(+a*g*v,ck)
    B=np.diag([np.exp(-eta/2),np.exp(eta/2)]); Gen=-w*np.eye(2)+B@M@np.linalg.inv(B)
    good5&=np.allclose(Gen,[[-a*np.exp(eta),a*np.exp(-eta)],[a*np.exp(eta),-a*np.exp(-eta)]]) and np.allclose(Gen.sum(axis=0),0)
report("4. decaying eigenvector of M_k is (1, -e^{-eta}); it is (1,-1)e^{-ax/c} boosted by -eta, slope e^{+ck t'}",good4)
report("5. -omega + B M_k B^-1 is the two-rate generator, columns sum to zero",good5)
# 6
x=1.4; good=np.isclose(np.log(1+np.exp(x)),x+sum((-1)**(n+1)*np.exp(-n*x)/n for n in range(1,200)))
sh=[np.exp(y)/(2+2*np.cosh(y)) for y in (1.4,5.6,14.0)]
report("6. log(1+e^{a tau}) = a tau + sum (-1)^(n+1) e^{-n a tau}/n ; share of the growing configuration",good and np.allclose(sh,[0.6435,0.99265,0.999998],atol=6e-5),"("+", ".join(f"{s:.6f}" for s in sh)+" at a tau = 1.4, 5.6, 14)")
# 7
good=True
for k in (0,0.7,2.0):
    w=np.sqrt(a*a+c*c*k*k); Gc=np.array([[-1j*c*k,a],[-a,1j*c*k]]); Pp=0.5*(np.eye(2)+1j*Gc/w); Pm=0.5*(np.eye(2)-1j*Gc/w)
    t=0.83; good&=np.allclose(expm(t*Gc),np.exp(-1j*w*t)*Pp+np.exp(1j*w*t)*Pm) and np.allclose(Pp+Pm,np.eye(2)) and np.allclose(Pp@Pm,0)
    # jumps: K(0+)-K(0-) = Pp+Pm ; K_F(0+)-K_F(0-) = Pp-(-Pm)
    good&=np.allclose(Pp-(-Pm),np.eye(2))
    # complex form of the real 4x4 projector: (u+ + i v+, u- + i v-)
    def cplx(Mat):   # 4x4 real commuting with J4 -> 2x2 complex
        return Mat[:2,:2]+1j*Mat[2:,:2]
    good&=np.allclose(cplx(P(k,+1)),Pp) and np.allclose(cplx(G(k)),Gc)
report("7. e^{tG} = e^{-iwt}Pi_+ + e^{+iwt}Pi_-; K and K_F both jump by 1 at the source; Pi_+- are P_+- with J written as i",good)
# 8
good=True; shown=[]
for k in (0.7,2.0):
    w=np.sqrt(a*a+c*c*k*k); Gc=np.array([[-1j*c*k,a],[-a,1j*c*k]]); Pp=0.5*(np.eye(2)+1j*Gc/w); Pm=0.5*(np.eye(2)-1j*Gc/w)
    psi=np.array([1,0]); wp=np.linalg.norm(Pp@psi)**2; wm=np.linalg.norm(Pm@psi)**2
    good&=np.isclose(wp+wm,1) and abs((wp-wm)-1)>1e-2; shown.append(f"k={k}: w+={wp:.3f}, w-={wm:.3f}")
report("8. one walker injected: both legs leaving carries q(w+ - w-) != q; leading leg leaving and the other arriving carries q(w+ + w-) = q",good,"("+"; ".join(shown)+")")
print("ALL PASS" if ok else "SOME FAIL")
