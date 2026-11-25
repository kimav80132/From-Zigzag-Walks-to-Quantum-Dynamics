"""zigzag_signs_and_twin_checks.py  --  recreated 9 Nov 2026 (the original is not on this disk).
Check numbers are the ones the master cites: 1, 4, 5a, 5b, 8, 9.  (The original's 2, 3, 6, 7 are cited nowhere
and are not recreated.)  Exact enumeration and symbolic algebra only.

 1   III.8c (4)-(6): each outward route with its twin (same turns, waiting times swapped pair by pair):
     they meet after every completed pair; the twin's weight is the outward weight with the launch sign;
     the signed sum of all routes launched backward is U^N T e_+.
 4   III.15a (A): the sign rule  w = (-1)^{n(+ -> -)} sin^m cos^(N-m), and the three-tick example.
 5a  III.8c (7)-(8): D rho = +-2aK, rho D theta = aW, symbolically, for any two columns.
 5b  III.8c (9)-(10): K = sqrt(rho+ rho-) sin phi, W = -sqrt(rho+ rho-) cos phi; one column alone has W = 0.
 8   the cut: in the walk's own clock the same start gives K = u+u- + v+v- >= 0 everywhere; on time slices both signs.
 9   III.15a: the event weight is the sum over ALL (out, back) pairs = Q^2; the cross pairs carry the interference.
"""
import itertools, numpy as np, sympy as sp
ok=True
def report(name,cond,extra=""):
    global ok; ok&=bool(cond); print(f"{name}: {'PASS' if cond else 'FAIL'}  {extra}")

th=0.37; c,s=np.cos(th),np.sin(th)
T=np.array([[0,1],[-1,0]])
def routes(N,launch):
    """all N-tick routes: at each tick move one site in the current sense, then keep (0) or flip (1)."""
    for flips in itertools.product((0,1),repeat=N):
        sense=launch; x=0; w=1.0; npm=0; pos=[]
        for f in flips:
            x+=sense; pos.append(x)
            if f:
                if sense==+1: w*=-s; npm+=1
                else: w*=+s
                sense=-sense
            else: w*=c
        yield flips,x,sense,w,npm,tuple(pos)
def U_column(N,launch_col,L=64):
    """the same walk as a matrix product on a ring big enough not to wrap: shift each sense, then the coin."""
    Q=np.zeros((L,2)); Q[0]=launch_col
    for _ in range(N):
        Qs=np.zeros_like(Q); Qs[:,0]=np.roll(Q[:,0],1); Qs[:,1]=np.roll(Q[:,1],-1)
        Q=np.stack([c*Qs[:,0]+s*Qs[:,1], -s*Qs[:,0]+c*Qs[:,1]],axis=1)
    return Q

# ---------------- 1. the twin
def stretches(flips,N):
    """stretch lengths of a route whose LAST tick is a flip-free end: lengths of constant-sense runs."""
    out=[]; run=0
    for f in flips:
        run+=1
        if f: out.append(run); run=0
    if run: out.append(run)
    return out
N=10; meet=True; wrel=True; nchecked=0; twins=set()
# (10 Nov 2026) the twin's own weight is now taken from the enumeration of the backward-launched routes, not from the sign-rule formula
back_w={}
for flips_b,x_b,sense_b,w_b,npm_b,pos_b in routes(N,-1):
    if not flips_b[-1]: back_w[tuple(stretches(flips_b,N))]=w_b
wenum=True
for flips,x,sense,w,npm,pos in routes(N,+1):
    if flips[-1]: continue                      # route ends inside a stretch
    d=stretches(flips,N)
    if len(d)%2: continue                       # completed pairs only
    # twin: launched backward, stretches swapped pair by pair
    dt=[]; 
    for i in range(0,len(d),2): dt+=[d[i+1],d[i]]
    # positions after each completed pair
    xo=np.cumsum([(+1 if i%2==0 else -1)*d[i] for i in range(len(d))])[1::2]
    xt=np.cumsum([(-1 if i%2==0 else +1)*dt[i] for i in range(len(dt))])[1::2]
    meet&=np.array_equal(xo,xt)
    # twin's own weight in the walk, launched in the backward sense
    m=len(d)-1; n_pm_twin=(len(d)//2)-1
    wt=(-1)**n_pm_twin*s**m*c**(N-m)
    wrel&=np.isclose(w,-wt)                     # "counted with the opposite sign": v(0) = -e_- g
    wenum&=(tuple(dt) in back_w) and np.isclose(back_w[tuple(dt)],wt)
    twins.add(tuple(dt)); nchecked+=1
report("1.  twin meets the outward route after every completed pair; weight = -(twin's own weight, read off the enumeration of backward launches)",meet and wrel and wenum,f"({nchecked} routes, N={N}; the twins are {len(twins)} distinct backward-launched routes)")
Qsum=np.zeros((64,2))
for flips,x,sense,w,npm,pos in routes(N,-1): Qsum[x%64,0 if sense==+1 else 1]+=-w     # signed: launch column is -e_-
report("1'. signed sum of all backward-launched routes = U^N T e_+",np.allclose(Qsum,U_column(N,T@np.array([1,0]))),"(III.8c (6))")

# ---------------- 4. the sign rule and the three-tick example
good=True
for N4 in (3,6,9):
    Q=np.zeros((64,2))
    for flips,x,sense,w,npm,pos in routes(N4,+1):
        m=sum(flips); good&=np.isclose(w,(-1)**npm*s**m*c**(N4-m)); Q[x%64,0 if sense==+1 else 1]+=w
    good&=np.allclose(Q,U_column(N4,np.array([1,0])))
report("4.  sign rule w=(-1)^{n(+->-)} sin^m cos^(N-m), and the path sum is the matrix product",good,"(N=3,6,9)")
ex=[(f,w,pos) for f,x,sense,w,npm,pos in routes(3,+1) if x==1 and sense==-1]
want={(0,1,0):(-c*c*s,(1,2,1)),(1,1,1):(s**3,(1,0,1))}
good=len(ex)==2 and all(np.isclose(w,want[f][0]) and pos==want[f][1] for f,w,pos in ex)
c4,s4=np.cos(np.pi/4),np.sin(np.pi/4)
report("4'. three ticks, end at x=1 in sense -: (keep,flip,keep) = -cos^2 sin through (1,2,1); (flip,flip,flip) = +sin^3 through (1,0,1); sum = 0 at theta=pi/4",good and abs(-c4*c4*s4+s4**3)<1e-15)

# ---------------- 5a. the fluid, symbolically
t,x,a,cc=sp.symbols('t x a c',real=True)
up,um,vp,vm=[sp.Function(n)(t,x) for n in ('u_p','u_m','v_p','v_m')]
Dp=lambda f: sp.diff(f,t)+cc*sp.diff(f,x); Dm=lambda f: sp.diff(f,t)-cc*sp.diff(f,x)
# the lines (7): D+ r+ = a r-,  D- r- = -a r+   ->  solve for the time derivatives and substitute
rules={sp.Derivative(up,t):a*um-cc*sp.diff(up,x), sp.Derivative(vp,t):a*vm-cc*sp.diff(vp,x),
       sp.Derivative(um,t):-a*up+cc*sp.diff(um,x), sp.Derivative(vm,t):-a*vp+cc*sp.diff(vm,x)}
K=up*um+vp*vm; W=up*vm-vp*um; rp=up**2+vp**2; rm=um**2+vm**2
e1=sp.simplify(Dp(rp).subs(rules)-2*a*K); e2=sp.simplify(Dm(rm).subs(rules)+2*a*K)
# rho D theta = u Dv - v Du
e3=sp.simplify((up*Dp(vp)-vp*Dp(up)).subs(rules)-a*W); e4=sp.simplify((um*Dm(vm)-vm*Dm(um)).subs(rules)-a*W)
report("5a. D+rho+ = 2aK, D-rho- = -2aK, rho+ D+theta+ = aW, rho- D-theta- = aW for any two columns",all(e==0 for e in (e1,e2,e3,e4)),"(symbolic)")

# ---------------- 5b. the angles
rng=np.random.default_rng(0); good=True
for _ in range(1000):
    Up,Um,Vp,Vm=rng.normal(size=4)
    rp_,rm_=Up*Up+Vp*Vp,Um*Um+Vm*Vm
    phi=(np.arctan2(Vm,Um)+np.pi/2)-np.arctan2(Vp,Up)          # (S_- - S_+)/hbar with the pi/2 of (9)
    good&=np.isclose(Up*Um+Vp*Vm,np.sqrt(rp_*rm_)*np.sin(phi)) and np.isclose(Up*Vm-Vp*Um,-np.sqrt(rp_*rm_)*np.cos(phi))
report("5b. K = sqrt(rho+rho-) sin(phi), W = -sqrt(rho+rho-) cos(phi) with S_- = hbar(theta_- + pi/2); one column alone (v=0) has W = 0",good and (W.subs({vp:0,vm:0})==0),"(1000 random pairs)")

# ---------------- 8. the cut
N8=40; L=128; g=np.exp(-0.5*((np.arange(L)-L//2)/4.0)**2)
# time slices: u launched forward with density g, v = U T u(0)
def evolve(Q,N,coin):
    for _ in range(N):
        Qs=np.zeros_like(Q); Qs[:,0]=np.roll(Q[:,0],1); Qs[:,1]=np.roll(Q[:,1],-1); Q=Qs@coin.T
    return Q
Rt=np.array([[c,s],[-s,c]]); u0=np.zeros((L,2)); u0[:,0]=g; v0=u0@T.T
u=evolve(u0.copy(),N8,Rt); v=evolve(v0.copy(),N8,Rt); Kt=u[:,0]*u[:,1]+v[:,0]*v[:,1]
# the walk's own clock: the index now runs along t, the coin is the exchange (all weights positive), the twins count positively
ch,sh=np.cosh(th),np.sinh(th); Rx=np.array([[ch,sh],[sh,ch]]); Sx=np.array([[0,1],[1,0]])
ux=evolve(u0.copy(),N8,Rx); vx=evolve(u0@Sx.T,N8,Rx); Kx=ux[:,0]*ux[:,1]+vx[:,0]*vx[:,1]
report("8.  walk's own clock: every crossing count >= 0 and K = u+u- + v+v- >= 0 everywhere; on time slices the same start gives both signs",
       ux.min()>=0 and vx.min()>=0 and Kx.min()>=0 and Kt.min()<-1e-6 and Kt.max()>1e-6, f"(clock: min K = {Kx.min():.1e}; time slice: K from {Kt.min():.3f} to {Kt.max():.3f})")

# ---------------- 9. the cross pairs
N9=8; by_end={}
for flips,x,sense,w,npm,pos in routes(N9,+1): by_end.setdefault((x,sense),[]).append(w)
good=True; biggest=0
for key,ws in by_end.items():
    ws=np.array(ws); allpairs=np.sum(np.outer(ws,ws)); diag=np.sum(ws**2); Q=ws.sum()
    good&=np.isclose(allpairs,Q*Q); biggest=max(biggest,abs(allpairs-diag))
w1,w2=-c*c*s,s**3
report("9.  event weight = sum over all (out, back) pairs = Q^2 at every endpoint; it is not the sum of routes paired with themselves",good and biggest>1e-3,f"(N={N9}; largest cross-pair total {biggest:.3f})")
report("9'. the example's cross pair: out by (keep,flip,keep), home by the reverse of (flip,flip,flip) = -cos^2 sin^4, either way round",np.isclose(w1*w2,-c**2*s**4) and np.isclose(w2*w1,-c**2*s**4))
print("ALL PASS" if ok else "SOME FAIL"); raise SystemExit(0 if ok else 1)
