"""Many-body form of 'a hole in the field theta is a walker above the sea in the field -theta'.
Brute-force Fock space (plain lists of configuration amplitudes, the determinant only in the one-tick step of each
configuration), lattice tick of III.15a on a ring of 4 sites, field = a different turn e^{i theta_x} of the pair on
each site at each tick (V.2).  9 Nov 2026.
Checks:  1 one hole  2 two holes  3 a pair made from the vacuum  4 the contrast (same-sign field) fails."""
import numpy as np, itertools
rng=np.random.default_rng(3)
N=4; th=0.6; c,s=np.cos(th),np.sin(th); D=2*N
U=np.zeros((D,D))
for y in range(N):
    yp=(y-1)%N; ym=(y+1)%N
    U[2*y,2*yp]+=c;  U[2*y,2*ym+1]+=s
    U[2*y+1,2*yp]+=-s; U[2*y+1,2*ym+1]+=c
ev,V=np.linalg.eig(U); eps=-np.angle(ev)
neg=[i for i in range(D) if eps[i]<-1e-9]; pos=[i for i in range(D) if eps[i]>1e-9]
Phi,_=np.linalg.qr(V[:,neg])
T=3; angles=rng.uniform(-0.5,0.5,size=(T,N))          # a field: site- and tick-dependent turn of the pair
def W(sign):
    M=np.eye(D,dtype=complex)
    for t in range(T):
        ph=np.repeat(np.exp(1j*sign*angles[t]),2)      # same turn on both senses of a site: a turn of the pair (J), not of the coin
        M=np.diag(ph)@U@M
    return M
def slater(cols):
    n=len(cols); out={}
    for idx in itertools.combinations(range(D),n):
        out[idx]=np.linalg.det(np.array([[cols[j][idx[i]] for j in range(n)] for i in range(n)]))
    return out
def remove(w,state):
    out={}
    for idx,amp in state.items():
        for p,X in enumerate(idx):
            rest=idx[:p]+idx[p+1:]; out[rest]=out.get(rest,0)+(-1)**p*np.conj(w[X])*amp
    return out
def add(v,state):
    out={}
    for idx,amp in state.items():
        for X in range(D):
            if X in idx: continue
            new=tuple(sorted(idx+(X,))); p=new.index(X); out[new]=out.get(new,0)+(-1)**p*v[X]*amp
    return out
def evolve(state,M):
    n=len(next(iter(state))); out={}
    for idx,amp in state.items():
        for jdx in itertools.combinations(range(D),n):
            out[jdx]=out.get(jdx,0)+np.linalg.det(M[np.ix_(jdx,idx)])*amp
    return out
def overlap(a,b): return sum(np.conj(a.get(k,0))*b.get(k,0) for k in set(a)|set(b))
e=lambda X: np.eye(D)[X]
sea=slater([Phi[:,j] for j in range(D//2)])
Wp,Wm=W(+1),W(-1)
vp=overlap(sea,evolve(sea,Wp)); vm=overlap(sea,evolve(sea,Wm))
print("0. vacuum factor: |<sea|W_theta|sea>| = %.6f, |<sea|W_-theta|sea>| = %.6f, equal moduli: %s"%(abs(vp),abs(vm),np.isclose(abs(vp),abs(vm))))
ok=True
# 1. one hole in theta  vs  one walker above the sea in -theta, every A, X
d1=0
for A in range(D):
    for X in range(D):
        h=overlap(remove(e(X),sea),evolve(remove(e(A),sea),Wp))/vp
        p=overlap(add(e(X),sea),evolve(add(e(A),sea),Wm))/vm
        d1=max(d1,abs(h-p))
print("1. one hole A->X in theta = one walker A->X in -theta, all 64 (A,X): max diff %.1e  %s"%(d1,"PASS" if d1<1e-12 else "FAIL")); ok&=d1<1e-12
# 2. two holes vs two walkers
d2=0; big2=0; pairs=[((0,3),(5,6)),((1,2),(4,7)),((2,5),(0,1)),((3,6),(2,7))]
for (A,B),(X,Y) in pairs:
    h=overlap(remove(e(X),remove(e(Y),sea)),evolve(remove(e(A),remove(e(B),sea)),Wp))/vp
    p=overlap(add(e(X),add(e(Y),sea)),evolve(add(e(A),add(e(B),sea)),Wm))/vm
    d2=max(d2,abs(h-p)); big2=max(big2,abs(h))
print("2. two holes (A,B)->(X,Y) in theta = two walkers in -theta: max diff %.1e (largest amplitude %.3f)  %s"%(d2,big2,"PASS" if d2<1e-12 else "FAIL")); ok&=d2<1e-12
# 3. a pair from the vacuum: walker at X and hole at Y in theta  vs  hole at X and walker at Y in -theta
d3=0; big=0
for X in range(D):
    for Y in range(D):
        a=overlap(add(e(X),remove(e(Y),sea)),evolve(sea,Wp))/vp
        b=overlap(remove(e(X),add(e(Y),sea)),evolve(sea,Wm))/vm
        d3=max(d3,abs(a-b)); big=max(big,abs(a))
print("3. pair from the vacuum (walker X, hole Y) in theta = (hole X, walker Y) in -theta, all 64: max diff %.1e (largest amplitude %.3f)  %s"%(d3,big,"PASS" if d3<1e-12 else "FAIL")); ok&=d3<1e-12
# 4. contrast: same-sign field does not match
d4=0
for A in range(D):
    for X in range(D):
        h=overlap(remove(e(X),sea),evolve(remove(e(A),sea),Wp))/vp
        p=overlap(add(e(X),sea),evolve(add(e(A),sea),Wp))/vp
        d4=max(d4,abs(h-p))
print("4. contrast: one hole in theta vs one walker in the SAME field differ: max diff %.2f  %s"%(d4,"PASS" if d4>1e-2 else "FAIL")); ok&=d4>1e-2
print("ALL PASS" if ok else "SOME FAIL")
raise SystemExit(0 if ok else 1)
