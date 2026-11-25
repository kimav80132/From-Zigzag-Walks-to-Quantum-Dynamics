import numpy as np, itertools
np.random.seed(1)
# III.15a tick on a ring: shift each sense, then the coin  R = cos 1 + sin T ; counters ordered (Q+(1),Q-(1),Q+(2),Q-(2),...)
N=4; th=0.6; c,s=np.cos(th),np.sin(th); D=2*N
U=np.zeros((D,D))
for y in range(N):
    yp=(y-1)%N; ym=(y+1)%N          # forward counter came from the left, backward from the right
    U[2*y,2*yp]+=c;  U[2*y,2*ym+1]+=s      # Q+'(y)= c Q+(y-1) + s Q-(y+1)
    U[2*y+1,2*yp]+=-s; U[2*y+1,2*ym+1]+=c  # Q-'(y)=-s Q+(y-1) + c Q-(y+1)
assert np.allclose(U.T@U,np.eye(D))
ev,V=np.linalg.eig(U); eps=-np.angle(ev)          # U phi = e^{-i eps} phi ; eps = "energy"
neg=[i for i in range(D) if eps[i]<-1e-9]; pos=[i for i in range(D) if eps[i]>1e-9]
assert len(neg)==len(pos)==D//2
Phi,_=np.linalg.qr(V[:,neg]); Lam,_=np.linalg.qr(V[:,pos]); V[:,neg]=Phi; V[:,pos]=Lam
Pm=Phi@Phi.conj().T; Pp=Lam@Lam.conj().T
print('0. Pi_- is a projector, Pi_-+Pi_+=1:',np.allclose(Pm@Pm,Pm),np.allclose(Pm+Pp,np.eye(D)))
print("1. conj(Pi_-) == Pi_+ (U real: the conjugate of a negative-energy column is a positive-energy one):",np.allclose(Pm.conj(),Pp))
# --- Fock space by hand: an n-walker state = antisymmetric n-index array (dict over sorted index tuples with sign)
def slater(cols):            # cols: list of D-vectors -> antisymmetric tensor as dict {sorted tuple: amplitude}
    n=len(cols); out={}
    for idx in itertools.combinations(range(D),n):
        M=np.array([[cols[j][idx[i]] for j in range(n)] for i in range(n)])
        out[idx]=np.linalg.det(M)
    return out
def remove(w,state):         # c_w : contract the first index with conj(w)
    out={}
    for idx,amp in state.items():
        for p,X in enumerate(idx):
            rest=idx[:p]+idx[p+1:]; sign=(-1)**p
            out[rest]=out.get(rest,0)+sign*np.conj(w[X])*amp
    return out
def add(v,state):            # c_v^dagger : wedge a new column in
    out={}
    for idx,amp in state.items():
        for X in range(D):
            if X in idx: continue
            new=tuple(sorted(idx+(X,))); p=new.index(X); sign=(-1)**p
            out[new]=out.get(new,0)+sign*v[X]*amp
    return out
def evolve(state,t):
    Ut=np.linalg.matrix_power(U,t); n=len(next(iter(state)))
    out={}
    for idx,amp in state.items():
        # evolve each index: sum over images with det structure -> apply U to each column of a Slater is the same as
        # transforming the tensor by U⊗...⊗U
        for jdx in itertools.combinations(range(D),n):
            M=Ut[np.ix_(jdx,idx)]; out[jdx]=out.get(jdx,0)+np.linalg.det(M)*amp
    return out
def overlap(a,b): return sum(np.conj(a.get(k,0))*b.get(k,0) for k in set(a)|set(b))
sea=slater([Phi[:,j] for j in range(D//2)])
t=3; Ut=np.linalg.matrix_power(U,t)
vac_amp=overlap(sea,evolve(sea,t))
e=lambda X: np.eye(D)[X]
A,X=1,4
holeA=remove(e(A),sea); holeX=remove(e(X),sea)
ratio=overlap(holeX,evolve(holeA,t))/vac_amp
print("2. hole propagator  <hole X|U^t|hole A>/<sea|U^t|sea> =",np.round(ratio,6))
print("   conj((U^t Pi_-)_{XA})                               =",np.round(np.conj((Ut@Pm)[X,A]),6),"  (VI.4 (7))")
print("   (U^t Pi_+)_{XA}  = the LEADING loop, forward, A->X   =",np.round((Ut@Pp)[X,A],6))
# energy cost: phase of hole state relative to vacuum vs phase of added leading loop
lam=V[:,pos[0]]; part=add(lam,sea)
r_part=overlap(part,evolve(part,t))/vac_amp
ell=Phi[:,0]; hole=remove(ell,sea)
r_hole=overlap(hole,evolve(hole,t))/vac_amp
print("3. add leading loop lam (eps=%.4f): phase/vac = e^{-i eps t}? "%eps[pos[0]], np.allclose(r_part,np.exp(-1j*eps[pos[0]]*t)), " |r|=%.12f"%abs(r_part))
print("   remove lagging loop ell (eps=%.4f): phase/vac = e^{+i eps t} = e^{-i|eps| t}? "%eps[neg[0]], np.allclose(r_hole,np.exp(+1j*eps[neg[0]]*t)), " |r|=%.12f"%abs(r_hole))
# copy swap (conjugation) maps the lagging loop to a leading loop
print("4. conj(ell) is a positive-energy column (U conj(ell) = e^{+i eps} conj(ell)):",np.allclose(U@ell.conj(),np.exp(+1j*eps[neg[0]])*ell.conj()))
# 10 Nov 2026: the lines above print booleans; here they are tested, and the three numbers of 2 compared.
c0 = np.allclose(Pm@Pm,Pm) and np.allclose(Pm+Pp,np.eye(D)); c1 = np.allclose(Pm.conj(),Pp)
d2 = max(abs(overlap(remove(e(X_),sea),evolve(remove(e(A_),sea),t))/vac_amp - np.conj((Ut@Pm)[X_,A_])) for A_ in range(D) for X_ in range(D))
d2b = abs(Ut@Pp - np.conj(Ut@Pm)).max()
c3 = np.allclose(r_part,np.exp(-1j*eps[pos[0]]*t)) and np.allclose(r_hole,np.exp(+1j*eps[neg[0]]*t)) and abs(abs(r_part)-1)<1e-12 and abs(abs(r_hole)-1)<1e-12
c4 = np.allclose(U@ell.conj(),np.exp(+1j*eps[neg[0]])*ell.conj())
ok_all = c0 and c1 and d2<1e-12 and d2b<1e-12 and c3 and c4
print("%s  0-4 tested: projectors %s; conj(Pi_-) = Pi_+ %s; hole amplitude = conj((U^t Pi_-)_XA) for all %d (A,X) to %.0e, and that = (U^t Pi_+)_XA to %.0e; "
      "phases of an added leading loop and of a removed lagging loop %s; conj(ell) a positive column %s"
      % ("PASS" if ok_all else "FAIL", c0, c1, D*D, d2, d2b, c3, c4))
raise SystemExit(0 if ok_all else 1)
