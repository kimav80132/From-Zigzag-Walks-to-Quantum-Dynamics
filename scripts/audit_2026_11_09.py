import numpy as np, itertools
from scipy.linalg import expm
np.set_printoptions(precision=6,suppress=True)
ok=lambda b: "PASS" if b else "FAIL"
I2=np.eye(2); sx=np.array([[0.,1.],[1.,0.]]); T=np.array([[0.,1.],[-1.,0.]]); Az=np.diag([1.,-1.]); J=np.array([[0.,-1.],[1.,0.]])
a=1.0;c=1.0
print("=== A. the dictionary at k != 0, real 4x4 (order u+,v+,u-,v-) ===")
for k in [0.0,0.7,2.0]:
    w=np.sqrt(a*a+c*c*k*k)
    G4=np.block([[-c*k*J,a*I2],[-a*I2,c*k*J]]); M4=np.kron(np.array([[-c*k,a],[a,c*k]]),I2); J4=np.kron(I2,J)
    Dm=np.block([[I2,0*I2],[0*I2,J]]); Di=np.linalg.inv(Dm)
    c1=np.allclose(Di@G4@Dm,J4@M4)                         # G -> J M_k under r~ = D^-1 r
    c2=np.allclose(expm(-0.6*J4@G4),Dm@expm(0.6*M4)@Di)     # "t -> -J tau"
    # III.9's own relabelling goes the other way: r~ = D r gives -J h_std, h_std = ck sigma_z + a sigma_x
    c3=np.allclose(Dm@G4@Di,-J4@np.kron(np.array([[c*k,a],[a,-c*k]]),I2))
    # growing eigvec of M_k (route only) -> D r : is it in the E<0 class  G r = +w J r ?
    ev,V=np.linalg.eig(np.array([[-c*k,a],[a,c*k]])); g=V[:,np.argmax(ev)]; d=V[:,np.argmin(ev)]
    rg=Dm@np.array([g[0],0,g[1],0]); rd=Dm@np.array([d[0],0,d[1],0])
    c4=np.allclose(G4@rg,+w*J4@rg) and np.allclose(G4@rd,-w*J4@rd)
    print(f"k={k}: G->J M_k {ok(c1)} | e^(-tau J G)=D e^(tau M) D^-1 {ok(c2)} | III.9 direction gives -J h_std {ok(c3)} | growing->E<0, decaying->E>0 {ok(c4)}")
print()
print("=== B. one real column: the non-local quarter turn J(k)=G(k)/omega exists (III.8c (15)) ===")
k=0.7; w=np.sqrt(1+k*k)
# real form of G(k) on one column needs cos/sin in x: 4 real numbers (Q+c,Q+s,Q-c,Q-s); d/dx on (cos,sin) pair = k*J
Gr=np.block([[-c*k*J,a*I2],[-a*I2,c*k*J]]); Jk=Gr/w
print("J(k) real, commutes with the step, squares to -1:",ok(np.allclose(Jk@Gr,Gr@Jk) and np.allclose(Jk@Jk,-np.eye(4))))
print("local (k-independent) real matrices commuting with A_z and T on one column: only multiples of 1:",
      ok(np.linalg.matrix_rank(np.vstack([np.kron(np.eye(2),X)-np.kron(X.T,np.eye(2)) for X in (Az,T)]))==3))
# twin rule v(0)=T u(0) vs exact partner v=J(k)u : equal only at k=0
u=np.array([1.,0.])
for kk in [0.0,0.7,2.0]:
    ww=np.sqrt(1+kk*kk); Jc=np.array([[-1j*kk,1],[-1,1j*kk]])/ww
    print(f"   k={kk}: |T u - J(k) u| = {np.linalg.norm(T@u-Jc@u):.4f}   (twin rule = exact quadrature partner only at rest; weight in wrong half {(ww-1)/(2*ww):.4f} = (w-a)/2w)")
print()
print("=== C. pair-creation matrix element: which projector order? ===")
k,kp=0.9,-0.4
def P(kk,sgn):
    ww=np.sqrt(1+kk*kk); G4=np.block([[-kk*J,I2],[-I2,kk*J]]); return 0.5*(np.eye(4)+sgn*np.kron(I2,J)@G4/ww)
J4=np.kron(I2,J)
print("P_+(k) J P_-(k') (lift a lagging loop to a leading one: creation from the vacuum) norm:",round(np.linalg.norm(P(k,+1)@J4@P(kp,-1)),4))
print("P_-(k) J P_+(k') (drop a leading loop into a lagging one)                          norm:",round(np.linalg.norm(P(k,-1)@J4@P(kp,+1)),4))
print("both nonzero, transposes of each other up to sign:",ok(np.allclose(P(k,+1)@J4@P(kp,-1),-(P(kp,-1)@J4@P(k,+1)).T)))
print("P_+ projects onto G r = -w J r (twin leads, E>0):",ok(np.allclose(np.block([[-k*J,I2],[-I2,k*J]])@P(k,+1),-np.sqrt(1+k*k)*J4@P(k,+1))))
print()
print("=== D. Dirac pair vs lone self-conjugate column: the count's normalization ===")
# rest mode, generator a*T. Dirac (complex column): H = c^dag (i a T) c on C^2 -> fill the negative eigenvalue
E_dirac=sum(e for e in np.linalg.eigvalsh(1j*a*T) if e<0)
# Majorana (one real column, no twin): two Majorana operators g1,g2 (Pauli X, Y), H = (i/4) sum g_a A_ab g_b, A = a T
X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]]); g=[X,Y]; A=a*T
H=sum(0.25j*A[i,j]*g[i]@g[j] for i in range(2) for j in range(2))
# check Heisenberg: d g/dt = i[H,g] reproduces  dg/dt = A g  (the walk's equation on the real column)
heis=all(np.allclose(1j*(H@g[i]-g[i]@H),sum(A[i,j]*g[j] for j in range(2))) for i in range(2))
E_maj=np.linalg.eigvalsh(H).min()
print(f"Dirac walker (route+twin): E0 = {E_dirac:.4f} = -omega ;  lone real column: E0 = {E_maj:.4f} = -omega/2 ; Heisenberg eq reproduces the walk's lines: {ok(heis)}")
print("count: det(1+e^{tau M}) per k has 4 configurations (0/1 in two slots) -> the Dirac count; -(1/tau)log Z -> -omega")
print()
print("=== E. the hole in a field: propagates as a particle of the OPPOSITE charge ===")
N=4; th=0.6; cc,ss=np.cos(th),np.sin(th); D=2*N
U=np.zeros((D,D))
for y in range(N):
    U[2*y,2*((y-1)%N)]+=cc;  U[2*y,2*((y+1)%N)+1]+=ss
    U[2*y+1,2*((y-1)%N)]+=-ss; U[2*y+1,2*((y+1)%N)+1]+=cc
ev,V=np.linalg.eig(U); eps=-np.angle(ev)
Phi,_=np.linalg.qr(V[:,eps<-1e-9]); Lam,_=np.linalg.qr(V[:,eps>1e-9])
theta=np.array([0.3,-0.2,0.5,0.1])                      # a field: V.2's turn, a phase per site per tick
def Wt(sign,t):
    ph=np.diag(np.exp(1j*sign*np.repeat(theta,2))); return np.linalg.matrix_power(ph@U,t)
t=3
Wp=Wt(+1,t); Wm=Wt(-1,t)
hole=Phi@np.linalg.inv(Phi.conj().T@Wp@Phi)@Phi.conj().T          # <hole X|W|hole A>/<sea|W|sea> = hole[A,X]
part_m=Wm-Wm@Phi@np.linalg.inv(Phi.conj().T@Wm@Phi)@Phi.conj().T@Wm   # particle above the sea, field -theta: VI.4 (11)
part_p=Wp-Wp@Phi@np.linalg.inv(Phi.conj().T@Wp@Phi)@Phi.conj().T@Wp   # same charge, for contrast
print("hole(theta)[A,X] == particle(-theta)[X,A]  (opposite charge):",ok(np.allclose(hole,part_m.T)),"  max diff %.1e"%np.abs(hole-part_m.T).max())
print("hole(theta)[A,X] differs from particle(+theta)[X,A] (same charge), as it should:",ok(not np.allclose(hole,part_p.T)),"  max diff %.1e"%np.abs(hole-part_p.T).max())
# brute-force Fock check of the determinant formula used for 'hole'
def slater(cols):
    n=len(cols); return {idx:np.linalg.det(np.array([[cols[j][idx[i]] for j in range(n)] for i in range(n)])) for idx in itertools.combinations(range(D),n)}
def remove(wv,st):
    out={}
    for idx,amp in st.items():
        for p,Xx in enumerate(idx):
            r=idx[:p]+idx[p+1:]; out[r]=out.get(r,0)+(-1)**p*np.conj(wv[Xx])*amp
    return out
def evolve(st,W):
    n=len(next(iter(st))); out={}
    for idx,amp in st.items():
        for jdx in itertools.combinations(range(D),n): out[jdx]=out.get(jdx,0)+np.linalg.det(W[np.ix_(jdx,idx)])*amp
    return out
ov=lambda p,q: sum(np.conj(p.get(kq,0))*q.get(kq,0) for kq in set(p)|set(q))
sea=slater([Phi[:,j] for j in range(N)]); A_,X_=1,4; e=np.eye(D)
brute=ov(remove(e[X_],sea),evolve(remove(e[A_],sea),Wp))/ov(sea,evolve(sea,Wp))
print("brute-force Fock hole amplitude vs determinant formula:",ok(np.isclose(brute,hole[A_,X_])),brute.round(6),hole[A_,X_].round(6))
print("vacuum persistence |<sea|W|sea>| =",abs(ov(sea,evolve(sea,Wp))).round(6),"(<1: the field makes pairs; the hole ratio is still exact)")
print()
print("=== F. the same in the order (u+,u-,v+,v-) used by the hole draft ===")
A4=np.diag([1.,-1.,1.,-1.]); S4=np.kron(I2,sx); T4=np.kron(I2,T); J4o=np.kron(J,I2)
Dd=np.array([[1.,0,0,0],[0,0,0,-1.],[0,0,1.,0],[0,1.,0,0]])      # (u-,v-) -> J(u-,v-) = (-v-,u-)
print("A S = T, S^2=+1, T^2=-1, J^2=-1:",ok(np.allclose(A4@S4,T4) and np.allclose(S4@S4,np.eye(4)) and np.allclose(T4@T4,-np.eye(4)) and np.allclose(J4o@J4o,-np.eye(4))))
print("J commutes with A,S,T; T anticommutes with A:",ok(all(np.allclose(J4o@X,X@J4o) for X in (A4,S4,T4)) and np.allclose(T4@A4,-A4@T4)))
for k in [0.0,0.7,2.0]:
    G=-c*k*A4@J4o+a*T4; M=-c*k*A4+a*S4
    print(f"k={k}: D^-1 G(k) D = J M_k:",ok(np.allclose(np.linalg.inv(Dd)@G@Dd,J4o@M)))
lead=np.array([1.,0,0,-1.]); lag=np.array([1.,0,0,1.])
print("leading (1,0,0,-1): T r = -J r:",ok(np.allclose(T4@lead,-J4o@lead)),"| lagging (1,0,0,1): T r = +J r:",ok(np.allclose(T4@lag,J4o@lag)))
print("D(1,1,0,0) = lagging:",ok(np.allclose(Dd@np.array([1.,1,0,0]),lag)),"| D(1,-1,0,0) = leading:",ok(np.allclose(Dd@np.array([1.,-1,0,0]),lead)))
print("copy swap v->-v turns lagging into leading:",ok(np.allclose(np.diag([1.,1,-1,-1])@lag,lead)))
# two-rate generator: -w + B M_k B^-1 has columns summing to zero
for eta in [0.5,1.2]:
    ck=a*np.sinh(eta); w=a*np.cosh(eta); B=np.diag([np.exp(-eta/2),np.exp(eta/2)]); Mk=np.array([[-ck,a],[a,ck]])
    Gen=-w*np.eye(2)+B@Mk@np.linalg.inv(B)
    print(f"eta={eta}: generator of the moving observer's two-rate walk {Gen.round(4).tolist()}  columns sum to 0:",ok(np.allclose(Gen.sum(axis=0),0)),"rates a e^(+-eta):",ok(np.isclose(Gen[1,0],a*np.exp(eta)) and np.isclose(Gen[0,1],a*np.exp(-eta))))
