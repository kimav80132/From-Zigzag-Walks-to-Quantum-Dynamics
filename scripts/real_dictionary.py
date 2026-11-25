import numpy as np
from scipy.linalg import expm
a=1.0
I2=np.eye(2); sx=np.array([[0.,1.],[1.,0.]]); T=np.array([[0.,1.],[-1.,0.]]); J=np.array([[0.,-1.],[1.,0.]])
# column r=(u_+,v_+,u_-,v_-): index = sense ⊗ copy
Sx=np.kron(sx,I2)          # x-clock exchange on both copies
G =np.kron(T,I2)*a         # t-clock real generator (III.5a (H) at rest), one T per copy
Jc=np.kron(I2,J)           # quarter turn between the copies (route<->twin)
Dm=np.block([[I2,np.zeros((2,2))],[np.zeros((2,2)),J]])   # quarter-turn the BACKWARD counter's copy pair
print("1. D^-1 (Jc G) D == -a Sx :",np.allclose(np.linalg.inv(Dm)@Jc@G@Dm,-a*Sx))
tau=0.7
print("   e^{-tau Jc G} == D e^{tau a Sx} D^-1 :",np.allclose(expm(-tau*Jc@G),Dm@expm(tau*a*Sx)@np.linalg.inv(Dm)))
# the two x-combinations, route only (v=0): (u_+,u_-)=(1,±1)
for name,u in [("positive count (1,1)",np.array([1.,1.])),("anti-count (1,-1)",np.array([1.,-1.]))]:
    r=np.array([u[0],0.,u[1],0.]); rt=Dm@r
    uu=rt[[0,2]]; vv=rt[[1,3]]
    lead = np.allclose(vv, T@uu); lag=np.allclose(vv,-T@uu)
    # energy sign in real terms: G r = -a Jc r  <=> E>0 ;  G r = +a Jc r <=> E<0
    Epos=np.allclose(G@rt,-a*Jc@rt); Eneg=np.allclose(G@rt,a*Jc@rt)
    print(f"2. {name}: after D -> u={uu}, v={vv};  v=+Tu (twin leads): {lead};  v=-Tu (twin lags): {lag};  E>0: {Epos}  E<0: {Eneg}")
# one real column alone: both combinations just rotate under e^{tG}, no distinction
for u in [np.array([1.,1.]),np.array([1.,-1.])]:
    q=expm(0.4*a*T)@u; print(f"3. one real column {u} after t=0.4: {q.round(4)}  (length {np.linalg.norm(q):.4f}, angle turned {np.degrees(np.arctan2(q[1],q[0])-np.arctan2(u[1],u[0])):.1f} deg)")
# (10 Nov 2026) the same statements as tested lines
good = np.allclose(np.linalg.inv(Dm)@Jc@G@Dm,-a*Sx) and np.allclose(expm(-tau*Jc@G),Dm@expm(tau*a*Sx)@np.linalg.inv(Dm))
res = {}
for name,u in [("grow",np.array([1.,1.])),("decay",np.array([1.,-1.]))]:
    rt=Dm@np.array([u[0],0.,u[1],0.]); uu=rt[[0,2]]; vv=rt[[1,3]]
    res[name]=(np.allclose(vv,T@uu),np.allclose(vv,-T@uu),np.allclose(G@rt,-a*Jc@rt),np.allclose(G@rt,a*Jc@rt))
good = good and res["grow"]==(False,True,False,True) and res["decay"]==(True,False,True,False)
print(("PASS  " if good else "FAIL  ")+"rest case of III.15e (5): the combination the count grows in, (1,1), is the lagging loop (v = -Tu, E < 0); the one that decays, (1,-1), is the leading loop (v = +Tu, E > 0)")
raise SystemExit(0 if good else 1)
