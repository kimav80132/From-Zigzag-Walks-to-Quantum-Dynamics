"""real_corners_against_feynmans_i.py
The walk read in t: a corner forward->backward weighs -a dt, backward->forward weighs +a dt (both real).
Feynman's checkerboard: every corner weighs the same imaginary number (here -i a dt; his sign convention differs).
Claim: the two sums over zigzags are the same thing, with (Feynman's backward amount) = i x (the walk's backward amount)."""
import numpy as np
a,c,T,N=1.3,1.0,2.0,400; dt=T/N; ok=True
def check(n,cnd):
    global ok; print(("PASS " if cnd else "FAIL ")+n); ok=ok and bool(cnd)
def run(wFB,wBF,dtype):
    F=np.zeros(2*N+3,dtype); B=np.zeros(2*N+3,dtype); F[N+1]=1
    for n in range(N):
        F=np.roll(F,1); B=np.roll(B,-1)
        F,B=F+wBF*B,B+wFB*F          # wFB: forward into backward; wBF: backward into forward
    return F,B
Fw,Bw=run(-a*dt,+a*dt,float)                   # the walk, real
Ff,Bf=run(-1j*a*dt,-1j*a*dt,complex)           # one imaginary weight for every corner
check("1 the walk's two counters stay real", np.isrealobj(Fw) and np.isrealobj(Bw))
check("2 forward amounts agree: Feynman's = the walk's (largest difference %.1e)"%np.max(np.abs(Ff-Fw)), np.allclose(Ff,Fw,atol=1e-12))
check("3 backward amounts: Feynman's = i x the walk's (largest difference %.1e)"%np.max(np.abs(Bf-1j*Bw)), np.allclose(Bf,1j*Bw,atol=1e-12))
check("4 two corners in a row weigh the same in both: (-a dt)(+a dt) = (-i a dt)^2 = -a^2 dt^2", np.isclose((-a*dt)*(a*dt),(-1j*a*dt)**2))
check("5 with ONE real weight for every corner (+a dt, +a dt) the sum is a different one (it grows: the telegraph's exchange)",
      not np.allclose(run(a*dt,a*dt,float)[0],Fw,atol=1e-6))
print("ALL PASS" if ok else "SOME FAILED")
