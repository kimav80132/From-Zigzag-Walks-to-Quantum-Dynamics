"""The Laplace weight e^{-ckt} of the x-clock count is III.15c's boost of the rest count.
Checks (exact linear algebra, no solving method):
 1. The rest count (1,1)e^{ax/c}, boosted by rapidity eta with the half-rapidity rescaling
    of the counters (III.15c (6)), solves the x-clock lines c d/dx' Q' = (-A_z d/dt' + a sigma_x) Q'
    with growth a cosh(eta) along the new clock and slope a sinh(eta) in the new time.
 2. The growing eigenvalue of M_k = [[-ck,a],[a,ck]] is a cosh(eta) and its eigenvector ratio
    Q_-/Q_+ is e^{eta}, with ck = a sinh(eta): the tilted count's slot is the boosted rest slot.
 3. hbar k = m c sinh(eta), hbar omega_k = m c^2 cosh(eta): the tilt is the momentum, the growth the energy.
 4. The boost does not conjugate M_0 into M_k (the clock changes too): eigenvalues differ by cosh(eta).
"""
import numpy as np
a=1.0; c=1.0; hbar=1.0; m=a*hbar/c**2
Az=np.diag([1.,-1.]); sx=np.array([[0.,1.],[1.,0.]])
npass=0
for eta in [0.0,0.3,0.5,1.2,2.0,-0.7]:
    g,gb=np.cosh(eta),np.sinh(eta); ck=a*gb; w=a*g
    v=np.array([np.exp(-eta/2),np.exp(eta/2)])
    ok1=np.allclose(w*v,(-Az*ck+a*sx)@v)                       # check 1
    Mk=np.array([[-ck,a],[a,ck]]); ev,V=np.linalg.eig(Mk); i=np.argmax(ev.real)
    u=V[:,i]/V[0,i]
    ok2=np.isclose(ev[i].real,a*g) and np.isclose(u[1].real,np.exp(eta))   # check 2
    ok3=np.isclose(hbar*ck/c,m*c*gb) and np.isclose(hbar*ev[i].real,m*c**2*g)  # check 3
    B=np.diag([np.exp(-eta/2),np.exp(eta/2)])
    ok4=(eta==0.0) or (not np.allclose(np.sort(np.linalg.eigvals(B@(a*sx)@np.linalg.inv(B)).real),np.sort(ev.real)))  # check 4
    print(f"eta={eta:+.2f}: ck={ck:+.4f} omega={w:.4f} Q-/Q+={u[1].real:.4f}=e^eta  checks:",ok1,ok2,ok3,ok4)
    npass+=sum([ok1,ok2,ok3,ok4])
print(f"{npass} of 24 PASS")
raise SystemExit(0 if npass==24 else 1)
