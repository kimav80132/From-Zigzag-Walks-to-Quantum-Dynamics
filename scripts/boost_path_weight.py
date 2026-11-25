import numpy as np
rng=np.random.default_rng(7); a=1.0; c=1.0
# a zigzag in the rest frame: senses eps_j, rest-clock segment lengths d_j (last segment open: survival only)
def weights(eps,d,eta):
    g,gb=np.cosh(eta),np.sinh(eta); ck=a*gb; w=a*g
    # primed clock lengths: forward segment shortened by e^{-eta}, backward lengthened by e^{+eta}
    dp=np.array([dj*np.exp(-ej*eta) for ej,dj in zip(eps,d)])
    S=dp.sum(); dt_p=(eps*dp).sum()
    # boosted process as the primed observer reads it: hazard a e^{eps eta} per unit primed clock in sense eps
    haz=a*np.exp(eps*eta)
    boosted=np.prod(haz[:-1])*np.exp(-(haz*dp).sum())
    symmetric=a**(len(d)-1)*np.exp(-a*S)
    formula=np.exp(-ck*dt_p)*np.exp(-(w-a)*S)*np.exp(eta*(eps[0]-eps[-1])/2)
    return boosted/symmetric, formula
ok=0
for trial in range(2000):
    m=rng.integers(1,9); e0=rng.choice([1,-1]); eps=np.array([e0*(-1)**j for j in range(m)])
    d=rng.exponential(1.0,size=m); eta=rng.uniform(-2,2)
    r,f=weights(eps,d,eta); ok+=np.isclose(r,f,rtol=1e-12)
print(ok,"of 2000 random zigzags: ratio == e^{-ck dt'} e^{-(w-a)S} e^{eta(eps_i-eps_f)/2}")
# the alternative sign e^{eta(eps_i+eps_f)/2} fails:
eps=np.array([1,-1,1]); d=np.array([0.7,1.1,0.4]); eta=0.9
r,f=weights(eps,d,eta); print("example: ratio",r," formula(e_i - e_f)",f," wrong-sign variant",f*np.exp(eta*eps[-1]))
# (10 Nov 2026) both statements as tested lines: the formula on all 2000, and the variant with eps_i + eps_f failing
rng=np.random.default_rng(8); bad=0; tot=0
for trial in range(2000):
    m=rng.integers(1,9); e0=rng.choice([1,-1]); eps=np.array([e0*(-1)**j for j in range(m)]); d=rng.exponential(1.0,size=m); eta=rng.uniform(-2,2)
    r,f=weights(eps,d,eta); tot+=1; bad+=not np.isclose(r,f*np.exp(eta*eps[-1]),rtol=1e-6)
good=(ok==2000) and bad>0.95*tot
print(("PASS  " if good else "FAIL  ")+f"(11) holds on 2000 of 2000 random zigzags to 1e-12; the variant with eps_i + eps_f fails on {bad} of {tot}")
raise SystemExit(0 if good else 1)
