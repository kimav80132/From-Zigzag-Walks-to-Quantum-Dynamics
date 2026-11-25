"""The vacuum on a time cut, from the document's continuum lines h(k) = ck sigma_z - a sigma_y with exact k on a ring (odd N, so k <-> -k
pairs exactly), compared with III.15a's lattice tick. Purpose: separate what is general from what is the tick's sublattice artefact."""
import numpy as np
def lines(N,a):
    ks=2*np.pi*np.fft.fftfreq(N); xs=np.arange(N); Pm=np.zeros((2*N,2*N),complex)
    for k in ks:
        w=np.sqrt(a*a+k*k); h=np.array([[k,1j*a],[-1j*a,-k]]); ph=np.exp(1j*k*xs)/np.sqrt(N)
        Pm+=np.kron(np.outer(ph,ph.conj()),0.5*(np.eye(2)-h/w))
    return Pm
def tick(N,th):
    c,s=np.cos(th),np.sin(th); D=2*N; U=np.zeros((D,D))
    for y in range(N):
        U[2*y,2*((y-1)%N)]+=c; U[2*y,2*((y+1)%N)+1]+=s; U[2*y+1,2*((y-1)%N)]+=-s; U[2*y+1,2*((y+1)%N)+1]+=c
    ev,V=np.linalg.eig(U); eps=-np.angle(ev); Phi,_=np.linalg.qr(V[:,eps<0]); return Phi@Phi.conj().T
def stats(Pm,x=3):
    pp,mm,pm=Pm[2*x,2*x].real,Pm[2*x+1,2*x+1].real,abs(Pm[2*x,2*x+1])**2
    nn=Pm[2*x,2*x].real*Pm[2*x+2,2*x+2].real-abs(Pm[2*x,2*x+2])**2
    return Pm.diagonal().real.min(),Pm.diagonal().real.max(),pp*mm-pm,nn
print("%-34s %-22s %-22s %-22s"%("","occupation of a state","P(both senses, 1 site)","P(same sense, neighbours)"))
for N,a in [(13,0.6),(49,0.6),(49,0.15),(193,0.15)]:
    lo,hi,both,nn=stats(lines(N,a)); print("%-34s %.6f..%.6f      %.4f                 %.4f"%("continuum lines  N=%d a=%.2f"%(N,a),lo,hi,both,nn))
for N,th in [(12,0.6),(48,0.6),(48,0.15)]:
    lo,hi,both,nn=stats(tick(N,th)); print("%-34s %.6f..%.6f      %.4f                 %.4f"%("lattice tick     N=%d theta=%.2f"%(N,th),lo,hi,both,nn))
print("independent walkers would give     0.5                    0.2500                 0.2500")
# 10 Nov 2026: tested.  (a) the two-state formula P(both) = n_1 n_2 - |<1|Pi|2>|^2 against the plain list of configuration probabilities |det|^2
# on a small ring; (b) the numbers the text quotes; (c) they do not depend on the ring's length.
import itertools
def brute(Pm, i, j):
    w, V = np.linalg.eigh(Pm); Phi = V[:, w > 0.5]; n = Phi.shape[1]; tot = 0.0; both = 0.0
    for idx in itertools.combinations(range(Pm.shape[0]), n):
        pr = abs(np.linalg.det(Phi[list(idx), :]))**2; tot += pr
        if i in idx and j in idx: both += pr
    return both, tot
Pb = lines(5, 0.6); x = 3
b_both, tot = brute(Pb, 2*x, 2*x + 1); b_nn, _ = brute(Pb, 2*x, 2*((x + 1) % 5))
f_both = Pb[2*x, 2*x].real*Pb[2*x + 1, 2*x + 1].real - abs(Pb[2*x, 2*x + 1])**2
f_nn = Pb[2*x, 2*x].real*Pb[2*((x + 1) % 5), 2*((x + 1) % 5)].real - abs(Pb[2*x, 2*((x + 1) % 5)])**2
ca = abs(tot - 1) < 1e-12 and abs(b_both - f_both) < 1e-12 and abs(b_nn - f_nn) < 1e-12
s13, s49, s49b, s193 = stats(lines(13, 0.6)), stats(lines(49, 0.6)), stats(lines(49, 0.15)), stats(lines(193, 0.15)); t12, t48 = stats(tick(12, 0.6)), stats(tick(48, 0.6))
cb = abs(s49[2] - 0.1993) < 5e-4 and abs(s193[2] - 0.2420) < 5e-4 and abs(t12[2] - 0.25) < 1e-12 and abs(t12[3] - 0.1965) < 5e-4 \
     and all(abs(v - 0.5) < 1e-12 for st in (s13, s49, s49b, s193, t12, t48) for v in st[:2])
cc = abs(s13[2] - s49[2]) < 1e-4 and abs(s49b[2] - s193[2]) < 1e-4 and abs(t12[3] - t48[3]) < 1e-4
ok_all = ca and cb and cc
print("%s  (a) on 5 sites the 252 configuration probabilities |det|^2 sum to 1 and give P(both senses) = %.6f, P(same sense, neighbours) = %.6f, the two-state formula "
      "%.6f, %.6f (%s); (b) continuum lines: both senses on a site 0.199 at a = 0.6 and 0.242 at a = 0.15, exactly one walker on a site %.2f and %.2f; "
      "the lattice tick 0.2500 exactly and 0.197 for neighbours (%s); (c) independent of the ring's length to 1e-4 (%s)"
      % ("PASS" if ok_all else "FAIL", b_both, b_nn, f_both, f_nn, ca, 1 - 2*s49[2], 1 - 2*s193[2], cb, cc))
raise SystemExit(0 if ok_all else 1)
