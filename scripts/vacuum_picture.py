import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
rng=np.random.default_rng(5)
N=12; th=0.6; cc,ss=np.cos(th),np.sin(th); D=2*N
U=np.zeros((D,D))
for y in range(N):
    U[2*y,2*((y-1)%N)]+=cc;  U[2*y,2*((y+1)%N)+1]+=ss
    U[2*y+1,2*((y-1)%N)]+=-ss; U[2*y+1,2*((y+1)%N)+1]+=cc
ev,V=np.linalg.eig(U); eps=-np.angle(ev); assert np.abs(eps).min()>1e-6
Phi,_=np.linalg.qr(V[:,eps<0]); Pm=Phi@Phi.conj().T
print("walkers in the vacuum:",Phi.shape[1],"on",N,"sites (",D,"states )")
print("occupation of every (site,sense) state: min %.6f max %.6f"%(Pm.diagonal().real.min(),Pm.diagonal().real.max()))
x=3; pp=Pm[2*x,2*x].real; mm=Pm[2*x+1,2*x+1].real; pmx=abs(Pm[2*x,2*x+1])**2
both=pp*mm-pmx; print("one site: P(both senses present)=%.4f  P(empty)=%.4f  P(exactly one)=%.4f   (independent would be 0.25,0.25,0.5)"%(both,both,1-2*both))
nn=Pm[2*x,2*x].real*Pm[2*(x+1),2*(x+1)].real-abs(Pm[2*x,2*(x+1)])**2
print("same sense on neighbouring sites: P=%.4f (independent 0.25)"%nn)
# sample configurations from the vacuum: |determinant|^2 (projection DPP), sequential algorithm
def sample(K):
    K=K.copy(); out=[]
    for i in range(K.shape[0]):
        p=K[i,i].real
        if rng.random()<p:
            out.append(i); K=K-np.outer(K[:,i],K[i,:])/K[i,i]
        else:
            K=K-np.outer(K[:,i],K[i,:])/(K[i,i]-1)
    return out
samples=[sample(Pm) for _ in range(6)]
print("walkers per sample:",[len(s) for s in samples])
fig,axs=plt.subplots(1,2,figsize=(16,6.5)); plt.subplots_adjust(wspace=0.22)
# panel A: beams (frames)
ax=axs[0]
cols=plt.cm.viridis(np.linspace(0.1,0.9,5))
for j,eta in enumerate([-0.9,-0.45,0.0,0.45,0.9]):
    for off in np.linspace(-3,3,4):
        xs=[0.0]; ts=[off]; e=1
        for d in [0.9,0.9]*4:
            dp=d*np.exp(-e*eta); xs.append(xs[-1]+dp); ts.append(ts[-1]+e*dp); e=-e
        xs=np.array(xs); ts=np.array(ts); keep=xs<=7.5
        ax.plot(xs[keep],ts[keep],color=cols[j],lw=1.6,alpha=0.85)
    ax.plot([],[],color=cols[j],lw=2,label="frame η = %+.2f"%eta)
ax.set_xlim(0,7.5); ax.set_ylim(-6,6); ax.set_xticks([]); ax.set_yticks([])
ax.set_xlabel("x (the clock)",fontsize=11); ax.set_ylabel("t",fontsize=11)
ax.legend(fontsize=9,loc="upper center",bbox_to_anchor=(0.5,-0.08),ncol=5,frameon=False)
ax.set_title("A. In the walk's clock: one SHEET of zigzags per frame, spread evenly over t.\nSheets of different frames pass through each other freely.\nRule: at most ONE unit in each sheet. Vacuum: every sheet holds its one unit.",fontsize=11)
# panel B: samples on a time cut
ax=axs[1]
for r,s in enumerate(samples):
    y=len(samples)-r
    ax.plot([-0.5,N-0.5],[y,y],color="lightgray",lw=0.8)
    for X in s:
        site,sense=divmod(X,2)
        if sense==0: ax.annotate("",xy=(site+0.28,y+0.33),xytext=(site-0.02,y+0.03),arrowprops=dict(arrowstyle="-|>",color="tab:blue",lw=2))
        else:        ax.annotate("",xy=(site-0.28,y+0.33),xytext=(site+0.02,y+0.03),arrowprops=dict(arrowstyle="-|>",color="tab:red",lw=2))
    ax.text(N-0.2,y+0.1,"%d walkers"%len(s),fontsize=9,va="center")
ax.set_xlim(-0.8,N+1.6); ax.set_ylim(0.5,len(samples)+0.9); ax.set_yticks([]); ax.set_xticks(range(N)); ax.set_xlabel("site x  (a cut at one time)",fontsize=11)
ax.set_title("B. On a time cut: six snapshots drawn from the same vacuum.\nBlue = moving right, red = moving left. Two senses may share a site (they cross);\nnever two of the same sense on a site. Always exactly %d walkers on %d sites."%(Phi.shape[1],N),fontsize=11)
for sp_ in ["top","right","left"]: ax.spines[sp_].set_visible(False)
plt.savefig("vacuum_picture.png",dpi=150,bbox_inches="tight")
# 10 Nov 2026: the numbers the text quotes, tested (the script only printed them).
dg=Pm.diagonal().real
c1=abs(dg-0.5).max()<1e-12; c2=all(len(s_)==N for s_ in samples) and all(len(set(s_))==len(s_) for s_ in samples)
c3=abs(nn-0.1965)<5e-4 and nn<0.25; c4=abs(both-0.25)<1e-12
c5=np.abs(Pm.conj()-(np.eye(D)-Pm)).max()<1e-12
ok_all=c1 and c2 and c3 and c4 and c5
print("%s  every (site, sense) state occupied with probability 1/2 (to %.0e), from conj(Pi_-) = 1 - Pi_- (%s); every sample has exactly %d walkers, no state twice (%s); "
      "same sense on neighbouring sites %.4f < 0.25 (%s); both senses on one site %.4f on this tick (%s)"
      %("PASS" if ok_all else "FAIL",abs(dg-0.5).max(),c5,N,c2,nn,c3,both,c4))
raise SystemExit(0 if ok_all else 1)
