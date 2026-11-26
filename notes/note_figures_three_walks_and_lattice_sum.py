"""make_figs.py -- figures for the note "From a random walk to the Dirac propagator".
Figure 1: three sampled walks (telegraph; an entwined pair built by Ord and Mann's rules, quant-ph/0208004 pp. 4-6; this note's walk).
Figure 2: the lattice sum over zigzags (step 9's weights) against -(a/2c) J_0(zeta) at a t = 12.
With --write it writes fig1.html, fig2.html (the figures as HTML fragments) and figdata.json; without it, it only runs its checks."""
import sys
WRITE="--write" in sys.argv
import numpy as np, json
from scipy.special import j0
W=160; OX,OY=30,196; VBW,VBH=204,226
X=lambda u: round(OX+u,1); Y=lambda v: round(OY-v,1)
def axes(hlabel, vlabel):
    s='<defs><marker id="ahM" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" class="fig-muted"/></marker></defs>'
    s+=f'<line x1="{OX-6}" y1="{OY}" x2="{OX+W+10}" y2="{OY}" class="fig-rule" stroke-width="1" marker-end="url(#ahM)"/>'
    s+=f'<line x1="{OX}" y1="{OY+6}" x2="{OX}" y2="{OY-W-10}" class="fig-rule" stroke-width="1" marker-end="url(#ahM)"/>'
    s+=f'<text x="{OX+W+8}" y="{OY+19}" text-anchor="end" font-size="13" font-style="italic" class="fig-muted">{hlabel}</text>'
    s+=f'<text x="{OX-8}" y="{OY-W-14}" text-anchor="start" font-size="13" font-style="italic" class="fig-muted">{vlabel}</text>'
    return s
def arrow(p,q,cls,size=7.5,frac=0.5):
    (x1,y1),(x2,y2)=p,q; mx,my=x1+(x2-x1)*frac,y1+(y2-y1)*frac
    d=np.array([x2-x1,y2-y1],float); d/=np.linalg.norm(d); n=np.array([-d[1],d[0]])
    tip=np.array([mx,my])+d*size*0.6; b1=np.array([mx,my])-d*size*0.5+n*size*0.5; b2=np.array([mx,my])-d*size*0.5-n*size*0.5
    return f'<path d="M{tip[0]:.1f},{tip[1]:.1f} L{b1[0]:.1f},{b1[1]:.1f} L{b2[0]:.1f},{b2[1]:.1f} z" class="{cls}"/>'
def poly(pts,cls,dash=None,w=2.2):
    d=' '.join(f'{x},{y}' for x,y in pts); da=f' stroke-dasharray="{dash}"' if dash else ''
    return f'<polyline points="{d}" fill="none" class="{cls}" stroke-width="{w}" stroke-linejoin="round" stroke-linecap="round"{da}/>'
dot=lambda x,y,cls='fig-fill',r=3.2: f'<circle cx="{x}" cy="{y}" r="{r}" class="{cls}"/>'

# ---------- sample for A (and C = A with x and ct exchanged)
L=150
def sampleA(seed):
    rng=np.random.default_rng(seed); u=[rng.uniform(55,85)]; v=[0.0]; d=1
    while True:
        tau=rng.exponential(42.0)
        if v[-1]+tau>=L: u.append(u[-1]+d*(L-v[-1])); v.append(L); break
        u.append(u[-1]+d*tau); v.append(v[-1]+tau); d=-d
    return np.array(u),np.array(v)
def goodA(u,v):
    seg=np.diff(v); nflip=len(u)-2
    if nflip not in (3,4): return None
    if seg.min()<22 or u.min()<10 or u.max()>150: return None
    # a level (value of u) crossed by at least 3 segments, 9 away from every corner
    best=None
    for lev in np.arange(20,141,1.0):
        cr=[i for i in range(len(u)-1) if (u[i]-lev)*(u[i+1]-lev)<0]
        if len(cr)==3 and np.min(np.abs(u-lev))>=11:
            m=np.min(np.abs(u-lev))
            if best is None or m>best[0]: best=(m,lev,cr)
    return best
for seedA in range(1,4000):
    u,v=sampleA(seedA); g=goodA(u,v)
    if g: break
_,lev,cr=g
# panel A: horizontal = x (u), vertical = t (v): telegraph
ptsA=[(X(a),Y(b)) for a,b in zip(u,v)]
sA=axes('x','t  (clock)')+poly(ptsA,'fig-ink')+''.join(dot(*p) for p in ptsA[1:-1])+dot(*ptsA[0],r=3.6)+arrow(ptsA[0],ptsA[1],'fig-fill')
# panel C: horizontal = x (= A's v), vertical = t (= A's u)
ptsC=[(X(b),Y(a)) for a,b in zip(u,v)]
sC=axes('x  (clock)','t')
ycross=Y(lev)
sC+=f'<line x1="{OX}" y1="{ycross}" x2="{OX+W+2}" y2="{ycross}" class="fig-acc" stroke-width="1.2" stroke-dasharray="5 5"/>'
for i in range(len(ptsC)-1):
    fwd=u[i+1]>u[i]
    sC+=poly([ptsC[i],ptsC[i+1]],'fig-fwd' if fwd else 'fig-bwd',None if fwd else '6 4')
for i in range(len(ptsC)-1):
    if i in (0,1): sC+=arrow(ptsC[i],ptsC[i+1],'fig-fwdfill' if u[i+1]>u[i] else 'fig-bwdfill')
sC+=''.join(dot(*p) for p in ptsC[1:-1])+dot(*ptsC[0],r=3.6)
for i in cr:
    f=(lev-u[i])/(u[i+1]-u[i]); xc=X(v[i]+f*(v[i+1]-v[i])); sC+=f'<circle cx="{xc}" cy="{ycross}" r="4.2" class="fig-accfill"/>'
sC+=f'<text x="{OX+W+2}" y="{ycross-6}" font-size="12" text-anchor="end" class="fig-accfill">one time</text>'

# ---------- sample for B: Ord and Mann's Eve (pp. 4-6): odd marks turn, even marks drop a marker; after t_R stop at the next marker
tR=112
def sampleB(seed):
    rng=np.random.default_rng(seed); u0=rng.uniform(60,100); anchors=[(u0,0.0)]; blue=[(u0,0.0)]; red=[(u0,0.0)]; d=1
    while True:
        t1,t2=rng.exponential(26.0,2)
        (ua,va)=anchors[-1]
        corner=(ua+d*t1,va+t1); nxt=(ua+d*t1-d*t2,va+t1+t2); rc=(ua-d*t2,va+t2)
        blue+= [corner,nxt]; red+=[rc,nxt]; anchors.append(nxt); d=-d
        if nxt[1]>=tR: break
    return anchors,blue,red
def goodB(an,bl,rd):
    if len(an)!=4: return False
    allp=np.array(bl+rd)
    if allp[:,0].min()<10 or allp[:,0].max()>150 or allp[:,1].max()>152: return False
    for path in (bl,rd):
        p=np.array(path); seg=np.abs(np.diff(p[:,1]))
        if seg.min()<15: return False
    # the two corners of a diamond well apart
    for k in range(len(an)-1):
        b=np.array(bl[2*k+1]); r=np.array(rd[2*k+1])
        if np.linalg.norm(b-r)<26: return False
    return True
for seedB in range(1,20000):
    an,bl,rd=sampleB(seedB)
    if goodB(an,bl,rd): break
P=lambda p:(X(p[0]),Y(p[1]))
sB=axes('x','t  (clock)')+poly([P(p) for p in rd],'fig-bwd','6 4')+poly([P(p) for p in bl],'fig-fwd')
sB+=arrow(P(bl[0]),P(bl[1]),'fig-fwdfill')+arrow(P(rd[-1]),P(rd[-2]),'fig-bwdfill')
sB+=''.join(dot(*P(p)) for p in an[1:])+dot(*P(an[0]),r=3.6)
ytr=Y(tR); sB+=f'<line x1="{OX-4}" y1="{ytr}" x2="{OX+4}" y2="{ytr}" class="fig-rule" stroke-width="1"/><text x="{OX+8}" y="{ytr+4}" font-size="11" class="fig-muted" font-style="italic">t</text><text x="{OX+13}" y="{ytr+7}" font-size="8" class="fig-muted">R</text>'

def panel(title,sub,svg,label):
    return (f'<div class="panel"><p class="ptitle">{title}</p><p class="psub">{sub}</p>'
            f'<svg viewBox="0 0 {VBW} {VBH}" role="img" aria-label="{label}">{svg}</svg></div>')
fig1='<figure>\n<div class="panels">\n'
fig1+=panel('<b>A</b>&ensp;Telegraph process','clock $t$ · the direction flips',sA,'A path that always advances in time and reverses its direction in space at random.')+'\n'
fig1+=panel('<b>B</b>&ensp;Entwined pair','clock $t$ · one turn back in time',sB,'A forward leg that zigzags up in time, and a return leg that comes back to the origin through the markers.')+'\n'
fig1+=panel('<b>C</b>&ensp;This walk','clock $x$ · the time sense flips',sC,'A path that always advances in x while its time runs forward, then backward, at random. A line of constant time crosses it three times.')+'\n'
fig1+='</div>\n<div class="legend"><span><svg viewBox="0 0 30 10"><line x1="1" y1="5" x2="29" y2="5" class="fig-fwd" stroke-width="2.4" stroke-linecap="round"/></svg> time forward</span><span><svg viewBox="0 0 30 10"><line x1="1" y1="5" x2="29" y2="5" class="fig-bwd" stroke-width="2.4" stroke-dasharray="6 4" stroke-linecap="round"/></svg> time backward</span><span><svg viewBox="0 0 12 10" style="width:12px"><circle cx="6" cy="5" r="3.2" class="fig-fill"/></svg> flip (in B: marker)</span></div>\n'
fig1+='<figcaption>Three walks at the speed $c$, each drawn from a random sample. <b>A</b> is the telegraph process of Goldstein and Kac. <b>B</b> is an entwined pair of Ord and Mann: a forward leg up to the first marker after the return time $t_R$, then a return leg back to the origin through the markers. <b>C</b> is the walk of this note: the same sample as A, with $x$ and $ct$ exchanged. Someone who reads C at one time meets the path three times.</figcaption>\n</figure>\n'
if WRITE: open('fig1.html','w').write(fig1)

# ---------- figure 2: lattice sum over zigzags against -(a/2c) J0
av,cv,T,N=1.0,1.0,12.0,12000; dt=T/N; h=cv*dt
Qp=np.zeros(2*N+3); Qm=np.zeros(2*N+3); Qp[N+1]=1.0
for n in range(N):
    Qp=np.roll(Qp,1); Qm=np.roll(Qm,-1)
    Qp,Qm=Qp+av*dt*Qm,Qm-av*dt*Qp
idx=np.arange(2*N+3)-(N+1); xx=idx*h
reach=np.where(((idx+N)%2==0)&(np.abs(xx)<=0.965*cv*T))[0]
pick=reach[np.round(np.linspace(0,len(reach)-1,45)).astype(int)]
lat=Qm[pick]/(2*h)/(av/(2*cv)); xs=xx[pick]/(cv*T)
form=-j0(av/cv*np.sqrt(cv**2*T**2-xx[pick]**2)); miss=float(np.max(np.abs(lat-form)))
PX0,PX1,PY0,PY1=58,540,24,196; ymin,ymax=-1.15,0.6
fx=lambda s: round(PX0+(s+1)/2*(PX1-PX0),1); fy=lambda q: round(PY1-(q-ymin)/(ymax-ymin)*(PY1-PY0),1)
s=''
for q in (-1,-0.5,0.5): s+=f'<line x1="{PX0}" y1="{fy(q)}" x2="{PX1}" y2="{fy(q)}" class="fig-grid" stroke-width="1"/>'
s+=f'<line x1="{PX0}" y1="{fy(0)}" x2="{PX1}" y2="{fy(0)}" class="fig-rule" stroke-width="1"/>'
for q,lab in ((-1,'−1'),(-0.5,'−0.5'),(0,'0'),(0.5,'0.5')): s+=f'<text x="{PX0-8}" y="{fy(q)+4}" text-anchor="end" font-size="12" class="fig-muted">{lab}</text>'
for sx,lab in ((-1,'−1'),(-0.5,'−0.5'),(0,'0'),(0.5,'0.5'),(1,'1')):
    s+=f'<line x1="{fx(sx)}" y1="{PY1}" x2="{fx(sx)}" y2="{PY1+5}" class="fig-rule" stroke-width="1"/><text x="{fx(sx)}" y="{PY1+19}" text-anchor="middle" font-size="12" class="fig-muted">{lab}</text>'
s+=f'<text x="{PX1}" y="{PY1+36}" text-anchor="end" font-size="13" font-style="italic" class="fig-muted">x / ct</text>'
s+=f'<text x="{PX0-44}" y="{PY0-8}" text-anchor="start" font-size="13" class="fig-muted">Q₋ in units of a/2c</text>'
sg=np.linspace(-1,1,1201); curve=-j0(av*T*np.sqrt(np.clip(1-sg**2,0,None)))
s+=poly([(fx(a),fy(b)) for a,b in zip(sg,curve)],'fig-ink',w=2)
for a,b in zip(xs,lat): s+=f'<circle cx="{fx(a)}" cy="{fy(b)}" r="3.4" class="fig-fwdfill fig-ring"/>'
fig2='<figure>\n<svg viewBox="0 0 560 240" role="img" aria-label="The backward counter at a time t, across the light cone: the sum over zigzags on a lattice lies on the Bessel function.">'+s+'</svg>\n'
fig2+='<div class="legend"><span><svg viewBox="0 0 30 10"><line x1="1" y1="5" x2="29" y2="5" class="fig-ink" stroke-width="2"/></svg> the formula, $-J_0(\\zeta)$</span><span><svg viewBox="0 0 12 10" style="width:12px"><circle cx="6" cy="5" r="3.4" class="fig-fwdfill"/></svg> the sum over all zigzags on a lattice</span></div>\n'
fig2+='<figcaption>The backward counter $Q_-$ at the time $t=12/a$, across the light cone. Each point adds up every zigzag that ends there, with the corner weights of step 9. The points lie on the Bessel function to FIGMISS.</figcaption>\n</figure>\n'.replace('FIGMISS',f'{miss:.3f}')
if WRITE: open('fig2.html','w').write(fig2)

# ---------- checks on what the figures show
okall=True
def check(name,cond):
    global okall
    print(("PASS " if cond else "FAIL ")+name); okall=okall and bool(cond)
du,dv=np.diff(u),np.diff(v)
check("1 panel A: every piece moves at the speed c (|dx| = c dt) and time always advances", np.allclose(np.abs(du),dv) and (dv>0).all())
check("2 panel A: the direction reverses at every flip", all(du[i]*du[i+1]<0 for i in range(len(du)-1)))
xc=[p[0] for p in ptsC]; tc=[OY-p[1] for p in ptsC]
check("3 panel C is panel A with x and ct exchanged: x always advances, the time sense reverses at every flip",
      (np.diff(xc)>0).all() and np.allclose(np.array(tc),u,atol=0.06) and np.allclose(np.array(xc)-OX,v,atol=0.06))
check("4 panel C: the line of one time meets the path three times", len(cr)==3)
bl_=np.array(bl); rd_=np.array(rd); an_=np.array(an)
light=lambda p: np.allclose(np.abs(np.diff(p[:,0])),np.abs(np.diff(p[:,1])))
check("5 panel B: every piece of the forward leg and of the return leg is on a light line", light(bl_) and light(rd_))
check("6 panel B: the two legs meet at the origin, at every marker and at the end, and nowhere else among their corners",
      np.allclose(bl_[::2],an_) and np.allclose(rd_[::2],an_) and all(np.linalg.norm(bl_[2*k+1]-rd_[2*k+1])>1 for k in range(len(an)-1)))
check("7 panel B: the forward leg turns once between consecutive markers and passes straight through each marker (Ord and Mann: odd marks turn, even marks drop a marker)",
      all(np.sign(bl_[2*k+1,0]-bl_[2*k,0])==-np.sign(bl_[2*k+2,0]-bl_[2*k+1,0]) for k in range(len(an)-1)) and
      all(np.sign(bl_[2*k,0]-bl_[2*k-1,0])==np.sign(bl_[2*k+1,0]-bl_[2*k,0]) for k in range(1,len(an)-1)))
check("8 panel B: the end is the first marker at or after the return time t_R", an_[-1,1]>=tR and (an_[1:-1,1]<tR).all())
check("9 panel B: at the end the return leg keeps the forward leg's direction in x and reverses in t",
      np.sign(bl_[-1,0]-bl_[-2,0])==np.sign(rd_[-2,0]-rd_[-1,0]))
check("10 figure 2: the lattice sum over zigzags (corner weights -a dt, +a dt) lies on -J_0(zeta) at a t = 12, 45 sites, largest miss %.4f in units of a/2c"%miss, miss<0.005)
print("ALL PASS" if okall else "SOME FAILED")
if WRITE:
    json.dump(dict(seedA=seedA,seedB=seedB,level=float(lev),flipsA=len(u)-2,markersB=len(an)-1,miss=miss,npoints=len(pick)),open('figdata.json','w'))
    print(json.load(open('figdata.json')))
