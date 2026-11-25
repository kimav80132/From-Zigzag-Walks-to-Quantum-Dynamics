"""
Checks of two drafts (5 Nov 2026):
 interaction draft ('Do different carrier responses on the two senses give gamma(Phi - vA)?'):
  1. the field terms of V.2 (E) turn each sense's copy pair: rho_+- D_+- theta_+- gets -b_+- rho_+-, b_+- = (q/hbar)(Phi -+ c A)
  2. density-weighted: (rho+ b+ + rho- b-)/rho = (q/hbar)(Phi - v A), v = c (rho+ - rho-)/rho
  3. per proper time: gamma (Phi - v A) = Phi', the scalar potential in the walker's rest frame (boost of the four-potential)
  4. unequal flip rates a +- delta: the extra term -delta [[1,1],[1,1]] changes the norm at rate -2 delta (Q+ + Q-)^2
 plates draft, section 5 and section 2:
  5. two surface penalties: det[I + g_surf diag(l1,l2)] / isolated factors = 1 - r1 r2 e^{-2 kappa L}, r_i = l_i/(2kappa + l_i)
  6. width alone does not fix the energy sign: f(p)(w_+ + w_-)/sqrt2 has equal weights however narrow f;
     a broad envelope times the rest-positive column (1,-i)/sqrt2 is one-sign to O((c k/a)^2)
"""
import sympy as sp, numpy as np
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all &= bool(cond); print(("PASS  " if cond else "FAIL  ") + msg)
t, x = sp.symbols('t x', real=True); c, a, q, hb = sp.symbols('c a q hbar', positive=True)
Phi, A = sp.symbols('Phi A', real=True)
up, um, vp, vm = [sp.Function(n)(t, x) for n in ('up', 'um', 'vp', 'vm')]
bp = q/hb*(Phi - c*A); bm = q/hb*(Phi + c*A)
rhs = {up: -c*sp.diff(up, x) + a*um + bp*vp, vp: -c*sp.diff(vp, x) + a*vm - bp*up,
       um: c*sp.diff(um, x) - a*up + bm*vm, vm: c*sp.diff(vm, x) - a*vp - bm*um}
sub = {sp.Derivative(f, t): r for f, r in rhs.items()}
Dp = lambda f: sp.diff(f, t) + c*sp.diff(f, x); Dm = lambda f: sp.diff(f, t) - c*sp.diff(f, x)
W = up*vm - vp*um; rp = up**2 + vp**2; rm = um**2 + vm**2
e1 = sp.simplify((up*Dp(vp) - vp*Dp(up)).subs(sub) - (a*W - bp*rp))
e2 = sp.simplify((um*Dm(vm) - vm*Dm(um)).subs(sub) - (a*W - bm*rm))
K = up*um + vp*vm
e3 = sp.simplify(Dp(rp).subs(sub) - 2*a*K); e4 = sp.simplify(Dm(rm).subs(sub) + 2*a*K)
ok("1. with the field: rho+ D+ theta+ = aW - b+ rho+, rho- D- theta- = aW - b- rho- (residuals %s, %s); the densities are untouched, D+-rho+- = +-2aK (%s, %s)" % (e1, e2, e3, e4), all(e == 0 for e in (e1, e2, e3, e4)))
Rp, Rm = sp.symbols('rho_p rho_m', positive=True)
v = c*(Rp - Rm)/(Rp + Rm)
ok("2. (rho+ b+ + rho- b-)/rho - (q/hbar)(Phi - v A) = %s" % sp.simplify((Rp*bp + Rm*bm)/(Rp + Rm) - q/hb*(Phi - v*A)), sp.simplify((Rp*bp + Rm*bm)/(Rp + Rm) - q/hb*(Phi - v*A)) == 0)
# 2'. (added 10 Nov 2026; the text cites it and the 5-Nov version did not have it) (10d) and (10e) on positive-energy columns at four speeds
import numpy as _np
def _pos_column(kk, av=1.0):                      # III.9's basis, h(k) = [[ck, i a], [-i a, -ck]], c = 1
    wv, Vv = _np.linalg.eigh(_np.array([[kk, 1j*av], [-1j*av, -kk]])); return Vv[:, 1], wv[1]
worst = 0.0; shown = []
for Vn in (0.0, 0.6, 0.866, 0.943):
    gam = 1/_np.sqrt(1 - Vn*Vn); kk = gam*Vn                     # p = gamma m V with a = m = 1
    col, En = _pos_column(kk); rp_, rm_ = abs(col[0])**2, abs(col[1])**2
    vcol = (rp_ - rm_)/(rp_ + rm_)                               # c (rho+ - rho-)/rho
    Ph, Aq = 0.37, -0.52; bp_, bm_ = (Ph - Aq), (Ph + Aq)        # q/hbar = 1
    weighted = (rp_*bp_ + rm_*bm_)/(rp_ + rm_)
    L = _np.array([[gam, -gam*Vn], [-gam*Vn, gam]]); Php = (L @ _np.array([Ph, Aq]))[0]     # the four-potential boosted to the walker's frame
    worst = max(worst, abs(vcol - Vn), abs(weighted - (Ph - Vn*Aq)), abs(gam*weighted - Php)); shown.append(f"V={Vn}: v={vcol:.4f}")
ok("2'. on positive-energy columns at V = 0, 0.6c, 0.866c, 0.943c: c(rho+ - rho-)/rho = V, the density-weighted turn rate is (q/hbar)(Phi - V A), and gamma times it is "
   "(q/hbar) Phi' with Phi' from the Lorentz boost of (Phi, A): largest error %.1e (%s)" % (worst, "; ".join(shown)), worst < 1e-12)
# 3. (rebuilt 10 Nov 2026: the earlier line compared gamma(Phi - V A) with itself)  the turn per step is a scalar, and in the walker's frame dx' = 0
V = sp.Symbol('V', real=True); g = 1/sp.sqrt(1 - V**2/c**2)
Lb = sp.Matrix([[g, -g*V/c], [-g*V/c, g]])                       # boost on (c t, x) and on (Phi/c, A)
dtt, dxx = sp.symbols('dt dx', real=True)
pot = Lb*sp.Matrix([Phi/c, A]); stp = Lb*sp.Matrix([c*dtt, dxx])
inv_lab = Phi*dtt - A*dxx; inv_w = c*pot[0]*stp[0]/c - pot[1]*stp[1]
along = {dxx: V*dtt}
ok("3. Phi dt - A dx is the same number in the walker's frame; along the walker dx = V dt it is (Phi - V A) dt, and there dx' = 0, dt' = dt/gamma, "
   "so it is Phi' dtau with Phi' = gamma (Phi - V A): the field turns the walker at q Phi'/hbar per unit of its own time",
   sp.simplify(inv_lab - inv_w) == 0 and sp.simplify(stp[1].subs(along)) == 0 and sp.simplify(stp[0].subs(along)/c - dtt/g) == 0
   and sp.simplify(c*pot[0] - g*(Phi - V*A)) == 0 and sp.simplify(inv_lab.subs(along) - c*pot[0]*(dtt/g)) == 0)
d = sp.Symbol('delta', real=True); Qp, Qm = sp.symbols('Q_p Q_m', real=True)
extra = -d*sp.Matrix([[1, 1], [1, 1]])*sp.Matrix([Qp, Qm])
rate = sp.expand(2*(sp.Matrix([[Qp, Qm]])*extra)[0])
ok("4. unequal flip rates: d/dt (Q+^2 + Q-^2) from the extra term = %s = -2 delta (Q+ + Q-)^2 (V.1, second row)" % sp.factor(rate), sp.simplify(rate + 2*d*(Qp + Qm)**2) == 0)
l1, l2, kap, qL = sp.symbols('lambda1 lambda2 kappa q_L', positive=True)
gs = sp.Matrix([[1, qL], [qL, 1]])/(2*kap)
dt = sp.simplify((sp.eye(2) + gs*sp.diag(l1, l2)).det()/((1 + l1/(2*kap))*(1 + l2/(2*kap))))
ok("5. surface penalties: normalized determinant = %s = 1 - r1 r2 q_L^2" % sp.factor(dt), sp.simplify(dt - (1 - l1/(2*kap + l1)*l2/(2*kap + l2)*qL**2)) == 0)
# 6
an = 1.0
def bands(k):
    H = np.array([[k, 1j*an], [-1j*an, -k]]); w, Vv = np.linalg.eigh(H); return Vv[:, 1], Vv[:, 0]
res = []
for sig in (2.0, 20.0):
    ks = np.linspace(-6/sig, 6/sig, 401); f = np.exp(-ks**2*sig**2/2)
    neg_mix = sum(abs(fk)**2*abs(np.vdot(bands(k)[1], (bands(k)[0] + bands(k)[1])/np.sqrt(2)))**2 for k, fk in zip(ks, f))/sum(abs(f)**2)
    e_rest = np.array([1, -1j])/np.sqrt(2)
    neg_rest = sum(abs(fk)**2*abs(np.vdot(bands(k)[1], e_rest))**2 for k, fk in zip(ks, f))/sum(abs(f)**2)
    res.append((sig, neg_mix, neg_rest))
ok("6. negative-energy weight: (w+ + w-)/sqrt2 times any envelope: %s; rest-positive column times an envelope of width sigma: %s - width alone does not fix the sign; the rest-positive column does, to O((c/a sigma)^2)"
   % (", ".join("sigma=%g: %.3f" % (s, n) for s, n, _ in res), ", ".join("sigma=%g: %.5f" % (s, r) for s, _, r in res)), all(abs(n - 0.5) < 1e-9 for _, n, _ in res) and res[1][2] < 0.01*res[0][2]*10)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
