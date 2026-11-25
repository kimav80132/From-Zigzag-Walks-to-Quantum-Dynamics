"""IV.7: the streaming term of one copy in the variables (rho, S, theta, phi).  The text says "multiply out; checked" and cites no script.
(New, 10 Nov 2026.)   psi = sqrt(rho) e^{iS/hbar} chi,  chi = (cos theta/2, e^{i phi} sin theta/2),  s = chi^dag sigma chi,  all four generic fields.
1. Re[chi^dag sigma_i d_i chi] = (1/2) div s;   Im[chi^dag sigma_i d_i chi] = sin^2(theta/2) s.grad(phi) + (1/2) s.(curl s).
2. psi^dag sigma_i (-i hbar d_i) psi = -(i hbar/2) s.grad(rho) + rho s.grad(S) - i hbar rho chi^dag sigma_i d_i chi.
3. real part = s.P with P = rho (grad S + hbar sin^2(theta/2) grad phi) + (hbar/2) curl(rho s);  imaginary part = -(hbar/2) div(rho s).
4. s.(grad rho x s) = 0, the identity used to write rho s.(curl s) = s.curl(rho s)."""
import sympy as sp
ok_all = True
def ok(n, b):
    global ok_all; ok_all &= bool(b); print(("PASS  " if b else "FAIL  ") + n)
x, y, z = X = sp.symbols('x y z', real=True); hbar = sp.symbols('hbar', positive=True)
rho = sp.Function('rho', positive=True)(*X); S, half, ph = [sp.Function(n, real=True)(*X) for n in ('S', 'u', 'phi')]
th = 2*half                                          # the half-angle u = theta/2 is the variable, so that only whole multiples of u appear
i = sp.I
sig = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -i], [i, 0]]), sp.diag(1, -1)]
chi = sp.Matrix([sp.cos(th/2), sp.exp(i*ph)*sp.sin(th/2)]); chid = sp.Matrix([[sp.cos(th/2), sp.exp(-i*ph)*sp.sin(th/2)]])   # the conjugate row, fields real
s = sp.Matrix([sp.sin(th)*sp.cos(ph), sp.sin(th)*sp.sin(ph), sp.cos(th)])
simp = lambda e: sp.simplify(sp.expand(sp.expand_trig(sp.expand(e.rewrite(sp.cos)))))
conj = lambda e: e.subs(i, -i)                      # every field is real, so conjugation is i -> -i
Re = lambda e: sp.expand((e + conj(e))/2); Im = lambda e: sp.expand((e - conj(e))/(2*i))
ok("0. chi^dag sigma chi = (sin theta cos phi, sin theta sin phi, cos theta)", all(simp((chid*sig[k]*chi)[0] - s[k]) == 0 for k in range(3)))
grad = lambda f: sp.Matrix([sp.diff(f, v) for v in X]); div = lambda V: sum(sp.diff(V[k], X[k]) for k in range(3))
curl = lambda V: sp.Matrix([sp.diff(V[2], y) - sp.diff(V[1], z), sp.diff(V[0], z) - sp.diff(V[2], x), sp.diff(V[1], x) - sp.diff(V[0], y)])
core = sum((chid*sig[k]*sp.diff(chi, X[k]))[0] for k in range(3))
re_core, im_core = Re(core), Im(core)
ok("1a. Re[chi^dag sigma_i d_i chi] = (1/2) div s", simp(re_core - div(s)/2) == 0)
ok("1b. Im[chi^dag sigma_i d_i chi] = sin^2(theta/2) s.grad(phi) + (1/2) s.(curl s)", simp(im_core - (sp.sin(th/2)**2*s.dot(grad(ph)) + s.dot(curl(s))/2)) == 0)
psi = sp.sqrt(rho)*sp.exp(i*S/hbar)*chi; psid = sp.sqrt(rho)*sp.exp(-i*S/hbar)*chid
full = sum((psid*sig[k]*(-i*hbar)*sp.diff(psi, X[k]))[0] for k in range(3))
ok("2. psi^dag sigma_i (-i hbar d_i) psi = -(i hbar/2) s.grad(rho) + rho s.grad(S) - i hbar rho chi^dag sigma_i d_i chi  (the product rule, with chi^dag sigma_i chi = s_i)",
   sp.expand(full - (-i*hbar/2*sum((chid*sig[k]*chi)[0]*sp.diff(rho, X[k]) for k in range(3)) + rho*sum((chid*sig[k]*chi)[0]*sp.diff(S, X[k]) for k in range(3))
                     - i*hbar*rho*core)) == 0)
P = rho*(grad(S) + hbar*sp.sin(th/2)**2*grad(ph)) + hbar/2*curl(rho*s)
# by 2. and 1.: Re = rho s.grad S + hbar rho Im[core],  Im = -(hbar/2) s.grad rho - hbar rho Re[core]
re_full = rho*s.dot(grad(S)) + hbar*rho*(sp.sin(th/2)**2*s.dot(grad(ph)) + s.dot(curl(s))/2)
im_full = -hbar/2*s.dot(grad(rho)) - hbar*rho*div(s)/2
ok("3a. Re[psi^dag sigma.p psi] = s.P with P = rho (grad S + hbar sin^2(theta/2) grad phi) + (hbar/2) curl(rho s)", simp(re_full - s.dot(P)) == 0)
ok("3b. Im[psi^dag sigma.p psi] = -(hbar/2) div(rho s), a divergence: it integrates to zero", simp(im_full + hbar/2*div(rho*s)) == 0)
ok("4. s.(grad rho x s) = 0, so rho s.(curl s) = s.curl(rho s)", simp(s.dot(grad(rho).cross(s))) == 0 and simp(rho*s.dot(curl(s)) - s.dot(curl(rho*s))) == 0)
# the same four statements at a point, on explicit fields, as numbers (a check on the symbolic simplifier itself)
num = {half: (sp.Rational(9, 10) + sp.Rational(3, 10)*x*y - z/5 + sp.sin(x + 2*z)/10)/2, ph: sp.Rational(2, 5)*x - sp.Rational(7, 10)*y*z + sp.cos(y)/5,
       rho: 1 + x**2/5 + sp.sin(y*z)/10, S: sp.Rational(3, 10)*x*z + y/2}
pt = {x: 0.37, y: -0.81, z: 0.52, hbar: 0.7}
val = lambda e: complex(e.subs(num).doit().subs(pt).evalf())
errs = [abs(val(core).real - val(div(s)/2).real), abs(val(core).imag - val(sp.sin(th/2)**2*s.dot(grad(ph)) + s.dot(curl(s))/2).real),
        abs(val(full).real - val(s.dot(P)).real), abs(val(full).imag - val(-hbar/2*div(rho*s)).real)]
ok(f"5. 1a, 1b, 3a, 3b directly from psi on explicit fields at a point, as numbers: errors {max(errs):.1e}", max(errs) < 1e-14)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
