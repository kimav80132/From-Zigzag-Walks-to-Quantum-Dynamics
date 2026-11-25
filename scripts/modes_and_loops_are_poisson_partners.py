"""Casimir document, "Modes and loops are Poisson partners".   Recreated 9 Nov 2026 (doc rev 68).  Units pi hbar c / L = 1.

f(n) = |n| e^{-eps |n|}, f^(k) = 2 (eps^2 - k^2)/(eps^2 + k^2)^2;  Poisson:  sum_n f(n) = sum_m f^(2 pi m).
  1. the finite part: sum_{n>=1} n e^{-eps n} = 1/eps^2 - 1/12 + O(eps^2).
  2. the images, term by term: the m = 0 image is the bulk 1/eps^2; the m-th is f^(2 pi m) -> -2/(2 pi m)^2 as eps -> 0, the loop that winds m
     times; their sum is -1/12; and Poisson's identity itself at finite eps.
  3. the half-integer modes of the turn-matrix wall: the m-th image is multiplied by (-1)^m; sum_{n>=0} (n + 1/2) e^{-eps (n + 1/2)}
     = 1/eps^2 + 1/24 + O(eps^2), and +1/24 = -2 sum_m (-1)^m/(2 pi m)^2.
  4. the Euler-Maclaurin endpoint: summing to n0 and correcting the tail with f(n0)/2 - f'(n0)/12 gives 1/eps^2 - 1/12 at n0 = 5, 20, 100."""
import sys
import mpmath as mp
mp.mp.dps = 40
fails = 0
def ok(name, cond):
    global fails
    print(("PASS  " if cond else "FAIL  ") + name); fails += (not cond)
fh = lambda k, e: 2*(e*e - k*k)/(e*e + k*k)**2
eps_list = [mp.mpf('0.1'), mp.mpf('0.01'), mp.mpf('0.001')]
# 1
fin = [mp.nsum(lambda n: n*mp.e**(-e*n), [1, mp.inf]) - 1/e**2 for e in eps_list]
ok("1. sum_{n>=1} n e^{-eps n} - 1/eps^2 = %s, %s, %s at eps = 0.1, 0.01, 0.001: -1/12 = %s, the remainder eps^2/240"
   % (mp.nstr(fin[0], 10), mp.nstr(fin[1], 10), mp.nstr(fin[2], 12), mp.nstr(-mp.mpf(1)/12, 12)),
   all(abs(f_ + mp.mpf(1)/12 - e*e/240) < e**4 for f_, e in zip(fin, eps_list)))
# 2
e = mp.mpf('0.001')
lhs = mp.nsum(lambda n: n*mp.e**(-e*n), [1, mp.inf])
rhs = (fh(0, e) + 2*mp.nsum(lambda m: fh(2*mp.pi*m, e), [1, mp.inf]))/2
img = [fh(2*mp.pi*m, e) for m in (1, 2, 3)]
ok("2. Poisson at eps = 0.001: sum_{n>=1} n e^{-eps n} = (1/2)[f^(0) + 2 sum_{m>=1} f^(2 pi m)] (relative %.0e); f^(0)/2 = 1/eps^2; the m = 1, 2, 3 "
   "images are %s, %s, %s against -2/(2 pi m)^2 = %s, %s, %s; sum over m of -2/(2 pi m)^2 = -2 zeta(2)/(4 pi^2) = %s"
   % (abs(lhs/rhs - 1), mp.nstr(img[0], 10), mp.nstr(img[1], 10), mp.nstr(img[2], 10), *[mp.nstr(-2/(2*mp.pi*m)**2, 10) for m in (1, 2, 3)],
      mp.nstr(-2*mp.zeta(2)/(4*mp.pi**2), 12)),
   abs(lhs/rhs - 1) < 1e-25 and all(abs(i_ + 2/(2*mp.pi*m)**2) < 1e-6 for i_, m in zip(img, (1, 2, 3))) and abs(-2*mp.zeta(2)/(4*mp.pi**2) + mp.mpf(1)/12) < 1e-30)
# 3
half = [mp.nsum(lambda n: (n + mp.mpf(1)/2)*mp.e**(-e_*(n + mp.mpf(1)/2)), [0, mp.inf]) - 1/e_**2 for e_ in eps_list]
alt = -2*mp.nsum(lambda m: (-1)**int(m)/(2*mp.pi*m)**2, [1, mp.inf])
rhs_h = (fh(0, e) + 2*mp.nsum(lambda m: (-1)**int(m)*fh(2*mp.pi*m, e), [1, mp.inf]))/2
lhs_h = mp.nsum(lambda n: (n + mp.mpf(1)/2)*mp.e**(-e*(n + mp.mpf(1)/2)), [0, mp.inf])
ok("3. half-integer modes: sum_{n>=0}(n + 1/2) e^{-eps(n + 1/2)} - 1/eps^2 = %s, %s, %s: +1/24 = %s; it is -2 sum_m (-1)^m/(2 pi m)^2 = %s; "
   "Poisson with the factor (-1)^m on the m-th image holds at eps = 0.001 (relative %.0e)"
   % (mp.nstr(half[0], 10), mp.nstr(half[1], 10), mp.nstr(half[2], 12), mp.nstr(mp.mpf(1)/24, 12), mp.nstr(alt, 12), abs(lhs_h/rhs_h - 1)),
   abs(half[2] - mp.mpf(1)/24) < 1e-6 and abs(alt - mp.mpf(1)/24) < 1e-30 and abs(lhs_h/rhs_h - 1) < 1e-25)
# 4
e = mp.mpf('0.01'); f = lambda n: n*mp.e**(-e*n); fp = lambda n: (1 - e*n)*mp.e**(-e*n)
vals = []
for n0 in (5, 20, 100):
    head = mp.fsum(f(n) for n in range(1, n0)); tail = mp.quad(f, [n0, mp.inf])
    vals.append(head + tail + f(n0)/2 - fp(n0)/12 - 1/e**2)
ok("4. eps = 0.01: sum to n0, the tail as an integral, plus f(n0)/2 - f'(n0)/12, minus 1/eps^2: %s, %s, %s at n0 = 5, 20, 100, against -1/12 = %s "
   "(what is left is the next Euler-Maclaurin term, f'''(n0)/720, of order eps^2)"
   % (mp.nstr(vals[0], 10), mp.nstr(vals[1], 10), mp.nstr(vals[2], 10), mp.nstr(-mp.mpf(1)/12, 10)), all(abs(v + mp.mpf(1)/12) < 1e-5 for v in vals))
print("ALL PASS" if not fails else "%d FAILED" % fails); sys.exit(1 if fails else 0)
