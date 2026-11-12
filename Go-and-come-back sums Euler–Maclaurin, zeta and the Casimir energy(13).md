# Go-and-come-back sums: Euler–Maclaurin, zeta and the Casimir energy

Oct 12, 2026 · @few4324

The go-to-infinity sum formula, rewritten with E = e^D = 1 + Δ, gives a Taylor series for S(x) whose coefficients are regularized sums of derivatives; a second route through decaying exponentials gives an exact integral. The same machinery gives zeta as a Taylor coefficient of the fractional sum, and the Casimir energy in one space dimension as a convergent sum over round trips.

## Way 1: derivatives

Factoring e^{xD} − 1 out of the whole bracket leaves something independent of x, and expanding that factor gives the Taylor series of S at 0.

```latex
\begin{aligned}
S(x)
&=\sum_{k=1}^{x} f(k)\\[6pt]
&=\lim_{n\to\infty}\left[\sum_{k=1}^{n-1}\big(f(k)-f(x+k)\big)+\sum_{k=0}^{x-1}f(n+k)\right]\\[6pt]
&=\lim_{n\to\infty}\left.\left[\sum_{k=1}^{n-1}\big(1-e^{xD}\big)e^{kD}+\frac{e^{xD}-1}{e^{D}-1}\,e^{nD}\right]f(t)\right|_{t=0}\\[6pt]
&=\lim_{n\to\infty}\left.\big(e^{xD}-1\big)\left[\frac{e^{nD}}{e^{D}-1}-\sum_{k=1}^{n-1}e^{kD}\right]f(t)\right|_{t=0}\\[6pt]
&=\lim_{n\to\infty}\left.\sum_{m=1}^{\infty}\frac{x^{m}D^{m}}{m!}\left[\sum_{j=0}^{\infty}\frac{B_j}{j!}D^{j-1}e^{nD}-\sum_{k=1}^{n-1}e^{kD}\right]f(t)\right|_{t=0}\\[6pt]
&=\sum_{m=1}^{\infty}\frac{x^{m}}{m!}\,\lim_{n\to\infty}\left[\sum_{j=0}^{\infty}\frac{B_j}{j!}\,f^{(m+j-1)}(n)-\sum_{k=1}^{n-1}f^{(m)}(k)\right]
\end{aligned}
```

Here B₁ = −½. Each Taylor coefficient is minus the sum of f^(m)(k) over all k, with the Bernoulli terms at n supplying what that sum needs to converge:

```latex
S^{(m)}(0)=\lim_{n\to\infty}\left[f^{(m-1)}(n)-\tfrac12 f^{(m)}(n)+\tfrac1{12}f^{(m+1)}(n)-\cdots-\sum_{k=1}^{n-1}f^{(m)}(k)\right]
```

The sum over j is asymptotic, not convergent, and the equalities above are to be read with it truncated at some J: the bracket then differs from S^(m)(0) by a remainder that goes to 0 as n → ∞ for fixed J (the omitted terms are O(f^(m+J)(n))), which is what the limit uses. It is only needed for the first few m; once f^(m) is summable the bracket is the plain sum. Stated as a definition rather than an identity: S^(m)(0) is the limit of the finite bracket, truncated at J, for the class of f whose derivatives f^(m+J)(n) → 0 — a finite expansion plus a remainder, not an infinite Bernoulli series substituted as an operator. The remainder bound itself is the open point listed at the end.

**Examples.** For f = ln the m = 1 bracket tends to −γ, and for f = t^a every coefficient is a zeta value:

```latex
\ln(x!)=-\gamma x+\sum_{m=2}^{\infty}\frac{(-1)^{m}\zeta(m)}{m}\,x^{m},
\qquad
\sum_{k=1}^{x}k^{a}=-\sum_{m=1}^{\infty}\binom{a}{m}\,\zeta(m-a)\,x^{m}
```

Both agree with the exact values at x = ½ to 25 digits (the second with a = ½ gives 0.4383320906).

## The same operator, without the Taylor step

Leaving the first sum alone and expanding only the jump gives closed forms in x that hold for any x, not just |x| < 1.

```latex
\begin{aligned}
S(x)
&=\lim_{n\to\infty}\left[\sum_{k=1}^{n-1}\big(f(k)-f(x+k)\big)+\left.\frac{E^{x}-1}{E-1}\,f(t)\right|_{t=n}\right]\\[6pt]
&=\lim_{n\to\infty}\left[\sum_{k=1}^{n-1}\big(f(k)-f(x+k)\big)+\sum_{j=1}^{\infty}\binom{x}{j}\left.\Delta^{j-1}f(t)\right|_{t=n}\right]\\[6pt]
&=\lim_{n\to\infty}\left[\sum_{k=1}^{n-1}\big(f(k)-f(x+k)\big)+\sum_{m=1}^{\infty}\frac{B_m(x)-B_m}{m!}\left.D^{m-1}f(t)\right|_{t=n}\right]\\[6pt]
&=\lim_{n\to\infty}\left[\sum_{k=1}^{n-1}\big(f(k)-f(x+k)\big)+\int_{n}^{n+x}f(t)\,dt+\sum_{m=1}^{\infty}\frac{B_m}{m!}\left.\big(e^{xD}-1\big)D^{m-1}f(t)\right|_{t=n}\right]
\end{aligned}
```

Line 2 is the video's formula (E = 1 + Δ). Lines 3 and 4 use E = e^D and are asymptotic in m. Line 4 is the textbook Euler–Maclaurin formula for the jump, from the factorization (e^{xD} − 1)/D · D/(e^D − 1).

| Order | Newton coefficient | Bernoulli-polynomial coefficient |
| --- | --- | --- |
| 1 | x | x |
| 2 | x(x−1)/2 | x(x−1)/2 |
| 3 | x(x−1)(x−2)/6 | x(x−1)(2x−1)/12 |
| 4 | x(x−1)(x−2)(x−3)/24 | x²(x−1)²/24 |

Summing the first geometric series in closed form as well cancels n and returns Euler–Maclaurin applied at t = 1, so that one sum has to stay discrete.

## Way 2: decaying exponentials

On a mode e^{−st} the shift is E = e^{−s} with |E| < 1, so every series converges and the result is an exact integral. Take f as a superposition of modes, minus the first p Newton terms of each mode about t = 1:

```latex
f(t)=\int_0^\infty\varphi(s)\left[e^{-st}-e^{-s}\sum_{i=0}^{p-1}\binom{t-1}{i}\big(e^{-s}-1\big)^{i}\right]ds
```

The subtracted part is a polynomial of degree p − 1 in t, and its sum is known exactly. Then:

```latex
\begin{aligned}
S(x)
&=\lim_{n\to\infty}\left[\sum_{k=1}^{n-1}\big(f(k)-f(x+k)\big)+\sum_{j=1}^{p}\binom{x}{j}\left.\Delta^{j-1}f(t)\right|_{t=n}\right]\\[6pt]
&=\lim_{n\to\infty}\int_0^\infty\varphi(s)\left[\sum_{k=1}^{n-1}\big(e^{-sk}-e^{-s(x+k)}\big)+\big(e^{-sn}-e^{-s}\big)\sum_{j=1}^{p}\binom{x}{j}\big(e^{-s}-1\big)^{j-1}\right]ds\\[6pt]
&=\lim_{n\to\infty}\int_0^\infty\varphi(s)\,\big(e^{-s}-e^{-sn}\big)\left[\frac{1-e^{-sx}}{1-e^{-s}}-\sum_{j=1}^{p}\binom{x}{j}\big(e^{-s}-1\big)^{j-1}\right]ds\\[6pt]
&=\lim_{n\to\infty}\int_0^\infty\varphi(s)\,\big(e^{-s}-e^{-sn}\big)\sum_{j=p+1}^{\infty}\binom{x}{j}\big(e^{-s}-1\big)^{j-1}\,ds\\[6pt]
&=\int_0^\infty\varphi(s)\,e^{-s}\sum_{j=p+1}^{\infty}\binom{x}{j}\big(e^{-s}-1\big)^{j-1}\,ds\\[6pt]
&=\int_0^\infty\varphi(s)\left[\frac{1-e^{-sx}}{e^{s}-1}-e^{-s}\sum_{j=1}^{p}\binom{x}{j}\big(e^{-s}-1\big)^{j-1}\right]ds
\end{aligned}
```

Line 4 holds at finite n, so the remainder is exact:

```latex
\Big[\ \cdot\ \Big]_{n}-S(x)=-\int_0^\infty\varphi(s)\,e^{-sn}\sum_{j=p+1}^{\infty}\binom{x}{j}\big(e^{-s}-1\big)^{j-1}\,ds=-\sum_{j=p+1}^{\infty}\binom{x}{j}\left.\Delta^{j-1}f(t)\right|_{t=n}
```

| p | φ(s) | f(t) | S(½) from the integral |
| --- | --- | --- | --- |
| 1 | −1/s | ln t | −0.120782237635 |
| 2 | 1/s² | t·ln t − 2·ln 2·(t − 1) | 0.045963787887 |
| 3 | −2/s³ | t²·ln t − quadratic | −0.044824227104 |

All three agree with the exact values to 30 digits, and line 4 at n = 7 for f = ln agrees with the direct bracket to 25 digits. For f = 1/t and f = ln the final line is the classical integral for the harmonic numbers and for ln Γ.

## Against plain Euler–Maclaurin

Both forms of Way 1 beat Euler–Maclaurin applied at the start of the sum by 10 to 12 orders of magnitude, and they match each other. That comparison is against the naive use; the standard practice — sum the first n terms directly and apply Euler–Maclaurin to the tail — is line 4 of the section above, and it reaches the same errors as the Taylor form. The gain is the shifted endpoint, nothing outside Euler–Maclaurin: the method is Euler–Maclaurin applied at n. What Way 1 adds is not accuracy over the shifted formula but the fractional sum itself: a value for non-integer x, from the same bracket.

Errors at x = ½. "Plain" is Euler–Maclaurin with endpoints 1 and 1 + x; the other two columns use n = 10 with 12 Bernoulli terms.

| f | Plain, best possible | Plain, 30 terms | Taylor form | Bernoulli-polynomial form |
| --- | --- | --- | --- | --- |
| ln t | 2.6 × 10⁻⁴ | 6.6 × 10⁵ | 2.9 × 10⁻¹⁶ | 1.2 × 10⁻¹⁵ |
| 1/t | 1.7 × 10⁻³ | 1.9 × 10⁷ | 3.9 × 10⁻¹⁶ | 1.6 × 10⁻¹⁵ |
| √t | 3.1 × 10⁻⁵ | 3.5 × 10⁴ | 7.1 × 10⁻¹⁷ | 3.1 × 10⁻¹⁶ |
| t^(3/2) | 6.9 × 10⁻⁶ | 1.9 × 10³ | 8.8 × 10⁻¹⁷ | 4.1 × 10⁻¹⁶ |

- **The gain comes from n.** With n = 1 the Taylor form hits the same floor as plain Euler–Maclaurin (2.6 × 10⁻⁴ for ln). With 6 Bernoulli terms it reaches 3.5 × 10⁻⁹ at n = 5 and 1.7 × 10⁻¹¹ at n = 10.
- **The Taylor form needs |x| < 1.** For ln at x = 0.9 it reaches only 1.3 × 10⁻⁶ after 80 powers of x, and at x = 2.5 it diverges. The Bernoulli-polynomial form gives 2.0 × 10⁻¹⁵ at x = 0.9.
- **For small x the Taylor form is the efficient one.** For ln at x = 0.1 it reaches 8.3 × 10⁻¹³ with 10 powers.
- **Beyond |x| < 1,** use the recursion S(x) = S(x − 1) + f(x) to bring x into (−½, ½\] first.

## Zeta from the fractional sum

Zeta is the first Taylor coefficient of the fractional sum of a power, and the video's formula and the derivative formula are one operator, D/Δ, expanded in Δ and in D.

**The video's method.** Differentiate the fractional sum of k^(1−s) in its upper limit z and set z = 0, using d/dz of C(z, k) at 0 equal to (−1)^(k−1)/k:

```latex
\begin{aligned}
\zeta(s)
&=\frac{1}{s-1}\left.\frac{d}{dz}\sum_{k=1}^{z}k^{1-s}\right|_{z=0}\\[6pt]
&=\frac{1}{s-1}\left.\frac{d}{dz}\lim_{n\to\infty}\left[\sum_{k=1}^{n-1}\big(k^{1-s}-(z+k)^{1-s}\big)+\sum_{k=1}^{m}\binom{z}{k}\Delta^{k-1}n^{1-s}\right]\right|_{z=0}\\[6pt]
&=\lim_{n\to\infty}\left[\sum_{k=1}^{n-1}\frac{1}{k^{s}}+\frac{1}{s-1}\sum_{k=1}^{m}\frac{(-1)^{k-1}}{k}\,\Delta^{k-1}n^{1-s}\right]
\end{aligned}
```

It holds for Re s > 1 − m, so raising m continues zeta strip by strip.

**The derivative version.** The same step applied to the operator gives both forms at once:

```latex
\left.\frac{d}{dx}\frac{E^{x}-1}{E-1}\right|_{x=0}=\frac{\ln E}{E-1}=\frac{D}{\Delta}
=\sum_{k=1}^{\infty}\frac{(-1)^{k-1}}{k}\Delta^{k-1}
=\sum_{j=0}^{\infty}\frac{B_j}{j!}D^{j}
```

```latex
\begin{aligned}
\zeta(s)
&=\lim_{n\to\infty}\left[\sum_{k=1}^{n-1}\frac{1}{k^{s}}+\frac{1}{s-1}\sum_{j=0}^{\infty}\frac{B_j}{j!}\,D^{j}n^{1-s}\right]\\[6pt]
&=\lim_{n\to\infty}\left[\sum_{k=1}^{n-1}\frac{1}{k^{s}}+\frac{n^{1-s}}{s-1}+\frac{1}{2n^{s}}+\frac{s}{12\,n^{s+1}}-\frac{s(s+1)(s+2)}{720\,n^{s+3}}+\cdots\right]
\end{aligned}
```

This is the m = 1 coefficient of the power series in Way 1; the higher coefficients give zeta at s + 1, s + 2, and so on.

**Same zeta, different speed.** Errors at n = 10 with the same number of correction terms:

| Correction terms | ζ(½), video | ζ(½), derivative | ζ(−½), video | ζ(−½), derivative |
| --- | --- | --- | --- | --- |
| 1 | 1.6 × 10⁻¹ | 1.6 × 10⁻¹ | 1.6 | 1.6 |
| 3 | 5.0 × 10⁻⁴ | 8.2 × 10⁻⁷ | 1.8 × 10⁻³ | 1.6 × 10⁻⁶ |
| 5 | 1.8 × 10⁻⁵ | 3.1 × 10⁻⁹ | 3.0 × 10⁻⁵ | 3.4 × 10⁻⁹ |
| 9 | 2.1 × 10⁻⁷ | 4.3 × 10⁻¹³ | 1.9 × 10⁻⁷ | 2.6 × 10⁻¹³ |

**Which corrections survive the limit.** Each term is one power of n smaller than the last, so the boundary is Re s = 1. The count below is for the consecutive Newton corrections of the video's form; in the Bernoulli derivative form the odd coefficients beyond B₁ vanish, so the number of surviving terms differs there:

| Range of s | Plain sum of 1/k^s | Corrections that do not vanish |
| --- | --- | --- |
| Re s > 1 | converges | none |
| 0 < Re s < 1 | diverges | the first |
| −1 < Re s < 0 | diverges | the first two |
| −2 < Re s < −1 | diverges | the first three |

Below Re s = 1 the sum and the surviving corrections each diverge, and zeta is the finite part left when they cancel. At s = 1 the first correction divides by zero, which is the pole.

**Negative integers.** At s = −q the function n^(q+1) is a polynomial, so both series stop after q + 2 terms and are exact at any n, including n = 1. Checked for s = −1, −2, −3 to 30 digits. For s = −1 at n = 1 the first sum is empty:

```latex
\zeta(-1)=\frac{1}{s-1}\Big[n^{2}-\tfrac12\cdot 2n+\tfrac1{12}\cdot 2\Big]_{n=1}=-\tfrac12\Big[1-1+\tfrac16\Big]=-\tfrac1{12}
```

In general ζ(−q) = −B\_(q+1)/(q+1) for q ≥ 1, with B₁ = −½ as above; at q = 0 the same bracket gives ζ(0) = −½ = B₁ (the formula with the minus sign would give +½ there, which is the usual B₁-convention trap). Either way it is the linear coefficient of the polynomial that sums k^(q+1), not a value of the divergent series — precisely, S′\_(q+1)(0) = −(q+1) ζ(−q), where S\_(q+1)(x) = Σ\_(k≤x) k^(q+1); the product of the vanishing binomial coefficient and the pole at s = 1 is taken as a limit, not set to zero first.

## Zeta defined by the fractional sum

The fractional sum of k^a is finite for every a, so its Taylor coefficients can be taken as the definition of zeta, with no analytic continuation in the definition. Fractional sums in this sense are Müller and Schleicher's ("How to add a non-integer number of terms, and how to produce unusual infinite summations", J. Comput. Appl. Math. 178, 347 (2005); "Fractional sums and Euler-like identities", Ramanujan J. 21, 123 (2010)), and the "master equation" below is the Taylor expansion of the Hurwitz zeta function, Σ\_{k=1}^{x} k^a = ζ(−a) − ζ(−a, 1+x): a theorem, not something to be checked at one value.

```latex
\sum_{k=1}^{x}k^{a}=-\sum_{m=1}^{\infty}\binom{a}{m}\,\zeta(m-a)\,x^{m}
\qquad\Longrightarrow\qquad
\zeta(m-a)=-\frac{1}{m!\binom{a}{m}}\left.\frac{d^{m}}{dx^{m}}\sum_{k=1}^{x}k^{a}\right|_{x=0}
```

With m = 1 and a = 1 − s this is the video's formula. Written out, it is one explicit limit for every s ≠ 1:

```latex
\zeta(s)=\lim_{n\to\infty}\left[\sum_{k=1}^{n-1}\frac{1}{k^{s}}+\frac{1}{s-1}\sum_{j=0}^{m}\frac{B_j}{j!}\,D^{j}n^{1-s}\right],\qquad m>1-\operatorname{Re}s
```

**Why it is legitimate.** For a < 0 nothing diverges and the identity is a Taylor expansion in which every zeta is an ordinary convergent sum:

```latex
\sum_{k=1}^{x}k^{a}=\sum_{k=1}^{\infty}\big(k^{a}-(x+k)^{a}\big)=-\sum_{k=1}^{\infty}\sum_{m=1}^{\infty}\binom{a}{m}x^{m}k^{a-m}=-\sum_{m=1}^{\infty}\binom{a}{m}\zeta(m-a)\,x^{m}
```

For a ≥ 0 the left side still exists because of the correction terms, and the right side has only one way to follow it.

**Self-consistency.** The same zeta value appears as different coefficients for different a, and they agree because

```latex
\frac{d}{dx}\sum_{k=1}^{x}k^{a}=-a\,\zeta(1-a)+a\sum_{k=1}^{x}k^{a-1}
```

**What this framing does and does not give.**

- **It explains the values.** A value such as −1/12 is the finite part left after the divergent sum and the surviving corrections cancel, and the formula shows which pieces those are.
- **It is known.** The explicit limit is the standard Euler–Maclaurin definition of zeta, used to compute it. What is found here is a route to it through the fractional sum.
- **Continuation becomes a theorem, not a definition.** Showing this function is the usual zeta still needs the uniqueness argument once, and needs the fractional sum to be analytic in a. One thing the recurrence does not do: S(x) − S(x−1) = f(x) with S(0) = 0 is also obeyed by S(x) + b sin(2πx), which agrees with S at every integer and changes S'(0) by 2πb. So the fractional sum — and with it the zeta coefficient — is fixed by the limit definition above (the finite-order corrected bracket, on its stated class of f), not by the recurrence and the integer values alone; that definition is what the uniqueness argument refers to.
- **s = 1 stays excluded.** There C(a, m) = 0 exactly where zeta has its pole; for a = 1, m = 2 the product has the finite limit that gives the ½x² in x(x+1)/2. That example, f(t) = t with S(x) = x(x+1)/2, is worth keeping in view: it is entire, so there is no universal radius of convergence of one for the Taylor form; the radius is set by the singularities of S, which for a polynomial do not exist and for 1/t sit at x = −1.
- **It is one definition among several.** The alternating series for Re s > 0 and the integral with 1/(e^t − 1) from Way 2 are others, and the functional equation comes more naturally from the integral forms.

The identity was checked numerically at a = ½ here; as the Hurwitz-zeta expansion it holds for all a ≠ −1 (the pole), so the step from a < 0 to all a is covered by that theorem rather than by the argument above.

## The path taken (read this first)

The Casimir part below contains three calculations, and until 1 November 2026 it did not say which one is the physics. It is this one.

**Two carriers, kept apart.** The coin's calculation below is the Dirac walker's own. III.14's kernel is P(x) = e^{−ax/c} δ(t − x/c) + (a/2) e^{−ax/c}\[I₀(ζ) + …\]: a term that never flips and goes at c, and the zigzags. Between two walls the no-flip leg crosses, is turned by the wall (the turn matrix, T per hit), crosses back, and after m round trips has length 2mL. That path is the loop of a *fermion* Casimir effect — the energy of the walker's own field between walls that turn it — and the massless calculations below are this with a → 0: the no-flip weight e^{−2maL/c} is 1 and the zigzag terms, which carry a factor a, are absent. It is not the force measured between two metal plates. That one is electromagnetic: what crosses the gap is the carrier of the main note's VI.2 (E), Part II's walker at m = 0 with a one-form index, and the walls are not turn matrices but plates of charged walkers, whose response to the carrier is what reflects it. "The plates, made of walkers" below does that calculation. The coin's sections stay as the fermion example and as the place where the loop machinery is shown to work. Three constructions, then: the *coin's* (alternating round trips, (r\_L r₀)^m = (−1)^m for the turn-matrix wall, the no-flip kernel), *Part II's scalar walker's* (Dirichlet images, no alternation, the −½ dτ/τ loop weight), and the *plates'* (the carrier reflected by what the plates' walkers do).

**The same −π/24L twice is not a coincidence.** The scalar's ½ Σ 1/m² = π²/12 and the coin's −Σ(−1)^m/m² = π²/12 agree because both fields have central charge c = 1: a Dirac fermion in 1+1 is one scalar (bosonization), and the energy of one chiral half of a conformal field on a circle of length ℓ is E = −(2πħc/ℓ)(c/24 − h), with h the weight of the sector (the whole field on a circle has both halves, c/12). A strip of width L between walls unfolds to one chiral half on ℓ = 2L, so E = (πħc/L)(h − c/24), and the coin's whole wall table is one line,

```latex
\frac{E\,L}{\pi\hbar c}=h-\frac1{24},\qquad h=\frac12\Big(\frac\varphi{2\pi}-\frac12\Big)^2,
```

which equals the round-trip sum (ħc/2πL) Σ\_{n≥1} cos(nφ)/n² = (πħc/2L) B₂(φ/2π) identically in φ. Row by row: h = 0 for the turn-matrix wall (φ = π, −π/24L — the untwisted sector, the scalar's Dirichlet number), h = 1/8 for the half-turn wall (φ = 0, +π/12L), h = 1/32 for the ±i wall (φ = π/2, −π/96L). The wall's phase is a twist of the boundary condition and the table is the weights of the twisted sectors; the smooth-cutoff mode sums of the scalar and of the coin's sea both land on −π/24L for the same reason. (`scripts/wall_twist_is_central_charge.py`, 1–4.)

**The displayed bridge: the coin's number from its own kernel.** Unfold the two movers into one. With walls at 0 and L, ψ₋(0) = r₀ψ₊(0) and ψ₊(L) = r\_Lψ₋(L); put φ(s) = ψ₊(s) on \[0, L) and φ(s) = r\_Lψ₋(2L − s) on \[L, 2L). Then φ is one right-mover on a circle of length 2L with φ(s + 2L) = r\_Lr₀φ(s), and r\_Lr₀ is the wall matrix applied twice: −1 for the turn matrix, +1 for the half-turn; the modes e^{2ikL} = r\_Lr₀ are the two rows of the table, (n+½)π/L and nπ/L. The no-flip kernel's two-point function is the chiral one, G₀(u) = 1/(2πiu), u = s − t, normalized by {ψ, ψ†} = δ — the walker's own normalization, not a choice. With the walls by images, G(u) = Σₘ (r\_Lr₀)^m G₀(u + 2mL) = (π/2L)/(2πi·sin(πu/2L)) for −1 and (π/2L)/(2πi·tan(πu/2L)) for +1. The square rule at the wall is the energy bilinear of a right-mover, −iψ†∂ₛψ, point-split on G: the m = 0 image is the bulk, and each other image gives exactly (r\_Lr₀)^m/(2π(2mL)²) — this is the 1/ℓ² of the round trip, with its normalization from G₀, and nothing from "line 2" of the coin's section is used. Summed: ⟨T⟩ = −π/(48L²) for the turn matrix and +π/(24L²) for the half-turn; over the circle of length 2L, E = −π/24L and +π/12L: the table's −0.130900 and +0.261799, from the kernel. (`scripts/coin_casimir_from_the_no_flip_leg.py`, 1–4.) The no-flip survival factor e^{−2maL/c} is only the mass factor; the 1/ℓ² and its normalization are the kernel's, as just shown.

**The three calculations are one.** (i) *Loops*: the sum above, now derived for the coin — the physics. For Part II's scalar walker the same step is the point-split kernel with Dirichlet images read with the square rule ("What the walk now supplies", below): ⟨T\_xx⟩ at the wall = −π/(24L²) = −dE/dL. (ii) *Modes*: Σ½ħωₙ with its finite part; its Poisson partner, image by image ("Modes and loops are Poisson partners"): the m-th image of the mode sum is the m-fold round trip, with (−1)^m for the half-integer modes of the turn-matrix wall. (iii) *Fractional sums*: on the loop side, x = δ/2L is the round trip that lands δ short of closing — the point-splitting made physical — and the Newton/Euler–Maclaurin expansion is the expansion in δ, with ζ(4), ζ(6) the response to pulling the two ends apart ("Loops that miss closing"); on the mode side, x is a fractional mode number, which means nothing, and the fractional sum there is a device that defines the finite part, validated by (ii). The extension of a sum by a small amount is a ballistic leg's shortfall, not a zigzag: the Newton expansion is binomial in the shift, the kernel's zigzags are Bessel in the flip count (main note III.15), two different smallnesses.

**What the force is.** On the route of (i) the wall feels the vacuum's momentum flux through it, a bilinear in the readings at the wall — what a wall made of walkers registers (main note V.5). Of V.5's three quantities the flux is the one a wall reads as a push, and it is what Lifshitz's formula computes for the plates below. V.5 (C)'s other reading, the shift of the wall's clock rate by the loops that touch it, with the force as the gradient of that shift, is the detector model; in 1+1 massless the two are one number (T\_xx = T₀₀), and for the plates the force is the flux.

**Inputs, after all this.** The wall factor (what the wall does to the leg), and the vacuum as the source-free kernel of III.16 with images (equivalently, on the mode side, the state that stays the same when the walker's added detectors are counted with the other sign). For the coin the main note now reduces the second to one thing taken: that the vacuum is the configuration the long-clock count of exclusive configurations selects (main note VI.4 (C), VI.5; 9 November 2026). Not inputs any more: the 1/τ and the ½ of the loop weight (counting), "the vacuum has an energy", "the vacuum is the lowest state". For the plates the inputs are two, stated in their section: the *material* (which walkers, how many, how confined) and the *state* (the lowest one: every column below E\_F filled, the sea included). Not inputs there: the wall factor, which is derived from the plates' walkers, and the field's own fluctuations, which the plates' own lines supply. **Open:** the massive case. The no-flip weight gives e^{−2maL/c} per round trip ("With mass", below), but a massive walker can also close a loop by flipping between the walls — the I₀, I₁ part of the kernel — and those closings are not computed here; the ballistic piece alone is the large-mL approximation. And the plates' far field: the carrier is used static there (the near field); retardation needs the moving end-to-end rule.

## Casimir energy as a sum over closed worldlines

The energy between the walls is a sum over round trips: in Part II's walk, closed worldlines that bounce off both walls; in the coin's version, paths that return to the same state after 2nL/c. One particle alone gives no attraction. (−½ Tr|H| is the energy of the filled Dirac sea written term by term; what the round-trip sum adds is the bookkeeping, not a different quantity.)

**One particle.** A single coin particle between walls has energy ħω with ω = (n + ½)πc/L. Pushing the walls in raises it, so it presses outward like any particle in a box.

**What the two sums count.** In Part II's walk a loop is a closed worldline, and by III.3 a worldline that returns to its starting point in spacetime runs backward in time somewhere, so the forward-time observer reads it as a pair: created, travelled, annihilated. In the coin's version nothing closes in spacetime — the walker always advances in t — and a "round trip" is a return to the same state, position and direction, after 2nL/c; the wall is the turn matrix applied once per hit, and two hits per round trip give T² = −1. That −1 is the whole sign of the coin's sum. It is not a time-sense flip at the wall, and no closed worldline with time flips at the walls exists here: the time-sense walker of III.1 never moves left, and a leg going left with Δt < 0 lies on the same null line as the leg before it and has ε = +1 (an earlier draft of this section drew such a loop; it was wrong). The two sums are different objects with different phases — +1 per loop in Part II's, −1 per round trip in the coin's — and give the same −πħc/24L; the coin's phase is an input there (the choice of wall), as the table below shows.

&#91;embedded content: one round trip between the walls, massless and with mass\]

With mass the legs are no longer straight: the time sense can also flip between the walls (III.1), each leg is III.14's zigzag, and the two members of the pair wander apart inside the wedge outside the light cone. A loop that must reach both walls then pays III.14's survival weight over its length, e^(−2mL), which is the K₁ of the massive formula below.

**The loops of Part II's walk.** Take the walk with τ summed (III.16), put in the walls by images, and count each loop with weight 1/τ. Loops that wind n times between the walls are displaced by 2nL (ħ = c = 1):

```latex
\begin{aligned}
E(L)&=-\tfrac12\,L\sum_{n\neq0}\int_0^\infty\frac{d\tau}{\tau}\,\frac{e^{-(2nL)^{2}/4\tau}}{4\pi\tau}\\[6pt]
&=-\tfrac12\,L\sum_{n\neq0}\frac{1}{4\pi\,n^{2}L^{2}}\\[6pt]
&=-\frac{1}{4\pi L}\,\zeta(2)=-\frac{\pi\hbar c}{24L}
\end{aligned}
```

The n = 0 loops touch no wall; they are proportional to L and are there without walls. Paths with an odd number of wall hits give a term that does not depend on L.

**With mass.** A factor e^(−m²τ) in the same integral gives

```latex
\int_0^\infty\frac{d\tau}{\tau}\,\frac{e^{-(2nL)^{2}/4\tau-m^{2}\tau}}{4\pi\tau}=\frac{m\,K_1(2m|n|L)}{2\pi\,|n|L},\qquad
E(L)=-\frac{m}{2\pi}\sum_{n=1}^{\infty}\frac{K_1(2nmL)}{n}\;\xrightarrow{\ m\to0\ }\;-\frac{\pi}{24L}
```

The loops die off as e^(−2mL): a pair that has to live long enough to reach both walls and return costs 2mL. For m = 1 the difference E(0.5) − E(1.5) is −0.103232 from this formula, and the smooth-cutoff mode sum gives −0.10307, −0.10319, −0.10322 at cutoffs 0.02, 0.01, 0.005, converging onto it. For a particle with the electron's mass the energy is negligible beyond about 10⁻¹² m; the long-range force needs a carrier that never flips.

**Inputs.**

- [ ] The loop weight, −½ and 1/τ, is the worldline form of "the vacuum has an energy". It is not derived from III.1 — but see "What the walk now supplies" below: the 1/τ and the ½ are counting, and the route through the kernel and the square rule needs neither "the vacuum has an energy" nor "the vacuum is the lowest state".
- [ ] These are loops of Part II's walk (Klein–Gordon). The pairs of III.3 belong to the coin (Dirac); the coin's own calculation is the next section and gives the same number for one choice of wall.
- [ ] Walls made of walkers whose fluctuations act on each other (the Lifshitz picture) is the calculation of "The plates, made of walkers" below.

**What the walk now supplies of these inputs (1 November 2026).** Two of the inputs above were "the vacuum has an energy", in two dresses: the loop weight −½ dτ/τ of the scalar calculation and E\_vac = −½ Tr|H| of the coin's. With Parts IV–VI of the main note in place, most of each is counting, and the rest can be replaced by a route that uses neither.

- *The 1/τ is counting.* ln det(1 − P) = −Σₙ Tr Pⁿ/n for the walk's step matrix P, and Tr Pⁿ/n is one count per closed loop divided by its symmetry: a loop of length n has n starting points. Checked on the ring (`casimir_from_the_walk.py` 1a, 1b).
- *The ½ is the undirected loop.* Part II's walker has no direction (VI.1 (F)): a loop and its reverse are one loop, and the directed trace counts both. For the coin the ½ is the symmetrization of the next item, not a return leg: III.8c's returns are the measurement's pairs, and an event's weight is already the product of two readings, not a half of anything.
- *The coin's |H| and the ½* are the sea and its symmetrization, and here one prescription stays. Exclusion in every state (VI.1 (D)) says at most one walker per state; it does not say which states the vacuum occupies — that is a choice of state. The choice the mode form uses is the symmetric one: when the walker's added detectors are counted with the other sign, J → −J, which sends e → −e and leaves the free generator (V.2 (H)), occupied and empty trade places, and the state symmetric under it has E₀ = ½(Σ\_{E<0}E − Σ\_{E>0}E) = −½Σ|E|. So: the ½ and the |E| follow from exclusion plus that symmetry; that the vacuum is the symmetric state is the prescription, the mode-side form of "the vacuum is the source-free kernel" below. The sign of a closed fermion loop relative to a scalar loop is one exchange, −1 (VI.1 (D)). The main note's VI.4–VI.5 (9 November 2026) put the prescription in one place: the vacuum is the configuration the long-clock count of exclusive configurations selects, and that configuration is the symmetric one; around it a closed history that winds n times carries (−1)^(n−1).
- *What was left*: that the loop sum is a free energy, and that the vacuum is the lowest state. Neither is needed, but one thing is stated in their place: which two-point function is the vacuum's. The energy and momentum densities are III.15a's square read at coincident points of the vacuum's two-point function, and the two-point function is taken to be the source-free kernel of III.16 with the walls by images — "no walker", read with III.16, in III.16's normalization (the source answer with unit source; the scale is the walk's ħ, input 17, not fixed here). That is the vacuum prescription on this route, stated once; deriving the kernel's equation does not choose it. For the massless scalar between walls at 0 and L: on the direct images (∂\_t∂\_t′ + ∂\_x∂\_x′)G₀ = −1/(πu²) with u = 2nL, on the reflected images exactly 0, so ⟨T₀₀⟩ = −π/(24L²), uniform, E(L) = −π/(24L) — the loop column and the mode column of the tables above, from the kernel and the square rule (`casimir_from_the_walk.py` 3). The force on the wall is the vacuum's momentum flux through it, ⟨T\_xx⟩ at the wall = −π/(24L²) = −dE/dL (ibid. 4): a bilinear in the readings at the wall, IV.4d (B)'s conserved pair, the thing a wall made of walkers registers (V.5). No sea is filled and nothing is minimized; the count of loops that reach both walls changes with L, and that change is the force.

What remains an input on this route: the wall factor (the turn matrix at the wall, as stated above) and the two-point function being the source-free kernel — the vacuum as "no walker", read with III.16. The mode sums with a smooth cutoff reproduce the table's −π/24 (turn-matrix wall) and +π/12 (half-turn wall) (ibid. 2), so the three routes — loops, modes, and the square on the kernel — are one number. (`scripts/casimir_from_the_walk.py`, 9 checks.)

## The coin's version: closed paths and the time shift

With energy as the generator of the time shift, a closed path of duration t carries energy ħ/(2πt²) times its wall factor, and summing over round trips gives −πħc/24L with no field in the derivation.

```latex
\begin{aligned}
U(t)&=e^{tD_t}=e^{-iHt/\hbar}\\[6pt]
|H|&=\frac{\hbar}{\pi}\int_{-\infty}^{\infty}\frac{1-U(t)}{t^{2}}\,dt\\[6pt]
E&=-\tfrac12\,\mathrm{Tr}\,|H|=\frac{\hbar}{2\pi}\int_{t\neq0}\frac{\mathrm{Tr}\,U(t)}{t^{2}}\,dt+\text{(bulk)}\\[6pt]
\mathrm{Tr}\,U(t)&=\frac{2L}{c}\sum_{n}e^{in\varphi}\,\delta\!\Big(t-\frac{2nL}{c}\Big)\\[6pt]
E(L)&=\frac{\hbar}{2\pi}\cdot\frac{2L}{c}\sum_{n\neq0}\frac{e^{in\varphi}}{(2nL/c)^{2}}=\frac{\hbar c}{2\pi L}\sum_{n=1}^{\infty}\frac{\cos n\varphi}{n^{2}}=\frac{\pi\hbar c}{2L}\,B_2\!\Big(\frac{\varphi}{2\pi}\Big)
\end{aligned}
```

| Line | What it is |
| --- | --- |
| 2 | The identity \|ω\| = (1/π)∫(1 − cos ωt)/t² dt applied to the shift: the size of the energy is a weighted sum of time shifts, with weight 1/t². The 1 is the zero shift. A check, not an input; the displayed bridge in "The path taken" does not use it. |
| 3 | Closed worldlines counted with the time-sense rule of IV.10. This is where the −½ and the absolute value come from: the walker's vacuum is the configuration the long-clock count of the main note selects — every growing slot occupied once, in every frame (main note VI.4 (C), VI.5) — and −½ Tr\|H\| is its energy. It is the energy of the walker read by all its detectors; a self-conjugate one, in the original detectors alone, would give half. What is taken is one thing: that the vacuum is what that limit selects. Feynman's rule and the source-free kernel then follow from exclusion. |
| 4 | Path counting. A massless particle between walls returns to its starting state only at t = 2nL/c. The factor 2L/c counts starting points and the two directions; e^{iφ} is the factor of the two wall flips in one round trip. |
| 5 | The sum over round trips. The n = 0 term is the zero-length path, proportional to L, so it is bulk. |

**The wall factor from the coin.** A wall is one of the two matrices of III.8 applied once. The current at a wall, c(|ψ₊|² − |ψ₋|²), vanishes for any relative phase ψ₋ = e^{iθ}ψ₊, so walls form a one-parameter family and the phase is a choice of boundary condition, not a consequence of the coin. The turn matrix (squares to −1, eigenvalues ±i) is the bag wall: each hit carries i, a round trip has two, so e^{iφ} = i·i = −1. The half-turn (eigenvalues ±1) is the other wall, φ = 0.

Stated without eigenvalues: a round trip applies the wall matrix twice, so its factor is the square of that matrix, −1 for the turn matrix and +1 for the telegraph's flip matrix. This assumes the same wall rule at both ends; in general the factor is the product of the two walls' phases.

```latex
E(L)=\frac{\hbar c}{2\pi L}\sum_{n=1}^{\infty}\frac{(-1)^{n}}{n^{2}}=\frac{\hbar c}{2\pi L}\cdot\Big(-\frac{\pi^{2}}{12}\Big)=-\frac{\pi\hbar c}{24L}
```

Attraction is the turn-matrix wall; the half-turn wall repels. The wall's phase is an input of the same kind as III.8's choice of the turn matrix for the mass flip. Values for L = 1, in units of ħc, with the same quantity computed independently from the modes using a smooth cutoff:

| Wall factor per round trip | Sum over closed paths | Mode sum with smooth cutoff, direct | Result |
| --- | --- | --- | --- |
| +1 (φ = 0) | +0.261799 | +0.261799 | repulsive |
| −1 (φ = π) | −0.130900 | −0.130899 | attractive, −π/24 |
| ±i (φ = π/2) | −0.032725 | −0.032725 | weakly attractive |

**What went in.**

- [ ] The energy of closed worldlines, taken as −½ Tr|H|. One particle alone pushes the walls apart; the attraction comes from loops, which the forward-time observer reads as pairs. (The |H| and the ½ now come from the main note — exclusion in every state, VI.1 (D), and one thing taken, that the vacuum is the configuration the long-clock count selects (VI.4 (C), VI.5), which counting the added detectors with the other sign (V.2 (H)) takes to itself — and the whole can be bypassed: "What the walk now supplies" in the closed-worldlines section.)
- [ ] The wall is the turn matrix of III.8 (phase i per hit), not the half-turn (phase ±1): a boundary condition chosen, not derived. The table shows the answer depends on it.
- [ ] Continuous time. On the lattice with a single step the weight vanishes for even path lengths (the top-of-band problem), so this is the continuum statement.
- [ ] Massless only. With mass, paths can also close by flipping between the walls, and line 4 is no longer a row of spikes.

## Loops that miss closing: finite separation and the master equation

Letting a loop end a distance δ from where it started gives a series in δ: the wall-induced energy at finite separation is the master equation at a = −2, and the Casimir energy is its zeroth coefficient. Every sum in this section converges.

**The paths.** A massless particle between walls a distance L apart moves at c and turns only at the walls. A zigzag from one point to another is labelled by its number of round trips n and has length 2nL + δ, where δ is the separation of the two points. The n = 0 path is the direct one and exists with or without walls.

**The weight of a path.** A loop of Part II's walk that ends a distance ℓ from its start, counted with the 1/τ weight, gives −½∫dτ/τ · e^(−ℓ²/4τ)/(4πτ) = −1/(2πℓ²). The same number is the mixed second derivative of the τ-summed kernel of III.16 at the two ends, written in terms of the real-time interval, with space and time parts contributing equally:

```latex
G=-\frac{1}{4\pi}\ln\!\big(\Delta x^{2}-c^{2}\Delta t^{2}\big),\qquad
\partial_x\partial_{x'}G\Big|_{\Delta t=0}=\frac{1}{c^{2}}\,\partial_t\partial_{t'}G\Big|_{\Delta t=0}=-\frac{1}{2\pi\,\ell^{2}},\qquad
w(\ell)=-\frac{\hbar c}{2\pi\,\ell^{2}}
```

The weight is linear in G, so images add with no cross terms.

**The sum.** With x = δ/2L:

```latex
\begin{aligned}
\mathcal H_{\text{walls}}(\delta)-\mathcal H_{\text{free}}(\delta)
&=-\frac{\hbar c}{2\pi}\sum_{n\neq0}\frac{1}{(2nL+\delta)^{2}}\\[6pt]
&=-\frac{\hbar c}{8\pi L^{2}}\Big[\sum_{n=1}^{\infty}\frac{1}{(n+x)^{2}}+\sum_{n=1}^{\infty}\frac{1}{(n-x)^{2}}\Big]\\[6pt]
&=-\frac{\hbar c}{8\pi L^{2}}\Big[2\zeta(2)-\sum_{k=1}^{x}\frac{1}{k^{2}}-\sum_{k=1}^{-x}\frac{1}{k^{2}}\Big]\\[6pt]
&=-\frac{\hbar c}{8\pi L^{2}}\Big[2\zeta(2)+\sum_{m=1}^{\infty}\binom{-2}{m}\zeta(m+2)\big(x^{m}+(-x)^{m}\big)\Big]\\[6pt]
&=-\frac{\hbar c}{4\pi L^{2}}\Big[\zeta(2)+3\,\zeta(4)\,x^{2}+5\,\zeta(6)\,x^{4}+7\,\zeta(8)\,x^{6}+\cdots\Big]
\end{aligned}
```

At δ = 0 only the first term is left:

```latex
\Delta\mathcal H=-\frac{\hbar c}{4\pi L^{2}}\,\zeta(2)=-\frac{\pi\hbar c}{24L^{2}},\qquad
E=L\,\Delta\mathcal H=-\frac{\pi\hbar c}{24L},\qquad
F=-\frac{dE}{dL}=-\frac{\pi\hbar c}{24L^{2}}
```

| Symbol in the master equation | Meaning here |
| --- | --- |
| k | number of round trips of the zigzag |
| a = −2 | minus the power of the path length in the weight |
| x | point separation in units of one round trip |
| sum of k^(−2) up to x | closed zigzags minus those ending a fraction x of a round trip away |
| coefficients ζ(4), ζ(6), ... | response of the energy to pulling the two points apart |

**Closed form and radius.** The sum over all n of 1/(n + x)² equals π²/sin²(πx), so the series converges for |x| < 1 and diverges at x = 1. That is δ = 2L: the two points are exactly one round trip apart.

**Checks.** At x = 0.2 the direct sum over zigzags gives 3.56685 for the bracket on line 2, and four terms of the series give 3.56681. As x goes to 0 the bracket tends to 2ζ(2) = 3.28987 (3.5669, 3.3062, 3.289875 at x = 0.2, 0.05, 0.001). The master equation at a = −2, x = 0.3 holds to 25 digits.

**Inputs and open points.**

- [ ] The loop weight, which fixes the π/24, is the input stated in the closed-worldline section.
- [ ] Paths with an odd number of wall hits are left out because their space and time parts cancel: on the reflected images (∂\_t∂\_t′ + ∂\_x∂\_x′)G is exactly 0 (the kernel route, `casimir_from_the_walk.py` 3). Reflecting and absorbing walls differ only in the sign of those paths, so the massless energy is the same for both.
- [ ] The image sum of G itself diverges; only the sum of its derivatives is used.
- [ ] These are loops of Part II's walk, with straight segments between wall hits in the massless case. The coin's version is its own section.
- [ ] 3+1 dimensions: weight 1/length⁴ would give ζ(4)/16π² = π²/1440 per unit area for this scalar loop weight (the electromagnetic value π²/720 is two polarizations of it); arithmetic checked, derivation not. The worldline route to the Casimir energy, with its 1/τ loop weight, is Gies, Langfeld and Moyaerts, JHEP 06 (2003) 018. The shifted-sum counterpart is Σ\_{n≥1} (n + x)^(−4) = ζ(4) − Σ\_{k=1}^{x} k^(−4); which separated-point observable a nonzero x is there has not been specified.
- [ ] The relation between this convergent side (ζ(2), ζ(4)) and the mode side (ζ(−1), ζ(−3)) is the reflection formula for zeta — which, for these sums, is Poisson summation, and is derived in "Modes and loops are Poisson partners" below: the m-th image of the mode sum is the loop that winds m times.
- [ ] Whether nudging the plates is equivalent to nudging the two points apart.

## Modes and loops are Poisson partners

The mode sum and the loop sum are not two regularizations that happen to agree; they are the two sides of one identity, Poisson summation, Σₙ f(n) = Σₘ f̂(2πm), and the m-th term on the right is the loop that winds m times between the walls. Take the mode sum with the exponential cutoff, f(n) = |n| e^{−ε|n|} (in units πħc/L = 1):

```latex
\hat f(k)=\int_{-\infty}^{\infty}|n|\,e^{-\varepsilon|n|}e^{-ikn}\,dn=\frac{2(\varepsilon^2-k^2)}{(\varepsilon^2+k^2)^2},\qquad \sum_{n\ge1}n\,e^{-\varepsilon n}=\frac12\Big[\underbrace{\hat f(0)}_{2/\varepsilon^2}+2\sum_{m\ge1}\hat f(2\pi m)\Big]=\frac1{\varepsilon^2}+\sum_{m\ge1}\frac{2(\varepsilon^2-4\pi^2m^2)}{(\varepsilon^2+4\pi^2m^2)^2}\ \xrightarrow{\ \varepsilon\to0\ }\ \frac1{\varepsilon^2}-\sum_{m\ge1}\frac{2}{4\pi^2m^2}=\frac1{\varepsilon^2}-\frac{2\zeta(2)}{4\pi^2}=\frac1{\varepsilon^2}-\frac1{12}.
```

Line by line: the m = 0 image is the bulk, 1/ε², the term proportional to L that is there without walls; every other image is finite as ε → 0 and is −2/(2πm)², the weight of the loop that winds m times (the closed-worldlines section's 1/(2mL)², in these units); their sum is −1/12. So ζ(−1) on the mode side *is* Σ 1/m² on the loop side, image by image, and the finite part that the fractional sum isolates ("Zeta from the fractional sum") is the loop sum. The Euler–Maclaurin corrections at a shifted endpoint carry the same number: summing to n₀ and correcting the tail with f(n₀)/2 − f′(n₀)/12 gives 1/ε² − 1/12 at n₀ = 5, 20 or 100 — the Bernoulli terms at the far end are the first Poisson images, which is what "go and come back" was doing all along — as a statement about sums. Two uses of the phrase are kept apart: shifting the endpoint of a sum to a far integer and correcting there is an identity of summation; a carrier crossing the gap and coming back is a propagation term. Poisson summation ties the one expansion to the other and does not make the endpoint shift a trajectory. What is physical about round trips is in the plates section.

The coin's turn-matrix wall has half-integer modes, f(n + ½): shifting n by ½ multiplies the m-th image by e^{iπm} = (−1)^m, and Σ\_{n≥0}(n+½)e^{−ε(n+½)} = 1/ε² + 1/24 with +1/24 = −2Σ\_{m≥1}(−1)^m/(4π²m²): the alternating images are the alternating round trips of the coin's table, wall factor −1 per hit. The half-turn wall (integer modes, (−1)^m absent) gives −1/12 and the sign of the force flips — the two rows of that table are one sign in the Poisson phase.

What this closes: the item in "Loops that miss closing" on how the convergent side (ζ(2), ζ(4)) is tied to the mode side (ζ(−1), ζ(−3)), the reflection formula. For these sums the reflection formula is Poisson summation, and Poisson summation is the statement that a sum over standing modes between two walls is a sum over round trips. (`scripts/modes_and_loops_are_poisson_partners.py`: 1 the finite part; 2 the images, term by term; 3 the half-integer case; 4 the Euler–Maclaurin endpoint at three n₀.)

## The plates, made of walkers: the Fabry–Pérot sum, closed

*The picture first.* Two plates of charged walkers face each other across a gap. Nothing is shone in. Each plate's charge jiggles on its own — its walkers step between the columns they fill and the columns they do not — and the carrier takes the jiggle across the gap, where it turns the other plate's walkers, which answer, and the answer comes back. Every round trip depends on the width of the gap; the plates' own jiggling does not. The force is the part that depends on the gap. Three things are derived below from the main note: what a plate does to the carrier (its reflection, from the Lorentz force on its walkers), that only round trips cross the gap (a determinant identity), and what the plates' own fluctuations are (their upward lines, from the filled sea). What is put in is the material and the state.

*Symbols.* κ is the carrier's decay rate across the gap at transverse wavenumber k and imaginary frequency ξ; r₁, r₂ the plates' reflections; L or d the gap; ω\_p the plasma frequency; χ a plate's charge response; v the carrier's kernel; ε a column's quasi-energy. Units ħ = c = 1 unless written. Nothing else is new.

**(A) Open paths and closed paths.** Classical optics adds up reflections too. A wave sent in from outside enters, bounces n times between the plates and leaves; summed,

```latex
t_{\rm total}=t_1t_2\sum_{n\ge0}\big(r_1r_2e^{2ik_zL}\big)^n=\frac{t_1t_2}{1-r_1r_2e^{2ik_zL}}\qquad(\text{Airy}),
```

an amplitude from a source to a detector, with no 1/n. The Casimir sum is the same reflections on paths that come back to where they started. A round trip wound n times is one loop with n starting points, so it is counted once in n:

```latex
\ln\big(1-r_1r_2e^{-2\kappa L}\big)=-\sum_{n\ge1}\frac1n\big(r_1r_2e^{-2\kappa L}\big)^n,
```

the "1/τ is counting" of the closed-worldlines section, with the same round-trip factor as Airy's. The force is the optics sum itself: differentiate in L and the 1/n cancels,

```latex
-\partial_L\ln\big(1-Re^{-2\kappa L}\big)=-2\kappa\sum_{n\ge1}\big(Re^{-2\kappa L}\big)^n=-\frac{2\kappa\,Re^{-2\kappa L}}{1-Re^{-2\kappa L}},\qquad R=r_1r_2,
```

each round trip weighted by the momentum 2κ it carries across the gap. And the modes are the zeros of the Airy denominator: 1 − r₁r₂e^{2ik\_zL} = 0 are the cavity's resonances, whose ½ħω summed and continued to imaginary frequency is the logarithm (van Kampen, Nijboer and Schram 1968). Mode sum: the zeros of the denominator; loop sum: the expansion of its logarithm — the Poisson partners of the previous section, in optics' words.

**(B) What the Gaussian supplies, in three lines.** The main note's VI.2 (F) gives the link field the weight e^{−½ a·Ka}, with K the free F² plus the plates' response. Three facts about any Gaussian, each one of the pictures above:

```latex
\begin{aligned}
\langle a\,a\rangle&=K^{-1}\equiv G\\[4pt]
\ln\det K&=\ln\det K_0+\operatorname{Tr}\ln(1+G_0V)=\ln\det K_0-\sum_n\frac{(-1)^n}n\operatorname{Tr}(G_0V)^n\\[4pt]
\partial_L\ln\det K&=\operatorname{Tr}\big(G\,\partial_LK\big)
\end{aligned}
```

Line 1: the covariance is the kernel — "the vacuum is the source-free kernel" and "the field has a Gaussian weight" are one sentence. Line 2: the determinant is the loops, free–plate–free–plate, counted with 1/n. Line 3: ∂\_L K lives only at the plate, so the force is the kernel read at the plate — the flux.

Shown on the smallest case, two surface penalties λᵢφ(zᵢ)²/2 a distance L apart. The free kernel restricted to the two surfaces at fixed κ, with q\_L = e^{−κL}, and the determinant it updates:

```latex
g_{\rm surf}=\frac1{2\kappa}\begin{pmatrix}1&q_L\\q_L&1\end{pmatrix},\qquad\det\Big[\begin{pmatrix}1&0\\0&1\end{pmatrix}+\frac1{2\kappa}\begin{pmatrix}1&q_L\\q_L&1\end{pmatrix}\begin{pmatrix}\lambda_1&0\\0&\lambda_2\end{pmatrix}\Big]=\Big(1+\frac{\lambda_1}{2\kappa}\Big)\Big(1+\frac{\lambda_2}{2\kappa}\Big)-\frac{\lambda_1\lambda_2q_L^2}{4\kappa^2}.
```

Divide by the two isolated-surface factors:

```latex
\frac{\det_{12}}{\det_1\det_2}=1-\frac{\lambda_1}{2\kappa+\lambda_1}\,\frac{\lambda_2}{2\kappa+\lambda_2}\,e^{-2\kappa L}=1-r_1r_2e^{-2\kappa L},\qquad r_i=\frac{\lambda_i}{2\kappa+\lambda_i}.
```

The gap enters only through round trips, and the reflection is what one surface does to the carrier. (`interaction_and_plates_checks.py` 5.) For plates of walkers rᵢ is derived next, and (G) shows the identity holds for the plates themselves, to 10⁻¹³.

**(C) The plate, as a material.** Take walkers of charge q, mass m, number density n₀, confined to a slab, with heavy opposite charges fixing neutrality and the slab's place. Two inputs, both stated: the material (these numbers, and that the walkers are slow, positive-energy and unpolarized, the leading local collisionless response), and the state (the Gaussian vacuum for the carrier, zero temperature). Nothing else is chosen; in particular no wall factor. The main note's VI.2 (D) gives, for a walker in the carrier's field, dπ/dt = q(E + v × B). Around a material at rest the induced velocity times B is second order, so to first order

```latex
m\,\partial_t\mathbf v=q\,\mathbf E,\qquad\partial_t\delta n+n_0\nabla\cdot\mathbf v=0\ \ (\text{the walkers' own conservation}),\qquad\mathbf j=qn_0\mathbf v,\qquad\partial_t\mathbf j=\frac{n_0q^2}m\,\mathbf E.
```

Both signs of charge add, through q². With the polarization ∂\_t 𝒫 = j and fields ∝ e^{−iωt}, −ω²𝒫 = (n₀q²/m)E, so

```latex
\epsilon(\omega)=1-\frac{\omega_p^2}{\omega^2},\qquad\omega_p^2=\frac{n_0q^2}m,\qquad\epsilon(i\xi)=1+\frac{\omega_p^2}{\xi^2}.
```

(The main note's rationalized units; in SI, ω\_p² = n₀q²/mε₀.) The retarded prescription is understood; the substitution ω → iξ evaluates the causal response where the force integral converges and changes no flip rate. This is the collisionless plasma, derived from the walker's Lorentz force under the stated approximations.

**(D) The wall rule, by matching.** Vacuum at z < 0, material at z > 0, wavevector k∥ along the surface; Maxwell from VI.2 (F) gives the normal wavenumbers k\_z² = ω² − k∥² and k\_mz² = ε(ω)ω² − k∥², which on the imaginary axis are iκ, iκ\_m with

```latex
\kappa=\sqrt{k_\parallel^2+\xi^2},\qquad\kappa_m=\sqrt{k_\parallel^2+\epsilon(i\xi)\,\xi^2}=\sqrt{\kappa^2+\omega_p^2}.
```

*Transverse electric* (electric field perpendicular to the plane of incidence; amplitudes in, reflected, transmitted = 1, r\_TE, t\_TE). Integrating Maxwell's equations across the interface makes the tangential electric and magnetic fields continuous; no singular surface sheet is added:

```latex
\begin{pmatrix}1&-1\\-k_z&-k_{mz}\end{pmatrix}\begin{pmatrix}r_{\rm TE}\\t_{\rm TE}\end{pmatrix}=\begin{pmatrix}-1\\-k_z\end{pmatrix}\qquad\Longrightarrow\qquad r_{\rm TE}(i\xi,k_\parallel)=\frac{\kappa-\kappa_m}{\kappa+\kappa_m}.
```

*Transverse magnetic* (the tangential magnetic amplitude; tangential-electric continuity brings in ε):

```latex
\begin{pmatrix}1&-1\\-k_z&-k_{mz}/\epsilon\end{pmatrix}\begin{pmatrix}r_{\rm TM}\\t_{\rm TM}\end{pmatrix}=\begin{pmatrix}-1\\-k_z\end{pmatrix}\qquad\Longrightarrow\qquad r_{\rm TM}(i\xi,k_\parallel)=\frac{\epsilon\kappa-\kappa_m}{\epsilon\kappa+\kappa_m}.
```

These are the wall rules the chosen walker material produces: a charge accelerated by the incoming field, its current feeding Maxwell, and the combined solution having a reflected part. No turn per wall hit is imposed. For a slab of thickness d with vacuum on both sides, the inside round trips sum to

```latex
r_{{\rm slab},p}=\frac{r_p\,(1-e^{-2\kappa_md})}{1-r_p^2e^{-2\kappa_md}},\qquad p={\rm TE},{\rm TM}:
```

the first reflection, plus transmission in, reflection at the back, and any number of extra inside round trips — Airy again, inside the plate. A thin or dilute slab transmits; the half-space rule returns when κ\_m d ≫ 1 over the modes that matter. For ω\_p → ∞ at fixed ξ, r\_TE → −1 and r\_TM → +1: two equal plates have round-trip product +1 in each polarization, the ideal rule, reached from the strong, optically thick response. (TE/TM by matching, and the slab: `plates_from_charged_walkers_checks.py` 1–2.)

**(E) The round trips, and the force.** Put the derived r's where the surface penalties were in (B). Eliminating the amplitudes in the gap gives, per polarization,

```latex
B_p=\begin{pmatrix}1&-r_{1p}e^{-\kappa L}\\-r_{2p}e^{-\kappa L}&1\end{pmatrix},\qquad\det B_p=1-r_{1p}r_{2p}e^{-2\kappa L},
```

the normalized round-trip factor once the single-plate factors are divided out. Integrated with the Gaussian's measure (zero temperature: (ħ/2π)∫₀^∞ dξ), the interaction energy per area and its L-derivative at fixed material are

```latex
\begin{aligned}
\frac E{\mathcal A}&=\frac\hbar{2\pi}\int_0^\infty d\xi\int\frac{d^2k_\parallel}{(2\pi)^2}\sum_{p}\ln\big(1-r_{1p}r_{2p}e^{-2\kappa L}\big)\\[6pt]
\frac F{\mathcal A}&=-\frac\hbar\pi\int_0^\infty d\xi\int\frac{d^2k_\parallel}{(2\pi)^2}\,\kappa\sum_p\frac{r_{1p}r_{2p}e^{-2\kappa L}}{1-r_{1p}r_{2p}e^{-2\kappa L}}
\end{aligned}
```

Lifshitz's formula, with the reflections fixed by the plates' walkers. Two identical passive plates have r₁ₚr₂ₚ = r\_p² ≥ 0: attraction. The field and the material modes are counted together in this one expression; a separate surface-mode energy added afterwards would count the same modes twice. (The use of scattering data in the determinant is Rahi, Emig, Graham, Jaffe and Kardar, *Phys. Rev. D* 80, 085021 (2009); the thin plasma sheet and its surface plasmons, Bordag, *J. Phys. A* 39, 6173 (2006).)

**(F) The ideal force emerges.** For identical thick plates, with α = ω\_p L, x = κL, μ = ξ/κ ∈ \[0, 1\]:

```latex
\frac{F\,L^4}{\mathcal A\,\hbar c}=-\frac1{2\pi^2}\int_0^\infty dx\,x^3\int_0^1d\mu\sum_p\frac{r_p(x,\mu)^2e^{-2x}}{1-r_p(x,\mu)^2e^{-2x}},
```

against the ideal pressure −π²ħc/240L⁴ (which the same integral gives with r\_p² = 1: `plates_from_charged_walkers_checks.py` 3):

| ω\_p L/c | F/F\_ideal |
| --- | --- |
| 0.1 | 0.0188044 |
| 1 | 0.1538280 |
| 5 | 0.4655710 |
| 10 | 0.6396879 |
| 50 | 0.9022738 |
| 200 | 0.9739226 |

Gauss–Legendre on both axes, 120 and 240 nodes, agree to the digits shown; 0 < x < 40, the tail exponentially small (ibid. 4; the series of the summation part of this note at x = 0.3, ibid. 5). This is the response integral and its approach to ideal reflection, not a simulation of an assembled slab. The next item is that simulation.

**(G) From the plates' own walkers.** (C)–(F) take the plate's response as the continuum plasma and the field's fluctuations from the Gaussian. The walk can do without the second: the fluctuations are the plates' own, and the sea makes them what a vacuum needs. Two runs on the 2+1 walk (a 1+1 neutral plate has no field outside it, by Gauss's law in one dimension, so two such plates do not feel each other; in 2+1 the carrier's static kernel at transverse momentum k is e^{−κ|z−z′|}/(2 sinh κ) with cosh κ = 2 − cos k, the walk's own occupation per transverse mode, `two_walker_plates_casimir_2plus1.py` 2).

*The plate is not a lamp.* Take one slab where the walkers' turn is weak (a mass well), fill every column below E\_F — the sea and a band of conduction walkers bound in the slab — and read the plate's charge at one site. The plate's charge jiggles by its *lines*, one from each filled column o to each empty column e, of rate Ω = ε\_e − ε\_o and weight |N\_eo|², with N\_eo the site overlap of the two columns (the sum over components of the conjugate of e times o). With the sea filled every line goes up, Ω > 0: the symmetrized correlation and the response are one object,

```latex
C(z,z')=-\frac1\pi\int_0^\infty d\xi\;\tfrac12\big[\chi(i\xi;z,z')+\overline{\chi(i\xi;z',z)}\big],\qquad\chi(i\xi)=\sum_{\rm lines}\frac{\overline{N(z)}\,N(z')}{i\xi-\Omega}\ \text{(and the reversed lines)},
```

to 7×10⁻⁹ (`two_walker_plates_casimir_2plus1.py` 3). That is the zero-temperature limit, coth(Ω/2T) → 1, of the plate's own sum-and-difference relation — "the ½ħω of the plate" with no field zero-point used. The mechanism is exchange: a lone conduction walker in the slab is a lamp, 93% of its line weight downward at Ω ≈ −2mc² (pair annihilation into the sea's empty columns); with the sea filled those columns are taken and its downward weight is exactly cancelled, line by line, by the exchange weights of the walkers below it (`plate_of_walkers_sum_and_difference.py` 3, 6). An excited plate (one walker lifted) has one downward line and is a lamp again; a hot plate obeys (W\_up + W\_down)/(W\_up − W\_down) = coth(Ω/2T) on 29,765 lines to 10⁻¹⁴ (ibid. 4, 5). So "not a lamp" is "no walker can step down", which is "the lowest columns are filled": the state input of (C), said in the walk's words.

*Two plates, near field.* Two such slabs, centre to centre d, coupled by the carrier's static kernel v = q²e^{−κ|z−z′|}/(2 sinh κ), which is the near-field end of the end-to-end rule (good while ad ≪ 1, a = q²D/2 with D the plate's Drude weight). The chain, each line checked:

```latex
\begin{aligned}
R&=\chi\,(1-v\chi)^{-1}\\[4pt]
r_R&=\frac{q^2}{2\sinh\kappa}\sum_{z,z'}e^{\kappa(z-z_c)}R(z,z')\,e^{\kappa(z'-z_c)}\\[4pt]
\frac{\det(1-v\chi)_{\rm both}}{\det(1-v\chi)_1\det(1-v\chi)_2}&=1-r_Rr_Le^{-2\kappa d}\\[4pt]
E(d)&=\int\frac{dk}{2\pi}\frac1{2\pi}\int_0^\infty d\xi\,\ln\big(1-r_Rr_Le^{-2\kappa d}\big)\\[4pt]
F&=-\frac{dE}{dd}=-\big\langle\partial_dV_{12}\big\rangle=\int\frac{dk}{2\pi}\sum_{z_1,z_2}\kappa\,v(z_1,z_2)\,C_{12}(z_1,z_2)
\end{aligned}
```

Line 1: one plate, screened. Line 2: what it sends back when e^{κz} arrives; r\_L is the same with e^{−κ(z−z\_c)}. Line 3 holds to 10⁻¹³ at every k, ξ, d: only round trips cross the gap — (B)'s identity, for the plates themselves. Line 5 holds to 10⁻¹⁶, with C₁₂ the gap-induced correlation of the two plates' charges, the plate-1–plate-2 block of −(1/π)∫dξ Herm\[χ(1 − vχ)⁻¹\]: negative (opposite charges face each other), falling as e^{−κd} (one crossing), while the force per k falls as e^{−2κd} (a round trip); the first round trip carries 86% of the force (`two_walker_plates_casimir_2plus1.py` 5, 6). Against Lifshitz's near-field formula with the material read off the walkers (their density profile, v\_F and D = 2v\_F/π from the walkers' own band; nothing fitted): agreement to 1.3% at d = 24 sites, 0.17% at 48, 0.04% from 96 on; the zero-thickness sheet is off by the slab's thickness against d (9% → 2%), and the plasma model of (C) is reached only as v\_F²/ad → 0, since it leaves out the walkers' own motion (ibid. 7; 7 checks, 2.5 min). The sea's own lines were left out: they are the vacuum's polarization, a dielectric constant 1 + 4×10⁻³ at this coupling — a renormalization of q², under 1% of the force.

So in the near field the whole force is the plates' own upward lines, carried across the gap and back. The carrier's own fluctuations — the Hadamard half of its real-time kernel, which a direct interaction between walkers does not have — were not used. That is the source-theory route to Casimir (Schwinger, DeRaad and Milton, *Ann. Phys.* 115, 1 (1978)), in the walk's words.

**(H) Inputs and open points.** Inputs: the material and the state of (C); in (G) the state is "every column below E\_F filled". Derived: the reflection (D), the round-trip structure (B), (E) and its exact lattice form (G), the force two ways (G), the plates' fluctuations as their response (G). Open: the far field — (G) uses the carrier static; Lifshitz's retarded κ = √(k² + ξ²/c²) needs the moving end-to-end rule between the plates' walkers, and whether the time-symmetric (Fokker) kernel gives it or needs the absorber condition to act as retarded is not decided here; 3+1 with two polarizations; collisions, binding and surface structure of a real metal.

(`scripts/plates_from_charged_walkers_checks.py` (D)–(F); `scripts/interaction_and_plates_checks.py` 5 for (B); `scripts/plate_of_walkers_sum_and_difference.py` and `scripts/two_walker_plates_casimir_2plus1.py` for (G).)

## Open points

- [ ] Way 1 is operator identities checked numerically, with no remainder bound.
- [ ] Way 2: which f have a representation with some φ. The three tested follow φ = (−1)^p (p−1)!/s^p for t^(p−1)·ln t.
- [ ] Way 2: a general f differs from the represented form by a polynomial of degree below p, summed separately by its closed form.
- [ ] The Taylor form's radius: it is not always set by a pole of S at x = −1 (for f(t) = t the sum x(x+1)/2 is entire); for which f the singularity at −1 is a pole, and what the radius is in general, is function-specific.
- [ ] The plates' far field: the carrier was used static (the near field, ad ≪ 1 in the plates section's units); retardation and the current–current term need the moving end-to-end rule, not yet run. Whether the plates' own fluctuations alone give the retarded Lifshitz force at all separations is decided there.
- [ ] The plates are 2+1 (a 1+1 neutral plate has no field outside it); the 3+1 electromagnetic run, with two polarizations, is not done.
