# Three random particles, three wave equations

**Part I.** A particle that takes random steps in space as time ticks → the Schrödinger equation (Nelson), every step written out, with the inputs listed in the ledger at the end of the Part. **Part II.** The same steps taken in spacetime — time included, with the sign relativity forces on it, which makes the time kick imaginary: a formal walk, real only in Euclidean signature (III.16), whose moments algebra is Part I’s → the Klein–Gordon equation. **Part III.** A particle that always advances at speed \\(c\\) while the sense in which *time* runs along its path flips at random (the time-sense particle) → the Dirac equation, by the same steps, with a real flip rate: nothing is continued in the rate — the sign that the continuation supplies elsewhere is here the choice of which coordinate is the clock (III.5a (J)), and the sum over paths is real (III.15a); the process’s own solutions are the kernel outside the cone, not the states a detector reads (III.15b). The walk read by \\(x\\) and read by \\(t\\) are one function on the two sides of the light cone (III.15); which of the two readings a detector makes is fixed by the causal cut (III.15b); that the record is the product of the two legs follows from the particle’s null tracks and frame independence (III.15c); its square is Part II. The telegraph particle (space direction random, time the clock) is kept as a remark: it needs its rate continued. Then what one such particle does, exactly: its kernel outside the light cone (III.14), the Dirac amplitude inside, which is the same function continued across the cone (III.15), and the source answer that Part II gives when its step label is summed (III.16). **Part IV.** Three dimensions: the same construction produces the vortex axis as the flatlander holds it, the two senses of circulation swapped by the mass, the spin density \\(\\frac\\hbar2\\rho\\hat{\\mathbf s}\\) as a canonical momentum, and the spin current \\(\\frac\\hbar2\\nabla\\times(\\rho\\hat{\\mathbf s})\\). Why the direction cannot be a classical label in three dimensions, with the mechanism (IV.2a); the propagators in three dimensions, on both sides of the light cone, checked against the three-coin automaton — a unitary lattice step, not a process with positive weights, which IV.2a rules out (IV.10). **Part V.** Where the electromagnetic field enters the walk: as a turn of the walker in its own plane that leaves its time sense alone, the only count-keeping knob a field can turn (V.1–V.2); the Pauli term and \\(g=2\\) from the step’s own noise (V.3); the effective flip rate (V.4); what a detector registers and why it is the square (V.5). And, in Part III, what an event is: the product of the walk’s two directions between an emitter and an absorber (III.15a). **Part VI.** Two walkers: a point of spacetime carries one process per direction, because the noise is the event’s (II.1); with one process per state as the input, identical walkers are odd under exchange — antisymmetry from exclusion in every state, by polarization; the free step keeps the sign and does not choose it (VI.1); a constraint on paths leaks and an elastic collision is a third statistics; Part II’s undirected walker meets no rule (bosons). The field from within: the angle of the clock hand is a connection — the coin does not read the hand, a local clock-zero is a potential with zero field, only loops read it, the force is the gradient of the turn angle per tick (the Lorentz force), the carrier is Part II’s walker at \\(m=0\\) with Coulomb’s law as its source answer, and the field’s own weight is \\(F^2\\), forced by no-preferred-frame and the relabelling (VI.2). Bell: two identical walkers in one packet are, by exclusion, the turn \\(T\\) in their axes, and the exchange term of VI.1 is the whole violation, \\(\\mathcal E=-\\cos(a-b)\\); two strangers scattered by the field give the same pair with \\(C=2|rt|/(|r|^2+|t|^2)\\) and \\(|S|=2\\sqrt{1+C^2}\\), quantum mechanics’ numbers; what the walk lacks of Bell’s three is outcome-independence, since an event has two endpoints (VI.3). Many walkers: a configuration on a cut is a determinant, adding and removing a walker are two operations on it, and the vacuum is the configuration that the count of exclusive configurations along the walk’s own clock selects; from it, Feynman’s rule, pair creation and the vacuum’s energy and polarization are computed (VI.4–VI.5); a missing walker is, in every amplitude, one walker of the opposite orientation and charge (VI.6); stopped at a finite length, the same count is a partition function, with the clock’s length as the inverse temperature, and the settled configuration seen from inside a region shares with the rest a third of a logarithm per edge (VI.8).

**Date:** 2026-11-09 (III.14–16, IV.2a, IV.10 added 10–11 October; III.15d–f, VI.1a, VI.2 (F′)–(G), VI.4–VI.6 added 9 November; v34: corrections from the review of the same day; v35: corrections from the scripts rebuilt on 10 November; v36: corrections from the re-verification of the scripts of 31 October, same day; v37: the same pass over the scripts written since, same day; v38: VI.7 added; v39: VI.7 in the document’s own names; v40: VI.7 (3) labelled; v41, v42: VI.7 (7a), same day; v43, 11 November: a remark in VI.1; v44, 12 November: the second column is read as detectors — III.8b, III.8c (B), III.15e, III.15f, IV.5a, V.2, V.4 rewritten, IV.5b new, input 11 restated; no equation of v43 changed except the sign of one term in V.4; v45, same day: the last section, “Next questions”, added; v46, same day: its item 6; v47, same day: its item 5, the thermodynamic route; v48, same day: VI.8, the count as a partition function and three entropies, from the draft of 11 November). Every numerical check named in the text is a script in `scripts/` next to this file. Every Part ends with a ledger of what it put in; the one list of all of them is the last section. Notation: \\(\\rho\\) density, \\(S\\) phase, \\(D=\\hbar/2m\\), \\(\\{\\rho(x),S(x')\\}=\\delta(x-x')\\), \\(\\psi=\\sqrt\\rho\\,e^{iS/\\hbar}\\). Part I in one space dimension, three dimensions noted where it differs; Part II in one space + one time dimension, signature \\((-,+)\\); Part III in one space + one time dimension with \\(x\\) as the clock. Every matrix is written out. Each step says what it uses. Two conventions: a field written without arguments is evaluated at the point \\((x,t)\\) that appears on the left-hand side of the same equation; and \\(\\langle\\,\\cdot\\,\\rangle\\) written with a condition, \\(\\langle\\,\\cdot\\,|\\,x_n=x\\rangle\\), is the average over the next step for a particle now at \\(x\\).

---

# Part I — The diffusing particle → Schrödinger

## I.1 The process

Time is cut into steps \\(\\Delta t\\). At each step the particle moves by a drift plus a random kick:

\\\[x\_{n+1}=x_n+b\_+(x_n,t_n)\\,\\Delta t+\\sqrt{2D\\,\\Delta t}\\;\\xi_n,\\\]

with \\(\\xi_n\\) a random number of mean \\(0\\) and variance \\(1\\), drawn fresh at each step, independent of \\(x_n\\) and of all earlier draws. \\(b\_+\\) is the *forward drift*: the mean velocity given where the particle is now. \\(D\\) is the diffusion constant. (In three dimensions: one such line per component, with independent \\(\\xi\\)’s.)

What is read off the rule directly, using only \\(\\langle\\xi\\rangle=0\\) and \\(\\langle\\xi^2\\rangle=1\\):

\\\[\\langle x\_{n+1}-x_n\\rangle=b\_+\\Delta t,\\qquad\\langle(x\_{n+1}-x_n)^2\\rangle=2D\\Delta t+(b\_+\\Delta t)^2 .\\\]

What is *not* in the rule as written: the shape of the distribution of \\(\\xi\\), and the probability \\(W(x,t+\\Delta t\\,|\\,x',t)\\) of being at \\(x\\) at time \\(t+\\Delta t\\) given that it was at \\(x'\\) at time \\(t\\) — read “\\(W\\) of \\(x\\) at \\(t+\\Delta t\\), given \\(x'\\) at \\(t\\)”. The time interval is part of the object: a step of \\(\\Delta t\\) and a step of \\(2\\Delta t\\) are different \\(W\\)’s. \\(W\\) is used for these one-step probabilities throughout, to keep them apart from densities. I.2 needs it. Both come out of building the process from a coin, next.

The kick is of size \\(\\sqrt{\\Delta t}\\), not \\(\\Delta t\\). Everything that follows hangs on that, so here is where it comes from.

### I.1a Where the update rule comes from: a coin

Here is the process as you would actually build it. A clock ticks every \\(\\tau\\). At each tick a fair coin is thrown and the particle steps \\(+\\ell\\) or \\(-\\ell\\). Two numbers, \\(\\ell\\) and \\(\\tau\\), both finite. Now compute what you see.

*One tick.* \\(x\_{\\mathrm{new}}=x\_{\\mathrm{old}}\\pm\\ell\\) by the coin, independent of \\(x\_{\\mathrm{old}}\\): for one tick, \\(W(x,t+\\tau\\,|\\,x',t)=\\tfrac12\[\\delta(x-x'-\\ell)+\\delta(x-x'+\\ell)\]\\). Density now = sum over where you were, times the probability of the step:

\\\[\\rho(x,t+\\tau)=\\int W(x,t+\\tau\\,|\\,x',t)\\,\\rho(x',t)\\,dx'=\\tfrac12\\!\\int\\delta(x-x'-\\ell)\\rho(x')\\,dx'+\\tfrac12\\!\\int\\delta(x-x'+\\ell)\\rho(x')\\,dx'=\\tfrac12\\rho(x-\\ell,t)+\\tfrac12\\rho(x+\\ell,t).\\\]

Expand, \\(\\rho(x\\pm\\ell)=\\rho\\pm\\ell\\partial_x\\rho+\\tfrac12\\ell^2\\partial_x^2\\rho+\\ldots\\), \\(\\rho(t+\\tau)=\\rho+\\tau\\partial_t\\rho+\\ldots\\); the \\(\\pm\\ell\\partial_x\\rho\\) cancel:

\\\[\\tau\\,\\partial_t\\rho=\\tfrac12\\ell^2\\partial_x^2\\rho\\qquad\\Longrightarrow\\qquad\\partial_t\\rho=\\frac{\\ell^2}{2\\tau}\\partial_x^2\\rho,\\qquad D=\\frac{\\ell^2}{2\\tau}.\\\]

*The next tick.* The same operation again, since the next coin is fresh: \\(\\rho(x,t+2\\tau)=\\int W(x,t+2\\tau|x'',t+\\tau)\\,\\rho(x'',t+\\tau)\\,dx''=\\int\\Big\[\\int W(x,t+2\\tau|x'',t+\\tau)\\,W(x'',t+\\tau|x',t)\\,dx''\\Big\]\\rho(x',t)\\,dx'\\), and the bracket is the two-tick probability, \\(W(x,t+2\\tau|x',t)\\):

\\\[W(x,t+2\\tau\\,|\\,x',t)=\\int W(x,t+2\\tau|x'',t+\\tau)\\,W(x'',t+\\tau|x',t)\\,dx''=\\tfrac14\\big\[\\delta(x-x'-2\\ell)+2\\delta(x-x')+\\delta(x-x'+2\\ell)\\big\]\\\]

— at \\(\\pm2\\ell\\) with probability \\(\\tfrac14\\) each, back at \\(0\\) with probability \\(\\tfrac12\\). After \\(N\\) ticks the same nested integral \\(N\\) times gives \\(W(x,t+N\\tau\\,|\\,x',t)\\): the binomial distribution on the points \\(x'+(N-2j)\\ell\\), Pascal’s triangle. **The elapsed time is the number of nested integrals**, \\(N=\\Delta t/\\tau\\). Nesting \\(N\\) integrals is tedious in \\(x\\) and a product in \\(k\\), which is the reason to transform.

*In Fourier space, to see the kicks.* \\(W(x,t+\\tau|x',t)\\) depends only on \\(x-x'\\), so the one-tick integral is a convolution and its transform is a product: with \\(\\hat\\rho(k)=\\int\\rho(x)e^{-ikx}dx\\) and \\(\\hat W(k;\\tau)=\\int W(x,t+\\tau|0,t)\\,e^{-ikx}dx=\\tfrac12(e^{-ik\\ell}+e^{ik\\ell})=\\cos k\\ell\\), one tick multiplies \\(\\hat\\rho\\) by \\(\\cos k\\ell\\), and \\(N\\) ticks multiply it by \\((\\cos k\\ell)^N\\) — that is, \\(\\hat W(k;N\\tau)=\\hat W(k;\\tau)^N\\): the elapsed time is the exponent. For \\(k\\ell\\ll1\\), \\(\\cos k\\ell=1-\\tfrac12k^2\\ell^2+\\ldots=e^{-k^2\\ell^2/2}\\) to this order, so

\\\[\\hat W(k;N\\tau)=(\\cos k\\ell)^N\\;\\to\\;e^{-Nk^2\\ell^2/2}\\qquad(k\\ell\\to0,\\ N\\to\\infty,\\ Nk^2\\ell^2\\ \\text{fixed}),\\\]

whose inverse transform is a Gaussian of variance \\(N\\ell^2\\). So after \\(N\\) ticks the particle has moved a Gaussian amount of width \\(\\ell\\sqrt N\\) — the width grows with the *square root* of the number of ticks because the transform is a power of something slightly less than one, and powers turn \\(1-\\epsilon\\) into \\(e^{-N\\epsilon}\\). Nothing was assumed beyond \\(k\\ell\\ll1\\) and many ticks — which is to say, this is the central limit theorem for this coin, as the next paragraph says.

*Now coarse-grain.* Look at the particle only every \\(\\Delta t=N\\tau\\), many ticks at a time. What you see between looks is the sum of \\(N\\) coin steps: a Gaussian kick of width

\\\[\\ell\\sqrt N=\\ell\\sqrt{\\Delta t/\\tau}=\\sqrt{2D\\,\\Delta t}.\\\]

That is the update rule of I.1. The \\(\\sqrt{\\Delta t}\\) is the \\(\\sqrt N\\) of the coin, with \\(N=\\Delta t/\\tau\\) ticks per look.

*Drift.* Load the coin: \\(+\\ell\\) with probability \\(p\\), \\(-\\ell\\) with \\(q\\), \\(p-q=b\\tau/\\ell\\). The transform of one tick is \\(p\\,e^{-ik\\ell}+q\\,e^{ik\\ell}=\\cos k\\ell-i(p-q)\\sin k\\ell=e^{-ikb\\tau-k^2\\ell^2/2}\\) to the same order, and \\(N\\) ticks give \\(e^{-ikb\\Delta t-k^2\\ell^2N/2}\\). Mean \\(\\propto\\Delta t\\), spread \\(\\propto\\sqrt{\\Delta t}\\), from the same coin: the bias adds up tick by tick, the spread only grows as \\(\\sqrt N\\).

*Transform back: the one-step probability.* After \\(N\\) ticks, \\(\\hat\\rho(k,t+\\Delta t)=e^{-ikb\\Delta t-Nk^2\\ell^2/2}\\,\\hat\\rho(k,t)\\). A product in \\(k\\) is a convolution in \\(x\\), so

\\\[\\rho(x,t+\\Delta t)=\\int W(x,t+\\Delta t\\,|\\,x',t)\\,\\rho(x',t)\\,dx',\\qquad W(x,t+\\Delta t\\,|\\,x',t)=\\frac1{2\\pi}\\int e^{-ikb\\Delta t-Nk^2\\ell^2/2}\\,e^{ik(x-x')}\\,dk,\\qquad\\Delta t=N\\tau,\\\]

\\(W(x,t+\\Delta t|x',t)\\) being the inverse transform of the factor: the probability that \\(N=\\Delta t/\\tau\\) coin steps carry the particle from \\(x'\\) at \\(t\\) to \\(x\\) at \\(t+\\Delta t\\). The transform of a Gaussian is a Gaussian — the standard pair is

\\\[\\frac1{2\\pi}\\int e^{-\\alpha k^2}e^{ikz}\\,dk=\\frac{1}{\\sqrt{4\\pi\\alpha}}\\,e^{-z^2/4\\alpha}\\\]

(a width \\(1/\\sqrt{2\\alpha}\\) in \\(k\\) becomes a width \\(\\sqrt{2\\alpha}\\) in \\(z\\); the factor \\(e^{-ikb\\Delta t}\\) shifts \\(z\\) by \\(b\\Delta t\\)). With \\(\\alpha=N\\ell^2/2\\) and \\(z=x-x'-b\\Delta t\\):

\\\[W(x,t+\\Delta t\\,|\\,x',t)=\\frac{1}{\\sqrt{2\\pi N\\ell^2}}\\exp\\Big\[-\\frac{(x-x'-b\\Delta t)^2}{2N\\ell^2}\\Big\]=\\frac{1}{\\sqrt{4\\pi D\\Delta t}}\\exp\\Big\[-\\frac{(x-x'-b\\Delta t)^2}{4D\\Delta t}\\Big\],\\qquad N\\ell^2=2D\\Delta t .\\\]

**This is the Gaussian.** It was not assumed: it is \\((\\cos k\\ell)^N\\) transformed back, Gaussian because a power of \\(1-\\epsilon\\) is \\(e^{-N\\epsilon}\\). (This computation *is* the central limit theorem, for this coin: the transform of a sum is the product of the transforms, and the product of many factors near \\(1\\) is \\(e^{-Nk^2\\ell^2/2}\\) whatever the fine shape of the factor.) The \\(\\xi\\) of I.1 is the sum of the \\(N\\) coin steps, minus its mean, in units of its width: \\(\\xi=(x\_{n+1}-x_n-b\\Delta t)/\\sqrt{2D\\Delta t}\\), with density \\(e^{-\\xi^2/2}/\\sqrt{2\\pi}\\). With \\(b\\to b\_+(x')\\) (the loading of the coin may depend on where you are), this \\(W(x,t+\\Delta t\\,|\\,x',t)\\) is what I.2 and I.3 use. The \\(\\Delta t\\) is on the right-hand side twice: as \\(N\\tau\\), the number of ticks nested, and inside the Gaussian’s width.

*What “\\(\\Delta t\\to0\\)” will mean in I.2:* \\(\\Delta t\\) small compared with the scales on which \\(\\rho\\) varies, but still many ticks. Two consequences used there: \\(\\Delta x^2=2D\\Delta t\\) is first order in \\(\\Delta t\\), so the second-order Taylor term survives; and \\(\\Delta x/\\Delta t=\\sqrt{2D/\\Delta t}\\) grows as you look more often, so the particle has a drift and a spread but no velocity (the \\(2Dd/\\Delta t\\) of I.5).

*The other coin.* Part III throws a different coin: at each tick the particle advances in \\(x\\) at speed \\(c\\), and with probability \\(a\\tau\\) the sense in which its time coordinate runs reverses. Then consecutive steps are not independent and the transform of \\(N\\) ticks is not a power of a number but a power of a \\(2\\times2\\) matrix; what comes out is in III.5.

## I.2 Forward Fokker–Planck

*Where the density comes from one step later.* To be at \\(x\\) now, the particle was at some \\(x'\\) one step ago and the step took it from \\(x'\\) to \\(x\\). Add up over every possible \\(x'\\), weighting by how likely it was to be there and how likely that step was:

\\\[\\rho(x,t+\\Delta t)=\\int W(x,t+\\Delta t\\,|\\,x',t)\\,\\rho(x',t)\\,dx'.\\\]

This is the one-tick line of I.1a, with \\(W\\) the Gaussian derived there for the interval \\(\\Delta t\\); the \\(\\Delta t\\) on the right sits inside \\(W\\).

*Why it can be expanded.* The step is small, of order \\(\\sqrt{\\Delta t}\\), so only \\(x'\\) close to \\(x\\) contribute. Write the step as \\(z=x-x'\\) and integrate over the step instead of the starting point:

\\\[\\rho(x,t+\\Delta t)=\\int W(x,t+\\Delta t\\,|\\,x-z,t)\\,\\rho(x-z,t)\\,dz.\\\]

Both factors are evaluated at the starting point \\(x-z\\): the density there, and the probability of the step \\(z\\) from there (through \\(b\_+(x-z)\\)). So the thing to expand in \\(z\\) is the product, as a function of the starting point. For a fixed step \\(z\\), let the starting point move from \\(x\\) to \\(x-z\\):

\\\[W(x,t+\\Delta t|x-z,t)\\,\\rho(x-z,t)=\\Big\[1-z\\,\\partial_x+\\tfrac12z^2\\partial_x^2-\\ldots\\Big\]\\Big\[W(x+z,t+\\Delta t\\,|\\,x,t)\\,\\rho(x)\\Big\],\\\]

the derivatives acting on everything to their right (the bracket is “density at the start, times probability of the step \\(z\\) from the start”, with the start now at \\(x\\)). Put this under the integral. The derivatives do not involve \\(z\\), so they come out of the \\(z\\)-integral, and what stays inside are averages over the step starting from \\(x\\) — the average of anything is \\(\\int(\\text{anything})\\,W(x+z,t+\\Delta t\\,|\\,x,t)\\,dz\\) — which I.1 already gave:

\\\[\\int W(x+z,t+\\Delta t\\,|\\,x,t)\\,dz=1,\\qquad\\int z\\,W(x+z,t+\\Delta t\\,|\\,x,t)\\,dz=\\langle x\_{n+1}-x_n\\rangle=b\_+(x)\\Delta t,\\qquad\\int z^2\\,W(x+z,t+\\Delta t\\,|\\,x,t)\\,dz=\\langle(x\_{n+1}-x_n)^2\\rangle=2D\\Delta t+O(\\Delta t^2).\\\]

Therefore, writing \\(\\rho\\) for \\(\\rho(x,t)\\) and \\(b\_+\\) for \\(b\_+(x,t)\\),

\\\[\\rho(x,t+\\Delta t)=\\rho-\\partial_x\\big(\\rho\\,b\_+\\big)\\Delta t+D\\Delta t\\,\\partial_x^2\\rho+O(\\Delta t^{3/2}),\\\]

the next term carrying \\(\\langle z^3\\rangle\\), which is \\(O(\\Delta t^{3/2})\\) because the step is \\(O(\\sqrt{\\Delta t})\\). Subtract \\(\\rho\\), divide by \\(\\Delta t\\), let \\(\\Delta t\\to0\\):

\\\[\\boxed{\\partial_t\\rho=-\\partial_x(\\rho\\,b\_+)+D\\,\\partial_x^2\\rho}\\qquad(\\mathrm{FP}\_+)\\\]

Drift carries the density along; the kicks spread it. For the coin of I.1a with no bias this is \\(\\partial_t\\rho=D\\partial_x^2\\rho\\) again, as it should be. In three dimensions nothing changes except \\(\\partial_x\\to\\nabla\\): \\(\\langle z_iz_j\\rangle=2D\\Delta t\\,\\delta\_{ij}\\) (independent components), so the spreading term is \\(D\\nabla^2\\rho\\) and the drift term is \\(-\\nabla\\cdot(\\rho\\,\\mathbf b\_+)\\). From here on the equations are written in three dimensions.

## I.3 The same process run backward: the backward drift

*The question.* Given that the particle is at \\(x\\) at time \\(t+\\Delta t\\), where was it at time \\(t\\)? (\\(W(x,t+\\Delta t|x',t)\\) is the Gaussian of I.1a.) The probability that it was at \\(x'\\) is proportional to (how likely it was to be at \\(x'\\)) \\(\\times\\) (how likely the step from \\(x'\\) to \\(x\\) was):

\\\[W(x',t\\,|\\,x,t+\\Delta t)=\\frac{\\rho(x',t)\\,W(x,t+\\Delta t|x',t)}{\\rho(x,t+\\Delta t)},\\\]

the denominator being the normalization (it does not depend on \\(x'\\)).

*Read it as a step backward from \\(x\\).* We want to know how this depends on the backward step \\(x'-x\\), which is of order \\(\\sqrt{\\Delta t}\\). Put the two \\(x'\\)-dependent factors in one exponent:

\\\[\\rho(x',t)\\,W(x,t+\\Delta t|x',t)\\propto\\exp\\Big\[-\\frac{\\big((x'-x)+b\_+\\Delta t\\big)^2}{4D\\Delta t}+\\ln\\rho(x',t)\\Big\]\\\]

(using \\(b\_+(x')\\Delta t=b\_+(x)\\Delta t+O(\\Delta t^{3/2})\\), negligible in the exponent). Expand \\(\\ln\\rho(x')=\\ln\\rho(x)+\\partial_x\\ln\\rho\\cdot(x'-x)+O(\\Delta t)\\) and complete the square in \\((x'-x)\\):

\\\[-\\frac{(x'-x)^2+2b\_+\\Delta t\\,(x'-x)}{4D\\Delta t}+\\partial_x\\ln\\rho\\,(x'-x)=-\\frac{(x'-x)^2+2\\big(b\_+-2D\\,\\partial_x\\ln\\rho\\big)\\Delta t\\,(x'-x)}{4D\\Delta t}=-\\frac{\\big((x'-x)+(b\_+-2D\\,\\partial_x\\ln\\rho)\\Delta t\\big)^2}{4D\\Delta t}+(\\text{no }x').\\\]

So, looked at backward, the step \\(x'-x\\) is again a Gaussian of variance \\(2D\\Delta t\\), now centred at \\(-(b\_+-2D\\partial_x\\ln\\rho)\\Delta t\\). A kick of the same size, a different centre: the time-reversed process is a process of the same kind as I.1. Write it in the form of I.1, going from \\(t+\\Delta t\\) to \\(t\\):

\\\[x_n=x\_{n+1}-b\_-(x\_{n+1})\\,\\Delta t+\\sqrt{2D\\Delta t}\\;\\xi,\\qquad\\boxed{b\_-=b\_+-2D\\,\\partial_x\\ln\\rho}\\\]

(\\(\\xi\\) Gaussian as in I.1, independent of \\(x\_{n+1}\\); the sign in front of it does not matter). \\(b\_-\\) is the *backward drift*: the mean velocity given where the particle is now, looking back.

*Backward Fokker–Planck.* Run I.2 on the reversed rule. It is the same calculation with the step going from \\(t+\\Delta t\\) to \\(t\\) and the centre \\(b\_+\\Delta t\\) replaced by \\(-b\_-\\Delta t\\):

\\\[\\rho(x,t)=\\rho(x,t+\\Delta t)+\\partial_x(\\rho\\,b\_-)\\Delta t+D\\Delta t\\,\\partial_x^2\\rho+O(\\Delta t^{3/2})\\qquad\\Longrightarrow\\qquad\\partial_t\\rho=-\\partial_x(\\rho\\,b\_-)-D\\,\\partial_x^2\\rho\\qquad(\\mathrm{FP}\_-)\\\]

The spreading term has the opposite sign: looking backward, the cloud contracts. In three dimensions, \\(\\partial_x\\to\\nabla\\) as in I.2.

It is FP\\(\_+\\) in disguise. Substitute \\(\\mathbf b\_-=\\mathbf b\_+-2D\\nabla\\ln\\rho\\):

\\\[-\\nabla\\cdot(\\rho\\mathbf b\_-)-D\\nabla^2\\rho=-\\nabla\\cdot(\\rho\\mathbf b\_+)+2D\\,\\nabla\\cdot(\\rho\\nabla\\ln\\rho)-D\\nabla^2\\rho=-\\nabla\\cdot(\\rho\\mathbf b\_+)+2D\\nabla^2\\rho-D\\nabla^2\\rho,\\\]

using \\(\\rho\\nabla\\ln\\rho=\\nabla\\rho\\). That is FP\\(\_+\\): one process, two descriptions.

*Two equations for one density.* FP\\(\_+\\) and FP\\(\_-\\) describe the same \\(\\rho\\); combine them. Add and halve — the \\(D\\nabla^2\\rho\\) terms cancel:

\\\[\\partial_t\\rho=-\\nabla\\cdot\\Big(\\rho\\;\\tfrac12(\\mathbf b\_++\\mathbf b\_-)\\Big).\\\]

This is a continuity equation: the density is carried along, with no spreading term at all, by the half-sum of the two drifts. So the half-sum is the velocity that actually transports probability — the *current velocity*,

\\\[\\boxed{\\partial_t\\rho+\\nabla\\cdot(\\rho\\,\\mathbf v)=0,\\qquad\\mathbf v=\\tfrac12(\\mathbf b\_++\\mathbf b\_-)}\\\]

Subtract and halve — the drifts’ difference is what I.3 computed:

\\\[\\tfrac12(\\mathbf b\_+-\\mathbf b\_-)=D\\nabla\\ln\\rho\\equiv\\mathbf u ,\\\]

the *osmotic velocity*. It is the part of the forward drift that is there only because the density is not flat: look at a Brownian cloud, mistake its spreading for a flow, and \\(\\mathbf u\\) is the flow you would write down. The two drifts are then \\(\\mathbf b\_\\pm=\\mathbf v\\pm\\mathbf u\\): the same transport, seen forward and backward, plus or minus the spreading.

## I.4 What I.1–I.3 give, and what they do not

So far: any density \\(\\rho\\), any drift \\(\\mathbf b\_+\\). From the process alone we got the backward drift, the current velocity \\(\\mathbf v\\) that transports \\(\\rho\\), and the osmotic velocity \\(\\mathbf u=D\\nabla\\ln\\rho\\). Nothing says how \\(\\mathbf v\\) changes in time. That is kinematics, and it is all a random walk can give.

*A check that the kinematics is right.* Take the one case where the answer is known: a Brownian particle in a fluid at temperature \\(T\\). Two facts about it come from thermal physics and are used here only for the check:

\\\[\\mathbf b\_+=-\\frac{\\nabla V}{m\\gamma}\\ \\ (\\text{drift = force / friction}),\\qquad\\rho\\propto e^{-V/k_BT}\\ \\ (\\text{Boltzmann}),\\ \\text{so}\\ \\ \\nabla\\ln\\rho=-\\frac{\\nabla V}{k_BT}.\\\]

In equilibrium the density is stationary and nothing flows: \\(\\mathbf v=0\\), which by I.3 means \\(\\mathbf b\_+=\\mathbf u=D\\nabla\\ln\\rho\\) — the drift current \\(\\rho\\mathbf b\_+\\) exactly cancels the diffusion current \\(-D\\nabla\\rho\\). Insert the two facts:

\\\[-\\frac{\\nabla V}{m\\gamma}=-D\\frac{\\nabla V}{k_BT}\\qquad\\Longrightarrow\\qquad D=\\frac{k_BT}{m\\gamma},\\\]

Einstein’s relation. So the osmotic velocity is the equilibrium drift, and I.3 is right where it can be tested.

*What is missing.* When \\(\\mathbf v\\neq0\\) the density flows, and to say how \\(\\mathbf v\\) evolves we need an energy and a law. The energy comes from the process (I.5); the law is one postulate (I.7); in between, I.6 shows that the continuity equation we already have is half of Hamilton’s equations for that energy, which fixes the bracket.

## I.5 Kinetic energy from the process

The kinetic energy is \\(\\frac m2\\times\\) (velocity)\\(^2\\), averaged over the noise, for a particle now at \\(\\mathbf x\\). From I.1 the velocity over the next step is \\((\\mathbf x\_{n+1}-\\mathbf x_n)/\\Delta t=\\mathbf b\_+(\\mathbf x)+\\sqrt{2D/\\Delta t}\\,\\boldsymbol\\xi_n\\). Square and average over \\(\\boldsymbol\\xi_n\\) at fixed \\(\\mathbf x_n=\\mathbf x\\); the cross term carries \\(\\langle\\boldsymbol\\xi\\rangle=0\\) and the square carries \\(\\langle|\\boldsymbol\\xi|^2\\rangle=d\\):

\\\[\\Big\\langle\\Big|\\frac{\\mathbf x\_{n+1}-\\mathbf x_n}{\\Delta t}\\Big|^2\\ \\Big|\\ \\mathbf x_n=\\mathbf x\\Big\\rangle=|\\mathbf b\_+(\\mathbf x)|^2+\\frac{2D\\,d}{\\Delta t}.\\\]

The second piece diverges as \\(\\Delta t\\to0\\), but it is the same number everywhere and at all times: an additive constant in the energy, which changes no equation of motion. Drop it. The same step read backward (I.3) gives \\(|\\mathbf b\_-|^2+2Dd/\\Delta t\\). The process read backward is the same process, so the energy must not depend on which way we read it: take the average of the two,

\\\[T=\\frac12\\Big\[\\frac m2|\\mathbf b\_+|^2+\\frac m2|\\mathbf b\_-|^2\\Big\]=\\frac m4\\big\[|\\mathbf v+\\mathbf u|^2+|\\mathbf v-\\mathbf u|^2\\big\]=\\boxed{\\frac m2\\big(v^2+u^2\\big)}\\\]

The osmotic velocity carries kinetic energy: a cloud that is not flat has energy \\(\\frac m2u^2\\) per particle even when nothing flows.

*This averaging is a rule, and here is what it decides.* Three combinations of \\(\\mathbf b\_+\\) and \\(\\mathbf b\_-\\) could be called “the kinetic energy”, and they are not equivalent:

- *Forward only,* \\(\\frac m2|\\mathbf b\_+|^2=\\frac m2|\\mathbf v+\\mathbf u|^2\\). This breaks the next step: with it, \\(\\delta H/\\delta S=-\\nabla\\cdot(\\rho\\,\\mathbf b\_+)\\), which is not \\(\\partial_t\\rho\\) (it is off by \\(D\\nabla^2\\rho\\)), so the continuity equation would not be Hamilton’s equation. Worse, it contains no quantum mechanics at all: \\(\\frac m2|\\mathbf v+\\mathbf u|^2=\\frac{|\\nabla S+mD\\nabla\\ln\\rho|^2}{2m}\\), and the shift \\(S'=S+mD\\ln\\rho\\) (canonical, since it leaves \\(\\{\\rho,S'\\}=\\delta\\)) turns it into \\(\\frac{|\\nabla S'|^2}{2m}\\): the classical fluid, with the osmotic velocity absorbed into the phase. The quantum potential exists only because \\(v^2+u^2\\) is *not* the square of a single gradient.
- *The product,* \\(\\frac m2\\,\\mathbf b\_+\\cdot\\mathbf b\_-=\\frac m2(v^2-u^2)\\). Time-symmetric too, and it also makes continuity Hamiltonian (the \\(u^2\\) term has no \\(S\\) in it either way). But it is not the mean square of any velocity — it is the product of two different ones — and it gives the quantum potential with the opposite sign, which does not linearize to Schrödinger (I.9 needs the \\(+\\)).
- *The average of the two energies,* \\(\\frac12\[\\frac m2|\\mathbf b\_+|^2+\\frac m2|\\mathbf b\_-|^2\]=\\frac m2(v^2+u^2)\\): the mean square velocity, averaged over the two readings of the process. This is the one used, and it is the one input of I.5. Time symmetry alone does not select it; “the energy is the average of the energies” does.

## I.6 The energy as a functional, and the bracket

The energy of the whole cloud is the density-weighted sum of kinetic and potential energy, \\(H=\\int\\rho\\,\[T+V\]\\,d^dx\\) with \\(T\\) from I.5. \\(\\mathbf u\\) is already a function of \\(\\rho\\). To make \\(H\\) a functional of fields we need \\(\\mathbf v\\) as a field too:

*Assumption (irrotational flow):* \\(\\mathbf v=\\nabla S/m\\), with \\(S\\) a scalar field (the velocity potential, in units of an action). Then

\\\[H\[\\rho,S\]=\\int\\rho\\Big\[\\frac m2(v^2+u^2)+V\\Big\]d^dx=\\int\\Big\[\\frac{\\rho|\\nabla S|^2}{2m}+\\frac{mD^2}{2}\\frac{|\\nabla\\rho|^2}{\\rho}+V\\rho\\Big\]d^dx,\\qquad\\frac{mD^2}2=\\frac{\\hbar^2}{8m}.\\\]

Now vary \\(H\\) with respect to \\(S\\). Integrating by parts once, \\(\\delta H/\\delta S=-\\nabla\\cdot(\\rho\\nabla S/m)=-\\nabla\\cdot(\\rho\\mathbf v)\\), and by I.3 that is \\(\\partial_t\\rho\\). So the continuity equation we already have reads

\\\[\\partial_t\\rho=\\frac{\\delta H}{\\delta S}.\\\]

That is Hamilton’s equation for \\(\\rho\\) when \\(\\rho\\) and \\(S\\) are a canonical pair,

\\\[\\boxed{\\{\\rho(\\mathbf x),S(\\mathbf y)\\}=\\delta(\\mathbf x-\\mathbf y)}\\\]

This bracket is the canonical one that makes the derived continuity equation Hamiltonian with the derived energy — a choice among Poisson structures, the standard one, not forced by continuity alone (one component of a vector field does not fix a bracket). \\(S\\) is the momentum conjugate to \\(\\rho\\).

## I.7 The dynamical postulate and Hamilton–Jacobi

A canonical pair has two Hamilton equations. One is the continuity equation, already derived. The other is the only thing in Part I that is put in:

*Assumption (the one dynamical input):* \\(\\dot S=-\\delta H/\\delta\\rho\\). In words: the energy \\(H\\) generates the time evolution of both fields.

Compute the derivative. For a density \\(F(\\rho,\\nabla\\rho)\\), \\(\\frac{\\delta}{\\delta\\rho}\\int F=\\frac{\\partial F}{\\partial\\rho}-\\nabla\\cdot\\frac{\\partial F}{\\partial(\\nabla\\rho)}\\) (the second term is the integration by parts that moves \\(\\nabla\\) off \\(\\delta\\rho\\)). For the osmotic term \\(F=|\\nabla\\rho|^2/\\rho\\):

\\\[\\frac{\\partial F}{\\partial\\rho}=-\\frac{|\\nabla\\rho|^2}{\\rho^2},\\qquad\\nabla\\cdot\\frac{\\partial F}{\\partial(\\nabla\\rho)}=\\nabla\\cdot\\frac{2\\nabla\\rho}{\\rho}=\\frac{2\\nabla^2\\rho}{\\rho}-\\frac{2|\\nabla\\rho|^2}{\\rho^2},\\\]

\\\[\\frac{\\delta}{\\delta\\rho}\\int\\frac{|\\nabla\\rho|^2}\\rho=\\frac{|\\nabla\\rho|^2}{\\rho^2}-\\frac{2\\nabla^2\\rho}{\\rho}=-4\\frac{\\nabla^2\\sqrt\\rho}{\\sqrt\\rho},\\\]

the last step by \\(\\nabla^2\\sqrt\\rho=\\nabla\\cdot\\frac{\\nabla\\rho}{2\\sqrt\\rho}=\\frac{\\nabla^2\\rho}{2\\sqrt\\rho}-\\frac{|\\nabla\\rho|^2}{4\\rho^{3/2}}\\). So

\\\[\\frac{\\delta H}{\\delta\\rho}=\\frac{|\\nabla S|^2}{2m}+V-2mD^2\\frac{\\nabla^2\\sqrt\\rho}{\\sqrt\\rho},\\qquad2mD^2=\\frac{\\hbar^2}{2m},\\\]

\\\[\\boxed{\\partial_tS+\\frac{|\\nabla S|^2}{2m}+V-\\frac{\\hbar^2}{2m}\\frac{\\nabla^2\\sqrt\\rho}{\\sqrt\\rho}=0}\\\]

Hamilton–Jacobi with one extra term, the quantum potential \\(Q=-\\frac{\\hbar^2}{2m}\\frac{\\nabla^2\\sqrt\\rho}{\\sqrt\\rho}\\). It is the osmotic kinetic energy \\(\\frac m2u^2\\) of I.5, after the variation with respect to \\(\\rho\\). Nothing quantum was put into this term: \\(Q\\) is there because the random walk carries kinetic energy in its spreading. (What is quantum, and is put in, is the single-valuedness of \\(\\psi\\) — Wallstrom’s point, in the ledger.)

## I.8 The same law as Nelson wrote it: mean acceleration

Nelson phrased the dynamical input as Newton’s law for a suitably defined acceleration. It is the same input; here is the dictionary. For any function \\(f(\\mathbf x,t)\\), the mean rate at which it changes over the next step, for a particle now at \\(\\mathbf x\\), is

\\\[D\_+f(\\mathbf x,t)\\equiv\\lim\_{\\Delta t\\to0}\\frac{\\big\\langle f(\\mathbf x\_{n+1},t+\\Delta t)-f(\\mathbf x,t)\\ \\big|\\ \\mathbf x_n=\\mathbf x\\big\\rangle}{\\Delta t}.\\\]

Compute it exactly as in I.2: Taylor-expand \\(f(\\mathbf x+\\mathbf z,t+\\Delta t)\\) in the step \\(\\mathbf z=\\mathbf x\_{n+1}-\\mathbf x_n\\) and in \\(\\Delta t\\), and use the two averages of I.1, \\(\\langle z_i\\rangle=b\_{+,i}\\Delta t\\) and \\(\\langle z_iz_j\\rangle=2D\\Delta t\\,\\delta\_{ij}\\):

\\\[\\langle f(\\mathbf x+\\mathbf z,t+\\Delta t)\\rangle=f+\\Delta t\\,\\partial_tf+\\langle z_i\\rangle\\partial_if+\\tfrac12\\langle z_iz_j\\rangle\\partial_i\\partial_jf+\\ldots\\qquad\\Longrightarrow\\qquad D\_+f=\\partial_tf+\\mathbf b\_+\\cdot\\nabla f+D\\nabla^2f .\\\]

The same over the previous step, read backward from \\(\\mathbf x\\) (I.3: the step \\(\\mathbf x\_{n-1}-\\mathbf x_n\\) has mean \\(-\\mathbf b\_-\\Delta t\\) and the same second moment),

\\\[D\_-f(\\mathbf x,t)\\equiv\\lim\_{\\Delta t\\to0}\\frac{\\big\\langle f(\\mathbf x,t)-f(\\mathbf x\_{n-1},t-\\Delta t)\\ \\big|\\ \\mathbf x_n=\\mathbf x\\big\\rangle}{\\Delta t}=\\partial_tf+\\mathbf b\_-\\cdot\\nabla f-D\\nabla^2f ,\\\]

the second-order term with the opposite sign because now the expansion is subtracted.

Check: \\(D\_+\\mathbf x=\\mathbf b\_+\\), \\(D\_-\\mathbf x=\\mathbf b\_-\\). Nelson’s time-symmetric acceleration is \\(\\mathbf a=\\tfrac12(D\_+D\_-+D\_-D\_+)\\mathbf x=\\tfrac12(D\_+\\mathbf b\_-+D\_-\\mathbf b\_+)\\). Expand with \\(\\mathbf b\_\\pm=\\mathbf v\\pm\\mathbf u\\):

\\\[\\mathbf a=\\partial_t\\mathbf v+\\tfrac12\\big\[(\\mathbf v+\\mathbf u)\\cdot\\nabla(\\mathbf v-\\mathbf u)+(\\mathbf v-\\mathbf u)\\cdot\\nabla(\\mathbf v+\\mathbf u)\\big\]+\\tfrac D2\\nabla^2(\\mathbf b\_--\\mathbf b\_+)=\\partial_t\\mathbf v+(\\mathbf v\\cdot\\nabla)\\mathbf v-(\\mathbf u\\cdot\\nabla)\\mathbf u-D\\nabla^2\\mathbf u.\\\]

With \\(\\mathbf v=\\nabla S/m\\) and \\(\\mathbf u=D\\nabla\\ln\\rho\\) both gradients, \\((\\mathbf v\\cdot\\nabla)\\mathbf v=\\nabla\\frac{v^2}2\\), \\((\\mathbf u\\cdot\\nabla)\\mathbf u=\\nabla\\frac{u^2}2\\), \\(\\nabla^2\\mathbf u=\\nabla(\\nabla\\cdot\\mathbf u)\\), and

\\\[\\frac{u^2}2+D\\,\\nabla\\cdot\\mathbf u=2D^2\\frac{\\nabla^2\\sqrt\\rho}{\\sqrt\\rho},\\\]

so \\(m\\mathbf a=\\nabla\\big\[\\frac{|\\nabla S|^2}{2m}+\\partial_tS\\big\]-\\nabla\\big\[2mD^2\\frac{\\nabla^2\\sqrt\\rho}{\\sqrt\\rho}\\big\]\\), and Newton’s law \\(m\\mathbf a=-\\nabla V\\) is the gradient of the boxed equation of I.7. Newton’s law for Nelson’s acceleration and the postulate \\(\\dot S=-\\delta H/\\delta\\rho\\) are the same statement.

## I.9 Madelung: the canonical transformation to \\(\\psi\\)

Set \\(\\psi=\\sqrt\\rho\\,e^{iS/\\hbar}\\), \\(\\psi^\*=\\sqrt\\rho\\,e^{-iS/\\hbar}\\). The bracket transforms as a Jacobian:

\\\[\\{\\psi(\\mathbf x),\\psi^\*(\\mathbf y)\\}=\\Big\[\\frac{\\partial\\psi}{\\partial\\rho}\\frac{\\partial\\psi^\*}{\\partial S}-\\frac{\\partial\\psi}{\\partial S}\\frac{\\partial\\psi^\*}{\\partial\\rho}\\Big\]\\{\\rho(\\mathbf x),S(\\mathbf y)\\}=\\Big\[\\frac\\psi{2\\rho}\\cdot\\frac{-i\\psi^\*}\\hbar-\\frac{i\\psi}\\hbar\\cdot\\frac{\\psi^\*}{2\\rho}\\Big\]\\delta(\\mathbf x-\\mathbf y)=\\boxed{-\\frac i\\hbar\\,\\delta(\\mathbf x-\\mathbf y)}\\\]

using \\(|\\psi|^2=\\rho\\); the square bracket is evaluated at \\(\\mathbf x\\), and the \\(\\delta\\) is the one of I.6. The Hamiltonian transforms by the identity \\(\\nabla\\psi=\\psi\\big(\\frac{\\nabla\\rho}{2\\rho}+\\frac i\\hbar\\nabla S\\big)\\), so \\(|\\nabla\\psi|^2=\\rho\\big\[\\frac{|\\nabla\\rho|^2}{4\\rho^2}+\\frac{|\\nabla S|^2}{\\hbar^2}\\big\]\\) and

\\\[\\frac{\\hbar^2}{2m}|\\nabla\\psi|^2=\\frac{\\hbar^2}{8m}\\frac{|\\nabla\\rho|^2}\\rho+\\frac{\\rho|\\nabla S|^2}{2m},\\qquad H=\\int\\Big\[\\frac{\\hbar^2}{2m}|\\nabla\\psi|^2+V|\\psi|^2\\Big\]d^dx.\\\]

The osmotic term and the kinetic term fuse into one \\(|\\nabla\\psi|^2\\). Hamilton’s equation:

\\\[\\dot\\psi=\\{\\psi,H\\}=\\int\\{\\psi(\\mathbf x),\\psi^\*(\\mathbf x')\\}\\frac{\\delta H}{\\delta\\psi^\*(\\mathbf x')}d^dx'=-\\frac i\\hbar\\frac{\\delta H}{\\delta\\psi^\*}=-\\frac i\\hbar\\Big\[-\\frac{\\hbar^2}{2m}\\nabla^2\\psi+V\\psi\\Big\],\\\]

\\\[\\boxed{i\\hbar\\,\\partial_t\\psi=-\\frac{\\hbar^2}{2m}\\nabla^2\\psi+V\\psi}\\\]

The \\(i\\) is the bracket. One comparison for later. The classical particle with no drift obeys FP\\(\_+\\) with \\(\\mathbf b\_+=0\\): \\(\\partial_t\\rho=D\\nabla^2\\rho\\), the heat equation. The free quantum particle obeys the equation just derived with \\(V=0\\): \\(\\partial_t\\psi=\\frac{i\\hbar}{2m}\\nabla^2\\psi=iD\\nabla^2\\psi\\). Same equation, \\(D\\) replaced by \\(iD\\). **Heat to Schrödinger is \\(D\\to iD\\).** Part II meets exactly this; Part III finds the process in which the \\(i\\) is already there.

### What went in (Part I)

| FP\\(\_\\pm\\), \\(\\mathbf u=D\\nabla\\ln\\rho\\), continuity | derived from I.1 |
| --- | --- |
| Einstein’s relation | check of I.3, using two thermal facts |
| \\(T=\\frac m2(v^2+u^2)\\) | from I.1, with one rule: the average of the forward and backward kinetic energies (I.5 says what the other choices do) |
| \\(\\{\\rho,S\\}=\\delta\\) | the canonical bracket that represents continuity with this \\(H\\); a choice (I.6) |
| \\(\\mathbf v=\\nabla S/m\\) | assumption |
| \\(\\dot S=-\\delta H/\\delta\\rho\\) (= Newton for Nelson’s acceleration) | the one dynamical input |
| \\(D=\\hbar/2m\\) | sets the scale |
| single-valued \\(\\psi\\) | needs circulation quantization (Wallstrom); an axiom |

---

# Part II — The random walk in spacetime → Klein–Gordon

## II.0 Why a second process

In Part I the clock ticks and the particle steps in space. Time is the clock, not a direction the particle can step in. Special relativity says time and space are directions of one spacetime and that no inertial observer is preferred. A relativistic version of Part I must therefore let the particle step in time as well as in space — and then the clock that counts the steps cannot be \\(t\\) any more. The clock is a new label, \\(\\tau\\), the number of steps taken; nothing observable is allowed to depend on it (II.7). Everything else is Part I, repeated with \\(x\\) replaced by the pair \\((t,x)\\), and one sign that relativity forces.

## II.1 The process

Coordinates \\(x^\\mu=(x^0,x^1)=(ct,x)\\). One step:

\\\[x^\\mu\_{n+1}=x^\\mu_n+b\_+^\\mu(x_n)\\,\\Delta\\tau+\\sqrt{2D\\,\\Delta\\tau}\\;\\xi^\\mu_n,\\qquad\\langle\\xi^\\mu\\rangle=0,\\qquad\\langle\\xi^\\mu\\xi^\\nu\\rangle=\\eta^{\\mu\\nu}=\\begin{pmatrix}-1&0\\\\0&1\\end{pmatrix}.\\\]

*Why \\(\\eta^{\\mu\\nu}\\).* In Part I the kick had the same variance in every space direction, so that no direction was preferred. Here no inertial frame may be preferred: the variance matrix \\(\\langle\\xi^\\mu\\xi^\\nu\\rangle\\) must look the same after any boost. A boost mixes \\(x^0\\) and \\(x^1\\) and preserves exactly one symmetric matrix up to a factor, \\(\\eta^{\\mu\\nu}\\) (the same matrix that makes \\(-c^2t^2+x^2\\) invariant). So the variance matrix is a multiple of \\(\\eta\\), and the factor is absorbed in \\(D\\).

*What it says.* \\(\\langle\\xi^1\\xi^1\\rangle=+1\\): the space kick is Part I’s. \\(\\langle\\xi^0\\xi^0\\rangle=-1\\): **the time kick has negative variance.** No real random number does that; \\(\\xi^0=i\\zeta\\) with \\(\\zeta\\) an ordinary Gaussian. This is the timelike noise. It means the process cannot be built with coins in the time direction — the time direction is Part I’s walk continued to imaginary steps, which is relativity’s way of saying time is not a space direction. Everything that follows is algebra on the two averages \\(\\langle z^\\mu\\rangle=b^\\mu\\Delta\\tau\\) and \\(\\langle z^\\mu z^\\nu\\rangle=2D\\Delta\\tau\\,\\eta^{\\mu\\nu}\\) of the step \\(z^\\mu=x^\\mu\_{n+1}-x^\\mu_n\\); I.2 showed that is all a Fokker–Planck derivation uses. Where the sign of \\(\\eta^{00}\\) matters, it is said.

*What is real in this walk.* The time kick is imaginary, so Part II as written is a formal walk. What is real is III.16’s object: the occupation density of a positive walk killed at the rate \\(m^2\\), in the signature where the kick is real, which is the source answer \\(K_0\\), \\(K_1\\). The Lorentzian walk is that density continued across the light cone — the same move III.15 makes for the time-sense walk’s kernel. Part II and Part III pay the Lorentzian sign once each, in the signature and in the clock.

Indices: \\(\\partial\_\\mu=(\\partial_0,\\partial_1)=(\\frac1c\\partial_t,\\partial_x)\\); raising with \\(\\eta\\), \\(\\partial^0=-\\partial_0\\), \\(\\partial^1=\\partial_1\\); and for any two vectors \\(a\\cdot b\\equiv\\eta\_{\\mu\\nu}a^\\mu b^\\nu=-a^0b^0+a^1b^1\\). The wave operator is \\(\\Box\\equiv\\partial^\\mu\\partial\_\\mu=-\\frac1{c^2}\\partial_t^2+\\partial_x^2\\).

## II.2 Fokker–Planck in spacetime

The density is now a density on spacetime at a given step: \\(\\rho(x^\\mu,\\tau)\\,d^2x\\) is the probability of finding the walker in the cell \\(d^2x=c\\,dt\\,dx\\) around \\(x^\\mu=(ct,x)\\) after \\(\\tau\\) worth of steps. (With the imaginary time kick of II.1 this is a weight rather than a probability in the \\(t\\) direction; II.9 returns to what is and is not positive.) I.2 word for word: to be at \\(x^\\mu\\) at step \\(n+1\\) you were at \\(x^\\mu-z^\\mu\\) and stepped by \\(z^\\mu\\); expand in the small step; the derivatives come out of the average; the two averages of II.1 go in:

\\\[\\rho(x,\\tau+\\Delta\\tau)=\\rho-\\partial\_\\mu\\big(\\rho\\,b\_+^\\mu\\big)\\Delta\\tau+\\tfrac12\\partial\_\\mu\\partial\_\\nu\\big(\\rho\\,2D\\Delta\\tau\\,\\eta^{\\mu\\nu}\\big)+\\ldots\\qquad\\Longrightarrow\\qquad\\boxed{\\partial\_\\tau\\rho=-\\partial\_\\mu(\\rho\\,b\_+^\\mu)+D\\,\\Box\\rho}\\qquad(\\mathrm{FP}\_+)\\\]

Written out, \\(D\\Box\\rho=-\\frac{D}{c^2}\\partial_t^2\\rho+D\\partial_x^2\\rho\\). In \\(x\\) the density spreads as in Part I. In \\(t\\) it *un*-spreads: the sign of \\(\\eta^{00}\\) has put a minus in front of \\(\\partial_t^2\\). The wave operator is on the page before anything quantum has been done; it is what the timelike noise produces.

*The same step, said once more.* “Where you were, times the step, then Taylor in the step” is the average of the shift \\(\\rho(x-z)=e^{-z\\cdot\\partial}\\rho(x)\\) over the step distribution. Its first moment is the drift, its second moment is \\(D\\Box\\): Klein–Gordon is the second moment of the shift. Part III will keep the first moment instead, with the direction of the step carried as a state rather than averaged — that is the whole difference between the two Parts, in one sentence.

## II.3 The same process run backward

I.3 again. Looking back from \\(x^\\mu\\) at step \\(n+1\\), the step \\(x^\\mu_n-x^\\mu\_{n+1}\\) has a Gaussian weight whose exponent is \\(-\\eta\_{\\mu\\nu}(z+b\_+\\Delta\\tau)^\\mu(z+b\_+\\Delta\\tau)^\\nu/4D\\Delta\\tau\\) (the inverse of the variance matrix \\(\\eta^{\\mu\\nu}\\) is \\(\\eta\_{\\mu\\nu}\\)) times \\(\\rho(x+z)\\approx\\rho(x)e^{+\\partial\_\\mu\\ln\\rho\\,z^\\mu}\\) (the previous position is \\(x+z\\) for this \\(z\\)). Complete the square exactly as in I.3 — the only change is that the shift of the centre picks up a raised index, because \\(\\eta\_{\\mu\\nu}w^\\mu w^\\nu+2g\_\\mu w^\\mu=\\eta\_{\\mu\\nu}(w+g^\\sharp)^\\mu(w+g^\\sharp)^\\nu-g\\cdot g\\) with \\((g^\\sharp)^\\mu=\\eta^{\\mu\\nu}g\_\\nu\\):

\\\[\\boxed{b\_-^\\mu=b\_+^\\mu-2D\\,\\partial^\\mu\\ln\\rho}\\\]

The reversed process is of the same kind, with FP\\(\_-\\): \\(\\partial\_\\tau\\rho=-\\partial\_\\mu(\\rho b\_-^\\mu)-D\\Box\\rho\\) — and substituting \\(b\_-\\) turns it back into FP\\(\_+\\) exactly as in I.3 (with \\(\\rho\\,\\partial^\\mu\\ln\\rho=\\partial^\\mu\\rho\\)). Add and halve, subtract and halve:

\\\[\\partial\_\\tau\\rho+\\partial\_\\mu(\\rho\\,v^\\mu)=0,\\qquad v^\\mu=\\tfrac12(b\_+^\\mu+b\_-^\\mu),\\qquad u^\\mu=\\tfrac12(b\_+^\\mu-b\_-^\\mu)=D\\,\\partial^\\mu\\ln\\rho .\\\]

Continuity in the step label \\(\\tau\\), with a current velocity \\(v^\\mu\\) and an osmotic velocity \\(u^\\mu\\) that now have a time component: \\(u^0=-D\\,\\partial_0\\ln\\rho\\). The osmotic velocity in the time direction is the un-spreading of II.2.

## II.4 Kinetic term from the process

I.5 with the dot product of II.1. The velocity over the next step is \\((x^\\mu\_{n+1}-x^\\mu_n)/\\Delta\\tau=b\_+^\\mu(x)+\\sqrt{2D/\\Delta\\tau}\\,\\xi^\\mu_n\\); its invariant square, averaged over \\(\\xi_n\\) at fixed \\(x_n=x\\), is

\\\[\\Big\\langle\\frac{x\_{n+1}-x_n}{\\Delta\\tau}\\cdot\\frac{x\_{n+1}-x_n}{\\Delta\\tau}\\ \\Big|\\ x_n=x\\Big\\rangle=b\_+\\cdot b\_++\\frac{2D}{\\Delta\\tau}\\,\\eta\_{\\mu\\nu}\\eta^{\\mu\\nu}=b\_+\\cdot b\_++\\frac{4D}{\\Delta\\tau},\\\]

the last term a state-independent constant, dropped. Average forward and backward as in I.5:

\\\[T=\\frac12\\Big\[\\frac m2\\,b\_+\\cdot b\_++\\frac m2\\,b\_-\\cdot b\_-\\Big\]=\\boxed{\\frac m2\\big(v\\cdot v+u\\cdot u\\big)}\\\]

The dot products carry the sign of \\(\\eta\\): \\(v\\cdot v=-(v^0)^2+(v^1)^2\\). For a slow particle \\(v^0\\approx c\\) and \\(v\\cdot v\\approx-c^2\\); this negative constant is where the rest energy will sit.

## II.5 The functional and the bracket

I.6. *Assumption (irrotational):* \\(v^\\mu=\\partial^\\mu S/m\\), with \\(S\\) a scalar field on spacetime. Then

\\\[H\[\\rho,S\]=\\int\\rho\\,\\frac m2\\big(v\\cdot v+u\\cdot u\\big)\\,d^2x=\\int\\Big\[\\frac{\\rho\\;\\partial^\\mu S\\,\\partial\_\\mu S}{2m}+\\frac{\\hbar^2}{8m}\\frac{\\partial^\\mu\\rho\\,\\partial\_\\mu\\rho}{\\rho}\\Big\]d^2x,\\qquad d^2x=c\\,dt\\,dx,\\\]

with \\(mD^2/2=\\hbar^2/8m\\) as before. Vary with respect to \\(S\\): \\(\\delta H/\\delta S=-\\partial\_\\mu(\\rho\\,\\partial^\\mu S/m)=-\\partial\_\\mu(\\rho v^\\mu)\\), which by II.3 is \\(\\partial\_\\tau\\rho\\). So continuity is Hamilton’s equation \\(\\partial\_\\tau\\rho=\\delta H/\\delta S\\) for

\\\[\\boxed{\\{\\rho(x),S(x')\\}=\\delta^2(x-x')}\\\]

the same bracket as I.6, now on spacetime, with \\(\\tau\\) as the evolution parameter. (No potential here: the free particle is the case treated.)

## II.6 The dynamical postulate

I.7: *Assumption (the one dynamical input):* \\(\\partial\_\\tau S=-\\delta H/\\delta\\rho\\). The functional derivative is I.7’s with \\(\\nabla\\to\\partial\_\\mu\\) and the \\(\\eta\\) signs carried through; the osmotic term gives \\(\\frac{\\delta}{\\delta\\rho}\\int\\frac{\\hbar^2}{8m}\\frac{\\partial\\rho\\cdot\\partial\\rho}{\\rho}=-\\frac{\\hbar^2}{2m}\\frac{\\Box\\sqrt\\rho}{\\sqrt\\rho}\\) (same algebra as I.7, checked). So

\\\[\\boxed{\\partial\_\\tau S+\\frac{\\partial^\\mu S\\,\\partial\_\\mu S}{2m}-\\frac{\\hbar^2}{2m}\\frac{\\Box\\sqrt\\rho}{\\sqrt\\rho}=0}\\\]

Hamilton–Jacobi in \\(\\tau\\), with the quantum potential built from the wave operator.

## II.7 The mass: nothing may depend on the step label

\\(\\tau\\) counts steps. It is not a clock anyone reads, so the physical solutions are the ones in which nothing observable depends on it: \\(\\rho=\\rho(x^\\mu)\\) only, and \\(S\\) advancing uniformly, \\(S(x,\\tau)=\\tilde S(x)-K\\tau\\) with \\(K\\) a constant. Put this into II.6 and II.3:

\\\[-K+\\frac{\\partial^\\mu\\tilde S\\,\\partial\_\\mu\\tilde S}{2m}-\\frac{\\hbar^2}{2m}\\frac{\\Box\\sqrt\\rho}{\\sqrt\\rho}=0,\\qquad\\partial\_\\mu(\\rho\\,\\partial^\\mu\\tilde S)=0.\\\]

What is \\(K\\)? Switch the noise off (\\(\\hbar\\to0\\)) and look at a particle with definite momentum: \\(\\partial\_\\mu\\tilde S=p\_\\mu=(-E/c,\\,p)\\) (the phase advances as \\(px-Et\\)), so \\(\\partial^\\mu\\tilde S\\,\\partial\_\\mu\\tilde S=p\\cdot p=-E^2/c^2+p^2\\). For a particle of mass \\(m\\) this is \\(-m^2c^2\\) — that is what “mass \\(m\\)” means in relativity, \\(E^2=p^2c^2+m^2c^4\\). So \\(2mK=-m^2c^2\\):

\\\[K=-\\frac{mc^2}{2},\\qquad\\boxed{\\partial^\\mu\\tilde S\\,\\partial\_\\mu\\tilde S+m^2c^2=\\hbar^2\\frac{\\Box\\sqrt\\rho}{\\sqrt\\rho},\\qquad\\partial\_\\mu(\\rho\\,\\partial^\\mu\\tilde S)=0}\\\]

The relativistic Hamilton–Jacobi equation with the quantum potential, and the conservation of the current \\(\\rho\\,\\partial^\\mu\\tilde S\\). **The mass entered here and only here:** as the value of the step-label energy \\(K\\) that makes the noiseless solutions obey \\(E^2=p^2c^2+m^2c^4\\). In Part I there was no such step because \\(t\\) was the clock; here the clock is \\(\\tau\\), and fixing what a unit of \\(\\tau\\) costs is fixing the mass.

## II.8 Madelung: Klein–Gordon

I.9: \\(\\psi=\\sqrt\\rho\\,e^{i\\tilde S/\\hbar}\\), bracket \\(\\{\\psi,\\psi^\*\\}=-i/\\hbar\\) by the same Jacobian. The identity of I.9 with \\(\\nabla\\to\\partial\_\\mu\\):

\\\[\\frac{\\Box\\psi}{\\psi}=\\frac{\\Box\\sqrt\\rho}{\\sqrt\\rho}-\\frac{\\partial^\\mu\\tilde S\\,\\partial\_\\mu\\tilde S}{\\hbar^2}+\\frac i\\hbar\\,\\frac{\\partial\_\\mu(\\rho\\,\\partial^\\mu\\tilde S)}{\\rho}\\\]

(multiply out \\(\\partial\_\\mu\\partial^\\mu(\\sqrt\\rho\\,e^{i\\tilde S/\\hbar})\\); checked). Insert the two boxed equations of II.7: the imaginary part is zero and the real part is \\(m^2c^2/\\hbar^2\\):

\\\[\\boxed{\\Big(\\Box-\\frac{m^2c^2}{\\hbar^2}\\Big)\\psi=0\\qquad\\Longleftrightarrow\\qquad\\partial_t^2\\psi=c^2\\partial_x^2\\psi-\\frac{m^2c^4}{\\hbar^2}\\psi}\\\]

The Klein–Gordon equation. Equivalently, before imposing II.7, the full \\(\\tau\\)-dependent field \\(\\Psi=\\sqrt\\rho\\,e^{iS/\\hbar}=e^{imc^2\\tau/2\\hbar}\\psi\\) obeys Hamilton’s equation with the \\(H\\) of II.5, \\(i\\hbar\\,\\partial\_\\tau\\Psi=-\\frac{\\hbar^2}{2m}\\Box\\Psi\\): Schrödinger’s equation with \\(\\tau\\) for time and \\(\\Box\\) for \\(\\nabla^2\\). II.7 is the statement that \\(\\Psi\\) is a stationary state of it with energy \\(K=-mc^2/2\\).

## II.9 What came out, and the slow limit

*The slow limit is Part I.* Write \\(\\psi=e^{-imc^2t/\\hbar}\\phi\\) with \\(\\phi\\) slowly varying (\\(|\\partial_t\\phi|\\ll mc^2|\\phi|/\\hbar\\)). Then \\(\\partial_t^2\\psi=e^{-imc^2t/\\hbar}\\big\[-\\frac{m^2c^4}{\\hbar^2}\\phi-\\frac{2imc^2}{\\hbar}\\partial_t\\phi+\\partial_t^2\\phi\\big\]\\); the first term cancels the mass term of Klein–Gordon, the last is dropped, and what is left is \\(i\\hbar\\,\\partial_t\\phi=-\\frac{\\hbar^2}{2m}\\partial_x^2\\phi\\): Part I’s result.

*Where the timelike noise is, in the end.* The only place \\(\\eta^{00}=-1\\) acted was to put the minus sign in front of \\(\\partial_t^2\\) in \\(\\Box\\), in II.2. Everything after that was Part I’s algebra with that operator. The mass came from II.7, not from the noise.

*What is and is not a probability.* \\(\\rho=|\\psi|^2\\) is positive. But the conserved quantity of II.7, \\(\\partial\_\\mu(\\rho\\,\\partial^\\mu\\tilde S)=0\\), has time component \\(\\rho\\,\\partial^0\\tilde S=\\rho E/c\\) — that is \\(m\\rho v^0\\), with \\(v^0=E/mc\\) — with \\(E\\) of either sign, and it is this, not \\(\\rho\\), that is conserved in \\(t\\). That is the Klein–Gordon “negative probability” problem, and the same Lorentzian sign that was in \\(\\langle\\xi^0\\xi^0\\rangle=-1\\) is behind it. It is left as it is here; Part III meets it again, as the conserved charge \\(P\_+-P\_-\\) of a classical zigzag — forward-time crossings minus backward-time crossings — which is where the indefinite sign gets its meaning (pair creation).

### What went in (Part II)

| steps in spacetime; variance matrix \\(\\propto\\eta^{\\mu\\nu}\\) | relativity: no preferred frame. This is the timelike noise. |
| --- | --- |
| FP\\(\_\\pm\\) with \\(\\Box\\); \\(b\_-\\); \\(u^\\mu=D\\partial^\\mu\\ln\\rho\\); continuity in \\(\\tau\\) | derived, as in I.2–I.3 |
| \\(T=\\frac m2(v\\cdot v+u\\cdot u)\\) | derived, as in I.5 |
| \\(\\{\\rho,S\\}=\\delta\\) | the canonical bracket, as in I.6 |
| \\(v^\\mu=\\partial^\\mu S/m\\) | assumption, as in I.6 |
| \\(\\partial\_\\tau S=-\\delta H/\\delta\\rho\\) | the one dynamical input, as in I.7 |
| \\(\\tau\\)-independence with \\(K=-mc^2/2\\) | the mass shell: this is where \\(m\\) enters |
| Klein–Gordon; Schrödinger as its slow limit | derived |

---

# Part III — The time-sense particle → Dirac (one dimension)

## III.0 Why a third process, and which one

Part II gives Klein–Gordon, which is second order in time. Dirac’s equation is first order. No process of the Part I/Part II kind can give a first-order equation: in both, the Hamiltonian came from \\(\\langle(\\text{velocity})^2\\rangle\\) and was quadratic in \\(\\partial S\\), and a quadratic Hamiltonian gives an equation of second order in the spatial derivative — and, once time is a direction the particle steps in (Part II), in \\(\\partial_t\\) as well. A first-order equation needs a process in which the velocity is a *state variable* rather than a derived quantity — a particle that carries its direction along (Dudley’s theorem says a relativistic Markov process must be of this kind).

There are two such particles in one space dimension, and they are not the same.

- **The telegraph particle (T).** \\(t\\) is the clock; the particle moves at \\(\\pm c\\) in \\(x\\) and reverses at random. Its worldline zigzags *in \\(x\\)*, going up: one position for each time.
- **The time-sense particle (S).** \\(x\\) is the clock; the particle always advances in \\(x\\) at \\(c\\), and what flips at random is the sense in which *time* runs along its path. Its worldline zigzags *in \\(t\\)*, going sideways: several positions for one time.

![One worldline of each](two_processes.png)

One worldline of each

Both have light-like segments and a random flip; the difference is which of the two coordinates carries the randomness. That difference decides the sign of the mass term (III.6): T gives the Klein–Gordon equation with the wrong sign and needs its flip rate continued to an imaginary value; **S gives Klein–Gordon with the right sign, and the Dirac equation, with a real flip rate and nothing continued.** That a walk which also runs backward in time reaches Dirac without the continuation is McKeon and Ord’s thesis (*Phys. Rev. Lett.* 69, 3 (1992): “when a particle always moves forward in time, an analytic continuation is required”); what is built here on that thesis is said in III.9’s remark. Part II already said what the time direction has to do: carry the Lorentzian sign. S is the process in which it does so with a coin instead of an imaginary variance. Part III is therefore built on S; T is kept as a remark wherever the two part ways, because it is the better-known process and because comparing them is what shows where the \\(i\\) of quantum mechanics comes from.

One more thing before starting. In one space dimension a particle can only go one way or the other; there is no plane for it to circulate in and no axis for a vortex. Part III is the *flattened* case; spin, the vortex, the plane seen as a flatlander, all need three dimensions and are Part IV.

## III.1 The process

State: a time coordinate \\(t\\) and a time sense \\(\\varepsilon=\\pm1\\). Clock: the position \\(x\\). In each \\(\\Delta x\\):

\\\[t\_{n+1}=t_n+\\varepsilon_n\\,\\frac{\\Delta x}{c},\\qquad\\varepsilon\_{n+1}=\\begin{cases}-\\varepsilon_n&\\text{with probability }\\dfrac{a}{c}\\,\\Delta x\\\\\[4pt\]\\ \\ \\varepsilon_n&\\text{with probability }1-\\dfrac ac\\,\\Delta x\\end{cases}\\\]

the flips independent of everything earlier. The particle always advances in \\(x\\) at speed \\(c\\); the one random thing is whether, as it does, \\(t\\) increases (\\(\\varepsilon=+1\\)) or decreases (\\(\\varepsilon=-1\\)). Two constants: \\(c\\) and the flip rate \\(a\\) — a rate per unit time along the path, \\(a/c\\) per unit \\(x\\). (Equivalently, with the path parameter \\(\\tau\\) of Part II: \\(dx/d\\tau=c\\), \\(dt/d\\tau=\\varepsilon\\), flips at rate \\(a\\) per unit \\(\\tau\\); since \\(x=c\\tau\\) along the path, \\(x\\) and \\(\\tau\\) are the same clock.)

*What a forward-time observer sees* (right panel of the figure). He reads the drawing upward. A climbing segment is a particle moving right; a descending segment, read upward, is a particle moving right *backward in time* — Feynman’s antiparticle. Where the zigzag turns, a pair appears (a trough) or disappears (a crest). At one time \\(t\\) the worldline is crossed several times: he sees several particles at once, though there is one path. That is not a defect of the description; it is the content of III.3.

*Remark: the telegraph particle T.* Make \\(t\\) the clock instead: \\(x\_{n+1}=x_n+\\sigma_nc\\Delta t\\), \\(\\sigma\\) flipping with probability \\(a\\Delta t\\). One position per time, right-movers and left-movers. Everything below can be redone for T by exchanging the roles of \\(x\\) and \\(t\\); where the result differs it is said.

Compare I.1 and I.1a: there the random part was a displacement of size \\(\\sqrt{\\Delta t}\\); here the displacement is exactly \\(\\varepsilon\\Delta x/c\\) and the random part is a sign that keeps its value until it flips. Compare Part II: there the time coordinate took random Gaussian steps with an imaginary variance, and \\(\\tau\\) was the clock; here the time coordinate takes steps \\(\\pm\\Delta x/c\\) of a definite sign, and the clock is \\(x\\). Part III is Part II with the time kick turned into a coin.

## III.2 The two densities, and what one step does to them

*What is being counted.* The state is \\((t,\\varepsilon)\\) and the clock is \\(x\\), so the probability is a density over \\(t\\) at a given \\(x\\): \\(P\_\\varepsilon(t;x)\\,dt\\) is the probability that the particle, when its clock reads \\(x\\), has time coordinate between \\(t\\) and \\(t+dt\\) and time sense \\(\\varepsilon\\). Two functions of \\((t,x)\\), with \\(\\int\[P\_++P\_-\]\\,dt=1\\) at every \\(x\\). \\(P\_+\\,dt\\) counts forward-time segments, \\(P\_-\\,dt\\) backward-time segments.

*The one-step probability.* Read off the rule: starting from \\((t',\\varepsilon')\\) at clock \\(x\\), one step lands at \\((t,\\varepsilon)\\) at clock \\(x+\\Delta x\\) with probability density

\\\[W(t,\\varepsilon,x+\\Delta x\\,|\\,t',\\varepsilon',x)=\\Big(1-\\frac ac\\Delta x\\Big)\\,\\delta\\Big(t-t'-\\varepsilon'\\frac{\\Delta x}c\\Big)\\,\\delta\_{\\varepsilon\\varepsilon'}+\\frac ac\\Delta x\\;\\delta(t-t')\\,\\delta\_{\\varepsilon,-\\varepsilon'} .\\\]

Two deltas, no spread — the coin of I.1a, now in the time direction.

*The step.* To be at \\((t,\\varepsilon)\\) when the clock reads \\(x+\\Delta x\\): you were at some \\((t',\\varepsilon')\\) when it read \\(x\\), and the step took you there. The line of I.2, with a sum over \\(\\varepsilon'\\) beside the integral over \\(t'\\):

\\\[P\_\\varepsilon(t;x+\\Delta x)=\\sum\_{\\varepsilon'=\\pm1}\\int W(t,\\varepsilon,x+\\Delta x\\,|\\,t',\\varepsilon',x)\\,P\_{\\varepsilon'}(t';x)\\,dt' .\\\]

The deltas do the integrals:

\\\[P\_\\varepsilon(t;x+\\Delta x)=\\Big(1-\\frac ac\\Delta x\\Big)P\_\\varepsilon\\Big(t-\\varepsilon\\frac{\\Delta x}c;x\\Big)+\\frac ac\\Delta x\\,P\_{-\\varepsilon}(t;x).\\\]

Taylor to first order in \\(\\Delta x\\) — \\(P\_\\varepsilon(t-\\varepsilon\\Delta x/c;x)=P\_\\varepsilon-\\varepsilon\\frac{\\Delta x}c\\partial_tP\_\\varepsilon\\), no second-order term survives because the step is \\(O(\\Delta x)\\); the flip term is written at \\(t\\) rather than at the advanced \\(t+\\varepsilon\\Delta x/c\\) of III.1’s “advance, then flip”, a difference of \\(O(\\Delta x^2)\\) that the limit does not see — cancel \\(P\_\\varepsilon\\), divide by \\(\\Delta x\\), multiply by \\(c\\):

\\\[\\boxed{c\\,\\partial_xP\_\\varepsilon+\\varepsilon\\,\\partial_tP\_\\varepsilon=a\\,(P\_{-\\varepsilon}-P\_\\varepsilon)}\\\]

This is the Fokker–Planck step of I.2 for this process, first order because the step is \\(O(\\Delta x)\\). *Remark (T):* the telegraph’s master equations are these with \\(x\\leftrightarrow t\\): \\(\\partial_tP\_\\sigma+\\sigma c\\,\\partial_xP\_\\sigma=a(P\_{-\\sigma}-P\_\\sigma)\\).

## III.3 What is conserved: the charge, and pair creation

Two equations for two densities; combine them. Add:

\\\[c\\,\\partial_x(P\_++P\_-)+\\partial_t(P\_+-P\_-)=0\\qquad\\Longleftrightarrow\\qquad\\boxed{\\partial_t\\,q+\\partial_x\\,j=0,\\quad q\\equiv P\_+-P\_-,\\quad j\\equiv c\\,(P\_++P\_-)}\\\]

Read in the clock \\(x\\): \\(\\partial_x(P\_++P\_-)=-\\frac1c\\partial_t(P\_+-P\_-)\\), the total probability over \\(t\\) is carried along \\(x\\) — the process’s own continuity. Read in \\(t\\), as the forward-time observer does: the quantity conserved in time is \\(q=P\_+-P\_-\\), **forward-time segments minus backward-time segments** — particles minus antiparticles — and every crossing carries it rightward at \\(c\\), so the current is \\(c(P\_++P\_-)\\). This is a conserved charge, not a positive count. Subtract:

\\\[\\partial_t(P\_++P\_-)+c\\,\\partial_x(P\_+-P\_-)=2a\\,(P\_--P\_+).\\\]

The *total* number of crossings is not conserved in \\(t\\): it changes at the rate \\(2a\\) at which the zigzag turns — pair creation at troughs, annihilation at crests. So the forward-time observer, who counts crossings at fixed \\(t\\), has a conserved charge and a non-conserved number. That is exactly the structure Klein–Gordon has (II.9: a conserved current whose time component is not positive), produced here by a classical zigzag before any quantum step. *Remark (T):* for the telegraph the conserved quantity in \\(t\\) is the positive total \\(P\_++P\_-\\) and the current is \\(c(P\_+-P\_-)\\); one particle, never more.

*Second-order form.* Eliminate \\(q\\): from the two equations, \\(c^2\\partial_x^2n+2ac\\,\\partial_xn=\\partial_t^2n\\) for \\(n=P\_++P\_-\\) (and the same for \\(q\\)) — a telegraph equation with \\(x\\) as the clock. *Remark (T):* the telegraph gives \\(\\partial_t^2\\rho+2a\\partial_t\\rho=c^2\\partial_x^2\\rho\\), \\(x\\leftrightarrow t\\).

## III.4 Running the clock backward

Run III.1 backward in \\(x\\), as I.3 did for the random walk. A segment with \\(\\varepsilon=+1\\) read with \\(x\\) decreasing has \\(t\\) decreasing: the two time senses trade places, and the advance is unchanged. The flips: the probability that a particle found with sense \\(\\varepsilon\\) at clock \\(x+\\Delta x\\) had sense \\(-\\varepsilon\\) at clock \\(x\\) is (how likely it was to be there) \\(\\times\\) (how likely that flip was), normalized,

\\\[\\frac{(a/c)\\Delta x\\,P\_{-\\varepsilon}(t;x)}{P\_\\varepsilon(t;x+\\Delta x)}=\\frac{\\tilde a\_\\varepsilon}{c}\\Delta x,\\qquad\\tilde a\_\\varepsilon=a\\,\\frac{P\_{-\\varepsilon}}{P\_\\varepsilon}.\\\]

The reversed process is of the same kind, with flip rates set by the ratio of the two densities. The geometric mean of the forward and reversed rates, times the density, \\(\\sqrt{a\\,\\tilde a\_\\varepsilon}\\,P\_\\varepsilon=a\\sqrt{P\_+P\_-}\\), is the same for both senses: the symmetric flux, the counterpart of I.3’s current velocity. As for the telegraph, and unlike I.3: the advance has no random part, so reversing the clock changes nothing about the motion — **no osmotic velocity appears** — and the information that reversal carries sits in \\(P\_-/P\_+\\), as it sat in \\(\\nabla\\ln\\rho\\) in Part I.

## III.5 The long-clock limit, and the two kinds of mixing

*Kac’s limit, in \\(x\\).* For solutions varying on an \\(x\\)-scale long compared with \\(c/a\\), the second \\(x\\)-derivative in III.3’s second-order form is negligible against the first, and

\\\[\\partial_xn=\\frac1{2ac}\\,\\partial_t^2n :\\\]

a heat equation in which the clock is \\(x\\) and the spreading is in \\(t\\). The time coordinate diffuses along the path, with coefficient \\(1/2a\\) per unit of \\(x/c\\). This is Part II’s walk restricted to its time direction, with a *positive* variance (compare II.1, where the variance had to be \\(-1\\)). *Remark (T):* the telegraph’s Kac limit is \\(\\partial_t\\rho=\\frac{c^2}{2a}\\partial_x^2\\rho\\), Part I’s heat equation, \\(D=c^2/2a\\).

*In Fourier space, as in I.1a.* Transform in \\(x\\) (\\(\\partial_x\\to ik\\)) and write the pair as a column. Take out the common loss first: both densities lose particles at rate \\(a/c\\) per unit \\(x\\), so put \\(P\_\\varepsilon=e^{-ax/c}Q\_\\varepsilon\\) (the steps, each an equality of the whole system, are written out in III.5a); then III.2 reads \\(c\\,\\partial_xQ\_\\varepsilon+\\varepsilon\\,\\partial_tQ\_\\varepsilon=a\\,Q\_{-\\varepsilon}\\), and solving for \\(\\partial_t\\),

\\\[\\partial_t\\begin{pmatrix}Q\_+\\\\Q\_-\\end{pmatrix}=-c\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}\\partial_x\\begin{pmatrix}Q\_+\\\\Q\_-\\end{pmatrix}+a\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}\\begin{pmatrix}Q\_+\\\\Q\_-\\end{pmatrix}.\\\]

(The factor \\(e^{-ax/c}\\) is an exponential of the *clock*, as \\(e^{-at}\\) will be for T; it is the same move in both cases.) In Fourier space \\(\\partial_t\\hat Q=G\\hat Q\\) with

\\\[G=\\begin{pmatrix}-ikc&a\\\\-a&ikc\\end{pmatrix},\\qquad G^2=-(k^2c^2+a^2)\\,1,\\qquad\\text{eigenvalues}\\quad\\pm i\\sqrt{a^2+k^2c^2}.\\\]

**Purely imaginary, for every \\(k\\): the pair rotates, it never decays.** For \\(kc\\ll a\\), \\(\\sqrt{a^2+k^2c^2}=a+\\frac{c^2k^2}{2a}+\\ldots\\): a fast rotation at the rate \\(a\\) times a slow one at \\(\\frac{c^2}{2a}k^2\\) — the dispersion of Schrödinger’s equation with \\(D=c^2/2a\\), from a real rate. *Remark (T):* for the telegraph, \\(G_T=\\begin{pmatrix}-ikc&a\\\\a&ikc\\end{pmatrix}\\), \\(G_T^2=(a^2-k^2c^2)1\\), eigenvalues \\(\\pm\\sqrt{a^2-k^2c^2}\\): real at long wavelength — one mode of the loss-removed \\(Q\\) grows and one decays, which with the loss \\(e^{-at}\\) put back is the *decay* of the \\(P\\) pair, heat — and only for \\(kc>a\\) does it oscillate. The two mixing matrices differ by one sign, and III.8a says what that sign is: \\(\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}\\) is the quarter-turn \\(R(-\\pi/2)\\) of the plane, \\(\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix}\\) is a half-turn seen from the plane. Rotation versus reflection; oscillation versus decay.

### III.5a The loss taken out, and the quarter-turn, line by line

III.5 took the loss out in one sentence. Here are the steps, each an equality of the whole system, so that nothing is done to a single term in words.

**(A) The system of III.2, solved for the \\(x\\)-derivative, stacked.** The two lines \\(c\\,\\partial_xP\_+=-\\partial_tP\_+-aP\_++aP\_-\\) and \\(c\\,\\partial_xP\_-=+\\partial_tP\_--aP\_-+aP\_+\\) are

\\\[c\\,\\partial_x\\begin{pmatrix}P\_+\\\\P\_-\\end{pmatrix}=\\begin{pmatrix}-\\partial_t-a & a\\\\ a & \\partial_t-a\\end{pmatrix}\\begin{pmatrix}P\_+\\\\P\_-\\end{pmatrix}.\\\]

(Top row times the column: \\((-\\partial_t-a)P\_++aP\_-\\), the first line. Bottom row: \\(aP\_++(\\partial_t-a)P\_-\\), the second.) The matrix splits by what each entry does:

\\\[\\begin{pmatrix}-\\partial_t-a & a\\\\ a & \\partial_t-a\\end{pmatrix}=\\underbrace{\\begin{pmatrix}-\\partial_t&0\\\\0&+\\partial_t\\end{pmatrix}}\_{\\text{advance}}+\\underbrace{\\begin{pmatrix}-a&0\\\\0&-a\\end{pmatrix}}\_{\\text{loss}}+\\underbrace{\\begin{pmatrix}0&a\\\\a&0\\end{pmatrix}}\_{\\text{exchange}} .\\\]

The loss is \\(-a\\) times the identity. Everything below uses only that.

**(B) Add \\(a\\binom{P\_+}{P\_-}\\) to both sides.** On the right it cancels the \\(-a\\) on the diagonal:

\\\[c\\,\\partial_x\\begin{pmatrix}P\_+\\\\P\_-\\end{pmatrix}+a\\begin{pmatrix}P\_+\\\\P\_-\\end{pmatrix}=\\begin{pmatrix}-\\partial_t & a\\\\ a & \\partial_t\\end{pmatrix}\\begin{pmatrix}P\_+\\\\P\_-\\end{pmatrix}.\\\]

**(C) An identity for the left side** (product rule, for any column \\(P\\)):

\\\[e^{-ax/c}\\;c\\,\\partial_x\\Big(e^{ax/c}P\\Big)=e^{-ax/c}\\Big(c\\cdot\\tfrac ac\\,e^{ax/c}P+e^{ax/c}\\,c\\,\\partial_xP\\Big)=aP+c\\,\\partial_xP .\\\]

The exponent \\(a/c\\) is the one value for which the right side is exactly the left side of (B); any other leaves a term \\((a-c\\lambda)P\\) over.

**(D) Put (C) into (B).**

\\\[e^{-ax/c}\\;c\\,\\partial_x\\Big(e^{ax/c}\\begin{pmatrix}P\_+\\\\P\_-\\end{pmatrix}\\Big)=\\begin{pmatrix}-\\partial_t & a\\\\ a & \\partial_t\\end{pmatrix}\\begin{pmatrix}P\_+\\\\P\_-\\end{pmatrix}.\\\]

**(E) Multiply both sides by \\(e^{ax/c}\\).** On the right it passes through the matrix, which contains only \\(\\partial_t\\) and constants while the factor depends only on \\(x\\):

\\\[c\\,\\partial_x\\Big(e^{ax/c}\\begin{pmatrix}P\_+\\\\P\_-\\end{pmatrix}\\Big)=\\begin{pmatrix}-\\partial_t & a\\\\ a & \\partial_t\\end{pmatrix}\\Big(e^{ax/c}\\begin{pmatrix}P\_+\\\\P\_-\\end{pmatrix}\\Big).\\\]

**(F) The same column stands on both sides.** It is \\(\\binom{Q\_+}{Q\_-}=e^{ax/c}\\binom{P\_+}{P\_-}\\), and (E) reads

\\\[c\\,\\partial_x\\begin{pmatrix}Q\_+\\\\Q\_-\\end{pmatrix}=\\begin{pmatrix}-\\partial_t & a\\\\ a & \\partial_t\\end{pmatrix}\\begin{pmatrix}Q\_+\\\\Q\_-\\end{pmatrix},\\qquad\\text{i.e.}\\qquad c\\,\\partial_xQ\_+=-\\partial_tQ\_++aQ\_-,\\quad c\\,\\partial_xQ\_-=+\\partial_tQ\_-+aQ\_+ .\\\]

The loss is gone; advance and exchange remain. This is the form the process evolves in (\\(x\\) is its clock) and the one III.14 solves.

**(G) Solve each line for \\(\\partial_t\\)** — the first by moving terms, the second by moving terms and multiplying the line by \\(-1\\):

\\\[\\partial_tQ\_+=-c\\,\\partial_xQ\_++aQ\_-,\\qquad \\partial_tQ\_-=+c\\,\\partial_xQ\_--aQ\_+ .\\\]

**(H) Stack (G), and (I) split the matrix:**

\\\[\\partial_t\\begin{pmatrix}Q\_+\\\\Q\_-\\end{pmatrix}=\\begin{pmatrix}-c\\,\\partial_x & a\\\\ -a & c\\,\\partial_x\\end{pmatrix}\\begin{pmatrix}Q\_+\\\\Q\_-\\end{pmatrix}=-c\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}\\partial_x\\begin{pmatrix}Q\_+\\\\Q\_-\\end{pmatrix}+a\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}\\begin{pmatrix}Q\_+\\\\Q\_-\\end{pmatrix}.\\\]

This is the form the forward-time observer reads (\\(t\\) is his clock), the one III.5 transforms and III.9 reads the Hamiltonian from. (F) and (H) are the same two equations solved for a different derivative: read sideways and read upward, as in III.0’s figure.

**(J) Where the \\(-1\\) in the lower-left came from.** Going from (F) to (H) multiplied the second line by \\(-1\\), which is left-multiplication by \\(\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}\\). Applied to the exchange matrix of (F):

\\\[\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix}=\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}.\\\]

The advance matrix times the exchange matrix is the turn matrix: in the \\(x\\)-clock the equation holds two half-turns, and reading it in \\(t\\) multiplies by one of them, and two half-turns make a rotation (III.8a). The minus is the backward sense running against \\(t\\), and nothing else. The exchange \\(\\sigma_x\\) is symmetric and the turn \\(T\\) antisymmetric; since \\(\\partial_t\\sum Q^2=Q\\cdot(\\mathcal L+\\mathcal L^{\\mathsf T})Q\\) for any generator \\(\\mathcal L\\), the clock \\(x\\) conserves a linear count and the clock \\(t\\) a length, and the one multiplication above is the whole difference (`symmetric_part_and_coincidence.py` C1–C3; III.15c).

**(K) And its square:**

\\\[\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}=\\begin{pmatrix}-1&0\\\\0&-1\\end{pmatrix},\\\]

the \\(-1\\) in front of \\(a^2\\) in III.6: at rest, \\(\\partial_t^2Q=-a^2Q\\), the pair oscillates at the frequency \\(a\\) instead of decaying at the rate \\(a\\) — mass as a frequency. (The telegraph’s exchange, \\(\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix}\\) in both clocks, squares to \\(+1\\): at rest it decays.)

*What a matrix is here.* One matrix per kind of thing that can happen in a step: advancing (the diagonal one, each sense carried by itself, in opposite directions) and flipping (the off-diagonal one, each sense fed into the other). They cannot be merged because one multiplies a derivative and the other multiplies the densities. In three dimensions the list of events is advance along \\(x\\), along \\(y\\), along \\(z\\), and flip: four matrices (IV.4a).

## III.6 Each component obeys Klein–Gordon, with the right sign

With \\(Q\_\\varepsilon=e^{ax/c}P\_\\varepsilon\\) the master equations are \\((c\\,\\partial_x+\\varepsilon\\,\\partial_t)Q\_\\varepsilon=aQ\_{-\\varepsilon}\\). Apply \\((c\\,\\partial_x-\\varepsilon\\,\\partial_t)\\) to the \\(\\varepsilon\\) equation and use the \\(-\\varepsilon\\) one:

\\\[(c\\,\\partial_x-\\varepsilon\\,\\partial_t)(c\\,\\partial_x+\\varepsilon\\,\\partial_t)Q\_\\varepsilon=a\\,(c\\,\\partial_x-\\varepsilon\\,\\partial_t)Q\_{-\\varepsilon}=a^2Q\_\\varepsilon,\\qquad(c\\,\\partial_x-\\varepsilon\\,\\partial_t)(c\\,\\partial_x+\\varepsilon\\,\\partial_t)=c^2\\partial_x^2-\\partial_t^2,\\\]

\\\[\\boxed{\\partial_t^2Q\_\\varepsilon=c^2\\partial_x^2Q\_\\varepsilon-a^2Q\_\\varepsilon}\\\]

**This is the Klein–Gordon equation of Part II, II.8, with \\(a=mc^2/\\hbar\\) — a real flip rate.** Nothing was continued. The rule is simple: a persistent-direction process with clock \\(u\\) and random-direction variable \\(w\\) always produces \\((\\partial_u^2-c^2\\partial_w^2)Q=+a^2Q\\), the clock entering with the plus sign. Klein–Gordon is \\((\\partial_t^2-c^2\\partial_x^2)Q=-a^2Q\\). So the clock has to be \\(x\\) and the random direction has to be \\(t\\). The randomness must be in the time direction — which is what \\(\\eta^{00}=-1\\) said in Part II, and what the coin in \\(t\\) does here.

*Remark (T).* For the telegraph the same computation gives \\(\\partial_t^2Q=c^2\\partial_x^2Q+a^2Q\\): the wrong sign. To reach Klein–Gordon one must set \\(a^2=-m^2c^4/\\hbar^2\\), i.e. continue the rate to \\(a=-imc^2/\\hbar\\). That continuation was the “one knob” of the earlier version of this Part. It is the price of putting the randomness in the space direction.

## III.7 Kinetic term from the process

From III.1, \\(\\Delta t/\\Delta x=\\varepsilon/c\\), so \\(\\langle(\\Delta t/\\Delta x)^2\\rangle=1/c^2\\): a number, independent of the state. In I.5 this step produced \\(\\frac m2(v^2+u^2)\\), the whole Hamiltonian including the osmotic term. Here it produces a constant, which drops. **The Hamiltonian of this particle contains no term built from the densities alone.** It must be found differently (III.10).

## III.8 Matrix form

III.5 already wrote the pair as a column:

\\\[\\partial_t\\begin{pmatrix}Q\_+\\\\Q\_-\\end{pmatrix}=-c\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}\\partial_x\\begin{pmatrix}Q\_+\\\\Q\_-\\end{pmatrix}+a\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}\\begin{pmatrix}Q\_+\\\\Q\_-\\end{pmatrix}.\\\]

Call \\(\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}\\) the **direction matrix** — its diagonal is the two time senses, which to the forward-time observer are the two directions of motion — and \\(\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}\\) the **turn matrix**. The direction matrix squares to \\(+1\\); the turn matrix squares to \\(-1\\):

\\\[\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}^2=\\begin{pmatrix}-1&0\\\\0&-1\\end{pmatrix}.\\\]

They anticommute:

\\\[\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}+\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}=\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix}+\\begin{pmatrix}0&-1\\\\-1&0\\end{pmatrix}=0.\\\]

That is why III.6 came out as a scalar equation with the right sign: squaring the operator on the right, the cross term cancels and the turn matrix contributes \\(-a^2\\). *Remark (T):* the telegraph’s second matrix is the flip \\(\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix}\\), which squares to \\(+1\\) and contributes \\(+a^2\\).

*The turn matrix from conservation alone.* The turn matrix need not be read off; it can be derived from one property of the pair in the \\(t\\)-reading (in the walk’s own clock the conserved form is the indefinite \\(\\int(Q\_+^2-Q\_-^2)dt\\), III.5a (J)’s remark). Multiply the column equation by \\(Q^{\\mathsf T}=(Q\_+,Q\_-)\\) and integrate over \\(x\\):

\\\[\\frac{d}{dt}\\int(Q\_+^2+Q\_-^2)\\,dx=-2c\\int Q^{\\mathsf T}\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}\\partial_xQ\\,dx+2a\\int Q^{\\mathsf T}\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}Q\\,dx=-c\\int\\partial_x(Q\_+^2-Q\_-^2)\\,dx+0=0 .\\\]

The first integral is a total derivative; the second vanishes because the turn matrix is antisymmetric, \\(Q^{\\mathsf T}BQ=Q\_+Q\_--Q\_-Q\_+=0\\). So the length \\(\\int(Q\_+^2+Q\_-^2)dx\\) is conserved in \\(t\\). Run it backward: for an exchange between two components to conserve that length, its matrix must satisfy \\(Q^{\\mathsf T}MQ=0\\) for every \\(Q\\), i.e. be antisymmetric, and a real antisymmetric \\(2\\times2\\) matrix is \\(\\lambda\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}\\) — one matrix, up to the scale \\(\\lambda=a\\). The quarter-turn is forced by “two components, length conserved”. The telegraph’s flip \\(\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix}\\) is symmetric: it conserves \\(\\int(P\_++P\_-)dx\\), a different length, and that is once more the whole difference between the two processes.

*The vector it turns.* With the loss out, the column is a vector in the plane of the two senses whose length cannot change, and a vector whose length cannot change and that moves a little per step can only turn. For a wave that does not move, \\(\\partial_x=0\\): \\(\\partial_tQ=aBQ\\), so \\(Q(t)=R(at)\\,Q(0)\\) — the pair rotates in its plane at the angular rate \\(a\\). That is III.5’s “the pair rotates, it never decays” as a picture, and III.12’s rest-energy phase is this rotation read as a phase. In three dimensions the vector behind the column is the axis \\(\\hat{\\mathbf s}=\\chi^\\dagger\\boldsymbol\\sigma\\chi\\), and it has its own rotation equation (IV.10): the matrices generate turns of the column, the column is the half-angle object, and the axis turns by the full angle.

### III.8a What the two matrices are: half-turns seen from the plane, and what a column is

The two matrices of III.8 have a geometric meaning, and it is worth one page now because the whole of Part IV is built from it.

*What the flatlander sees of a rotation.* Take the \\(xy\\)-plane and look at the shadows of three-dimensional rotations on it. Nothing is added to the rotations; we only ask what each one looks like from the plane. There are three cases.

A rotation about \\(z\\) by \\(\\alpha\\) casts a rotation of the plane, \\(R(\\alpha)=\\begin{pmatrix}\\cos\\alpha&-\\sin\\alpha\\\\\\sin\\alpha&\\cos\\alpha\\end{pmatrix}\\).

A rotation about an axis *in* the plane by a partial angle \\(\\gamma\\) casts a squeeze: about the \\(x\\)-axis, \\((x,y,z)\\to(x,\\,y\\cos\\gamma-z\\sin\\gamma,\\ldots)\\), shadow \\((x,y)\\to(x,\\,y\\cos\\gamma)\\). Not a rotation, not invertible at \\(\\gamma=\\pi/2\\); nothing in this Part uses them.

A **half-turn** about an axis in the plane is the case that matters. The half-turn about the unit axis \\(\\mathbf n\\) is the matrix \\(2\\mathbf n\\mathbf n^{\\mathsf T}-1\\) (it keeps \\(\\mathbf n\\) and reverses everything perpendicular to it). For \\(\\mathbf n=(\\cos\\beta,\\sin\\beta,0)\\) its \\(xy\\)-block is

\\\[2\\begin{pmatrix}\\cos^2\\beta&\\cos\\beta\\sin\\beta\\\\\\cos\\beta\\sin\\beta&\\sin^2\\beta\\end{pmatrix}-\\begin{pmatrix}1&0\\\\0&1\\end{pmatrix}=\\begin{pmatrix}\\cos2\\beta&\\sin2\\beta\\\\\\sin2\\beta&-\\cos2\\beta\\end{pmatrix}\\equiv m(\\beta).\\\]

Determinant \\(-1\\), square \\(1\\): seen from the plane the half-turn acts as a reflection, fixing the line of \\(\\mathbf n\\) (check: \\(m(\\beta)(\\cos\\beta,\\sin\\beta)^{\\mathsf T}=(\\cos\\beta,\\sin\\beta)^{\\mathsf T}\\), the half-turn keeps its own axis) and reversing the line perpendicular to it — the line where the plane of the half-turn cuts the flatlander’s plane. Two of them span all the others:

\\\[m(\\beta)=\\cos2\\beta\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}+\\sin2\\beta\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix}:\\\]

\\(\\beta=0\\) is the half-turn about the \\(x\\)-axis and gives the direction matrix; \\(\\beta=\\pi/4\\) is the half-turn about the diagonal and gives the flip matrix. **The two matrices of III.8 are what half-turns about two in-plane axes look like from the plane.** The anticommutation found in III.8 is a fact about those two half-turns, computed below.

*The label on a plane-axis half-turn is a vector at twice the angle.* The half-turn about \\(\\mathbf n\\) and about \\(-\\mathbf n\\) are the same half-turn, and indeed \\(m(\\beta+\\pi)=m(\\beta)\\). So the half-turns are labelled not by the oriented axis (angle \\(\\beta\\), period \\(2\\pi\\)) but by the *line* of the axis (period \\(\\pi\\)), and the natural way to write a line as a vector is to double its angle: \\(\\hat{\\mathbf v}=(\\cos2\\beta,\\sin2\\beta)\\). Write \\(M(\\hat{\\mathbf v})\\equiv m(\\beta)\\), and for any vector

\\\[M(\\mathbf v)=v_1\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}+v_2\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix}=\\begin{pmatrix}v_1&v_2\\\\v_2&-v_1\\end{pmatrix},\\\]

so the direction matrix is \\(M(\\hat{\\mathbf x})\\) and the flip matrix is \\(M(\\hat{\\mathbf y})\\). Multiply out:

\\\[M(\\mathbf v)^2=|\\mathbf v|^2\\,1,\\qquad M(\\mathbf v)M(\\mathbf w)+M(\\mathbf w)M(\\mathbf v)=2\\,(\\mathbf v\\cdot\\mathbf w)\\,1 .\\\]

Two plane-axis half-turns anticommute exactly when their labels are perpendicular — their axes \\(45^\\circ\\) apart, as the \\(x\\)-axis and the diagonal are. That is the anticommutation of III.8.

*Two half-turns make a rotation.* Half-turns about axes \\(\\beta_1\\) and \\(\\beta_2\\), one after the other:

\\\[m(\\beta_2)m(\\beta_1)=\\begin{pmatrix}\\cos2\\beta_2\\cos2\\beta_1+\\sin2\\beta_2\\sin2\\beta_1&\\cos2\\beta_2\\sin2\\beta_1-\\sin2\\beta_2\\cos2\\beta_1\\\\ \\sin2\\beta_2\\cos2\\beta_1-\\cos2\\beta_2\\sin2\\beta_1&\\sin2\\beta_2\\sin2\\beta_1+\\cos2\\beta_2\\cos2\\beta_1\\end{pmatrix}=R\\big(2(\\beta_2-\\beta_1)\\big).\\\]

Two half-turns about axes \\(\\gamma\\) apart compose to the rotation by \\(2\\gamma\\) about the normal (true in three dimensions, and the shadow says so). **To rotate the plane by \\(\\alpha\\), you apply two plane-axis half-turns whose axes are \\(\\alpha/2\\) apart.** In the labels, \\(M(\\hat{\\mathbf w})M(\\hat{\\mathbf v})=R(\\alpha)\\) when \\(\\hat{\\mathbf w}\\) is \\(\\hat{\\mathbf v}\\) turned by \\(\\alpha\\).

*Once or twice.* \\(R(\\alpha)\\) can hit two kinds of things. It hits a column once, \\(\\chi\\to R(\\alpha)\\chi\\), and the column turns by \\(\\alpha\\). It hits a half-turn by conjugation, \\(m(\\beta)\\to R(\\alpha)\\,m(\\beta)\\,R(\\alpha)^{-1}\\), which is the half-turn about the rotated axis, \\(m(\\beta+\\alpha)\\) (if \\(m\\) keeps the column \\(\\chi\\), \\(RmR^{-1}\\) keeps \\(R\\chi\\)). The axis turned by \\(\\alpha\\) — but the *label*, at twice the angle of the axis, turned by \\(2\\alpha\\):

\\\[R(\\alpha)\\,M(\\mathbf v)\\,R(\\alpha)^{-1}=M\\big(R(2\\alpha)\\,\\mathbf v\\big).\\\]

The same matrix turns a column by \\(\\alpha\\) and a half-turn’s label by \\(2\\alpha\\). A full turn of the label is a half turn of the column, \\(R(\\pi)=-1\\): the column comes back reversed, because the axis \\(\\mathbf n\\) came back as \\(-\\mathbf n\\), which is the same half-turn.

*So:* a **spinor** is a column — an oriented axis, as the flatlander holds it, the thing one plane-axis half-turn hits. A **vector** is the half-turn about the axis — it has forgotten the orientation, and it is hit from both sides. The same rotations act on both; the spinor turns at half the rate of the vector’s label, and \\(\\chi\\) and \\(-\\chi\\) are the same vector. The spinor *of* a unit vector \\(\\hat{\\mathbf v}\\) is the axis its half-turn keeps, \\((\\cos\\beta,\\sin\\beta)\\) with \\(2\\beta\\) the angle of \\(\\hat{\\mathbf v}\\). Nothing here is quantum: it is what rotations look like from a plane.

*Reading III.8 with this.* The column \\((Q\_+,Q\_-)\\) of III.8 is a spinor in a plane whose two axes are “streaming” and “turning”. The direction matrix is the plane-axis half-turn of the first axis: it keeps the pure forward-time column \\((1,0)\\) and reverses the pure backward-time one. The turn matrix is \\(R(-\\pi/2)=m(0)\\,m(\\pi/4)\\), a quarter-turn of the plane — the product of the two half-turns, the thing that squares to \\(-1\\). In the standard basis of III.9 (lower component multiplied by \\(i\\)) it becomes the half-turn of the second axis, \\(\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix}\\), which keeps the equal mixture \\((1,1)/\\sqrt2\\); the two descriptions differ by a relabelling of the relative phase and nothing else.

### III.8b The complex pairing: one rule, and its three uses in this document

*The picture first.* Everywhere this document turns something real into something complex it makes the same move, and the move is not “add an imaginary unit”. It is: find a quarter-turn that is already there, and call it \\(i\\). Two real numbers that a quarter-turn mixes are one complex number. Said once as a rule, it tells in advance what stays linear and what does not, and it says what the \\(i\\) of III.9 is in the walk’s own words: a turn of the walker’s clock hand that leaves its time sense alone. III.1’s two detectors cannot read that turn; the completed detectors of III.8c can.

*Symbols.* \\(J\\) is any real matrix with \\(J^2=-1\\), a quarter-turn. \\(z\\) is the complex number a pair becomes. Nothing else is new.

**(A) The rule.** A real space carrying a quarter-turn \\(J\\) is a complex space with \\(i=J\\). For a pair \\((x,y)\\) and the quarter-turn of the plane,

\\\[J=\\begin{pmatrix}0&-1\\\\1&0\\end{pmatrix},\\qquad z=x+iy,\\qquad Jz=-y+ix=i\\,z .\\\]

The eigenvectors of \\(J\\) are the basis \\((z,\\bar z)\\):

\\\[U^{-1}=\\frac1{\\sqrt2}\\begin{pmatrix}1&i\\\\1&-i\\end{pmatrix},\\qquad U^{-1}\\begin{pmatrix}x\\\\y\\end{pmatrix}=\\frac1{\\sqrt2}\\begin{pmatrix}z\\\\\\bar z\\end{pmatrix},\\qquad U^{-1}JU=\\begin{pmatrix}i&0\\\\0&-i\\end{pmatrix}.\\\]

A real matrix \\(M\\) on the pair becomes, on \\(z\\), one of two things:

\\\[MJ=JM:\\quad Mz=(\\text{a complex number})\\,z\\qquad\\text{linear; a rotation }e^{\\theta J}\\text{ is the number }e^{i\\theta};\\\]

\\\[MJ=-JM:\\quad Mz=(\\text{a complex number})\\,\\bar z\\qquad\\text{antilinear; a half-turn reverses the quarter-turn, so it conjugates.}\\\]

Every matrix in this document either commutes or anticommutes with the \\(J\\) in question, so there is no third case.

**(B) The canonical layer.** If the pair is canonical the quarter-turn is the bracket:

\\\[\\{x,y\\}=1,\\quad z=\\frac{x+iy}{\\sqrt2}\\qquad\\Longrightarrow\\qquad\\{z,\\bar z\\}=-i ,\\\]

and with the polar form \\(\\psi=\\sqrt\\rho\\,e^{iS/\\hbar}\\) and \\(\\{\\rho,S\\}=1\\) the same computation gives \\(\\{\\psi,\\psi^\*\\}=-i/\\hbar\\), which is I.9. The \\(i\\) of quantum mechanics is a canonical pair read as one complex number.

**(C) The three uses.**

1. *I.9, fields.* \\((\\rho,S)\\to\\psi\\). The quarter-turn is the bracket, and Hamilton’s two equations for the pair become one,

\\\[i\\hbar\\,\\partial_t\\psi=\\frac{\\delta H}{\\delta\\psi^\*}.\\\]

2. *III.9 and V.2, the column: the completed detectors.* The four detectors of III.8c (B), \\((Q,Q')\\) — III.1’s two and the two added — and \\(J\\) the walker’s own turn, which carries each \\(Q\_\\varepsilon\\) onto \\(Q'\_\\varepsilon\\) and leaves the time sense alone. \\(J\\) commutes with \\(A\\) and with \\(T\\), so by (A) both stay linear and the pair obeys the real column’s equation written once:

\\\[\\psi=Q+iQ',\\qquad i\\hbar\\,\\partial_t\\psi=\\big\[cA\\hat p+i\\hbar aT\\big\]\\psi .\\\]

Three things about this use. *It needs no complex number:* V.2 runs the four real counters with \\(J\\) as a matrix and the field as a real angle, and the real pair and the complex column agree on the lattice to \\(10^{-15}\\) after 140 ticks (`two_copy_walk_em.py` 5a). *It is the only use in which the dimension doubles, and the rule says why:* on one real column no real \\(J\\) with \\(J^2=-1\\) commutes with \\(A\\) — what commutes with \\(A\\) is diagonal — and the one such \\(J\\) that is a turn is \\(\\pm T\\), with which \\(A\\) anticommutes, so pairing \\((Q\_+,Q\_-)\\) with \\(i=T\\) would make the direction matrix antilinear — a complex reading of one column in which streaming conjugates the phase, which is not Dirac (`complex_pairing_tool.py` D2, D2′). The \\(i\\) needed the added detectors. (On the rest column of positive energy, \\((1,-i)^{\\mathsf T}/\\sqrt2\\) in III.9’s basis, \\(T\\) acts as \\(-i\\): that is this pairing on the one solution where it does no harm, an eigenvector of \\(T\\).) *And what the added detectors are:* III.8c (B) — detectors for the two combinations of time sense and hand position that III.1 does not read. For a walker at rest they hold what the original two hold a quarter of the Compton clock later, \\(v=Tu\\); in motion that is no longer so, and \\(v=Tu\\) is not the one-sign partner either (III.8c (E)). The walker’s own turn, \\(J\\), and the time flip, \\(T\\), are different turns and they commute. The complex column’s length, \\(Q^2+Q'^2\\), is what a detector registers (V.5); without a field the original and the added detectors close separately (III.15a (F)).

3. *The flatlander’s matrices: Cartesian basis to eigenbasis.* \\((x,y)\\to z\\) with \\(J=B_3=-T\\), the quarter-turn of the plane (III.8a’s matrices under the names of the older spin notes: \\(B_2=m(0)\\) the direction matrix, \\(B_1=m(\\pi/4)\\) the flip, \\(B_3=B_1B_2=R(\\pi/2)\\)). Then

\\\[U^{-1}B_3U=i\\sigma_z,\\qquad U^{-1}B_2U=\\sigma_x,\\qquad U^{-1}B_1U=-\\sigma_y ,\\\]

\\\[B_3z=iz,\\qquad B_2z=\\bar z,\\qquad B_1z=i\\bar z .\\\]

The rotation generator becomes a number; the two in-plane half-turns, which anticommute with it, become antilinear. Read backward, \\(\\sigma_x\\to B_2\\) and \\(\\sigma_y\\to-B_1\\) are real — the in-plane spin matrices are the half-turns of III.8a — and only \\(\\sigma_z\\to-iB_3\\), the normal, needs the \\(i\\): III.8a’s “the half-turn about an axis out of his plane”. In 2+1 these three are all a flatlander has, and exponentiated with a real angle \\(B_1,B_2\\) are the runs of IV.4d, not tilts of a coin; the tilts need the second pair of counters (IV.7a). This is the map between the Cartesian basis of the plane, where a rotation is the real \\(B_3\\), and the eigenbasis of the rotation generator, where spin along the normal is diagonal. (`B_to_A_dictionary.py` 2a–2d.)

**(D) Three quarter-turns in three dimensions, one job each** (IV.4b’s real matrices, `B_to_A_dictionary.py` 4a–4d).

\\\[A_xA_yA_z:\\ \\text{commutes with }A_x,A_y,A_z;\\ \\text{anticommutes with }T .\\\]

The \\(i\\) of the massless column: the real four-counter column is one complex two-counter column, IV.3’s, exactly. For that reading the mass is antilinear, so a real column with mass is not a complex two-column — V.2’s “real is neutral” as algebra.

\\\[A_xA_y:\\ \\text{commutes with }A_z\\text{ and }T .\\\]

The \\(i\\) of the massive column moving along \\(z\\): use 2’s complex column read inside the real walker, with \\(i\\) the turn about the axis of motion.

\\\[\\text{no real matrix but a multiple of }1\\text{ commutes with all of }A_x,A_y,A_z,T\\ \\text{— in particular no real }J\\text{ with }J^2=-1 .\\\]

So the charge \\(i\\) of use 2 cannot be a matrix on these four counters. It needs the four added detectors of IV.5a, and there it is the real \\(8\\times8\\) matrix \\(J\\): the walker’s own turn.

*What went in.* Nothing new: (A) is linear algebra, (B) is I.9’s computation, (C)–(D) are matrices and checks already on the page. What it adds is the rule that decides, before any computation, which operators survive a complex reading as linear ones, and the reading of the one doubling the document makes.

(`scripts/complex_pairing_tool.py` (A)–(C); `scripts/B_to_A_dictionary.py` (C)3, (D); `scripts/two_copy_walk_em.py`, `scripts/copy_is_return_leg.py` (C)2.)

### III.8c The measurement, and the completed detectors

*The picture first.* A walker leaves the source and is absorbed at a detector. Two things are read off the walk here, and they are kept apart.

- **The measurement.** Between the source and the detector there are routes: routes that go out from the source, and routes that come back from the detector to the source. Take every outward route and every return route that agree on the endpoint — the same detector, the same sense. Each pair contributes the product of its two signed weights. The sum over all agreeing pairs is the weight of the event, \\(q_X^2\\) (III.15a).
- **The completed detectors.** III.1’s two detectors read the time sense. The walker is also a clock, and a flip moves its hand a quarter of a turn; the two detectors read the hand only through the time sense. Detectors for the two missing combinations — forward with the hand a quarter on, backward with the hand at the start — complete the set to four. The added two are a second column of the walk. They are what lets the walker be at rest, and what gives the phase, and with it the fluid of III.10.

The return routes of the measurement are the walk run backward from the detector: the one process in the document that runs backward. The added detectors are not a way back; they read the same walker. Both are read on time slices, with the source and the detector each sitting on one (III.15b).

*Symbols.* The walk read in \\(t\\) (III.5a (H)), with the loss taken out, on a column of two real counters \\(q=(q\_+,q\_-)^{\\mathsf T}\\):

\\\[\\partial_tq=\\mathcal L\\,q,\\qquad\\mathcal L=-c\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}\\partial_x+a\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}=\\begin{pmatrix}-c\\,\\partial_x&a\\\\-a&c\\,\\partial_x\\end{pmatrix}.\\\]

Three \\(2\\times2\\) matrices, kept apart by name: the direction matrix \\(A_z=\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}\\) (the forward counter moves right, the backward one left), the turn \\(T=\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}\\) (one flip, with its sign), and the sense swap \\(S=\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix}\\) (exchange the two counters; not in the equation). \\(e\_+=(1,0)^{\\mathsf T}\\), \\(e\_-=(0,1)^{\\mathsf T}\\). \\(U(t_2,t_1)\\) is the evolution from \\(t_1\\) to \\(t_2\\); it is orthogonal, \\(U^{\\mathsf T}U=1\\). One tick of length \\(\\Delta t\\) is a shift, then the coin:

\\\[U\_{\\rm tick}=\\begin{pmatrix}\\cos a\\Delta t&\\sin a\\Delta t\\\\-\\sin a\\Delta t&\\cos a\\Delta t\\end{pmatrix}\\begin{pmatrix}L\_-&0\\\\0&L\_+\\end{pmatrix},\\qquad L\_\\pm g(x)=g(x\\pm c\\Delta t).\\\]

**(A) The measurement: every outward route with every return route that agrees on the endpoint.** Let \\(f\\) be the source column, \\(q(t)=U(t,0)f\\) the outward column, \\(U\\equiv U(t_f,0)\\), and \\(e_b\\) the column concentrated at one detector endpoint \\(b\\) (a site and a sense). The returns from \\(b\\) are the walk run backward from \\(e_b\\):

\\\[r_b(t)=U(t_f,t)^{\\mathsf T}e_b ,\\qquad r_b(0)=U^{\\mathsf T}e_b .\\tag{1}\\\]

1. the coin undone, then the shift undone, tick by tick (III.15a (B)); its weight on a route is the product of the same signed coin entries, read in the other order.

\\\[\\alpha_b\\equiv e_b^{\\mathsf T}Uf=f^{\\mathsf T}U^{\\mathsf T}e_b=f^{\\mathsf T}r_b(0),\\qquad W_b=\\sum\_{\\text{out to }b}\\ \\sum\_{\\text{back from }b}w\_{\\rm out}\\,w\_{\\rm back}=\\alpha_b\\cdot\\alpha_b=\\alpha_b^2 .\\tag{2}\\\]

2. the outward sum arriving at \\(b\\) equals the return sum arriving at the source, because \\(U^{\\mathsf T}=U^{-1}\\); every outward route that reaches \\(b\\) pairs with every return route that leaves \\(b\\), so the double sum is the product of the two sums. Single pairs can be negative; the sum per endpoint is a square. \\(\\sum_bW_b=\\|Uf\\|^2=\\|f\\|^2\\). This is III.15a (C), with the event named as what it is: *an outward route and a return route that agree on the endpoint*.

\\\[r_b(t)^{\\mathsf T}q(t)=e_b^{\\mathsf T}U(t_f,t)\\,U(t,0)f=\\alpha_b\\ \\ \\text{at every }t,\\qquad\\sum_b\\alpha_b\\,r_b(t)=q(t).\\tag{3}\\\]

3. the overlap of the two legs is the same on every time slice of the event; all the returns, weighted by their outward amplitudes, rebuild the outward column. (`second_copy_is_the_hilbert_transform.py` (e1)–(e3): \\(10^{-16}\\).)

**(B) The completed detectors.** With the loss out the walker at rest is a clock, \\(u(t)=(\\cos at,\\,-\\sin at)^{\\mathsf T}\\) ((11) below), and a flip is a quarter of a turn of its hand (\\(T^2=-1\\): two flips, half a turn). So III.1’s two detectors are not only “forward” and “backward”: the forward one reads the hand at the start, the backward one reads it a quarter on. The time sense and the position of the hand are one label. A set of detectors that reads both has four:

|  | hand at the start | hand a quarter on |
| --- | --- | --- |
| forward | \\(+u\_+\\) | \\(+v\_+\\) |
| backward | \\(-v\_-\\) | \\(+u\_-\\) |

4. \\(u\_\\pm\\) are III.1’s two; \\(v\_\\pm\\) are the two added. \\(v\_\\varepsilon\\) has the time sense of \\(u\_\\varepsilon\\) and the hand a quarter further on (for \\(u\_-\\): from a quarter on to half a turn, which is the start counted with a minus — the sign on \\(v\_-\\)). The rules do not ask where the hand is: advance by time sense; a flip reverses the time sense and moves the hand a quarter on. So the added two obey the lines of the original two, \\(\\partial_tv=\\mathcal Lv\\). The walker at rest in all four:

\\\[u(t)=\\begin{pmatrix}\\cos at\\\\-\\sin at\\end{pmatrix},\\qquad v(t)=T\\,u(t)=\\begin{pmatrix}-\\sin at\\\\-\\cos at\\end{pmatrix}:\\\]

\\\[\\underbrace{(u\_+,\\,v\_+)}\_{\\text{forward}}=\\underbrace{(-v\_-,\\,u\_-)}\_{\\text{backward}}=(\\cos at,\\,-\\sin at).\\tag{4}\\\]

4a. the same amount in the forward and the backward detector at each position of the hand. Its flow, \\(c\\,(u\_+^2+v\_+^2-u\_-^2-v\_-^2)\\), is zero at every moment: the walker is at rest. The original two alone flow at \\(c\\cos2at\\) ((14) below): at rest only on the average over a turn.

\\\[v(0)=-\\,e\_-\\,g=\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}\\begin{pmatrix}1\\\\0\\end{pmatrix}g=T\\,u(0),\\qquad u(0)=e\_+\\,g .\\tag{5}\\\]

5. the launch at rest: the forward detector with the hand at the start holds the density \\(g(x)\\), and the backward detector with the hand at the start holds the same, \\(-v\_-=g\\). So the added two start as the original two turned by one flip.

\\\[v(t)=U(t,0)\\,T\\,u(0),\\qquad\\partial_tv=\\mathcal Lv .\\tag{6}\\\]

6. the added two are a solution of the same lines, with their own starting data. The free lines never move an amount between the original two and the added two, so the lines do not force the added two: they are taken (input 11), for what (4), (C) and V.2 need of them — rest, the phase as a field, and something for a field to turn. (Read in the walk’s own clock \\(x\\), with data on the line \\(x=0\\), the same four have \\(u\_+u\_-+v\_+v\_-\\ge0\\) everywhere, III.15b; on time slices the same start gives both signs. `zigzag_signs_and_twin_checks.py` 8.)

**(C) What the completed detectors give: the fluid.** The original two, \\(u\\), and the added two, \\(v\\), each obeying (6)’s lines. At each point and time sense form the pair of the two hand positions, and the light-line derivatives:

\\\[r\_\\pm=\\begin{pmatrix}u\_\\pm\\\\v\_\\pm\\end{pmatrix},\\qquad D\_\\pm=\\partial_t\\pm c\\,\\partial_x,\\qquad D\_+r\_+=a\\,r\_-,\\quad D\_-r\_-=-a\\,r\_+ .\\tag{7}\\\]

7. the two lines of \\(\\mathcal L\\), written per sense: each sense advances along its light line and is fed by the other.

\\\[\\rho\_\\pm=u\_\\pm^2+v\_\\pm^2,\\quad\\theta\_\\pm=\\operatorname{atan2}(v\_\\pm,u\_\\pm),\\qquad K=u\_+u\_-+v\_+v\_-,\\quad W=u\_+v\_--v\_+u\_-:\\qquad D\_\\pm\\rho\_\\pm=\\pm2aK,\\qquad\\rho\_\\pm D\_\\pm\\theta\_\\pm=aW .\\tag{8}\\\]

8. differentiate and use (7). The first pair says what each sense gains from the other; the second says how fast each sense’s angle turns.

\\\[S\_+=\\hbar\\theta\_+,\\quad S\_-=\\hbar\\big(\\theta\_-+\\tfrac\\pi2\\big),\\quad\\varphi=\\frac{S\_--S\_+}\\hbar:\\qquad K=\\sqrt{\\rho\_+\\rho\_-}\\,\\sin\\varphi,\\qquad W=-\\sqrt{\\rho\_+\\rho\_-}\\,\\cos\\varphi .\\tag{9}\\\]

9. the constant \\(\\pi/2\\) on the backward sense is III.9’s relabelling to the standard basis (\\(\\psi\_-\\to i\\psi\_-\\)). With \\(a=mc^2/\\hbar\\), (8) becomes

\\\[D\_\\pm\\rho\_\\pm=\\pm\\frac{2mc^2}\\hbar\\sqrt{\\rho\_+\\rho\_-}\\,\\sin\\varphi,\\qquad D\_\\pm S\_\\pm=-mc^2\\sqrt{\\frac{\\rho\_\\mp}{\\rho\_\\pm}}\\,\\cos\\varphi ,\\tag{10}\\\]

10. III.10’s four equations, derived in real variables, with no complex number. They hold for *any* two columns of the walk. What the added detectors buy is that \\(S\_\\pm\\) is a variable at all: with the original two alone the angles are frozen at \\(0\\) or \\(\\pi\\) and the fluid has densities and an imbalance current but no phase gradients. (`zigzag_signs_and_twin_checks.py` 5a, 5b.)

*With the field.* V.2 turns the walker in its own plane — the pair \\((u\_\\varepsilon,v\_\\varepsilon)\\) of each time sense — at each step. In the \\(t\\)-reading (V.2 (E)) the lines (7) become

\\\[D\_+r\_+=a\\,r\_--b\_+\\,Jr\_+,\\qquad D\_-r\_-=-a\\,r\_+-b\_-\\,Jr\_-,\\qquad J=\\begin{pmatrix}0&-1\\\\1&0\\end{pmatrix},\\qquad b\_\\pm=\\frac e\\hbar\\big(\\Phi\\mp cA_x\\big),\\tag{10a}\\\]

10a. \\(Jr\_\\varepsilon=(-v\_\\varepsilon,u\_\\varepsilon)^{\\mathsf T}\\), so the field term is \\(b\_\\varepsilon(v\_\\varepsilon,-u\_\\varepsilon)^{\\mathsf T}\\): V.2 (E)’s two lines per sense. \\(\\Phi\\) turns both senses alike, \\(A_x\\) oppositely.

\\\[D\_\\pm\\rho\_\\pm=\\pm2aK\\quad(\\text{unchanged}),\\qquad\\rho\_\\pm D\_\\pm\\theta\_\\pm=aW-b\_\\pm\\rho\_\\pm .\\tag{10b}\\\]

10b. differentiate as in (8): \\(u(b v)+v(-b u)=0\\) in the density, \\(u(-bu)-v(bv)=-b\\rho\\) in the angle. A turn of the pair changes no density and turns every angle at the rate \\(b\_\\varepsilon\\). With (9),

\\\[D\_\\pm S\_\\pm=-mc^2\\sqrt{\\frac{\\rho\_\\mp}{\\rho\_\\pm}}\\cos\\varphi\\;-\\;e\\big(\\Phi\\mp cA_x\\big):\\tag{10c}\\\]

10c. III.10’s phase equations with minimal coupling, \\(\\hat p\\to\\hat p-eA_x\\) and \\(+e\\Phi\\), i.e. V.2 (F) in fluid variables. The density-weighted rate of the two phases carries

\\\[\\frac{\\rho\_+b\_++\\rho\_-b\_-}{\\rho}=\\frac e\\hbar\\Big\[\\Phi-cA_x\\,\\frac{\\rho\_+-\\rho\_-}{\\rho}\\Big\]=\\frac e\\hbar\\big(\\Phi-v\\,A_x\\big),\\qquad v=\\frac j\\rho=c\\,\\frac{\\rho\_+-\\rho\_-}\\rho,\\tag{10d}\\\]

10d. the \\(vA_x\\) is the imbalance of the two senses, III.11’s \\(c\\cos\\theta\\). Along the packet, \\(dt=\\gamma\\,d\\tau\\):

\\\[\\frac{d\\theta\_{\\rm int}}{d\\tau}=-\\frac e\\hbar\\,\\gamma\\big(\\Phi-vA_x\\big)=-\\frac e\\hbar\\,\\Phi',\\tag{10e}\\\]

10e. \\(\\Phi'=\\gamma(\\Phi-vA_x)\\) is the scalar potential in the walker’s own frame: the field turns the walker’s clock at \\(e\\Phi'/\\hbar\\) per unit of its proper time, with the \\(\\gamma\\) from the change of clock, not fitted. (`interaction_and_plates_checks.py` 1, 2, 2\\('\\), 3: (10b) symbolically; (10d) and \\(\\gamma(\\Phi-VA_x)=\\Phi'\\) on positive-energy columns at \\(V=0,\\,0.6c,\\,0.866c,\\,0.943c\\) to roundoff.) What (10a) assumes is V.2’s rule — the field as the walker’s own turn; VI.2 (G) is where that turn is produced by another walker.

**(D) What the original two already have: the clock, and de Broglie.** At rest the original two turn at the rate \\(a\\):

\\\[u(t)=\\begin{pmatrix}\\cos at\\\\-\\sin at\\end{pmatrix},\\qquad e^{-at}\\,e^{aTt}=\\sum\_{n\\ge0}e^{-at}\\frac{(at)^n}{n!}\\,T^n .\\tag{11}\\\]

11. the weight passes from the forward counter to the backward one and back; the right side is the Poisson sum over the number of flips, \\(T^n\\) for \\(n\\) flips. The angle of the summed rotation is \\(at\\): the flip rate times the time, through a signed sum. (It equals the mean number of flips; it is not an accumulated turn, which would be \\(\\tfrac\\pi2at\\).)

\\\[u(t,x)=\\begin{pmatrix}e^{\\eta/2}&0\\\\0&e^{-\\eta/2}\\end{pmatrix}\\begin{pmatrix}\\cos\\theta\\\\-\\sin\\theta\\end{pmatrix},\\qquad\\theta=a\\gamma\\Big(t-\\frac{Vx}{c^2}\\Big),\\qquad\\tanh\\eta=\\frac Vc .\\tag{12}\\\]

12. the same clock, moving: the senses are rescaled by half the rapidity (III.15c (B)) and the angle is read at the walker’s own time. The angle is read back from the two counters by undoing the rescaling, \\(\\theta=\\operatorname{atan2}\\big(-e^{\\eta/2}u\_-,\\ e^{-\\eta/2}u\_+\\big)\\).

\\\[S=-\\hbar\\theta:\\qquad p=\\partial_xS=\\gamma mV,\\qquad E=-\\partial_tS=\\gamma mc^2,\\qquad\\partial^\\mu S\\,\\partial\_\\mu S=-\\frac{E^2}{c^2}+p^2=-m^2c^2 .\\tag{13}\\\]

13. the momentum is how fast the clock winds along \\(x\\); the energy is how fast it turns in \\(t\\); together they obey II.7’s mass shell. De Broglie’s relations, from the original two detectors.

\\\[\\rho_u=\\cosh\\eta+\\sinh\\eta\\cos2\\theta,\\qquad j_u=c\\,(\\sinh\\eta+\\cosh\\eta\\cos2\\theta).\\tag{14}\\\]

14. what the original two lack: their density and current beat at twice the clock (`second_copy_is_the_hilbert_transform.py` 8). A real column is half a positive-energy mode and half its conjugate, which moves at the same velocity; the beat is their interference. And they have no momentum *density*: the twist \\(\\hbar\\sum_j(u_j\\partial_xv_j-v_j\\partial_xu_j)\\) needs the added detectors (IV.5a (H)).

**(E) The one-sign part, and how wide the source is.** Klein–Gordon (III.6) says \\(\\mathcal L^2=-\\Omega^2\\) with \\(\\Omega=\\sqrt{a^2-c^2\\partial_x^2}\\). So \\(\\mathcal L\\) is a frequency times a quarter turn:

\\\[\\mathcal L=\\Omega\\,J,\\qquad J=\\mathcal L\\,\\Omega^{-1},\\qquad J(k)=\\frac1{\\omega_k}\\begin{pmatrix}-ick&a\\\\-a&ick\\end{pmatrix},\\quad\\omega_k=\\sqrt{a^2+c^2k^2},\\qquad J^2=-1 .\\tag{15}\\\]

15. on a mode of wavenumber \\(k\\); at rest, \\(k=0\\), \\(J=T\\): the quarter turn is the flip itself. Moving, it is the flip read with Klein–Gordon’s frequency, reaching over a Compton length \\(c/a\\).

\\\[u+iJu=2u\_{(+)},\\qquad u-iJu=2u\_{(-)}\\qquad(\\text{every real column }u).\\tag{16}\\\]

16. \\(u\_{(\\pm)}\\) the positive- and negative-energy parts of \\(u\\). The partner \\(v=Ju\\) makes a state of one sign of energy, for any \\(u\\). A walker launched at rest has \\(v(0)=T\\,u(0)\\) ((5)), which is \\(J\\) at \\(k=0\\); for a source of width \\(\\sigma\\) the difference is first order in \\(c/(a\\sigma)\\) and the negative-energy weight second order:

| source width \\(\\sigma\\), in units \\(c/a\\) | \\(\\lVert Tu_0-Ju_0\\rVert/\\lVert u_0\\rVert\\) | negative-energy weight |
| --- | --- | --- |
| point | \\(1.09\\) | \\(0.30\\) |
| \\(5\\) | \\(0.14\\) | \\(0.005\\) |
| \\(15\\) | \\(0.047\\) | \\(0.0006\\) |

So a walker launched at rest, \\(v(0)=Tu(0)\\) of (5), over more than a Compton length is one-sign to order \\((c/a\\sigma)^2\\), and a sharper launch brings its own negative-energy part. Width alone does not do it: the column \\((w\_++w\_-)/\\sqrt2\\) is half negative-energy at every width (`interaction_and_plates_checks.py` 6); it is the launch at rest, spread wide, that is one-sign. That is standard Dirac physics, not a defect: a spinor localized below the Compton wavelength cannot be purely positive-energy (Newton and Wigner, *Rev. Mod. Phys.* 21, 400 (1949)). The document uses one sign of energy in three places: III.11’s slow limit, which works to order \\((v/c)^2\\) and is met by any source wider than \\(c/a\\); Feynman’s rule (IV.10, bottom row), which is not an input of its own but follows from exclusion once the vacuum is taken (VI.4 (C), VI.5); and the split \\(\\Pi\_\\pm\\) on each line of VI.4–VI.6, which is one walker’s own property, the sense in which its clock turns (III.15e (C)). (`ord_entwining_in_the_walk_1plus1.py` 3, 4; `second_copy_is_the_hilbert_transform.py` 9.)

At rest three columns coincide: the added two of the walker launched at rest, \\(v=Tu\\); the one-sign partner \\(Ju\\); and the original two read backward in time with the senses swapped, \\(S\\,u(-t)=-Tu\\). Moving, they part: the one-sign partner \\(J(k)u\\) is a different column from \\(Tu\\), and \\(v\\) is no longer \\(u\\) a quarter clock later. One identity holds at any speed for a launch in the forward sense with \\(v(0)=Tu(0)\\): \\(S\\,u(-t)=-v\\), because \\(S\\mathcal LS=-\\mathcal L\\) (\\(SA_zS=-A_z\\), \\(STS=-T\\)) and \\(S\\,u(0)=-T\\,u(0)\\) (ibid. (c), (c′): \\(3\\times10^{-16}\\) on a moving packet).

*What this settles.* The measurement keeps every outward route with every return route that agrees on the endpoint, and gives \\(q_X^2\\). The completed detectors are four, III.1’s two and the two added; they give the walker at rest, the second column, the phase and III.10’s fluid. Both are read on time slices. What the original two already carry is the clock and its winding, \\(p=\\partial_xS\\); what the added two bring is the phase as a field, the momentum density and a record without the beat. One sign of energy is a property of how wide the source is.

*What went in.* (A) is orthogonality and III.15a’s event, stated as *an outward route and a return route that agree on the endpoint* (input 9). (B) is one thing taken: the detectors read the time sense and the position of the clock hand separately, four instead of two (input 11); the lines do not force it. (C)–(E) are equalities. Not derived: that the weight of an event is its frequency (input 10). (An earlier version built the second column from pairs of routes with the same turns — Ord’s entwining, Ord and Mann, *Phys. Rev. A* 67, 022105 (2003). Here it is a pair of detectors, and no pairing of routes is used outside the measurement.)

(`scripts/second_copy_is_the_hilbert_transform.py`: (e1)–(e3) for (1)–(3); 5, 6, 8, 5d for (11)–(14); 9 for (15)–(16); (c), (c′) for the retraced walk. `scripts/zigzag_signs_and_twin_checks.py`: 5a–5b for (7)–(10), 8 for the cut. `scripts/ord_entwining_in_the_walk_1plus1.py`: 3, 4 for the table. `scripts/detectors_in_master_names.py`: 5 for the table of (B), (4) and the two loops.)

## III.9 The quantum step, the same move as I.9 — and what it no longer needs

Each of the two counters of III.8 becomes a Madelung pair, with \\(\\rho\_\\pm\\) a new variable, not III.1’s \\(P\_\\pm\\):

\\\[\\psi\_\\pm=\\sqrt{\\rho\_\\pm}\\,e^{iS\_\\pm/\\hbar},\\qquad\\{\\rho\_\\pm(x),S\_\\pm(x')\\}=\\delta(x-x'),\\qquad\\text{so}\\qquad\\{\\psi\_\\pm,\\psi\_\\pm^\*\\}=-\\tfrac i\\hbar\\delta(x-x')\\\]

by the Jacobian computation of I.9, the two pairs having zero bracket with each other. Here \\(\\rho\_\\pm\\) is written as a definition; III.15a derives it as the weight of an event, the product of the walk’s two directions between an emitter and an absorber. It is the \\(t\\)-reading’s square, conserved as \\(\\rho\_++\\rho\_-\\) (III.8); III.1’s \\(P\_\\pm\\), conserved as \\(P\_+-P\_-\\) in \\(t\\) (III.3), are the other reading and do not appear from here to III.13.

*Assumption (the one dynamical input, the same as I.7):* the two Madelung fields evolve by Hamilton’s equations, \\(\\dot\\psi\_\\pm=\\{\\psi\_\\pm,H\\}=-\\frac i\\hbar\\,\\delta H/\\delta\\psi\_\\pm^\*\\), for a real energy \\(H\\), and the equation they obey is the equation of III.8 — the classical equation of motion of the two densities, now read for the column \\((\\psi\_+,\\psi\_-)\\). For \\(H=\\int\\Psi^\\dagger\\hat H\\Psi\\) with \\(\\hat H\\) a matrix of operators, Hamilton’s equation is \\(i\\hbar\\dot\\Psi=\\hat H\\Psi\\), and \\(H\\) is real for every \\(\\Psi\\) exactly when \\(\\hat H\\) is Hermitian. Multiply III.8 by \\(i\\hbar\\) to read off \\(\\hat H\\):

\\\[i\\hbar\\,\\partial_t\\begin{pmatrix}\\psi\_+\\\\\\psi\_-\\end{pmatrix}=\\Big\[c\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}\\hat p+i\\hbar a\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}\\Big\]\\begin{pmatrix}\\psi\_+\\\\\\psi\_-\\end{pmatrix},\\qquad\\hat p=-i\\hbar\\partial_x .\\\]

The first operator is Hermitian (\\(-i\\hbar\\partial_x\\) is, and the direction matrix is real symmetric). The second is \\(i\\hbar a\\) times a real *antisymmetric* matrix: \\((i\\hbar a\\,T)^\\dagger=-i\\hbar a\\,T^{\\mathsf T}=i\\hbar a\\,T\\). **It is Hermitian for real \\(a\\).** No condition on the rate; nothing to continue. Name the rate by the mass it produces in III.6, \\(a=mc^2/\\hbar\\), and note \\(i\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}=\\begin{pmatrix}0&i\\\\-i&0\\end{pmatrix}=-\\sigma_y\\):

\\\[\\boxed{a=\\frac{mc^2}\\hbar\\ \\ (\\text{real})}\\qquad\\boxed{i\\hbar\\,\\partial_t\\begin{pmatrix}\\psi\_+\\\\\\psi\_-\\end{pmatrix}=\\Big\[c\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}\\hat p+mc^2\\begin{pmatrix}0&i\\\\-i&0\\end{pmatrix}\\Big\]\\begin{pmatrix}\\psi\_+\\\\\\psi\_-\\end{pmatrix}}\\\]

**The Dirac equation in one dimension**, with the mass matrix \\(-\\sigma_y\\). The matrix form of III.8 with real \\(Q\\) is this same equation on its real solutions. The *lines* are shared; the *solutions* are not: the process evolves data on \\(x=\\) const and its kernel lives outside the cone with the growth \\(e^{ax/c}\\) (III.14), while a normalizable Dirac state is data on \\(t=\\) const inside the cone (III.15b) — the two readings of one equation, and only the second is what a detector reads.

*The standard basis.* Multiply the lower component by a phase, \\(\\psi\_-\\to i\\psi\_-\\) (the unitary \\(U=\\mathrm{diag}(1,i)\\), which leaves the direction matrix alone). Then \\(-\\sigma_y\\to\\sigma_x\\) and the equation takes the form in which the rest of this Part is written:

\\\[i\\hbar\\,\\partial_t\\Psi=\\Big\[c\\,\\sigma_z\\hat p+mc^2\\sigma_x\\Big\]\\Psi,\\qquad\\text{components:}\\quad i\\hbar\\,\\partial_t\\psi\_\\pm=\\mp i\\hbar c\\,\\partial_x\\psi\_\\pm+mc^2\\psi\_\\mp .\\\]

The classical real solutions sit in this basis at relative phase \\(\\varphi=(S\_--S\_+)/\\hbar=\\pi/2\\) modulo \\(\\pi\\) (the sign of \\(Q\_+Q\_-\\) picks the branch; undefined where a reading vanishes). Everything from III.10 on uses this basis; the only thing to remember is that the mass couples the two time senses through a quarter-turn, and the \\(\\sigma_x\\) is that quarter-turn after a relabelling of the phase.

Remarks, all exact.

1. *Where the \\(i\\) comes from.* In Part I the \\(i\\) of Schrödinger’s equation came from the knob \\(D\\to iD\\). Here the classical equation of III.8 already mixes the two time senses by a quarter-turn — a matrix that squares to \\(-1\\), which is all that “\\(i\\)” means — and the Hamiltonian form needs no further input. III.5 showed the consequence before any quantum step: the pair rotates instead of decaying, with Schrödinger’s dispersion \\(\\frac{c^2}{2a}k^2\\) at long wavelength. **The \\(i\\) of quantum mechanics is the quarter-turn by which a particle’s time sense mixes with its opposite.** Two \\(i\\)’s, to keep apart: the one in the *equation* is \\(T\\), a real matrix with \\(T^2=-1\\), and this is Ord’s sentence (quant-ph/0206016, “what we need is the alternating pattern of sign”); the one in the *solutions* — the freedom for the two components to be complex, the phases \\(S\_\\pm\\) — cannot be \\(T\\), because \\(T\\) anticommutes with the direction matrix, and needs the added detectors: III.8b (C)2 and III.8c (B), the completed detectors. In three dimensions this is written out as an \\(8\\times8\\) real system in IV.5a: the \\(i\\) is the walker’s own turn \\(J\\), which carries the original four counters onto the added four, every entry \\(0\\) or \\(\\pm1\\).
2. *The mass is the reversal rate*, \\(a=mc^2/\\hbar\\), the Compton frequency — the rate at which the particle reverses its time sense. A massless particle never reverses.
3. *Klein–Gordon.* III.6 already gave \\(\\partial_t^2\\psi=c^2\\partial_x^2\\psi-\\frac{m^2c^4}{\\hbar^2}\\psi\\) for each component with the real rate: the equation of Part II (II.8). Every solution of Dirac solves Klein–Gordon (III.13 does this in both languages).
4. *The dispersion relation, from III.5.* The eigenvalues \\(\\pm i\\sqrt{a^2+k^2c^2}\\) are \\(\\pm iE/\\hbar\\) with \\(E=\\sqrt{m^2c^4+\\hbar^2k^2c^2}\\): the relativistic energy–momentum relation, from the time-sense particle with its real rate.
5. *The Hamiltonian is a plane-axis half-turn.* In the standard basis, III.8a’s notation, \\(\\hat H=cp\\,M(\\hat{\\mathbf x})+mc^2M(\\hat{\\mathbf y})=M(\\mathbf v)\\) with \\(\\mathbf v=(cp,\\,mc^2)\\): the half-turn of the vector whose two components are the streaming and the mixing. \\(\\hat H^2=|\\mathbf v|^2=c^2p^2+m^2c^4\\) — the energy is the length of that vector — and the positive-energy column is the one \\(M(\\mathbf v)\\) keeps, at *half* the angle of \\((cp,mc^2)\\). For \\(p\\to+\\infty\\) the column is \\((1,0)\\); for \\(p\\to-\\infty\\) it is \\((0,1)\\); for \\(p=0\\) it is \\((1,1)/\\sqrt2\\). A particle at rest is the equal in-phase mixture of the two time senses: the column is steady, the particle going backward in time as often as forward (III.9a). Zitterbewegung, in the standard sense used from here on, is the other thing: the beat between the positive- and negative-energy parts of a column, which a single real column shows (III.8c (14)) and a steady one does not.

*Remark (T).* For the telegraph the second operator is \\(i\\hbar a\\) times the real *symmetric* flip matrix, Hermitian only if \\(i\\hbar a\\) is real: \\(a=-imc^2/\\hbar\\). The telegraph reaches the same Dirac equation, but only by continuing its rate (Gaveau, Jacobson, Kac and Schulman, *Phys. Rev. Lett.* 53, 419 (1984)), after which its \\(P\_\\pm\\) are no longer probabilities. (McKeon and Ord’s 1993 *Ann. Phys.* paper is the \\(3{+}1\\) version of their time-reversal route, not of this one.) For S the rate stays real; its densities solve the same lines, on the other cut (III.15b). The thesis that time reversal removes the continuation is McKeon and Ord’s, *Phys. Rev. Lett.* 69, 3 (1992); their construction has \\(t\\) for its clock, four signed charge densities (forward and return portions of entwined paths), and the sign built in; here the clock is \\(x\\), the two states are positive, and the sign appears when the lines are read in \\(t\\) (Ord, “The Dirac equation in classical statistical mechanics”, quant-ph/0206016, says the same thing about the \\(i\\): “we do not actually need \\(i\\), what we need is the alternating pattern of sign”).

### III.9a Where the Lorentzian sign sits, in S and in T

Part II’s relativity sat in the noise: \\(\\langle\\xi^0\\xi^0\\rangle=-1\\). S has no Gaussian noise; where is its Lorentzian sign? Look at the mixing alone (drop the streaming for a moment), in the \\(t\\)-reading — III.5a (J): in the clock \\(x\\) the mixing matrix is the symmetric \\(\\sigma_x\\), a boost of the pair; it is the clock \\(t\\) that carries \\(T\\). One tick of the \\(t\\)-reading on the column \\((Q\_+,Q\_-)\\) after the loss is removed is, to first order in \\(\\Delta t\\),

\\\[1+a\\,\\Delta t\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}=\\exp\\Big\[a\\,\\Delta t\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}\\Big\]+O(\\Delta t^2),\\qquad\\exp\\big\[\\theta\\,T\\big\]=\\cos\\theta\\,1+\\sin\\theta\\,T,\\\]

a **rotation** of the \\((Q\_+,Q\_-)\\) plane by the angle \\(a\\,\\Delta t\\) — the mixing of the two time senses is a circular rotation, and the full flow conserves a definite form: from III.8,

\\\[\\frac{d}{dt}\\int\\big(Q\_+^2+Q\_-^2\\big)dx=\\int\\big\[2Q\_+(-c\\partial_xQ\_++aQ\_-)+2Q\_-(c\\partial_xQ\_-\\!-aQ\_+)\\big\]dx=-c\\int\\partial_x\\big(Q\_+^2-Q\_-^2\\big)dx=0 .\\\]

So where is the sign? In the *clock*: the quantity the process itself conserves is the total over \\(t\\) at fixed \\(x\\) (III.3), and what is conserved in \\(t\\) is the signed charge \\(P\_+-P\_-\\). The Lorentzian sign of S is the sign of the time sense, \\(\\varepsilon\\), which the forward-time observer reads as particle versus antiparticle.

*Remark (T).* For the telegraph the mixing is \\(\\exp\[a\\Delta t\\,\\sigma_x\]=\\cosh+\\sigma_x\\sinh\\), a **hyperbolic** rotation — a boost of the pair — and the real-rate flow conserves the indefinite \\(\\int(Q\_+^2-Q\_-^2)dx\\). The continuation \\(a\\to-imc^2/\\hbar\\) turns \\(\\cosh+\\sigma_x\\sinh\\) into \\(\\cos-i\\sigma_x\\sin\\): a boost into a rotation, the indefinite form into the definite one. That is what the knob did. S has the rotation from the start, because the coin is in the time direction; its indefinite object is not a quadratic form on the pair but the charge \\(P\_+-P\_-\\), and that one is physical.

*Forward and backward in time, twice.* The turn matrix is the operation \\(\\varepsilon\\to-\\varepsilon\\) composed with a sign. In the standard basis the mass matrix \\(\\sigma_x\\) has eigenvectors \\((1,\\pm1)/\\sqrt2\\) with eigenvalues \\(\\pm1\\), which at rest are the states of energy \\(\\pm mc^2\\); the negative-energy one is the one Feynman reads as the antiparticle. Two layers, then: \\(\\varepsilon=-1\\) is a property of a *segment* of the classical worldline — a backward-time crossing, an antiparticle for the forward-time observer of III.1 — and negative energy is a property of the *wave*, the out-of-phase mixture of both \\(\\varepsilon\\). At rest a positive-energy particle is the equal in-phase mixture of forward and backward segments: the particle goes backward in time as often as forward — the in-phase mixture, not the beat.

## III.10 Dirac in \\((\\rho,S)\\), and the Hamiltonian read off

(Standard basis of III.9, mass matrix \\(\\sigma_x\\).) Multiply the \\(\\psi\_+\\) equation by \\(\\psi\_+^\*\\), using \\(\\psi\_+^\*\\partial\\psi\_+=\\frac12\\partial\\rho\_++\\frac i\\hbar\\rho\_+\\partial S\_+\\) (the identity of I.9) and \\(\\psi\_+^\*\\psi\_-=\\sqrt{\\rho\_+\\rho\_-}\\,e^{i\\varphi}\\) with \\(\\varphi\\equiv(S\_--S\_+)/\\hbar\\):

\\\[\\frac{i\\hbar}2\\partial_t\\rho\_+-\\rho\_+\\partial_tS\_+=-\\frac{i\\hbar c}2\\partial_x\\rho\_++c\\rho\_+\\partial_xS\_++mc^2\\sqrt{\\rho\_+\\rho\_-}(\\cos\\varphi+i\\sin\\varphi).\\\]

Imaginary part, divided by \\(\\hbar/2\\); real part, divided by \\(-\\rho\_+\\); and the same for the left-movers:

\\\[\\boxed{\\partial_t\\rho\_\\pm\\pm c\\,\\partial_x\\rho\_\\pm=\\pm\\frac{2mc^2}\\hbar\\sqrt{\\rho\_+\\rho\_-}\\,\\sin\\varphi,\\qquad\\partial_tS\_\\pm\\pm c\\,\\partial_xS\_\\pm=-mc^2\\sqrt{\\frac{\\rho\_\\mp}{\\rho\_\\pm}}\\cos\\varphi}\\\]

III.15a shows that the first pair are the equations of the event weights of III.1, with \\(\\rho\_\\pm\\) derived rather than defined.

Set them beside III.2. The streaming is unchanged. The exchange has become \\(\\frac{2mc^2}\\hbar\\sqrt{\\rho\_+\\rho\_-}\\sin\\varphi\\): its magnitude \\(2a\\sqrt{\\rho\_+\\rho\_-}\\) is the product of the two readings — for the original two detectors \\(2a\\,Q\_+Q\_-\\) with its sign, for all four \\(2a\\,(Q\_+Q\_-+Q\_+'Q\_-')\\) — — the cross term of IV.4d (B), the exchange of III.15a (G) — not III.4’s geometric mean \\(a\\sqrt{P\_+P\_-}\\) of the densities, which it only resembles; its direction is set by the relative phase. The classical process is the case \\(\\varphi=\\pi/2\\) modulo \\(\\pi\\) (real solutions in the original basis, III.9), where the exchange is \\(\\frac{2mc^2}\\hbar\\sqrt{\\rho\_+\\rho\_-}\\): the two time senses feed each other at the full rate. Letting \\(\\varphi\\) vary is the entire content of the quantum step in fluid language.

Now the Hamiltonian, by requiring \\(\\partial_t\\rho\_\\pm=\\delta H/\\delta S\_\\pm\\) and \\(\\partial_tS\_\\pm=-\\delta H/\\delta\\rho\_\\pm\\) (I.6 run backwards: there the bracket was found from \\(H\\); here \\(H\\) from the bracket):

- \\(\\delta H/\\delta S\_+=-c\\,\\partial_x\\rho\_+\\) requires \\(H\\supset\\int c\\rho\_+\\partial_xS\_+\\), since \\(\\frac\\delta{\\delta S\_+}\\int c\\rho\_+\\partial_xS\_+=-\\partial_x(c\\rho\_+)\\). Likewise \\(H\\supset-\\int c\\rho\_-\\partial_xS\_-\\).
- \\(\\delta H/\\delta S\_+=\\frac{2mc^2}\\hbar\\sqrt{\\rho\_+\\rho\_-}\\sin\\varphi\\) requires \\(H\\supset\\int2mc^2\\sqrt{\\rho\_+\\rho\_-}\\cos\\varphi\\), since \\(\\partial\_{S\_+}\\cos\\frac{S\_--S\_+}\\hbar=\\frac1\\hbar\\sin\\frac{S\_--S\_+}\\hbar\\).

\\\[\\boxed{H=\\int\\Big\[c\\,\\rho\_+\\partial_xS\_+-c\\,\\rho\_-\\partial_xS\_-+2mc^2\\sqrt{\\rho\_+\\rho\_-}\\cos\\frac{S\_--S\_+}\\hbar\\Big\]dx}\\\]

Check the \\(S\\) equations: \\(-\\delta H/\\delta\\rho\_+=-c\\,\\partial_xS\_+-mc^2\\sqrt{\\rho\_-/\\rho\_+}\\cos\\varphi\\) ✓. No further term is needed, as III.7 said. This \\(H\\) is \\(\\mathrm{Re}\\int\\Psi^\\dagger\\hat H\\Psi\\) with \\(\\hat H\\) from III.9 (the imaginary part, \\(-\\frac{\\hbar c}2\\partial_x(\\rho\_+-\\rho\_-)\\), integrates to zero). Reading: the right-movers carry momentum density \\(\\rho\_+\\partial_xS\_+\\) at speed \\(+c\\), the left-movers at \\(-c\\), and the mass couples them with weight \\(\\mathrm{Re}(\\psi\_+^\*\\psi\_-)\\). Compare I.6: kinetic energy linear in \\(\\partial S\\) instead of quadratic, and no osmotic term.

## III.11 The slow limit: Part I comes out

Part I had one density and one phase. Here there are two of each. For a particle moving slowly compared with \\(c\\), Part I should come back, so the question is: which combinations of \\((\\rho\_\\pm,S\_\\pm)\\) are Part I’s \\((\\rho,S)\\), and what do the other two do?

Two of the four combinations are already on the page. The total density \\(\\rho=\\rho\_++\\rho\_-\\) appeared in III.3, and the current there was \\(j=c(\\rho\_+-\\rho\_-)\\), so the velocity in units of \\(c\\) is the *imbalance* \\((\\rho\_+-\\rho\_-)/\\rho\\) — a number between \\(-1\\) and \\(1\\), which we write as a cosine. The relative phase \\(\\varphi=(S\_--S\_+)/\\hbar\\) appeared in III.10 as the thing that steers the exchange. That leaves the mean phase. So:

\\\[\\rho=\\rho\_++\\rho\_-,\\qquad\\bar S=\\tfrac12(S\_++S\_-),\\qquad\\cos\\theta=\\frac{\\rho\_+-\\rho\_-}\\rho,\\qquad\\varphi=\\frac{S\_--S\_+}\\hbar,\\\]

inverted: \\(\\rho\_\\pm=\\rho\\cos^2\\frac\\theta2,\\ \\rho\\sin^2\\frac\\theta2\\) (so that \\(\\rho\_+-\\rho\_-=\\rho\\cos\\theta\\)), \\(S\_\\pm=\\bar S\\mp\\frac\\hbar2\\varphi\\), \\(2\\sqrt{\\rho\_+\\rho\_-}=\\rho\\sin\\theta\\), and \\(\\rho\_+\\partial S\_+-\\rho\_-\\partial S\_-=\\rho\\cos\\theta\\,\\partial\\bar S-\\frac{\\hbar\\rho}2\\partial\\varphi\\). In these variables the Hamiltonian of III.10 is

\\\[H=\\int\\Big\[c\\,\\rho\\cos\\theta\\,\\partial_x\\bar S-\\frac{\\hbar c\\rho}2\\partial_x\\varphi+mc^2\\rho\\sin\\theta\\cos\\varphi\\Big\]dx.\\\]

To use Hamilton’s equations in the new variables we need their brackets. From III.9, \\(\\{\\rho\_+(x),S\_+(x')\\}=\\{\\rho\_-(x),S\_-(x')\\}=\\delta(x-x')\\) and all other brackets zero. Every bracket below carries the same \\(\\delta(x-x')\\), so it is written once on the right:

\\\[\\{\\rho(x),\\bar S(x')\\}=\\{\\rho\_++\\rho\_-,\\tfrac12(S\_++S\_-)\\}=\\tfrac12(1+1)\\,\\delta(x-x'),\\qquad\\{\\varphi(x),P\_\\varphi(x')\\}=\\Big\\{\\frac{S\_--S\_+}\\hbar,\\frac\\hbar2(\\rho\_+-\\rho\_-)\\Big\\}=\\tfrac12\\big(-\\{S\_-,\\rho\_-\\}-\\{S\_+,\\rho\_+\\}\\big)=\\tfrac12(1+1)\\,\\delta(x-x'),\\\]

where the partner of \\(\\varphi\\) turned out to be \\(P\_\\varphi=\\frac\\hbar2(\\rho\_+-\\rho\_-)=\\frac\\hbar2\\rho\\cos\\theta\\) — \\(\\hbar/2\\) times the imbalance density — and the cross brackets \\(\\{\\rho,\\varphi\\},\\{\\rho,P\_\\varphi\\},\\{\\bar S,\\varphi\\},\\{\\bar S,P\_\\varphi\\}\\) vanish by the same count. So \\((\\rho,\\bar S)\\) and \\((\\varphi,P\_\\varphi)\\) are two canonical pairs, each with the bracket \\(\\delta(x-x')\\): Part I’s pair, and a new one made of the relative phase and the imbalance.

### The fast pair

The new pair carries the mass term, so it is the one that moves fast. Write \\(H\\) in the canonical variables: \\(\\cos\\theta=2P\_\\varphi/\\hbar\\rho\\) and \\(\\sin\\theta=\\sqrt{1-(2P\_\\varphi/\\hbar\\rho)^2}\\), so \\(\\partial\\sin\\theta/\\partial P\_\\varphi=-\\frac{2}{\\hbar\\rho}\\frac{\\cos\\theta}{\\sin\\theta}\\). Hamilton’s equations for \\((\\varphi,P\_\\varphi)\\) at fixed \\((\\rho,\\bar S)\\) (the \\(\\partial_x\\varphi\\) term integrated by parts for \\(\\delta/\\delta\\varphi\\)):

\\\[\\dot P\_\\varphi=-\\frac{\\delta H}{\\delta\\varphi}=-\\frac{\\hbar c}2\\partial_x\\rho+mc^2\\rho\\sin\\theta\\sin\\varphi,\\qquad\\dot\\varphi=\\frac{\\delta H}{\\delta P\_\\varphi}=\\frac{2c}\\hbar\\Big\[\\partial_x\\bar S-mc\\cos\\varphi\\cot\\theta\\Big\].\\\]

The \\(mc^2\\) terms drive \\((\\theta,\\varphi)\\) at frequency \\(\\sim mc^2/\\hbar\\), the Compton frequency; \\((\\rho,\\bar S)\\) change on the much slower scale set by \\(v\\) and \\(u\\). Keep the slow solutions, exactly as III.5 separated the slow rotation at \\(\\frac{c^2}{2a}k^2\\) from the fast one at \\(a\\) (the fast, Compton-frequency oscillations are the Zitterbewegung and the antiparticle sector). On it, the right-hand sides vanish to leading order in \\(v/c\\) and \\(u/c\\):

\\\[\\cot\\theta=\\frac{\\partial_x\\bar S}{mc\\cos\\varphi},\\qquad\\sin\\theta\\sin\\varphi=\\frac\\hbar{2mc}\\frac{\\partial_x\\rho}\\rho,\\\]

and with \\(\\theta=\\frac\\pi2+O(1/c)\\), \\(\\varphi=O(1/c)\\), to leading order

\\\[\\boxed{c\\cos\\theta=\\frac{\\partial_x\\bar S}m=v,\\qquad c\\,\\varphi=\\frac\\hbar{2m}\\partial_x\\ln\\rho=u}\\\]

with corrections of relative order \\((v/c)^2,(u/c)^2\\). Nelson’s current velocity is the imbalance between right- and left-movers; Nelson’s osmotic velocity is their relative phase. In I.3 the relation \\(u=D\\partial_x\\ln\\rho\\) came from time-reversing the diffusion; here it comes from the mass term. And III.5 now carries a number: \\(D=c^2/2a=\\hbar/2m\\).

### The energy

Insert into \\(H\\). Write \\(\\cos\\theta=\\delta\\), \\(\\delta=\\partial_x\\bar S/mc\\); \\(\\sin\\theta=1-\\frac{\\delta^2}2+O(\\delta^4)\\), \\(\\cos\\varphi=1-\\frac{\\varphi^2}2+O(\\varphi^4)\\).

| term of \\(H\\) | after substitution |
| --- | --- |
| \\(c\\rho\\cos\\theta\\,\\partial_x\\bar S\\) | \\(\\dfrac{\\rho(\\partial_x\\bar S)^2}m\\) |
| \\(mc^2\\rho\\sin\\theta\\cos\\varphi\\) | \\(mc^2\\rho-\\dfrac{\\rho(\\partial_x\\bar S)^2}{2m}-\\dfrac{\\hbar^2}{8m}\\dfrac{(\\partial_x\\rho)^2}\\rho\\) |
| \\(-\\dfrac{\\hbar c\\rho}2\\partial_x\\varphi=-\\dfrac{\\hbar^2}{4m}\\rho\\,\\partial_x^2\\ln\\rho\\) | \\(+\\dfrac{\\hbar^2}{4m}\\dfrac{(\\partial_x\\rho)^2}\\rho-\\dfrac{\\hbar^2}{4m}\\partial_x^2\\rho\\) |

(third row: \\(\\rho\\,\\partial_x^2\\ln\\rho=\\partial_x^2\\rho-(\\partial_x\\rho)^2/\\rho\\)). Sum: the \\((\\partial_x\\bar S)^2\\) terms give \\(\\frac1m-\\frac1{2m}=\\frac1{2m}\\); the \\((\\partial_x\\rho)^2/\\rho\\) terms give \\(\\frac{\\hbar^2}{4m}-\\frac{\\hbar^2}{8m}=\\frac{\\hbar^2}{8m}\\); the last piece integrates to zero.

\\\[\\boxed{H\\to\\int\\Big\[mc^2\\rho+\\frac{\\rho(\\partial_x\\bar S)^2}{2m}+\\frac{\\hbar^2}{8m}\\frac{(\\partial_x\\rho)^2}\\rho\\Big\]dx}\\\]

Rest energy plus the Hamiltonian of I.6, osmotic term and coefficient included. In Part I the osmotic term came from I.5; here III.7 gave nothing, and the osmotic term came from the relative phase of right- and left-movers: the streaming term supplies \\(m\\rho u^2\\) and the mass term subtracts half, the same bookkeeping as for the kinetic term.

## III.11a The exact relation: the current is the twist, weighted by the mass term

*The picture first.* III.11 found, in the slow limit, that the imbalance of the two senses is the slope of the mean phase, \\(c\\cos\\theta=\\partial_x\\bar S/m\\). Relativity gives the exact form. The current is the slope of the phase, weighted by how strongly the mass term couples the two senses, plus two terms that are derivatives of what one column already has.

\\\[m\\,j=\\rho\\sin\\theta\\cos\\varphi\\;\\partial_x\\bar S\\;+\\;\\frac\\hbar2\\,\\rho\\sin\\varphi\\;\\partial_x\\theta\\;-\\;\\frac\\hbar{2c}\\,\\partial_t\\big(\\rho\\sin\\theta\\sin\\varphi\\big),\\qquad j=c\\,\\rho\\cos\\theta .\\tag{1}\\\]

1. III.11’s variables. \\(\\rho\\sin\\theta\\cos\\varphi=2\\sqrt{\\rho\_+\\rho\_-}\\cos\\varphi\\) is the weight of III.10’s mass term; \\(\\rho\\sin\\theta\\sin\\varphi=2\\sqrt{\\rho\_+\\rho\_-}\\sin\\varphi\\) is its exchange. Gordon’s decomposition (1928) in the walk’s variables.

\\\[\\text{one real column:}\\qquad m\\,j=\\hbar\\,\\big(u\_+\\partial_xu\_--u\_-\\partial_xu\_+\\big)-\\frac\\hbar c\\,\\partial_t\\big(u\_+u\_-\\big).\\tag{2}\\\]

2. the same identity for a single column, in the walk’s own counters. Its current is the twist between the two *senses* plus the rate of change of their product (at rest, \\(j=c\\cos2at\\) comes entirely from the second term: the beat). The twist of the pairs \\((u\_\\varepsilon,v\_\\varepsilon)\\), \\(\\rho\\,\\partial_x\\bar S\\) in (1), appears only with the added detectors.

\\\[\\text{slow:}\\ \\ \\theta\\to\\tfrac\\pi2,\\ \\varphi\\to0:\\quad mc\\cos\\theta=\\partial_x\\bar S;\\qquad\\text{plane wave:}\\ \\ \\cos\\theta=\\frac{pc}E,\\ \\sin\\theta=\\frac{mc^2}E,\\quad m\\,j=\\frac{\\rho\\,p}{\\gamma}.\\tag{3}\\\]

3. the slow limit is III.11’s box; for a plane wave at any speed the current is the momentum density divided by \\(\\gamma\\). (`scripts/current_is_the_weighted_twist.py` 1–5.)

## III.12 The same limit on the wave equation

Each component obeys Klein–Gordon with the real rate (III.6). Put \\(a=mc^2/\\hbar\\) and write \\(e^{-iat}=e^{-imc^2t/\\hbar}\\) for the rest-energy phase. Put \\(\\psi=e^{-iat}\\tilde\\phi\\): \\(\\partial_t^2\\psi=e^{-iat}(-a^2\\tilde\\phi-2ia\\partial_t\\tilde\\phi+\\partial_t^2\\tilde\\phi)\\), so \\(\\tilde\\phi\\) obeys \\(\\partial_t^2\\tilde\\phi-2ia\\,\\partial_t\\tilde\\phi=c^2\\partial_x^2\\tilde\\phi\\) — the telegraph equation of the T remark in III.3, with an imaginary rate. Drop \\(\\partial_t^2\\) under the condition of III.5 (\\(|\\partial_t\\tilde\\phi|\\ll a|\\tilde\\phi|\\)):

\\\[\\partial_t\\tilde\\phi=\\frac{ic^2}{2a}\\partial_x^2\\tilde\\phi=\\frac{i\\hbar}{2m}\\partial_x^2\\tilde\\phi\\qquad\\Longleftrightarrow\\qquad i\\hbar\\,\\partial_t\\tilde\\phi=-\\frac{\\hbar^2}{2m}\\partial_x^2\\tilde\\phi.\\\]

Schrödinger, with \\(D=c^2/2a=\\hbar/2m\\) and the \\(i\\) supplied by the quarter-turn of III.8 (it is the \\(i\\) in \\(e^{-iat}\\)). *Remark (T):* for the telegraph the same two lines with a real rate give the heat equation, \\(\\partial_t\\tilde\\phi=\\frac{c^2}{2a}\\partial_x^2\\tilde\\phi\\), and Schrödinger only after \\(a\\to-ia\\). The Madelung form of this equation is III.11.

## III.13 The square root: Dirac and Klein–Gordon in both languages

Part II gave Klein–Gordon; III.9 gave Dirac; III.9 remark 3 said every component of Dirac solves Klein–Gordon. Here is that relation written out, first in the \\(\\psi\\) language, then in \\((\\rho,S)\\), so that “Dirac is the square root of Klein–Gordon” is a computation in both and not a phrase.

*In \\(\\psi\\).* Apply \\(i\\hbar\\partial_t\\) twice to \\(i\\hbar\\partial_t\\Psi=\\hat H\\Psi\\):

\\\[-\\hbar^2\\partial_t^2\\Psi=\\hat H^2\\Psi,\\qquad\\hat H^2=\\big(c\\,\\sigma_z\\hat p+mc^2\\sigma_x\\big)^2=c^2\\hat p^2+m^2c^4+c\\,mc^2\\,\\hat p\\,(\\sigma_z\\sigma_x+\\sigma_x\\sigma_z)=\\big(c^2\\hat p^2+m^2c^4\\big)1,\\\]

the cross term vanishing because the two plane-axis half-turns anticommute (III.8). So each component obeys \\(\\partial_t^2\\psi=c^2\\partial_x^2\\psi-\\frac{m^2c^4}{\\hbar^2}\\psi\\): II.8. The same thing covariantly: multiply the Dirac equation \\(i\\hbar\\partial_t\\Psi=-i\\hbar c\\sigma_z\\partial_x\\Psi+mc^2\\sigma_x\\Psi\\) on the left by \\(\\sigma_x/c\\) and write \\(\\partial_0=\\partial_t/c\\), \\(\\partial_1=\\partial_x\\):

\\\[i\\hbar\\big(\\gamma^0\\partial_0+\\gamma^1\\partial_1\\big)\\Psi=mc\\,\\Psi,\\qquad\\gamma^0=\\sigma_x=\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix},\\quad\\gamma^1=\\sigma_x\\sigma_z=\\begin{pmatrix}0&-1\\\\1&0\\end{pmatrix}.\\\]

Look at what the two matrices are in III.8a’s plane (standard basis): \\(\\gamma^0\\) is the mass matrix \\(\\sigma_x\\), a plane-axis half-turn (squares to \\(+1\\)); \\(\\gamma^1\\) is the *quarter-turn* \\(R(\\pi/2)\\) of the plane, the product of the two half-turns (squares to \\(-1\\)). (In the original basis of III.8 the roles are exchanged by the phase \\(U\\) of III.9; the algebra is the same.) Time is a half-turn, space a quarter-turn: that is the Lorentzian signature, in rotation language. Then \\((\\gamma^0\\partial_0+\\gamma^1\\partial_1)^2=\\partial_0^2-\\partial_1^2\\) (the cross term cancels, \\(\\gamma^0\\gamma^1=-\\gamma^1\\gamma^0\\)), and

\\\[\\big(i\\hbar\\gamma\\!\\cdot\\!\\partial-mc\\big)\\big(i\\hbar\\gamma\\!\\cdot\\!\\partial+mc\\big)=-\\hbar^2(\\partial_0^2-\\partial_1^2)-m^2c^2=\\hbar^2\\Box-m^2c^2 ,\\\]

Klein–Gordon’s operator, in Part II’s signature. **The factorization is nothing but \\(M(\\mathbf v)^2=\\mathbf v\\cdot\\mathbf v\\) of III.8a, with the Lorentzian dot product**: represent the gradient vector by the plane-axis rotation it generates, and its square is its length. Klein–Gordon is about the length of the four-momentum; Dirac says the four-momentum’s half-turn keeps the column.

*In \\((\\rho,S)\\).* The Dirac fluid of III.10 is first order: it gives \\(\\partial_t\\rho\_\\pm\\) and \\(\\partial_tS\_\\pm\\) in terms of \\(x\\)-derivatives and the other component. Use it twice — differentiate it once more in \\(t\\), substituting it again for the time derivatives that appear — and look at the right-movers’ pair \\((\\rho\_+,S\_+)\\) alone. What comes out (checked symbolically) is Part II’s pair, II.7, for that one component:

\\\[\\partial\_\\mu\\big(\\rho\_+\\partial^\\mu S\_+\\big)=0,\\qquad\\partial^\\mu S\_+\\,\\partial\_\\mu S\_++m^2c^2=\\hbar^2\\frac{\\Box\\sqrt{\\rho\_+}}{\\sqrt{\\rho\_+}},\\\]

and the same for \\((\\rho\_-,S\_-)\\). So each component of the Dirac fluid is a Klein–Gordon fluid of Part II, by itself. What Dirac adds is the other component and the first-order tie between the two: the exchange \\(\\frac{2mc^2}\\hbar\\sqrt{\\rho\_+\\rho\_-}\\sin\\varphi\\) and the phase coupling \\(-mc^2\\sqrt{\\rho\_\\mp/\\rho\_\\pm}\\cos\\varphi\\), both built from the Madelung product \\(\\psi\_+^\*\\psi\_-=\\sqrt{\\rho\_+\\rho\_-}e^{i\\varphi}\\). Part II has one density and one phase; Part III has two of each, and the second pair is the time sense of III.1 — which way along the worldline time runs — carried as a field. Klein–Gordon is the same particle with that information averaged out, and that is why it is second order and has twice as many solutions per component.

## III.14 The kernel: what one zigzag does, exactly

*The picture first.* Start the particle at \\(x=0\\), \\(t=0\\), with its time running forward. Let its clock advance to \\(x\\). What does the forward-time observer find? Three things, and they are worth naming in words before any formula.

1. With some probability the particle never flipped. Then it is exactly on the light line \\(t=x/c\\): one point, no spread. The probability of “no flip up to clock \\(x\\)” is the survival probability of a Poisson event with rate \\(a/c\\) per unit \\(x\\), which is \\(e^{-ax/c}\\).
2. If it flipped an odd number of times it is now running backward in time — the forward-time observer calls it an antiparticle. It cannot be on the light line, because at least one segment went the other way.
3. If it flipped an even number of times (at least two) it is running forward again, but it has lost ground in \\(t\\): it is somewhere strictly inside the wedge \\(|t|\<x/c\\).

Everything inside the wedge is *outside* the light cone of the starting point as a forward-time observer draws it — a zigzag in \\(t\\) with \\(x\\) as clock can only reach \\(|t|\\le x/c\\). So the whole kernel lives on and outside the cone. That is the opposite of the telegraph particle, whose kernel lives on and inside it, and it is the reason the two processes continue into each other across the light line (III.15).

III.5 kept only the long-clock limit of III.2, which throws away the light-line point (item 1) and smooths the rest. Here III.2 is solved exactly, so that all three items are visible.

*Take out the loss, as in III.5.* With \\(P\_\\varepsilon=e^{-ax/c}Q\_\\varepsilon\\) the master equations are \\((c\\,\\partial_x+\\varepsilon\\,\\partial_t)Q\_\\varepsilon=aQ\_{-\\varepsilon}\\), and each \\(Q\_\\varepsilon\\) obeys III.6, \\(c^2\\partial_x^2Q\_\\varepsilon-\\partial_t^2Q\_\\varepsilon=a^2Q\_\\varepsilon\\).

*The data at \\(x=0\\).* The clock is \\(x\\), so what is given is the pair at \\(x=0\\) and its first \\(x\\)-derivative there: \\(Q\_+(t;0)=\\delta(t)\\), \\(Q\_-(t;0)=0\\), and from the master equations at \\(x=0\\),

\\\[c\\,\\partial_xQ\_+(t;0)=-\\partial_t\\delta(t),\\qquad c\\,\\partial_xQ\_-(t;0)=a\\,\\delta(t).\\\]

*Fourier in \\(t\\)* — the variable that is not the clock, \\(\\partial_t\\to-i\\omega\\), \\(\\hat Q(\\omega;x)=\\int Q(t;x)e^{i\\omega t}dt\\). III.6 becomes \\(c^2\\partial_x^2\\hat Q\_\\varepsilon=(a^2-\\omega^2)\\hat Q\_\\varepsilon\\), an ordinary equation in \\(x\\) whose two solutions are \\(\\cosh(\\kappa x/c)\\) and \\(\\sinh(\\kappa x/c)/\\kappa\\) with \\(\\kappa=\\sqrt{a^2-\\omega^2}\\). The data pick the combination:

\\\[\\hat Q\_+(\\omega;x)=\\cosh\\frac{\\kappa x}c+i\\omega\\,\\frac{\\sinh(\\kappa x/c)}{\\kappa},\\qquad\\hat Q\_-(\\omega;x)=a\\,\\frac{\\sinh(\\kappa x/c)}{\\kappa}.\\\]

*Back to \\(t\\).* One transform does everything, because \\(\\cosh(\\kappa x/c)=c\\,\\partial_x\[\\sinh(\\kappa x/c)/\\kappa\]\\) and \\(i\\omega\\leftrightarrow-\\partial_t\\). The one transform is a standard pair:

\\\[\\frac1{2\\pi}\\int\_{-\\infty}^{\\infty}\\frac{\\sinh(\\kappa x/c)}{\\kappa}\\,e^{-i\\omega t}\\,d\\omega=\\frac12\\,I_0(\\zeta)\\,\\Theta(x-c|t|),\\qquad\\zeta=\\frac ac\\sqrt{x^2-c^2t^2},\\\]

\\(I_0\\) the modified Bessel function. Then \\(Q\_-=\\frac a2I_0(\\zeta)\\Theta\\) directly, and \\(Q\_+=(c\\,\\partial_x-\\partial_t)\[\\frac12I_0\\Theta\]\\): the derivatives acting on \\(\\Theta\\) give \\(\\frac12\[\\delta(t-x/c)+\\delta(t+x/c)\]+\\frac12\[\\delta(t-x/c)-\\delta(t+x/c)\]=\\delta(t-x/c)\\), and acting on \\(I_0\\) give \\(\\frac a2\\frac{x+ct}{\\sqrt{x^2-c^2t^2}}I_1(\\zeta)\\). Restoring \\(e^{-ax/c}\\):

\\\[\\boxed{\\begin{aligned}P\_+(t;x)&=e^{-ax/c}\\Big\[\\delta\\big(t-\\tfrac xc\\big)+\\frac a2\\,\\frac{x+ct}{\\sqrt{x^2-c^2t^2}}\\,I_1(\\zeta)\\,\\Theta(x-c|t|)\\Big\]\\\\\[4pt\]P\_-(t;x)&=e^{-ax/c}\\;\\frac a2\\,I_0(\\zeta)\\,\\Theta(x-c|t|)\\end{aligned}}\\\]

(Monte Carlo of III.1 at \\(c=2\\): both lines, including the coefficient \\(a/2\\) — not \\(a/2c\\) — and the mass \\(e^{-ax/c}\\) of the \\(\\delta\\).)

*Reading it against the picture.* The \\(\\delta\\) is item 1, the never-flipped path, on the light line with weight \\(e^{-ax/c}\\), in \\(P\_+\\) only. \\(P\_-\\) is item 2, the antiparticle, starting at order \\(a\\) (one flip at least). The \\(I_1\\) term of \\(P\_+\\) is item 3. The bulk of \\(P\_-\\) depends on \\(x\\) and \\(t\\) only through \\(x^2-c^2t^2\\) — through \\(\\zeta=a\\sqrt{-t_1t_2}\\), the flip rate times the geometric mean of the event’s two null track coordinates (III.15c (A)); the bulk of \\(P\_+\\) carries besides the null prefactor \\((x+ct)/\\sqrt{x^2-c^2t^2}\\), one track coordinate over the other, which is not a function of the interval: it is the \\(+\\) sense’s own direction, III.15c (B)’s \\(k^{\\mp1/2}\\) — times the loss factor — the loss factor is the only part that knows the frame.

*One more thing the equation allows.* Both \\(I_0(\\zeta)\\) and \\(K_0(\\zeta)\\) solve \\(c^2\\partial_x^2Q-\\partial_t^2Q=a^2Q\\) outside the cone. The kernel is \\(I_0\\) because the question was “given the particle at \\(x=0\\), where is it at clock \\(x\\)” — data on the line \\(x=0\\), regular at \\(\\zeta=0\\). \\(K_0\\) blows up at \\(\\zeta=0\\) and decays far away: it answers “a source sits at the origin; what field does it make”, a different question about the same equation (III.16).

## III.15 What a particle is, and what the amplitude sums over

*The picture first.* The forward-time observer of III.1 counts crossings at fixed \\(t\\). Between two flips he sees one thing moving at \\(c\\): a segment of the worldline with no flip in it. At a flip he sees a pair appear or disappear (III.3). So if asked “what is a particle”, he has a definite answer: **a no-flip segment**. That is classical — it is III.1 read by someone who refuses to read \\(t\\) backward — and its kernel is III.14’s \\(\\delta\\)-term: light-like, with survival weight \\(e^{-ax/c}\\).

After the quantum step of III.9 there is a second object, the column \\(\\psi\\) obeying Dirac. Here a flip is not an event at all: it is the mass term, the quarter-turn acting on the column, and it acts at every instant on every component. Asking “how many times did this wave flip” is asking the wrong kind of question — the wave is the sum over *all* numbers of flips at once. Two layers, then: on the path a flip ends a particle; in the wave the flips, summed to all orders, are what give the particle its mass. The two do not contradict each other, and the rest of this section makes each statement a computation.

*The amplitude.* With \\(\\hat H\\) of III.9 and \\(\\hat H^2=E^2\\), \\(E=\\sqrt{p^2c^2+m^2c^4}\\),

\\\[e^{-i\\hat Ht/\\hbar}=\\cos\\frac{Et}\\hbar-i\\,\\frac{\\hat H}{E}\\sin\\frac{Et}\\hbar ,\\\]

and the amplitude to go from \\((0,0)\\) with sense \\(+\\) to \\((t,x)\\) is \\(K(t,x)=\\int\\frac{dp}{2\\pi\\hbar}\\,e^{ipx/\\hbar}e^{-i\\hat H(p)t/\\hbar}\\). Doing the \\(p\\)-integral, with \\(\\zeta'=\\frac ac\\sqrt{c^2t^2-x^2}\\) and \\(a=mc^2/\\hbar\\),

\\\[K\_{++}(t,x)=\\delta(x-ct)-\\frac a{2c}\\,\\frac{ct+x}{\\sqrt{c^2t^2-x^2}}\\,J_1(\\zeta')\\,\\Theta(ct-|x|),\\qquad K\_{-+}(t,x)=-\\frac{ia}{2c}\\,J_0(\\zeta')\\,\\Theta(ct-|x|).\\\]

*Three computations of one thing.* The following were run against each other numerically and agree: (i) the \\(p\\)-integral above; (ii) the lattice sum over *every* zigzag in \\(x\\) with \\(t\\) as clock, each flip weighted \\(-ia\\,\\Delta t\\) in place of the probability \\(a\\,\\Delta t\\) (Feynman’s checkerboard); (iii) the Bessel formulas. So:

1. **The amplitude is the sum over all zigzags, all numbers of flips.** Calling a no-flip segment a particle throws nothing away; the flipped paths are all in \\(K\\).
2. **The same amplitude is a superposition of straight paths.** Written in \\(p\\), \\(K=\\int dp\\,e^{ipx/\\hbar}\\big\[e^{-iEt/\\hbar}\\Pi\_+(p)+e^{+iEt/\\hbar}\\Pi\_-(p)\\big\]\\) with \\(\\Pi\_\\pm=\\frac12(1\\pm\\hat H/E)\\): each term is a plane wave of definite momentum and energy — a particle that never flips, on the mass shell. The sum over zigzags *is* this sum, resummed. This is III.9a’s statement that the mixing is a rotation, carried out to all orders.

*One function on both sides of the light cone.* Set III.14 beside the amplitude:

|  | clock | support | no-flip term | bulk |
| --- | --- | --- | --- | --- |
| III.14 (S, classical) | \\(x\\) | \\(x\\ge c\\lvert t\\rvert\\) | \\(\\delta(t-x/c)\\,e^{-ax/c}\\) | \\(\\frac a2\\,I_0(\\zeta)\\), \\(\\frac a2\\frac{x+ct}{\\sqrt{x^2-c^2t^2}}I_1(\\zeta)\\) |
| Dirac amplitude | \\(t\\) | \\(ct\\ge\\lvert x\\rvert\\) | \\(\\delta(x-ct)\\) | \\(-\\frac{ia}{2c}J_0(\\zeta')\\), \\(-\\frac a{2c}\\frac{ct+x}{\\sqrt{c^2t^2-x^2}}J_1(\\zeta')\\) |

\\(I_0(\\zeta)=J_0(i\\zeta)\\) and \\(\\sqrt{x^2-c^2t^2}=i\\sqrt{c^2t^2-x^2}\\). So the bulk of the classical kernel outside the cone and the bulk of the amplitude inside the cone are **the same analytic function of the invariant interval**, continued across the light line. What changes in crossing is the weight on a flip — \\(\\frac a2\\) per unit \\(t\\) outside, \\(-\\frac{ia}{2c}\\) per unit \\(x\\) inside — and which coordinate is the clock. That is the \\(i\\) of III.9a as a formula: the real-rate classical process lives outside the cone, the amplitude inside, and the quarter-turn takes one into the other.

*The same statement from the lines, without the word “continued”.* \\(c=1\\); \\(\\sigma\\equiv\\frac{a^2}4(x^2-t^2)\\), so \\(\\zeta=2\\sqrt\\sigma\\) outside and \\(\\zeta'=2\\sqrt{-\\sigma}\\) inside; \\(f\\), \\(h\\) functions of \\(\\sigma\\); \\(c_0\\) a number. Every line is checked (`one_function_from_the_lines.py`).

\\\[(\\partial_x-\\varepsilon\\partial_t)(\\partial_x+\\varepsilon\\partial_t)Q\_\\varepsilon=a(\\partial_x-\\varepsilon\\partial_t)Q\_{-\\varepsilon}=a^2Q\_\\varepsilon\\quad\\Longrightarrow\\quad(\\partial_x^2-\\partial_t^2)Q\_\\varepsilon=a^2Q\_\\varepsilon\\tag{1}\\\]

1. III.5a (F) applied twice; the same whether the lines are later solved for \\(\\partial_x\\) or for \\(\\partial_t\\).

\\\[Q\_-=f(\\sigma):\\qquad(\\partial_x^2-\\partial_t^2)f=f''\\tfrac{a^4}4(x^2-t^2)+f'\\tfrac{a^2}2+f'\\tfrac{a^2}2=a^2(\\sigma f''+f')\\tag{2}\\\]

2. chain rule, \\(\\partial_xf=f'\\frac{a^2x}2\\), \\(\\partial_tf=-f'\\frac{a^2t}2\\).

\\\[\\sigma f''+f'=f\\tag{3}\\\]

3. 2. into (1). It contains \\(\\sigma\\) and nothing else — not \\(x\\), not \\(t\\), not the sign of \\(\\sigma\\).

\\\[f=\\sum_nc_n\\sigma^n:\\quad\\sum_n(n+1)^2c\_{n+1}\\sigma^n=\\sum_nc_n\\sigma^n\\quad\\Longrightarrow\\quad f=c_0\\,g_0(\\sigma),\\qquad g_0(\\sigma)=\\sum_n\\frac{\\sigma^n}{(n!)^2}\\tag{4}\\\]

4. converges for every \\(\\sigma\\) of either sign; the second solution of (3) is singular at \\(\\sigma=0\\) (it is III.16’s \\(K_0\\)) and a kernel starting from a point on the cone cannot contain it.

\\\[(\\partial_x+\\partial_t)Q\_+=aQ\_-=ac_0g_0(\\sigma),\\qquad Q\_+^{\\rm bulk}=(x+t)\\,h(\\sigma):\\qquad(\\partial_x+\\partial_t)\\big\[(x+t)h\\big\]=2h+(x+t)\\tfrac{a^2}2(x-t)h'=2h+2\\sigma h'\\tag{5}\\\]

5. the \\(+\\) line with (4); the null factor \\(x+t\\) turns it into an equation in \\(\\sigma\\) again.

\\\[h+\\sigma h'=\\tfrac{ac_0}2g_0(\\sigma)\\quad\\Longrightarrow\\quad h=\\tfrac{ac_0}2\\,g_1(\\sigma),\\qquad g_1(\\sigma)=\\sum_n\\frac{\\sigma^n}{n!\\,(n+1)!}\\tag{6}\\\]

6. by the same series, \\((n+1)d_n=1/(n!)^2\\); every \\(\\sigma\\).

\\\[c_0=Q\_-\\big|\_{\\sigma=0}=\\begin{cases}+\\frac a2&\\text{walk's clock: the flip into }-\\text{ weighs }+a\\ \\ (\\text{exchange }\\sigma_x,\\ \\text{III.5a (A)})\\\\\[2pt\]-\\frac a2&t\\text{-reading: the flip into }-\\text{ weighs }-a\\ \\ (\\text{turn }T,\\ \\text{III.5a (J)})\\end{cases}\\tag{7}\\\]

7. the weight the no-flip \\(\\delta\\) hands to the first flip at the cone: the one number that differs.

\\\[Q\_-=c_0\\,g_0(\\sigma),\\qquad Q\_+^{\\rm bulk}=\\tfrac{ac_0}2\\,(x+t)\\,g_1(\\sigma),\\qquad g_0(\\sigma)=I_0(2\\sqrt\\sigma)=J_0(2\\sqrt{-\\sigma})\\tag{8}\\\]

8. with \\(c_0=+\\frac a2\\) this is III.14’s box (\\(\\frac{a}{2}\\frac{x+t}{\\sqrt{x^2-t^2}}I_1(\\zeta)=\\frac{a^2}4(x+t)g_1\\)); with \\(c_0=-\\frac a2\\) it is the real column of the amplitude above (the \\(-i\\) is this \\(-\\) after III.9’s relabelling of the lower component). Same \\(g_0\\), same \\(g_1\\), same \\(x+t\\) on both sides of the cone. The two columns of the table are one regular solution of (3) and (6) with one number different; the equations never knew which side they were on.

### III.15a The transaction: the square as the product of the two directions

*The picture first.* III.15 called a no-flip segment a particle and showed that the amplitude is the sum over every zigzag, all numbers of flips. This section asks a different question of the same paths: not “how does the reading get from the source to a point”, but “what is an event”. The answer it finds is that an event has two endpoints. The walker is emitted at the source and absorbed at a point; between the two endpoints it goes out along every path and comes back along every path; the weight of the event is the sum over the outbound paths times the sum over the return paths. Each sum is the reading \\(Q\\). Their product is \\(Q^2\\). That is where the square comes from, and it needs no complex number and no added detectors: one real column, read once in each direction between two endpoints.

*Symbols.* \\(w_i\\) is the signed weight of one path \\(i\\) of the walk read in \\(t\\), the product of its step factors. Nothing else is new.

**(A) The walk read in \\(t\\), step by step.** III.5a (H) with the loss out: in one tick \\(dt\\), the forward counter moves right by \\(c\\,dt\\) and the backward counter moves left, and then the coin acts,

\\\[\\begin{pmatrix}Q\_+\\\\Q\_-\\end{pmatrix}\\to\\begin{pmatrix}\\cos a\\,dt&\\sin a\\,dt\\\\-\\sin a\\,dt&\\cos a\\,dt\\end{pmatrix}\\begin{pmatrix}Q\_+\\\\Q\_-\\end{pmatrix},\\\]

which is \\(\\cos(a\\,dt)\\,1+\\sin(a\\,dt)\\,T\\). A path is a sequence of ticks with, at each tick, “keep the sense” or “flip”. Its weight \\(w_i\\) is \\(\\cos a\\,dt\\) per keep, \\(+\\sin a\\,dt\\) per flip into the forward sense, \\(-\\sin a\\,dt\\) per flip into the backward sense. The reading at a point is the sum of the weights of the paths that end there with that sense: \\(Q\_\\varepsilon(t,x)=\\sum_iw_i\\). (This is III.15’s “sum over every zigzag”; the \\(-i\\) of Feynman’s rule is the sign of \\(T\\), and the sum is real. `scripts/checkerboard_is_the_walk.py`: the real signed path sum, column kept real at every tick, lands on III.15’s \\(-\\frac a{2c}J_0(\\zeta')\\) to \\(0.7\\%\\) at \\(1600\\) ticks, site for site. The \\(0.7\\%\\) is a one-site offset of the backward counter, which is written after the shift: read one site over, it agrees to \\(8\\times10^{-6}\\), and the forward counter agrees with its \\(J_1\\) form to \\(5\\times10^{-6}\\).)

*Follow the sign at each turn.* Start a sign at \\(+1\\). A keep leaves it alone; a flip from \\(+\\) to \\(-\\) reverses it; a flip from \\(-\\) to \\(+\\) leaves it alone — the signs of the three entries of the coin. With \\(\\theta=a\\,dt\\), a route of \\(N\\) ticks with \\(m\\) flips weighs

\\\[w=(-1)^{n\_{+\\to-}}\\,(\\sin\\theta)^m(\\cos\\theta)^{N-m} .\\\]

For example, release the walker in \\(+\\) for three ticks, distance in units of \\(c\\,dt\\). Two routes end at \\(x=1\\) in sense \\(-\\): *keep, flip, keep*, through positions \\((1,2,1)\\), weight \\(-\\cos^2\\theta\\sin\\theta\\); and *flip, flip, flip*, through \\((1,0,1)\\), weight \\(+\\sin^3\\theta\\). At \\(\\theta=\\pi/4\\) the detector is reachable by both and its reading is exactly zero. Go out by the first and come home by the reverse of the second: the pair weighs \\(-\\cos^2\\theta\\sin^4\\theta\\), and so does the reverse pairing. The two routes have different numbers of turns; pairs like these are where the interference sits, and they are why the event weight is not a sum over routes paired with themselves. (`zigzag_signs_and_twin_checks.py` 4, 9.)

**(B) The backward walk is the same walk with the senses swapped.** One tick is a shift of each sense followed by the coin \\(R=\\cos(a\\,dt)\\,1+\\sin(a\\,dt)\\,T\\). Both are orthogonal, so the evolution \\(U\\) over any number of ticks is orthogonal, \\(U^{\\mathsf T}U=1\\), and its inverse is its transpose. The inverse tick is \\(R^{\\mathsf T}=\\cos(a\\,dt)\\,1-\\sin(a\\,dt)\\,T\\) followed by the reversed shift — the order reverses with the inverse, \\((RS)^{-1}=S^{-1}R^{\\mathsf T}\\) — and \\(R^{\\mathsf T}\\) is the same coin with the flip signs reversed. Reversing the flip signs is swapping the senses, \\(A_xR^{\\mathsf T}A_x=R\\). So a return path from a point to the source is a path of the same walk, read with the senses swapped, and its weight is the product of the elements of \\(R^{\\mathsf T}\\) it selects. (`scripts/transaction_lines.py`, checks 1, 4; `scripts/copy_is_return_leg.py`, check 1: \\(A_xQ(-t,x)\\) solves III.1.)

**(C) The event, line by line.** Source \\(A=(0,0,+)\\), point \\(X=(t,x,\\varepsilon)\\).

1. *Outbound sum.* The sum over paths from \\(A\\) to \\(X\\) of the product of the selected matrix elements is the matrix element of the product: \\(\\displaystyle\\sum_iw_i=U\_{XA}=Q\_\\varepsilon(t,x)\\). This is III.15’s “the amplitude is the sum over all zigzags” (check 2).
2. *Return sum.* The sum over paths from \\(X\\) back to \\(A\\) along the backward walk is the matrix element of the inverse evolution, and the inverse is the transpose: \\(\\displaystyle\\sum_jw_j=(U^{-1})\_{AX}=(U^{\\mathsf T})\_{AX}=U\_{XA}=Q\_\\varepsilon(t,x)\\). Transposition swaps the endpoints; the weight to come back from a point is the weight to reach it (check 3).
3. *Event weight.* The event “emitted at \\(A\\), absorbed at \\(X\\)” is realized by every outbound path paired with every return path, with no condition linking the two. Its weight is the sum over all pairs, and the double sum factorizes by distributivity alone: \\\[\\sum\_{i,j}w_iw_j=\\Big(\\sum_iw_i\\Big)\\Big(\\sum_jw_j\\Big)=Q\_\\varepsilon(t,x)\\cdot Q\_\\varepsilon(t,x)=Q\_\\varepsilon^2 .\\\] (check 5; by enumeration on an 8-tick lattice, 256 paths each way, equal to \\(Q\_\\varepsilon^2\\) at every point and sense to \\(2\\times10^{-16}\\), `scripts/loop_closure.py`, check 1.)
4. *Over the sense and over space.* \\(\\rho(t,x)=Q\_+^2+Q\_-^2\\), and \\(\\sum_XU\_{XA}^2=(U^{\\mathsf T}U)\_{AA}=1\\) (check 6). The same holds when the tick matrices vary from tick to tick, so long as each is real and orthogonal: \\(U=U_N\\cdots U_1\\), \\(W_X=(e_X^{\\mathsf T}Uf)(f^{\\mathsf T}U^{\\mathsf T}e_X)=q_X^2\\), \\(\\sum_XW_X=\\|f\\|^2\\) (`transaction_lines.py`, check 7: nine different orthogonal ticks, a random real start).

**(D) Three things that follow.**

1. *The count.* Summed over the sense, the event weight through \\(x\\) is \\(Q\_+^2+Q\_-^2=\\rho(x)\\): the count of IV.4d (B) is the weight of being absorbed at \\(x\\) (`loop_closure.py`, check 2).
2. *Conservation.* Summed over \\(x\\), the event weights add to \\(1\\) at every \\(t\\) (ibid., check 4): the walker is absorbed somewhere, once. This is the orthogonality of the step, read as a statement about events.
3. *Exclusivity.* An event has one absorber. Two detectors are two candidate endpoints with weights \\(Q^2\\) at each; one event happens. The amplitude picture, which gives every detector a shift, never had this; the event picture has it by definition of an event.

**(E) Where the signs are.** The weights \\(w_i\\) carry signs, but no event does.

- A pair that returns along its own outbound path has weight \\(w_i^2>0\\). These pairs alone do not give the event weight: at one point of the lattice run they add to \\(0.046\\) against \\(0.320\\) (`loop_closure.py`, check 3). The rest comes from pairs that go out one way and come back another.
- Among all pairs, about half are negative (3304 of 6864 in the run; positives \\(+2.43\\), negatives \\(-1.43\\), total \\(1.00\\); `scripts/negative_pushes.py`, check c). Every negative pair has \\(j\\neq i\\): the two directions disagree in the signs of their U-turns. Those pairs are the interference, and they are inside the sums, never on the event.
- A backward-in-time crossing is not a negative weight: the event weight through \\((x,-)\\) is \\(Q\_-^2\\ge0\\) (check a).
- The two directions push a detector opposite ways: the outbound sum and the return sum are equal, so anything linear in “out minus back” is zero (check b). This is the push-then-pull-back of the gauge-coupled detector (`scripts/detector_as_force.py`), with zero net momentum; what survives is the product, the displacement. A detector cannot read \\(Q\\) linearly because a closed event has no linear part.

**(F) With the added detectors.** Without a field the original two and the added two do not mix, so events close within each two, and the event weight through \\((x,\\varepsilon)\\) is \\(Q\_\\varepsilon^2+Q'^2\_\\varepsilon=|\\psi\_\\varepsilon|^2\\). Closing within each two and adding gives \\(Q\_\\varepsilon^2+Q'^2\_\\varepsilon\\) for any preparation of the pair; the closure imposes no phase relation between the original two and the added two. What the added two are is III.8c (B).

**(G) III.10 is the bookkeeping of where events close.** Differentiate the event weights along the light lines, using III.1’s equation (`scripts/III10_from_closure.py`):

\\\[(\\partial_t+c\\,\\partial_x)\\rho\_+=+\\frac{2mc^2}\\hbar\\sqrt{\\rho\_+\\rho\_-}\\,\\sin\\varphi,\\qquad(\\partial_t-c\\,\\partial_x)\\rho\_-=-\\frac{2mc^2}\\hbar\\sqrt{\\rho\_+\\rho\_-}\\,\\sin\\varphi,\\\]

with \\(\\rho\_\\pm=|\\psi\_\\pm|^2\\) and \\(\\sqrt{\\rho\_+\\rho\_-}\\sin\\varphi=Q\_+Q\_-+Q'\_+Q'\_-\\) in the basis of III.9. These are the equations of III.10, exchange term included. So the walk of III.10 — right-movers and left-movers at \\(\\pm c\\), forward in \\(t\\), one place per run, converting at a rate steered by \\(\\varphi\\) — is the forward-in-\\(t\\) bookkeeping of where the events of III.1 close, and its density is no longer a definition: it is the event weight. This is the link between the two walks of Part III, and between Part I, where \\(\\rho=|\\psi|^2\\) is Madelung’s definition, and Part III, where it is a count of the walk’s own paths in two directions. The phases close the same way. With \\(\\theta\_\\pm=\\operatorname{atan2}(Q'\_\\pm,Q\_\\pm)\\) and \\(W=Q\_+Q'\_--Q'\_+Q\_-\\), \\(\\rho\_\\pm(\\partial_t\\pm c\\,\\partial_x)\\theta\_\\pm=aW\\), which with III.9’s relabelling is III.10’s phase equation (III.8c (C); `III10_from_closure.py` 4, 5).

**(H) The detector, once more.** The clock shift of a detector coupled through the walker’s own turn (`scripts/go_and_come_back_detector.py`) was the overlap of the joint state with itself after the coupling, \\(1-(1-\\cos g)(Q^2+Q'^2)\\) at the detector: going out through the detector and coming back. That is this section’s event weight with the detector as the absorber endpoint, computed from the detector’s side. The two computations are one computation.

*What went in (III.15a).* The walk of III.1, read in \\(t\\); its step is a rotation, so it runs backward as the same walk with the senses swapped (derived). That an event has two endpoints, an emitter and an absorber, and is realized by all outbound and all return paths at once: this is the definition of an event, the walk’s own form of Wheeler–Feynman’s and Cramer’s picture, with everything real. Ord’s entwined pairs (quant-ph/0206016) also pair an outbound path with a return path, and count their net charge, forward minus return; here the two legs *multiply*, and that product is the square. III.15c fixes the degree: a boost stretches a track by \\(k\\) and a reading by \\(\\sqrt k\\), so a frame-independent record is two readings, and (C) says which two. The weight of an event is then \\(Q^2\\) (derived), the weights are positive, exclusive and add to one (derived), and the equations of III.10 follow (derived). That the weight of an event is its frequency is a reading, of the same kind, and at the same place, as III.1’s reading of the step factors as probabilities: one reading, made once. Why the physical absorber reads in \\(t\\) rather than in the walk’s own clock \\(x\\) — in the clock \\(x\\) a detector at a place is passed by every path exactly once and its record is \\(\\int(P\_++P\_-)\\,dt=1\\), linear, with no square — is III.15b: a record is data on a spacelike cut, and the clock cut is not one.

(`scripts/transaction_lines.py` (B)–(C) as matrix identities, `scripts/loop_closure.py` (C)–(E) by enumeration, `scripts/negative_pushes.py` (E), `scripts/copy_is_return_leg.py` (B), `scripts/zigzag_signs_and_twin_checks.py` the sign rule, the three-tick example and the cross pairs, `scripts/III10_from_closure.py` (G), `scripts/checkerboard_is_the_walk.py` (A).)

### III.15b Which reading is physical: the cone and the cut

*The picture first.* The walk can be read in its own clock \\(x\\) or in \\(t\\), and the two readings are one function on the two sides of the light cone (III.15). III.15a and V.5 found that the square, and with it everything a detector does, belongs to the \\(t\\)-reading, and left owed the sentence that says why the \\(t\\)-reading is the physical one. Here is that sentence. It is an input, of the same kind and standing as II.1’s “no preferred frame”: the noise of the walk is the metric, the metric fixes the light cone, and a record must be data on a cut the cone calls spacelike. The walk’s own clock cut is not one.

*Symbols.* None new.

**(A) The noise is the metric, read from either side.** II.1 put it in: the only variance matrix no boost prefers is \\(\\eta^{\\mu\\nu}\\), with the time kick carrying the Lorentzian sign. The Dirac walk has no Gaussian kick, but its two-step moment plays the same part (V.3), and it is the same metric in both readings. In the clock \\(t\\) (IV.4b (H)) the matrices on the three space directions anticommute and square to \\((+1,+1,+1)\\):

\\\[\\partial_t^2Q=c^2\\big(\\partial_x^2+\\partial_y^2+\\partial_z^2\\big)Q-a^2Q .\\\]

In the clock \\(z\\) (IV.4b (F)) the matrices on \\((t,x,y)\\) are \\((A_z,\\,A_xA_z,\\,A_yA_z)\\), anticommuting and squaring to \\((+1,-1,-1)\\):

\\\[c^2\\partial_z^2Q=\\partial_t^2Q-c^2\\partial_x^2Q-c^2\\partial_y^2Q+a^2Q .\\\]

The same equation, rearranged: in the clock \\(z\\), \\(t\\) spreads with the sign a space direction had in the clock \\(t\\), and \\(x,y\\) with the opposite sign. The cone, \\(c^2t^2=x^2+y^2+z^2\\), is the same on both lines; what moves is which coordinate is the clock, and the clock’s sign goes with it. That is the Lorentzian sign of III.9a, once more.

**(B) Two cuts, two kernels.** A reading evolves data off a cut. The \\(t\\)-reading evolves data given on \\(t=\\text{const}\\), a spacelike surface; the walk’s own clock evolves data given on \\(x=\\text{const}\\), a timelike surface — a worldline. III.14 and III.15 computed both kernels from one source at \\((0,0)\\):

\\\[x\\text{-reading (III.14):}\\quad\\text{support }x\\ge c|t|,\\ \\ x^2-c^2t^2\\ge0:\\ \\text{spacelike or null from the source, including }t\<0;\\\]

\\\[t\\text{-reading (III.15):}\\quad\\text{support }ct\\ge|x|,\\ \\ c^2t^2-x^2\\ge0,\\ t>0:\\ \\text{the future cone.}\\\]

Concretely, a detector sitting at a place \\(x_D\\) and reading the walk in its own clock is crossed by every path exactly once, \\(\\int(P\_++P\_-)\\,dt=1\\) (V.5 (F)), but a finite part of that record sits at \\(t\<0\\): for \\(x_D=c/a\\) it is \\(27\\%\\), for \\(x_D=3c/a\\) it is \\(38\\%\\) (`which_reading_is_physical.py` 1). Those crossings are the backward-in-\\(t\\) segments of III.1 — the forward-time observer’s antiparticles. A detector that read the walk in the walk’s own clock would therefore click, a third of the time, **in the past of the trigger**: the source fires at \\(t=0\\) and the counter at \\(x_D\\) has already registered. No counter does that. The \\(t\\)-reading has no weight there at all.

**(C) The input.** *A record is data on a spacelike cut: what a detector registers at an event depends on the source only through the event’s past cone.* Data on a timelike cut is not a causal problem — the \\(x\\)-reading would correlate the detector with the source at spacelike separation, the click before the trigger of (B). With this sentence the \\(t\\)-reading is the physical one, and III.15a follows: the signed weights of \\(T\\), the two legs, the square, and V.5’s detector that can register nothing else. The cut is used three times in the document, and it is one statement: the source and the detector each sit on a time slice; the measurement pairs the outward and return routes that agree there (III.8c (A)); the completed detectors read there (III.8c (B)); and a record is data there (this section). In the walk’s own clock the same lock gives only \\(\\sin\\varphi\\ge0\\) (III.8c (B), item 6): the cut matters for the state, not only for the counter.

**(D) What it does and does not do.** It does not derive the cone from the coin; the cone is \\(\\eta\\), an input, and this is the same input read once more: the metric says which cuts are spacelike, and records live on those. It answers the question as it was asked — not “why do detectors wait in time”, but “why is the signed reading the one with the physics” — and the answer is that the other reading’s support is spacelike. One thing to keep apart: the *source* answer (III.16, IV.10’s bottom row) does not vanish outside the cone, it decays there; the statement here is about the evolution of a state from a cut, which is what a detection is (V.5 evolves the joint state in \\(t\\)), and the source problem is a different question, treated as one.

*What went in.* (A), (B): the matrices of IV.4b and the kernels of III.14–III.15, checked. (C): the one input of this section, the causal cut, of the same standing as \\(\\eta^{\\mu\\nu}\\) in II.1. It replaces the posit in the Part III ledger — “the signed \\(t\\)-reading is the physical one” — with the sentence that selects it.

(`scripts/which_reading_is_physical.py`: 1, 2, 3a, 3b.)

### III.15c Why the record is a square: light-cone tracks and frame independence

(Primes in this section mark the boosted readings, \\(Q'\_\\pm\\); they are not V.2’s added detectors \\(Q'\\).)

*The picture first.* III.15a found the weight of an event to be the product of its two legs, \\(Q^2\\); III.15b said which cut a record lives on. Here is why the record on that cut has to be a square. The particle moves at \\(c\\): its two senses are its two null tracks, and an event’s position is the pair of track coordinates it sits on. A boost stretches one track coordinate by \\(k\\) and shrinks the other by \\(1/k\\); it stretches the two readings by \\(\\sqrt k\\) and \\(1/\\sqrt k\\). A count of crossings has to go as a track coordinate, so a count is a product of two readings — and III.15a already names them: the out reading and the return reading, which are one number. \\(c=1\\).

*Symbols.* \\(u=t-x\\), \\(v=t+x\\); \\(t_1=u(X)\\), \\(t_2=v(X)\\) for an event \\(X\\); \\(k=e^\\eta\\). Nothing else is new.

**(A) The two senses are the two null tracks.** III.1: \\(dt=\\varepsilon\\,dx\\), so

\\\[du=(\\varepsilon-1)\\,dx,\\qquad dv=(\\varepsilon+1)\\,dx :\\tag{1}\\\]

1. a \\(+\\) segment keeps \\(u\\) and advances in \\(v\\); a \\(-\\) segment keeps \\(v\\) and advances in \\(u\\). An event is where a \\(+\\) track and a \\(-\\) track cross, and its two track coordinates give its position and its interval:

\\\[t=\\tfrac{t_1+t_2}2,\\qquad x=\\tfrac{t_2-t_1}2,\\qquad t^2-x^2=t_1t_2,\\qquad\\zeta=a\\sqrt{x^2-t^2}=a\\sqrt{-t_1t_2}.\\tag{2}\\\]

2. III.14’s support \\(x\\ge|t|\\) is \\(t_1\\le0\\le t_2\\); III.15’s is \\(t_1,t_2\\ge0\\). The argument of both kernels is the flip rate times the geometric mean of the two track coordinates. (`bondi_and_the_two_senses.py` 1–2, `radar_and_the_two_legs.py` 2.)

**(B) A boost stretches a track by \\(k\\) and a reading by \\(\\sqrt k\\).** The lines with the loss out, read in \\(t\\) (III.5a (H)), in null coordinates, \\(\\partial_t\\pm\\partial_x=2\\partial\_{v,u}\\):

\\\[2\\,\\partial_vQ\_+=a\\,Q\_-,\\qquad2\\,\\partial_uQ\_-=-a\\,Q\_+ .\\tag{3}\\\]

3. each sense advances along its own track and is fed by the other. A boost, \\(t'=\\gamma(t-\\beta x)\\), \\(x'=\\gamma(x-\\beta t)\\):

\\\[u'=\\gamma(1+\\beta)\\,u=k\\,u,\\qquad v'=\\gamma(1-\\beta)\\,v=\\frac vk,\\qquad t_1't_2'=t_1t_2 .\\tag{4}\\\]

4. Bondi’s \\(k\\)-calculus, cited for the coordinates and the factor. (`bondi_and_the_two_senses.py` 3.) Put \\(\\partial_v=\\frac1k\\partial\_{v'}\\), \\(\\partial_u=k\\,\\partial\_{u'}\\) and allow the readings a factor each, \\(Q\_+=\\lambda Q'\_+\\), \\(Q\_-=\\mu Q'\_-\\):

\\\[2\\,\\partial\_{v'}Q'\_+=a\\,\\frac{k\\mu}\\lambda\\,Q'\_-,\\qquad2\\,\\partial\_{u'}Q'\_-=-a\\,\\frac\\lambda{k\\mu}\\,Q'\_+ .\\tag{5}\\\]

5. 4. into (3), divided through. These are (3) again iff \\(\\lambda/\\mu=k\\). The scale is fixed by what a boost is: two boosts compose by adding rapidities, so \\(\\lambda=e^{\\alpha\\eta}\\); the mirror \\(x\\to-x\\) swaps the senses and reverses the boost, so \\(\\mu(\\eta)=\\lambda(-\\eta)\\) and \\(\\lambda\\mu=1\\). Then \\(e^{2\\alpha\\eta}=e^\\eta\\):

\\\[Q'\_+=\\frac{Q\_+}{\\sqrt k},\\qquad Q'\_-=\\sqrt k\\,Q\_-,\\qquad Q'\_+Q'\_-=Q\_+Q\_- .\\tag{6}\\\]

6. half the rapidity on each reading; the product of the two senses is invariant, as \\(t_1t_2\\) is — it is the mass term’s weight (III.10). (`half_rapidity_from_the_lines.py` 1–4; IV.4d (D) reaches (6) from the content.)

**(C) A count goes as a track coordinate, so a record is two readings.** IV.4d (B)’s conserved content in null components, with \\(\\partial_t=\\partial_u+\\partial_v\\), \\(\\partial_x=\\partial_v-\\partial_u\\):

\\\[j^u=j^0-j^1=2Q\_-^2\\ \\to\\ k\\,(2Q\_-^2),\\qquad j^v=j^0+j^1=2Q\_+^2\\ \\to\\ \\frac{2Q\_+^2}k .\\tag{7}\\\]

7. by (6). A record on a stretch of \\(-\\) track (constant \\(v\\), length \\(du\\)) is the number of \\(+\\) segments crossing it, \\(j^v\\) times the surface element \\(\\tfrac12du\\) (since \\(dt\\,dx=\\tfrac12du\\,dv\\)):

\\\[\\int Q\_+^2\\,du\\ \\to\\ \\int\\frac{Q\_+^2}k\\,k\\,du=\\int Q\_+^2\\,du :\\tag{8}\\\]

8. the same number in every frame. That is what frame independence of a count is. Two things fix the degree, and neither does it alone. First, a record is a *conserved local density* — a count on a cut that the next cut agrees with — and on the lines, which are linear and first order, the only densities without derivatives whose flux closes are the quadratic forms. This is a theorem, three lines long. Ask for \\(f(Q\_+,Q\_-)\\) and a flux \\(g(Q\_+,Q\_-)\\), no derivatives of the readings and no explicit \\(x,t\\), with \\(\\partial_tf+\\partial_xg=0\\) on every solution of the lines \\(\\partial_tQ\_\\pm=\\mp c\\,\\partial_xQ\_\\pm\\pm aQ\_\\mp\\). Substituting, \\((g_u-cf_u)\\,\\partial_xQ\_++(g_v+cf_v)\\,\\partial_xQ\_-+a\\,(vf_u-uf_v)=0\\) with \\(u=Q\_+\\), \\(v=Q\_-\\); the values and the slopes at a point are independent, so all three brackets vanish. The first two give \\(cf\_{uv}=-cf\_{uv}\\), so \\(f=F(u)+G(v)\\); the third then reads \\(vF'(u)=uG'(v)\\), so \\(F'(u)/u=G'(v)/v=k\\): \\(f=\\tfrac k2(Q\_+^2+Q\_-^2)\\) up to a constant, \\(g=\\tfrac{ck}2(Q\_+^2-Q\_-^2)\\) (`quadratic_density_theorem.py`, symbolically to degree four). With \\(a=0\\) the third bracket is absent and any \\(F(Q\_+)+G(Q\_-)\\) is conserved: it is the mass that pairs the two readings into one square. (A quartic such as \\(Q\_+^3Q\_-\\) scales as \\(1/k\\) exactly like \\(Q\_+^2\\); the boost alone would not exclude it, the theorem does.) Second, among the quadratic forms the boost says which one pairs with which null length: \\(Q\_+^2\\sim1/k\\) goes with \\(du\\), \\(Q\_-^2\\sim k\\) with \\(dv\\), and \\(Q\_+Q\_-\\sim1\\) with neither — it is the mass density, the one scalar two readings make (`counting_in_spacetime.py`). So a record is bilinear in readings. The two readings it is made of are III.15a (C)1–2 — the out sum \\(U\_{XA}\\) and the return sum \\((U^{\\mathsf T})\_{AX}\\), one number because the step is orthogonal:

\\\[W(X)=U\_{XA}\\,(U^{\\mathsf T})\_{AX}=U\_{XA}^2=Q\_\\varepsilon(X)^2,\\qquad\\sum_XW(X)=(U^{\\mathsf T}U)\_{AA}=1 .\\tag{9}\\\]

9. III.15a (C)3–4, with the degree now forced by (8). In one sentence: a reading is the square root of a track coordinate’s worth of crossings; an event has two track coordinates, so it has two readings, and III.15a says which two.

*What this leaves.* The square is the walk’s agreement with its own return: for a coin that agreement decays (the walk forgets), for this walk it is held at one on every cut, because the orthogonal step retraces (`symmetric_part_and_coincidence.py` G1–G3). Why the return leg exists on the physical cut and not on the walk’s own clock is III.15b. What remains an input, the only one left in III.15a: that the weight of an event is its frequency — a number read as a probability, at the place III.1 read the step factors that way.

*What went in.* No new input. (A) is III.1 in null coordinates; (B) is the Lorentz transformation on the lines, with the scale fixed by composition and the mirror; (C) is IV.4d (B) in null components and III.15a (C). Not derived: count = frequency.

(`scripts/bondi_and_the_two_senses.py` (A), (4); `scripts/radar_and_the_two_legs.py` (A); `scripts/quadratic_density_theorem.py` (8), the theorem; `scripts/half_rapidity_from_the_lines.py`, `scripts/sqrt_k_as_an_eigenvalue_problem.py` (B); `scripts/frame_independence_forces_the_square.py`, `scripts/counting_in_spacetime.py` (C); `scripts/reversal_makes_the_return_leg.py`, `scripts/symmetric_part_and_coincidence.py` for the return leg and the agreement.)

### III.15d A count along the walk’s own clock, and why its weight is a boost

*The picture first.* III.14’s kernel says where one zigzag is, in \\(t\\), when the clock reads \\(x\\). Ask a counting question of it: of the zigzags that have advanced to \\(x\\), how many have come back to the sense they started in? With the loss out the answer grows along the clock, at the rate \\(a\\). Count the same zigzags with a weight on how far each has netted in \\(t\\) and the rate becomes \\(\\omega_k=\\sqrt{a^2+c^2k^2}\\). This section shows that the weighted count is not a new object: it is the plain count seen from a moving frame. The rest zigzag spends equal clock-weight in both senses and grows as \\(e^{ax/c}\\); a boost shears that growth partly into the new clock and partly into a slope in the new time, and stretches the two null tracks oppositely, so the moving observer sees one sense heavier than the other. The slope is the weight, \\(k\\) is the momentum of the frame, and the growth rate is the energy over \\(\\hbar\\).

*Symbols.* \\(\\tau=x/c\\) the walk’s clock; \\(Q\_{\\varepsilon'\\varepsilon}(t;x)\\) the kernel of III.14 with the loss out, started in sense \\(\\varepsilon\\) at \\(t=0\\) and read in sense \\(\\varepsilon'\\); \\(k\\) a real number; \\(\\eta\\) a rapidity. Nothing else is new.

**(A) III.14’s kernel, written out.** Start at \\(x=0\\), \\(t=0\\) in the forward sense and advance the clock to \\(x\\). III.14 gives where the walker is in \\(t\\), and in which sense:

\\\[Q\_{++}(t;x)=\\delta(t-\\tau)+\\frac a2\\sqrt{\\frac{\\tau+t}{\\tau-t}}\\;I_1\\!\\big(a\\sqrt{\\tau^2-t^2}\\big),\\qquad Q\_{-+}(t;x)=\\frac a2\\,I_0\\!\\big(a\\sqrt{\\tau^2-t^2}\\big),\\qquad |t|\<\\tau,\\tag{1}\\\]

1. the \\(\\delta\\) is the path that never flipped; \\(I_1\\) collects the even-flip paths, back in the forward sense having lost ground in \\(t\\); \\(I_0\\) the odd-flip paths, now running backward. By the mirror, \\(Q\_{--}(t;x)=\\delta(t+\\tau)+\\frac a2\\sqrt{\\frac{\\tau-t}{\\tau+t}}\\,I_1\\) and \\(Q\_{+-}=Q\_{-+}\\). Every weight is positive, and with the loss put back \\(e^{-a\\tau}\\int dt\\,\[Q\_{++}+Q\_{-+}\]=1\\): it is the distribution of where one zigzag is at clock \\(x\\) (`sea_in_coordinates_check.py` 1).

**(B) The count, and the weighted count.** The number of zigzags back in their starting sense, both starts:

\\\[C_0(\\tau)=\\int dt\\,\\big\[Q\_{++}(t;x)+Q\_{--}(t;x)\\big\]=2\\cosh(a\\tau).\\tag{2}\\\]

2. it grows. The reading in \\(x\\) conserves a linear count, not a length (III.5a (J), III.15b); with the survival factor out the number of alternatives grows as \\(e^{a\\tau}\\), which is III.5a (F)’s \\(Q=e^{ax/c}P\\). The reading in \\(t\\) cannot do this: it conserves \\(\\int(Q\_+^2+Q\_-^2)\\) and only turns. Now weight each zigzag by \\(e^{-ck\\,\\Delta t}\\), \\(\\Delta t\\) its net displacement in \\(t\\) (ups minus downs):

\\\[L_k(\\tau)\_{\\varepsilon'\\varepsilon}=\\int dt\\,e^{-ckt}\\,Q\_{\\varepsilon'\\varepsilon}(t;x),\\\]

\\\[L_k(\\tau)=\\begin{pmatrix}\\cosh\\omega_k\\tau-\\dfrac{ck}{\\omega_k}\\sinh\\omega_k\\tau & \\dfrac a{\\omega_k}\\sinh\\omega_k\\tau\\\\\[6pt\] \\dfrac a{\\omega_k}\\sinh\\omega_k\\tau & \\cosh\\omega_k\\tau+\\dfrac{ck}{\\omega_k}\\sinh\\omega_k\\tau\\end{pmatrix}=e^{\\tau M_k},\\qquad M_k=\\begin{pmatrix}-ck&a\\\\a&ck\\end{pmatrix},\\qquad M_k^2=\\omega_k^2\\begin{pmatrix}1&0\\\\0&1\\end{pmatrix}.\\tag{3}\\\]

3. a Laplace transform in \\(t\\) of a real kernel with bounded support, \\(|t|\\le\\tau\\), so it converges for every \\(k\\) of either sign; no continuation is made and no complex number appears (ibid. 2, the integrals against the closed form to \\(10^{-15}\\)). It is III.14’s own \\(\\hat Q(\\omega;x)\\), \\(\\kappa=\\sqrt{a^2-\\omega^2}\\), read where \\(\\kappa=\\omega_k\\). The same-sense count is

\\\[C_k(\\tau)=\\operatorname{tr}L_k(\\tau)=2\\cosh(\\omega_k\\tau),\\qquad\\frac1\\tau\\log C_k(\\tau)\\to\\omega_k .\\tag{4}\\\]

4. *\\(\\omega_k\\), the frequency of the \\(t\\)-reading, is the growth rate along the walk’s own clock of the count of returning zigzags weighted by \\(e^{-ck\\Delta t}\\)* (ibid. 3). That the trace is the sum over where the turns are placed takes three lines: a history with \\(2m\\) turns is an ordered run of \\(m+1\\) stretches weighted \\(e^{\\mp cku}\\) and \\(m\\) weighted \\(e^{\\pm cku}\\); the Laplace transform in \\(\\tau\\) of an ordered convolution is a product, so the \\(2m\\)-turn term is \\(a^{2m}\\big\[(p+ck)^{-(m+1)}(p-ck)^{-m}+(p-ck)^{-(m+1)}(p+ck)^{-m}\\big\]=\\frac{2p}{p^2-c^2k^2}\\big(\\frac{a^2}{p^2-c^2k^2}\\big)^m\\), and the geometric sum is \\(\\frac{2p}{p^2-\\omega_k^2}\\), the transform of \\(2\\cosh\\omega_k\\tau\\) (`audit_gluing_draft_2026_11_09.py` 2). With the survival factor left in, the same number is an average under III.1 itself:

\\\[\\tfrac12e^{-a\\tau}\\,C_k(\\tau)=\\mathbb E_P\\!\\left\[e^{-ck\\,\\Delta t}\\,\\mathbf 1\_{\\varepsilon(\\tau)=\\varepsilon(0)}\\right\],\\qquad \\omega_k=a+\\lim\_{\\tau\\to\\infty}\\frac1\\tau\\log\\mathbb E_P\\!\\left\[e^{-ck\\,\\Delta t}\\,\\mathbf 1\_{\\varepsilon(\\tau)=\\varepsilon(0)}\\right\].\\tag{5}\\\]

5. launch in either sense, run the clock to \\(\\tau\\), keep the histories that return to their starting sense, average \\(e^{-ck\\Delta t}\\). Near rest,

\\\[\\omega_k-a=\\frac{c^2k^2}{2a}+O(k^4)=\\frac{\\hbar k^2}{2m}+\\ldots,\\qquad\\frac{c^2}{2a}=\\frac\\hbar{2m}:\\tag{6}\\\]

6. the Schrödinger kinetic term, with Part I’s diffusion constant \\(D=\\hbar/2m\\) as the coefficient: the growth rate of the exponential moment of the net time advance is, near rest, the diffusive \\(Dk^2\\) of the zigzag’s wandering in \\(t\\) along its own clock. The closure in (4) is in *sense*; no path that always advances in \\(x\\) returns to its \\(x\\).

**(C) The weight is a boost of the plain count.** The lines in the walk’s clock, \\(c\\,\\partial_xQ=(-A_z\\partial_t+a\\sigma_x)Q\\) with \\(A_z=\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}\\) and \\(\\sigma_x=\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix}\\), have at \\(k=0\\) the growing solution

\\\[Q(t,x)=\\begin{pmatrix}1\\\\1\\end{pmatrix}e^{ax/c}.\\tag{7}\\\]

7. the plain count of (2): both senses equal, no dependence on \\(t\\). Boost it with III.15c: rapidity \\(\\eta\\), \\(v=c\\tanh\\eta\\), \\(x=\\gamma(x'+vt')\\), \\(t=\\gamma(t'+vx'/c^2)\\), and III.15c (6), half the rapidity on each reading:

\\\[\\begin{pmatrix}Q'\_+\\\\Q'\_-\\end{pmatrix}=\\begin{pmatrix}e^{-\\eta/2}&0\\\\0&e^{+\\eta/2}\\end{pmatrix}\\begin{pmatrix}Q\_+\\\\Q\_-\\end{pmatrix},\\\]

\\\[Q'(t',x')=\\begin{pmatrix}e^{-\\eta/2}\\\\e^{+\\eta/2}\\end{pmatrix}\\exp\\Big\[\\frac ac\\gamma(x'+vt')\\Big\]=\\begin{pmatrix}e^{-\\eta/2}\\\\e^{+\\eta/2}\\end{pmatrix}e^{\\omega x'/c}\\,e^{ck\\,t'},\\qquad\\omega=a\\cosh\\eta,\\quad ck=a\\sinh\\eta.\\tag{8}\\\]

8. this column solves the primed lines exactly (`tilt_is_a_boost.py` 1). The invariant exponent \\(ax/c\\) has become growth along the new clock at \\(\\omega\\) and a slope in the new time at \\(ck\\); the column is one-sided, \\(Q'\_-/Q'\_+=e^{\\eta}\\). Now compare with (3). A solution that slopes as \\(e^{+ckt}\\) obeys \\(c\\,\\partial_xq=M_kq\\), the same matrix the Laplace weight evolves by, and \\(M_k\\) has the growing eigenvalue \\(\\omega_k\\) with eigenvector

\\\[\\begin{pmatrix}1\\\\(\\omega_k+ck)/a\\end{pmatrix}=\\begin{pmatrix}1\\\\e^{\\eta}\\end{pmatrix}\\qquad\\text{when}\\quad ck=a\\sinh\\eta,\\ \\ \\omega_k=a\\cosh\\eta:\\tag{9}\\\]

9. the boosted column of (8), entry for entry (ibid. 2). With \\(a=mc^2/\\hbar\\),

\\\[\\hbar k=mc\\sinh\\eta=p,\\qquad\\hbar\\omega_k=mc^2\\cosh\\eta=E,\\qquad\\omega_k^2-c^2k^2=a^2\\ \\Longleftrightarrow\\ \\cosh^2\\eta-\\sinh^2\\eta=1 .\\tag{10}\\\]

10. the tilt is the momentum and the growth the energy of the rest configuration in motion (ibid. 3); the imbalance is III.15c’s, the sense running against the observer’s motion stretched and the other compressed. The decaying eigenvector, \\((1,-e^{-\\eta})\\) with rate \\(-\\omega_k\\), is the same statement for the rest solution \\((1,-1)e^{-ax/c}\\) seen from the frame \\(-\\eta\\). One caveat, kept exact: the boost does not turn \\(M_0\\) into \\(M_k\\) — their eigenvalues differ by \\(\\cosh\\eta\\) — because the clock changes with the frame. It is the solutions that map, not the step matrices (ibid. 4).

**(D) Path by path.** A forward segment of rest clock \\(d\\) has primed clock \\(d\\,e^{-\\eta}\\), a backward one \\(d\\,e^{+\\eta}\\) (III.15c (A)–(B)), so the moving observer sees the flip rate \\(a\\,e^{+\\eta}\\) while the walker runs forward and \\(a\\,e^{-\\eta}\\) while it runs backward. Against the one-rate walk he would write in his own clock, every zigzag’s weight differs by

\\\[\\frac{\\text{boosted path weight}}{\\text{rate-}a\\text{ path weight}}=e^{-ck\\,\\Delta t'}\\;e^{-(\\omega_k-a)S}\\;e^{\\eta(\\varepsilon_i-\\varepsilon_f)/2},\\qquad S=S\_++S\_-,\\ \\ \\Delta t'=S\_+-S\_-,\\tag{11}\\\]

11. from \\(-a(e^{\\eta}-1)S\_+-a(e^{-\\eta}-1)S\_-=-(\\omega_k-a)S-ck\\,\\Delta t'\\), with the flip factors \\(e^{\\pm\\eta}\\) per completed segment collapsing to the two ends. \\(S\_\\pm\\) is the primed clock spent in each sense, \\(\\varepsilon_i,\\varepsilon_f\\) the first and last sense. So the weight of (B) is the ratio of two path measures, on every path and not only on the growing solution; the uniform \\(e^{-(\\omega_k-a)S}\\) is the loss the moving observer must take out (\\(\\omega_k=\\gamma a\\), V.4’s rate at a fixed position); the end factors are III.15c (6). In matrices, with \\(B=\\begin{pmatrix}e^{-\\eta/2}&0\\\\0&e^{+\\eta/2}\\end{pmatrix}\\): the boosted kernel is \\(e^{-(\\omega_k-a)S}\\,B\\,L_k\\,B^{-1}\\), which is why \\(M_k\\)’s growing eigenvector is \\(B\\,(1,1)^{\\mathsf T}\\); and the weighted count of the one-rate walk is the plain count of the two-rate walk the moving observer sees, whose generator

\\\[-\\omega_k\\begin{pmatrix}1&0\\\\0&1\\end{pmatrix}+B\\,M_k\\,B^{-1}=\\begin{pmatrix}-a\\,e^{\\eta}&a\\,e^{-\\eta}\\\\a\\,e^{\\eta}&-a\\,e^{-\\eta}\\end{pmatrix}\\tag{12}\\\]

12. has columns that sum to zero: a process, with probabilities. (`boost_path_weight.py`: 2000 random zigzags, the ratio equals (11) to \\(10^{-12}\\), and the variant with \\(\\varepsilon_i+\\varepsilon_f\\) fails; `audit_2026_11_09.py` F for (12).)

*What went in.* No new input. (A)–(B) are III.14’s kernel and integrals of it; (C)–(D) are III.15c’s boost applied to the plain count. What (C) gives a name to — “energy” for the growth rate — is III.8c (D)’s definition of energy as \\(\\hbar\\) times a turn rate in \\(t\\); the two rates are one number because the \\(t\\)-generator is a quarter turn times the clock’s (III.15e (4)). The definition itself, with \\(\\hbar\\) as the scale (input 17), is not derived. VI.5 uses (3)–(4) for many walkers.

(`scripts/sea_in_coordinates_check.py` 1–3 for (A)–(B); `scripts/tilt_is_a_boost.py`, 24 checks at six rapidities, for (C), and `scripts/master_v33_new_text_checks.py` 4–5 for the decaying slot and (12); `scripts/boost_path_weight.py` for (11); `scripts/audit_gluing_draft_2026_11_09.py` 2 for the turn-placement sum; `scripts/audit_2026_11_09.py` F for (12).)

### III.15e One system, isolated both ways: the original two detectors, and the four

*The picture first.* One zigzag read by III.1’s two detectors is a real column of two counters. Isolating the clock derivative gives a count that grows; isolating the time derivative gives a pair of counters that turns. The original two turn one way only, so they have a frequency and no sign of frequency: read by them alone the walker is its own antiparticle. Read by all four detectors of III.8c (B), the walker’s clock can turn one way or the other: the same amount in the forward and the backward detector at each position of the hand, or opposite amounts. These are the two loops, leading and lagging. This section writes the four numbers as one real column, isolates the two derivatives, and shows which combination of the walk’s own count corresponds to which sense of the clock. No complex number is used: the only quarter turn besides the walk’s own \\(T\\) is the walker’s own turn, \\(\\mathbb J\\).

*Symbols.* Four real numbers, the original two first, \\(r=(u\_+,u\_-,v\_+,v\_-)^{\\mathsf T}\\): III.1’s two detectors \\(u\\) and the two added ones \\(v\\) (III.8c (B)), each a forward and a backward counter. \\(a=mc^2/\\hbar\\), \\(\\omega_k=\\sqrt{a^2+c^2k^2}\\). Four fixed matrices:

\\\[\\mathbb A=\\begin{pmatrix}1&0&0&0\\\\0&-1&0&0\\\\0&0&1&0\\\\0&0&0&-1\\end{pmatrix},\\quad\\mathbb S=\\begin{pmatrix}0&1&0&0\\\\1&0&0&0\\\\0&0&0&1\\\\0&0&1&0\\end{pmatrix},\\quad\\mathbb T=\\begin{pmatrix}0&1&0&0\\\\-1&0&0&0\\\\0&0&0&1\\\\0&0&-1&0\\end{pmatrix},\\quad\\mathbb J=\\begin{pmatrix}0&0&-1&0\\\\0&0&0&-1\\\\1&0&0&0\\\\0&1&0&0\\end{pmatrix}.\\\]

\\(\\mathbb A\\) is the advance (the backward sense runs down in \\(t\\)), \\(\\mathbb S\\) the exchange and \\(\\mathbb T\\) the turn, each the same \\(2\\times2\\) matrix on the original two and on the added two; \\(\\mathbb J\\) is the walker’s own turn, each \\(u\_\\varepsilon\\) onto \\(v\_\\varepsilon\\) and each \\(v\_\\varepsilon\\) onto \\(-u\_\\varepsilon\\). \\(\\mathbb 1\\) is the \\(4\\times4\\) unit matrix.

**(A) The system, isolated both ways.** The walk’s clock, loss out (III.5a (F)), once for the original two and once for the added two:

\\\[c\\,\\partial_x\\!\\begin{pmatrix}u\_+\\\\u\_-\\\\v\_+\\\\v\_-\\end{pmatrix}=\\begin{pmatrix}-\\partial_t&a&0&0\\\\a&\\partial_t&0&0\\\\0&0&-\\partial_t&a\\\\0&0&a&\\partial_t\\end{pmatrix}\\!\\begin{pmatrix}u\_+\\\\u\_-\\\\v\_+\\\\v\_-\\end{pmatrix},\\qquad\\text{i.e.}\\quad c\\,\\partial_xr=-\\mathbb A\\,\\partial_tr+a\\,\\mathbb S\\,r .\\tag{1}\\\]

1. III.5a (A)’s two lines, written twice. Move the \\(\\partial_t\\) term to one side, multiply by \\(\\mathbb A\\) (\\(\\mathbb A^2=\\mathbb 1\\)), and do the one product \\(\\mathbb A\\mathbb S=\\mathbb T\\):

\\\[\\partial_tr=-c\\,\\mathbb A\\,\\partial_xr+a\\,\\mathbb T\\,r,\\qquad\\partial_t\\!\\begin{pmatrix}u\_+\\\\u\_-\\\\v\_+\\\\v\_-\\end{pmatrix}=\\begin{pmatrix}-c\\partial_x&a&0&0\\\\-a&c\\partial_x&0&0\\\\0&0&-c\\partial_x&a\\\\0&0&-a&c\\partial_x\\end{pmatrix}\\!\\begin{pmatrix}u\_+\\\\u\_-\\\\v\_+\\\\v\_-\\end{pmatrix}.\\tag{2}\\\]

2. the same two rows for each two; one row has changed sign (III.5a (J)). What that does to the turns is III.15a (A)’s rule:

| turn | walk’s clock | \\(t\\)-reading |
| --- | --- | --- |
| \\(+\\to-\\) | \\(+a\\) | \\(-a\\) |
| \\(-\\to+\\) | \\(+a\\) | \\(+a\\) |

The sign sits on one orientation of turn, not on time spent in the backward sense. Which orientation carries it is a convention (\\(Q\_-\\to-Q\_-\\) moves it); the product over a pair of turns is not:

\\\[\\mathbb S^2=+\\mathbb 1\\ \\Rightarrow\\ \\cosh,\\ \\sinh\\ \\ (\\text{the count grows}),\\qquad\\mathbb T^2=-\\mathbb 1\\ \\Rightarrow\\ \\cos,\\ \\sin\\ \\ (\\text{the counters turn}).\\\]

And \\(\\mathbb J^2=-\\mathbb 1\\); it commutes with \\(\\mathbb A\\), \\(\\mathbb S\\) and \\(\\mathbb T\\), because it moves whole blocks and they act alike inside each. \\(\\mathbb T\\) anticommutes with \\(\\mathbb A\\). Two matrices square to \\(-\\mathbb 1\\); only \\(\\mathbb J\\) commutes with the whole step.

**(B) The original two alone: a frequency, no sign, half a walker.** Set \\(v=0\\). The \\(t\\)-reading turns the column in the sense plane at \\(\\omega_k\\): the dispersion relation, and \\(|E|=\\hbar\\omega_k\\) (III.5 remark 4, III.8c (D)). Every real column turns the same way, so nothing on the original two can carry a sign of energy. Three statements:

- No *local* real matrix on the original two commutes with \\(A_z\\) and \\(T\\) except multiples of the unit matrix (III.8b).
- A *non-local* one exists, III.8c (15): \\(J(k)=G(k)/\\omega_k\\), real, commuting with the step, \\(J(k)^2=-1\\). Measured against it the original two have one sign of frequency only, \\(\\partial_tu=\\omega_kJ(k)u\\) for every \\(u\\).
- So a walker read by the original two alone is self-conjugate — its own antiparticle — and it is half a walker: filled to its lowest configuration it has, at rest, energy \\(-\\tfrac12\\hbar a\\); read by all four, \\(-\\hbar a\\) (VI.4 (A), VI.6 (A)).

A walker in the original two alone can still be acted on, through its coin: V.1’s scalar coupling \\(a\\to a+V/\\hbar\\). What it cannot be is turned.

**(C) All four: the two senses of the clock.** With all four there is something to compare. At rest (2) is \\(\\partial_tr=a\\,\\mathbb T\\,r\\), and since \\(\\mathbb T\\) and \\(\\mathbb J\\) commute and both square to \\(-\\mathbb 1\\), the four-number space splits into two classes:

\\\[\\begin{aligned}\\mathbb T\\,r&=-\\mathbb J\\,r:&v&=+Tu&&(\\text{forward and backward the same}),\\\\\\mathbb T\\,r&=+\\mathbb J\\,r:&v&=-Tu&&(\\text{forward and backward opposite}).\\end{aligned}\\tag{3}\\\]

3. the same amount, or opposite amounts, in the forward and the backward detector of each hand position. For a walker launched forward with the hand at the start, \\(u=(1,0)\\), the two are

\\\[\\text{leading loop }\\lambda=\\begin{pmatrix}1\\\\0\\\\0\\\\-1\\end{pmatrix},\\qquad\\text{lagging loop }\\ell=\\begin{pmatrix}1\\\\0\\\\0\\\\1\\end{pmatrix},\\qquad\\ell=\\begin{pmatrix}1&0&0&0\\\\0&1&0&0\\\\0&0&-1&0\\\\0&0&0&-1\\end{pmatrix}\\lambda\\ \\ (\\text{the added two counted with the other sign, }v\\to-v).\\\]

In the detectors of III.8c (B): the exchange of the forward and the backward detector of one hand position is \\(\\mathbb J\\,\\mathbb T\\), with square \\(+1\\), and the two classes are \\(\\mathbb J\\mathbb T\\,r=+r\\) and \\(\\mathbb J\\mathbb T\\,r=-r\\). With the half-sums \\(\\tfrac12(u\_+-v\_-)\\), \\(\\tfrac12(v\_++u\_-)\\) and the half-differences \\(\\tfrac12(u\_++v\_-)\\), \\(\\tfrac12(v\_+-u\_-)\\) of the two detectors of each hand position, the leading loop has sums only and the lagging loop differences only (`detectors_in_master_names.py` 5).

Written out in time, with \\(u(t)=(\\cos at,\\,-\\sin at)\\) for both (III.5a (K)):

\\\[\\lambda:\\quad v(t)=\\begin{pmatrix}-\\sin at\\\\-\\cos at\\end{pmatrix}=u\\Big(t+\\frac\\pi{2a}\\Big),\\qquad\\qquad\\ell:\\quad v(t)=\\begin{pmatrix}\\sin at\\\\\\cos at\\end{pmatrix}=u\\Big(t-\\frac\\pi{2a}\\Big).\\\]

The added two hold what the original two hold a quarter period on, or a quarter period back. The original two alone, \\(v=0\\), are half of each: \\((1,0,0,0)=\\tfrac12\\lambda+\\tfrac12\\ell\\), the original two adding and the added two cancelling. In the first class \\(\\partial_tr=-a\\,\\mathbb J\\,r\\): each pair \\((u\_\\varepsilon,v\_\\varepsilon)\\) turns as \\(e^{-at\\mathbb J}\\), which with \\(\\mathbb J\\) written as \\(i\\) is \\(e^{-iat}\\), \\(E=+\\hbar a\\). In the second it turns the other way, \\(E=-\\hbar a\\). *Energy with its sign is the sense in which the walker turns in its own plane: the sense of its clock* (IV.5a (H)). Which class is called positive is the convention \\(\\psi\\propto e^{-iEt/\\hbar}\\) together with the orientation of the record’s cut (VI.4 (F)). The walker launched at rest, III.8c (5), is the leading class.

In motion the two classes are \\(v=\\pm J(k)u\\) with III.8c (15)’s \\(J(k)\\), and the projectors onto them are

\\\[\\mathbb P\_\\pm(k)=\\tfrac12\\Big(\\mathbb 1\\pm\\frac{\\mathbb J\\,\\mathbb G(k)}{\\omega_k}\\Big),\\qquad\\mathbb G(k)=-ck\\,\\mathbb A\\,\\mathbb J+a\\,\\mathbb T,\\\]

\\(\\mathbb P\_+\\) the leading class (\\(\\mathbb G(k)\\) is (2) on a mode, next paragraph). The launch \\(v(0)=Tu(0)\\) is the leading class exactly at rest; at wavenumber \\(k\\) it leaves the weight

\\\[\\frac{\\omega_k-a}{2\\omega_k}=\\frac{E-mc^2}{2E}\\\]

in the other class: \\(0\\), \\(0.090\\), \\(0.276\\) at \\(ck/a=0,\\,0.7,\\,2\\). The added two are the original two’s quadrature partner — the same signal a quarter period on — exactly at rest, and to order \\((c/a\\sigma)^2\\) for a packet of width \\(\\sigma\\) (III.8c (E)).

**(D) The two clocks on the four.** Put in a mode. In the \\(t\\)-reading the pair winds along \\(x\\), \\(\\partial_x\\to k\\,\\mathbb J\\); in the walk’s clock the boosted rest count slopes in \\(t\\), \\(\\partial_t\\to ck\\) (III.15d (8)):

\\\[\\partial_tr=\\mathbb G(k)\\,r,\\quad\\mathbb G(k)=-ck\\,\\mathbb A\\,\\mathbb J+a\\,\\mathbb T;\\qquad\\qquad c\\,\\partial_xr=\\mathbb M_k\\,r,\\quad\\mathbb M_k=-ck\\,\\mathbb A+a\\,\\mathbb S .\\\]

\\(\\mathbb M_k\\) is III.15d’s \\(M_k\\) on the original two and on the added two. Turn the two backward detectors into each other, \\((u\_-,v\_-)\\to(-v\_-,u\_-)\\):

\\\[\\mathcal D=\\begin{pmatrix}1&0&0&0\\\\0&0&0&-1\\\\0&0&1&0\\\\0&1&0&0\\end{pmatrix},\\qquad\\mathcal D^{-1}\\,\\mathbb G(k)\\,\\mathcal D=\\mathbb J\\,\\mathbb M_k .\\tag{4}\\\]

4. the \\(t\\)-generator is \\(\\mathbb J\\) times the clock’s generator. So in each eigen-direction of \\(\\mathbb M_k\\) the walk’s clock scales and the \\(t\\)-reading turns the pair:

\\\[\\mathbb M_k\\,r=+\\omega_k\\,r:\\ \\ e^{+\\omega_k\\tau}\\ \\text{along }x,\\ \\ e^{+\\omega_kt\\mathbb J}\\ \\text{in }t\\ \\ (E\<0);\\qquad\\mathbb M_k\\,r=-\\omega_k\\,r:\\ \\ e^{-\\omega_k\\tau},\\ \\ e^{-\\omega_kt\\mathbb J}\\ \\ (E>0).\\tag{5}\\\]

5. at rest, in the original two only: \\(\\mathcal D\\,(1,1,0,0)^{\\mathsf T}=(1,0,0,1)^{\\mathsf T}=\\ell\\) and \\(\\mathcal D\\,(1,-1,0,0)^{\\mathsf T}=(1,0,0,-1)^{\\mathsf T}=\\lambda\\). *The combination the count grows in is the lagging loop; the one that decays is the leading loop.* It holds at every \\(k\\). \\(\\mathcal D\\) is a relabelling of III.9’s kind, in the opposite sense to III.9’s own; with \\(\\mathbb J\\) written as \\(i\\) it is \\(\\begin{pmatrix}1&0\\\\0&i\\end{pmatrix}\\), and (4) is \\(iG(k)=-D\\,M_k\\,D^{-1}\\): the growth rate along the clock and the turn rate in \\(t\\) are one number, \\(\\omega_k\\).

*What went in.* No new input: III.5a’s lines written for the original two and the added two, and products of the four matrices. (B) says what the original two have and lack; (C)–(D) need the four detectors of III.8c (B), which are taken there. What V.1–V.2 show is that a walker a field can *turn* must be read by all four; one read by the original two alone can still be coupled through its coin. The names “positive” and “negative” are a convention and an orientation; the correspondence growing ↔︎ lagging, decaying ↔︎ leading is not.

(`scripts/audit_2026_11_09.py`: A for (4)–(5) at \\(k=0,\\,0.7,\\,2\\) and III.9’s direction; B for \\(J(k)\\), the local commutant and the weight \\((\\omega_k-a)/2\\omega_k\\); C for \\(\\mathbb P\_\\pm\\); D for all four against the original two alone; F for (1)–(3). `scripts/real_dictionary.py`: the rest case of (5). `scripts/master_v33_new_text_checks.py` 1–3: the two loops written out in time, the original two alone as half of each, the weight in the other class.)

### III.15f What turns the clock hand: the flip rate between the two time senses

*The picture first.* The two forward detectors of III.8c (B) read the clock hand at its two positions, and so do the two backward ones. Each pair is a hand: a length and an angle. What turns it? Not the advance, which only carries it along its light line. The flips do: each flip feeds one time sense from the other. So how fast the forward hand turns depends on how much the backward detectors hold at the same positions of the hand. The rate is the flip rate \\(a\\) times a ratio, and the ratio is the overlap of the forward and the backward detectors of the same hand position, over what the forward ones hold. At rest with the same amount in both, the ratio is one and the hand turns at \\(a\\). In motion the ratio is not one: the flip rate is deformed. With a field, the field’s own turn is added: the sum is V.4’s effective flip rate.

*Symbols.* \\(r\_\\varepsilon=(u\_\\varepsilon,v\_\\varepsilon)^{\\mathsf T}\\), \\(\\rho\_\\varepsilon=u\_\\varepsilon^2+v\_\\varepsilon^2\\), \\(\\theta\_\\varepsilon=\\operatorname{atan2}(v\_\\varepsilon,u\_\\varepsilon)\\), \\(D\_\\pm=\\partial_t\\pm c\\,\\partial_x\\), \\(W=u\_+v\_--v\_+u\_-\\) and \\(b\_\\pm=\\frac e\\hbar(\\Phi\\mp cA_x)\\), all of III.8c (C). Nothing is new.

**(A) The rate of each hand.** III.8c (10a), written out for the two time senses:

\\\[\\begin{aligned}D\_+u\_+&=a\\,u\_-+b\_+\\,v\_+,&D\_+v\_+&=a\\,v\_--b\_+\\,u\_+,\\\\D\_-u\_-&=-a\\,u\_++b\_-\\,v\_-,&D\_-v\_-&=-a\\,v\_+-b\_-\\,u\_- .\\end{aligned}\\\]

The angle of a pair turns at \\(\\rho\\,D\\theta=u\\,Dv-v\\,Du\\). For the forward hand and for the backward one:

\\\[\\rho\_+D\_+\\theta\_+=u\_+\\big(a\\,v\_--b\_+\\,u\_+\\big)-v\_+\\big(a\\,u\_-+b\_+\\,v\_+\\big)=a\\,W-b\_+\\,\\rho\_+,\\\]

\\\[\\rho\_-D\_-\\theta\_-=u\_-\\big(-a\\,v\_+-b\_-\\,u\_-\\big)-v\_-\\big(-a\\,u\_++b\_-\\,v\_-\\big)=a\\,W-b\_-\\,\\rho\_- .\\tag{1}\\\]

1. III.8c (8) and (10b). The same \\(W\\) turns both hands.

**(B) \\(W\\) is an overlap.** With the table of III.8c (B) — forward: \\(+u\_+\\) at the start, \\(+v\_+\\) a quarter on; backward: \\(-v\_-\\) at the start, \\(+u\_-\\) a quarter on:

\\\[-W=\\underbrace{u\_+\\,(-v\_-)}\_{\\text{hand at the start}}+\\underbrace{v\_+\\,u\_-}\_{\\text{hand a quarter on}}\\qquad\\big(\\text{forward}\\times\\text{backward in each}\\big).\\tag{2}\\\]

2. so (1) says: the hand of each time sense turns at

\\\[D\_\\pm\\theta\_\\pm=-\\Big(a\\;\\frac{-W}{\\rho\_\\pm}+b\_\\pm\\Big):\\tag{3}\\\]

3. the flip rate, times the overlap of the two time senses at the same hand position over the record of that time sense; plus the field’s rate. Four cases (`detectors_in_master_names.py` 5–7):

| the walker | \\(-W/\\rho\_+\\), \\(-W/\\rho\_-\\) | the hands turn at |
| --- | --- | --- |
| in the original two alone, \\(v=0\\) | \\(0\\), \\(0\\) | they do not turn: the angles are frozen at \\(0\\) or \\(\\pi\\) (III.8c (C) 10) |
| at rest, the same amount forward and backward, \\(v=Tu\\) | \\(1\\), \\(1\\) | \\(a+b\\): V.4’s effective flip rate |
| at rest, opposite amounts, \\(v=-Tu\\) | \\(-1\\), \\(-1\\) | \\(-a+b\\): the clock turns the other way, and the same field slows it — the opposite charge |
| one clock run to rapidity \\(\\eta\\) | \\(e^{-\\eta}\\), \\(e^{+\\eta}\\) | \\(a\\,e^{-\\eta}+b\_+\\) and \\(a\\,e^{+\\eta}+b\_-\\), each along its own light line |

4. in the last row \\(\\rho\_\\pm=e^{\\pm\\eta}\\) and \\(W=-1\\) (III.8c (12) with its added two), so \\(-W/\\rho\_\\pm=\\sqrt{\\rho\_\\mp/\\rho\_\\pm}\\): III.10’s \\(D\_\\pm S\_\\pm=-mc^2\\sqrt{\\rho\_\\mp/\\rho\_\\pm}\\cos\\varphi\\), read as a rate. The time sense that holds more is fed less, in proportion, by the one that holds less: that is the deformation. The two rates \\(a\\,e^{\\mp\\eta}\\) are the one clock read along the two light lines.

**(C) What an earlier version of this section asked.** It asked whether the figure closed by an outward route and “its twin” was Feynman’s loop with its \\(-1\\). There is no such figure: the added detectors are not a route. What that section found stays, where it belongs. Two routes of a measurement that make their flips in opposite senses multiply to a negative pair, \\((-\\sin\\theta)(+\\sin\\theta)\\): the two mass corners (III.15a (E)). Out from \\((S,+)\\) and home to \\((S,-)\\) sums to \\(0\\) over all far points; home to \\((S,+)\\) sums to \\(1\\) (VI.1 (A)). The statistics minus is in neither: it is exclusion blocking a Z into an occupied state (VI.1a (C)) and the exchange sign on closed histories of the configuration count (VI.5 (A)). (`scripts/route_twin_loop.py`: the exact enumeration behind these three statements; what it calls the twin is a route launched in the backward sense.)

*What went in.* Nothing: III.8c (7), (8), (10a), (10b), and the table of III.8c (B).

## III.16 The other question about the same equation: Part II is the source problem

*The picture first.* There are two different questions one can ask of a wave equation, and they have different answers even though the equation is the same.

- *The evolution question.* “The particle was here; where is it now?” You give the state on one line (the clock reading \\(x=0\\)) and the equation carries it forward. The answer is regular everywhere: nothing is singular about “being here”. That is III.14, and the answer is \\(I_0\\).
- *The source question.* “Something sits at the origin and keeps emitting; what field does it make?” The answer must be singular at the origin — the source is there — and must die away far from it. That is \\(K_0\\).

Part III asks the first question, because a process with a clock is a thing that evolves. Part II, it turns out, asks the second. The reason is II.7: Part II’s step label \\(\\tau\\) is not a physical time and nothing may depend on it. So the only physical object Part II can produce is its kernel *summed over all \\(\\tau\\)* — “the particle passed through here at some step, any step” — which is a source sitting at the origin for all \\(\\tau\\). That sum is the Green’s function.

*The computation.* Part II’s kernel after \\(\\tau\\) steps, in its own clock, is the Gaussian of II.2, \\((4\\pi D\\tau)^{-d/2}e^{-s^2/4D\\tau}\\) in \\(d\\) dimensions, times the weight that II.7’s constant \\(K\\) attaches to each unit of \\(\\tau\\); in the signature where the walk is real (II.1’s time kick continued, so that \\(s^2=x^2+c^2t_E^2\\)) and in units \\(\\hbar=c=1\\), with the clock rescaled to \\(D=1\\), \\(\\tau=D\\tau\_{\\text{p}}=\\tau\_{\\text{p}}/2m\\) (Schwinger’s proper-time parameter), II.7’s weight \\(e^{-mc^2\\tau\_{\\text{p}}/2\\hbar}\\) per unit of proper time is \\(e^{-m^2\\tau}\\). This is the proper-time representation of the propagator (Feynman, *Phys. Rev.* 80, 440 (1950); Schwinger, *Phys. Rev.* 82, 664 (1951)). Summing over \\(\\tau\\) (the Lorentzian statement is the continuation of III.15):

\\\[\\int_0^\\infty\\frac{d\\tau}{4\\pi\\tau}\\,e^{-s^2/4\\tau}\\,e^{-m^2\\tau}=\\frac1{2\\pi}K_0(ms)\\qquad(1{+}1),\\qquad\\int_0^\\infty\\frac{d\\tau}{(4\\pi\\tau)^2}\\,e^{-s^2/4\\tau}\\,e^{-m^2\\tau}=\\frac{m}{4\\pi^2s}K_1(ms)\\qquad(3{+}1),\\\]

both checked by quadrature (`scripts/brownian_killed_4d.py` 1). The \\(3{+}1\\) line is dimension-blind in its derivation: Part II never used the number of space dimensions.

*The same integral read as a process.* \\(e^{-m^2\\tau}\\) is a survival probability and \\((4\\pi\\tau)^{-d/2}e^{-s^2/4\\tau}\\) is where a Brownian particle is after \\(\\tau\\); the integral over \\(\\tau\\) is therefore the *expected time spent* at \\(x\\) by a Brownian particle that is killed at the rate \\(m^2\\) — its occupation density. Nothing in it is negative or complex. Run as a Monte Carlo in four Euclidean dimensions (200 000 paths, time in each radial shell divided by the shell’s volume), the occupation density matches \\(\\frac{m}{4\\pi^2s}K_1(ms)\\), averaged over the same shell, to \\(0.6\\%\\) in every shell from \\(s=0.3/m\\) to \\(3/m\\), within \\(1.7\\) standard errors (`scripts/brownian_killed_4d.py` 2). So the \\(3{+}1\\) source answer *is* produced by a positive process — Part II’s own walk, with the Lorentzian sign paid once in II.1’s \\(\\eta^{00}\\) by working in the signature where the walk is real. The sign is in the signature, not in the weights; that is why there is no sign problem here and why there is one in IV.2a.

*What this settles.* The same \\(K_0\\) is reached by a second road: the on-shell plane waves of III.15 item 2, written with the frame variable \\(\\eta\\) (\\(p=mc\\sinh\\eta\\), \\(dp/E=d\\eta/c\\)) and continued to imaginary time, \\(\\int d\\eta\\,e^{-\\mu x\\cosh\\eta}=2K_0(\\mu x)\\) — every frame weighted equally. Where does “every frame equally” come from? Not from III.1: III.1 has a clock, hence a frame, and its loss factor \\(e^{-ax/c}\\) is not invariant. It comes from the same place as \\(\\eta^{\\mu\\nu}\\) in II.1 — Lorentz invariance taken as an input. Part II pays for that input with an imaginary time-kick and gets the source answer; Part III pays with a chosen frame and gets the evolution answer; the frame sum is the exchange rate between the two payments. In the ledger: *Lorentz invariance, as in II.1*.

### What went in (Part III)

| the clock is \\(x\\) and the coin is in \\(t\\) (the time-sense particle S) | the choice of process; the telegraph T is the other choice, kept as remarks |
| --- | --- |
| master equations, conserved charge \\(P\_+-P\_-\\), pair creation, telegraph-in-\\(x\\), Kac’s limit in \\(x\\) | derived from III.1 |
| clock reversal: no osmotic velocity; information in \\(P\_-/P\_+\\) | derived |
| \\(\\langle(\\Delta t/\\Delta x)^2\\rangle=1/c^2\\): no \\(\\rho\\)-only term in \\(H\\) | derived |
| direction matrix and turn matrix anticommute; turn matrix squares to \\(-1\\) ⇒ Klein–Gordon with the right sign, real rate | derived (III.6, III.8) |
| Madelung pairs + Hamiltonian flow | the walker is read by four detectors, III.1’s two and the two added (III.8c (B)); III.10’s four equations then follow for any two columns (III.8c (C)), and the Hamiltonian and bracket of III.9–III.10 are a packaging of them, a reading. It imposes no condition on \\(a\\). What it posits differs from Part I: there \\(H\\) came out of the process; here the posit was that the signed \\(t\\)-reading is the physical one and that an event is its square; III.15b and III.15c reduce both to II.1 — a record is data on a spacelike cut, and on that cut a reading goes as the square root of a null track coordinate while a count goes as the coordinate — leaving one reading: the weight of an event is its frequency |
| \\(a=mc^2/\\hbar\\) | the rate named by the mass it produces |
| the complex pairing: \\(i\\) = a quarter-turn already present; commuting matrices stay linear, anticommuting ones go antilinear; III.9’s \\(i\\) needs the added detectors because \\(A\\) anticommutes with \\(T\\), and it is the walker’s own turn (III.8c (B)) | derived (III.8b), from the matrices on the page |
| the measurement: every outward route with every return route that agrees on the endpoint; the completed detectors: time sense and position of the clock hand read separately, four instead of two; the walker at rest, \\(v=Tu\\); III.10’s fluid for any two columns; one sign of energy to \\(O((c/a\\sigma)^2)\\) for a source of width \\(\\sigma\\) | III.8c: derived, the added detectors taken (input 11); the lines do not force them |
| the current is the phase slope weighted by the mass term, exactly; III.11’s box is its slow limit | derived (III.11a, Gordon’s decomposition) |
| Dirac; \\(H\\); Nelson’s \\(v,u\\), \\(D=\\hbar/2m\\), the osmotic term, Schrödinger | derived |
| \\(c,\\hbar,m\\) | constants |
| what the forward-time observer may call a probability at fixed \\(t\\) | open: the charge \\(P\_+-P\_-\\) is conserved and indefinite, the count \\(P\_++P\_-\\) is positive and not conserved (III.3) |
| the exact kernel: \\(\\delta\\) on the light line with mass \\(e^{-ax/c}\\), \\(I_0\\), \\(I_1\\) bulk outside the cone | derived from III.2 (III.14), Monte Carlo checked |
| a particle is a no-flip segment; the amplitude is the sum over all zigzags \\(=\\) the on-shell superposition; classical kernel and amplitude are one function across the light cone | derived (III.15), three computations checked against each other |
| Part II with \\(\\tau\\) summed is the source answer, \\(K_0\\); the frame sum \\(\\int d\\eta\\) reaches the same \\(K_0\\) | derived (III.16); the equal weight of frames is Lorentz invariance as an input, as \\(\\eta^{\\mu\\nu}\\) in II.1 |
| the count of returning zigzags along the walk’s clock, weighted by \\(e^{-ck\\Delta t}\\), grows at \\(\\omega_k\\); the weight is III.15c’s boost of the plain count: \\(\\hbar k=mc\\sinh\\eta\\), \\(\\hbar\\omega_k=mc^2\\cosh\\eta\\), \\(Q\_-/Q\_+=e^\\eta\\); near rest the rate is \\(\\hbar k^2/2m\\), Part I’s \\(D\\) | derived (III.15d), from III.14’s kernel and III.15c; “energy” for the growth rate is III.8c (D)’s definition |
| a walker in the original two detectors has a frequency and no sign of it, and is half a walker; read by all four it has two classes, its clock turning one way (\\(E>0\\)) or the other (\\(E\<0\\)); the combination the walk’s count grows in is the lagging one, \\(\\mathcal D^{-1}\\mathbb G(k)\\mathcal D=\\mathbb J\\mathbb M_k\\) | derived (III.15e), given the four detectors of III.8c (B); the names of the classes are a convention and the cut’s orientation |
| the clock hand of each time sense turns at the flip rate times the overlap of the forward and the backward detectors of the same hand position, over that sense’s record, plus the field’s rate: \\(D\_\\pm\\theta\_\\pm=-(a\\,(-W)/\\rho\_\\pm+b\_\\pm)\\) | derived (III.15f), from III.8c (8), (10b) |

---

# Part IV — Three dimensions: the vortex, its axis, and the two time senses

## IV.0 The physical objects

*The picture first.* In one dimension the walker had a time sense and nothing else. In three dimensions it has one more thing, an **axis**: the direction it moves along between flips. Everything spin-like in this Part is that axis and what the walk does with it. Three objects, each derived later in the Part:

- the **axis**, a unit vector \\(\\hat{\\mathbf s}(\\mathbf x,t)\\). It is not a label the walker can carry (IV.2); it is held as a two-component column, the half-angle object of III.8a (IV.3). It is at once the direction the walker streams along and the direction its angular momentum points, \\(\\tfrac\\hbar2\\) per particle (IV.6, IV.7a).
- the **circulation about the axis**. A cloud of walkers with an axis is not still: its density gradient drives a flow around the axis, \\(\\mathbf u\\times\\hat{\\mathbf s}\\), the osmotic velocity turned a quarter turn (IV.8 (C)). It carries no density, it is regular where the cloud is fullest, and for a walker at rest in the original four detectors alone it is the walker itself going round at \\(c\\) with radius \\(\\hbar/2mc\\) (IV.8 (D), `true_spin.py` 3b). This is the vortex, with no medium: the circulating thing is the walker’s own cloud.
- the **sense** of that circulation. Read forward in time it goes one way; read backward, the other, because the circulation reverses with time while the osmotic velocity does not (IV.8 (C)). The two senses are the two time senses of Part III (IV.4): a left-handed particle is a right-handed one with its time run the other way. One axis, two senses, one per time direction.

The flatlander of III.8a is the plane perpendicular to the axis. What he sees of a rotation of the axis is III.8a; the thing one plane-axis rotation hits is the column; and what his matrices become once the column is written in \\(\\rho\\), \\(S\\) and \\(\\hat{\\mathbf s}\\) is IV.7a — the axis itself.

This Part shows that the time-sense construction of Part III, carried to three dimensions with nothing added but “one speed \\(c\\) in every direction”, produces exactly these objects and nothing else: the column that *is* the axis as the flatlander holds it (IV.3); the doubling into the two senses, which the mass swaps at the Compton rate (IV.4–IV.5); the spin density \\(\\frac\\hbar2\\rho\\,\\hat{\\mathbf s}\\) as a canonical momentum, with the rotation algebra (IV.6); the momentum density carrying \\(\\frac\\hbar2\\nabla\\times(\\rho\\hat{\\mathbf s})\\) (IV.7); and the fluid of the axis — its bending energy, its circulation, its slow limit with \\(g=2\\), and the vortex worked out — in IV.8. What is *not* here, and was in earlier versions of this picture: a medium circulating around a vortex line, a circulation \\(\\hbar/m\\) attached to that line, and kicks from the medium that tilt the axis. The vortex with a hole at its centre is the orbital one, a phase winding; the spin needs no hole (IV.8 (E)).

## IV.1 Setting up

The rule of III.1 with an axis: the particle moves at \\(c\\) along its axis \\(\\hat{\\mathbf s}\\), and its time sense flips,

\\\[\\frac{d\\mathbf x}{d\\tau}=c\\,\\hat{\\mathbf s},\\qquad\\frac{dt}{d\\tau}=\\varepsilon=\\pm1,\\qquad\\varepsilon\\to-\\varepsilon\\ \\text{at rate }a .\\\]

To the forward-time observer the \\(\\varepsilon=-1\\) particle moves along \\(-\\hat{\\mathbf s}\\). What “the axis \\(\\hat{\\mathbf s}\\)” is as a state is not assumed; IV.2 shows it cannot be a classical label and IV.3 shows what it is. Keep the structure of III.8: a column of densities, one per possible state of motion, one direction matrix per axis, one turn matrix \\(T\\) for the reversal of the time sense, and the first-order equation

\\\[\\partial_t(\\text{column})=-c\\big(A_x\\partial_x+A_y\\partial_y+A_z\\partial_z\\big)(\\text{column})+a\\,T\\,(\\text{column}).\\\]

How many components and which matrices is not assumed. The physical demand, which held automatically in III.6: the particle moves at \\(c\\) with no preferred direction, i.e. squaring the operator gives \\(c^2(\\partial_x^2+\\partial_y^2+\\partial_z^2)\\) on every component plus a constant. Squaring,

\\\[(A_x\\partial_x+A_y\\partial_y+A_z\\partial_z)^2=A_x^2\\partial_x^2+A_y^2\\partial_y^2+A_z^2\\partial_z^2+(A_xA_y+A_yA_x)\\partial_x\\partial_y+(A_yA_z+A_zA_y)\\partial_y\\partial_z+(A_zA_x+A_xA_z)\\partial_z\\partial_x,\\\]

so we need \\(A_k^2=1\\), \\(A_kA_l+A_lA_k=0\\) for \\(k\\ne l\\), and likewise \\(A_kT+TA_k=0\\) and, for the constant to be \\(-a^2\\) as in III.6, \\(T^2=-1\\). In the words of III.8a: three mutually perpendicular plane-axis half-turns for the axes, and a quarter-turn for the time sense that anticommutes with all three. The question is how many components a column needs to be hit by all four.

## IV.2 The components cannot be labelled by “which way I am going”

Suppose they were — a sign for each axis — so that every \\(A_k\\) is diagonal with entries \\(\\pm1\\). Diagonal matrices commute, so \\(A_xA_y+A_yA_x=2A_xA_y\\) is diagonal with entries \\(\\pm2\\), never zero: the square has \\(\\partial_x\\partial_y\\) cross terms, and the particle does not move at \\(c\\) in every direction. **In three dimensions the components cannot be labelled by a direction of motion.** The direction matrices for different axes must fail to commute. This is the one place where the three-dimensional construction cannot be a classical process of the type III.1 — and it is also the first sign that what the column holds is not a direction of motion but an *axis*, the thing of IV.0.

## IV.2a Why the direction cannot be a classical label: the mechanism

*The picture first.* In one dimension a particle moving at \\(c\\) has two possible directions, and III.8 found that the free step is the pair of projectors \\(\\frac12(1\\pm\\sigma_z)\\): a right-mover keeps going right, a left-mover keeps going left, nothing else happens between flips. Try to do the same in three dimensions. The particle should move at \\(c\\) in some direction \\(\\hat n\\), and the free step should be “move along \\(\\hat n\\) with some weight \\(W(\\hat n)\\) attached to the column”. Here is the mechanism that breaks this, in words.

Take a column whose axis is \\(\\hat m\\) (IV.3: the column *is* its axis). It should move, on average, along \\(\\hat m\\) at speed \\(c\\). But a column with axis \\(\\hat m\\) is not orthogonal to a column with axis \\(\\hat n\\) for \\(\\hat n\\ne\\hat m\\): the overlap is \\(\\cos^2(\\theta/2)\\). So a move “along \\(\\hat n\\) with the projector onto axis \\(\\hat n\\)” picks up part of the \\(\\hat m\\) column for every \\(\\hat n\\) — the \\(\\hat m\\) particle is sent off in all directions at once, with weights that are largest near \\(\\hat m\\) but nowhere zero. Average the displacement over those weights and you get only one third of \\(c\\,\\hat m\\) (on a sphere, \\(\\langle n_in_j\\rangle=\\frac13\\delta\_{ij}\\)). To recover the full \\(c\\,\\hat m\\) you must triple the weights, and then they no longer add up to one unless you subtract something — and what you subtract is negative weight. There is no way to spread a column over all directions, move each piece at \\(c\\), and still have the whole thing advance at \\(c\\), with every piece counted positively. In one dimension the “sphere” is two points, the overlap between the two axes is zero, and the average is the full \\(c\\): nothing to fix.

*The computation.* A step of the massless column is \\(\\psi(\\mathbf x,t+\\epsilon)=\\int\\frac{d\\Omega}{4\\pi}\\,W(\\hat n)\\,\\psi(\\mathbf x-\\epsilon c\\hat n,t)\\), and it has to reproduce \\(\\partial_t\\psi=-c\\,\\boldsymbol\\sigma\\cdot\\nabla\\psi\\) (IV.3). Expanding to first order in \\(\\epsilon\\): \\(\\int\\frac{d\\Omega}{4\\pi}W=1\\) and \\(\\int\\frac{d\\Omega}{4\\pi}W(\\hat n)\\,\\hat n_j=\\sigma_j\\). Write \\(W=w_0(\\hat n)+\\mathbf w(\\hat n)\\cdot\\boldsymbol\\sigma\\). The second condition, component by component, is \\(\\int\\frac{d\\Omega}{4\\pi}w_k(\\hat n)\\,\\hat n_j=\\delta\_{kj}\\); set \\(k=j\\) and sum over the three directions: \\(\\int\\frac{d\\Omega}{4\\pi}\\mathbf w\\cdot\\hat n=3\\). For \\(W\\) to be a positive weight on columns, \\(w_0\\ge|\\mathbf w|\\). Then

\\\[3=\\int\\frac{d\\Omega}{4\\pi}\\,\\mathbf w\\cdot\\hat n\\;\\le\\;\\int\\frac{d\\Omega}{4\\pi}\\,|\\mathbf w|\\;\\le\\;\\int\\frac{d\\Omega}{4\\pi}\\,w_0=1 .\\\]

Impossible. The isotropic solution is unique, \\(W=1+3\\,\\boldsymbol\\sigma\\cdot\\hat n\\), and its eigenvalues are \\(4\\) and \\(-2\\). In one dimension the same chain reads \\(1\\le1\\): saturated, by the projectors.

*Where the minus sign sits.* It is not a negative probability of existing. Split the weight by helicity, \\(P(\\hat n)=\\frac12(1+\\boldsymbol\\sigma\\cdot\\hat n)\\):

\\\[1+3\\,\\boldsymbol\\sigma\\cdot\\hat n=4\\,P(\\hat n)-2\\,P(-\\hat n):\\\]

weight \\(4\\) on the column whose spin points *along* the motion, weight \\(-2\\) on the column whose spin points *against* it. For the right-handed column of IV.3, spin-against-motion is the negative-energy column. III.9a keeps that apart from the backward time sense of a segment: the link between the two is Feynman’s rule (III.15, IV.10), which carries the negative-energy half backward in \\(t\\), and under that rule the sign is the sign of \\(q=P\_+-P\_-\\) in III.3, the forward-time observer’s signed count. The same sign in another dress: multiply the projectors along a path and the product is a positive magnitude times \\(e^{-i\\Omega/2}\\), \\(\\Omega\\) the solid angle swept by the axis — the circulation of the axis, with the sign of its handedness. This is a holonomy, but on the sphere of directions, not in spacetime: in flat space a translation turns nothing (\\(\\partial\_\\mu\\) commute), while a change in the direction of motion turns the column by twice the angle between the two steps (III.8a). Carried around a closed loop of directions, the column comes back with the phase \\(\\Omega/2\\) of the solid angle enclosed (`scripts/direction_weights_and_holonomy.py` 4: 200 random spherical triangles, to \\(10^{-14}\\); 1–3 there for the average over directions, the weight \\(1+3\\,\\boldsymbol\\sigma\\cdot\\hat n\\) and the inequality).

*What is different from one dimension.* In III.1 the massless particle never flips; flips come with the mass, at the finite rate \\(a\\). Here the massless step already needs the backward sense with weight \\(-2\\) against \\(4\\) at *every* step — a flip rate that is infinite per unit time. The forward-time observer’s \\(q\\) is then a small difference of two large positive counts. That is the reason the direction cannot be a classical label in three dimensions, and why Part IV carries it as the column instead. Read physically: the column’s velocity components are each \\(\\pm c\\) at all times while the axis advances at \\(c\\) on average — circulation at \\(c\\) about an axis that drifts. The vortex of IV.0 is not a picture added to the equation; it is what the equation does with the freedom the inequality above forbids.

*The same thing seen on the light cone.* In three dimensions a massless kernel from a point is not \\(\\delta(ct-r)\\) on the expanding sphere; for data in \\(\\psi\\) it is \\(\\partial_t\[\\delta(ct-r)/4\\pi r\]\\), which contains \\(\\delta'(ct-r)\\): a positive skin just outside the shell and a negative skin just inside. In III.3’s language that is a forward-time shell with a backward-time shell trailing it by nothing — a pair on the cone, a dipole layer of \\(q\\). The two positive shells are each infinite; only their difference is finite. Same mechanism, same place.

*A test that was run.* The obvious attempt — run III.1 with the radius \\(r\\) as the clock and divide the result by the spherical Jacobian \\(4\\pi r^2\\) — was run against \\(\\frac{m}{4\\pi^2s}K_1(ms)\\) (`scripts/radial_test.py`). It fails, and it fails before any formula: a Green’s function depends on \\(s=\\sqrt{r^2-c^2t^2}\\) alone, and the radial kernel over \\(4\\pi r^2\\) does not — along a fixed \\(r\\) it grows with \\(s\\) while the Green’s function falls. The three reasons are the ones above: a one-clock process answers the evolution question, not the source question; \\(u/r\\) is the s-wave map and \\(u/r^2\\) is nothing; and the loss factor \\(e^{-ar/c}\\) knows the frame. The source answer in three dimensions comes from Part II (III.16), not from a radial clock.

## IV.3 Massless: two components, and the column is the axis as the flatlander holds it

The \\(2\\times2\\) solution of IV.1 (unique up to a change of basis and one choice of orientation: flipping the sign of one axis changes \\(A_xA_yA_z=\\pm i\\), which no similarity undoes; the document takes \\(+i\\)):

\\\[A_x=\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix},\\qquad A_y=\\begin{pmatrix}0&-i\\\\i&0\\end{pmatrix},\\qquad A_z=\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}.\\\]

\\(A_y^2=\\begin{pmatrix}(-i)(i)&0\\\\0&(i)(-i)\\end{pmatrix}=1\\); \\(A_xA_y+A_yA_x=\\begin{pmatrix}i&0\\\\0&-i\\end{pmatrix}+\\begin{pmatrix}-i&0\\\\0&i\\end{pmatrix}=0\\); the others likewise. The operator, with \\(\\hat p_k=-i\\hbar\\partial_k\\),

\\\[c\\big(A_x\\hat p_x+A_y\\hat p_y+A_z\\hat p_z\\big)=c\\begin{pmatrix}\\hat p_z&\\hat p_x-i\\hat p_y\\\\\\hat p_x+i\\hat p_y&-\\hat p_z\\end{pmatrix},\\\]

squares to \\(c^2(\\hat p_x^2+\\hat p_y^2+\\hat p_z^2)\\) times the identity: top-left \\(\\hat p_z^2+(\\hat p_x-i\\hat p_y)(\\hat p_x+i\\hat p_y)\\), off-diagonal \\(\\hat p_z(\\hat p_x-i\\hat p_y)-(\\hat p_x-i\\hat p_y)\\hat p_z=0\\). Row by row,

\\\[i\\hbar\\,\\partial_t\\psi_1=-i\\hbar c\\big\[\\partial_z\\psi_1+(\\partial_x-i\\partial_y)\\psi_2\\big\],\\qquad i\\hbar\\,\\partial_t\\psi_2=-i\\hbar c\\big\[(\\partial_x+i\\partial_y)\\psi_1-\\partial_z\\psi_2\\big\].\\\]

Component 1 streams along \\(+z\\) and component 2 along \\(-z\\), as in Part III. Motion along \\(x\\) or \\(y\\) is not a property of a component; it is a *transition between the two*, with the factor \\(\\partial_x\\mp i\\partial_y=e^{\\mp i\\phi}(\\partial_r\\mp\\frac ir\\partial\_\\phi)\\) in polar coordinates of the \\(xy\\)-plane: a sideways step carries one unit of rotation about \\(z\\). This is the \\(i\\) of III.8b, met as the only way to move sideways at speed \\(c\\) with two components. (The massless two-component Dirac equation, the Weyl equation.)

*These are the three plane-axis rotations of space, and the column is the axis.* For a vector \\(\\mathbf n=(n_x,n_y,n_z)\\) put \\(M(\\mathbf n)=n_xA_x+n_yA_y+n_zA_z=\\begin{pmatrix}n_z&n_x-in_y\\\\n_x+in_y&-n_z\\end{pmatrix}\\), the three-component version of III.8a. Multiplying out, \\(M(\\mathbf n)^2=|\\mathbf n|^2\\,1\\) and \\(M(\\mathbf n)M(\\mathbf n')+M(\\mathbf n')M(\\mathbf n)=2(\\mathbf n\\cdot\\mathbf n')\\,1\\) — the same two identities, and they *are* the conditions of IV.1. And as in III.8a a vector is its own half-turn: multiplying out, \\(M(\\mathbf n)\\,M(\\mathbf v)\\,M(\\mathbf n)=M\\big(2(\\mathbf n\\cdot\\mathbf v)\\mathbf n-\\mathbf v\\big)\\) for unit \\(\\mathbf n\\), the half-turn of \\(\\mathbf v\\) about \\(\\mathbf n\\). Here no projection is needed: \\(M(\\mathbf n)\\) hits \\(M(\\mathbf v)\\) from both sides and turns it; it hits a column once. The third one, \\(A_y\\), is the one the flatlander of III.8a could not write with real numbers: the half-turn about an axis out of his plane.

For the vortex axis \\(\\hat{\\mathbf s}=(\\sin\\theta\\cos\\phi,\\sin\\theta\\sin\\phi,\\cos\\theta)\\) the column its own half-turn keeps is

\\\[M(\\hat{\\mathbf s})\\begin{pmatrix}\\cos\\frac\\theta2\\\\ e^{i\\phi}\\sin\\frac\\theta2\\end{pmatrix}=\\begin{pmatrix}\\cos\\theta\\cos\\frac\\theta2+\\sin\\theta e^{-i\\phi}\\,e^{i\\phi}\\sin\\frac\\theta2\\\\ \\sin\\theta e^{i\\phi}\\cos\\frac\\theta2-\\cos\\theta\\,e^{i\\phi}\\sin\\frac\\theta2\\end{pmatrix}=\\begin{pmatrix}\\cos(\\theta-\\frac\\theta2)\\\\ e^{i\\phi}\\sin(\\theta-\\frac\\theta2)\\end{pmatrix}=\\begin{pmatrix}\\cos\\frac\\theta2\\\\ e^{i\\phi}\\sin\\frac\\theta2\\end{pmatrix}\\equiv\\chi(\\hat{\\mathbf s}).\\\]

Half angles, as in III.8a. **This column is the vortex axis as the flatlander holds it**: the spin doublet of IV.6–IV.8, with its half angles now explained, and with \\(\\hat{\\mathbf s}=\\chi^\\dagger\\boldsymbol\\sigma\\chi\\) recovering the axis from the column. Conversely \\(\\chi\\) and \\(-\\chi\\) give the same axis — the sign that a full turn of the axis puts on the column, and the half-quantum of circulation of IV.8 (B).

The massless operator is \\(c\\,M(\\mathbf p)=c|\\mathbf p|\\,M(\\hat{\\mathbf p})\\), so the positive-energy column of momentum \\(\\mathbf p\\) is \\(\\chi(\\hat{\\mathbf p})\\): the axis lies along the motion. In the vortex picture: a massless particle’s circulation plane is perpendicular to its velocity, and “the direction in which you turn is always the same” — helicity.

## IV.4 Mass needs four components: the two senses

A fourth matrix anticommuting with \\(A_x,A_y,A_z\\) does not exist in \\(2\\times2\\), whatever its square. The general \\(2\\times2\\) matrix is \\(F=\\begin{pmatrix}f_0+f_3&f_1-if_2\\\\f_1+if_2&f_0-f_3\\end{pmatrix}\\); \\(A_zF+FA_z=\\begin{pmatrix}2(f_0+f_3)&0\\\\0&2(f_0-f_3)\\end{pmatrix}\\) forces \\(f_0=f_3=0\\), then \\(A_xF+FA_x=2f_1\\cdot1\\) forces \\(f_1=0\\) and \\(A_yF+FA_y\\) forces \\(f_2=0\\). So \\(F=0\\). In the words of IV.3: a \\(2\\times2\\) matrix is \\(f_0\\,1+M(\\mathbf f)\\); the identity part commutes with everything, and the vector part is the half-turn about a vector in a space that has only three perpendicular directions, all taken. A fourth perpendicular half-turn needs a bigger column. Four components are the minimum (the number must be even, since \\(TA_kT^{-1}=-A_k\\) makes \\(A_k\\) traceless), and the solution is two blocks, each IV.3, with the turn matrix acting between the blocks (a block is half of the column of four counters; it is not the added four detectors of IV.5a):

\\\[A_x=\\begin{pmatrix}0&1&0&0\\\\1&0&0&0\\\\0&0&0&-1\\\\0&0&-1&0\\end{pmatrix},\\quad A_y=\\begin{pmatrix}0&-i&0&0\\\\i&0&0&0\\\\0&0&0&i\\\\0&0&-i&0\\end{pmatrix},\\quad A_z=\\begin{pmatrix}1&0&0&0\\\\0&-1&0&0\\\\0&0&-1&0\\\\0&0&0&1\\end{pmatrix},\\quad T=\\begin{pmatrix}0&0&1&0\\\\0&0&0&1\\\\-1&0&0&0\\\\0&-1&0&0\\end{pmatrix}.\\\]

The three \\(A_k\\) square to \\(1\\), \\(T\\) squares to \\(-1\\), and all six pairs anticommute (checked). \\(T\\) is III.8’s turn matrix acting on the pair of blocks: the quarter-turn of the time-sense plane. Components 1,2 are one block of IV.3 and components 3,4 are the other block *with the opposite sign of the direction matrices*: the second block streams the opposite way along its axis. So the column is

\\\[\\Psi=\\begin{pmatrix}\\psi_R\\\\ \\psi_L\\end{pmatrix},\\qquad\\psi_R,\\psi_L\\ \\text{each a two-component column of the IV.3 kind, with its own axis }\\hat{\\mathbf s}\_R,\\hat{\\mathbf s}\_L .\\\]

**Reading, from the rule of IV.1.** The two blocks are the two time senses \\(\\varepsilon=\\pm1\\). The lower block carries the direction matrices with the opposite sign because \\(dt/d\\tau=-1\\) with \\(d\\mathbf x/d\\tau=c\\hat{\\mathbf s}\\) unchanged: the forward-time observer sees it move along \\(-\\hat{\\mathbf s}\\). Nothing was identified; the rule says it. \\(T\\) is the reversal of \\(\\varepsilon\\), a quarter-turn between the blocks exactly as in III.8. Now the handedness. For the \\(\\varepsilon=+1\\) block, IV.3: the positive-energy column of momentum \\(\\mathbf p\\) is \\(\\chi(\\hat{\\mathbf p})\\), axis along the motion — right-handed. For the \\(\\varepsilon=-1\\) block the observer sees motion along \\(-\\hat{\\mathbf s}\\), so the axis is against the motion — left-handed. **A left-handed particle is a right-handed one with its time sense reversed** — the forward-time observer’s reading of a backward segment, velocity reversed and axis kept (Feynman–Stueckelberg), which is not Wigner’s time reversal, under which the spin reverses too (IV.8 (C) uses that one). That is “two handednesses from the two time directions”, here as a consequence of the parametrization, and \\(T\\) — the reversal of \\(\\varepsilon\\) — is the operation that turns one into the other. III.9a adds the second layer: at rest the eigencombinations of the mass term (in the standard basis, \\(\\psi_R\\pm\\psi_L\\)) are the positive- and negative-energy states, so a particle of definite energy sign at rest is the equal mixture of the two time senses, in phase or out of phase.

### IV.4a The extension is forced: perpendicular mirrors, and why four

*The picture first.* In one dimension the process handed us two matrices and said what each one does: the direction matrix squares to \\(+1\\) — applied twice it is nothing, the mark of a reflection — and the time-sense turn squares to \\(-1\\). Read them as mirrors: one mirror for the one spatial direction, and a turn for the time sense that is a mirror of the other kind. The question of this Part is what the three-dimensional version of that pair has to be, and the answer is that there is no freedom: it is fixed by the two things IV.1 asks for, and nothing else.

*One speed \\(c\\) in every direction forces perpendicular mirrors.* Move along a unit direction \\(\\hat{\\mathbf n}\\); the direction matrix for that motion is \\(\\hat{\\mathbf n}\\cdot\\mathbf A\\), and “speed \\(c\\) along \\(\\hat{\\mathbf n}\\), no leakage sideways” is \\((\\hat{\\mathbf n}\\cdot\\mathbf A)^2=1\\) for every \\(\\hat{\\mathbf n}\\). Expand the square:

\\\[(\\hat{\\mathbf n}\\cdot\\mathbf A)^2=\\sum_j n_j^2A_j^2+\\sum\_{j\<k}n_jn_k\\,(A_jA_k+A_kA_j).\\\]

For this to equal \\(\\hat{\\mathbf n}\\cdot\\hat{\\mathbf n}=1\\) for *all* \\(\\hat{\\mathbf n}\\), each coefficient must match: \\(A_j^2=1\\) and \\(A_jA_k+A_kA_j=0\\) for \\(j\\neq k\\). Three mirrors, each squaring to one, mutually perpendicular — perpendicular is what anticommuting mirrors are (III.8a: two half-turns about perpendicular axes anticommute). No other choice reproduces isotropy.

*Klein–Gordon with the right sign in every direction forces the turn to be perpendicular to all three.* Repeat III.6 along \\(\\hat{\\mathbf n}\\): the square of \\((\\hat{\\mathbf n}\\cdot\\mathbf A)\\,\\partial+a\\,T\\) has the cross term \\(a\\,\\partial\\,(\\hat{\\mathbf n}\\cdot\\mathbf A\\,T+T\\,\\hat{\\mathbf n}\\cdot\\mathbf A)\\), which must vanish for every \\(\\hat{\\mathbf n}\\) — so \\(T\\) anticommutes with each \\(A_j\\) — and the constant is \\(a^2T^2\\), which must be \\(-a^2\\): \\(T^2=-1\\), as in III.6. A fourth mirror, perpendicular to the three, of the other kind.

*How many components.* Four mutually perpendicular mirrors cannot live in an odd number of dimensions (from \\(AT=-TA\\) with both invertible, \\(\\det A\\det T=(-1)^d\\det T\\det A\\), so \\(d\\) is even) and cannot live in two: the three Pauli matrices already fill the \\(2\\times2\\) anticommuting set, and nothing anticommutes with all three (the only solution of \\(X\\sigma_j+\\sigma_jX=0\\), \\(j=x,y,z\\), is \\(X=0\\)). The first dimension that admits four is \\(4\\). So the column has four components — two time senses times the two components of an axis, as IV.4 counted it — and the count is a theorem, not a convention.

*And there is only one answer.* Two explicit solutions are on the page (the basis of IV.4 and the standard one, \\(\\gamma^0=F\\), \\(\\gamma^k=FA_k\\) of IV.9); they are related by an invertible change of basis \\(S\\), found by solving \\(S\\gamma^\\mu=\\gamma'^{\\,\\mu}S\\) (one-dimensional solution space, checked). That any two sets of four \\(4\\times4\\) matrices obeying these relations are so related is Pauli’s theorem (1936), cited here, not re-proved. Everything in the rest of Part IV is therefore independent of which basis is written down. (`scripts/gamma_uniqueness.py`: the expansion in \\(\\hat{\\mathbf n}\\) symbolically, the determinant argument, the \\(2\\times2\\) obstruction as a rank computation, the \\(4\\times4\\) solutions (IV.4’s set as printed, IV.4b’s real set and the standard one), the right-sign Klein–Gordon along a random direction, and \\(S\\).)

*Said with mirrors.* A rotation of space is two of the three spatial mirrors in succession, by twice the angle between them; a boost is a spatial mirror followed by the time mirror. The column is what one mirror hits once (III.8a). So “the matrices of Part III extend to three dimensions” means exactly: add two more mirrors, perpendicular to the two you have, keep their kinds; the smallest column that all four can hit has four entries, and there is one such column up to how you name its entries.

### IV.4b The same four matrices from the process’s clock: the four lines, line by line

*The picture first.* III.5a began with the two master equations of the one-dimensional walk, written in the walk’s own clock \\(x\\), and reached the Dirac form by equalities only: add the loss back, take it out with one exponential, name the column, solve for \\(\\partial_t\\), stack, split. The same chain runs in three dimensions, with one axis as the clock, and it shows where the two new matrices of IV.4 sit before anyone reads the system in \\(t\\). Take \\(z\\) as the clock. There are four lines instead of two. Along \\(z\\) they are III.2 twice over: two time-sense walks side by side, each with its advance, its loss and its flip. Across \\(z\\) there are two more matrices, one for \\(\\partial_x\\) and one for \\(\\partial_y\\), and both are antisymmetric: turns, in the clock itself. In one dimension the turn appeared only when the \\(x\\)-clock equation was read in \\(t\\) (III.5a (J)); in three dimensions moving sideways is already a turn in the process’s own clock, and reading in \\(t\\) swaps the roles, flip for turn and turns for advances, by the same single multiplication as in (J).

*Where the four lines come from.* Only part of (A) is read off a one-step probability. With \\(\\partial_x=\\partial_y=0\\) the lines are III.2 along \\(z\\), twice; that part is a walk. The \\(\\partial_x\\) and \\(\\partial_y\\) terms are not, and IV.2a says they cannot be: a weight on columns that moves each piece at \\(c\\) and advances the whole at \\(c\\) has negative entries. They are what IV.4a’s two demands fix once the three single-axis walks share one column: speed \\(c\\) in every direction and Klein–Gordon with the right sign, and IV.4c gives them their physical source. That this is so is checked at the end, in (K), and nowhere before: the chain below uses the four lines only as a system of equations, and reads off what the matrices are at the last step, not the first.

**(A) The system.** Four real functions \\(P_1,P_2,P_3,P_4\\) of \\((t,x,y,z)\\), and \\(z\\) the clock:

\\\[c\\,\\partial_zP_1=-\\partial_tP_1-aP_1+aP_3\\;-\\;c\\,\\partial_xP_2\\;-\\;c\\,\\partial_yP_3\\tag{1}\\\] \\\[c\\,\\partial_zP_2=+\\partial_tP_2-aP_2-aP_4\\;+\\;c\\,\\partial_xP_1\\;+\\;c\\,\\partial_yP_4\\tag{2}\\\] \\\[c\\,\\partial_zP_3=+\\partial_tP_3-aP_3+aP_1\\;-\\;c\\,\\partial_xP_4\\;+\\;c\\,\\partial_yP_1\\tag{3}\\\] \\\[c\\,\\partial_zP_4=-\\partial_tP_4-aP_4-aP_2\\;+\\;c\\,\\partial_xP_3\\;-\\;c\\,\\partial_yP_2\\tag{4}\\\]

The same four lines as a matrix acting on the column \\(P=(P_1,P_2,P_3,P_4)^{\\mathsf T}\\). Each row reads “what multiplies \\(P_1\\), what multiplies \\(P_2\\), what multiplies \\(P_3\\), what multiplies \\(P_4\\)”:

\\\[c\\,\\partial_z\\begin{pmatrix}P_1\\\\P_2\\\\P_3\\\\P_4\\end{pmatrix}=\\begin{pmatrix}-\\partial_t-a & -c\\,\\partial_x & a-c\\,\\partial_y & 0\\\\ c\\,\\partial_x & +\\partial_t-a & 0 & -a+c\\,\\partial_y\\\\ a+c\\,\\partial_y & 0 & +\\partial_t-a & -c\\,\\partial_x\\\\ 0 & -a-c\\,\\partial_y & c\\,\\partial_x & -\\partial_t-a\\end{pmatrix}\\begin{pmatrix}P_1\\\\P_2\\\\P_3\\\\P_4\\end{pmatrix}.\\\]

Check row 1 times the column: \\((-\\partial_t-a)P_1-c\\,\\partial_xP_2+(a-c\\,\\partial_y)P_3+0\\), which is (1). Row 4: \\(0+(-a-c\\,\\partial_y)P_2+c\\,\\partial_xP_3+(-\\partial_t-a)P_4\\), which is (4). Rows 2 and 3 the same way.

Split the matrix by what each entry does. Five parts, where III.5a had three:

\\\[\\underbrace{\\begin{pmatrix}-\\partial_t&0&0&0\\\\0&+\\partial_t&0&0\\\\0&0&+\\partial_t&0\\\\0&0&0&-\\partial_t\\end{pmatrix}}\_{\\text{advance}}+\\underbrace{\\begin{pmatrix}-a&0&0&0\\\\0&-a&0&0\\\\0&0&-a&0\\\\0&0&0&-a\\end{pmatrix}}\_{\\text{loss}}+\\underbrace{\\begin{pmatrix}0&0&a&0\\\\0&0&0&-a\\\\a&0&0&0\\\\0&-a&0&0\\end{pmatrix}}\_{\\text{exchange}}\\\]

\\\[+\\underbrace{c\\,\\partial_x\\begin{pmatrix}0&-1&0&0\\\\1&0&0&0\\\\0&0&0&-1\\\\0&0&1&0\\end{pmatrix}}\_{x\\text{-turn}}+\\underbrace{c\\,\\partial_y\\begin{pmatrix}0&0&-1&0\\\\0&0&0&1\\\\1&0&0&0\\\\0&-1&0&0\\end{pmatrix}}\_{y\\text{-turn}} .\\\]

Look at the shapes before going on. The advance is diagonal, \\(-\\partial_t\\) on \\(P_1,P_4\\) and \\(+\\partial_t\\) on \\(P_2,P_3\\): forward and backward in \\(t\\), as in III.5a. The loss is \\(-a\\) times the identity, as there. The exchange is symmetric and couples \\(P_1\\) with \\(P_3\\) and \\(P_2\\) with \\(P_4\\); it squares to the identity, a flip. (The \\(-a\\) of the second pair is a sign convention on \\(P_4\\) and nothing else; \\(P_4\\to-P_4\\) makes it \\(+a\\).) The two new matrices are antisymmetric, entry \\((i,j)\\) minus entry \\((j,i)\\), and each squares to minus the identity: quarter-turns, III.8’s turn matrix, one between the pairs \\((P_1,P_2)\\) and \\((P_3,P_4)\\) driven by \\(\\partial_x\\), one between \\((P_1,P_3)\\) and \\((P_2,P_4)\\) driven by \\(\\partial_y\\).

**(B) Add \\(a\\,P\\) to both sides.** On the right it cancels the loss:

\\\[c\\,\\partial_zP+aP=\\big\[\\text{advance}+\\text{exchange}+x\\text{-turn}+y\\text{-turn}\\big\]P .\\\]

**(C) The identity for the left side.** For any column, product rule:

\\\[e^{-az/c}\\,c\\,\\partial_z\\big(e^{az/c}P\\big)=e^{-az/c}\\Big(c\\cdot\\tfrac ac\\,e^{az/c}P+e^{az/c}\\,c\\,\\partial_zP\\Big)=aP+c\\,\\partial_zP .\\\]

The exponent \\(a/c\\) is the one value for which the right side is exactly the left side of (B), as in III.5a (C).

**(D) Put (C) into (B).**

\\\[e^{-az/c}\\,c\\,\\partial_z\\big(e^{az/c}P\\big)=\\big\[\\text{advance}+\\text{exchange}+x\\text{-turn}+y\\text{-turn}\\big\]P .\\\]

**(E) Multiply both sides by \\(e^{az/c}\\).** On the right it passes through the bracket, because the bracket contains only \\(\\partial_t\\), \\(\\partial_x\\), \\(\\partial_y\\) and constants, and \\(e^{az/c}\\) depends on \\(z\\) alone. This is the one step that three dimensions could have spoiled and does not: no \\(\\partial_z\\) sits inside the bracket.

\\\[c\\,\\partial_z\\big(e^{az/c}P\\big)=\\big\[\\text{advance}+\\text{exchange}+x\\text{-turn}+y\\text{-turn}\\big\]\\big(e^{az/c}P\\big).\\\]

**(F) The same column stands on both sides.** It is \\(Q=e^{az/c}P\\), four components \\(Q_1,\\ldots,Q_4\\), and (E) reads

\\\[c\\,\\partial_z\\begin{pmatrix}Q_1\\\\Q_2\\\\Q_3\\\\Q_4\\end{pmatrix}=\\begin{pmatrix}-\\partial_t & -c\\,\\partial_x & a-c\\,\\partial_y & 0\\\\ c\\,\\partial_x & +\\partial_t & 0 & -a+c\\,\\partial_y\\\\ a+c\\,\\partial_y & 0 & +\\partial_t & -c\\,\\partial_x\\\\ 0 & -a-c\\,\\partial_y & c\\,\\partial_x & -\\partial_t\\end{pmatrix}\\begin{pmatrix}Q_1\\\\Q_2\\\\Q_3\\\\Q_4\\end{pmatrix},\\\]

i.e.

\\\[c\\,\\partial_zQ_1=-\\partial_tQ_1+aQ_3-c\\,\\partial_xQ_2-c\\,\\partial_yQ_3,\\qquad c\\,\\partial_zQ_2=+\\partial_tQ_2-aQ_4+c\\,\\partial_xQ_1+c\\,\\partial_yQ_4,\\\] \\\[c\\,\\partial_zQ_3=+\\partial_tQ_3+aQ_1-c\\,\\partial_xQ_4+c\\,\\partial_yQ_1,\\qquad c\\,\\partial_zQ_4=-\\partial_tQ_4-aQ_2+c\\,\\partial_xQ_3-c\\,\\partial_yQ_2 .\\\]

The loss is gone; advance, exchange and the two turns remain. This is the form the process evolves in.

**(G) Solve each line for \\(\\partial_t\\).** Lines 1 and 4: move \\(-\\partial_tQ\\) to the left and \\(c\\,\\partial_zQ\\) to the right. Lines 2 and 3: the same, then multiply the whole line by \\(-1\\):

\\\[\\partial_tQ_1=-c\\,\\partial_zQ_1-c\\,\\partial_xQ_2-c\\,\\partial_yQ_3+aQ_3,\\qquad \\partial_tQ_2=+c\\,\\partial_zQ_2-c\\,\\partial_xQ_1-c\\,\\partial_yQ_4+aQ_4,\\\] \\\[\\partial_tQ_3=+c\\,\\partial_zQ_3+c\\,\\partial_xQ_4-c\\,\\partial_yQ_1-aQ_1,\\qquad \\partial_tQ_4=-c\\,\\partial_zQ_4+c\\,\\partial_xQ_3-c\\,\\partial_yQ_2-aQ_2 .\\\]

**(H) Stack (G).**

\\\[\\partial_t\\begin{pmatrix}Q_1\\\\Q_2\\\\Q_3\\\\Q_4\\end{pmatrix}=\\begin{pmatrix}-c\\,\\partial_z & -c\\,\\partial_x & a-c\\,\\partial_y & 0\\\\ -c\\,\\partial_x & +c\\,\\partial_z & 0 & a-c\\,\\partial_y\\\\ -a-c\\,\\partial_y & 0 & +c\\,\\partial_z & c\\,\\partial_x\\\\ 0 & -a-c\\,\\partial_y & c\\,\\partial_x & -c\\,\\partial_z\\end{pmatrix}\\begin{pmatrix}Q_1\\\\Q_2\\\\Q_3\\\\Q_4\\end{pmatrix}.\\\]

**(I) Split that matrix into its four parts** (what multiplies \\(\\partial_z\\), \\(\\partial_x\\), \\(\\partial_y\\), and \\(a\\)):

\\\[=-c\\begin{pmatrix}1&0&0&0\\\\0&-1&0&0\\\\0&0&-1&0\\\\0&0&0&1\\end{pmatrix}\\partial_z-c\\begin{pmatrix}0&1&0&0\\\\1&0&0&0\\\\0&0&0&-1\\\\0&0&-1&0\\end{pmatrix}\\partial_x-c\\begin{pmatrix}0&0&1&0\\\\0&0&0&1\\\\1&0&0&0\\\\0&1&0&0\\end{pmatrix}\\partial_y+a\\begin{pmatrix}0&0&1&0\\\\0&0&0&1\\\\-1&0&0&0\\\\0&-1&0&0\\end{pmatrix}.\\\]

The four matrices are, in this order, \\(A_z\\), \\(A_x\\), \\(A_y\\), \\(T\\): the \\(A_z\\), \\(A_x\\) and \\(T\\) of IV.4 exactly, and an \\(A_y\\) that is real (the real \\(4\\times4\\) representation is Majorana’s, *Nuovo Cimento* 14, 171 (1937)). One label to keep straight: in the clock \\(z\\) the time sense is a property of each *counter* — \\(1,4\\) forward, \\(2,3\\) backward, the sign of \\(dt/dz\\), and the mass couples \\(1\\leftrightarrow3\\), \\(2\\leftrightarrow4\\), forward to backward; IV.4’s “the two blocks are the two time senses” is the \\(t\\)-reader’s name for the same pairing, where the block is the chirality and the backward counter of a block is its antiparticle. Same matrices, two names for the index. Now the three \\(A\\)’s are symmetric and \\(T\\) is antisymmetric. In (A) it was the other way round: the exchange was symmetric and the two transverse matrices antisymmetric. (J) says how the roles swapped.

**(J) Where the signs came from, as a matrix identity.** Going from (F) to (H) multiplied lines 2 and 3 by \\(-1\\): left-multiplication by the direction matrix \\(A_z=\\mathrm{diag}(1,-1,-1,1)\\). Apply it to each matrix of (F). To the exchange:

\\\[\\begin{pmatrix}1&0&0&0\\\\0&-1&0&0\\\\0&0&-1&0\\\\0&0&0&1\\end{pmatrix}\\begin{pmatrix}0&0&1&0\\\\0&0&0&-1\\\\1&0&0&0\\\\0&-1&0&0\\end{pmatrix}=\\begin{pmatrix}0&0&1&0\\\\0&0&0&1\\\\-1&0&0&0\\\\0&-1&0&0\\end{pmatrix}=T .\\\]

Symmetric in, antisymmetric out: the flip of the clock is the turn of the \\(t\\)-reader, III.5a (J) twice over. To the \\(x\\)-turn:

\\\[\\begin{pmatrix}1&0&0&0\\\\0&-1&0&0\\\\0&0&-1&0\\\\0&0&0&1\\end{pmatrix}\\begin{pmatrix}0&-1&0&0\\\\1&0&0&0\\\\0&0&0&-1\\\\0&0&1&0\\end{pmatrix}=\\begin{pmatrix}0&-1&0&0\\\\-1&0&0&0\\\\0&0&0&1\\\\0&0&1&0\\end{pmatrix}=-A_x .\\\]

Antisymmetric in, symmetric out: the turn of the clock is an advance for the \\(t\\)-reader. The \\(y\\)-turn gives \\(-A_y\\) the same way. One multiplication by the direction matrix does the whole swap. In the clock, moving sideways is turning; read in \\(t\\), it is just moving. So the process form holds the quarter-turns and the Dirac form the half-turns — the rotation matrices and the spin matrices are the same terms in two clocks (IV.7a).

**(K) The squares.**

\\\[T^2=\\begin{pmatrix}0&0&1&0\\\\0&0&0&1\\\\-1&0&0&0\\\\0&-1&0&0\\end{pmatrix}^2=\\begin{pmatrix}-1&0&0&0\\\\0&-1&0&0\\\\0&0&-1&0\\\\0&0&0&-1\\end{pmatrix},\\\]

the \\(-a^2\\) of Klein–Gordon, and \\(A_x^2=A_y^2=A_z^2=1\\) with all six pairs anticommuting (multiply any two), so that applying the operator of (H) twice gives \\(\\partial_t^2Q_k=c^2\\nabla^2Q_k-a^2Q_k\\) for each of the four components, IV.4a’s two demands met. This is the check that (A) was the right system.

*The real \\(A_y\\).* IV.3 said the flatlander could not write the third half-turn with real numbers: in \\(2\\times2\\) the only matrix anticommuting with \\(\\sigma_x\\) and \\(\\sigma_z\\) is \\(\\sigma_y\\), imaginary. In \\(4\\times4\\) he can. The \\(A_y\\) of (I) is real, and it is IV.4’s imaginary \\(A_y\\) after the unitary change of basis \\(U=(1+A_yA_y^{\\mathrm{IV.4}})/\\sqrt2\\), which leaves \\(A_x\\), \\(A_z\\) and \\(T\\) alone (checked). Nothing in (A)–(K) needed an \\(i\\): the four lines are real, their solutions are real, and III.9’s sentence holds in three dimensions as it did in one — the classical system is the real slice of the quantum one, and what the quantum step adds is the freedom for the components to be complex.

*What the turns do.* In position space the \\(x\\)-turn of (F) reads, for a step \\(dz\\) of clock, \\(Q_1(z+dz)\\ni-dz\\,\\partial_xQ_2\\approx\\tfrac12\\big\[Q_2(x-dz)-Q_2(x+dz)\\big\]\\): what arrives in component 1 from component 2 counts with a plus if it came from the left and a minus if it came from the right. A step through \\(x\\) carries a sign. The same term seen by the \\(t\\)-reader, after (J), is \\(-cA_x\\partial_x\\) with \\(A_x\\) symmetric; on a plane wave \\(e^{ik_xx}\\) it turns the column at the rate \\(ck_x\\) in \\(t\\), which is IV.10’s precession. Seen in the clock, before (J), the same plane wave meets \\(ik_x\\) times an antisymmetric matrix, whose eigenvalues are real, \\(\\pm k_x\\): the pair grows or shrinks along \\(z\\). One term, a rotation in \\(t\\) and a growth in \\(z\\); that is what “a turn in the clock” means for a mode. And the \\(x\\)-turn pairs \\(Q_1\\) with \\(Q_2\\), a forward-in-\\(t\\) component with a backward one, a \\(+z\\) mover with a \\(-z\\) mover: to move in \\(x\\) the particle is turning between \\(+z\\) and \\(-z\\), which is IV.2a’s “direction is not a label” as a mechanism.

*What the \\(P\\)’s are.* Probabilities exactly as long as nothing depends on \\(x\\) and \\(y\\): then (A) is two walks of III.1 along \\(z\\) — the pair \\((P_1,P_3)\\) as written, the pair \\((P_2,P_4)\\) after the relabelling \\(P_4\\to-P_4\\), since (2) and (4) carry the exchange with \\(-a\\) where III.2 has \\(+a\\); the sign is a convention of the counters, fixed here before the \\(P\\)’s are called probabilities. With a transverse gradient the turns put signs on the hops, and the growth just found is the \\(z\\)-clock’s growing mode, \\(-a+\\sqrt{a^2+k\_\\perp^2}>0\\) (the square of the clock operator is \\(a^2+k\_\\perp^2-\\omega^2\\) times the identity), which no probability process has. So in three dimensions the clock-axis walk is a walk and the transverse motion is not; it is what three walks sharing one column cost, and the cost is paid in turns. The signs on the paths of IV.10’s three-coin step are these. (`scripts/four_lines_clock.py`: the matrix form against the four lines, the shapes in (A), the identity (C), (J) and (K) as matrix products, Klein–Gordon symbolically, the unitary to IV.4’s \\(A_y\\), and the reduction to III.2 at \\(\\partial_x=\\partial_y=0\\).)

### IV.4c Where the signs come from: the one-dimensional walk seen by a runner

*The picture first.* IV.4b found the four lines and said that their \\(\\partial_x\\) and \\(\\partial_y\\) parts are not a walk: they carry signs, and IV.4a’s two demands fix them. This section gives those signs a physical source.

Start from the one-dimensional walk along \\(z\\) and nothing else: two time senses, each with its advance, loss and flip, at fixed \\(x\\) and \\(y\\). Watch it from someone running along \\(x\\). The runner’s time is tilted into \\(x\\). What was a pure move forward in \\(t\\) is, for him, also a step toward \\(+x\\), and a move backward in \\(t\\) is also a step toward \\(-x\\).

So the flip of the walk looks different to him. In the flip, a path that was going forward in \\(t\\) goes backward. The \\(t\\)-reader sees this as a \\(+z\\) mover and a \\(-z\\) mover meeting and annihilating. The runner sees a path that drifts toward \\(+x\\) and then back toward \\(-x\\): a sideways U-turn. In the advance, backward moves count \\(-1\\), so the runner’s sideways U-turns inherit that \\(-1\\).

Asking that the runner see the same kind of walk fixes the sideways part completely. A second runner, along \\(y\\), does the same thing in \\(y\\). The two runners together leave a turn in the \\((x,y)\\) plane. That turn is why two paths that make their \\(x\\) and \\(y\\) U-turns in opposite orders cancel. Two counters have no room for that turn, and so there are four.

*Symbols introduced here:* - \\(v\\) is the runner’s speed, small, so everything is to first order in \\(v/c\\). - \\(K_x\\) and \\(K_y\\) say what a run along \\(x\\) or along \\(y\\) does to the counters. - \\(W_x\\) and \\(W_y\\) are the matrices that multiply \\(c\\,\\partial_x\\) and \\(c\\,\\partial_y\\) in the clock form, i.e. IV.4b’s \\(x\\)-turn and \\(y\\)-turn. - \\(\\varepsilon_i\\) is IV.4’s time sense of counter \\(i\\): \\(+1\\) forward, \\(-1\\) backward.

**(A) What is given.** The one-dimensional walk along \\(z\\), with the loss taken out as in IV.4b (B)–(F), is

\\\[c\\,\\partial_zQ=-A_z\\,\\partial_tQ+a\\,(A_zT)\\,Q,\\qquad A_z=\\mathrm{diag}(\\varepsilon_1,\\varepsilon_2,\\ldots).\\\]

\\(-A_z\\partial_t\\) is the advance (forward \\(-\\partial_t\\), backward \\(+\\partial_t\\)), and \\(A_zT\\) is the exchange. For one walk the column has two counters, with \\(A_z=\\mathrm{diag}(1,-1)\\) and \\(T=\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}\\), so the exchange is \\(\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix}\\). Now allow motion in \\(x\\) and \\(y\\), with matrices that are not known:

\\\[c\\,\\partial_zQ=-A_z\\,\\partial_tQ+a\\,(A_zT)\\,Q+c\\,W_x\\,\\partial_xQ+c\\,W_y\\,\\partial_yQ .\\\]

Nothing is assumed about \\(W_x\\) and \\(W_y\\).

**(B) The runner along \\(x\\).** To first order, his coordinates are

\\\[t'=t+\\frac{v}{c^2}\\,x,\\qquad x'=x+v\\,t .\\\]

A forward time move, \\(\\Delta t=+\\Delta z/c\\) with \\(\\Delta x=0\\), is for him also \\(\\Delta x'=+v\\,\\Delta z/c\\). A backward time move is \\(\\Delta x'=-v\\,\\Delta z/c\\). So a flip, which turns forward into backward, is for him a drift toward \\(+x'\\) followed by a drift toward \\(-x'\\), which is a sideways U-turn. In the other direction, a sideways step \\(\\Delta x\\) is for him also a time step \\(\\Delta t'=v\\,\\Delta x/c^2\\).

The derivatives change accordingly:

\\\[\\partial_t=\\partial\_{t'}+v\\,\\partial\_{x'},\\qquad \\partial_x=\\partial\_{x'}+\\frac{v}{c^2}\\,\\partial\_{t'} .\\\]

**(C) What the runner writes down.** Put (B) into (A):

\\\[c\\,\\partial_zQ=-A_z\\partial\_{t'}Q\\;\\underbrace{-\\,v\\,A_z\\,\\partial\_{x'}Q}\_{\\text{time moves seen sideways}}\\;+\\;c\\,W_x\\partial\_{x'}Q\\;\\underbrace{+\\,\\frac vc\\,W_x\\,\\partial\_{t'}Q}\_{\\text{sideways moves seen in time}}\\;+\\;c\\,W_y\\partial_yQ+a\\,(A_zT)\\,Q .\\\]

His readings are \\(Q'=(1+\\tfrac vcK_x)\\,Q\\). Multiply the line by \\((1+\\tfrac vcK_x)\\) and write it in \\(Q'\\). Call the bracket of (A) \\(O\\). To first order this adds \\(\\tfrac vc(K_xO-OK_x)\\). Collect what multiplies each derivative:

\\\[\\partial\_{t'}:\\quad \\frac vc\\Big\[\\,W_x-(K_xA_z-A_zK_x)\\,\\Big\]\\\] \\\[\\partial\_{x'}:\\quad v\\Big\[\\,K_xW_x-W_xK_x-A_z\\,\\Big\]\\\] \\\[\\partial\_{y}:\\quad v\\Big\[\\,K_xW_y-W_yK_x\\,\\Big\]\\\] \\\[\\text{exchange}:\\quad \\frac vc\\,a\\Big\[\\,K_x(A_zT)-(A_zT)K_x\\,\\Big\]\\\]

The runner must see the same four kinds of term with the same matrices: the same advance, the same exchange at the same rate \\(a\\), and the same turns. This is universality, and it says that every bracket vanishes:

\\\[\\text{(C1)}\\ \\ W_x=K_xA_z-A_zK_x\\qquad\\text{(C2)}\\ \\ K_xW_x-W_xK_x=A_z\\qquad\\text{(C3)}\\ \\ K_x(A_zT)=(A_zT)K_x\\qquad\\text{(C4)}\\ \\ K_xW_y=W_yK_x\\\]

**(D) What (C1) says: the sign is the time sense.** Since \\(A_z\\) is diagonal, (C1) can be read entry by entry:

\\\[(W_x)\_{ij}=(K_x)\_{ij}\\,\\varepsilon_j-\\varepsilon_i\\,(K_x)\_{ij}=(K_x)\_{ij}\\,(\\varepsilon_j-\\varepsilon_i).\\\]

- **The diagonal is zero.** The runner can produce no carried sideways move; every sideways move he sees is a switch. An ordinary walk, in which each counter keeps its own \\(x\\)-direction, cannot come out.
- **Same time sense gives zero.** Between two counters of the same time sense the entry is also zero. Sideways moves connect a forward counter to a backward one, and nothing else.
- **Opposite time senses give \\(\\pm2(K_x)\_{ij}\\).** The sign of the entry is the time sense of the counter the switch leaves, minus the time sense of the counter it lands on.

Read it in one clock step \\(dz\\), as in IV.4b’s *What the turns do*:

\\\[Q_i(z+dz)\\ni dz\\,(W_x)\_{ij}\\,\\partial_xQ_j\\approx\\tfrac12(W_x)\_{ij}\\big\[Q_j(x+dz)-Q_j(x-dz)\\big\].\\\]

So a U-turn from \\(j\\) to \\(i\\) that arrives from the left weighs \\(-\\tfrac12(W_x)\_{ij}=(K_x)\_{ij}\\,\\varepsilon_i\\). Its sign is the time sense of the counter it lands on, times the sign of \\((K_x)\_{ij}\\). A U-turn landing on a backward counter takes the \\(-1\\) that backward moves carry in the advance.

**(E) One runner fits in two counters: the 2+1 walk.** Take the two counters of one walk.

- (C3) with the exchange \\(\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix}\\) forces \\(K_x=\\begin{pmatrix}k&m\\\\m&k\\end{pmatrix}\\).
- Then (C1), from (D), gives \\((W_x)\_{+-}=m(-1-1)=-2m\\) and \\((W_x)\_{-+}=m(1-(-1))=2m\\), so \\(W_x=\\begin{pmatrix}0&-2m\\\\2m&0\\end{pmatrix}\\).
- For (C2), multiply:

\\\[K_xW_x=\\begin{pmatrix}2m^2&-2km\\\\2km&-2m^2\\end{pmatrix},\\qquad W_xK_x=\\begin{pmatrix}-2m^2&-2mk\\\\2mk&2m^2\\end{pmatrix},\\qquad K_xW_x-W_xK_x=\\begin{pmatrix}4m^2&0\\\\0&-4m^2\\end{pmatrix}=A_z .\\\]

So \\(m^2=\\tfrac14\\). The part \\(k\\,1\\) multiplies both counters by the same factor and does nothing. With \\(m=\\tfrac12\\):

\\\[K_x=\\tfrac12\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix},\\qquad W_x=\\begin{pmatrix}0&-1\\\\1&0\\end{pmatrix}.\\\]

The choice \\(m=-\\tfrac12\\) only exchanges \\(+x\\) and \\(-x\\). The lines are

\\\[c\\,\\partial_zQ\_+=-\\partial_tQ\_++aQ\_--c\\,\\partial_xQ\_-,\\qquad c\\,\\partial_zQ\_-=+\\partial_tQ\_-+aQ\_++c\\,\\partial_xQ\_+ .\\\]

These are the 2+1 lines, obtained here rather than assumed. As a walk, the coins of the sideways move are:

| arrives from | lands on \\(+\\) (forward) | lands on \\(-\\) (backward) |
| --- | --- | --- |
| the left (\\(x-\\ell\\)) | \\(+\\tfrac12\\) | \\(-\\tfrac12\\) |
| the right (\\(x+\\ell\\)) | \\(-\\tfrac12\\) | \\(+\\tfrac12\\) |

Each counter also keeps \\(\\tfrac12\\) from each side. One clock step of this walk reproduces the 2+1 \\(P\\) lines exactly.

**(F) A second runner, along \\(y\\).** He gives the same four conditions with \\(y\\) in place of \\(x\\). In particular, he tilts time into \\(y\\) and leaves \\(x\\) alone, so he must see the \\(x\\)-turn unchanged:

\\\[\\text{(C1y)}\\ \\ W_y=K_yA_z-A_zK_y,\\quad\\text{(C2y)}\\ \\ K_yW_y-W_yK_y=A_z,\\quad\\text{(C3y)}\\ \\ K_y(A_zT)=(A_zT)K_y,\\quad\\text{(C4y)}\\ \\ K_yW_x=W_xK_y .\\\]

**(G) Two runners leave a turn in the \\((x,y)\\) plane.**

*In directions.* The \\(x\\)-runner’s walk turns in the \\((z,x)\\) plane: an \\(x\\)-U-turn turns a \\(+z\\) mover into a \\(-z\\) mover while it steps in \\(x\\). The \\(y\\)-runner’s walk turns in the \\((z,y)\\) plane. Make a small \\((z,x)\\) turn, then a small \\((z,y)\\) turn, then undo the first, then undo the second. You do not come back to where you started: what is left over is a small turn in the \\((x,y)\\) plane. For the matrices that turn the directions \\((ct,x,y,z)\\),

\\\[\[\\text{turn}\_{zx},\\,\\text{turn}\_{zy}\]=-\\,\\text{turn}\_{xy}.\\\]

*In the walk.* An \\(x\\)-U-turn followed by a \\(y\\)-U-turn is \\(W_yW_x\\), and the other order is \\(W_xW_y\\). Every solution of (C1)–(C4) together with (C1y)–(C4y) has

\\\[W_yW_x=2\\,(K_xK_y-K_yK_x)=-\\,W_xW_y,\\qquad W_x^2=W_y^2=-1 ,\\\]

and this is derived, not found: (C1) gives \\(\\{A_z,W_x\\}=\\{A_z,\[K_x,A_z\]\\}=\[K_x,A_z^2\]=0\\); commuting that with \\(K_x\\) and using (C2), \\(0=\\{\[K_x,A_z\],W_x\\}+\\{A_z,\[K_x,W_x\]\\}=2W_x^2+2\\); commuting \\(\\{A_z,W_y\\}=0\\) with \\(K_x\\) and using (C4), \\(\\{W_x,W_y\\}=0\\). Then \\(A_x\\equiv-A_zW_x\\) and \\(A_y\\equiv-A_zW_y\\) square to \\(1\\) and anticommute with each other and with \\(A_z\\); \\(K_x-A_x/2\\) commutes with all three, so it lies in the span of \\(1\\) and \\(\\omega\\equiv A_zA_xA_y\\) (the commutant of three anticommuting \\(4\\times4\\) matrices of square one; \\(\\omega^2=-1\\)), both of which commute with \\(A_x\\) and \\(A_y\\): \\(K_x=A_x/2+k\\,1+\\beta\\,\\omega\\). (C3) then says two things: \\(\\{A_x,T\\}=0\\) — the mass turn anticommutes with the new direction, as IV.4 requires — and \\(\\beta=0\\). In steps: a part of \\(T\\) that commuted with \\(\\omega\\) would lie in the algebra of the three \\(A\\)’s, anticommute with \\(A_z\\) and, by (C3) for the two runners, with \\(A_x\\) and \\(A_y\\), and nothing in that algebra does; so \\(T\\omega=-\\omega T\\), and (C3) reads \\(\\tfrac12\\{A_x,T\\}=2\\beta\\,\\omega T\\). Multiply by \\(T^{-1}\\): \\(\\tfrac12(A_x+TA_xT^{-1})=2\\beta\\,\\omega\\); conjugate by \\(T\\): the left side is unchanged and the right side changes sign, so \\(\\beta=0\\) and \\(\\{A_x,T\\}=0\\). (\\(A_zT\\), the exchange, is symmetric with square \\(+1\\) and anticommutes with \\(\\omega\\); the two are not the same matrix.) Either way \\(2\[K_x,K_y\]=\\tfrac12\[A_x,A_y\]=A_xA_y=W_yW_x\\), since \\(\\omega\\) commutes with \\(A_x\\) and \\(A_y\\) (`runner_identities_proof.py`: the operator identities symbolically, the rest on 147 solutions from random starts; the derived \\(A_x\\), \\(A_y\\) need not be symmetric — a solution is IV.4b’s up to a change of basis that keeps \\(A_z\\) and \\(A_zT\\)). The random-start search of `runner_section_3plus1.py` (f) confirms what the four lines force.

\\(K_xK_y-K_yK_x\\) is the turn of the \\((x,y)\\) plane that the two runs produce on the counters. So the two orders give the two senses of that turn, and they differ by a full sign.

**(H) Mirror cancellation.** Two paths with the same moves taken in opposite order, one \\(x\\)-U-turn and one \\(y\\)-U-turn from the same sides, arrive with \\(W_yW_x\\) and \\(W_xW_y=-W_yW_x\\). The side factors \\(\\pm\\tfrac12\\) are the same for both. So the two paths cancel. The cancellation is the \\((x,y)\\) turn of (G), and that turn was produced by the two runners together. Nothing is chosen.

**(I) Four counters are forced.** In two counters, (C1) with (C3) makes every sideways matrix antisymmetric. There is only one antisymmetric \\(2\\times2\\) matrix up to a factor, \\(\\begin{pmatrix}0&-1\\\\1&0\\end{pmatrix}\\). So \\(W_x\\) and \\(W_y\\) would be multiples of it, they would commute, and the two orders in (H) would be equal instead of opposite. Solving (C1y)–(C4y) directly in two counters, with \\(W_x\\) from (E), gives no solution.

The next column is two one-dimensional walks side by side. This is IV.4b (A) at \\(\\partial_x=\\partial_y=0\\), “III.2 twice over”, with advance \\(\\mathrm{diag}(1,-1,-1,1)\\) and exchange \\(A_zT\\).

**(J) The four-counter solution is IV.4b’s four lines.** Take \\(K_x=\\tfrac12A_x\\) and \\(K_y=\\tfrac12A_y\\), with the real matrices of IV.4b (I). Since \\(A_x\\) anticommutes with \\(A_z\\), (C1) gives

\\\[W_x=\\tfrac12A_xA_z-\\tfrac12A_zA_x=A_xA_z=\\begin{pmatrix}0&-1&0&0\\\\1&0&0&0\\\\0&0&0&-1\\\\0&0&1&0\\end{pmatrix},\\\]

which is the \\(x\\)-turn of IV.4b (A). In the same way,

\\\[W_y=A_yA_z=\\begin{pmatrix}0&0&-1&0\\\\0&0&0&1\\\\1&0&0&0\\\\0&-1&0&0\\end{pmatrix},\\\]

which is the \\(y\\)-turn. The other conditions hold:

- **(C2):** \\(K_xW_x-W_xK_x=\\tfrac12A_xA_xA_z-\\tfrac12A_xA_zA_x=\\tfrac12A_z+\\tfrac12A_z=A_z\\).
- **(C3):** \\(A_x\\) anticommutes with both \\(A_z\\) and \\(T\\), so it commutes with \\(A_zT\\).
- **(C4y):** \\(K_yW_x=\\tfrac12A_yA_xA_z=\\tfrac12A_xA_zA_y=W_xK_y\\), moving \\(A_y\\) past two matrices it anticommutes with.

The \\(y\\) conditions work the same way. With the loss put back, the lines are (1)–(4) of IV.4b (A), signs included.

In (D)’s terms: for the pair \\((1,2)\\), \\((K_x)\_{12}=+\\tfrac12\\), so an \\(x\\)-U-turn from the left counts \\(+\\) when it lands forward. For the pair \\((3,4)\\), \\((K_x)\_{34}=-\\tfrac12\\), and the rule is reversed. This is IV.4’s second block, which streams the opposite way.

**(K) Two U-turns in both orders, worked out.** Move \\(A_z\\) through using \\(A_zA_x=-A_xA_z\\) and \\(A_z^2=1\\):

\\\[W_yW_x=A_yA_zA_xA_z=-A_yA_xA_zA_z=-A_yA_x=A_xA_y,\\qquad W_xW_y=A_xA_zA_yA_z=-A_xA_y .\\\]

\\\[A_xA_y=\\begin{pmatrix}0&0&0&1\\\\0&0&1&0\\\\0&-1&0&0\\\\-1&0&0&0\\end{pmatrix},\\qquad (A_xA_y)^2=-A_xA_xA_yA_y=-1 .\\\]

So \\(\\tfrac12A_xA_y\\) is the \\((x,y)\\) turn on the counters. It works at half the angle: a full turn of the \\((x,y)\\) plane multiplies every counter by \\(-1\\).

One walker, starting on counter 4, with entries read from \\(W_x\\) and \\(W_y\\) above (row = lands on, column = leaves):

| order | first U-turn | second U-turn | arrives on 1 with |
| --- | --- | --- | --- |
| \\(x\\) then \\(y\\) | \\(4\\to3\\): \\(-1\\) | \\(3\\to1\\): \\(-1\\) | \\(+1\\) |
| \\(y\\) then \\(x\\) | \\(4\\to2\\): \\(+1\\) | \\(2\\to1\\): \\(-1\\) | \\(-1\\) |

**(L) What was used.** - **Universality:** the runner sees the same lines (C), and so does the second runner (F). - **Additivity:** the moves within one clock step add, since the line is a sum of advance, exchange, \\(x\\)-turn and \\(y\\)-turn; and the first-order changes add.

Nothing else was used. The signs in IV.4b (A) are not chosen: they are the backward time sense of the one-dimensional walk, tilted into \\(x\\) and \\(y\\) by runners. Finite runs, with \\(\\cosh\\tfrac\\eta2+\\sinh\\tfrac\\eta2\\,A_x\\) on the counters and the full Lorentz change on the coordinates, leave the four lines unchanged. The same holds along \\(y\\).

*What this does not show.* Every solution of the runner conditions satisfies (G) (checked). That every solution is IV.4b’s up to relabelling the counters is not proved here; for the matrices themselves that is IV.4a’s uniqueness statement.

*The constructive side.* The runner’s rule can be run forward. Take a \\(1{+}1\\) walk along \\(z\\) in its own clock, positive (the branching reading, \\(Q=e^{az}P\\): an expected count of a walk that keeps going and spawns the other sense at rate \\(a\\), not a conserved probability), and read it with a runner of rapidity \\(\\eta\\) along \\(x\\): \\(Q\_\\eta(t,x,z)=e^{\\eta A_x/2}\\,Q_0\\big(\\gamma(t-\\beta x),z\\big)\\) solves the \\(2{+}1\\) lines exactly (every \\(1{+}1\\) plane wave, both branches). The piece is a positive process; its reading is not: \\(e^{\\eta A_x/2}=\\cosh\\frac\\eta2+\\sinh\\frac\\eta2A_x\\) has a negative entry for \\(\\eta\<0\\), and the runner reads \\(Q\_-\<0\\) wherever \\(Q\_{0+}>\\coth\\frac{|\\eta|}2\\,Q\_{0-}\\) — on \\(29\\%\\) of the sites in the run, from paths that are all positive. A general \\(2{+}1\\) column is a sum of such sheets with coefficients of both signs, and the unitary \\(t\\)-clock step evolves that sum as it should (one sheet, a sum, a difference: \\(2\\)–\\(3\\%\\) at \\(h=0.025\\), length conserved to \\(10^{-9}\\)). What cannot be done is to run the sum as a positive process in the clock \\(z\\): the transverse turn has real growth rates \\(\\pm k\\) (IV.4b, *What the turns do*), and the \\(2{+}1\\) \\(z\\)-clock walk on transverse data diverges as the spacing shrinks (\\(0.3\\to20\\to5\\times10^5\\) over two halvings). That is not a numerical accident: a surface \\(z=\\) const is timelike in \\(2{+}1\\), data on it pose an ill-posed Cauchy problem, modes with \\(k>\\omega\\) are evanescent in \\(z\\) and grow as \\(e^{kz}\\) with \\(k\\) up to \\(\\pi/h\\) on the lattice; in \\(1{+}1\\) there is no transverse \\(k\\), which is why the walk’s own clock works there. So the \\(2{+}1\\) column is boosted positive \\(1{+}1\\) walks plus the signs special relativity puts on the runners’ readings and on the coefficients — IV.2a’s \\(-2\\) is those signs, and its “small difference of two large counts” is the growth just seen. (`scripts/boosted_sheets_make_2plus1.py`, 1–5.)

(`scripts/runner_section_3plus1.py`: (B) and the first-order collection in (C) for general \\(4\\times4\\) \\(K_x,W_x,W_y\\); (D) entry by entry; (E) solved, and no second runner in two counters; (G) over all solutions found from random starts; (J) against IV.4b (A) line by line; (K) as matrix products and the walker table; finite runs along \\(x\\) and \\(y\\). `scripts/runner_identities_proof.py`: (G)’s identities derived. `scripts/walk_1d_to_2d_by_runs.py` W8: the 2+1 walk step against the 2+1 \\(P\\) lines.)

IV.4d says what the runner of this section is looking at: the streams, the no-flip segments of one path sorted by direction, whose readings add and whose contents do not.

### IV.4d The streams: what the runner of IV.4c is looking at

*The picture first.* Between two flips the walker moves in a straight line at \\(c\\); III.15 called that segment a particle. Collect, at a point, every segment passing through it in the same direction. That collection is a *stream*. In one dimension there are two: the forward segments, which are \\(Q\_+\\), and the backward segments, which are \\(Q\_-\\) and which the forward-time observer calls antiparticles. In two dimensions there is a stream for every direction, and the direction matrices sort the segments into them. There is still one walker. A stream is one path cut at its flips and binned by direction, and the readings are the weights of the bins with the U-turn signs carried along. This section says what the readings, the squares and the runs of IV.4c mean in those terms. Everything is \\(2{+}1\\) in the \\(t\\)-reading, with \\(A_z=\\mathrm{diag}(1,-1)\\) and \\(A_x=\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix}\\), and \\(a=0\\) until (H).

**(A) A stream is a column that a direction matrix gives back.** With \\(a=0\\) the lines are \\(\\partial_tQ=-cA_z\\partial_zQ-cA_x\\partial_xQ\\). Four solutions move at \\(c\\) in one direction without changing shape:

\\\[\\begin{pmatrix}1\\\\0\\end{pmatrix}f\\big(t-\\tfrac zc\\big),\\qquad\\begin{pmatrix}0\\\\1\\end{pmatrix}f\\big(t+\\tfrac zc\\big),\\qquad\\begin{pmatrix}1\\\\1\\end{pmatrix}f\\big(t-\\tfrac xc\\big),\\qquad\\begin{pmatrix}1\\\\-1\\end{pmatrix}f\\big(t+\\tfrac xc\\big).\\\]

Check the third: \\(A_x\\binom11=\\binom11\\), so \\(-cA_x\\partial_x\\) gives \\(+f'\\) on both lines, and \\(\\partial_t\\) gives \\(f'\\). Along \\(z\\) the streams are the single counters; along \\(x\\) they are the sum and the difference of the counters. “The speeds are the eigenvalues” (IV.4a) means exactly this: a stream along a direction is a column that the direction’s matrix returns multiplied by \\(+1\\) (moving at \\(+c\\)) or \\(-1\\) (at \\(-c\\)).

**(B) How much of a stream there is: the content.** Multiply line \\(+\\) by \\(2Q\_+\\), line \\(-\\) by \\(2Q\_-\\), add; the \\(x\\) terms combine into one derivative:

\\\[\\partial_t\\big(Q\_+^2+Q\_-^2\\big)+c\\,\\partial_z\\big(Q\_+^2-Q\_-^2\\big)+c\\,\\partial_x\\big(2Q\_+Q\_-\\big)=0 .\\\]

So the count is \\(Q\_+^2+Q\_-^2\\), the flow along \\(z\\) is \\(Q\_+^2-Q\_-^2\\) and the flow along \\(x\\) is \\(2Q\_+Q\_-\\). Call the three together the *content* of the column, \\(j=(Q\_+^2+Q\_-^2,\\ 2Q\_+Q\_-,\\ Q\_+^2-Q\_-^2)\\) in the order \\((ct,x,z)\\). The \\(+z\\) stream \\((1,0)\\) has content \\((1,0,1)\\): count one, all of it moving along \\(z\\). The \\(+x\\) stream \\((1,1)\\) has content \\((2,2,0)\\).

**(C) Readings add; contents do not.** The lines are linear, so the sum of two streams is a solution, and its reading is the sum of the readings. Its content is not the sum of the contents, because a square has a cross term. The same column is therefore several streams in several sortings:

\\\[Q\_+^2+Q\_-^2=\\tfrac12(Q\_++Q\_-)^2+\\tfrac12(Q\_+-Q\_-)^2 .\\\]

Sorted by \\(z\\) the column is “\\(Q\_+^2\\) forward and \\(Q\_-^2\\) backward”; sorted by \\(x\\) it is “\\(\\tfrac12(Q\_++Q\_-)^2\\) to the right and \\(\\tfrac12(Q\_+-Q\_-)^2\\) to the left”. Neither sorting is more real than the other. This is where the streams of the walk differ from streams of dust: dust adds in the density, these add in the reading, and the density is read off afterwards.

**(D) The runner, from the content.** IV.4c found the half-rapidity on the counters from the lines. Here it comes from the content, with one input: a runner along \\(z\\) with rapidity \\(\\eta\\) sees the content as a flow, so its light-cone parts \\(j^0\\pm j^z\\) scale by \\(e^{\\pm\\eta}\\) and \\(j^x\\) is unchanged. But \\(j^0+j^z=2Q\_+^2\\) and \\(j^0-j^z=2Q\_-^2\\). So

\\\[Q\_+^2\\to e^{\\eta}Q\_+^2,\\qquad Q\_-^2\\to e^{-\\eta}Q\_-^2,\\qquad\\text{hence}\\qquad Q\_+\\to e^{\\eta/2}Q\_+,\\quad Q\_-\\to e^{-\\eta/2}Q\_- .\\\]

The reading gets the square root of the factor its content gets. That is the whole origin of the half: a reading is half of a stream, in the sense that two readings multiplied make one content. (IV.4c reaches the same factors from the lines alone, with no input about the content.)

**(E) The same along \\(x\\): universality.** The \\(x\\)-streams \\(Q\_+\\pm Q\_-\\) obey the same one-dimensional law along \\(x\\) that \\(Q\_\\pm\\) obey along \\(z\\). So a runner along \\(x\\) does to them what the runner along \\(z\\) did to \\(Q\_\\pm\\): \\((Q\_++Q\_-)\\to e^{\\eta/2}(Q\_++Q\_-)\\) and \\((Q\_+-Q\_-)\\to e^{-\\eta/2}(Q\_+-Q\_-)\\). Back in the counters,

\\\[\\begin{pmatrix}Q\_+\\\\Q\_-\\end{pmatrix}\\to\\begin{pmatrix}\\cosh\\frac\\eta2&\\sinh\\frac\\eta2\\\\ \\sinh\\frac\\eta2&\\cosh\\frac\\eta2\\end{pmatrix}\\begin{pmatrix}Q\_+\\\\Q\_-\\end{pmatrix},\\\]

and the content again scales as a vector, now along \\(x\\). Both runs leave the lines unchanged.

**(F) Two runs leave the U-turn coin.** Run a little along \\(z\\), then along \\(x\\), then back along \\(z\\), then back along \\(x\\). On the directions what is left is a small turn of the \\((z,x)\\) plane. On the readings the same four steps leave \\(1+\\varepsilon\\delta\\,J\\) with

\\\[J=K_zK_x-K_xK_z=\\begin{pmatrix}0&\\frac12\\\\-\\frac12&0\\end{pmatrix},\\qquad K_z=\\tfrac12A_z,\\quad K_x=\\tfrac12A_x .\\\]

\\(J\\) turns the readings by half the angle the directions turn, and \\(2J\\) is the sideways coin of IV.4b’s \\(z\\)-clock lines: \\(+\\tfrac12\\) into the forward counter, \\(-\\tfrac12\\) into the backward one. This is IV.4c (G) read on the streams: the signed U-turn is what two runs in different directions leave behind on the readings.

**(G) One stream in a direction is a pair of streams along \\(z\\).** Take the reading \\((\\cos\\frac\\theta2,\\sin\\frac\\theta2)\\). The matrix \\(\\cos\\theta\\,A_z+\\sin\\theta\\,A_x\\) returns it unchanged, so it is the stream moving at \\(+c\\) in the direction \\(\\theta\\) from \\(z\\). Sorted by \\(z\\) the same reading is \\(\\cos^2\\frac\\theta2\\) of a forward stream and \\(\\sin^2\\frac\\theta2\\) of a backward one, and its cross term is \\(2Q\_+Q\_-=\\sin\\theta\\): exactly its flow along \\(x\\). So “one stream in the direction \\(\\theta\\)” and “forward segments mixed with backward segments along \\(z\\)” are the same column in two sortings, and the interference of the forward and backward streams is the sideways motion. The forward-time observer reads the backward part as pairs (III.3); this is the sense in which the streams are a picture of virtual pairs: the pairs are the other sorting.

**(H) The flips move content between streams.** Restore \\(a\\), in one dimension:

\\\[(\\partial_t+c\\,\\partial_z)\\,Q\_+^2=+2aQ\_+Q\_-,\\qquad(\\partial_t-c\\,\\partial_z)\\,Q\_-^2=-2aQ\_+Q\_- .\\\]

What the forward stream gains the backward one loses. With \\(a=0\\) each stream runs on its own forever; the mass is the rate at which segments leave one stream and join the other. This is III.3’s pair creation and annihilation, seen from the streams.

**(I) What the streams are.** Every step above is a real rescaling or a real turn of a single column: one direction of traversal, sorted by direction. What the streams do not contain is the second direction of III.15a; there, the content \\(Q\_+^2+Q\_-^2\\) is the weight of an event between two endpoints, and the streams are what one leg of that event looks like, binned by the direction of its segments.

(`scripts/streams_section.py`: (A)–(H) as written. `scripts/boosts_streams.py`, `scripts/elements_from_boosts.py`, `scripts/two_step_covariance.py`, `scripts/universality_additivity.py`: the runs, the loop, and the content under a boost.)

## IV.5 The quantum step and the Dirac equation

III.9 word for word: four Madelung pairs, multiply by \\(i\\hbar\\). The turn term \\(i\\hbar a\\,T\\) is \\(i\\hbar a\\) times a real antisymmetric matrix, Hermitian for real \\(a\\): nothing to continue, \\(a=mc^2/\\hbar\\). As in III.9, multiply the lower block by \\(i\\) (\\(U=\\mathrm{diag}(1,1,i,i)\\) — in the chiral basis of IV.4, to which one returns first from IV.4b’s real counters: the phase leaves the chiral \\(A_k\\) alone but not the real \\(A_y\\) of IV.4b — which leaves the \\(A_k\\) alone) to reach the standard form with the symmetric \\(F=\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix}\\otimes1\\) (checked):

\\\[i\\hbar\\,\\partial_t\\begin{pmatrix}\\psi_1\\\\\\psi_2\\\\\\psi_3\\\\\\psi_4\\end{pmatrix}=\\Big\[c\\big(A_x\\hat p_x+A_y\\hat p_y+A_z\\hat p_z\\big)+mc^2F\\Big\]\\begin{pmatrix}\\psi_1\\\\\\psi_2\\\\\\psi_3\\\\\\psi_4\\end{pmatrix},\\qquad F=\\begin{pmatrix}0&0&1&0\\\\0&0&0&1\\\\1&0&0&0\\\\0&1&0&0\\end{pmatrix},\\\]

\\\[i\\hbar\\,\\partial_t\\psi_1=-i\\hbar c\\big\[\\partial_z\\psi_1+(\\partial_x-i\\partial_y)\\psi_2\\big\]+mc^2\\psi_3\\\] \\\[i\\hbar\\,\\partial_t\\psi_2=-i\\hbar c\\big\[(\\partial_x+i\\partial_y)\\psi_1-\\partial_z\\psi_2\\big\]+mc^2\\psi_4\\\] \\\[i\\hbar\\,\\partial_t\\psi_3=+i\\hbar c\\big\[\\partial_z\\psi_3+(\\partial_x-i\\partial_y)\\psi_4\\big\]+mc^2\\psi_1\\\] \\\[i\\hbar\\,\\partial_t\\psi_4=+i\\hbar c\\big\[(\\partial_x+i\\partial_y)\\psi_3-\\partial_z\\psi_4\\big\]+mc^2\\psi_2\\\]

**The Dirac equation.** The mass couples \\(1\\leftrightarrow3\\), \\(2\\leftrightarrow4\\): it reverses the time sense \\(\\varepsilon\\) at the real rate \\(mc^2/\\hbar\\), the Compton frequency — “the mass is the reversal rate” of III.9 — and with it, by IV.4, the handedness of the circulation. Squaring gives \\(c^2\\hat p^2+m^2c^4\\) on every component: Klein–Gordon for each, as in Part II. As in III.9 remark 5, \\(\\hat H\\) is the half-turn of the four-component vector \\((c\\hat p_x,c\\hat p_y,c\\hat p_z,mc^2)\\) and the energy is its length.

### IV.5a Eight detectors: a time sense and an orientation; and the \\(i\\) is the walker’s own turn

*The picture first.* IV.4b ran the argument of III.5a in three dimensions with four counters. Look at what the four are. Each is a detector for one time sense and one orientation of the walker: forward, as launched; forward, turned in the plane \\((x,y)\\); backward, turned in \\((y,z)\\); backward, turned in \\((z,x)\\). A counter there is forward exactly when its orientation keeps the direction of the clock axis. Time sense and orientation are one label. They are two different things: a walker can be turned over and still go forward in time. A set of detectors that reads both has one for each time sense and each orientation, eight. The rules are IV.4b’s three, unchanged, and the four added detectors obey IV.4b’s four lines again. There is one walker throughout: no second walker, and no process run backward. The only process in the document that runs backward is the back-leg of a measurement (III.15a).

What the four added detectors change. With four, the walker moves at \\(c\\) at every moment. With eight, the walker at rest is the same amount in the forward and the backward detector of each orientation. And the \\(i\\) of IV.5 is a real thing: a turn of the walker in its own plane that leaves its time sense alone, which carries each original detector onto an added one. Written as a real \\(8\\times8\\) matrix it is \\(J\\), every entry \\(0\\) or \\(\\pm1\\).

*Names.* A *detector* is a counter. An *orientation* is how the walker is turned in space relative to the walker as launched (the *reference*, the walker in counter 1 alone). The walker’s *own plane* is the plane in which its clock turns; for the reference it is \\((y,z)\\).

*Symbols.* \\(A_x,A_y,A_z,T\\) are IV.4b’s real matrices. \\(u=(u_1,u_2,u_3,u_4)^{\\mathsf T}\\) are IV.4b’s four counters, \\(v=(v_1,v_2,v_3,v_4)^{\\mathsf T}\\) the four added ones, \\(X=(u,v)\\) the column of eight. \\(\\mathcal A_k\\), \\(\\mathcal T\\) and \\(J\\) are the \\(8\\times8\\) matrices of (B). \\(R\_{yz},R\_{zx},R\_{xy}\\) are written out in (1).

**(A) What the four counters are.** The products of two direction matrices:

\\\[R\_{yz}=A_yA_z=\\begin{pmatrix}0&0&-1&0\\\\0&0&0&1\\\\1&0&0&0\\\\0&-1&0&0\\end{pmatrix},\\qquad R\_{zx}=A_zA_x=\\begin{pmatrix}0&1&0&0\\\\-1&0&0&0\\\\0&0&0&1\\\\0&0&-1&0\\end{pmatrix},\\\]

\\\[R\_{xy}=A_xA_y=\\begin{pmatrix}0&0&0&1\\\\0&0&1&0\\\\0&-1&0&0\\\\-1&0&0&0\\end{pmatrix}.\\tag{1}\\\]

1. each is a turn of the walker in the plane named. Multiply out: \\(R\_{yz}\\) commutes with \\(A_x\\) and anticommutes with \\(A_y\\) and \\(A_z\\), so \\(R\_{yz}^{-1}A_xR\_{yz}=A_x\\), \\(R\_{yz}^{-1}A_yR\_{yz}=-A_y\\), \\(R\_{yz}^{-1}A_zR\_{yz}=-A_z\\): it reverses the two directions of its plane and keeps the third. The other two the same way. Their products:

\\\[R\_{yz}^2=R\_{zx}^2=R\_{xy}^2=-1,\\qquad\\begin{aligned}R\_{zx}R\_{yz}&=+R\_{xy},&R\_{yz}R\_{zx}&=-R\_{xy},\\\\R\_{xy}R\_{zx}&=+R\_{yz},&R\_{zx}R\_{xy}&=-R\_{yz},\\\\R\_{yz}R\_{xy}&=+R\_{zx},&R\_{xy}R\_{yz}&=-R\_{zx}.\\end{aligned}\\tag{2}\\\]

2. a general orientation is a combination \\(q_0+q\_{yz}R\_{yz}+q\_{zx}R\_{zx}+q\_{xy}R\_{xy}\\): by (2) a product of two such is again one, and a turn by \\(\\theta\\) in the plane \\((j,k)\\) is \\(\\cos\\tfrac\\theta2+\\sin\\tfrac\\theta2R\_{jk}\\) (IV.4c (K)). An orientation is four numbers. Turn the reference:

\\\[R\_{yz}\\begin{pmatrix}1\\\\0\\\\0\\\\0\\end{pmatrix}=\\begin{pmatrix}0\\\\0\\\\1\\\\0\\end{pmatrix},\\qquad R\_{zx}\\begin{pmatrix}1\\\\0\\\\0\\\\0\\end{pmatrix}=\\begin{pmatrix}0\\\\-1\\\\0\\\\0\\end{pmatrix},\\qquad R\_{xy}\\begin{pmatrix}1\\\\0\\\\0\\\\0\\end{pmatrix}=\\begin{pmatrix}0\\\\0\\\\0\\\\-1\\end{pmatrix}.\\tag{3}\\\]

3. (the first column of each matrix of (1).) So each counter of IV.4b is the reference in one orientation, with the time sense IV.4b (I) gave it:

| counter | time sense | orientation | keeps the direction of \\(z\\)? | amount in the detector |
| --- | --- | --- | --- | --- |
| \\(u_1\\) | forward | as launched | yes | \\(+u_1\\) |
| \\(u_4\\) | forward | turned in \\((x,y)\\) | yes | \\(-u_4\\) |
| \\(u_3\\) | backward | turned in \\((y,z)\\) | no, reversed | \\(+u_3\\) |
| \\(u_2\\) | backward | turned in \\((z,x)\\) | no, reversed | \\(-u_2\\) |

4. the minus signs are the convention of the counters (IV.4b (A): “the \\(-a\\) of the second pair is a sign convention on \\(P_4\\)”). There is no detector for “turned in \\((y,z)\\) and forward”, and none for “as launched and backward”. Read the three rules off IV.4b (G), as sentences about a time sense and an orientation:

- *Advance.* A forward detector is carried along \\(+z\\), a backward one along \\(-z\\).
- *Time flip.* It connects \\(u_1\\) with \\(u_3\\) and \\(u_4\\) with \\(u_2\\): forward with backward, and the orientations differ by \\(R\_{yz}\\) applied *on the right* (\\(1\\cdot R\_{yz}=R\_{yz}\\), \\(R\_{zx}\\cdot R\_{yz}=R\_{xy}\\)). A turn on the right is made before the walker’s own orientation: a turn in the walker’s own plane.
- *Sideways step.* Along \\(x\\) it connects \\(u_1\\) with \\(u_2\\) and \\(u_4\\) with \\(u_3\\): forward with backward, the orientations differing by \\(R\_{zx}\\) *on the left* (\\(R\_{zx}\\cdot1=R\_{zx}\\), \\(R\_{zx}\\cdot R\_{xy}=-R\_{yz}\\)). Along \\(y\\) it connects \\(u_1\\) with \\(u_3\\) and \\(u_4\\) with \\(u_2\\): \\(R\_{yz}\\) on the left. A turn on the left is made after the walker’s orientation: a turn in the plane as the detectors see it, the plane of the clock axis and the direction of the step (IV.4b: “to move in \\(x\\) the particle is turning between \\(+z\\) and \\(-z\\)”).

5. none of the three sentences asks which orientation the walker has.

**(B) Eight detectors, the same three rules.** One detector for each time sense and each orientation. The four that IV.4b does not have are \\(v_1,\\dots,v_4\\), and \\(v_j\\) is the detector with the time sense of \\(u_j\\) and the orientation of \\(u_j\\) turned in the walker’s own plane (\\(R\_{yz}\\) on the right):

|  | as launched | turned in \\((y,z)\\) | turned in \\((z,x)\\) | turned in \\((x,y)\\) |
| --- | --- | --- | --- | --- |
| forward | \\(+u_1\\) | \\(+v_1\\) | \\(+v_4\\) | \\(-u_4\\) |
| backward | \\(-v_3\\) | \\(+u_3\\) | \\(-u_2\\) | \\(-v_2\\) |

6. (\\(u_1\\): \\(1\\cdot R\_{yz}=R\_{yz}\\), so \\(v_1\\) is forward, turned in \\((y,z)\\). \\(u_2\\): \\(R\_{zx}R\_{yz}=R\_{xy}\\). \\(u_3\\): \\(R\_{yz}R\_{yz}=-1\\), as launched with a minus. \\(u_4\\): \\(R\_{xy}R\_{yz}=-R\_{zx}\\). `detectors_in_master_names.py`, check 1.) The rules of (A) turn the orientation on the left (a sideways step) or on the right by \\(R\_{yz}\\) (a flip); both commute with the one turn on the right that takes \\(u_j\\) to \\(v_j\\). So the added four obey the lines the original four obey, IV.4b (H):

\\\[\\partial_tu=\\mathcal Lu,\\qquad\\partial_tv=\\mathcal Lv,\\qquad\\mathcal L=-c\\,(A_x\\partial_x+A_y\\partial_y+A_z\\partial_z)+a\\,T.\\tag{4}\\\]

7. stacked, with \\(X=(u,v)\\) the column of eight:

\\\[X=\\begin{pmatrix}u_1\\\\u_2\\\\u_3\\\\u_4\\\\v_1\\\\v_2\\\\v_3\\\\v_4\\end{pmatrix},\\qquad\\partial_tX=\\begin{pmatrix}\\mathcal L&0\\\\0&\\mathcal L\\end{pmatrix}X=\\begin{pmatrix}\\mathcal L\_{11}&\\mathcal L\_{12}&\\mathcal L\_{13}&\\mathcal L\_{14}&0&0&0&0\\\\\\mathcal L\_{21}&\\mathcal L\_{22}&\\mathcal L\_{23}&\\mathcal L\_{24}&0&0&0&0\\\\\\mathcal L\_{31}&\\mathcal L\_{32}&\\mathcal L\_{33}&\\mathcal L\_{34}&0&0&0&0\\\\\\mathcal L\_{41}&\\mathcal L\_{42}&\\mathcal L\_{43}&\\mathcal L\_{44}&0&0&0&0\\\\0&0&0&0&\\mathcal L\_{11}&\\mathcal L\_{12}&\\mathcal L\_{13}&\\mathcal L\_{14}\\\\0&0&0&0&\\mathcal L\_{21}&\\mathcal L\_{22}&\\mathcal L\_{23}&\\mathcal L\_{24}\\\\0&0&0&0&\\mathcal L\_{31}&\\mathcal L\_{32}&\\mathcal L\_{33}&\\mathcal L\_{34}\\\\0&0&0&0&\\mathcal L\_{41}&\\mathcal L\_{42}&\\mathcal L\_{43}&\\mathcal L\_{44}\\end{pmatrix}X ,\\tag{5}\\\]

eight real equations, the operator \\(8\\times8\\), block diagonal: the original four and the added four do not exchange anything. Each \\(\\mathcal L\_{ij}\\) is the entry of (4), a real first-order operator, \\(\\mathcal L=-cA_x\\partial_x-cA_y\\partial_y-cA_z\\partial_z+aT\\) with

\\\[A_x=\\begin{pmatrix}0&1&0&0\\\\1&0&0&0\\\\0&0&0&-1\\\\0&0&-1&0\\end{pmatrix},\\quad A_y=\\begin{pmatrix}0&0&1&0\\\\0&0&0&1\\\\1&0&0&0\\\\0&1&0&0\\end{pmatrix},\\quad A_z=\\begin{pmatrix}1&0&0&0\\\\0&-1&0&0\\\\0&0&-1&0\\\\0&0&0&1\\end{pmatrix},\\quad T=\\begin{pmatrix}0&0&1&0\\\\0&0&0&1\\\\-1&0&0&0\\\\0&-1&0&0\\end{pmatrix}.\\\]

\\(\\mathcal L\\) is real and, with the integration by parts, antisymmetric, so \\(\\sum u^2\\) and \\(\\sum v^2\\) are each conserved (III.8). Define the real matrix

\\\[J=\\begin{pmatrix}0&-1_4\\\\1_4&0\\end{pmatrix}=\\begin{pmatrix}0&0&0&0&-1&0&0&0\\\\0&0&0&0&0&-1&0&0\\\\0&0&0&0&0&0&-1&0\\\\0&0&0&0&0&0&0&-1\\\\1&0&0&0&0&0&0&0\\\\0&1&0&0&0&0&0&0\\\\0&0&1&0&0&0&0&0\\\\0&0&0&1&0&0&0&0\\end{pmatrix},\\qquad J\\begin{pmatrix}u_1\\\\u_2\\\\u_3\\\\u_4\\\\v_1\\\\v_2\\\\v_3\\\\v_4\\end{pmatrix}=\\begin{pmatrix}-v_1\\\\-v_2\\\\-v_3\\\\-v_4\\\\u_1\\\\u_2\\\\u_3\\\\u_4\\end{pmatrix},\\qquad J^2=-1_8 ,\\tag{6}\\\]

the walker’s own turn: each \\(u_j\\) onto \\(v_j\\) and each \\(v_j\\) onto \\(-u_j\\), the time sense untouched; every entry \\(0\\) or \\(\\pm1\\). With \\(\\mathcal A_k=\\mathrm{diag}(A_k,A_k)\\) and \\(\\mathcal T=\\mathrm{diag}(T,T)\\) — for instance

\\\[\\mathcal T=\\begin{pmatrix}T&0\\\\0&T\\end{pmatrix}=\\begin{pmatrix}0&0&1&0&0&0&0&0\\\\0&0&0&1&0&0&0&0\\\\-1&0&0&0&0&0&0&0\\\\0&-1&0&0&0&0&0&0\\\\0&0&0&0&0&0&1&0\\\\0&0&0&0&0&0&0&1\\\\0&0&0&0&-1&0&0&0\\\\0&0&0&0&0&-1&0&0\\end{pmatrix},\\qquad J\\mathcal T=\\mathcal TJ=\\begin{pmatrix}0&-T\\\\T&0\\end{pmatrix}\\\]

— the two quarter-turns side by side: \\(\\mathcal T\\) turns entries \\(1\\leftrightarrow3\\), \\(2\\leftrightarrow4\\) within each half, \\(J\\) turns the upper half into the lower —

\\\[J\\mathcal A_k=\\mathcal A_kJ,\\qquad J\\mathcal T=\\mathcal TJ ,\\tag{7}\\\]

\\(J\\) commutes with every walker matrix, because it acts between the two sets of four and they act within each.

8. on a \\(7^3\\) periodic lattice with random real \\(u,v\\) the \\(8\\times8\\) real evolution and the complex evolution of \\(z=u+iv\\) agree to \\(3\\times10^{-14}\\) after \\(t=1.3\\), and \\(\\sum(u^2+v^2)\\) is conserved to \\(2\\times10^{-15}\\) (`eight_real_counters.py` 1–3a). Squaring: \\(\\mathcal A_k^2=1\\), \\(\\mathcal T^2=-1\\), all pairs anticommute, so each of the eight obeys \\(\\partial_t^2X_j=c^2\\nabla^2X_j-a^2X_j\\): IV.4a’s two demands hold on the eight.

**(C) What the lines do not force.** In (4) the lines of \\(u\\) contain only \\(u\\) and the lines of \\(v\\) only \\(v\\):

9. the original four are a closed set and so are the added four. The free lines never move an amount from one set to the other. A walker put into the original four stays there, and nothing in (4) says the added four are needed. They are taken, for three reasons: rest ((D)), the two numbers that say how far below \\(c\\) the walker moves (IV.5b (A)), and charge ((E), V.2).
10. with nothing depending on \\(x\\) and \\(y\\) the eight are four two-counter walks of III.5a (H) side by side: \\((u_1,u_3)\\), \\((u_4,u_2)\\), \\((v_1,v_3)\\), \\((v_4,v_2)\\), each a forward detector with the backward detector whose orientation is turned in the walker’s own plane (`eight_detectors.py`, check 11).

**(D) Rest.** The record and the flow of the eight are the sums over the two sets, \\(\\rho=u^{\\mathsf T}u+v^{\\mathsf T}v\\), \\(j_k=c\\,(u^{\\mathsf T}A_ku+v^{\\mathsf T}A_kv)\\). Exchanging the forward and the backward detector of every orientation is one matrix, the product of the two quarter-turns of (B):

\\\[J\\,\\mathcal T=\\begin{pmatrix}0&-T\\\\T&0\\end{pmatrix},\\qquad J\\mathcal T\\begin{pmatrix}u\\\\v\\end{pmatrix}=\\begin{pmatrix}-Tv\\\\Tu\\end{pmatrix},\\qquad(J\\mathcal T)^2=1.\\tag{8}\\\]

11. it carries \\(u_1\\) onto \\(v_3\\) and back, \\(v_1\\) onto \\(u_3\\), \\(v_4\\) onto \\(u_2\\), \\(u_4\\) onto \\(v_2\\), each with the signs of the table: the two detectors of one orientation (`detectors_in_master_names.py`, check 2). It is the time flip followed by the walker’s own turn, which undoes the quarter-turn the flip made in the orientation.
12. *With eight.* Let nothing depend on \\(x,y,z\\). Then \\(\\partial_tX=a\\,\\mathcal T\\,X\\). The walker with the same amount in the forward and the backward detector of every orientation is \\(J\\mathcal T\\,X=X\\), that is \\(v=Tu\\). For it \\(j_k=0\\) at all times, and it stays so; \\(J\\mathcal T\\,X=-X\\), \\(v=-Tu\\), the same. These are the only two ways: writing a walker at rest as first set \\(u+x\\), second set \\(T(u-x)\\), the flow vanishes at all times only if \\(x^{\\mathsf T}A_ku=0\\) and \\(x^{\\mathsf T}TA_ku=0\\) for \\(k=x,y,z\\), and the six vectors \\(A_ku\\), \\(TA_ku\\) fill the four counters whenever \\(u\\neq0\\) (the determinant of their table of scalar products is \\(4(u^{\\mathsf T}u)^4\\)): \\(x=0\\), or \\(u=0\\) (`walker_at_rest_needs_the_pair.py`, check 5).
13. *With the original four.* No orientation has both a forward and a backward detector there. Example: \\(u(t)=(\\cos at,\\,0,\\,-\\sin at,\\,0)^{\\mathsf T}\\) solves (4) at rest, with \\(\\rho=1\\) and

\\\[j_x=0,\\qquad j_y=-c\\,\\sin2at,\\qquad j_z=c\\,\\cos2at:\\\]

14. speed \\(c\\) at every moment, going round in the plane \\((y,z)\\) at the rate \\(2a\\) (IV.8 (D): the walker going round its axis at \\(c\\) with radius \\(\\hbar/2mc\\)). For any walker in the original four alone, \\(j_x^2+j_y^2+j_z^2=c^2\\rho^2\\) identically (`walker_at_rest_needs_the_pair.py`, check 2). At rest only on the average over a turn.
15. *The same walker in all eight:* \\(u=(\\cos at,0,-\\sin at,0)^{\\mathsf T}\\), \\(v=Tu=(-\\sin at,0,-\\cos at,0)^{\\mathsf T}\\). The added four carry the flow \\(c\\,(0,+\\sin2at,-\\cos2at)\\), the opposite one at every moment. One walker, \\(\\rho=2\\), at rest. By (4) its two sets turn as \\(\\cos at+\\sin at\\,T\\); said on the eight, \\(\\partial_tX=a\\,\\mathcal T\\,X=-a\\,J\\,X\\) when \\(J\\mathcal T\\,X=X\\): the walker turns in its own plane at the rate \\(a\\). This is its clock. For \\(J\\mathcal T\\,X=-X\\) it turns the other way.

**(E) The \\(i\\) is the walker’s own turn.** Multiply (5) by \\(\\hbar J\\):

\\\[\\hbar\\,J\\,\\partial_tX=\\Big\[-c\\hbar\\sum_kJ\\mathcal A_k\\,\\partial_k+\\hbar a\\,J\\mathcal T\\Big\]X ,\\tag{9}\\\]

16. every matrix and every field real. Now abbreviate: write \\(z=u+iv\\) for the pair. Then \\(JX\\) is written \\(-v+iu=iz\\), and (9) *is* IV.5’s equation, \\(i\\hbar\\,\\partial_tz=\\big\[c\\sum_kA_k(-i\\hbar\\partial_k)+i\\hbar a\\,T\\big\]z\\), letter for letter: the \\(i\\) in front of \\(\\hbar\\partial_t\\) is \\(J\\), the \\(i\\) in \\(\\hat p=-i\\hbar\\partial\\) is \\(J\\), the \\(i\\) in \\(i\\hbar aT\\) is \\(J\\).

**(F) What \\(J\\) is not.** It is not the mass turn. \\(T\\) turns the counters *within* each set of four — the two time senses into each other, at the rate \\(a\\) — and anticommutes with the \\(A_k\\); \\(J\\) turns the walker in its own plane and leaves its time sense alone: it carries each of the original four detectors onto one of the added four, and commutes with everything. That is why \\(J\\) can stand in for a scalar \\(i\\) and \\(T\\) cannot: the original four counters have no real matrix on them, other than multiples of \\(1\\), that commutes with all four of \\(A_x,A_y,A_z,T\\) (III.8b), so it cannot be read as “two complex components”. The \\(i\\) needs the added four detectors. The two turns commute with each other, \\(J\\mathcal T=\\mathcal TJ\\): the walker’s own turn and the time flip are different turns (III.8b again). Two rotations, kept apart by what they change. Turning the counters of one set, \\(q\\to e^{\\alpha T}q\\), changes its current: from \\(q=(1,0)^{\\mathsf T}\\) the current \\(c(q\_+^2-q\_-^2)\\) becomes \\(c\\cos2\\alpha\\) — that rotation is the mass flip itself. Turning the walker in its own plane, \\(X\\to e^{\\alpha J}X\\), changes nothing a detector reads, \\(\\rho\\) and \\(j\\) alike, and shifts the momentum density by \\(\\hbar\\rho\\,\\partial_x\\alpha\\) (VI.2 (B), (D)). The original four admit no rotation of the second kind; the eight admit exactly one.

**(G) The slice and the general state.** With \\(v=0\\) at \\(t=0\\), \\(v=0\\) for all \\(t\\): the real four-counter solution is a closed slice of (5), and the complex evolution of real initial data stays real (`eight_real_counters.py` 3b). It is the walker in the original four detectors alone, (C) 9: it moves at \\(c\\) at every moment and no field can turn it. It is not the general state: a general column of IV.5 is eight real functions, and a real column is the special case in which the second is zero — or the walker launched at rest, \\(v(0)=Tu(0)\\), (D) 12, of one sign of energy when the source is wider than a Compton length (III.8c (E)). A change of complex basis in IV.5 (\\(U=\\mathrm{diag}(1,1,i,i)\\), the chiral map of IV.4) is a real \\(8\\times8\\) rotation that mixes \\(u\\) and \\(v\\) componentwise; it need not carry a real slice to a real slice, which is what IV.5’s sentence about the chiral basis is about.

**(H) The phase rotates each pair of counters.** Since \\(u\\) and \\(v\\) are two real fields on the same lattice, the “phase” of the complex column is a real thing at each point — but not one angle between the two four-vectors. Counterexample: \\(z=e^{i\\theta}(1,i)/\\sqrt2\\) has \\(u\\cdot v=0\\) for every \\(\\theta\\); the two real vectors stay at \\(90^\\circ\\) while the phase turns. What the phase does is rotate every pair \\((u_j,v_j)\\) by the same angle: \\(\\varphi_j=\\operatorname{atan2}(v_j,u_j)\\) is \\(\\theta+\\arg\\chi_j\\) component by component, and a common phase is the real rotation

\\\[e^{\\theta J}=\\begin{pmatrix}\\cos\\theta\\,1_4&-\\sin\\theta\\,1_4\\\\\\sin\\theta\\,1_4&\\cos\\theta\\,1_4\\end{pmatrix},\\qquad z=\\sqrt\\rho\\,e^{i\\theta}\\chi,\\ \\chi=r+is:\\quad u=\\sqrt\\rho\\,(r\\cos\\theta-s\\sin\\theta),\\ \\ v=\\sqrt\\rho\\,(r\\sin\\theta+s\\cos\\theta),\\\]

with \\(r,s\\) the real and imaginary parts of the axis column (IV.3’s columns are complex; \\(u=\\sqrt\\rho\\cos\\theta\\,\\chi\\) is the case \\(s=0\\)). The momentum has an exact real-field form that needs no single angle:

\\\[\\mathcal P_x=\\hbar\\sum\_{j=1}^4\\big(u_j\\,\\partial_xv_j-v_j\\,\\partial_xu_j\\big)=\\mathrm{Re}\\big\[z^\\dagger(-i\\hbar\\partial_x)z\\big\]=\\rho\\,\\partial_xS+\\hbar\\rho\\,\\mathrm{Im}\\big(\\chi^\\dagger\\partial_x\\chi\\big),\\\]

VI.2 (D)’s “the momentum density is the twist of the pair”, one pair per component: with \\(u_j=R_j\\cos\\varphi_j\\), \\(v_j=R_j\\sin\\varphi_j\\) each term is \\(R_j^2\\,\\partial_x\\varphi_j\\), the twist of that pair weighted by its density, and the last form says why a varying axis cannot be dropped — its own twist, \\(\\hbar\\,\\mathrm{Im}(\\chi^\\dagger\\partial_x\\chi)\\), is IV.6’s spin term. Integrated, \\(\\mathcal P_x\\) is \\(\\langle p\\rangle\\) (\\(0.30000\\) against \\(0.30000\\) for a \\(1{+}1\\) packet). So a walker moving in \\(+x\\) is not “\\(u\\) at \\(0^\\circ\\) and \\(v\\) at \\(90^\\circ\\)”: it is every pair’s angle *winding* along \\(+x\\); moving in \\(-x\\), winding the other way, and \\(\\mathcal P_x\\) changes sign. Each pair’s angle turns in \\(t\\) at \\(E/\\hbar=a\_{\\text{eff}}\\) (V.4). For a plane wave of one sign of energy \\(v_j(x)=u_j(x-\\lambda/4)\\) for every component: the second walk is the first displaced by a quarter wavelength in \\(x\\) (and a quarter period in \\(t\\)). For a packet only \\(\\mathcal P_x\\) and the rate are the walker’s momentum and energy. (`phase_is_the_angle_between_two_walks.py`: 1 the counterexample and the per-pair angle; 2 \\(e^{\\theta J}\\) and the \\(r,s\\) form; 3 \\(\\mathcal P_x\\) against \\(\\mathrm{Re}\[z^\\dagger\\hat pz\]\\) site by site and against \\(\\rho\\partial_xS+\\hbar\\rho\\,\\mathrm{Im}(\\chi^\\dagger\\partial_x\\chi)\\); 4 the quarter wavelength, the sign with the direction, the rate \\(E/\\hbar\\).)

One consequence, for the vortex of IV.8 (E) and for input 4. The original four and the added four are two real sheets; each has its zero set; a vortex is where the two zero sets cross, and the circulation around it is the number of times the angle of the pairs \\((u_j,v_j)\\) turns on the way round. With the original four alone there are zeros but nothing winds. And because each component’s phase is the angle of a pair of single-valued real fields, it is defined only modulo \\(2\\pi\\), and its winding around any loop that avoids that pair’s zeros is \\(2\\pi n\\) by itself: Wallstrom’s condition (input 4) holds for each pair and is not an extra condition; the circulation of the whole current adds the axis term of IV.8 (B), \\((h/m)\\,k\\sin^2\\tfrac{\\theta_0}2\\). In Part I, with one density and a real \\(S\\) from Hamilton–Jacobi, it is an input; from III.8c on, where the walker is read by all its detectors (input 11), it follows.

*What this changes.* Nothing in the equations; one thing in the words. “The quantum step is the classical step with \\(i\\)” (III.9) now reads: the quantum step is the classical step on all the detectors, with the walker’s own turn \\(J\\) written as \\(i\\) — and \\(J\\) is what V.2’s field turns. (The measurement’s returns are a different thing: every return route that agrees on the endpoint, III.8c (A).)

*What went in.* IV.4b’s lines and its three rules, unchanged. One thing is taken: the detectors read the orientation and the time sense separately, eight instead of four (input 11). The lines do not force it ((C) 9). It is what rest and charge need.

(`scripts/detectors_in_master_names.py`: 1 the table; 2 the exchange \\(J\\mathcal T\\) and rest. `scripts/eight_detectors.py`: 2–4, 11 the rules on the eight, the two closed sets, the four walks along the clock. `scripts/walker_at_rest_needs_the_pair.py`: 2, 5. `scripts/eight_real_counters.py`: 0 the matrices; 1 (6)–(7); 2 the two evolutions agree; 3a conservation; 3b the slice; 4 \\(T\\) is not a complex structure. `scripts/phase_is_the_angle_between_two_walks.py`: (H).)

### IV.5b Sums and differences: the eight sorted by axis, the planes as operations, and the frame

*The picture first.* Add and subtract the forward and the backward detector of one orientation. The eight sums and differences sort by axis: one sum and one difference for each of \\(t,x,y,z\\). In that sorting three things become plain. The flow along \\(x\\) is made of the pairs \\((t,x)\\) and \\((y,z)\\), and the eight lines are a divergence, a gradient and a curl. The six planes of spacetime are *operations*: each moves the amount in every detector onto one other detector; the three rotations keep sums as sums, the three runs exchange sums and differences. And the eight amounts are a *frame*: the eight products of the three runs, applied to one detector, land on the eight, so each amount says how much the walker is turned in one way. For a walker at rest relative to the detectors the sums hold the walker whose clock turns one way and the differences the walker whose clock turns the other way; a walker in motion is in both, the differences being \\(\\tanh(\\eta/2)\\) times the sums.

*Symbols.* \\(s_t,s_x,s_y,s_z\\) and \\(d_t,d_x,d_y,d_z\\) are the sums and differences of (1); \\(\\mathbf s=(s_x,s_y,s_z)\\), \\(\\mathbf d=(d_x,d_y,d_z)\\); \\(n=(s_t,s_x,s_y,s_z\\,|\\,d_t,d_x,d_y,d_z)^{\\mathsf T}\\). \\(\\mathcal A_k\\), \\(\\mathcal T\\), \\(J\\) are IV.5a’s matrices, written here in the order of \\(n\\). \\(M\\) is the frame, a \\(4\\times4\\) matrix; \\(\\tilde M\\) its mirror image; \\(\\overline M\\) is \\(M\\) with its six plane turns reversed; \\(\\omega=A_xA_yA_z\\). \\(\\eta\\) is a rapidity and \\(\\mathbf e\\) a unit direction.

**(A) The eight sorted by axis.** For each orientation, half the sum and half the difference of the amounts in its forward and its backward detector (the table of IV.5a (B)):

\\\[\\begin{aligned}s_t&=\\tfrac12(u_1-v_3),&d_z&=\\tfrac12(u_1+v_3)&&\\text{(as launched: }+u_1,\\ -v_3),\\\\s_x&=\\tfrac12(v_1+u_3),&d_y&=-\\tfrac12(v_1-u_3)&&\\text{(turned in }(y,z)\\text{: }+v_1,\\ +u_3),\\\\s_y&=\\tfrac12(v_4-u_2),&d_x&=\\tfrac12(v_4+u_2)&&\\text{(turned in }(z,x)\\text{: }+v_4,\\ -u_2),\\\\s_z&=-\\tfrac12(u_4+v_2),&d_t&=\\tfrac12(v_2-u_4)&&\\text{(turned in }(x,y)\\text{: }-u_4,\\ -v_2).\\end{aligned}\\tag{1}\\\]

As a matrix on \\(X=(u_1,u_2,u_3,u_4\\,|\\,v_1,v_2,v_3,v_4)^{\\mathsf T}\\):

\\\[n={\\scriptsize\\left(\\begin{array}{rrrr|rrrr}\\tfrac12&0&0&0&0&0&-\\tfrac12&0\\\\0&0&\\tfrac12&0&\\tfrac12&0&0&0\\\\0&-\\tfrac12&0&0&0&0&0&\\tfrac12\\\\0&0&0&-\\tfrac12&0&-\\tfrac12&0&0\\\\\\hline 0&0&0&-\\tfrac12&0&\\tfrac12&0&0\\\\0&\\tfrac12&0&0&0&0&0&\\tfrac12\\\\0&0&\\tfrac12&0&-\\tfrac12&0&0&0\\\\\\tfrac12&0&0&0&0&0&\\tfrac12&0\\end{array}\\right)}\\,X.\\\]

1. the exchange \\(J\\mathcal T\\) of IV.5a (8) keeps the four sums and reverses the four differences: in this order it is \\(\\mathrm{diag}(1,1,1,1\\,|\\,-1,-1,-1,-1)\\). So the walker at rest of IV.5a (D) 12 has sums only (\\(v=Tu\\)) or differences only (\\(v=-Tu\\)). Why the axis names: putting (1) into the record and the flow,

\\\[\\rho=2\\,\\big(s_t^2+|\\mathbf s|^2+d_t^2+|\\mathbf d|^2\\big),\\qquad\\mathbf j=4c\\,\\big(s_t\\,\\mathbf d+d_t\\,\\mathbf s-\\mathbf s\\times\\mathbf d\\big).\\tag{2}\\\]

2. the flow along \\(x\\) is made by the \\(t\\)-detector of one kind with the \\(x\\)-detector of the other kind, and by the \\(y\\)- and \\(z\\)-detectors crossed. The sum \\(s_x\\) is the amount turned about the \\(x\\) axis; the difference \\(d_x\\) is what a walker as launched (\\(s_t\\) only) needs in order to flow along \\(x\\), \\(j_x=4c\\,s_td_x\\). The eight lines of IV.5a (4), added and subtracted by (1):

\\\[\\begin{aligned}\\partial_ts_t&=-c\\,\\partial_zd_z-c\\,\\partial_xd_x-c\\,\\partial_yd_y+a\\,s_x\\\\\\partial_ts_x&=c\\,\\partial_zd_y-c\\,\\partial_xd_t-c\\,\\partial_yd_z-a\\,s_t\\\\\\partial_ts_y&=-c\\,\\partial_zd_x+c\\,\\partial_xd_z-c\\,\\partial_yd_t+a\\,s_z\\\\\\partial_ts_z&=-c\\,\\partial_zd_t-c\\,\\partial_xd_y+c\\,\\partial_yd_x-a\\,s_y\\\\\\partial_td_t&=-c\\,\\partial_zs_z-c\\,\\partial_xs_x-c\\,\\partial_ys_y+a\\,d_x\\\\\\partial_td_x&=-c\\,\\partial_zs_y-c\\,\\partial_xs_t+c\\,\\partial_ys_z-a\\,d_t\\\\\\partial_td_y&=c\\,\\partial_zs_x-c\\,\\partial_xs_z-c\\,\\partial_ys_t-a\\,d_z\\\\\\partial_td_z&=-c\\,\\partial_zs_t+c\\,\\partial_xs_y-c\\,\\partial_ys_x+a\\,d_y\\end{aligned}\\tag{3}\\\]

3. with \\(\\nabla=(\\partial_x,\\partial_y,\\partial_z)\\) they are a divergence, a gradient and a curl:

\\\[\\begin{aligned}\\partial_ts_t&=-c\\,\\nabla\\!\\cdot\\mathbf d+a\\,s_x,&\\partial_t\\mathbf s&=-c\\,\\nabla d_t-c\\,\\nabla\\times\\mathbf d+a\\,(-s_t,\\;s_z,\\;-s_y),\\\\\\partial_td_t&=-c\\,\\nabla\\!\\cdot\\mathbf s+a\\,d_x,&\\partial_t\\mathbf d&=-c\\,\\nabla s_t+c\\,\\nabla\\times\\mathbf s+a\\,(-d_t,\\;-d_z,\\;d_y).\\end{aligned}\\tag{4}\\\]

4. (an observation, not used: with \\(a=0\\) and \\(s_t=d_t=0\\), (4) has the shape of Maxwell’s equations in empty space.) Two more numbers, besides the record and the flow:

\\\[\\begin{gathered}\\alpha=-u^{\\mathsf T}Tv=\\big(s_t^2+|\\mathbf s|^2\\big)-\\big(d_t^2+|\\mathbf d|^2\\big),\\\\\\gamma=u^{\\mathsf T}T\\omega\\,v=2\\,\\big(s_td_t-\\mathbf s\\cdot\\mathbf d\\big),\\qquad\\rho^2-\\frac{|\\mathbf j|^2}{c^2}=4\\,\\big(\\alpha^2+\\gamma^2\\big).\\end{gathered}\\tag{5}\\\]

5. \\(\\alpha\\) is the overlap of the forward and the backward detectors of the same orientation. The walker moves below \\(c\\) exactly as far as \\(\\alpha\\) and \\(\\gamma\\) are not zero. On the original four alone, or the added four alone, both are zero: a walker in four detectors moves at \\(c\\) (IV.5a (D) 14).

**(B) The six planes are operations.** In the order of \\(n\\) the three direction matrices \\(\\mathcal A_k=\\mathrm{diag}(A_k,A_k)\\) are

\\\[\\mathcal A_x={\\scriptsize\\left(\\begin{array}{rrrr|rrrr}0&0&0&0&0&1&0&0\\\\0&0&0&0&1&0&0&0\\\\0&0&0&0&0&0&0&-1\\\\0&0&0&0&0&0&1&0\\\\\\hline 0&1&0&0&0&0&0&0\\\\1&0&0&0&0&0&0&0\\\\0&0&0&1&0&0&0&0\\\\0&0&-1&0&0&0&0&0\\end{array}\\right)},\\qquad\\mathcal A_y={\\scriptsize\\left(\\begin{array}{rrrr|rrrr}0&0&0&0&0&0&1&0\\\\0&0&0&0&0&0&0&1\\\\0&0&0&0&1&0&0&0\\\\0&0&0&0&0&-1&0&0\\\\\\hline 0&0&1&0&0&0&0&0\\\\0&0&0&-1&0&0&0&0\\\\1&0&0&0&0&0&0&0\\\\0&1&0&0&0&0&0&0\\end{array}\\right)},\\\]

\\\[\\mathcal A_z={\\scriptsize\\left(\\begin{array}{rrrr|rrrr}0&0&0&0&0&0&0&1\\\\0&0&0&0&0&0&-1&0\\\\0&0&0&0&0&1&0&0\\\\0&0&0&0&1&0&0&0\\\\\\hline 0&0&0&1&0&0&0&0\\\\0&0&1&0&0&0&0&0\\\\0&-1&0&0&0&0&0&0\\\\1&0&0&0&0&0&0&0\\end{array}\\right)}.\\tag{6}\\\]

6. each has one entry \\(\\pm1\\) in every row and column: it moves the amount in every detector onto one other detector, with a sign. \\(\\mathcal A_k\\) is a run along \\(k\\): when the walker runs along \\(+k\\) with rapidity \\(\\eta\\) relative to the detectors they read \\(\\big(\\cosh\\tfrac\\eta2+\\sinh\\tfrac\\eta2\\,\\mathcal A_k\\big)n\\) in place of \\(n\\) (IV.4c (L)). The product of two different runs is a rotation:

\\\[\\mathcal A_x\\mathcal A_y={\\scriptsize\\left(\\begin{array}{rrrr|rrrr}0&0&0&-1&0&0&0&0\\\\0&0&1&0&0&0&0&0\\\\0&-1&0&0&0&0&0&0\\\\1&0&0&0&0&0&0&0\\\\\\hline 0&0&0&0&0&0&0&1\\\\0&0&0&0&0&0&1&0\\\\0&0&0&0&0&-1&0&0\\\\0&0&0&0&-1&0&0&0\\end{array}\\right)}=-\\mathcal A_y\\mathcal A_x,\\\]

7. the rotation in \\((x,y)\\), with square \\(-1\\): IV.4c (G), two runs in different directions leave a turn. The products of the three runs are eight operations and no more: the identity, three runs, three rotations, and the product of all three. Where each sends the amount of each detector (column: the detector the amount starts in; entry: where it lands, with its sign):

| operation | \\(s_t\\) | \\(s_x\\) | \\(s_y\\) | \\(s_z\\) | \\(d_t\\) | \\(d_x\\) | \\(d_y\\) | \\(d_z\\) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rotation in \\((x,y)\\) | \\(s_z\\) | \\(-s_y\\) | \\(s_x\\) | \\(-s_t\\) | \\(-d_z\\) | \\(-d_y\\) | \\(d_x\\) | \\(d_t\\) |
| rotation in \\((y,z)\\) | \\(s_x\\) | \\(-s_t\\) | \\(-s_z\\) | \\(s_y\\) | \\(-d_x\\) | \\(d_t\\) | \\(-d_z\\) | \\(d_y\\) |
| rotation in \\((z,x)\\) | \\(s_y\\) | \\(s_z\\) | \\(-s_t\\) | \\(-s_x\\) | \\(-d_y\\) | \\(d_z\\) | \\(d_t\\) | \\(-d_x\\) |
| run along \\(x\\) | \\(d_x\\) | \\(d_t\\) | \\(-d_z\\) | \\(d_y\\) | \\(s_x\\) | \\(s_t\\) | \\(s_z\\) | \\(-s_y\\) |
| run along \\(y\\) | \\(d_y\\) | \\(d_z\\) | \\(d_t\\) | \\(-d_x\\) | \\(s_y\\) | \\(-s_z\\) | \\(s_t\\) | \\(s_x\\) |
| run along \\(z\\) | \\(d_z\\) | \\(-d_y\\) | \\(d_x\\) | \\(d_t\\) | \\(s_z\\) | \\(s_y\\) | \\(-s_x\\) | \\(s_t\\) |
| all three runs | \\(d_t\\) | \\(-d_x\\) | \\(-d_y\\) | \\(-d_z\\) | \\(-s_t\\) | \\(s_x\\) | \\(s_y\\) | \\(s_z\\) |
| quarter turn in the walker’s own plane | \\(s_x\\) | \\(-s_t\\) | \\(s_z\\) | \\(-s_y\\) | \\(-d_x\\) | \\(d_t\\) | \\(d_z\\) | \\(-d_y\\) |
| time flip | \\(-s_x\\) | \\(s_t\\) | \\(-s_z\\) | \\(s_y\\) | \\(-d_x\\) | \\(d_t\\) | \\(d_z\\) | \\(-d_y\\) |

8. the three rotations keep sums as sums and differences as differences; the three runs and the product of all three exchange them. A rotation in \\((x,y)\\) turns \\(x\\) into \\(y\\) and, with it, \\(t\\) into \\(z\\); a run along \\(x\\) exchanges \\(t\\) with \\(x\\) and \\(y\\) with \\(z\\). The last two rows are not among the eight: the walker’s own turn \\(J\\) acts on the other side and commutes with all eight, and the time flip is \\(\\mathcal T=-J\\,(J\\mathcal T)\\), the walker’s own turn after the exchange of forward and backward.

**(C) The walker in motion: a clock at rest, read by detectors that run.** In the order of \\(n\\) the lines are \\(\\partial_tn=-c\\,(\\mathcal A_x\\partial_x+\\mathcal A_y\\partial_y+\\mathcal A_z\\partial_z)\\,n+a\\,\\mathcal T\\,n\\), with \\(J\\mathcal A_k=\\mathcal A_kJ\\), \\((J\\mathcal T)\\mathcal A_k=-\\mathcal A_k(J\\mathcal T)\\), \\(J(J\\mathcal T)=(J\\mathcal T)J\\). Take the walker at rest with sums only, \\(J\\mathcal T\\,n_0=n_0\\): by IV.5a (D) 15 it is a clock, \\(n(t)=(\\cos at-\\sin at\\,J)\\,n_0\\). Let it run with rapidity \\(\\eta\\) along \\(\\mathbf e\\) relative to the detectors, \\(N=e_x\\mathcal A_x+e_y\\mathcal A_y+e_z\\mathcal A_z\\), \\(N^2=1\\). Its own time at the event \\((t,\\mathbf x)\\) is \\(\\tau=t\\cosh\\eta-(\\mathbf e\\cdot\\mathbf x/c)\\sinh\\eta\\) (III.8c (12)), and the detectors read

\\\[n(t,\\mathbf x)=\\Big(\\cosh\\tfrac\\eta2+\\sinh\\tfrac\\eta2\\,N\\Big)\\,\\big(\\cos a\\tau-\\sin a\\tau\\;J\\big)\\,n_0.\\tag{7}\\\]

9. differentiate: only \\(\\tau\\) depends on \\(t\\) and \\(\\mathbf x\\), and \\(J\\) commutes with \\(N\\),

\\\[\\partial_tn=-a\\cosh\\eta\\;J\\,n,\\qquad c\\,\\partial_kn=+a\\sinh\\eta\\;e_k\\,J\\,n,\\\]

\\\[\\partial_tn+c\\,\\mathcal A_k\\partial_kn=-a\\,J\\,(\\cosh\\eta-\\sinh\\eta\\,N)\\,n.\\\]

10. with \\(\\cosh\\eta-\\sinh\\eta\\,N=\\big(\\cosh\\tfrac\\eta2-\\sinh\\tfrac\\eta2N\\big)^2\\) and \\(\\big(\\cosh\\tfrac\\eta2-\\sinh\\tfrac\\eta2N\\big)\\big(\\cosh\\tfrac\\eta2+\\sinh\\tfrac\\eta2N\\big)=1\\):

\\\[\\begin{aligned}(\\cosh\\eta-\\sinh\\eta\\,N)\\,n&=\\Big(\\cosh\\tfrac\\eta2-\\sinh\\tfrac\\eta2\\,N\\Big)\\big(\\cos a\\tau-\\sin a\\tau\\,J\\big)\\,n_0\\\\&=\\Big(\\cosh\\tfrac\\eta2-\\sinh\\tfrac\\eta2\\,N\\Big)\\,(J\\mathcal T)\\,\\big(\\cos a\\tau-\\sin a\\tau\\,J\\big)\\,n_0\\\\&=(J\\mathcal T)\\,\\Big(\\cosh\\tfrac\\eta2+\\sinh\\tfrac\\eta2\\,N\\Big)\\big(\\cos a\\tau-\\sin a\\tau\\,J\\big)\\,n_0=(J\\mathcal T)\\,n,\\end{aligned}\\\]

11. the middle step because \\(J\\mathcal T\\,n_0=n_0\\) and \\(J\\mathcal T\\) commutes with \\(J\\), the last because it anticommutes with \\(N\\). Hence \\(\\partial_tn+c\\,\\mathcal A_k\\partial_kn=-a\\,J\\,(J\\mathcal T)\\,n=a\\,\\mathcal T\\,n\\): the lines. *The eight lines are what a clock at rest looks like to detectors that run.* The walker with differences only at rest goes the same way with its clock turning the other way. From (7) at \\(\\tau=0\\), \\(n=\\cosh\\tfrac\\eta2\\,n_0+\\sinh\\tfrac\\eta2\\,N\\,n_0\\), sums and then differences:

\\\[\\text{differences}=\\tanh\\tfrac\\eta2\\;\\times\\;N\\,(\\text{sums}),\\qquad\\frac{|\\mathbf j|}{c\\,\\rho}=\\tanh\\eta .\\tag{8}\\\]

12. a run does not reverse a clock: it moves amounts between the sums and the differences. For a walker that is one clock, however it is run or rotated, the sign of \\(\\alpha\\) says which clock it is: \\(\\alpha>0\\) if at rest it has sums only (the leading loop), \\(\\alpha\<0\\) if differences only (the lagging loop). A packet is not one clock, and for it \\(\\alpha\\) does not decide: a walker launched with \\(v=Tu\\) and a winding of wavenumber \\(k\\) has \\(\\alpha>0\\) everywhere and still has the share \\((\\omega_k-a)/2\\omega_k\\) of its weight in the other clock (III.15e (C)).

**(D) The frame: the eight as products of the three runs.** Read the column \\(s_t\\) of (B)’s table:

| product of runs | on the detectors | lands on | the plane | on IV.4b’s four counters |
| --- | --- | --- | --- | --- |
| none | \\(1\\) | \\(s_t\\) | not turned | \\(1\\) |
| \\(x\\) | \\(\\mathcal A_x\\) | \\(d_x\\) | time and \\(x\\) | \\(A_x\\) |
| \\(y\\) | \\(\\mathcal A_y\\) | \\(d_y\\) | time and \\(y\\) | \\(A_y\\) |
| \\(z\\) | \\(\\mathcal A_z\\) | \\(d_z\\) | time and \\(z\\) | \\(A_z\\) |
| \\(y\\!\\cdot\\!z\\) | \\(\\mathcal A_y\\mathcal A_z\\) | \\(s_x\\) | \\((y,z)\\) | \\(R\_{yz}\\) |
| \\(z\\!\\cdot\\!x\\) | \\(\\mathcal A_z\\mathcal A_x\\) | \\(s_y\\) | \\((z,x)\\) | \\(R\_{zx}\\) |
| \\(x\\!\\cdot\\!y\\) | \\(\\mathcal A_x\\mathcal A_y\\) | \\(s_z\\) | \\((x,y)\\) | \\(R\_{xy}\\) |
| \\(x\\!\\cdot\\!y\\!\\cdot\\!z\\) | \\(\\mathcal A_x\\mathcal A_y\\mathcal A_z\\) | \\(d_t\\) | — | \\(\\omega\\) |

13. the eight products land on the eight detectors, each once. Give each amount its product:

\\\[M=s_t+s_x\\,R\_{yz}+s_y\\,R\_{zx}+s_z\\,R\_{xy}+d_x\\,A_x+d_y\\,A_y+d_z\\,A_z+d_t\\,\\omega,\\\]

\\\[\\omega=A_xA_yA_z=\\begin{pmatrix}0&0&0&1\\\\0&0&-1&0\\\\0&1&0&0\\\\-1&0&0&0\\end{pmatrix},\\tag{9}\\\]

\\\[M=\\begin{pmatrix}s_t+d_z&s_y+d_x&-s_x+d_y&s_z+d_t\\\\-s_y+d_x&s_t-d_z&s_z-d_t&s_x+d_y\\\\s_x+d_y&-s_z+d_t&s_t-d_z&s_y-d_x\\\\-s_z-d_t&-s_x+d_y&-s_y-d_x&s_t+d_z\\end{pmatrix}\\\]

\\\[\\phantom{M}=\\begin{pmatrix}u_1&v_4&-v_1&-u_4\\\\u_2&-v_3&-v_2&u_3\\\\u_3&v_2&-v_3&-u_2\\\\u_4&-v_1&-v_4&u_1\\end{pmatrix}.\\tag{10}\\\]

14. the second form is the first with (1) put in: the first column of \\(M\\) is \\(u\\) and its third column is \\(-v\\); the forward detectors fill rows \\(1\\) and \\(4\\), the backward ones rows \\(2\\) and \\(3\\). *The change of variable, step by step.* \\(\\tfrac12(1+A_z)=\\mathrm{diag}(1,0,0,1)\\) keeps rows \\(1,4\\) and \\(\\tfrac12(1-A_z)=\\mathrm{diag}(0,1,1,0)\\) keeps rows \\(2,3\\). With the forward and the backward row of the table of IV.5a (B),

\\\[\\begin{aligned}M&=\\tfrac12(1+A_z)\\big(u_1+v_1R\_{yz}+v_4R\_{zx}-u_4R\_{xy}\\big)\\\\&\\quad+\\tfrac12(1-A_z)\\big(-v_3+u_3R\_{yz}-u_2R\_{zx}-v_2R\_{xy}\\big),\\end{aligned}\\tag{11}\\\]

15. which is half the sum of the two brackets plus \\(A_z\\) times half their difference. The half-sum is the four sums on \\(1,R\_{yz},R\_{zx},R\_{xy}\\). In the other term three one-line products turn the four into the runs and \\(\\omega\\),

\\\[A_zR\_{yz}=A_zA_yA_z=-A_y,\\qquad A_zR\_{zx}=A_zA_zA_x=A_x,\\qquad A_zR\_{xy}=A_zA_xA_y=\\omega,\\\]

16. and the result is \\(d_z\\), \\(d_y\\) (with its minus), \\(d_x\\) and \\(d_t\\) of (1). *The rules on the frame.* Each is one multiplication of \\(M\\):

- a run, \\(\\mathcal A_k\\), is \\(A_kM\\) (on the left);
- the walker’s own turn, \\(J\\), is \\(M\\,R\_{yz}\\) (on the right);
- the exchange of forward and backward, \\(J\\mathcal T\\), is the mirror image \\(\\tilde M=TMT^{-1}=s_t+s_xR\_{yz}+s_yR\_{zx}+s_zR\_{xy}-d_xA_x-d_yA_y-d_zA_z-d_t\\,\\omega\\): every run and \\(\\omega\\) reversed, every rotation kept;
- the time flip, \\(\\mathcal T\\), is \\(-\\tilde M\\,R\_{yz}\\);
- counting the added four with the other sign, VI.7’s \\(\\mathbb K\\), is \\(M\\,A_z\\) (on the right).

Put into the lines:

\\\[\\boxed{\\;\\partial_tM=-c\\,\\big(A_x\\,\\partial_xM+A_y\\,\\partial_yM+A_z\\,\\partial_zM\\big)-a\\,\\tilde M\\,R\_{yz}\\;}\\tag{12}\\\]

17. moving is a run applied to the frame; a time flip replaces the frame by its mirror image and turns it in its own plane. For the reference at rest \\(M(t)=\\cos at-\\sin at\\,R\_{yz}\\), and in motion \\(M(t,\\mathbf x)=\\big(\\cosh\\tfrac\\eta2+\\sinh\\tfrac\\eta2\\,(e_xA_x+e_yA_y+e_zA_z)\\big)\\big(\\cos a\\tau-\\sin a\\tau\\,R\_{yz}\\big)\\). *What no run and no rotation changes:*

\\\[M\\,\\overline M=\\alpha+\\gamma\\,\\omega,\\qquad\\rho=\\tfrac12\\operatorname{Tr}\\big(M^{\\mathsf T}M\\big),\\tag{13}\\\]

18. with \\(\\overline M=s_t-s_xR\_{yz}-s_yR\_{zx}-s_zR\_{xy}-d_xA_x-d_yA_y-d_zA_z+d_t\\,\\omega\\) and the \\(\\alpha,\\gamma\\) of (5). Writing \\(\\alpha+\\gamma\\,\\omega=\\sigma(\\cos\\beta+\\sin\\beta\\,\\omega)\\), every frame with \\(\\sigma>0\\) is a size \\(\\sqrt\\sigma\\), an angle \\(\\beta\\) and a product of runs and rotations: *the eight numbers are six plane turns, a size and an angle.* For one clock \\(\\beta=0\\) or \\(\\pi\\).
19. *Where “product” is literal.* The reference run along \\(z\\), then \\(y\\), then \\(x\\), with \\(t_k=\\tanh\\tfrac{\\eta_k}2\\), is \\(M=\\big(\\cosh\\tfrac{\\eta_1}2+\\sinh\\tfrac{\\eta_1}2A_x\\big)\\big(\\cosh\\tfrac{\\eta_2}2+\\sinh\\tfrac{\\eta_2}2A_y\\big)\\big(\\cosh\\tfrac{\\eta_3}2+\\sinh\\tfrac{\\eta_3}2A_z\\big)\\), and multiplied out

\\\[\\frac{d_x}{s_t}=t_1,\\qquad\\frac{d_y}{s_t}=t_2,\\qquad\\frac{d_z}{s_t}=t_3,\\\]

\\\[\\frac{s_z}{s_t}=t_1t_2,\\qquad\\frac{s_x}{s_t}=t_2t_3,\\qquad\\frac{s_y}{s_t}=-t_1t_3,\\qquad\\frac{d_t}{s_t}=t_1t_2t_3 :\\\]

20. the amount turned in \\((x,y)\\) is the amount run along \\(x\\) times the amount run along \\(y\\). For a general walker the eight amounts are free and the products are products of operations, as in the table.

*What went in.* Nothing beyond IV.5a: (1) is a change of variable, (B) and (D) are products of IV.4b’s matrices, (C) is III.15c’s run applied to the clock at rest. The form (12) of the equation on such a frame is in the literature (Hestenes, *J. Math. Phys.* 8, 798 (1967)).

(`scripts/detectors_in_master_names.py`: 3 for (1), (2), (5); 4 for \\(\\mathbb K\\). `scripts/eight_detectors.py`: 8 for (B); 11 for (2), (5); 12 for (C) and the invariance of \\(\\alpha,\\gamma\\); 13 for (3)–(4); 14 for (D). `scripts/one_walker_one_frame.py`: 3, 4, 8 for (12), (13). `scripts/boosting_the_rest_walker.py`: (C) by another road.)

## IV.6 The spin density is a canonical momentum, and it obeys the rotation algebra

(IV.7a says what becomes of the matrices in the variables of this section: every one of them becomes the axis.)

Take one block and write it in Madelung form with its axis: \\(\\psi=\\sqrt\\rho\\,e^{iS/\\hbar}\\chi(\\hat{\\mathbf s})\\), i.e.

\\\[\\psi_1=\\sqrt\\rho\\,\\cos\\tfrac\\theta2\\,e^{iS/\\hbar},\\qquad\\psi_2=\\sqrt\\rho\\,\\sin\\tfrac\\theta2\\,e^{i(S+\\hbar\\phi)/\\hbar}.\\\]

Reverse-Madelung each component: \\(\\rho_1=\\rho\\cos^2\\frac\\theta2\\), \\(\\rho_2=\\rho\\sin^2\\frac\\theta2\\), and for the phases take III.11’s mean, \\(S_1=\\bar S-\\frac\\hbar2\\phi\\), \\(S_2=\\bar S+\\frac\\hbar2\\phi\\) (the column above has an overall phase \\(e^{-i\\phi/2}\\) relative to this, which changes no density and no bracket; with \\(S_1=S\\) itself one would get \\(\\{S,P\_\\phi\\}=-\\hbar/2\\neq0\\)). This is III.11’s change of variables with \\(\\varphi\\to\\phi\\), and its bracket computation applies unchanged: \\((\\rho,\\bar S)\\) is a canonical pair and so is \\((\\phi,P\_\\phi)\\) with

\\\[P\_\\phi=\\frac\\hbar2(\\rho_1-\\rho_2)=\\frac\\hbar2\\,\\rho\\cos\\theta=\\frac\\hbar2\\,\\rho\\,\\hat s_z .\\\]

**The momentum conjugate to the azimuth of the axis is \\(\\hbar/2\\) times the density times the axis component.** That is the spin density, \\(\\frac\\hbar2\\rho\\hat{\\mathbf s}\\), obtained as a canonical momentum. Its three components, written in the canonical variables (\\(\\cos\\theta=2P\_\\phi/\\hbar\\rho\\)),

\\\[s_z=P\_\\phi,\\qquad s_x=\\sqrt{(\\tfrac{\\hbar\\rho}2)^2-P\_\\phi^2}\\,\\cos\\phi,\\qquad s_y=\\sqrt{(\\tfrac{\\hbar\\rho}2)^2-P\_\\phi^2}\\,\\sin\\phi,\\\]

have the brackets, from \\(\\{\\phi(\\mathbf x),P\_\\phi(\\mathbf x')\\}=\\delta(\\mathbf x-\\mathbf x')\\) alone (the \\(\\delta\\) factors out of every line below and is dropped),

\\\[\\{s_x,s_y\\}=\\frac{\\partial s_x}{\\partial\\phi}\\frac{\\partial s_y}{\\partial P\_\\phi}-\\frac{\\partial s_x}{\\partial P\_\\phi}\\frac{\\partial s_y}{\\partial\\phi}=(-A\\sin\\phi)\\Big(-\\frac{P\_\\phi}A\\sin\\phi\\Big)-\\Big(-\\frac{P\_\\phi}A\\cos\\phi\\Big)(A\\cos\\phi)=P\_\\phi=s_z,\\\]

\\(A\\equiv\\sqrt{(\\hbar\\rho/2)^2-P\_\\phi^2}\\), and cyclically \\(\\{s_y,s_z\\}=s_x\\), \\(\\{s_z,s_x\\}=s_y\\). **The spin density obeys the rotation algebra**, as a Poisson algebra, with no operator anywhere: the flow it generates on \\(\\hat{\\mathbf s}\\) is the ordinary \\(3\\times3\\) rotation (III.8a, “twice”), and on the column it is the half-rate one. This is the angular momentum of the plane of circulation: the direction you would turn in, carried as a canonical momentum of size \\(\\hbar/2\\) per particle.

## IV.7 The fluid form: the momentum density carries \\(\\frac\\hbar2\\nabla\\times(\\rho\\hat{\\mathbf s})\\)

(The Pauli equation in fluid variables is Takabayasi’s and Bohm, Schiller and Tiomno’s, *Suppl. Nuovo Cimento* 1, 48 and 2, 1 (1955); the curl of the spin density in the current is Gordon’s decomposition (1928). What is this document’s is where the terms come from in the walk.)

Now the streaming term of one block, \\(\\mathrm{Re}\\big\[\\psi^\\dagger\\,c\\,M(\\hat{\\mathbf p})\\,\\psi\\big\]=c\\,\\mathrm{Re}\\big\[\\psi^\\dagger\\sigma_i(-i\\hbar\\partial_i)\\psi\\big\]\\), in the variables of IV.6. Differentiate \\(\\psi=\\sqrt\\rho\\,e^{iS/\\hbar}\\chi\\):

\\\[\\partial_i\\psi=\\psi\\Big(\\frac{\\partial_i\\rho}{2\\rho}+\\frac i\\hbar\\partial_iS\\Big)+\\sqrt\\rho\\,e^{iS/\\hbar}\\,\\partial_i\\chi,\\\]

so, with \\(\\chi^\\dagger\\sigma_i\\chi=\\hat s_i\\),

\\\[\\psi^\\dagger\\sigma_i(-i\\hbar\\partial_i)\\psi=-\\frac{i\\hbar}2\\,\\hat{\\mathbf s}\\cdot\\nabla\\rho+\\rho\\,\\hat{\\mathbf s}\\cdot\\nabla S-i\\hbar\\,\\rho\\,\\chi^\\dagger\\sigma_i\\partial_i\\chi .\\\]

Everything hangs on the last term, the derivative of the axis itself. Compute it for \\(\\chi=(\\cos\\frac\\theta2,e^{i\\phi}\\sin\\frac\\theta2)\\) (multiply out; `scripts/IV7_streaming_term_in_axis_variables.py` 1a, 1b, and 3a, 3b for the box below):

\\\[\\mathrm{Re}\\big\[\\chi^\\dagger\\sigma_i\\partial_i\\chi\\big\]=\\tfrac12\\nabla\\cdot\\hat{\\mathbf s},\\qquad\\mathrm{Im}\\big\[\\chi^\\dagger\\sigma_i\\partial_i\\chi\\big\]=\\sin^2\\tfrac\\theta2\\;\\hat{\\mathbf s}\\cdot\\nabla\\phi+\\tfrac12\\,\\hat{\\mathbf s}\\cdot(\\nabla\\times\\hat{\\mathbf s}).\\\]

Therefore the real part of the streaming term is

\\\[\\mathrm{Re}\\big\[\\psi^\\dagger\\sigma_i(-i\\hbar\\partial_i)\\psi\\big\]=\\rho\\,\\hat{\\mathbf s}\\cdot\\nabla S+\\hbar\\rho\\sin^2\\tfrac\\theta2\\;\\hat{\\mathbf s}\\cdot\\nabla\\phi+\\frac\\hbar2\\,\\rho\\,\\hat{\\mathbf s}\\cdot(\\nabla\\times\\hat{\\mathbf s}),\\\]

and the imaginary part, \\(-\\frac\\hbar2\\big\[\\hat{\\mathbf s}\\cdot\\nabla\\rho+\\rho\\nabla\\cdot\\hat{\\mathbf s}\\big\]=-\\frac\\hbar2\\nabla\\cdot(\\rho\\hat{\\mathbf s})\\), integrates to zero. Use \\(\\hat{\\mathbf s}\\cdot(\\nabla\\rho\\times\\hat{\\mathbf s})=0\\) to write \\(\\rho\\,\\hat{\\mathbf s}\\cdot(\\nabla\\times\\hat{\\mathbf s})=\\hat{\\mathbf s}\\cdot\\nabla\\times(\\rho\\hat{\\mathbf s})\\). Then

\\\[\\boxed{\\mathrm{Re}\\big\[\\psi^\\dagger\\,c\\,M(\\hat{\\mathbf p})\\,\\psi\\big\]=c\\,\\hat{\\mathbf s}\\cdot\\mathbf P,\\qquad\\mathbf P=\\rho\\big(\\nabla S+\\hbar\\sin^2\\tfrac\\theta2\\nabla\\phi\\big)+\\frac\\hbar2\\,\\nabla\\times(\\rho\\,\\hat{\\mathbf s})}\\\]

**Reading.** The block streams at speed \\(c\\) along its axis, carrying the momentum density \\(\\mathbf P\\). \\(\\mathbf P\\) has two parts. The first is \\(m\\) times the current velocity, with the phase of the second component (\\(\\hbar\\sin^2\\frac\\theta2\\nabla\\phi\\)) added to \\(\\nabla S\\) — the completed square of IV.8 (A). The second is **the curl of the spin density**, \\(\\frac\\hbar2\\nabla\\times(\\rho\\hat{\\mathbf s})\\): the circulation of the vortex, which is zero in the bulk of an irrotational flow but not where \\(\\rho\\) or \\(\\hat{\\mathbf s}\\) change. This is the term the vortex picture obtains by accounting for the change of density across the circulation plane (IV.8 (C), (E)); here its component along \\(\\hat{\\mathbf s}\\) comes out of \\(\\chi^\\dagger\\sigma_i\\partial_i\\chi\\), the derivative of the axis. The streaming term fixes only \\(\\hat{\\mathbf s}\\cdot\\mathbf P\\) — any vector perpendicular to \\(\\hat{\\mathbf s}\\) could be added to \\(\\mathbf P\\) without changing it, which is why the identity \\(\\hat{\\mathbf s}\\cdot(\\nabla\\rho\\times\\hat{\\mathbf s})=0\\) was needed — so the vector form in the box is fixed by the slow limit of the current, IV.8 (C), not by this section alone. The piece \\(\\frac{\\hbar c}2\\rho\\,\\hat{\\mathbf s}\\cdot(\\nabla\\times\\hat{\\mathbf s})\\), the twist of the axis field along itself, is what the three axis-streamings of IV.3, each read in its own basis, add up to.

The full Hamiltonian is the sum over the two blocks with opposite signs of \\(c\\) (the second block streams the other way), plus the mass coupling between them, \\(2mc^2\\sqrt{\\rho_R\\rho_L}\\,\\mathrm{Re}\\big\[e^{i(S_R-S_L)/\\hbar}\\chi_L^\\dagger\\chi_R\\big\]\\) with \\(|\\chi_L^\\dagger\\chi_R|^2=\\frac12(1+\\hat{\\mathbf s}\_L\\cdot\\hat{\\mathbf s}\_R)\\) — the two senses couple most when their axes agree:

\\\[H=\\int\\Big\\{c\\,\\hat{\\mathbf s}\_R\\cdot\\mathbf P_R-c\\,\\hat{\\mathbf s}\_L\\cdot\\mathbf P_L+2mc^2\\sqrt{\\rho_R\\rho_L}\\,\\mathrm{Re}\\big\[e^{i(S_R-S_L)/\\hbar}\\chi_L^\\dagger\\chi_R\\big\]\\Big\\}d^3x ,\\\]

whose Hamilton equations with the four Madelung brackets are IV.5 (checked in the earlier notes this section grew from; not re-run here), and whose current is \\(\\mathbf j=c(\\rho_R\\hat{\\mathbf s}\_R-\\rho_L\\hat{\\mathbf s}\_L)\\).

## IV.7a The same matrices three times: in the clock, in \\(t\\), and in the Madelung variables

*The picture first.* Three questions kept coming back in the writing of this Part and they are one question. Are the walk’s matrices the Pauli matrices, or the rotations a flatlander sees? Why does the Dirac form hold spin matrices while the process holds turns? And what becomes of all of them once the column is written in \\(\\rho\\), \\(S\\) and the axis? The answer is that the same four matrices are read three times, in three languages, and in the last one they stop being matrices at all.

*Symbols.* \\(\\vartheta\\) is an angle of rotation (the polar angle \\(\\theta\\) of the axis is taken, IV.6). Nothing else is new.

**(A) In the walk’s clock: turns.** IV.4b (F), with \\(z\\) the clock,

\\\[c\\,\\partial_zQ=-A_z\\,\\partial_tQ+a\\,(A_zT)\\,Q+c\\,(A_xA_z)\\,\\partial_xQ+c\\,(A_yA_z)\\,\\partial_yQ .\\\]

The matrix on \\(\\partial_x\\) is \\(A_xA_z\\) and the one on \\(\\partial_y\\) is \\(A_yA_z\\): the rotation matrices of the two planes that contain the clock axis, antisymmetric, squaring to \\(-1\\) — IV.4c’s \\(W_x\\), \\(W_y\\), the U-turn coins. The mass is the symmetric exchange \\(A_zT\\), a half-turn. In the process’s own clock, moving sideways is turning and flipping is swapping.

**(B) Read in \\(t\\): spin matrices.** Multiply once by \\(A_z\\) (IV.4b (J)):

\\\[A_z\\,(A_xA_z)=-A_x,\\qquad A_z\\,(A_yA_z)=-A_y,\\qquad A_z\\,(A_zT)=T,\\\]

\\\[\\partial_tQ=-c\\big(A_x\\partial_x+A_y\\partial_y+A_z\\partial_z\\big)Q+aTQ .\\\]

Antisymmetric in, symmetric out, and the other way for the mass. The sideways terms are now the half-turns \\(A_x,A_y\\) — the Pauli matrices, in their eigenbasis of direction — and the mass is the quarter-turn \\(T\\). **The Dirac form holds the half-turns, the process form the quarter-turns, and they are one multiplication apart.** Which you see depends on which coordinate is the clock. The third rotation, \\(A_xA_y\\), the plane perpendicular to the clock, is written in neither; it is what the other two leave behind, \\(W_yW_x=A_xA_y=-W_xW_y\\) (IV.4c (G), (K)), and it does its work there: the mirror cancellation, \\(g=2\\) (V.3), the precession (IV.10). The flatlander’s matrices of III.8b (C)3 are the \\(2\\times2\\) blocks of (B): \\(B_2=\\sigma_z\\), \\(B_1=\\sigma_x\\) the two directions, \\(B_3=-T\\) the one turn; in 2+1 that is all there is, and exponentiated with a real angle \\(B_1,B_2\\) give the runs of IV.4d, not tilts. The tilts are the two products of (A) that need the second pair of counters.

**(C) After Madelung: the axis.** Write one block as \\(\\psi=\\sqrt\\rho\\,e^{iS/\\hbar}\\chi\\), \\(\\chi=(\\cos\\tfrac\\theta2,\\,e^{i\\phi}\\sin\\tfrac\\theta2)\\), \\(\\hat{\\mathbf s}=\\chi^\\dagger\\boldsymbol\\sigma\\chi\\) (IV.6). A matrix enters the equations for \\(\\rho,S,\\theta,\\phi\\) only sandwiched, and for any \\(M=m_0+\\mathbf m\\cdot\\boldsymbol\\sigma\\)

\\\[\\psi^\\dagger M\\psi=\\rho\\,\\big(m_0+\\mathbf m\\cdot\\hat{\\mathbf s}\\big).\\\]

So a matrix does not become another matrix; it becomes a number at each point, and every Pauli matrix becomes a component of the axis: \\(\\sigma_x\\to\\rho\\hat s_x\\), \\(\\sigma_y\\to\\rho\\hat s_y\\), \\(\\sigma_z\\to\\rho\\hat s_z\\). The imaginary \\(\\sigma_y\\) becomes the real \\(\\rho\\sin\\theta\\sin\\phi\\): the \\(i\\) is gone. The column carries \\(\\theta/2\\) and the axis \\(\\theta\\): the half-angle is gone. Each of the three readings of a matrix lands on the same object:

- *as a direction:* the streaming term \\(c\\,\\mathrm{Re}\[\\psi^\\dagger\\boldsymbol\\sigma\\cdot\\hat{\\mathbf p}\\,\\psi\]=c\\,\\hat{\\mathbf s}\\cdot\\mathbf P\\) (IV.7) — the block moves at \\(c\\) along its axis;
- *as a generator:* \\(e^{-i\\vartheta\\sigma_k/2}\\chi\\) has axis \\(R_k(\\vartheta)\\,\\hat{\\mathbf s}\\), the full angle; and in the canonical variables the generator is the spin density \\(\\mathbf s=\\tfrac\\hbar2\\rho\\hat{\\mathbf s}\\), with \\(\\{\\phi,s_z\\}=1\\) and \\(\\{s_x,s_y\\}=s_z\\) (IV.6) — the commutators have become brackets, and \\(s_z\\) rotates the axis about \\(z\\) by shifting its azimuth;
- *as the mass:* \\(\\mathrm{Re}\[\\psi_L^\\dagger\\psi_R\]=\\sqrt{\\rho_L\\rho_R}\\,|\\chi_L^\\dagger\\chi_R|\\cos\\big(\\tfrac{S_R-S_L}\\hbar+\\arg\\chi_L^\\dagger\\chi_R\\big)\\) — the overlap’s own phase stays in the cosine; a local phase choice can move it into \\(S\\), together with the connection terms — with \\(|\\chi_L^\\dagger\\chi_R|^2=\\tfrac12(1+\\hat{\\mathbf s}\_L\\cdot\\hat{\\mathbf s}\_R)\\) (IV.7) — the turn between the senses, weighted by how well their axes agree.

**(D) What this settles.** The distinction of (A)–(B), half-turn against quarter-turn, direction against spin, collapses in (C) into one vector: \\(\\hat{\\mathbf s}\\) is at once the direction the block streams along and the direction its angular momentum points, and its components, as densities, generate its own rotations. That is the vortex of IV.0: an axis that is both the velocity and the spin, \\(\\tfrac\\hbar2\\) per particle, with circulation \\(\\tfrac\\hbar2\\nabla\\times(\\rho\\hat{\\mathbf s})\\) wherever \\(\\rho\\) or \\(\\hat{\\mathbf s}\\) change, and, for a walker at rest in the original four detectors alone, the walker going round it at \\(c\\) with radius \\(\\hbar/2mc\\) (`true_spin.py` 3b). The flatlander of III.8a and III.8b is the column: his real vector is an axis lying in his plane at twice his angle, his in-plane rotation turns it by the double angle, and his tilt lifts it out of the plane with the lifted part in the other pair of counters (`madelung_of_the_matrices.py` 5). The Madelung variables are the axis itself. The one is the half-angle, complex reading of the other; nothing in the walk is a coin in a medium.

*What went in.* Nothing: (A) and (B) are IV.4b’s lines, (C) is IV.6–7’s change of variables, each line a check already on the page.

(`scripts/madelung_of_the_matrices.py` (C), (D); `scripts/B_to_A_dictionary.py` (B); `scripts/four_lines_clock.py` (A)–(B).)

## IV.8 The fluid of the axis: the square, the circulation, the current, and the slow limit

*The picture first.* IV.6 made the axis a canonical variable and IV.7 found \\(\\tfrac\\hbar2\\nabla\\times(\\rho\\hat{\\mathbf s})\\) in the streaming term, but three things were still only pointed to: the energy of a bent axis, the circulation it carries, and the slow limit in which Part I comes back with the spin. Here they are derived on the page, each line a check, with Nelson’s velocity \\(\\mathbf v=\\nabla S/m\\) kept irrotational throughout, as I.6 assumed it. What comes out is the vortex of IV.0 with no medium: in an irrotational fluid the spin is a circulation of the cloud about its axis driven by the change of density — the osmotic velocity turned a quarter turn about the axis — and the vortex with a hole at its centre is a different object, the orbital one.

*Symbols.* \\(\\mathbf u=D\\nabla\\ln\\rho\\) is I.3’s osmotic velocity, \\(D=\\hbar/2m\\). \\(\\eta\\) is the small component of the column in the slow limit. Nothing else is new.

**(A) The completed square.** One block, \\(\\Psi=\\sqrt\\rho\\,e^{iS/\\hbar}\\chi\\) with \\(\\chi=(\\cos\\tfrac\\theta2,e^{i\\phi}\\sin\\tfrac\\theta2)\\) and \\(\\hat{\\mathbf s}=\\chi^\\dagger\\boldsymbol\\sigma\\chi\\). Two facts about the column alone:

\\\[\\chi^\\dagger\\nabla\\chi=i\\sin^2\\tfrac\\theta2\\,\\nabla\\phi,\\qquad|\\nabla\\chi|^2=\\tfrac14|\\nabla\\hat{\\mathbf s}|^2+\\sin^4\\tfrac\\theta2\\,(\\nabla\\phi)^2 .\\\]

The connection of the column is purely imaginary, and its kinetic density is a quarter of the axis’s plus a square. Expand \\(|\\nabla\\Psi|^2\\): the \\(\\sqrt\\rho\\)–\\(S\\) and \\(\\sqrt\\rho\\)–\\(\\chi\\) cross terms vanish, the \\(S\\)–\\(\\chi\\) cross term survives and completes a square with the \\(\\sin^4\\) piece:

\\\[\\frac{\\hbar^2}{2m}|\\nabla\\Psi|^2=\\frac{\\hbar^2}{2m}(\\nabla\\sqrt\\rho)^2+\\frac m2\\rho\\,\\big|\\mathbf v+\\mathbf v\_\\phi\\big|^2+\\frac{\\hbar^2}{8m}\\rho\\,|\\nabla\\hat{\\mathbf s}|^2,\\qquad\\mathbf v=\\frac{\\nabla S}m,\\quad\\mathbf v\_\\phi=\\frac\\hbar m\\sin^2\\tfrac\\theta2\\,\\nabla\\phi .\\\]

Three energies: Part I’s osmotic term, the kinetic energy of Nelson’s velocity plus the axis’s own phase velocity \\(\\mathbf v\_\\phi\\) (the phase of the second component, which the axis carries with it), and the energy of bending the axis, \\(\\tfrac{\\hbar^2}{8m}|\\nabla\\hat{\\mathbf s}|^2\\) per particle. The current of the column is \\(\\tfrac\\hbar m\\mathrm{Im}\[\\Psi^\\dagger\\nabla\\Psi\]=\\rho(\\mathbf v+\\mathbf v\_\\phi)\\). For a uniform axis \\(\\mathbf v\_\\phi=0\\) and this is I.6.

**(B) The circulation of the phase.** \\(\\mathbf v\\) is a gradient; \\(\\mathbf v\_\\phi\\) is not:

\\\[\\nabla\\times(\\mathbf v+\\mathbf v\_\\phi)=\\frac\\hbar{2m}\\,\\sin\\theta\\,\\nabla\\theta\\times\\nabla\\phi .\\\]

The flow of the phase is irrotational except where the axis twists. On a loop along which \\(\\theta=\\theta_0\\) and the azimuth winds \\(k\\) times while \\(S\\) winds \\(n\\) times,

\\\[\\oint(\\mathbf v+\\mathbf v\_\\phi)\\cdot d\\mathbf l=\\frac hm\\Big\[n+k\\sin^2\\tfrac{\\theta_0}2\\Big\]:\\\]

the circulation is quantized in \\(h/m\\) plus a fraction fixed by the tilt of the axis, one half for an axis in the plane. That is the half-quantum of the column: a full turn of an axis lying in the plane puts \\(\\pi\\hbar/m\\) of circulation on the phase (for a tilted axis, \\(\\sin^2\\tfrac{\\theta_0}2\\) of a quantum, as the formula says), because the column turned by half a turn (III.8a).

**(C) The slow limit of Dirac, with the axis.** III.11 did this in one dimension: the fast pair follows the slow one, and \\(c\\cos\\theta=v\\), \\(c\\varphi=u\\). The same step in three dimensions, in the basis where the mass is diagonal: write the column as a large part \\(\\phi=\\Psi\\) and a small part \\(\\eta\\); to leading order in \\(v/c\\) and \\(u/c\\) the small part follows the large one, \\(\\eta=\\boldsymbol\\sigma\\cdot\\hat{\\mathbf p}\\,\\phi/2mc\\) — the stationary point of the fast pair, as in III.11. Then:

*The current.* The Dirac current \\(c\\,\\Psi^\\dagger\\boldsymbol\\alpha\\Psi=2c\\,\\mathrm{Re}\[\\phi^\\dagger\\boldsymbol\\sigma\\eta\]\\) becomes, using \\(\\sigma_i\\sigma_j=\\delta\_{ij}+i\\epsilon\_{ijk}\\sigma_k\\) once,

\\\[\\mathbf j=\\rho\\,(\\mathbf v+\\mathbf v\_\\phi)+\\frac\\hbar{2m}\\nabla\\times(\\rho\\hat{\\mathbf s}),\\qquad\\frac\\hbar{2m\\rho}\\nabla\\times(\\rho\\hat{\\mathbf s})=\\mathbf u\\times\\hat{\\mathbf s}+D\\,\\nabla\\times\\hat{\\mathbf s}.\\\]

This is the vector form of IV.7’s \\(\\mathbf P=m\\mathbf j\\): the streaming term fixed only \\(\\hat{\\mathbf s}\\cdot\\mathbf P\\), and the slow limit fixes the rest. The curl term is the spin current, and for a uniform axis it is exactly

\\\[\\mathbf j=\\rho\\,\\mathbf v+\\rho\\,(\\mathbf u\\times\\hat{\\mathbf s}),\\qquad\\nabla\\times\\mathbf v=0,\\qquad\\nabla\\cdot\\big\[\\rho\\,(\\mathbf u\\times\\hat{\\mathbf s})\\big\]=0 .\\\]

Call \\(\\mathbf u\\times\\hat{\\mathbf s}\\) the *osmotic circulation*. It is built from the osmotic velocity and vanishes with it where the cloud is flat; it is perpendicular to \\(\\mathbf u\\) and to the axis; it carries no density, as a vortex in an irrotational fluid carries none; and it reverses under time reversal, since \\(\\mathbf u\\) is even and \\(\\hat{\\mathbf s}\\) odd — read backward in time it circulates the other way, which is IV.0’s sentence. So it enters the two drifts with one sign, \\(\\mathbf b\_\\pm=\\mathbf v+\\mathbf u\\times\\hat{\\mathbf s}\\pm\\mathbf u\\): a transport velocity made of the spreading. The density gradient, which in Part I only spread the cloud, here also circulates it around the axis. Nelson’s \\(\\mathbf v\\) stays irrotational.

*The current, coin by coin.* The current is exact before any limit: \\(j_k=c\\,\\psi^\\dagger A_k\\psi\\) is the net flow of the \\(k\\)-th coin, the weight moving \\(+c\\) along \\(k\\) minus the weight moving \\(-c\\), and IV.4d (B) wrote it in the walk’s words — along the clock axis \\(Q\_+^2-Q\_-^2\\), the imbalance of the two streams, and along \\(x\\) the cross term \\(2Q\_+Q\_-\\) of the forward and backward \\(z\\)-streams. Take the slow limit in the \\(2{+}1\\) walk, in III.9’s own basis, \\(\\hat H=c\\,\\sigma_z\\hat p_z+c\\,\\sigma_x\\hat p_x+i\\hbar a\\,T\\). The mass matrix is \\(iT=-\\sigma_y\\); its eigenvectors are

\\\[e\_+=\\frac1{\\sqrt2}\\begin{pmatrix}1\\\\-i\\end{pmatrix}\\ (+1,\\ \\text{positive energy at rest}),\\qquad e\_-=\\frac1{\\sqrt2}\\begin{pmatrix}1\\\\ i\\end{pmatrix}\\ (-1),\\\]

and the axis of \\(e\_+\\) is the normal, \\(e\_+^\\dagger\\boldsymbol\\sigma e\_+=(0,-1,0)\\): in \\(2{+}1\\) the spin can only point out of the plane of motion. Write \\(\\psi=\\phi\\,e\_++\\eta\\,e\_-\\) with \\(\\phi=\\sqrt\\rho\\,e^{iS/\\hbar}\\) slow; the streaming term is off-diagonal in this basis, with matrix element

\\\[e\_-^\\dagger(\\sigma_z\\hat p_z+\\sigma_x\\hat p_x)e\_+=\\hat p_z-i\\hat p_x\\qquad\\big(e\_-^\\dagger\\sigma_ze\_+=\\tfrac12\[1\\cdot1+(-i)(i)\]=1,\\ \\ e\_-^\\dagger\\sigma_xe\_+=\\tfrac12\[1\\cdot(-i)+(-i)\\cdot1\]=-i\\big),\\\]

so the small component follows the large one as in III.11, \\(\\eta=(\\hat p_z-i\\hat p_x)\\phi/2mc\\). The flows need \\(e\_+^\\dagger\\sigma_ze\_+=e\_+^\\dagger\\sigma_xe\_+=0\\) and \\(e\_+^\\dagger\\sigma_ze\_-=1\\), \\(e\_+^\\dagger\\sigma_xe\_-=i\\), so to first order in \\(\\eta\\)

\\\[j_z=c\\,\\psi^\\dagger\\sigma_z\\psi=2c\\,\\mathrm{Re}(\\phi^\*\\eta),\\qquad j_x=c\\,\\psi^\\dagger\\sigma_x\\psi=2c\\,\\mathrm{Re}(i\\phi^\*\\eta)=-2c\\,\\mathrm{Im}(\\phi^\*\\eta).\\\]

Now \\(\\phi^\*\\partial\\phi=\\tfrac12\\partial\\rho+\\tfrac i\\hbar\\rho\\,\\partial S\\) (I.9), so

\\\[2mc\\,\\phi^\*\\eta=\\phi^\*(-i\\hbar\\partial_z\\phi)-i\\,\\phi^\*(-i\\hbar\\partial_x\\phi)=\\big(\\rho\\,\\partial_zS-\\tfrac\\hbar2\\partial_x\\rho\\big)-i\\big(\\rho\\,\\partial_xS+\\tfrac\\hbar2\\partial_z\\rho\\big),\\\]

\\\[j_z=\\frac{\\rho\\,\\partial_zS-\\frac\\hbar2\\partial_x\\rho}m=\\rho\\,v_z-\\rho\\,u_x,\\qquad j_x=\\frac{\\rho\\,\\partial_xS+\\frac\\hbar2\\partial_z\\rho}m=\\rho\\,v_x+\\rho\\,u_z ,\\\]

which is \\(\\mathbf j=\\rho(\\mathbf v+\\mathbf u\\times\\hat{\\mathbf s})\\) with \\(\\hat{\\mathbf s}=-\\hat{\\mathbf y}\\). Each coin’s net flow is Nelson’s velocity along that coin plus the osmotic velocity of the *perpendicular* direction, turned a quarter about the axis. At rest \\(j_x=\\rho u_z\\) exactly; and since \\(j_x=2c\\,\\mathrm{Re}(\\psi\_+^\*\\psi\_-)\\) is the cross term of the two \\(z\\)-streams, this says that the relative phase of the two time senses — III.11’s \\(c\\varphi=u\\) — is read sideways by \\(A_x\\) as a flow (`spin_current_by_direction.py`). That is the spin’s flow decomposed by direction: three sortings of one column into streams (IV.4d (C)), and in each sorting the osmotic velocity of that walk feeds the flow across it.

*The mechanism, as a sum over arrivals.* Why a density gradient along \\(z\\) makes a flow along \\(x\\): a walker arrives from one side of a step along \\(z\\), is turned a quarter about the axis by the sideways U-turn of IV.4c, with the sign its time sense gives it ((D) there), and leaves along \\(x\\). Write the sum over the six arrival directions with the walk’s own weights — the coin’s \\(\\pm\\tfrac12\\) per side (IV.4c (E)), the step \\(\\lambda=\\hbar/2mc=c/2a\\) (the orbit radius of `true_spin.py` 3b — half the mean distance \\(c/a\\) between flips, not one flip’s step; the sum over arrivals is the mechanism read on the slow limit, (C)’s derivation is the slow limit itself), and the quarter-turn rule \\(\\hat{\\mathbf e}\_j\\to\\hat{\\mathbf s}\\times\\hat{\\mathbf e}\_j\\) — the sense that makes the circulation of (C) and (E) — i.e. the factor \\(\\epsilon\_{ikj}\\hat s_k=-\\epsilon\_{ijk}\\hat s_k\\) for leaving along \\(i\\):

\\\[j^{\\rm spin}\_i=-c\\sum\_{j}\\Big\[\\tfrac12\\rho(\\mathbf x-\\lambda\\hat{\\mathbf e}\_j)-\\tfrac12\\rho(\\mathbf x+\\lambda\\hat{\\mathbf e}\_j)\\Big\]\\epsilon\_{ijk}\\hat s_k .\\\]

Expand each side, \\(\\rho(\\mathbf x\\mp\\lambda\\hat{\\mathbf e}\_j)=\\rho\\mp\\lambda\\,\\partial_j\\rho+\\tfrac12\\lambda^2\\partial_j^2\\rho\\mp\\ldots\\): the even orders cancel between the two sides and the bracket is \\(-\\lambda\\,\\partial_j\\rho+O(\\lambda^3)\\), so

\\\[j^{\\rm spin}\_i=+c\\lambda\\,\\epsilon\_{ijk}\\,\\partial_j\\rho\\,\\hat s_k=+\\frac\\hbar{2m}\\big\[\\nabla\\rho\\times\\hat{\\mathbf s}\\big\]\_i=+\\frac\\hbar{2m}\\big\[\\nabla\\times(\\rho\\hat{\\mathbf s})\\big\]\_i\\qquad(c\\lambda=\\tfrac\\hbar{2m};\\ \\text{uniform axis}),\\\]

the spin current, with nothing free: the density difference between the two sides of one step, times the turn. A sphere of arrival directions instead of six gives the same: with \\(\\int\\frac{d\\Omega}{4\\pi}n_jn_l=\\tfrac13\\delta\_{jl}\\),

\\\[\\int\\frac{d\\Omega}{4\\pi}\\,3\\,(\\hat{\\mathbf n}\\cdot\\nabla\\rho)\\,(\\hat{\\mathbf n}\\times\\hat{\\mathbf s})\_i=3\\cdot\\tfrac13\\,\\epsilon\_{ijk}\\,\\partial_j\\rho\\,\\hat s_k=\\big\[\\nabla\\rho\\times\\hat{\\mathbf s}\\big\]\_i\\\]

(`kick_integral_done_right.py`). The older vortex notes reached this sum with a background flow for the handedness, a free step length and a factor one half put in by hand; here the handedness is the time sense, the step is the orbit radius, and the half is the coin’s.

*The two senses.* The imbalance and the relative phase of the two blocks of IV.4 are III.11’s, projected on the axis, each with a term from the axis field itself:

\\\[\\rho_R-\\rho_L=\\frac\\rho c\\Big\[\\hat{\\mathbf s}\\cdot(\\mathbf v+\\mathbf v\_\\phi)+\\frac\\hbar{2m}\\hat{\\mathbf s}\\cdot(\\nabla\\times\\hat{\\mathbf s})\\Big\],\\qquad\\frac{S_L-S_R}\\hbar=\\frac1c\\Big\[\\hat{\\mathbf s}\\cdot\\mathbf u+\\frac\\hbar{2m}\\nabla\\cdot\\hat{\\mathbf s}\\Big\].\\\]

The pair \\(\\mathbf v\\pm\\mathbf u\\) is the pair of time senses projected on the axis: the twist of the axis along itself goes with the imbalance, its divergence with the phase. The osmotic circulation is in neither, being perpendicular to the axis: the streaming along the axis does not see it, and only the current does.

*The energy.* To this order the energy is \\(mc^2\\rho+\\tfrac1{2m}|\\boldsymbol\\sigma\\cdot\\boldsymbol\\pi\\,\\phi|^2\\) with \\(\\boldsymbol\\pi=\\hat{\\mathbf p}-\\tfrac ec\\mathbf A\\), and

\\\[(\\boldsymbol\\sigma\\cdot\\boldsymbol\\pi)^2=\\boldsymbol\\pi^2-\\frac{e\\hbar}c\\,\\boldsymbol\\sigma\\cdot\\mathbf B :\\\]

the Pauli term with \\(g=2\\), from the square of the direction matrices — V.3’s statement, reached from the slow side. With (A) the energy in fluid variables is Part I’s Hamiltonian, plus the bending energy of the axis, plus the Zeeman term: the Pauli fluid, derived.

**(D) Three angular momenta, and the factor two.** Take a localized cloud with a uniform axis. The curl term of \\(\\mathbf P\\) carries, by \\(\\int\\mathbf r\\times(\\nabla\\times\\mathbf A)\\,d^3x=2\\int\\mathbf A\\,d^3x\\),

\\\[\\int\\mathbf r\\times\\frac\\hbar2\\nabla\\times(\\rho\\hat{\\mathbf s})\\,d^3x=\\hbar\\,\\hat{\\mathbf s}:\\\]

twice the spin density’s \\(\\tfrac\\hbar2\\) of IV.6. That is not an error; it is where \\(g=2\\) lives. The momentum of the *energy flow* — the symmetrized stress tensor of the column, which in this limit is half the canonical momentum plus half of \\(m\\) times the current — has for its curl part half the Gordon curl, \\(\\tfrac\\hbar4\\nabla\\times(\\rho\\hat{\\mathbf s})\\), and that carries exactly \\(\\tfrac\\hbar2\\hat{\\mathbf s}\\): the spin. The *charge* flows with the full current, and its moment is

\\\[\\frac e{2c}\\int\\mathbf r\\times\\mathbf j\_{\\rm spin}\\,d^3x=\\frac{e\\hbar}{2mc}\\hat{\\mathbf s}=\\frac e{mc}\\cdot\\frac\\hbar2\\hat{\\mathbf s}:\\\]

the Bohr magneton on a spin of one half, \\(g=2\\), as the ratio of the charge’s circulation to the energy’s. The third number is the walker itself: a walker at rest in the original four detectors alone goes round its axis at \\(c\\) with radius \\(\\hbar/2mc\\) (`true_spin.py` 3b), angular momentum \\(m\\cdot\\tfrac\\hbar{2mc}\\cdot c=\\tfrac\\hbar2\\), the energy flow’s number.

**(E) The vortex, worked.** A cloud with its axis along \\(z\\) and no orbital phase: \\(\\mathbf v=0\\) and the whole flow is the osmotic circulation. For a Gaussian cloud of width \\(a\\),

\\\[\\mathbf u\\times\\hat{\\mathbf s}=\\frac\\hbar{ma^2}\\,(-y,\\,x,\\,0)=\\omega\\,\\hat{\\mathbf z}\\times\\mathbf r,\\qquad\\omega=\\frac\\hbar{ma^2}:\\\]

a rigid rotation of the cloud about its axis, fastest where the density falls fastest, zero at the centre, where the cloud is fullest. On a circle of radius \\(r\\) the circulation is \\(2\\pi\\hbar r^2/ma^2\\), growing as a rigid body’s, and it reaches \\(h/m\\) at the cloud’s radius \\(r=a\\); the energy flow carries half of it, \\(h/2m\\) there — the half-quantum of (B), at the one radius the cloud has. The angular momentum per particle is \\(\\hbar\\) in the charge flow and \\(\\tfrac\\hbar2\\) in the energy flow, for *any* axisymmetric profile — the azimuthal speed is \\(-D\\,\\partial_r\\ln\\rho\\) and one integration by parts gives \\(2mD=\\hbar\\) whatever the width — and the moment is \\(e\\hbar/2mc\\).

*The vortex has the hole; the spin does not.* The spin’s flow needs no hole, because its driver is \\(\\nabla\\rho\\) and not a phase. The vortex with a hole at its centre is the other one: a phase winding \\(S=\\hbar\\varphi\_{\\rm az}\\), \\(\\mathbf v=\\hbar/mr\\), singular at the axis so that \\(\\rho\\sim r^2\\) there — the orbital vortex, angular momentum \\(\\hbar\\) per particle, independent of the density. The old vortex notes modelled the spin as this holed vortex with \\(\\ell=\\hbar/m\\), and that is why they kept needing a factor of one half. When both are present, with \\(\\rho\\sim r^2e^{-r^2/a^2}\\), the osmotic circulation itself becomes a \\(1/r\\) flow at the core, \\(w=-\\hbar/mr+\\hbar r/ma^2\\), of circulation \\(-h/m\\) and in the sense *opposite* to the outer rotation; with the orbital sense equal to the spin’s the two \\(1/r\\) parts cancel and the total current is regular at the core, and the angular momenta add to \\(l+s=\\tfrac32\\hbar\\) in the energy flow and \\(l+2s=2\\hbar\\) in the charge flow, the Landé combination. With the senses opposite the singularities double.

**(F) What this closes.** The axis \\(\\hat{\\mathbf s}\\) is: the direction the block streams along (IV.7); the direction its spin points, \\(\\tfrac\\hbar2\\) per particle (IV.6); the thing whose bending costs \\(\\tfrac{\\hbar^2}{8m}|\\nabla\\hat{\\mathbf s}|^2\\) (A); the carrier of a half-quantum of phase circulation per turn (B); the axis about which the density gradient circulates the cloud, with \\(g=2\\) as the ratio of charge flow to energy flow (C)–(E). That is the whole of the vortex picture, with no medium: the particle of the medium going round a line, the circulation \\(\\hbar/m\\) attached to it, the kicks that tilt its axis, are replaced by the walker’s own column, whose axis does all of the above because the column is the half-angle object of III.8a. In an irrotational fluid the spin is carried by the change of density, and the two senses of circulation are the two time senses, because the osmotic circulation reverses with time while the osmotic velocity does not.

*What went in.* (A), (B), (D), (E): identities of the column and one vector identity, checked. (C): the slow limit, with one step — the small component follows the large one to leading order in \\(v/c\\), the same step as III.11’s fast pair — in the standard Dirac representation. The symmetrized momentum of (D) is the symmetric stress tensor of the Dirac column, computed in the same limit (`IV8_axis_fluid.py` E2′); the reading of its curl part as the spin is Ohanian, *Am. J. Phys.* 54, 500 (1986). Extra inputs in three dimensions, unchanged: one speed \\(c\\) in every direction for every component (IV.1), and minimal coupling for the field.

(`scripts/IV8_axis_fluid.py`: 0, A1–A3, B1–B2, C1–C5, D, E1–E4, E2′. `scripts/spin_current_by_direction.py`, `scripts/kick_integral_done_right.py`: (C), coin by coin and the sum over arrivals. `scripts/vortex_worked.py` and `scripts/vortex_core.py`: (E).)

---

## IV.9 The square root in three dimensions

III.13 again, with four components. In \\(\\psi\\): \\(\\hat H^2=\\big(c\\,\\mathbf A\\!\\cdot\\!\\hat{\\mathbf p}+mc^2F\\big)^2=c^2\\hat p^2+m^2c^4\\) because all six pairs anticommute (IV.4; in the standard basis \\(F^2=1\\), in the original basis \\((i\\hbar aT)^2=+\\hbar^2a^2\\) because \\(T^2=-1\\)), so every one of the four components obeys Klein–Gordon. Covariantly, with \\(\\gamma^0=F\\) and \\(\\gamma^k=F A_k\\): \\((\\gamma^0)^2=1\\), \\((\\gamma^k)^2=FA_kFA_k=-F^2A_k^2=-1\\), all anticommuting, and \\((\\gamma\\!\\cdot\\!\\partial)^2=\\partial_0^2-\\nabla^2\\) — the same factorization, four half-turns instead of two. In \\((\\rho,S)\\): each of the four Madelung pairs \\((\\rho_a,S_a)\\), \\(a=1..4\\), obeys II.7 by itself, by the same argument as III.13 (each \\(\\psi_a\\) solves Klein–Gordon, and II.8’s identity turns that into II.7 for its own pair). What Dirac adds beyond four Klein–Gordon fluids is the first-order tie between them: within a block, the sideways motion as a transition between the two components of the axis (IV.3), carried in the fluid by \\(\\frac\\hbar2\\nabla\\times(\\rho\\hat{\\mathbf s})\\) (IV.7); between the blocks, the reversal of the time sense by the mass, carried by the Madelung product \\(\\psi_L^\\dagger\\psi_R\\) (IV.7). Klein–Gordon is the same particle with the axis and the time sense averaged out.

## IV.10 The propagators in three dimensions

*The picture first.* Part IV gave the Dirac equation in \\(3{+}1\\) (IV.5) and its square root structure (IV.9). What has not been written is what a point does — the kernels. The lesson of IV.2a is that there is no classical one-clock process behind the \\(3{+}1\\) column, so the two objects available are the two that do not need one: the *source* answer, which Part II supplies in any dimension (III.16), and the *amplitude*, which the wave equation supplies by the on-shell sum of III.15. Both are derived here from inputs already in the document; neither is a new process.

*The scalar source answer.* III.16, \\(3{+}1\\) line: \\(\\Delta(s)=\\frac{m}{4\\pi^2 s}K_1(ms)\\), \\(s\\) the Euclidean distance.

*The Dirac source answer: the factorization doing work.* IV.9 wrote Dirac’s operator as a square root of Klein–Gordon’s: \\((i\\gamma\\cdot\\partial+m)(i\\gamma\\cdot\\partial-m)=-(\\gamma\\cdot\\partial)^2-m^2=\\Box-m^2\\) in Part II’s signature (\\(\\Box=-\\partial_t^2+\\nabla^2\\); the gammas square to \\(+,-,-,-\\), IV.9, so \\((\\gamma\\cdot\\partial)^2=-\\Box\\)). So the Dirac source answer is the Klein–Gordon one with the other factor applied, \\(S=(i\\gamma\\cdot\\partial+m)\\Delta\\), and that needs only \\(\\partial\_\\mu\\Delta\\). With \\(\\partial\_\\mu s=x\_\\mu/s\\) and \\(\\frac{d}{ds}\\big\[\\frac{K_1(ms)}s\\big\]=-\\frac{m}{s}K_2(ms)\\) (`three_coins_symbolic.py` E),

\\\[\\partial\_\\mu\\Delta=-\\frac{m^2}{4\\pi^2}\\,\\frac{x\_\\mu}{s^2}\\,K_2(ms),\\qquad S=\\underbrace{\\frac{m^2}{4\\pi^2 s}K_1(ms)}\_{\\text{scalar channel}}\\;-\\;i\\,\\underbrace{\\frac{m^2}{4\\pi^2 s^2}K_2(ms)\\,\\gamma\\cdot x}\_{\\text{vector channel}} .\\\]

Two channels, as in one dimension (\\(K_0\\) and \\(K_1\\) there): a part that does not know the direction to the source and a part along \\(\\gamma\\cdot x\\) that does. The relative \\(i\\) between them is the quarter-turn of III.8, now between the mass half-turn \\(\\gamma^0\\) and the direction.

*The amplitude.* From IV.5, exactly as in III.15: \\(e^{-i\\hat Ht/\\hbar}=\\cos\\frac{Et}\\hbar-i\\frac{\\hat H}{E}\\sin\\frac{Et}\\hbar\\) with \\(\\hat H=c\\,\\boldsymbol\\alpha\\cdot\\hat{\\mathbf p}+\\beta mc^2\\). Every term is a derivative of one scalar function,

\\\[D(t,\\mathbf x)=\\int\\frac{d^3p}{(2\\pi\\hbar)^3}\\,e^{i\\mathbf p\\cdot\\mathbf x/\\hbar}\\,\\frac{\\sin(Et/\\hbar)}{E/\\hbar}=\\frac1{2\\pi}\\,\\mathrm{sgn}(t)\\Big\[\\delta(s^2)-\\frac{m}{2s}J_1(ms)\\,\\Theta(s^2)\\Big\],\\qquad s^2=c^2t^2-r^2\\\]

(units \\(c=\\hbar=1\\) in the bracket; the bulk by a radial integral at three points, to \\(1.5\\times10^{-3}\\), `three_coins_symbolic.py` H), and

\\\[K(t,\\mathbf x)=\\Big(\\partial_t-c\\,\\boldsymbol\\alpha\\cdot\\nabla-i\\,\\frac{mc^2}\\hbar\\,\\beta\\Big)D(t,\\mathbf x).\\\]

The \\(\\delta(s^2)\\) differentiated is the \\(\\delta'\\) of IV.2a — the sharp front with its two skins; the \\(J_1/s\\) bulk is the \\(3{+}1\\) version of the \\(J_0,J_1\\) of III.15, and \\(J_1(ms)/s\\) continued across the cone is \\(I_1(ms')/s'\\), the same function on both sides, as before. What is *not* available here is the first column of III.15’s table: a classical kernel outside the cone that this continues into. IV.2a says why.

*The process in three dimensions, from the rules already on the page.* IV.2a forbids one coin for the direction; it does not forbid three. Three rules, none new: along any one axis the column decides a coin of \\(\\pm c\\) by its own projectors \\(\\frac12(1\\pm\\alpha_j)\\) — III.1 applied to that axis; the same rule on all three axes, since no axis is preferred (IV.1); and the coins taken *in turn*, since IV.2a says a joint rule for three non-commuting projectors has no positive form. With III.8’s turn matrix for the mass, one step of size \\(\\epsilon\\) is forced: an \\(x\\)-coin, a \\(y\\)-coin, a \\(z\\)-coin, then \\(\\cos(mc^2\\epsilon/\\hbar)-i\\beta\\sin(mc^2\\epsilon/\\hbar)\\). Each sub-step is unitary, so the step is, and nothing blows up. One step moves the particle by \\((\\pm c\\epsilon,\\pm c\\epsilon,\\pm c\\epsilon)\\): every component \\(\\pm c\\) at every instant, \\(c\\) on average along the axis — a full tick moves it by \\(\\sqrt3\\,c\\epsilon\\), so the three sub-steps are not one trajectory at speed \\(c\\) within the tick; it is the continuum generator, not each lattice path, that carries the relativistic limit.

*The limit, taken as in I.2.* Write each coin as a shift, Taylor in the step, divide by \\(\\epsilon\\) (`scripts/three_coins_symbolic.py`, symbolic):

\\\[\\frac{\\psi(t+\\epsilon)-\\psi(t)}{\\epsilon}\\;\\longrightarrow\\;-c\\,\\boldsymbol\\alpha\\cdot\\nabla\\psi-\\frac{i\\,mc^2}{\\hbar}\\,\\beta\\,\\psi,\\qquad\\text{i.e.}\\qquad i\\hbar\\,\\partial_t\\psi=\\big(-i\\hbar c\\,\\boldsymbol\\alpha\\cdot\\nabla+\\beta mc^2\\big)\\psi ,\\\]

IV.5’s Dirac equation. The order in which the three coins are taken changes the step only at \\(O(\\epsilon^2)\\), and that change is a sum of commutators \\(\[\\alpha_jp_j,\\alpha_kp_k\]\\) (checked), so the order is not an input: every order has the same limit.

*The flips are the rotation.* Each coin is a half-turn structure, and two half-turns about different axes compose to a rotation (III.8a). So between coins the column is turned, not reweighted. For the massless column with momentum \\(\\mathbf p\\) the axis \\(\\hat{\\mathbf s}=\\chi^\\dagger\\boldsymbol\\sigma\\chi\\) obeys, after one step (symbolic, and numerically to \\(O(\\epsilon^2)\\)),

\\\[\\hat{\\mathbf s}(t+\\epsilon)=\\hat{\\mathbf s}+\\epsilon\\,\\frac{2c}{\\hbar}\\,\\mathbf p\\times\\hat{\\mathbf s}+O(\\epsilon^2):\\\]

the axis precesses about \\(\\hat{\\mathbf p}\\) by the angle \\(2|\\mathbf p|c\\,\\epsilon/\\hbar\\) per step. “The direction flips at every step” and “the axis rotates” are the same sentence; the flips are what a rotation looks like when each component may only be \\(\\pm c\\). That is IV.2a’s circulation, derived. A path is a sequence of axis choices and its weight is the product of the projectors it passed through — a product of overlaps \\(\\langle\\chi(\\hat a_k)|\\chi(\\hat a\_{k+1})\\rangle\\), positive magnitudes times the solid-angle phase; the phase appears as soon as two different axes are composed, \\(P_x^+P_y^+=\\frac14(1+\\sigma_x+\\sigma_y+i\\sigma_z)\\), and that is where IV.2a’s \\(-2\\) lives.

Checked numerically as well (`scripts/dirac3d_momentum.py`): the step is translation-invariant, so per momentum it is a \\(4\\times4\\) matrix \\(U(\\mathbf p)\\); \\(U(\\mathbf p)^N\\) and \\(e^{-i\\hat H(\\mathbf p)t}\\) were computed on the same grid, applied to the same smooth start, and compared point by point inside the cone at \\(t=3.2/m\\). With \\(N=32\\) and \\(N=64\\) steps the automaton converges onto the Dirac column (the gaps shrink by \\(1.6\\) to \\(6.5\\) times entry by entry, \\(4.6\\) times in sum), and the Dirac column agrees with the bulk formulas above — \\(i\\,m^2J_1(ms)/4\\pi s\\) in the lower block, \\(m^2J_2(ms)(t+\\boldsymbol\\sigma\\cdot\\mathbf x)/4\\pi s^2\\) in the upper — up to what smoothing a start of width \\(0.35/m\\) produces — \\(+10\\%\\) in the lower block, \\(2\\)–\\(5\\%\\) in the upper — which is \\(\\frac{w^2}2\\nabla^2\\) of the formula, to \\(2.7\\%\\). Near the cone the \\(\\delta'\\) skin dominates, as it should: at the point \\((0,0,2.5)\\), \\(0.7/m\\) inside the cone, the smoothed continuum column has the opposite sign from the unsmoothed bulk formula (\\(-0.030\\) against \\(+0.040\\)) because the smoothing reaches the skin, and the automaton agrees with the smoothed column there (\\(-0.032\\)), not with the formula. A first attempt in position space with a point-like start was swamped by lattice-scale modes; that is an artefact of the lattice’s dispersion at short wavelength, not of the physics, and the momentum-space evaluation removes it.

*Two questions, two sides: the table in three dimensions.* III.14–15 had one function on both sides of the light cone because both sides answered the *evolution* question. In three dimensions that row has one cell proven empty (IV.2a), so it is worth laying out all four cells and saying which process fills which:

|  | outside the cone, \\(s'^2=r^2-c^2t^2\\) | inside the cone, \\(s^2=c^2t^2-r^2\\) |
| --- | --- | --- |
| evolution question | \\(I_1(ms')/s'\\): no one-clock process (IV.2a) | \\(J_1(ms)/s\\), \\(J_2\\): the three-coin automaton |
| source question | \\(\\frac{m}{4\\pi^2s'}K_1(ms')\\): Brownian motion killed at \\(m^2\\) (III.16) | \\(\\frac{im}{8\\pi s}H_1^{(2)}(ms)\\): the automaton with Feynman’s rule |

The bottom row is *one* function, \\(\\Delta_F(t,\\mathbf x)=\\int\\frac{d^3p}{(2\\pi)^3}e^{i\\mathbf p\\cdot\\mathbf x}\\,\\frac{e^{-iE|t|}}{2E}\\): at \\(t=0\\) it is the Brownian \\(K_1\\) (ratio checked, \\(1.001\\)), inside the cone it is \\(\\frac{im}{8\\pi s}H_1^{(2)}(ms)\\) (the ratio \\(\\Delta_F\\,s/H_1^{(2)}\\) is the same constant at three unrelated points), and the two are each other’s continuation — \\(K_1\\) across the cone is the Hankel function, both Bessel functions \\(J_1\\) and \\(Y_1\\) with the \\(i\\) between them, not \\(J_1\\) alone. What turns the evolution kernel into this is Feynman’s rule, which in this document’s words is a rule about the *time sense*: carry the positive-energy half of the column forward in \\(t\\) and the negative-energy half backward, \\(K_F=e^{-i\\hat Ht}\\Pi\_+\\) for \\(t>0\\) — the \\(\\Pi\_\\pm\\) of III.15 item 2, each given its sense. Applied to the automaton with the automaton’s *own* projector (built from the eigenphases of its unitary step), the lower block converges onto \\(m\\Delta_F\\): at \\(N=32,64\\) steps the gap to the continuum column closes, and the continuum column with the smoothing narrowed from \\(0.35/m\\) to \\(0.15/m\\) lands on \\(\\frac{im^2}{8\\pi s}H_1^{(2)}(ms)\\) to \\(2\\%\\) at six interior points (`scripts/feynman_bottom_row.py`). Compared with the evolution kernel’s lower block \\(\\frac{im^2}{4\\pi s}J_1\\): Feynman’s rule keeps half of the \\(J_1\\) and adds the \\(Y_1\\) — the part that remembers there is a source. So in three dimensions the two sides of the cone are joined on the source row by the same two processes, the Euclidean walk outside and the three coins inside; on the evolution row only the inside exists — and the constructive form of the empty cell is IV.4c’s last paragraph: boosted positive \\(1{+}1\\) walks, summed with the signs relativity puts on the readings.

*What went in (IV.10).* \\(\\Delta\\): Part II with \\(\\tau\\) summed (III.16), checked as a killed Brownian occupation density. \\(S\\): IV.9’s factorization. \\(K\\): IV.5’s Hamiltonian and the on-shell sum. The three-coin process: III.1 applied axis by axis, in turn (IV.1, IV.2a) — no new input; its limit is IV.5, taken symbolically, and its kernel matches \\(K\\) numerically. No new input. Feynman’s rule as the time-sense rule for the two halves of the column. It is not an input of its own. Given the vacuum, it is exclusion: a walker added to the vacuum propagates by \\(U^t\\Pi\_+\\), and a walker removed from it by the conjugate of \\(U^t\\Pi\_-\\), which is \\(U^t\\Pi\_+\\) again — a leading loop, of the opposite charge when a field is on (VI.4 (C), VI.6). What is taken is one thing: that the vacuum is the configuration the long-clock count selects (VI.4 (C), (H); VI.5). The carrier’s state is a separate choice (VI.2 (F′)). What the three-coin automaton is: a process in the sense of the \\(1{+}1\\) checkerboard — each sub-step a classical coin on the column’s axis — and not a process with positive scalar weights, which IV.2a says cannot exist. As a unitary lattice step it is the Dirac cellular automaton of Bialynicki-Birula, *Phys. Rev. D* 49, 6920 (1994).

---

# Part V — The electromagnetic field: where it enters the walk

## V.0 The question, and the knobs

*The picture first.* Parts III and IV built a particle from one walk: a time sense that flips at a fixed rate, carried at speed \\(c\\). The walk has very few parts. In one clock step \\(dx\\) of III.2,

\\\[P\_\\varepsilon(t;x+dx)=\\underbrace{\\Big(1-\\tfrac ac\\,dx\\Big)}\_{\\text{survive}}\\,P\_\\varepsilon\\Big(\\underbrace{t-\\varepsilon\\tfrac{dx}c}\_{\\text{advance}};x\\Big)+\\underbrace{\\tfrac ac\\,dx}\_{\\text{flip}}\\,P\_{-\\varepsilon}(t;x),\\\]

there is a coin, a speed and a survival weight, and each can be changed either by the same amount for both senses or by opposite amounts. That is the whole list of places a force could enter. This Part turns each knob and reads off what comes out, finds that the electromagnetic field fits none of them as a real weight, and builds the one step that does carry it. Everything below is \\(1{+}1\\) unless said otherwise; the \\(3{+}1\\) statements use the matrices of IV.4b.

*Symbols introduced in this Part.* \\(\\Phi(t,x)\\) and \\(A_x(t,x)\\) are the electric and vector potentials and \\(e\\) the charge; \\(F\_{\\mu\\nu}\\) is the field, with \\(\\mathbf E=-\\nabla\\Phi-\\partial_t\\mathbf A\\) and \\(\\mathbf B=\\nabla\\times\\mathbf A\\). \\(D\_\\mu=\\partial\_\\mu-\\frac{ie}{\\hbar}A\_\\mu\\) is a move in the presence of the potential. A prime on a counter, \\(Q'\_\\varepsilon\\), marks the added detector of III.8c (B): \\(Q\_\\varepsilon\\) is \\(u\_\\varepsilon\\) and \\(Q'\_\\varepsilon\\) is \\(v\_\\varepsilon\\). (\\(A_x\\) is a number, the component of the potential; Part IV’s direction matrices \\(A_k\\) are a different object and do not appear in this Part’s formulas.)

## V.1 Turning each knob

Read in \\(t\\) with the loss taken out, III.5a (H) is \\(\\partial_tQ=-cA\\,\\partial_xQ+aTQ\\). Each knob adds a term to it. (`scripts/process_knobs.py`.)

| knob | the step changed by | what is added | keeps \\(Q\_+^2+Q\_-^2\\) | what it is |
| --- | --- | --- | --- | --- |
| coin, same for both senses | \\(a\\to a+v\\) | \\(v\\,TQ\\) (and a common loss, taken out by \\(e^{-vx/c}\\)) | yes | **a change of mass:** \\(mc^2\\to mc^2+\\hbar v\\) |
| coin, opposite for the two senses | out of \\(+\\) at \\(a+v\\), out of \\(-\\) at \\(a-v\\) | \\(-v\\,(Q\_++Q\_-)\\) in both lines | no: rate \\(-2v(Q\_++Q\_-)^2\\) | a real rescaling of the counts |
| speed, by sense | \\(t\\to t+\\varepsilon(1+k\_\\varepsilon)\\,dx/c\\) | each line divided by \\(1+k\_\\varepsilon\\): per unit \\(t\\), speed \\(c/(1+k\_\\varepsilon)\\) and feed \\(a/(1+k\_\\varepsilon)\\) | the weighted \\((1+k\_+)Q\_+^2+(1+k\_-)Q\_-^2\\); the plain length only if \\(k\_+=k\_-\\) | **a change of the light cone:** a medium, or gravity |
| survival, by sense | extra loss \\(l\_\\varepsilon\\,dx\\) | \\(-c\\,\\varepsilon l\_\\varepsilon\\,Q\_\\varepsilon\\) on the diagonal | no | the slots where \\(\\Phi\\) (\\(l\_\\varepsilon=\\varepsilon l\\), by sense) and \\(A_x\\) (common \\(l\\)) sit in Dirac in the clock reading — the \\(t\\)-reading swaps the two, V.2 (E) — filled with a **real** number |

(The table is for constant knobs. A knob that varies with \\(x\\) is not just the same term with a coefficient: a position-dependent speed, \\(-c\\,\\mathrm{diag}(\\ldots)\\partial_x\\) with the diagonal a function of \\(x\\), conserves \\(\\sum Q^2\\) only with the symmetrized form \\(-\\tfrac c2\\{\\mathrm{diag}(\\ldots),\\partial_x\\}\\), a derivative correction; and a loss factor that varies in \\(x\\) adds terms of its own. The “yes” entries are statements about constant coefficients.) The last row is the one that matters. Put the potentials into Dirac and solve for the clock derivative:

\\\[c\\,\\partial_xQ\_\\varepsilon=-\\varepsilon\\,\\partial_tQ\_\\varepsilon+aQ\_{-\\varepsilon}\\ -\\ i\\varepsilon\\frac{e\\Phi}{\\hbar}\\,Q\_\\varepsilon\\ +\\ ic\\frac{eA_x}{\\hbar}\\,Q\_\\varepsilon .\\\]

The field sits on the diagonal, as a survival weight, exactly where the fourth knob acts: \\(\\Phi\\) with opposite signs for the two senses, \\(A_x\\) with the same sign. But the number in the slot is imaginary. A real survival skew rescales the counts from point to point; that was Weyl’s gauge theory of 1918, and it failed for this reason until London put in the \\(i\\) in 1927. An imaginary survival weight is a phase: over one step of the path, displacing by \\(dt=\\varepsilon\\,dx/c\\) and \\(dx\\), it is

\\\[e^{\\frac{ie}{\\hbar}(-\\Phi\\,dt+A_x\\,dx)} .\\\]

So the field enters the walk as a phase on each step, and nowhere else: not in the coin, not in the speed, not in the spread of the steps. The rest of this Part is what that one statement produces.

*Why no real knob can do it.* Among all first-order real changes to the two-counter walk, the only ones that keep the length \\(Q\_+^2+Q\_-^2\\) are the coin changes \\(\\mu T\\), and \\(T\\) anticommutes with \\(A\\): it changes the sense. There is no real, length-keeping change that leaves the sense alone (checked). Now take the four detectors of III.8c (B), which read the time sense and the position of the clock hand separately, and ask again for real, length-keeping changes that leave both the sense and the flips alone: there is exactly one, a turn of the hand — the walker’s own turn (checked; with \\(n\\) sets of detectors there would be the \\(\\mathfrak{so}(n)\\) of turns among them, which is the door to non-abelian fields and is not opened here). That is the knob electromagnetism uses.

## V.2 The step that carries the field: the walker’s own turn

*The picture first.* The walker is a clock (III.8c (D)), and the four detectors of III.8c (B) read where its hand stands separately from its time sense. A flip moves the hand a quarter on and reverses the time sense. The field does the first without the second: it turns the hand, bit by bit, as the walk moves, and leaves the time sense alone. On the detectors that turn carries each original detector onto the added one of the same time sense, \\(Q\_\\varepsilon\\) onto \\(Q'\_\\varepsilon\\). A walker read by the original two alone cannot be turned this way, because the turn takes it out of them: it is neutral. A charged particle needs all four. This is the same statement as “a real Dirac field is neutral (Majorana) and a charged one is complex”: the original two and the added two are the real and imaginary parts.

**(A) The step.** Counters \\(P\_\\varepsilon\\) (the original two) and \\(P'\_\\varepsilon\\) (the added two), before the loss is taken out. The turn in item 3 can make a counter negative: the step is real and keeps \\(P\_\\varepsilon^2+P\_\\varepsilon'^2\\), but it carries signed weights, as the \\(t\\)-reading of III.5a already does through \\(T\\); it is not a probability process. In each clock step \\(dx\\):

1. *advance:* \\(t\\to t+\\varepsilon\\,dx/c\\);
2. *flip:* with probability \\(a\\,dx/c\\) the sense reverses and the hand moves a quarter on: \\(P\_\\varepsilon\\) feeds \\(P\_{-\\varepsilon}\\) and \\(P'\_\\varepsilon\\) feeds \\(P'\_{-\\varepsilon}\\), as in III.8c (B);
3. *turn:* the hand is turned and the sense kept: the pair \\((P\_\\varepsilon,P'\_\\varepsilon)\\) of each sense is turned by the angle

\\\[\\theta\_\\varepsilon=-\\frac e\\hbar\\Big(\\varepsilon\\,\\Phi\\,\\frac{dx}c-A_x\\,dx\\Big),\\\]

which is the phase of V.1 read as an angle: \\(\\Phi\\) turns the two senses in opposite directions, \\(A_x\\) turns them the same way. The turn keeps \\(P\_\\varepsilon^2+P\_\\varepsilon'^2\\).

**(B) One step.** For the original two (the added two with \\(\\sin\\to-\\sin\\) on the first bracket and \\(\\cos\\) on the second),

\\\[P\_\\varepsilon(t;x+dx)=\\cos\\theta\_\\varepsilon\\Big\[(1-\\tfrac{a\\,dx}c)P\_\\varepsilon(t-\\varepsilon\\tfrac{dx}c)+\\tfrac{a\\,dx}cP\_{-\\varepsilon}(t)\\Big\]-\\sin\\theta\_\\varepsilon\\Big\[(1-\\tfrac{a\\,dx}c)P'\_\\varepsilon(t-\\varepsilon\\tfrac{dx}c)+\\tfrac{a\\,dx}cP'\_{-\\varepsilon}(t)\\Big\].\\\]

**(C) Taylor to first order in \\(dx\\)** (\\(\\cos\\theta=1\\), \\(\\sin\\theta=\\theta\\), \\(P(t-\\varepsilon\\,dx/c)=P-\\varepsilon\\frac{dx}c\\partial_tP\\)), subtract \\(P\_\\varepsilon\\), divide by \\(dx\\), multiply by \\(c\\):

\\\[c\\,\\partial_xP\_\\varepsilon=-\\varepsilon\\,\\partial_tP\_\\varepsilon-aP\_\\varepsilon+aP\_{-\\varepsilon}+\\frac e\\hbar\\big(\\varepsilon\\Phi-cA_x\\big)P'\_\\varepsilon,\\qquad c\\,\\partial_xP'\_\\varepsilon=-\\varepsilon\\,\\partial_tP'\_\\varepsilon-aP'\_\\varepsilon+aP'\_{-\\varepsilon}-\\frac e\\hbar\\big(\\varepsilon\\Phi-cA_x\\big)P\_\\varepsilon .\\\]

III.2 twice over, once for the original two and once for the added two, plus the turn. The turn is the only place the field appears.

**(D) Loss out,** \\(P=e^{-ax/c}Q\\), as in III.5a (C): the \\(-aQ\\) cancels the loss and the lines keep their form without it.

**(E) Read in \\(t\\)** (multiply each line by \\(\\varepsilon\\)):

\\\[\\partial_tQ\_\\varepsilon=-\\varepsilon c\\,\\partial_xQ\_\\varepsilon+\\varepsilon a\\,Q\_{-\\varepsilon}+\\frac e\\hbar\\big(\\Phi-\\varepsilon cA_x\\big)Q'\_\\varepsilon,\\qquad \\partial_tQ'\_\\varepsilon=-\\varepsilon c\\,\\partial_xQ'\_\\varepsilon+\\varepsilon a\\,Q'\_{-\\varepsilon}-\\frac e\\hbar\\big(\\Phi-\\varepsilon cA_x\\big)Q\_\\varepsilon .\\\]

In the \\(t\\)-reading \\(\\Phi\\) acts the same way on both senses and \\(A_x\\) oppositely: the two clocks swap the roles, as in III.5a (J).

**(F) Write each pair as one complex number,** \\(\\psi\_\\varepsilon=Q\_\\varepsilon+iQ'\_\\varepsilon\\). Line plus \\(i\\) times primed line, with \\(Q'\_\\varepsilon-iQ\_\\varepsilon=-i\\psi\_\\varepsilon\\):

\\\[\\boxed{\\partial_t\\psi\_\\varepsilon=-\\varepsilon c\\Big(\\partial_x-\\frac{ie}\\hbar A_x\\Big)\\psi\_\\varepsilon+\\varepsilon a\\,\\psi\_{-\\varepsilon}-\\frac{ie\\Phi}\\hbar\\,\\psi\_\\varepsilon}\\qquad\\text{i.e.}\\qquad i\\hbar\\,\\partial_t\\psi=\\Big\[c\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}\\big(\\hat p-eA_x\\big)+i\\hbar a\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}+e\\Phi\\Big\]\\psi .\\\]

III.9’s Dirac equation with \\(\\hat p\\to\\hat p-eA_x\\) and \\(+e\\Phi\\): minimal coupling, from the step. The walker’s own turn is the \\(i\\) that multiplies the field. In III.9’s words, the added detectors are “the freedom for the two components to be complex”: the quantum step of III.9 and the completed detectors are the same thing. The complex column is the four detectors written as two complex numbers; III.8c (B) says what the added two are, and III.15a (C) how an event at a detector is weighed.

**(G) Gauge freedom is where the clock hand is set.** With \\(\\Phi\\to\\Phi-\\partial_t\\chi\\), \\(A_x\\to A_x+\\partial_x\\chi\\) and the walker turned in its own plane by \\(e\\chi/\\hbar\\) at each point, the boxed equation is unchanged. Where the hand of the walker’s clock is said to stand is a choice made point by point, and the potentials are the bookkeeping of that choice.

**(H) Why the antiparticle has the opposite charge.** The angle \\(\\theta\_\\varepsilon\\) carries \\(\\varepsilon\\Phi\\): a backward-in-time segment is turned the other way by the same \\(\\Phi\\). The statement on the detectors: conjugating the complex form of (F) is counting the added two with the other sign, \\(Q'\\to-Q'\\) (the pair read with \\(J\\to-J\\)). It leaves the free real generator as it is and sends \\(e\\to-e\\): a solution with its added two counted with the other sign is a solution with the opposite charge. It is also the walker whose clock turns the other way (III.15e (C), V.4 (2)). The scalar coupling of V.1’s first row, \\(a\\to a+V/\\hbar\\), has no \\(\\varepsilon\\): it is a phase per unit of the walker’s *proper* time, \\(V\\,d\\tau/\\hbar\\), and acts alike on both senses. Electromagnetism is a phase per unit of the *coordinate* time the walker moves through, \\(e\\Phi\\,dt/\\hbar\\) with \\(dt=\\varepsilon\\,d\\tau\\). The two agree at rest, \\(E=mc^2+V\\), and part company for a moving particle and on the antiparticle branch: \\(E\_-=-(mc^2+V)\\) for the scalar, \\(E\_-=-mc^2+V\\) for the field.

(`scripts/two_copy_walk_em.py`: (B)–(G) symbolically, and a lattice run of the four real detectors against the complex column with a potential step: agreement to \\(10^{-15}\\), total length conserved, partial reflection. `scripts/imaginary_skew.py`: the slot of the field in the clock reading, and the form on the four detectors. `scripts/charge_needs_doubling.py`: no real \\(J\\) on two counters.)

## V.3 Where the field is in the noise: an output, not an input

*The picture first.* Part II put the metric into the walk as the spread of the steps, \\(\\langle\\xi^\\mu\\xi^\\nu\\rangle=\\eta^{\\mu\\nu}\\), an input. In the Dirac process the spread is not an input: it is computed from the step matrices, as the two-step moment, and that is where the field shows up.

*Why the spread cannot carry the field as an input.* Any noise matrix enters the equation as \\(N^{\\mu\\nu}\\partial\_\\mu\\partial\_\\nu\\), and \\(\\partial\_\\mu\\partial\_\\nu\\) is symmetric. So an antisymmetric addition to the noise drops out identically, \\((\\eta^{\\mu\\nu}+\\kappa F^{\\mu\\nu})\\partial\_\\mu\\partial\_\\nu=\\Box\\), and a symmetric one of the same signature is \\(\\Lambda\\eta\\Lambda^{\\mathsf T}\\) for some linear \\(\\Lambda\\) — a change of coordinates, not a moving observer: a Lorentz \\(\\Lambda\\) leaves \\(\\eta\\) itself unchanged, \\(\\Lambda\\eta\\Lambda^{\\mathsf T}=\\eta\\). A rotation of the noise vector by the field, \\(\\xi\\to(1+\\frac emF\\Delta\\tau)\\xi\\), gives back \\(\\langle\\xi'\\xi'\\rangle=\\Lambda\\eta\\Lambda^{\\mathsf T}=\\eta\\): nothing. And no change of a second moment can produce the terms the field needs, which have fewer than two derivatives. The field goes into the first moment (the phase of V.1, which is the drift shift \\(m v^\\mu=\\partial^\\mu S-eA^\\mu\\) of Part II with the field), and the spread stays as it was. (`scripts/em_noise_factorization.py` A, `scripts/lorentz_on_noise.py`, `scripts/skewed_covariance.py`.)

*The two-step moment of the Dirac process.* With \\(\\gamma^0=\\beta=iT\\) and \\(\\gamma^k=\\beta A_k\\) (IV.4b’s real matrices), the step is \\(\\gamma^\\mu D\_\\mu\\) and two steps give

\\\[\\gamma^\\mu\\gamma^\\nu D\_\\mu D\_\\nu=\\underbrace{-\\eta^{\\mu\\nu}D\_\\mu D\_\\nu}\_{\\text{symmetric: the metric}}+\\underbrace{\\tfrac12\[\\gamma^\\mu,\\gamma^\\nu\]}\_{\\text{antisymmetric}}\\cdot\\tfrac12\[D\_\\mu,D\_\\nu\],\\qquad \[D\_\\mu,D\_\\nu\]=-\\frac{ie}\\hbar F\_{\\mu\\nu}.\\\]

The antisymmetric part is made of the six turns: \\(\\tfrac12\[\\gamma^0,\\gamma^k\]=A_k\\), the time-space turns (IV.4c’s boosts), and \\(\\tfrac12\[\\gamma^j,\\gamma^k\]=-A_jA_k\\), the space-space turns (IV.4b’s U-turn planes). Without a field, two moves commute and the antisymmetric part is invisible: this is the noise matrix of the free process, \\(\\mathrm{diag}(1,-1,-1,1)\\) in the clock \\(z\\) and \\(\\mathrm{diag}(1,1,1,-1)\\) in the clock \\(t\\) (`scripts/noise_matrix_dirac.py`). With a field, the two orders of two moves differ by \\(F\\), and the turns show. Computed out, for general \\(\\Phi(t,\\mathbf x)\\) and \\(\\mathbf A(t,\\mathbf x)\\):

\\\[(i\\hbar\\partial_t-e\\Phi)^2\\psi=\\Big\[c^2(\\hat{\\mathbf p}-e\\mathbf A)^2+m^2c^4\\ -\\ e\\hbar c^2\\,\\boldsymbol\\Sigma\\cdot\\mathbf B\\ +\\ ie\\hbar c\\,\\boldsymbol\\alpha\\cdot\\mathbf E\\Big\]\\psi,\\qquad\\Sigma_z=-iA_xA_y\\ \\text{(cyclic)}.\\\]

Klein–Gordon with the field, plus the two antisymmetric pieces: the space-space turn times \\(\\mathbf B\\) and the time-space turn times \\(\\mathbf E\\). In the slow limit \\(\\sqrt{H^2}=mc^2+\\frac{\\pi^2}{2m}-\\frac{e\\hbar}{2m}\\boldsymbol\\Sigma\\cdot\\mathbf B+\\ldots\\): the magnetic moment is \\(\\frac em\\cdot\\frac\\hbar2\\boldsymbol\\Sigma\\), \\(g=2\\), without being put in. In the walk’s words: the magnetic moment is the \\((x,y)\\) turn left behind by two U-turns taken in the two orders (IV.4c (G)), and the field couples to it because, in a field, the two orders also differ by \\(B_z\\). (`scripts/em_noise_factorization.py` B–D.)

*The six turns are one object.* Under a boost along \\(x\\) (`scripts/turns_under_boost.py`): \\(A_z\\to\\cosh\\eta\\,A_z-\\sinh\\eta\\,A_xA_z\\) and \\(A_xA_z\\to\\cosh\\eta\\,A_xA_z-\\sinh\\eta\\,A_z\\), while \\(A_x\\) and \\(A_yA_z\\) are unchanged. A time-space turn picks up a space-space part and the other way round, exactly as \\((\\mathbf E,\\mathbf B)\\) mix, and the set of six is closed. The field \\(F\_{\\mu\\nu}\\) is a six-component object of the same kind, which is why it couples to the turns and to nothing else: the exchange \\(A_zT\\) (the mass) is unchanged by every boost, and none of the six turns has a component along it.

*What a tilt of the step matrices gives instead.* Tilting the matrices from point to point, \\(\\gamma'^\\mu(x)=S^{-1}(x)\\gamma^\\mu S(x)\\), does produce a first-order coupling, \\(-\\gamma^\\mu(\\partial\_\\mu S)S^{-1}\\), but it is valued in the generators of \\(S\\): a local boost gives a real (Weyl) vector coupling, a local rotation gives \\(\\tfrac12(\\partial\_\\mu\\theta)A_xA_y\\), the spin connection of a rotating frame. That is the coupling to a local frame, gravity, not electromagnetism. The electromagnetic generator must commute with every step matrix, and among real matrices on the original four counters only the identity does. (`scripts/tilted_matrix_noise.py`.)

## V.4 The effective flip rate

*The picture first.* The coin turns the column at the rate \\(a=mc^2/\\hbar\\) per unit of proper time in the frame where the column’s flow vanishes (\\(\\cos\\theta=0\\) in III.11; there \\(\\bar S=-\\hbar a\\tau\\) measures that time), in every frame and every gauge: that is the mass, and the field does not change it. The statement is about the column: a segment of the path is light-like and has no proper time, and the path’s clock \\(x\\) is frame-bound (III.16). What an observer sees turning in time is more than the coin. On the detectors the mass and the potential are the same operation, a turn of the walker in its own plane. The mass makes that turn together with the exchange of the forward and the backward detectors; the potential makes it without the exchange. For a walker at rest with the same amount forward and backward the exchange does nothing, the two turns are one, and its clock turns at \\(a+e\\Phi/\\hbar\\): for that walker a potential cannot be told apart from a change of the flip rate. With opposite amounts the same potential slows the clock: the opposite charge. In motion the two part, and the field is seen to be what it is: a change of the effective flip rate, not of the flip rate.

*Symbols.* \\(r=(Q\_+,Q\_-,Q'\_+,Q'\_-)^{\\mathsf T}\\), the four detectors in III.15e’s order, with III.15e’s matrices \\(\\mathbb T\\) (the time flip) and \\(\\mathbb J\\) (the walker’s own turn), written out in (1). \\(M\\), \\(\\tilde M\\), \\(R\_{yz}\\), \\(d_x,\\dots,d_t\\) and \\(\\omega\\) are IV.5b’s.

**The flip rate** is the coin: \\(a=mc^2/\\hbar\\), the rate of time flips, i.e. of pair creation and annihilation (III.3), per unit of the column’s proper time as above. It is a Lorentz scalar as the mass is, and gauge-invariant, and the field leaves it alone: in the square above, \\(m^2c^4\\) stands untouched, with no term of the form (mass × field), and in V.2’s step the flip probability is \\(a\\,dx/c\\) with or without \\(\\Phi\\). The mass shell with the field, \\((E-e\\Phi)^2-c^2(\\mathbf p-e\\mathbf A)^2=m^2c^4\\), says the same thing: the kinetic momentum keeps the length \\(m\\).

**The effective flip rate** is the total turn rate of the column in time, \\(a\_{\\text{eff}}=E/\\hbar\\).

**(A) At rest: one turn, with and without the exchange.** V.2 (E) with nothing depending on \\(x\\) and a scalar potential only, the four lines as one matrix equation:

\\\[\\partial_t\\begin{pmatrix}Q\_+\\\\Q\_-\\\\Q'\_+\\\\Q'\_-\\end{pmatrix}=\\Bigg\[\\,a\\underbrace{\\begin{pmatrix}0&1&0&0\\\\-1&0&0&0\\\\0&0&0&1\\\\0&0&-1&0\\end{pmatrix}}\_{\\mathbb T}-\\frac{e\\Phi}\\hbar\\underbrace{\\begin{pmatrix}0&0&-1&0\\\\0&0&0&-1\\\\1&0&0&0\\\\0&1&0&0\\end{pmatrix}}\_{\\mathbb J}\\Bigg\]\\begin{pmatrix}Q\_+\\\\Q\_-\\\\Q'\_+\\\\Q'\_-\\end{pmatrix}.\\tag{1}\\\]

1. the generator is \\(a\\,\\mathbb T-\\frac{e\\Phi}\\hbar\\,\\mathbb J\\). \\(\\mathbb T\\) and \\(\\mathbb J\\) commute and both square to \\(-1\\). Their product is the exchange of the forward and the backward detector of each hand position (III.8c (B)’s table: \\(Q\_+\\leftrightarrow-Q'\_-\\), \\(Q'\_+\\leftrightarrow Q\_-\\)):

\\\[\\mathbb J\\,\\mathbb T=\\begin{pmatrix}0&0&0&-1\\\\0&0&1&0\\\\0&1&0&0\\\\-1&0&0&0\\end{pmatrix},\\qquad(\\mathbb J\\mathbb T)^2=1,\\qquad a\\,\\mathbb T=-a\\,\\mathbb J\\,(\\mathbb J\\mathbb T).\\\]

2. so the mass term is the turn \\(\\mathbb J\\) after the exchange, and the potential’s term is the same turn without it. Split \\(r\\) into the part the exchange keeps and the part it reverses, \\(r_s=\\tfrac12(1+\\mathbb J\\mathbb T)\\,r\\) (the same amount forward and backward) and \\(r_d=\\tfrac12(1-\\mathbb J\\mathbb T)\\,r\\) (opposite amounts). Then (1) is

\\\[\\partial_t\\,r_s=-\\Big(a+\\frac{e\\Phi}\\hbar\\Big)\\,\\mathbb J\\,r_s,\\qquad\\partial_t\\,r_d=+\\Big(a-\\frac{e\\Phi}\\hbar\\Big)\\,\\mathbb J\\,r_d .\\tag{2}\\\]

3. three cases.

- *The same amount forward and backward* (\\(r_d=0\\): the walker launched at rest, \\(Q'=TQ\\)). Its clock turns at \\\[a\_{\\text{eff}}=a+\\frac{e\\Phi}\\hbar=\\frac{mc^2+e\\Phi}\\hbar .\\\] One turn, in one plane, at one rate. For this walker a potential and a change of the flip rate are the same thing.
- *Opposite amounts* (\\(r_s=0\\)). The clock turns the other way, at \\(a-e\\Phi/\\hbar\\). The potential that speeds the first clock slows this one: the opposite charge (V.2 (H)).
- *In motion.* Split the potential’s term with \\(1=\\mathbb J\\mathbb T+(1-\\mathbb J\\mathbb T)\\):

\\\[\\frac{e\\Phi}\\hbar\\,\\mathbb J=\\underbrace{\\frac{e\\Phi}\\hbar\\,\\mathbb J\\,(\\mathbb J\\mathbb T)}\_{\\text{the form of the mass term}}+\\underbrace{\\frac{e\\Phi}\\hbar\\,\\mathbb J\\,(1-\\mathbb J\\mathbb T)}\_{\\text{acts on }r_d\\text{ only}}.\\tag{3}\\\]

4. the first term is a flip rate changed from \\(a\\) to \\(a+e\\Phi/\\hbar\\). The second is zero on \\(r_s\\) and twice the turn on \\(r_d\\). For one clock run to rapidity \\(\\eta\\), \\(r_d\\) is \\(\\tanh\\tfrac\\eta2\\) times \\(r_s\\) (IV.5b (8)). So to first order in the speed a potential is a change of the flip rate, and beyond that it is not.

**(B) The same in three dimensions, on the frame.** With a scalar potential, IV.5b (12) is (the \\(A_k\\) here are IV.4b’s direction matrices)

\\\[\\partial_tM=-c\\,\\big(A_x\\partial_xM+A_y\\partial_yM+A_z\\partial_zM\\big)\\;\\underbrace{-\\,a\\,\\tilde M\\,R\_{yz}}\_{\\text{mass}}\\;\\underbrace{-\\,\\frac{e\\Phi}\\hbar\\,M\\,R\_{yz}}\_{\\text{potential}} .\\tag{4}\\\]

5. both terms turn the frame in its own plane; the mass turns the mirror image \\(\\tilde M\\), the potential turns \\(M\\) itself. Sums only, \\(\\tilde M=M\\): \\(\\partial_tM=-(a+e\\Phi/\\hbar)\\,M\\,R\_{yz}\\). Differences only, \\(\\tilde M=-M\\): the other way at \\(a-e\\Phi/\\hbar\\). In general, with \\(M-\\tilde M=2\\,(d_xA_x+d_yA_y+d_zA_z+d_t\\,\\omega)\\),

\\\[\\frac{e\\Phi}\\hbar\\,M\\,R\_{yz}=\\frac{e\\Phi}\\hbar\\,\\tilde M\\,R\_{yz}+\\frac{2e\\Phi}\\hbar\\,\\big(d_xA_x+d_yA_y+d_zA_z+d_t\\,\\omega\\big)\\,R\_{yz}:\\tag{5}\\\]

6. a change of the flip rate, plus a term on the differences only. This is (3).

**(C) What the field is not.** A change of the flip rate itself, V.1’s first row, is \\(a\\to a+V/\\hbar\\) in the mass term alone. It has no second term in (3). It changes the mass, and it turns both clocks faster:

\\\[\\partial_t\\,r_s=-\\Big(a+\\frac V\\hbar\\Big)\\,\\mathbb J\\,r_s,\\qquad\\partial_t\\,r_d=+\\Big(a+\\frac V\\hbar\\Big)\\,\\mathbb J\\,r_d ,\\\]

7. and gives \\((mc^2+V)^2\\) in the square. The field leaves the coin at \\(a\\) and shifts the two clocks in opposite directions. The two agree on a walker at rest with the same amount forward and backward, and nowhere else. So: the field deforms the *effective* flip rate, the total rate at which the walker turns in its own plane, and leaves the flip rate alone. Per time sense the same statement is III.15f (3): the forward hand turns at \\(a\\,(-W)/\\rho\_++b\_+\\) and the backward one at \\(a\\,(-W)/\\rho\_-+b\_-\\), the flip rate times the overlap ratio, plus the field’s rate.

**(D) In motion: three readings of one clock.** For a moving particle the one clock of III.8c (D), \\(\\theta=a\\gamma(t-Vx/c^2)\\), is read three ways. At a fixed position it turns at \\(\\gamma a\\), up with the motion, like \\(\\gamma m\\); sampled along the motion, at \\((E-pV)/\\hbar=a/\\gamma\\), the dilated rate; per unit of the walker’s own time, at \\(a\\) — the coin, this section’s opening sentence. With the field, III.8c (10d)–(10e) add the same three readings of one term: \\(e\\Phi/\\hbar\\) at a fixed position, \\(e(\\Phi-VA_x)/\\hbar\\) along the motion, \\(e\\Phi'/\\hbar\\) per proper time with \\(\\Phi'=\\gamma(\\Phi-VA_x)\\). So

\\\[\\text{per proper time:}\\qquad a\_{\\rm eff}=a+\\frac{e\\Phi'}\\hbar ,\\\]

and the lab’s \\(a\_{\\rm eff}=\\gamma a+e\\Phi/\\hbar\\) below is the fixed-position reading of the same clock:

\\\[a\_{\\text{eff}}=\\underbrace{\\gamma a}\_{\\text{frame-dependent, like }\\gamma m}+\\underbrace{\\frac{e\\Phi}\\hbar}\_{\\text{gauge-dependent, like potential energy}} ,\\\]

and the coin \\(a\\) is the invariant part, like the rest mass. A constant \\(\\Phi\\) raises \\(a\_{\\text{eff}}\\) and changes nothing observable (the densities are identical; the column differs by an overall phase: checked).

*What went in.* Nothing: V.2 (E)’s lines and products of III.15e’s matrices. (v43 printed the generator at rest as \\(aT+\\frac{e\\Phi}\\hbar J\\). With V.2 (E) and the \\(J\\) of IV.5a and III.15e the sign of the second term is minus, as in (1); the rate \\(a+e\\Phi/\\hbar\\) was right.)

(`scripts/detectors_in_master_names.py`: 6 and 8 for (1)–(5) and (C). `scripts/eight_detectors.py`: 15, the same on the eight detectors. `scripts/field_and_coin.py`, `scripts/energy_vs_coin.py`, `scripts/relativistic_flip_rate.py`: (D), as before.)

## V.5 Detection: what a detector registers, and why it is the square

*The picture first.* A detector is made of the same stuff as the particle: it is a walker too, with a flip rate of its own, and it can be coupled to the particle only through the knob the field uses, the walker’s own turn (V.1, V.2). This section builds the smallest such detector, a bit with two states, each a pair \\((Q,Q')\\), turned by an angle \\(g\\) while the particle is on its site, and asks three questions: what it receives, what it registers, and what it does. The answers: it receives the amplitude, phase and all; it registers a return — going out through the detector and coming back — which is III.15a’s event weight with the detector as the absorber endpoint, so it registers \\(Q^2+Q'^2\\) and nothing linear; and what it does is move, pushed by the gradient of that weight, which is the Casimir mechanism applied to one particle. Everything is \\(1{+}1\\) with the four detectors of V.2, written as real pairs \\((Q,Q')\\); nothing is complex.

*Symbols introduced in this section.* \\(x_D\\) is the detector’s site; \\(g\\) is the angle by which the detector’s bit is turned, in its \\((0,1)\\) plane, while the particle is at \\(x_D\\). Nothing else is new.

**(A) The coupling, and what is handed over.** The joint state of particle and detector is one real vector: at each site, a pair \\((Q,Q')\\) for the particle with the bit in state \\(0\\), and a pair with the bit in state \\(1\\). Before the coupling the bit is in \\(0\\) everywhere. While the particle is at \\(x_D\\) the bit is turned by \\(g\\), the same turn for \\(Q\\) and for \\(Q'\\). After one turn, state \\(1\\) at \\(x_D\\) carries (Here \\(g\\) is the detector’s turn angle, dimensionless; it is not the \\(g=2\\) of V.3, and (C)’s contact coupling is \\(\\kappa_D\\).)

\\\[\\sin g\\cdot\\big(Q,Q'\\big)(x_D),\\\]

linear in the amplitude, phase included (`scripts/go_and_come_back_detector.py`, check 1). So the detector does receive the full column, with its phase. Nothing is squared at this step.

**(B) What the detector registers is a return.** The detector’s own change — its effective flip rate shifting, V.4 — is how far the joint state has been turned away from itself: the overlap of the state before with the state after. At \\(x_D\\) the overlap pairs the unturned part with itself and the handed-over branch (bit in \\(1\\)) with the bit’s original state-\\(1\\) amplitude, which is zero:

\\\[Q\\cdot(\\cos g\\,Q)+Q'\\cdot(\\cos g\\,Q')+0\\cdot(\\sin g\\,Q)+0\\cdot(\\sin g\\,Q')=\\cos g\\,\\big(Q^2+Q'^2\\big),\\\]

and over the whole state the overlap is \\(1-(1-\\cos g)\\,(Q^2+Q'^2)(x_D)\\), a real number between \\(0\\) and \\(1\\) for \\(g\\le\\pi/2\\) (check 2; the absorber of (D) has \\(g=\\pi/2\\); for \\(g>\\pi/2\\) it can be negative — it is an overlap, not a probability). Three different numbers come out of this one operation and should be kept apart: this overlap; the weight of the clicked branch, \\(\\sin^2g\\,(Q^2+Q'^2)(x_D)\\), which is a probability; and the clock shift of (C), which needs the contact coupling there, not the angle \\(g\\). Going out hands over \\((Q,Q')\\); coming back compares with \\((Q,Q')\\); a comparison of a column with itself is its length. This is III.8c (A)’s measurement read from the detector’s side: every outward route to the detector with every return route from it that agrees on the endpoint, and their sum is \\(\\alpha_b^2\\) — the weight of the clicked branch, \\(\\sin^2g\\,(Q^2+Q'^2)\\), is that event weight times the coupling. The overlap and the clock shift of (C) are the two other bilinears a detector carries; the event weight is the second of the three. Turn the particle in its own plane at \\(x_D\\) by any angle: the handed-over branch changes, the overlap does not (check 3). The phase is in the detector, in the branch it was handed, and it would show if that branch were later brought back together with the rest of the particle; it is not in the detector’s clock, because a return pairs every component with itself.

The algebra agrees with the mechanism: under the walker’s own turn there is no invariant linear in \\((Q\_+,Q'\_+,Q\_-,Q'\_-)\\), and the quadratic invariants are a four-dimensional space spanned by \\(Q\_+^2+Q\_+'^2\\), \\(Q\_-^2+Q\_-'^2\\), \\(Q\_+Q\_-+Q'\_+Q'\_-\\), \\(Q\_+Q'\_--Q'\_+Q\_-\\) (`scripts/gauge_detector.py`, A1, A2). Anything coupled through the turn can respond only to lengths and cross-lengths of the column.

**(C) The clock shift and the force: the Casimir mechanism for one particle.** This is a second model of the detector, not the bit of (A): write the coupling as a contact term \\(V=\\kappa_D\\,\\delta(x-x_D)\\) on the detector’s effective flip rate, \\(\\kappa_D\\) a coupling with the units of energy times length (the angle \\(g\\) of (A) is dimensionless; the two are related by the dwell, \\(g\\sim\\kappa_D\\times\\) dwell\\(/\\hbar\\), for a sharply localized heavy pointer). To first order the shift is \\(\\kappa_D\\,(Q^2+Q'^2)(x_D)\\), independent of where the particle’s clock hand stands; to second order every matrix element pairs the particle’s column at \\(x_D\\) with an intermediate state, so the shift is a bilinear in the column again, \\(\\kappa_D^2\\,\\psi^\\dagger(x_D)\\,G\\,\\psi(x_D)\\). For a detector whose own level spacing is large against the particle’s energies the sum over the particle’s states closes, \\(G\\) is a multiple of \\(1\\), and the shift is \\(\\kappa_D^2(Q^2+Q'^2)(x_D)\\) times a detector-only number (\\(0.6\\%\\) off at a spacing of \\(50\\) against a band of width \\(3\\)); otherwise it is a combination of the four invariants of (B) (`scripts/detector_as_force.py`, 1a, 1b, 1b′, by exact diagonalization). The shift depends on where the detector is, and its gradient is a force, exactly (Ehrenfest):

\\\[\\frac{d\\langle p_D\\rangle}{dt}=-\\frac{\\partial\\,\\Delta E}{\\partial x_D}=-\\kappa_D\\,\\partial_x\\big(Q^2+Q'^2\\big)\\Big|\_{x_D}.\\\]

The detector is pushed down the slope of \\(Q^2+Q'^2\\), as an atom is pushed by the gradient of a light intensity (check 2). On a packet passing a heavy pointer on the lattice: the pointer is pushed as the packet arrives and pulled back as it leaves, net momentum \\(\\sim10^{-7}\\), but it is displaced, by \\(3\\times10^{-4}\\) lattice units in the run, and the displacement is a functional of \\(Q^2+Q'^2\\) only (3b); the force law itself, from the joint evolution of the particle and a quantum pointer, holds tick by tick to \\(1.3\\%\\) (3d). Its clock ran fast while the particle was there: the extra turns it accumulated are \\((\\kappa_D/\\hbar)\\int(Q^2+Q'^2)(x_D,t)\\,dt=\\kappa_D\\times1.016\\), the time-integrated event weight at the detector (3c′). In III.15a’s words the push and the pull-back are the two directions of the event, equal and opposite, so nothing linear survives; the displacement is the product.

This is the mechanism of the Casimir document, applied to one particle instead of the vacuum. There, the walls feel a force because the sum over closed paths between them depends on where the walls are, and the force is the gradient of that sum; “one particle alone pushes the walls apart; the attraction comes from loops”: one particle in mode \\(n\\) has \\(dE_n/dL=-n\\pi\\hbar c/L^2\\), pushing out, while the loop sum gives \\(+\\pi\\hbar c/24L^2\\), pulling in (check 4). The detection force is the one-particle term; the Casimir force is the vacuum term; both are gradients of a closed-loop sum with respect to a position.

**(D) The perfect absorber conserves the total.** Set \\(g=\\pi/2\\) and take the flipped branch out at every tick, counting it: the detector absorbs everything that reaches it and reflects nothing. Then at every tick

\\\[\\text{absorbed}+\\text{escaped through the far edge}+\\text{still in the box}=1,\\\]

to \\(4\\times10^{-15}\\) (`scripts/perfect_absorber.py`, A): in the run \\(0.478+0.514+0.008\\). The particle’s own norm falls to \\(0.008\\) and nothing is lost: the missing weight is in the record. The record per tick is \\((Q\_-^2+Q\_-'^2)(\\text{edge})\\), the flux into the wall, tick for tick (B); a partial absorber reflecting \\(30\\%\\) gives absorbed \\(0.334\\) with the same total of \\(1\\) (C). This is the absorbing-boundary rule of the standard theory (Tumulka, “Detection time distribution for Dirac particles”, arXiv:1601.04571 (2016)), here as the \\(g=\\pi/2\\) limit of the detector of (A). “Perfect” means the record is the whole flux.

**(E) A screen of detectors: the weights, and the one that happens.** Put a detector at every site of a screen and send the packet in (`scripts/screen_of_detectors.py`).

- Weak coupling: all seventy detectors are shifted, each by \\(\\kappa_D\\) times its own dwell; the shifts form the smooth profile of the time-integrated \\(Q^2+Q'^2\\), and their sum is \\(\\kappa_D\\times47.7\\): one unit of presence spread continuously over the screen. No detector clicks; every one moves a little.
- Perfect absorbers: the joint state becomes a sum of branches “absorbed at \\(D\\)” with squared lengths adding to \\(1\\) (\\(0.9846\\) absorbed plus \\(0.0154\\) that never arrived). The squared lengths are the Born weights of the outcomes. The dynamics keeps all the branches; it never picks one.

So the amplitude picture, which is the walk averaged over its paths, gives the distribution of outcomes exactly and never a single outcome. The single outcome is in III.15a: an event has one absorber endpoint; the detectors of the screen are the candidate endpoints, the branch weights above are the weights of the candidate events, and one event happens. This is the same division of labour as III.15a’s: the readings decide the weights, the event decides which. In the walk’s own clock the same thing is visible path by path: a screen at a fixed \\(x\\) with a detector at every \\(t\\) is hit by each of \\(200{,}000\\) paths of III.1 exactly once, and the hits are distributed as \\(P\_++P\_-\\) (3).

**(F) The clock the detector reads.** In III.1’s own clock, \\(x\\), a detector at a place is passed by every path exactly once, and its record is \\(\\int(P\_++P\_-)\\,dt=1\\): linear, positive, with no square and no residue (`scripts/perfect_absorber.py`, E). The square appears when the detector lives in \\(t\\), where real detectors live, and where the walk read in \\(t\\) carries the signs of \\(T\\) and an event needs both directions. Why the physical detector reads in \\(t\\) rather than in the walk’s clock is III.15b: read in the walk’s clock, a counter would click in the past of the trigger.

*What went in (V.5).* The detector as a walker of the same kind, coupled through the walker’s own turn — V.2’s knob, the only one a field can turn — turned by \\(g\\) while the particle is on its site: this is the model, and it is the one input of the section. Derived from it: that the detector receives the amplitude linearly and registers a return, so that its clock shift, its force and its record are \\(Q^2+Q'^2\\) and never linear; the force as the gradient of the event weight (the Casimir mechanism for one particle); conservation of the total with a perfect absorber; the Born weights of the outcomes on a screen. Taken from III.15a: that one event happens at one endpoint, with those weights. Why detectors read in \\(t\\): III.15b.

(`scripts/go_and_come_back_detector.py` (A), (B); `scripts/gauge_detector.py` (B), the invariants; `scripts/detector_as_force.py` (C); `scripts/perfect_absorber.py` (D), (F); `scripts/screen_of_detectors.py` (E); the Casimir document “Go-and-come-back sums: Euler–Maclaurin, zeta and the Casimir energy” for the vacuum term of (C).)

### What went in (Part V)

| the four knobs of the step: coin, speed, survival, each common or by sense | derived: mass, light cone, Weyl rescaling; none is the field (V.1) |
| --- | --- |
| the field as a phase \\(e^{ie(-\\Phi\\,dt+A_x\\,dx)/\\hbar}\\) per step | **the one input of this Part**: the angle by which the walker is turned in its own plane, its time sense untouched; its form is minimal coupling read per step (V.2) |
| the walker’s own turn | on the four detectors of III.8c (B), the only real, length-keeping knob that leaves sense and flips alone (real, not positive: the counters carry signed weights); the four detectors are III.9’s complex components. The added detectors themselves are taken (input 11) |
| Dirac with minimal coupling; gauge freedom as where the clock hand is set; antiparticle charge from \\(\\varepsilon\\) in the angle | derived (V.2) |
| the spread of the steps unchanged; the field in the antisymmetric two-step moment; Pauli term, \\(g=2\\) | derived (V.3) |
| a tilt of the step matrices is a local frame (spin connection), not the field | derived (V.3) |
| flip rate (the coin, invariant) vs effective flip rate \\(E/\\hbar\\) | derived (V.4): the mass and the potential are one turn, with and without the exchange of forward and backward; \\(a\\pm e\\Phi/\\hbar\\) at rest for the two clocks |
| why the field produces that turn angle | relocated by VI.2: (B), (D) show a non-removable link angle is a field and acts as the Lorentz force exactly, on every column; why there is one is VI.2 (F)’s single input |
| a detector as a walker coupled through the walker’s own turn, turned by \\(g\\) on its site | **the model of V.5**: derived from it, the detector registers a return, \\(Q^2+Q'^2\\) and nothing linear; its clock shift, its force (the gradient of the event weight: the Casimir mechanism for one particle), a perfect absorber’s conservation of the total, and the Born weights on a screen (V.5) |
| one event at one endpoint, with those weights | taken from III.15a |
| why detectors read in \\(t\\) rather than in the walk’s clock \\(x\\) | III.15b: a record is data on a spacelike cut (one input, of II.1’s kind); in the walk’s clock a counter clicks before the trigger |
| the field’s own dynamics (Maxwell), and the detector as a charged absorber | Maxwell in VI.2 (F), from one input of II.1’s kind (the same input that answers the row above); the charged absorber not in this Part |

---

# Part VI — Two walkers

## VI.1 What two particles are: one future per point and direction

*The picture first.* Parts III–V say what one particle is: an event with an emitter and an absorber, an out-leg and a return leg, carried by a half-angle column (III.15a, III.8c, III.8a). They do not say what two particles are, and the test of the definition is whether two identical particles exclude each other. By the rules of III.15a alone they do not: two walkers are Bose-like. What makes them fermions is not a new rule but an old one read once more: the noise of the walk belongs to the event, not to the walker (II.1, III.15b). Then a point of spacetime has one future in each direction, two walkers in the same state at the same point are one process, and the exchange of two identical walkers is a move the evolution can never make. The sign under exchange then follows from the rule read for every state (D), not from conservation, which keeps both signs. This section runs the test, shows two natural readings of the rule that fail and why, and gives the one that holds in any dimension.

*Symbols.* \\(D\\) and \\(E\\) are the two pairings of two emitters with two absorbers, \\(D=U\_{X_1A_1}U\_{X_2A_2}\\) and \\(E=U\_{X_2A_1}U\_{X_1A_2}\\), with \\(U\_{XA}\\) the reading of III.15a (C)1. Nothing else is new.

**(A) Two independent walkers are Bose-like.** The joint step of two walkers is \\(U\\otimes U\\). It commutes with the exchange of the two: the step is blind to which walker is which, and no sign lives in it. For identical walkers both pairings are events, and III.15a’s rule — every out-path with every return path, no condition — gives

\\\[W=(D+E)^2 ,\\\]

which at a common absorber \\(X_1=X_2\\) is \\(4D^2>0\\) at every reachable site (10 of them on the lattice of `pauli_test.py`). Exclusion would be \\((D-E)^2=0\\). The closed loop of one walker weighs \\((U^{\\mathsf T}U)\_{AA}=+1\\). The \\(-1\\) that the Dirac vacuum’s closed loops need is not on this loop: it is (D)’s exchange sign on a closed cycle of the configuration count, \\((-1)^{n-1}\\) for a cycle of \\(n\\) — \\(+1\\) for one walker returning to itself, as here, \\(-1\\) for two exchanging places (VI.5 (A), VI.4 (D)).

**(B) One future per point and direction.** The walker’s phase space is the pair (point, direction): the position and the light-like momentum. The column is a density on it, and the step is a measure-preserving map on it — \\(U^{\\mathsf T}U=1\\), the determinant of the step is one, which is III.15a’s conservation. Now ask what evolves two walkers that sit in the same state at the same event. The coin. And the coin is the event’s, not the walker’s: the noise has the covariance \\(\\eta\\) because it is a property of spacetime (II.1), and the metric it defines is read from either clock (III.15b (A)). So two walkers in the same state at the same event draw the same coin there, and every coin after; their futures coincide (`pauli_test.py` 4: identical at every step in 300 of 300 runs). They are not two processes but one, with twice the weight. Two walkers at the same event in *different* directions draw the same coin too, but one outcome sends them to different states; they stay two (separate at once, 300 of 300). This is a reading of II.1, not a consequence of it: II.1 fixes the covariance of the noise, and “the coin is drawn at the event, once” is one more thing said about it, of the same standing as III.15b’s cut. It is the input of this section. Hence:

\\\[\\text{a point of spacetime carries at most one process in each direction.}\\\]

In \\(1{+}1\\) a direction is the sense; in \\(3{+}1\\) it is the column — the sense and the axis. Two walkers in different directions may share a point; two in the same direction may not. This is uniqueness of the flow on the walker’s phase space, given that reading.

The configuration space of two identical walkers is the set of pairs of *distinct* phase points, with the two orderings identified (in \\(2{+}1\\) the identification can be braided, VI.1 (E)).

**(C) Two readings of the rule that fail, and what they teach.** The rule says what a point cannot hold; it still has to be read into the step. Two obvious readings:

*As a constraint on paths — “no shared point”.* Keep only path pairs that never share a point. Exclusion at once, and on a planar lattice the sum over such pairs is the determinant \\(D-E\\) (Lindström–Gessel–Viennot): it equals \\(D-E\\) at 80 of 104 endpoint pairs. It fails at the other 24, because the lemma needs the swap of two paths’ continuations at their first meeting to keep the weight, and the walk’s weight at a point is the coin — \\(\\cos\\) for a keep, \\(\\pm\\sin\\) for a flip — which depends on the sense; every first meeting is head-on (928 of 928), and the swap changes the weight in 368 of them. The event weights then sum to \\(0.827\\): **the constraint leaks.** A rule on paths is not consistent with III.15a’s closure.

*As an elastic collision.* Let two walkers that would land on the same point bounce instead, both senses reversed, no coin charged. This is a permutation of the non-coincident configurations, so the joint step stays orthogonal — weight conserved, \\(1\\) to \\(10^{-12}\\) — and no two walkers ever share a point. But its distribution over absorber pairs is neither the fermion’s \\((D-E)^2\\) nor the boson’s \\((D+E)^2\\): total variation \\(0.26\\) from the first and \\(0.27\\) from the second (\\(0.21\\) and \\(0.22\\) in positions alone). It is a third statistics, hard-core walkers. The reason is (B) exactly: **the rule forbids two processes in one state, not two at one point.** The free fermion has weight \\(0.045\\) on two walkers at one point with *opposite* senses — the head-on moment — and zero on the same point with the same sense. The bounce forbids the head-on moment, which the rule allows.

**(D) The sign, from exclusion in every state.** Keep the free step — a walker’s step does not change because another is near; away from a coincidence the two draw independent coins, and (B)’s shared coin is used only to say what a coincidence would be. Write the two-walker column as \\(C\_{IJ}\\), \\(I=(x,\\text{sense})\\) (the four-counter index in \\(3{+}1\\)): the free step is \\(C\\to UCU^{\\mathsf T}\\), one \\(U\\) per walker, the transpose because both indices are walker indices; exchange is \\(C\\to C^{\\mathsf T}\\); the step keeps both parities and conserves the weight of each (both sum to exactly \\(1\\); `pauli_test.py` 2c, `pauli_from_one_process_per_point.py` 5, `pauli_exclusion_every_state.py` 3) — so neither the step nor conservation chooses the sign. Nor does exclusion read in the site–sense basis alone. It is true that the coin at a head-on point produces same-state components from the pair \\((x,+;\\,x,-)\\), which the even combination carries and the odd one cancels (`pauli_exclusion_every_state.py` 4); but a symmetric column with nothing on any same-site pair at all — \\(C=(\\mathbf 1\\mathbf 1^{\\mathsf T}-1)\\otimes 1\\), two walkers in the same sense, never on one site, spread uniformly — is kept by the free step exactly, tick after tick (ibid. 1). “Free step plus no two walkers on one site and sense” admits an even sector; an earlier version of this section claimed otherwise, and the claim was wrong. What settles the sign is (B) read as it is meant: the coin is drawn once per *state*, and a state of the walker is any column \\(v\\). That the rule cannot be stated for sites and senses only is one line of the free step. With \\(C\\to UCU^{\\mathsf T}\\),

\\\[\\big(UCU^{\\mathsf T}\\big)\_{II}=\\sum\_{J,K}U\_{IJ}\\,C\_{JK}\\,U\_{IK},\\\]

so the same-state entries after one tick are made of the entries \\(C\_{JK}\\), \\(J\\neq K\\), before it: the coin sends \\(e_I\\) to a combination of several \\(e_J\\), and “no weight on any same site and sense” is a statement about one basis that the next tick reads in another. On the lattice of III.15a a symmetric column with every \\(C\_{II}=0\\) has same-state entries of \\(0.26\\) after one tick, and read in the basis the tick makes, \\((Ue_I)^{\\mathsf T}C(Ue_I)\\), the same column already has \\(0.32\\); the antisymmetric column keeps \\(C\_{II}=0\\) in every basis for every tick, to \\(10^{-16}\\) (`exclusion_must_be_basis_free.py` 1–3). The only form of the rule that every basis agrees on is the one that holds for every column:

Exclusion in every state is \\(v^{\\mathsf T}Cv=0\\) for every column \\(v\\): two walkers are never in one state. Then

\\\[v=e_I:\\ C\_{II}=0;\\qquad v=e_I+e_J:\\ 0=C\_{II}+C\_{IJ}+C\_{JI}+C\_{JJ}=C\_{IJ}+C\_{JI}\\qquad\\Longrightarrow\\qquad C^{\\mathsf T}=-C ,\\\]

two lines of polarization (for complex columns the same, with the bilinear \\(v^{\\mathsf T}Cv\\); ibid. 2), and the free step keeps it, \\((UCU^{\\mathsf T})^{\\mathsf T}=UC^{\\mathsf T}U^{\\mathsf T}=-UCU^{\\mathsf T}\\). So

\\\[\\text{one process per state, every state}\\quad\\Longrightarrow\\quad C^{\\mathsf T}=-C,\\qquad W=(D-E)^2 :\\\]

exclusion exact in every basis, conservation exact, the exchange of identical walkers carries \\(-1\\), and no same-state pair ever carries weight. The reading of (B) made the two orderings one configuration; the step cannot choose the sign and does not; exclusion stated for every state chooses odd. Said plainly: antisymmetry is derived from exclusion alone, by polarization; exclusion, for every state, is the input. That it holds for every state and not only for sites and senses is forced by the step, as just shown: a site-and-sense rule is not a rule the tick can keep. Exclusion is the input; its basis-free form is derived.

*Remark: what determinism gives, and what it does not.* (B)’s step is a volume-preserving map of the walker’s phase space, and that alone is already a classical exclusion: an incompressible flow never drives two elements of phase fluid onto each other, so an occupation that is \\(0\\) or \\(1\\) stays \\(0\\) or \\(1\\) (Lynden-Bell’s argument for collisionless stars, 1967, which ends in a Fermi-type distribution; “deterministic and reversible” read as “volume-preserving” is Carcassi and Aidala’s). It *keeps* an exclusion; it does not make one — the same flow keeps an occupation of \\(2\\). Three statements are needed, one line each. (i) *The unit*: a walker is identified by its state, so two in one state are one walker — (B)’s rule, Lynden-Bell’s “one fixed density”. (ii) *It is kept*: the step preserves volume, \\(\\det U=1\\). (iii) *The sign*: the rule holds in every state, (D). Said with volumes the three are one object: \\(N\\) walkers are the volume element spanned by their columns,

\\\[u_1\\wedge u_2\\wedge\\dots\\wedge u_N,\\qquad u\\wedge u=0,\\qquad u\\wedge v=-\\,v\\wedge u,\\qquad(Uu_1)\\wedge\\dots\\wedge(Uu_N)=\\big\[\\text{the }N\\times N\\text{ minors of }U\\big\]\\ \\big(u_1\\wedge\\dots\\wedge u_N\\big).\\\]

Two equal columns span no volume: exclusion. The orientation of the volume: the exchange sign. The step acting on volumes: the many-walker step of VI.5 (1). The walk’s own clock shows where the criterion bites. III.2’s probability process, loss in, shrinks the area of the pair of counters as \\(e^{-2a\\tau}\\) — the sense is being forgotten (VI.7 (3)); with the loss out the same clock is a squeeze of determinant \\(1\\) (\\(\\operatorname{Tr}M_k=0\\)), and the reading in \\(t\\) is a rotation. The volume-preserving readings are the two the document works in. Not obtained this way: the size of the unit (here one slot per site and sense; Lynden-Bell’s density is whatever the start had), and spin–statistics as a theorem. (`scripts/exclusion_as_volume.py` 1–4.)

**(E) Three dimensions.** Nothing in (B) or (D) used a planar lattice, non-crossing paths or a half-turn: (B) uses that the noise is the event’s, (D) that \\(U\\otimes U\\) commutes with exchange and that odd columns vanish on same-state configurations, which holds for the four-counter column of IV.4b as for the two-counter one. In \\(2{+}1\\) the identification of the two orderings by a reflection is not the only one available (the exchange can be braided), so what follows is for \\(3{+}1\\) and \\(1{+}1\\). So in \\(3{+}1\\): two walkers at one point with opposite spin along an axis are two processes and may coexist; the same spin may not; a forward and a backward sense at one point — particle and antiparticle — may. That is Pauli as it is used.

**(F) Bosons.** The rule is about processes that *have* a direction at a point. Part II’s walker has none: its step is a kick with a drift and a spread and no velocity (I.5, II.4), no sense, no column. For it “in each direction” is empty, and (B) imposes nothing: unconstrained paths, symmetric counting. The Casimir document already counts Part II’s loops with the plus sign and needs the minus for the Dirac loops; here that is (B) reading the two kinds of walker. Statistics, in this document, follows from whether the process carries a direction at a point: the directed walker is the fermion, the undirected walker the boson; the link angle of Part V is a relation between points, not a process at one, and what carries it from event to event is Part II’s walker at \\(m=0\\) with a one-form index — a boson (VI.2 (E)). This is a reading, consistent with everything on the page, and not a theorem: why the half-angle column must be the directed one is the spin–statistics sentence this Part owes.

*What went in.* One input: the coin is the event’s — II.1’s noise read once more, as III.15b read it; from it, the rule “one process per state” (B). Derived from the rule read for every state, by polarization: the odd combination, the \\(-1\\) under exchange (D) — the free step keeps the sign and does not choose it, and read in one basis the rule admits an even sector (D), in \\(1{+}1\\) and \\(3{+}1\\) (E); and that the two natural readings of the rule on the step are not the free walker’s statistics (C), with the numbers. Not here: spin–statistics as a theorem; bosons beyond Part II’s walker.

(`scripts/pauli_test.py` (A), (B); `scripts/pauli_from_one_process_per_point.py` (C) first part, (D); `scripts/pauli_exclusion_every_state.py` (D), the symmetric column and the polarization; `scripts/kinematic_collision_vs_fermion.py` (C) second part.)

### VI.1a One zigzag, one cut, a second ribbon — by hand

*The picture first.* One path of III.1 with two flips, drawn on the page; a cut at one time; then a second ribbon through the same point. Read on the cut, the path is a particle, a pair appearing and a pair disappearing — the Z — with a weight that has no sign in it. The second ribbon is where exclusion touches the Z: a pair created into a state already occupied has weight zero, the direct term cancelled by its exchange. Every line is one step of the rules already on the page.

**What is used, and nothing else.**

- III.1, the process: clock \\(x\\), state \\((t,\\varepsilon)\\). In each step \\(\\Delta x\\): advance \\(t\\to t+\\varepsilon\\,\\Delta x/c\\), then a coin: flip \\(\\varepsilon\\to-\\varepsilon\\) with probability \\(\\frac ac\\Delta x\\), keep with \\(1-\\frac ac\\Delta x\\).
- III.5a (F), the loss taken out: on \\(Q=e^{ax/c}P\\) a keep weighs \\(1\\) and a flip weighs \\(+\\frac ac\\Delta x\\), **both** flips, \\(+\\to-\\) and \\(-\\to+\\) (the exchange matrix \\(\\begin{pmatrix}0&a\\\\a&0\\end{pmatrix}\\) has the same sign in both corners). To the order the document keeps (III.2: “a difference of \\(O(\\Delta x^2)\\) that the limit does not see”).
- III.3: on a time cut the forward-time observer counts crossings; a forward crossing carries \\(+1\\), a backward crossing \\(-1\\); the sum (the charge) is conserved in \\(t\\), the number is not; a trough is a pair appearing, a crest a pair disappearing.
- III.8c (A) and III.15a: the measurement. Routes of the walk *read in \\(t\\)*; weight of the event at endpoint \\(b\\) = (sum of outward routes to \\(b\\)) × (sum of return routes from \\(b\\)) \\(=\\alpha_b^2\\); the return is the same walk run backward, \\(U^{\\mathsf T}\\), the same coin entries read in the other order.
- III.15b: a record is data on a time cut. “In the clock \\(x\\) a detector at a place is passed by every path exactly once and its record is \\(\\int(P\_++P\_-)\\,dt=1\\), linear, with no square.”
- VI.1 (B): a point of spacetime carries at most one process in each direction. VI.1 (D): two walkers are a column \\(C\_{IJ}\\) over pairs of states, \\(C^{\\mathsf T}=-C\\), free step \\(C\\to MCM^{\\mathsf T}\\), weight \\((D-E)^2\\). VI.1 (A): the closed loop of one walker weighs \\((U^{\\mathsf T}U)\_{AA}=+1\\).

Units: \\(\\Delta x=1\\), \\(c=1\\), so one step moves \\(x\\) by \\(1\\) and \\(t\\) by \\(\\pm1\\). Write \\(\\beta\\equiv\\frac ac\\Delta x\\) for the flip weight.

**(A) One zigzag between two chosen points.**

**Picture.** Source \\(A=(x,t,\\varepsilon)=(0,0,+)\\). Endpoint \\(B=(5,3,+)\\). The path: senses \\(\\varepsilon_0\\ldots\\varepsilon_4=(+,+,-,+,+)\\) for the five steps.

```
 t
 3 |                          B
 2 |          o           o
 1 |      o       o
 0 |  A
   +--0---1---2---3---4---5--- x
```

Step by step, each line one step of III.1 (advance by the current sense, then the coin that sets the next sense):

| step | from \\((x,t)\\) | sense carried | to \\((x,t)\\) | coin after the step | weight (loss out) |
| --- | --- | --- | --- | --- | --- |
| 1 | \\((0,0)\\) | \\(+\\) | \\((1,1)\\) | keep \\(+\\) | \\(1\\) |
| 2 | \\((1,1)\\) | \\(+\\) | \\((2,2)\\) | **flip** to \\(-\\) | \\(\\beta\\) |
| 3 | \\((2,2)\\) | \\(-\\) | \\((3,1)\\) | **flip** to \\(+\\) | \\(\\beta\\) |
| 4 | \\((3,1)\\) | \\(+\\) | \\((4,2)\\) | keep \\(+\\) | \\(1\\) |
| 5 | \\((4,2)\\) | \\(+\\) | \\((5,3)\\) | keep \\(+\\) (arrival sense at \\(B\\)) | \\(1\\) |

**Its signed weight.**

\\\[w_1=\\underbrace{1}\_{\\text{keep}}\\cdot\\underbrace{\\beta}\_{+\\to-}\\cdot\\underbrace{\\beta}\_{-\\to+}\\cdot\\underbrace{1}\_{\\text{keep}}\\cdot\\underbrace{1}\_{\\text{keep}}=+\\beta^2 .\\\]

Every factor is positive. In the walk’s own clock there is no sign anywhere: III.5a (A)’s exchange matrix is \\(+a\\) in both corners, and III.1’s probabilities are probabilities. (The minus of III.5a (J) appears only when the same equations are *read in \\(t\\)*: it is “the backward sense running against \\(t\\), and nothing else”.)

**Where \\(B\\) sits.** \\(t_B=3\<x_B=5\\): \\(B\\) is outside the forward light cone of \\(A\\) — spacelike from it. Every path of III.1 has this property (III.14: “a zigzag in \\(t\\) with \\(x\\) as clock can only reach \\(|t|\\le x/c\\)”). The space-clock paths are the walk’s *spacelike* connections.

**The cut.** Cut at \\(t=\\tfrac32\\). The path crosses it three times:

| crossing | between | \\(x\\) at the cut | direction of \\(t\\) | III.3’s charge |
| --- | --- | --- | --- | --- |
| 1 | steps 1–2, i.e. on step 2 | \\(1.5\\) | up | \\(+1\\) |
| 2 | on step 3 | \\(2.5\\) | **down** | \\(-1\\) |
| 3 | on step 4 | \\(3.5\\) | up | \\(+1\\) |

Charge on the cut: \\(+1-1+1=+1\\) = the charge of \\(A\\) (III.3, conserved). Number of crossings: \\(3\\) (III.3, not conserved).

**The Z, and which steps made it.** Read upward in \\(t\\), as the forward-time observer reads:

- at \\(t=1\\), \\(x=3\\) (the **trough**, where step 3 ends and step 4 begins — the coin “flip to \\(+\\)” after step 3): a pair appears; one member goes right and up (step 4, 5: a particle), the other is step 3 read upward: \\(x\\) decreasing as \\(t\\) increases — a left-moving crossing with charge \\(-1\\), the antiparticle;
- at \\(t=2\\), \\(x=2\\) (the **crest**, where step 2 ends — the coin “flip to \\(-\\)” after step 2): the particle of steps 1–2, moving right, meets that antiparticle, moving left; they disappear.

So one path of III.1 with two flips, read on time cuts, is: particle in from \\(A\\), pair created at \\((3,1)\\), annihilation at \\((2,2)\\), particle out to \\(B\\). That is the Z. Its weight is the two flips, \\(\\beta^2\\), no sign; and *in \\(t\\) the creation comes before the annihilation* (\\(t=1\\) before \\(t=2\\)), the reverse of their order along the path (\\(x=2\\) before \\(x=3\\)).

**(B) Its measurement return: what the document defines, exactly.**

The measurement of III.8c (A)/III.15a is built on **routes read in \\(t\\)**: one position per time, the walker’s sense is a direction of motion, out from the source and back from the detector by the same coin in the other order, weight \\(\\alpha_b^2\\). The space-clock path of (A) is not such a route (it has three positions at \\(t=\\tfrac32\\)). So, taking the document as written:

**(a) In the walk’s own clock there is no return.** III.15b says it: a detector at a place \\(x_d\\) is crossed by every path exactly once (the path always advances in \\(x\\)), its record is \\(\\int(P\_++P\_-)dt=1\\), linear, no square — nothing to pair, nothing to return. The space-clock zigzag has a weight (\\(\\beta^2\\)), and it has no measurement return of its own.

**(b) On a time cut, each crossing is an endpoint for the \\(t\\)-reading’s measurement — and the \\(t\\)-reading’s routes do not reach two of the three.** A route read in \\(t\\) from \\(A=(0,0)\\) moves one site per tick, so at \\(t=\\tfrac32\\) it can only be at \\(|x|\\le\\tfrac32\\). Crossing 1 at \\(x=1.5\\) is on the cone; crossings 2 and 3, at \\(x=2.5\\) and \\(3.5\\), are outside it. Written with III.8c (A)’s symbols, for \\(b=(2.5,\\tfrac32,-)\\) or \\((3.5,\\tfrac32,+)\\):

\\\[\\alpha_b=e_b^{\\mathsf T}U(\\tfrac32,0)\\,e_A=0\\qquad\\Rightarrow\\qquad W_b=\\alpha_b^2=0 .\\\]

The \\(t\\)-reading’s outward sum is zero there, so its measurement rule, out-and-return, gives those two crossings weight zero from the source \\(A\\), at any finite step.

**(c) Where the two readings meet.** III.15 states that the \\(t\\)-reading’s sum over routes (inside the cone) and the space-clock sum over zigzags (outside the cone) are *one function across the light line*, as functions — the amplitude’s Bessel \\(J\\) inside continuing into the kernel’s Bessel \\(I\\) outside. That continuation is a statement about the two solved kernels, not a step rule that carries a route across the line. So the Z’s two outer crossings enter the measurement only through that continuation.

**That is the exact standing of “the return as the thing that makes pairs”:** the document’s return is a route read in \\(t\\), and the Z is a zigzag read in \\(x\\); the document has no step that turns one into the other; it has III.15’s identity of the two kernels. No approximation was made to say this; it is read off III.8c (A), III.15, III.15b.

**(C) A second ribbon, and every crossing on the cut.**

**Picture.** Keep path 1. Add ribbon 2 from \\(A_2=(2,0,+)\\) with senses \\((+,+,+)\\): \\((2,0)\\to(3,1)\\to(4,2)\\to(5,3)\\).

```
 t
 3 |                          B  (both paths arrive here: see (iii))
 2 |          o           o
 1 |      o       o(3,1)      <- path 1's trough; ribbon 2 passes here going up
 0 |  A           A2
   +--0---1---2---3---4---5--- x
```

Ribbon 2’s weight: \\(w_2=1\\cdot1\\cdot1=1\\) (three keeps).

**Crossings on the cut \\(t=\\tfrac32\\)**, both paths together:

| \\(x\\) | from | direction | charge |
| --- | --- | --- | --- |
| \\(1.5\\) | path 1, step 2 | up | \\(+1\\) |
| \\(2.5\\) | path 1, step 3 | down | \\(-1\\) |
| \\(3.5\\) | path 1, step 4 | up | \\(+1\\) |
| \\(3.5\\) | ribbon 2, step 2 | up | \\(+1\\) |

Two crossings at the same place, same direction: \\(x=3.5\\), up. Follow them back: both come from the point \\((3,1)\\) with sense \\(+\\) — path 1 after its coin at \\((3,1)\\) (flip \\(-\\to+\\)), ribbon 2 after its coin at \\((3,1)\\) (keep \\(+\\)). Same point, same direction. VI.1 (B): **not allowed as two processes.** This is the place where exclusion touches the Z: the trough’s created particle and a ribbon already there.

**How exclusion combines the alternatives (VI.1 (D)), written out.** Look at the clock value \\(x=3\\), where both arrive. Just before the coin at \\(x=3\\), the two walkers are in the states

\\\[I=(t{=}1,\\ \\varepsilon{=}-)\\ \\ \\text{(path 1, arriving on step 3)},\\qquad J=(t{=}1,\\ \\varepsilon{=}+)\\ \\ \\text{(ribbon 2, arriving on step 1)} .\\\]

Different states (opposite senses), so both may be there (VI.1 (E): “a forward and a backward sense at one point — particle and antiparticle — may”). The two-walker column has, on this pair of states, the antisymmetric entries

\\\[C\_{IJ}=+w,\\qquad C\_{JI}=-w,\\qquad C\_{II}=C\_{JJ}=0,\\qquad w=\\underbrace{\\beta}\_{\\text{path 1's first flip}}\\cdot\\underbrace{1}\_{\\text{ribbon 2's keep}} .\\\]

The coin at \\(x=3\\) is one \\(2\\times2\\) matrix per walker on the sense, with the loss out (III.5a (F)):

\\\[M=\\begin{pmatrix}M\_{++}&M\_{+-}\\\\M\_{-+}&M\_{--}\\end{pmatrix}=\\begin{pmatrix}1&\\beta\\\\\\beta&1\\end{pmatrix}\\qquad(\\text{rows: new sense; columns: old sense}),\\\]

and the free two-walker step is \\(C\\to MCM^{\\mathsf T}\\) (VI.1 (D)), i.e. \\(C'\_{I'J'}=\\sum\_{IJ}M\_{I'I}\\,C\_{IJ}\\,M\_{J'J}\\). Ask for the entry on the **same state twice**, \\(I'=J'=(1,+)\\) — the configuration “path 1 flipped into \\(+\\) and ribbon 2 kept \\(+\\)”:

\\\[\\begin{aligned} C'\_{(1,+),(1,+)} &=M\_{+,-}\\,C\_{IJ}\\,M\_{+,+}\\;+\\;M\_{+,+}\\,C\_{JI}\\,M\_{+,-}\\;+\\;M\_{+,-}\\,C\_{II}\\,M\_{+,-}\\;+\\;M\_{+,+}\\,C\_{JJ}\\,M\_{+,+}\\\\ &=\\beta\\cdot(+w)\\cdot1\\;+\\;1\\cdot(-w)\\cdot\\beta\\;+\\;0\\;+\\;0\\\\ &=\\beta w-\\beta w=0 . \\end{aligned}\\\]

Line by line: the first term is \\(D\\) — slot 1 (path 1) flips into \\(+\\), slot 2 (ribbon 2) keeps \\(+\\); the second is \\(E\\) — the same two moves with the slots exchanged, carrying \\(C\_{JI}=-C\_{IJ}\\). They are the same history counted once in each order, and exclusion’s sign makes the pair vanish. **The history “the Z creates a particle into the state ribbon 2 occupies” has weight zero.** Its would-be weight, each ordering, was \\(\\beta w=\\beta^2\\cdot1\\): the Z’s two flips times ribbon 2’s keep. What exclusion subtracts is exactly \\(-\\beta^2\\), attached to the creation vertex.

The entries that survive the same coin:

\\\[C'\_{(1,-),(1,+)}=M\_{--}C\_{IJ}M\_{++}+M\_{-+}C\_{JI}M\_{+-}=w-\\beta^2w=(1-\\beta^2)\\,w,\\qquad C'\_{(1,+),(1,-)}=-(1-\\beta^2)\\,w .\\\]

Read as pictures: \\(w\\) is “both keep” — path 1 continues its backward segment past \\(x=3\\) (no trough at \\((3,1)\\); the trough moves elsewhere), ribbon 2 goes on up. \\(-\\beta^2w\\) is “both flip” — path 1’s trough at \\((3,1)\\) *and* ribbon 2 turning down at \\((3,1)\\) (read upward: an antiparticle from \\((4,0)\\) annihilates ribbon 2 at \\((3,1)\\) at the same point where path 1’s pair appears): the two lines cross at a point, and exclusion counts that crossing with the exchange sign. After \\(x=3\\) the configuration always holds one \\(+\\) and one \\(-\\) at \\(t=1\\) — never two of the same sense.

**(iii) The endpoint.** Both paths as drawn arrive at \\(B=(5,3,+)\\): same state again, and by the same computation the two-walker column has zero weight there. So the pair of paths of the picture is, as a two-walker configuration, weight zero — twice over. The alternatives that carry weight are those in which the two walkers are in different states at every \\(x\\): ribbon 2 in the other sense at \\((3,1)\\), or at another \\(t\\), or path 1 without the trough there.

**A near miss, for contrast.** Ribbon 2 from \\((2,2,+)\\) instead, arriving at \\((3,3)\\): it never shares a state with path 1; \\(D\\) and \\(E\\) are different histories with different weights, \\((D-E)^2\\ne0\\), and the Z of path 1 is untouched. Exclusion acts only on the coincidence.

**(D) A closed loop?**

None appears, and here is why from the steps. A closed worldline would have to come back in \\(x\\) and in \\(t\\). A path of III.1 never comes back in \\(x\\) (every step is \\(+\\Delta x\\)). A route read in \\(t\\) (III.15a) never comes back in \\(t\\). So no path of either reading closes, and no combination of the two is a step of either. The one closed object in the document is the **event**: out-leg and return leg, III.15a; it closes at the source and weighs \\((U^{\\mathsf T}U)\_{AA}=+1\\) (VI.1 (A)), with no sign on it.

What (C) *does* produce is a minus sign on a Z: when the trough creates into an occupied state, the exchange term \\(E\\) carries \\(-1\\) and cancels the direct term. That is VI.1’s sign at a creation vertex — the first and only place in the original rules where a minus attaches to pair creation.

**(E) What this settles and what it does not.**

Against the three rules that make Feynman’s counting:

| needed | in the original rules? | where |
| --- | --- | --- |
| backward segment carries \\(-q\\) | **yes** | III.3: a backward crossing on a time cut counts \\(-1\\); the charge is forward minus backward |
| which segments run backward (the Z) | **yes, as space-clock zigzags**, weight \\(+\\beta^2\\) per trough–crest pair, living outside the source’s cone | III.1, III.14; (A) |
| the Z inside the *measurement* (out-and-return) | **no step does it**; the \\(t\\)-reading’s routes stop at the cone, the Z lives outside; the two meet only in III.15’s identity of the two kernels | (B) |
| \\(-1\\) when a Z meets an occupied state | **yes, derived**: \\(D-E=0\\) at the coincidence, \\(-\\beta^2\\) subtracted at the creation vertex | VI.1 (D); (C) |
| \\(-1\\) per closed loop | no closed *path* exists in either reading (the only closed path-object is the event, weight \\(+1\\)); the \\(-1\\) per closed *history* lives in the configuration count: the cyclic exchange sign \\((-1)^{n-1}\\) of VI.1 (D) (VI.5 (A)) | VI.1 (A), (D); (D) |

So the status of the idea “the return makes the pairs” is this: the pairs are made by the space-clock flips (A), the return is a \\(t\\)-route (B), and the document has no step that joins them. The sign that exclusion gives (C) is derived, but it is a sign on a blocked creation, not on a loop; Feynman’s \\(-1\\) per loop is the same exchange sign on closed histories of the configuration count, and that count is also what selects which states are occupied when no ribbon is present (VI.5 (A), (C)).

(`scripts/one_zigzag_one_cut_two_ribbons.py`, exact count, not a solving method: it enumerates the \\(2^5\\) sense sequences of III.1 from \\(A\\) with \\(\\Delta x=1\\) and checks: path 1’s weight is \\(\\beta^2\\) and its three crossings carry \\(+1,-1,+1\\); every path from \\(A\\) that ends above the cut has charge \\(+1\\) on it and an odd number of crossings, every path that ends below has charge \\(0\\) and an even number (III.3: the charge on a cut is what has passed it); the two-walker column through the coin at \\(x=3\\) has the same-state entry \\(0\\) and the surviving entries \\(\\pm(1-\\beta^2)w\\), with the \\(D\\) and \\(E\\) terms \\(\\pm\\beta w\\) printed separately. Nothing is solved, summed to a limit, or approximated.)

## VI.2 The field from within: the angle of the pair is a connection, and only loops read it

*The picture first.* Part V found where the field enters the walk — as a turn of the walker in its own plane, the pair \\((Q\_\\varepsilon,Q'\_\\varepsilon)\\) of each time sense (III.8c (B)) — and left two lines open in its ledger: why the field produces that turn angle, and the field’s own law. This section starts from V.2’s step on the four detectors and from one rule, III.1’s own: the coin is drawn at the event and depends on nothing else. Written as an equation on the step, that rule says the position of the clock hand is a spectator of the coin, so the pair can be turned at every event without changing the walk (A). A turn that differs from event to event changes the lines by exactly a potential, and a potential that is nothing but such a turn has zero field (B). An event reads the pair only through closed loops, and a loop reads only the enclosed angle (C). The force on the mean step is the gradient of the turn angle per tick, which in three dimensions is the Lorentz force with the velocity the direction matrices (D). The carrier between two walkers is Part II’s walker at \\(m=0\\) — no flips, no direction, VI.1 (F)’s boson — and its source answer is Coulomb’s law (E). If the angle at the event is itself a noise of the event with a local quadratic weight, that weight is \\(F^2\\), forced by the two invariances already on the page (F). The rule between two walkers is then Maxwell solved for one walker’s current and read by the other’s segments, with the carrier’s state left as a named choice (F′); and the chain is run as one process on a lattice (G).

*Symbols.* \\(\\mathbf P\_\\varepsilon=(P\_\\varepsilon,P'\_\\varepsilon)^{\\mathsf T}\\) is the pair of counters of one sense, the original and the added detector (III.8c (B), Part V). \\(J=\\begin{pmatrix}0&-1\\\\1&0\\end{pmatrix}\\) is the quarter-turn of the pair and \\(R(\\theta)=\\cos\\theta\\,1+\\sin\\theta\\,J\\) its turn, with \\(dR/d\\theta=JR=RJ\\). \\(\\theta_i\\) is the angle accumulated along a path \\(i\\). Nothing else is new.

**(A) The coin does not read the clock hand.** V.2 (A)–(B): in one clock step \\(dx\\) the pair of sense \\(\\varepsilon\\) advances, flips with probability \\(a\\,dx/c\\), and is turned by \\(\\theta\_\\varepsilon\\):

\\\[\\mathbf P\_\\varepsilon(t;x+dx)=R(\\theta\_\\varepsilon)\\Big\[\\big(1-\\tfrac{a\\,dx}c\\big)\\mathbf P\_\\varepsilon\\big(t-\\varepsilon\\tfrac{dx}c;x\\big)+\\tfrac{a\\,dx}c\\,\\mathbf P\_{-\\varepsilon}(t;x)\\Big\].\\\]

III.1’s rule is that the flip probability is \\(a\\,dx/c\\) at the event, independent of everything earlier; where the hand stands is not among the event’s data, so the flip multiplies the pair by a number, not by a matrix on the pair. That is the whole content of “the noise belongs to the event” for one walker, and it is visible in the bracket: every term is a scalar times a pair. A scalar commutes with \\(J\\), so a constant turn \\(R(\\theta_0)\\) passes through the bracket:

\\\[R(\\theta_0)\\,\\mathbf P\_\\varepsilon(t;x+dx)=R(\\theta\_\\varepsilon)\\Big\[\\big(1-\\tfrac{a\\,dx}c\\big)R(\\theta_0)\\mathbf P\_\\varepsilon\\big(t-\\varepsilon\\tfrac{dx}c;x\\big)+\\tfrac{a\\,dx}c\\,R(\\theta_0)\\mathbf P\_{-\\varepsilon}(t;x)\\Big\]:\\\]

\\(R(\\theta_0)\\mathbf P\\) obeys the step \\(\\mathbf P\\) obeys. The walk never reads the absolute angle of the pair. On the lattice of III.15a, with \\(U=(R\_{\\rm coin}\\otimes1)\\,S\\), this is \\(JU=UJ\\) (`pair_angle_is_a_connection.py` 1). V.2 (G) observed this freedom; here it is III.1’s rule.

**(B) A turn that varies from event to event is a potential, and a potential that is only a turn has no field.** Taylor the step as in V.2 (C): with \\(\\theta\_\\varepsilon=-\\frac e\\hbar(\\varepsilon\\Phi\\frac{dx}c-A_x\\,dx)\\) and \\(R(\\theta)=1+\\theta J+O(\\theta^2)\\), the lines for the pair are

\\\[c\\,\\partial_x\\mathbf P\_\\varepsilon=-\\varepsilon\\,\\partial_t\\mathbf P\_\\varepsilon-a\\mathbf P\_\\varepsilon+a\\mathbf P\_{-\\varepsilon}-\\frac e\\hbar\\big(\\varepsilon\\Phi-cA_x\\big)J\\mathbf P\_\\varepsilon .\\tag{L}\\\]

(Written out, \\(-\\frac e\\hbar(\\varepsilon\\Phi-cA_x)J\\mathbf P\_\\varepsilon=\\frac e\\hbar(\\varepsilon\\Phi-cA_x)\\,(P'\_\\varepsilon,\\,-P\_\\varepsilon)^{\\mathsf T}\\), which is V.2 (C)’s pair of lines.) Now turn the pair by an angle that depends on the event, \\(\\tilde{\\mathbf P}\_\\varepsilon=R(\\theta)\\mathbf P\_\\varepsilon\\) with \\(\\theta=\\theta(t,x)\\), the same for both senses. Then \\(\\mathbf P\_\\varepsilon=R(-\\theta)\\tilde{\\mathbf P}\_\\varepsilon\\), and by the product rule with \\(dR/d\\theta=JR\\),

\\\[\\partial_t\\mathbf P\_\\varepsilon=R(-\\theta)\\,\\partial_t\\tilde{\\mathbf P}\_\\varepsilon-(\\partial_t\\theta)\\,J\\,R(-\\theta)\\tilde{\\mathbf P}\_\\varepsilon,\\qquad\\partial_x\\mathbf P\_\\varepsilon=R(-\\theta)\\,\\partial_x\\tilde{\\mathbf P}\_\\varepsilon-(\\partial_x\\theta)\\,J\\,R(-\\theta)\\tilde{\\mathbf P}\_\\varepsilon .\\\]

Put these into (L) and multiply the whole line on the left by \\(R(\\theta)\\); \\(R(\\theta)\\) commutes with \\(J\\) and with every scalar, and \\(R(\\theta)R(-\\theta)=1\\):

\\\[c\\,\\partial_x\\tilde{\\mathbf P}\_\\varepsilon-c\\,(\\partial_x\\theta)J\\tilde{\\mathbf P}\_\\varepsilon=-\\varepsilon\\,\\partial_t\\tilde{\\mathbf P}\_\\varepsilon+\\varepsilon(\\partial_t\\theta)J\\tilde{\\mathbf P}\_\\varepsilon-a\\tilde{\\mathbf P}\_\\varepsilon+a\\tilde{\\mathbf P}\_{-\\varepsilon}-\\frac e\\hbar\\big(\\varepsilon\\Phi-cA_x\\big)J\\tilde{\\mathbf P}\_\\varepsilon .\\\]

Move the two \\(\\theta\\) terms to the right and collect what multiplies \\(J\\tilde{\\mathbf P}\_\\varepsilon\\):

\\\[c\\,\\partial_x\\tilde{\\mathbf P}\_\\varepsilon=-\\varepsilon\\,\\partial_t\\tilde{\\mathbf P}\_\\varepsilon-a\\tilde{\\mathbf P}\_\\varepsilon+a\\tilde{\\mathbf P}\_{-\\varepsilon}-\\frac e\\hbar\\Big\[\\varepsilon\\Big(\\Phi-\\frac\\hbar e\\partial_t\\theta\\Big)-c\\Big(A_x+\\frac\\hbar e\\partial_x\\theta\\Big)\\Big\]J\\tilde{\\mathbf P}\_\\varepsilon .\\\]

This is (L) again, with

\\\[\\tilde\\Phi=\\Phi-\\frac\\hbar e\\,\\partial_t\\theta,\\qquad\\tilde A_x=A_x+\\frac\\hbar e\\,\\partial_x\\theta .\\\]

Two readings of one equality. Read with \\(\\Phi,A_x\\) given: V.2 (G), the potentials change by a gradient when the clock-zero is moved, now line by line. Read with \\(\\Phi=A_x=0\\): a walker whose pair is turned by \\(\\theta(t,x)\\) obeys the lines of a walker in the potentials \\(\\tilde\\Phi=-\\frac\\hbar e\\partial_t\\theta\\), \\(\\tilde A_x=\\frac\\hbar e\\partial_x\\theta\\), and for those

\\\[\\tilde E=-\\partial_x\\tilde\\Phi-\\partial_t\\tilde A_x=\\frac\\hbar e\\,\\partial_x\\partial_t\\theta-\\frac\\hbar e\\,\\partial_t\\partial_x\\theta=0 .\\\]

So the lines fix where the field sits — it can only be the angle on the link, nothing else in (L) acts on the pair — and they fix what is not a field: any angle that is the difference of clock-zeros at the two ends of the link. A field is the part of the link angle that no choice of clock-zeros removes. On the lattice: \\(U\[\\theta\_{\\rm link}+\\theta\_{x+e}-\\theta_x\]=R(\\theta)\\,U\[\\theta\_{\\rm link}\]\\,R(\\theta)^{\\mathsf T}\\) exactly, and the event weights \\(Q^2+Q'^2\\) are unchanged (`pair_angle_is_a_connection.py` 2, 3, 3\\('\\)). In the clock reading \\(\\theta\_\\varepsilon\\) carries \\(\\varepsilon\\), so the backward sense is turned the other way by the same \\(\\Phi\\): V.2 (H).

**(C) Only loops read it.** III.15a (A): the reading at a point is the sum over paths of their signed weights, \\(Q\_\\varepsilon=\\sum_iw_i\\). With link angles the pair arriving along path \\(i\\) is turned by the sum of the link angles along \\(i\\), \\(\\theta_i\\); starting from the pair \\((1,0)^{\\mathsf T}\\) at the emitter,

\\\[\\begin{pmatrix}Q\_\\varepsilon\\\\Q'\_\\varepsilon\\end{pmatrix}=\\sum_iw_i\\,R(\\theta_i)\\begin{pmatrix}1\\\\0\\end{pmatrix}=\\sum_iw_i\\begin{pmatrix}\\cos\\theta_i\\\\\\sin\\theta_i\\end{pmatrix}.\\\]

The event weight is the length of the pair (V.5 (B)):

\\\[W=Q\_\\varepsilon^2+Q'^2\_\\varepsilon=\\Big(\\sum_iw_i\\cos\\theta_i\\Big)^2+\\Big(\\sum_iw_i\\sin\\theta_i\\Big)^2=\\sum\_{i,j}w_iw_j\\big(\\cos\\theta_i\\cos\\theta_j+\\sin\\theta_i\\sin\\theta_j\\big)=\\sum\_{i,j}w_iw_j\\cos(\\theta_i-\\theta_j).\\\]

Under the relabelling of (B) every path from the emitter \\(A\\) to the absorber \\(X\\) changes by the same amount, \\(\\theta_i\\to\\theta_i+\\theta(X)-\\theta(A)\\), so every difference \\(\\theta_i-\\theta_j\\) is unchanged. And \\(\\theta_i-\\theta_j\\) is the sum of the link angles out along \\(i\\) and back along \\(j\\) — the angle around the loop the two paths close. The return leg is the transposed step, which carries \\(-\\theta\_{\\rm link}\\) on each link (III.15a (B) with the turn: \\(R(\\theta)^{\\mathsf T}=R(-\\theta)\\)); out and back along one path is the free loop, \\((U^{\\mathsf T}U)\_{AA}=1\\) for any link field (5a). A link field whose loop sums all vanish gives the free weights to machine precision; one link carrying an angle changes them; a line that both crossings pick up with the same sign is a flux, and the weights are periodic in it with period \\(2\\pi\\) (4, 5b): the Aharonov–Bohm effect as the statement that an event is a loop. With the pairing \\(\\psi=Q+iQ'\\) of V.2 (F) the same lines read \\(\\psi=\\sum_iw_ie^{i\\theta_i}\\), \\(W=|\\psi|^2\\), and \\(\\theta\_{\\rm link}=\\frac e\\hbar A\_\\mu dx^\\mu\\) is the phase added to the walker’s own phase per step: Part II’s drift shift \\(mv^\\mu=\\partial^\\mu S-eA^\\mu\\) (V.3) as a one-form.

**(D) The force is the gradient of the turn angle per tick.** In the \\(t\\)-reading one tick turns the pair at \\(x\\) by \\(\\theta(x)=-\\frac e\\hbar\\Phi(x)\\,dt\\) (V.2 (E)). Take the momentum density of the pair, the cross term \\(\\hbar\\,(Q\\,\\partial_xQ'-Q'\\partial_xQ)=-\\hbar\\,\\mathbf Q^{\\mathsf T}J\\,\\partial_x\\mathbf Q\\) (in V.2 (F)’s pairing, \\(\\hbar\\,\\mathrm{Im}\\,\\psi^\*\\partial_x\\psi\\)). After the turn, with \\(\\partial_x(R\\mathbf Q)=R\\,\\partial_x\\mathbf Q+(\\partial_x\\theta)JR\\,\\mathbf Q\\), \\(R^{\\mathsf T}JR=J\\) and \\(J^2=-1\\):

\\\[(R\\mathbf Q)^{\\mathsf T}J\\,\\partial_x(R\\mathbf Q)=\\mathbf Q^{\\mathsf T}R^{\\mathsf T}JR\\,\\partial_x\\mathbf Q+(\\partial_x\\theta)\\,\\mathbf Q^{\\mathsf T}R^{\\mathsf T}JJR\\,\\mathbf Q=\\mathbf Q^{\\mathsf T}J\\,\\partial_x\\mathbf Q-(\\partial_x\\theta)\\,|\\mathbf Q|^2 ,\\\]

so the momentum density gains \\(\\hbar\\,(\\partial_x\\theta)\\,|\\mathbf Q|^2\\). Integrated over \\(x\\), with \\(\\theta=-\\frac e\\hbar\\Phi\\,dt\\), the mean momentum changes per tick by

\\\[\\Delta\\langle p\\rangle=\\hbar\\int(\\partial_x\\theta)\\,|\\mathbf Q|^2dx=e\\int(-\\partial_x\\Phi)\\,|\\mathbf Q|^2dx\\;dt=e\\,\\langle E\\rangle\\,dt .\\\]

The force is the gradient of the angle the event turns the pair by, per unit time. In three dimensions the same statement with the direction matrices: with \\(\\boldsymbol\\pi=\\hat{\\mathbf p}-e\\mathbf A\\), \\(H=c\\,\\boldsymbol\\alpha\\cdot\\boldsymbol\\pi+\\beta mc^2+e\\Phi\\),

\\\[\\frac{d\\pi_i}{dt}=\\frac i\\hbar\[H,\\pi_i\]+\\frac{\\partial\\pi_i}{\\partial t},\\qquad\[\\pi_j,\\pi_i\]=i\\hbar e\\,(\\partial_jA_i-\\partial_iA_j)=i\\hbar e\\,\\epsilon\_{jik}B_k,\\qquad\[e\\Phi,\\pi_i\]=i\\hbar e\\,\\partial_i\\Phi ,\\\]

\\\[\\frac i\\hbar\\,c\\,\\alpha_j\[\\pi_j,\\pi_i\]=\\frac i\\hbar\\,c\\,\\alpha_j\\,i\\hbar e\\,\\epsilon\_{jik}B_k=e\\,c\\,\\epsilon\_{ijk}\\alpha_jB_k=e\\,(c\\,\\boldsymbol\\alpha\\times\\mathbf B)\_i,\\qquad\\frac i\\hbar\[e\\Phi,\\pi_i\]-e\\,\\partial_tA_i=-e\\,\\partial_i\\Phi-e\\,\\partial_tA_i=e\\,E_i ,\\\]

\\\[\\boxed{\\frac{d\\boldsymbol\\pi}{dt}=e\\,\\big(\\mathbf E+c\\,\\boldsymbol\\alpha\\times\\mathbf B\\big)}\\\]

exactly, for any \\(\\Phi(t,\\mathbf x)\\), \\(\\mathbf A(t,\\mathbf x)\\), on every column (`lorentz_force_turn_angle.py` 1, 1\\('\\), 1\\(''\\)). \\(\[\\pi_j,\\pi_i\]\\) is the angle of a plaquette, the two orders of two moves (V.3); \\(\[e\\Phi,\\pi_i\]\\) is the gradient of the time-link angle. On the lattice, with the two real legs, a packet in a uniform \\(E\\) gains momentum at the rate \\(eE\\) to \\(0.3\\%\\) (2). This closes Part V’s ledger line: the turn angle is not taken from a phase; a phase per step is what a connection is, and its gradient per tick is the force.

*How a curvature becomes a push, and the three forces are one shape.* Two things the chain above says without saying. First, the momentum of a walker is the twist of its pair along \\(x\\) — the cross term \\(\\mathbf Q^{\\mathsf T}J\\partial_x\\mathbf Q\\) is how much the angle of the pair changes from site to site — and “the turn per tick differs between \\(x\\) and \\(x+dx\\)” is a loop angle: a step at \\(t\\), a tick at \\(x+dx\\), the step back at \\(t+dt\\), the tick back at \\(x\\),

\\\[\\theta\_{\\square}=\\theta_x(t)+\\theta_t(x+dx)-\\theta_x(t+dt)-\\theta_t(x)=\\frac e\\hbar\\big(-\\partial_tA_x-\\partial_x\\Phi\\big)\\,dt\\,dx=\\frac e\\hbar\\,E_x\\,dt\\,dx ,\\\]

the curvature of the pair angle on the time–space plaquette, read exactly as \\(B\\) is read on a space plaquette (V.3, `lorentz_force_turn_angle.py` 1\\('\\)). The force of (D) is this loop angle per tick, divided by \\(dx\\): the pair comes back from the rectangle twisted, and a twist is a momentum. Second, V.4’s effective flip rate makes the same statement in one line. At rest in a potential the pair at \\(x\\) advances by \\(a\_{\\text{eff}}(x)\\,dt=\\big(a+\\tfrac e\\hbar\\Phi(x)\\big)dt\\) per tick, and

\\\[F=-\\partial_x\\big(\\hbar\\,a\_{\\text{eff}}\\big)=-e\\,\\partial_x\\Phi=eE :\\\]

the force is the gradient of the total turn rate, and the coin \\(a\\), being the same everywhere, drops out of the gradient — the field pushes without touching the flip rate (V.4). On the lattice the two are measured separately: the phase the pair advances per tick, site by site, is \\(a+\\tfrac e\\hbar\\Phi(x)\\) with slope \\(-eE\\) to \\(2\\times10^{-4}\\) for a packet \\(e^{-x^2/2w^2}\\) fifty Compton lengths wide — the difference is the packet’s own \\(\\tfrac12(c/aw)^2\\), four times larger at half the width — and its gradient equals the measured \\(d\\langle p\\rangle/dt\\) to \\(0.02\\%\\), the force itself being \\(eE\\) to \\(3\\times10^{-7}\\). The local rate has to be read over whole periods of the beat at \\(2a\\), or on a packet of one sign of energy: over \\(40\\) ticks of an unprojected packet it reads \\(80\\%\\) off, the beat’s phase sloping with \\(x\\) (`force_is_a_turn_rate_gradient.py` 1, 2). This is the shape the document’s other two forces have. V.5 (C): the detector is pushed by \\(-\\partial\_{x_D}(\\hbar\\,\\Delta a_D)\\), the gradient of its own clock shift \\(g(Q^2+Q'^2)(x_D)\\). The Casimir document: the walls are pushed by \\(-\\partial_LE\_{\\text{vac}}\\), \\(E\_{\\text{vac}}=-\\tfrac12\\operatorname{Tr}|H|\\) of the carrier, the gradient of the vacuum’s total turn rate. Three forces, one shape — a total turn rate differentiated with respect to a position — and three things that set the rate: for the Lorentz force a second walker’s charge, carried to \\(x\\) by the carrier of (E); for the detector, the particle on its site; for Casimir, a wall. No new process: the link angle of (F) is the only thing added, and (D) is III.1’s rule on it.

**(E) The carrier.** Two walkers feel each other through the link angles, and the angles must be carried from one event to another. The Casimir document: “the long-range force needs a carrier that never flips”. Part II’s walker with \\(m=0\\) has \\(K=0\\) in II.7, so its phase does not advance with the step label — it never flips — and it has no direction at a point, so by VI.1 (F) it is a boson. What it carries is a link angle, \\(a\_\\mu dx^\\mu\\), so the walker carries a one-form index: its kernel is \\(\\eta\_{\\mu\\nu}\\) times the scalar one, and the time component comes with \\(\\eta\_{00}=-1\\) — II.1’s sign, once more. That sign decides the sense of the force: a scalar carrier makes like sources attract, a one-form carrier makes them repel (static field energy on a grid, both signs, against \\(\\mp1/4\\pi r\\): `carrier_sign_and_current.py` 1a–1c). The scalar part of the kernel is III.16’s \\(\\tau\\)-sum of II.2’s Gaussian with \\(e^{-m^2\\tau}\\to1\\). With \\(u=s^2/4\\tau\\), \\(d\\tau=-\\frac{s^2}{4u^2}du\\), \\(\\frac1{(4\\pi\\tau)^2}=\\frac{u^2}{\\pi^2s^4}\\):

\\\[\\int_0^\\infty\\frac{d\\tau}{(4\\pi\\tau)^2}\\,e^{-s^2/4\\tau}=\\int_0^\\infty\\frac{u^2}{\\pi^2s^4}\\cdot\\frac{s^2}{4u^2}\\,e^{-u}\\,du=\\frac1{4\\pi^2s^2}\\int_0^\\infty e^{-u}du=\\frac1{4\\pi^2s^2},\\\]

the \\(m\\to0\\) limit of \\(\\frac{m}{4\\pi^2s}K_1(ms)\\) (`carrier_and_field_action.py` 1). For a source that sits still, \\(s^2=r^2+t^2\\) and the sum over its time is

\\\[\\int\_{-\\infty}^{\\infty}\\frac{dt}{4\\pi^2(r^2+t^2)}=\\frac1{4\\pi^2}\\Big\[\\frac1r\\arctan\\frac tr\\Big\]\_{-\\infty}^{\\infty}=\\frac1{4\\pi^2}\\cdot\\frac\\pi r=\\frac1{4\\pi r}:\\\]

Coulomb’s law (2), with the sign of the one-form: like walkers repel. In \\(1{+}1\\) the same sum is \\(-\\frac1{2\\pi}\\ln s\\) and its static field is \\(\\int dt\\,\\partial_xG=-\\tfrac12\\,\\mathrm{sgn}\\,x\\), a potential \\(|x|/2\\) that grows without bound (3). And \\(1/(4\\pi^2s^2)\\) is the transform of \\(1/k^2\\) in four dimensions (4): the carrier’s kernel, with its \\(\\eta\_{\\mu\\nu}\\), is the propagator of (F) (Feynman gauge), so “the carrier is Part II’s walker at \\(m=0\\) with a one-form index” and “the field’s weight is \\(F^2\\)” are one statement at the free level. One carrier is not another: the companion document’s coin calculation uses the walker’s *own* no-flip leg between walls, which is the fermion field’s Casimir effect; the electromagnetic force between charged walkers is carried by this walker, Part II’s at \\(m=0\\), with no sense and a one-form index.

**(F) The field’s own law.** Let the link angle be a noise of the event, as Part II’s step was: \\(a\_\\mu(x)\\) drawn with a weight \\(e^{-\\int\\mathcal L\[a\]}\\), \\(\\mathcal L\\) local and quadratic. II.1 fixed the walker’s variance by one demand; here there are two, both already on the page:

\\\[\\mathcal L\[\\Lambda a\]=\\mathcal L\[a\]\\quad(\\text{no preferred frame, II.1}),\\qquad\\mathcal L\[a+\\partial\\chi\]=\\mathcal L\[a\]+\\partial\_\\mu(\\ldots)\\quad(\\text{the relabelling of (B)}).\\\]

With two derivatives the frame-invariant scalars are

\\\[T_0=a\_\\mu a^\\mu,\\qquad T_1=\\partial\_\\mu a\_\\nu\\,\\partial^\\mu a^\\nu,\\qquad T_2=\\partial\_\\mu a\_\\nu\\,\\partial^\\nu a^\\mu,\\qquad T_3=(\\partial\_\\mu a^\\mu)^2 ,\\\]

and \\(\\mathcal L=c_0T_0+c_1T_1+c_2T_2+c_3T_3\\). Its equations of motion are linear in \\(a\\); under \\(a\\to a+\\partial\\chi\\) they shift by terms in \\(\\chi\\) alone, and with \\(\\chi=e^{ik\\cdot x}\\) the shift is a polynomial in \\(k\\) times \\(e^{ik\\cdot x}\\) that must vanish for every \\(k\\):

\\\[c_0=0,\\qquad c_1+c_2+c_3=0\\\]

(5). \\(T_2\\) and \\(T_3\\) give the same equations, since

\\\[\\partial\_\\mu a\_\\nu\\,\\partial^\\nu a^\\mu=\\partial\_\\mu\\big(a\_\\nu\\partial^\\nu a^\\mu\\big)-\\partial^\\nu\\big(a\_\\nu\\partial\_\\mu a^\\mu\\big)+(\\partial\_\\mu a^\\mu)^2 ,\\\]

two total derivatives apart. So \\(\\mathcal L=c_1(T_1-T_2)\\) up to total derivatives, and

\\\[T_1-T_2=\\partial\_\\mu a\_\\nu\\,\\partial^\\mu a^\\nu-\\partial\_\\mu a\_\\nu\\,\\partial^\\nu a^\\mu=\\tfrac12\\big(\\partial\_\\mu a\_\\nu-\\partial\_\\nu a\_\\mu\\big)\\big(\\partial^\\mu a^\\nu-\\partial^\\nu a^\\mu\\big)=\\tfrac12F\_{\\mu\\nu}F^{\\mu\\nu}\\\]

(5\\('\\)): Maxwell’s free weight, unique up to its scale. This step is textbook — a Lorentz-invariant, gauge-invariant, quadratic weight for a one-form is \\(F^2\\) (Weinberg, *The Quantum Theory of Fields* I, ch. 8; the massless spin-one argument) — and is cited, not claimed — though the step is short enough to show: a translation-invariant quadratic weight for \\(a\_\\mu\\) with at most two derivatives has the Fourier kernel \\(K^{\\mu\\nu}(k)=c_0\\eta^{\\mu\\nu}+c_1k^2\\eta^{\\mu\\nu}+c_2k^\\mu k^\\nu\\); the relabelling \\(a\_\\mu\\to a\_\\mu+\\partial\_\\mu\\lambda\\) of (A) requires \\(K^{\\mu\\nu}k\_\\nu=0\\) for every \\(k\\), so \\(c_0=0\\) and \\(c_2=-c_1\\), which is \\(F\_{\\mu\\nu}F^{\\mu\\nu}\\) up to scale and total derivatives (\\(F\\tilde F\\) is one). The two-derivative restriction is the premise; with more derivatives there are more invariants. What is this document’s is the origin of the gauge invariance, (A): the coin does not read the clock hand. Two things to keep straight about the weight. It is a weight in the signature where Part II’s walk is real (III.16), where \\(F^2\\) is a sum of squares; in Lorentzian signature \\(F^2=2(B^2-E^2)\\) is indefinite and \\(e^{-\\int\\mathcal L}\\) is not a weight. And its source must be conserved, \\(\\partial\_\\nu j^\\nu=0\\), or \\(\\partial\_\\mu F^{\\mu\\nu}=j^\\nu\\) is inconsistent (\\(\\partial\_\\nu\\partial\_\\mu F^{\\mu\\nu}=0\\)): the source is the \\(t\\)-reading’s content, \\(\\rho=\\psi^\\dagger\\psi\\), \\(j=c\\,\\psi^\\dagger\\alpha\\psi\\), and \\(\\partial_t\\rho+\\partial_xj=0\\) holds identically for any \\(\\Phi(t,x)\\), \\(A_x(t,x)\\) in V.2 (F)’s equation (`carrier_sign_and_current.py` 2a; on the lattice, V.2’s run on the four detectors with a potential step keeps the total length to \\(10^{-15}\\)). This is IV.4d (B)’s conservation with the field on. III.3’s conserved pair \\((P\_+-P\_-,\\ c(P\_++P\_-))\\) is the walk’s own clock reading — a different object, also conserved; the source of Maxwell is the square, as everything a detector reads is. With the walker’s turn angle as the source term, \\(-\\frac14F^2+j^\\mu a\_\\mu\\) — the plus sign belongs to this document’s signature \\((-,+)\\), where \\(j^\\mu a\_\\mu=-\\rho\\Phi+\\mathbf j\\cdot\\mathbf A\\); with \\(-j^\\mu a\_\\mu\\) one gets \\(\\nabla\\cdot\\mathbf E=-\\rho\\) and like charges attract — the equations of motion are

\\\[\\partial\_\\mu F^{\\nu\\mu}=j^\\nu ,\\qquad\\text{i.e.}\\qquad\\nabla\\cdot\\mathbf E=\\rho,\\quad\\nabla\\times\\mathbf B-\\partial_t\\mathbf E=\\mathbf j,\\\]

with \\(j^\\nu\\) what multiplies \\(a\_\\nu\\) in the walker’s energy — V.2 (F)’s \\(e\\Phi\\,\\psi^\\dagger\\psi-ec\\,\\psi^\\dagger\\boldsymbol\\alpha\\cdot\\mathbf A\\,\\psi\\), i.e. \\(j^\\nu=e\\,(\\rho,\\ c\\,\\psi^\\dagger\\boldsymbol\\alpha\\psi)\\), the current of IV.7 (6; the check is strict about the sign, \\(\\nabla\\cdot\\mathbf E=+\\rho\\)).

**(F′) The rule between two walkers is Maxwell, solved and read.** *The picture first.* Walker 1 moves; every one of its segments sends the carrier out; the carrier reaches a segment of walker 2 and turns its pair. Summed over walker 1’s segments that is a potential at walker 2, and (D) says what the potential does. This paragraph writes the rule and shows it is (F)’s Maxwell with nothing added, then says the one thing Maxwell leaves open.

*The rule.* With \\(c=1\\), \\(\\eta=\\mathrm{diag}(-1,+1,+1,+1)\\), a carrier path from a segment \\(dx_1^\\mu\\) of walker 1 to a segment \\(dx_2^\\nu\\) of walker 2 turns walker 2’s pair by

\\\[\\delta\\theta_2=\\frac{q_1q_2}\\hbar\\,\\eta\_{\\mu\\nu}\\,dx_1^\\mu\\,dx_2^\\nu\\;\\mathcal G(x_2,x_1),\\qquad\\eta\_{\\mu\\nu}dx_1^\\mu dx_2^\\nu=-dt_1dt_2+d\\mathbf x_1\\!\\cdot d\\mathbf x_2 ,\\tag{1}\\\]

1. \\(\\mathcal G\\) the carrier’s kernel: the sum over all carrier paths from \\(x_1\\) to \\(x_2\\) of their weights. Walker 1’s segments, weighted by its column, are its time segments with density \\(\\rho_1\\) and its moves with current \\(\\mathbf j_1=\\psi_1^\\dagger\\boldsymbol\\alpha\\psi_1\\) (IV.4d (B)), so

\\\[\\Phi(x_2)=q_1\\!\\int\\!d^4x_1\\,\\mathcal G(x_2,x_1)\\,\\rho_1(x_1),\\qquad\\mathbf A(x_2)=q_1\\!\\int\\!d^4x_1\\,\\mathcal G(x_2,x_1)\\,\\mathbf j_1(x_1),\\qquad\\delta\\theta_2=\\frac{q_2}\\hbar\\big(\\mathbf A\\cdot d\\mathbf x_2-\\Phi\\,dt_2\\big),\\tag{2}\\\]

2. V.2’s rule, with the minus of \\(-\\Phi\\,dt_2\\) supplied by \\(\\eta\_{00}\\): that sign is what makes like charges repel and parallel currents attract ((G), check 3–4).

*It is Maxwell.* (F) gives \\(\\partial\_\\mu F^{\\nu\\mu}=j^\\nu\\) with \\(j^\\nu=q(\\rho,\\psi^\\dagger\\boldsymbol\\alpha\\psi)\\). Write \\(F\_{\\mu\\nu}=\\partial\_\\mu A\_\\nu-\\partial\_\\nu A\_\\mu\\) in the Lorenz gauge \\(\\partial\_\\mu A^\\mu=0\\):

\\\[\\Box A^\\nu=-j^\\nu,\\qquad\\Box=-\\partial_t^2+\\nabla^2,\\qquad A^\\nu(x)=\\int d^4x'\\,\\mathcal G(x-x')\\,j^\\nu(x'),\\qquad\\Box\\,\\mathcal G=-\\delta^4 .\\tag{3}\\\]

3. put this \\(A\\) into V.2’s turn and (2) is reproduced letter for letter: the rule is Maxwell solved for walker 1’s current and read by walker 2’s segments, not a rule on top of (F). The gauge is free: (F)’s source is conserved, \\(\\partial\_\\nu j^\\nu=0\\), so \\(\\partial\_\\nu A^\\nu=\\int\\mathcal G\\,\\partial'\_\\nu j^\\nu=0\\) — the \\(A\\) the rule builds is always in Lorenz gauge. And \\(\\mathcal G\\) is the same object twice: Part II’s step has generator \\(D\\Box\\) and summed over its clock is \\(\\Box^{-1}\\), \\(\\frac1{4\\pi^2s^2}\\) in the signature where the walk is real ((E)); Maxwell’s inverse is \\(\\Box^{-1}\\) in real time. Static check for both, \\(\\int dt\\,\\mathcal G=\\frac1{4\\pi r}\\), Coulomb.

*What Maxwell does not fix.* \\(\\Box\\mathcal G=-\\delta^4\\) has many solutions,

\\\[\\mathcal G\_{\\rm ret}=\\frac{\\delta(t-r)}{4\\pi r},\\qquad\\mathcal G\_{\\rm adv}=\\frac{\\delta(t+r)}{4\\pi r},\\qquad\\mathcal G\_{\\rm sym}=\\tfrac12(\\mathcal G\_{\\rm ret}+\\mathcal G\_{\\rm adv})=\\frac{\\delta(t^2-r^2)}{4\\pi},\\qquad\\mathcal G_F=\\mathcal G\_{\\rm sym}+\\tfrac i2\\mathcal G^{(1)},\\tag{4}\\\]

4. differing by solutions of \\(\\Box f=0\\): free radiation present besides what the walkers make. Choosing \\(\\mathcal G\\) is a boundary condition, a statement about the state. Retarded: none comes in. Time-symmetric (Fokker, Wheeler–Feynman): what passes between two walkers, half out and half back on the light cone, which the absorber condition makes act as retarded. Feynman: the vacuum’s, adding the Hadamard function \\(\\mathcal G^{(1)}\\), the field’s own fluctuations, which a direct interaction between walkers does not have. On the lattice the three agree for steady sources ((G) uses the static sum \\(\\sum_nP^n=(1-P)^{-1}\\), which cannot tell them apart). This is, for the carrier, the question VI.4 answers for the walker — which state: *the walk and Maxwell fix how a current makes a field; what neither fixes is how much free field there is with no source.* For the walker the long-clock count selects the state (VI.4 (C)); for the carrier it is left here as a choice. With the Fokker choice there is none, and the fluctuations that a vacuum calculation needs have to come from the walkers themselves — which the companion document’s plates section (G) shows they do, in the near field: the plates’ own upward lines, carried by \\(\\mathcal G\_{\\rm sym}\\)’s static end, give Lifshitz’s force with no \\(\\mathcal G^{(1)}\\).

*What went in.* (1)–(3): no new input — V.2’s rule, (F)’s Maxwell and III.16’s kernel. (4): the choice of inverse is the state of the carrier — a choice, named.

**(G) The force, run as one process.** (A)–(F′) were checked piece by piece. Here the chain is run on a \\(3{+}1\\) lattice with nothing solved by hand: a walk gives the potential, a walk is pushed; then walker 1 is set moving and the magnetic sign is read.

*Coulomb.* Walker 1 sits at the centre of a \\(96^3\\) box; its density is the source. The carrier is Part II’s walker at \\(m=0\\) ((E)): a simple random walk from the source, absorbed at the box faces, whose total occupation \\(G(\\mathbf x)=\\sum_nP^n\\delta_0\\) — the kernel summed over its clock, III.16, i.e. \\((6-A)G=6\\delta_0\\) with \\(A\\) the neighbour matrix — is the potential, \\(\\Phi=eG\\). The link angle at each site is V.2’s time-link angle, \\(\\theta_t=-e\\Phi\\,dt\\). Walker 2 is IV.10’s three-coin automaton with the mass turn, a packet at rest at distance \\(d\\), positive-energy projected with the automaton’s own eigenvectors, turned by \\(\\theta_t\\) each tick; its momentum is read from its own lattice modes. (`force_from_two_walks.py`, 5 checks.) The carrier’s occupation, computed exactly and as \\(2\\times10^5\\) Monte Carlo walks, agree to \\(0.3\\%\\) shell by shell. Coulomb’s law from the walk: \\(\[G(r)-G(16)\]/\[\\frac6{4\\pi}(\\frac1r-\\frac1{16})\]=1.026,\\,1.012,\\,1.008,\\,1.006\\) at \\(r=4,6,8,12\\) (the box shifts \\(G\\) by a nearly uniform amount, which the differences remove); the gradient follows \\(\\frac6{4\\pi}/d^2\\) to \\(1\\)–\\(3\\%\\) at \\(d=8\\)–\\(16\\). The push on walker 2 is \\(-e\\langle\\partial_x\\Phi\\rangle\\), tick by tick, to \\(0.8\\%\\) — (D)’s statement with the turn rate made by another walk — away from walker 1 for like charges, toward it with walker 2’s charge reversed. And the walker is a cloud: \\(Fd^2/(e^2\\cdot\\frac6{4\\pi})=0.516,\\,0.865,\\,0.986\\) against the cloud’s weight within \\(d\\) of its centre, \\(0.511,\\,0.859,\\,0.979\\), at \\(d=8,12,16\\): a point source’s field averaged over a round cloud is the field at the centre on the part inside \\(d\\) and zero on the shells that enclose the source — the mean-value property of a harmonic field, which the walk’s occupation has.

*Magnetism, and the sign of \\(\\eta\\).* Now the rule (F′)(1) with walker 1 moving: steady line currents along \\(z\\) (plane waves \\(e^{ik_zz}\\) in \\(z\\), Gaussian in \\(x,y\\), \\(72^2\\times32\\) sites, \\(m=0.5\\)), so the static \\(\\mathcal G\\) suffices and the sources are walker 1’s time-segment density \\(\\rho_1\\) and its \\(z\\)-moves \\(j_1=\\psi_1^\\dagger\\alpha_z\\psi_1\\), with \\(A_t=\\kappa q_1\\eta\_{tt}(G\*\\rho_1)\\), \\(A_z=\\kappa q_1(G\*j_1)\\) — the one-form index of (E). Walker 2’s time turn per tick is \\(e^{iq_2A_t}\\) and each of its \\(z\\)-hops carries \\(e^{\\pm iq_2A_z}\\). (`vertex_rule_magnetism.py`, 4 checks.) Walker 1’s \\(z\\)-segments average to its current: \\(j_1/\\rho_1=0.7727\\) against the automaton’s group velocity \\(d\\omega/dk_z\\) averaged over the packet’s transverse modes, \\(0.7736\\). Walker 2’s push per tick is \\(q_2\\langle-\\partial_x\\Phi+\\alpha_z\\partial_xA_z\\rangle\\), the Lorentz force of (D), to \\(1.3\\%\\) on every run. With \\(\\eta\\): like static charges repel (\\(F_0>0\\)); parallel currents give \\(F/F_0=0.401\\) against \\(1-v_1v_2=0.403\\), antiparallel \\(1.617\\) against \\(1+v_1^2=1.597\\) — parallel currents attract magnetically. With a Euclidean \\(\\delta\\) in place of \\(\\eta\\) (time–time entry \\(+1\\), and \\(\\kappa\\) flipped so that like static charges still repel), parallel currents give \\(1.617\\) against \\(1+v_1v_2\\): the wrong magnetism. The relative sign of electricity and magnetism is \\(\\eta\_{00}\\)’s, II.1’s sign once more.

*What went in ((G)).* V.2’s step, (F′)’s rule with the sign of the time link that (E) ties to \\(\\eta\_{00}\\), and the carrier as a walk ((E), a reading). The carrier is static here — the walk summed over its clock — so this is Coulomb and the magnetostatic force, not radiation; (F)’s own dynamics of the link angle and the light-cone kernel of (F′)(4) are not exercised. In three space dimensions that kernel is supported on the light cone (Huygens): a walk at \\(c\\) that never turns — the companion document’s “carrier that never flips” — with no sense and a one-form index.

*What went in (the section).* (A)–(D): no new input — III.1’s rule written on V.2’s step, and equalities. (E): the carrier as Part II’s walker at \\(m=0\\) is the Casimir document’s sentence read with VI.1 (F); the kernel and Coulomb’s law are then III.16, derived. (F): one input — the link angle is a noise of the event with a local quadratic weight, of the same standing as II.1’s variance matrix; given it, \\(F^2\\) is forced and the source is the walker’s current. (F′): no new input, and one choice named, the carrier’s state. (G): runs. Not derived: that input, and the scale \\(e\\). Part V’s two open lines, “why the field produces that turn angle” and “the field’s own dynamics”, are not closed but relocated: (B) and (D) show that a link angle no clock-zero removes *is* a field and acts as the Lorentz force; why there is such an angle at all is (F)’s one input. Two open lines became one input, named.

(`scripts/pair_angle_is_a_connection.py` (A)–(C); `scripts/lorentz_force_turn_angle.py` (D); `scripts/force_is_a_turn_rate_gradient.py` (D), the three forces; `scripts/carrier_and_field_action.py` (E)–(F); `scripts/carrier_sign_and_current.py` the sign of the carrier and the conservation of the source; `scripts/interaction_and_plates_checks.py` (F′); `scripts/force_from_two_walks.py`, `scripts/vertex_rule_magnetism.py` (G).)

## VI.3 Bell: where the entanglement is, and which of Bell’s three the walk does not have

*The picture first.* Bell’s inequality needs two particles whose outcomes are correlated beyond what any assignment of outcomes at the source allows. The question for this document is not whether the walk reproduces the correlations — its column obeys Dirac, so it must — but where the entangled event comes from with nothing written down by hand, and which of Bell’s assumptions the walk fails to have. The first answer: the entanglement is the exchange term VI.1 was forced to keep — the pairing \\(E\\) in \\((D-E)^2\\) — when the two walkers share one packet and differ only in their axes; the free step carries it to the two detectors. The second: the walk has local steps and setting-independent coins; what it does not have is an outcome fixed at the source, because its weights belong to closed events whose return legs start at the detectors. One more thing the run showed: in \\(1{+}1\\) there is no violation, because the only two-valued thing a walker carries there is its direction of motion, which the geometry has already read by the time it reaches a detector. Bell needs the axis of IV.3; it is a \\(3{+}1\\) statement.

*Symbols.* \\(\\chi(\\hat{\\mathbf n})\\) the half-angle column whose axis is \\(\\hat{\\mathbf n}\\) (III.8a, IV.3), \\(P(\\hat{\\mathbf n})=\\chi(\\hat{\\mathbf n})\\chi(\\hat{\\mathbf n})^{\\mathsf T}=\\tfrac12(1+M(\\hat{\\mathbf n}))\\) the sorting of a column by that axis (IV.2a). \\(\\varphi_1,\\varphi_2\\) the two walkers’ columns; the two-walker column \\(C(x,y)=\\varphi_1(x)\\varphi_2(y)^{\\mathsf T}-\\varphi_2(x)\\varphi_1(y)^{\\mathsf T}\\) of VI.1 (D), a matrix in the two walkers’ axis indices at each pair of sites, carried by the free step as \\(C\\to UCU^{\\mathsf T}\\) and read at the detectors as \\(P_aC(x,y)P_b^{\\mathsf T}\\); \\(D\\), \\(E\\) its two terms, VI.1’s two pairings of emitters with absorbers; \\(L\\), \\(R\\) the two detector regions; \\(a\\), \\(b\\) the angles of the two detector axes in the plane perpendicular to the line joining them, measured from \\(\\hat{\\mathbf x}\\) (the axis angle here, not the flip rate; and the correlation \\(\\mathcal E(a,b)\\) of (6), not the pairing). The weights \\(W\\) below are the unnormalized weights of III.15a; the probabilities of the selected trials — one walker at each detector — are \\(W/\\sum W\\), which is where (6) divides. One basis note: in IV.3 the column \\((1,0)^{\\mathsf T}\\) was the axis \\(+\\hat{\\mathbf z}\\); here the real columns of III.8a are used with \\((1,0)^{\\mathsf T}=\\chi(\\hat{\\mathbf x})\\), the two related by the fixed rotation that takes \\(\\sigma_x\\to\\sigma_z\\), \\(\\sigma_y\\to\\sigma_x\\), under which \\(M(\\hat{\\mathbf n}\_a)=\\cos a\\,\\sigma_z+\\sin a\\,\\sigma_x\\) is real. Nothing else is new.

**(A) The ingredients, and how the state is prepared.** Five ingredients, all on the page. (1) Two identical walkers in one packet — one place, one motion, two axes — which VI.1 (E) allows: two processes at one point in different directions. (2) One process per state, hence the odd combination \\(D-E\\), VI.1 (D). (3) The free step on each walker, IV.10’s three coins and the mass turn. (4) A setting: the detector sorts the column it absorbs by an axis, \\(P(\\pm\\hat{\\mathbf n}\_a)\\) — the two streams of the half-turn \\(M(\\hat{\\mathbf n}\_a)\\) (IV.4d (A): a stream is a column a direction matrix gives back), one absorbed as “\\(+\\)”, the other as “\\(-\\)” (V.5 (D), by component). (5) III.15a’s weight on the closed two-walker event. No interaction, no state written in.

*The preparation is (1) and (2) together, and it needs no vertex.* Write the two walkers’ columns as the packet times their axis columns, \\(\\varphi_i=f\\chi_i\\). For any two 2-vectors,

\\\[\\chi_1\\chi_2^{\\mathsf T}-\\chi_2\\chi_1^{\\mathsf T}=\\big(\\chi\_{1,1}\\chi\_{2,2}-\\chi\_{1,2}\\chi\_{2,1}\\big)\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}=\\big(\\chi_1^{\\mathsf T}T\\chi_2\\big)\\,T ,\\tag{0}\\\]

0. so the odd combination of two identical walkers sharing one packet has, for its axis part, the turn matrix \\(T\\) of III.8 — now a quarter-turn between the two walkers’ axes — times a number, *whatever two axes they were launched with*, so long as they differ. \\(T\\) is the column of two axes with no net axis: the pair has no direction, not because angular momentum was conserved at a vertex, but because exclusion leaves one antisymmetric column and this is it. On the lattice, walkers launched with axes \\(180^\\circ\\), \\(60^\\circ\\) and \\(135^\\circ\\) apart give the same \\(S\\) to four digits (\\(2.7853\\), on the \\(1/(mw)=0.08\\) row of (D)), only the weight of the pair changing, as \\(\\sin^2\\tfrac{\\theta\_{12}}2\\) — the number in (0) (`bell_two_walkers_3plus1.py` with the second axis as argument). So nothing has to be prepared beyond putting two identical walkers in one packet; the free step then spreads the packet to both detectors, and the post-selection of one walker at each is what the detectors do. Below, the axes are taken as \\(\\pm\\hat{\\mathbf x}\\) for the arithmetic; (0) says the choice does not matter.

*The state from an interaction: two independent walkers made to share a packet.* The experiment can also start from two walkers that have nothing to do with each other, coming in from opposite sides, and let them meet. Two facts from Part V and VI.2 do the work. The field of VI.2 couples to the charge, not to the axis — V.4: the coin is untouched; V.3: the only axis coupling is the Pauli term, \\(O(\\hbar/m)\\) — so the two-walker amplitude is a spatial part \\(F(x_1,x_2)\\) times the two axis columns, and VI.1’s exchange makes it \\(F(1,2)\\chi_1\\chi_2^{\\mathsf T}-F(2,1)\\chi_2\\chi_1^{\\mathsf T}\\). At a pair of detectors \\((L,R)\\) the axis part is \\(F(L,R)\\chi_1\\chi_2^{\\mathsf T}-F(R,L)\\chi_2\\chi_1^{\\mathsf T}\\): “walker 1 bounced back to \\(L\\)” against “walker 1 went through to \\(R\\)” — (0)’s two orderings, now carrying the *bounce* and *pass* amplitudes \\(r\\) and \\(t\\). Without an interaction the walkers pass through each other, \\(F(L,R)=0\\), and there is nothing at the detectors to entangle. What the detectors read from this pair is computed in closed form in (B′) below, after the shared-packet case; the lattice (two \\(1{+}1\\) walkers, the rotation coin, the other walker’s potential entering as V.2’s phase per step, a Gaussian repulsion of range \\(6\\) and height \\(V_0\\); coin angle \\(0.5\\) per tick, packets of width \\(14\\) with wavenumbers \\(\\pm0.35\\), \\(300\\) ticks; the axis carried as the spectator it is) gives: free, \\(w\_{\\rm b}=0.000\\), \\(S=2.000\\); \\(V_0=0.2\\), bounce \\(0.39\\), pass \\(0.61\\), \\(S=2.49\\), matching (B′) at every \\(V_0\\). So the two preparations are one: a bound pair let go, or two strangers made to bounce off each other half the time — either way, exclusion with one shared packet, and nothing written in.

**(B) The weight of the two-walker event, and the correlation.** The slow walkers of IV.6: each column is the packet times its axis column, \\(\\varphi_1=f\\,\\chi(\\hat{\\mathbf x})\\), \\(\\varphi_2=f\\,\\chi(-\\hat{\\mathbf x})\\), and the free step keeps this form — the coins move \\(f\\), the mass turn acts on the senses, neither turns the axis of a packet at rest (IV.10: the precession is \\(\\mathbf p\\times\\hat{\\mathbf s}\\), zero on average), and both walkers get the same \\(f\\). In the plane of the detector axes the columns are III.8a’s real ones: \\(\\chi(\\hat{\\mathbf x})=(1,0)^{\\mathsf T}\\), \\(\\chi(-\\hat{\\mathbf x})=(0,1)^{\\mathsf T}\\), \\(\\chi(\\hat{\\mathbf n}\_a)=(\\cos\\tfrac a2,\\sin\\tfrac a2)^{\\mathsf T}\\), half the angle of the axis. Absorber endpoints: \\((x\\in L,\\ \\pm_a)\\) and \\((y\\in R,\\ \\pm_b)\\).

\\\[D=\\big\[P(\\hat{\\mathbf n}\_a)\\varphi_1(x)\\big\]\\big\[P(\\hat{\\mathbf n}\_b)\\varphi_2(y)\\big\]^{\\mathsf T},\\qquad E=\\big\[P(\\hat{\\mathbf n}\_a)\\varphi_2(x)\\big\]\\big\[P(\\hat{\\mathbf n}\_b)\\varphi_1(y)\\big\]^{\\mathsf T}\\tag{1}\\\]

1. VI.1’s two pairings, with the detectors’ sorting on the endpoints: walker 1 absorbed at \\((x,+\_a)\\) and walker 2 at \\((y,+\_b)\\), or the other way round. Each is an outer product of two sorted columns.

\\\[W(+\_a,+\_b)=\\sum\_{x\\in L}\\sum\_{y\\in R}\\big\\|D-E\\big\\|^2\\tag{2}\\\]

2. VI.1 (D): the odd combination, squared, III.15a’s rule on a two-walker event; the other ordering is the same configuration.

\\\[\\|D-E\\|^2=|P_a\\varphi_1(x)|^2|P_b\\varphi_2(y)|^2+|P_a\\varphi_2(x)|^2|P_b\\varphi_1(y)|^2-2\\,\\big(P_a\\varphi_1(x)\\big)\\!\\cdot\\!\\big(P_a\\varphi_2(x)\\big)\\,\\big(P_b\\varphi_2(y)\\big)\\!\\cdot\\!\\big(P_b\\varphi_1(y)\\big)\\tag{3}\\\]

3. for outer products, \\(\\|uv^{\\mathsf T}-u'v'^{\\mathsf T}\\|^2=|u|^2|v|^2+|u'|^2|v'|^2-2(u\\cdot u')(v\\cdot v')\\). The first two terms are the two pairings each with itself — what two *distinguishable* walkers give. The third is \\(D\\cdot E\\): the two orderings interfering, III.15a (E)’s pairs that go out one way and come back the other, now across the two detectors. It is the only term that ties the two readings beyond a product, and it exists because VI.1 (B) made the two orderings one configuration.

\\\[\\varphi_i=f\\chi_i:\\qquad|P_a\\varphi_1(x)|^2=f(x)^2\\cos^2\\tfrac a2,\\quad|P_a\\varphi_2(x)|^2=f(x)^2\\sin^2\\tfrac a2,\\quad\\big(P_a\\varphi_1(x)\\big)\\!\\cdot\\!\\big(P_a\\varphi_2(x)\\big)=f(x)^2\\,\\chi(\\hat{\\mathbf x})^{\\mathsf T}P(\\hat{\\mathbf n}\_a)\\chi(-\\hat{\\mathbf x})=f(x)^2\\cos\\tfrac a2\\sin\\tfrac a2\\tag{4}\\\]

4. \\(P(\\hat{\\mathbf n}\_a)\\chi(\\hat{\\mathbf x})=\\chi(\\hat{\\mathbf n}\_a)\\,(\\chi(\\hat{\\mathbf n}\_a)^{\\mathsf T}\\chi(\\hat{\\mathbf x}))=\\cos\\tfrac a2\\,\\chi(\\hat{\\mathbf n}\_a)\\); the squares are IV.7’s \\(\\tfrac12(1\\pm\\hat{\\mathbf n}\_a\\cdot\\hat{\\mathbf x})=\\tfrac12(1\\pm\\cos a)\\) — how much of the \\(\\hat{\\mathbf x}\\) column is in the \\(\\hat{\\mathbf n}\_a\\) stream — and the cross term is the product of the two overlaps, \\(\\cos\\tfrac a2\\cdot\\sin\\tfrac a2=\\tfrac12\\sin a\\). Same at \\(R\\) with \\(b\\) and the roles of \\(\\chi(\\pm\\hat{\\mathbf x})\\) exchanged.

\\\[W(+\_a,+\_b)=\\rho_L\\rho_R\\Big\[\\cos^2\\tfrac a2\\sin^2\\tfrac b2+\\sin^2\\tfrac a2\\cos^2\\tfrac b2-2\\cdot\\tfrac12\\sin a\\cdot\\tfrac12\\sin b\\Big\]=\\tfrac12\\rho_L\\rho_R\\big\[1-\\cos a\\cos b-\\sin a\\sin b\\big\]=\\tfrac12\\rho_L\\rho_R\\big\[1-\\cos(a-b)\\big\]\\tag{5}\\\]

5. 4. into (3) into (2), \\(\\rho_L=\\sum_Lf^2\\), \\(\\rho_R=\\sum_Rf^2\\); \\(\\cos^2\\tfrac a2\\sin^2\\tfrac b2+\\sin^2\\tfrac a2\\cos^2\\tfrac b2=\\tfrac12(1-\\cos a\\cos b)\\). The \\(-\\sin a\\sin b\\) is the exchange term. Likewise \\(W(-\_a,-\_b)=W(+\_a,+\_b)\\) and \\(W(+\_a,-\_b)=W(-\_a,+\_b)=\\tfrac12\\rho_L\\rho_R\[1+\\cos(a-b)\]\\); the four add to \\(2\\rho_L\\rho_R\\).

\\\[\\mathcal E(a,b)=\\frac{W(++)+W(--)-W(+-)-W(-+)}{W(++)+W(--)+W(+-)+W(-+)}=-\\cos(a-b)\\tag{6}\\\]

6. the correlation of the two readings. It depends only on the angle between the two detector axes, as the pair’s \\(J=0\\) demands, and it is \\(-1\\) when the axes agree: the walker read \\(+\\) at one detector is read \\(-\\) at the other, which is (1): the two walkers’ axes are opposite.

\\\[S=\\mathcal E(0,\\tfrac\\pi4)-\\mathcal E(0,\\tfrac{3\\pi}4)+\\mathcal E(\\tfrac\\pi2,\\tfrac\\pi4)+\\mathcal E(\\tfrac\\pi2,\\tfrac{3\\pi}4)=-\\tfrac{\\sqrt2}2-\\tfrac{\\sqrt2}2-\\tfrac{\\sqrt2}2-\\tfrac{\\sqrt2}2=-2\\sqrt2\\tag{7}\\\]

7. the CHSH sum at the standard angles; \\(|S|=2\\sqrt2>2\\).

\\\[\\text{distinguishable walkers (first two terms of (3) only):}\\qquad W(+\_a,+\_b)=\\tfrac12\\rho_L\\rho_R\\big\[1-\\cos a\\cos b\\big\],\\qquad \\mathcal E=-\\cos a\\cos b,\\qquad|S|\\le2\\tag{8}\\\]

8. without the exchange term the correlation factorizes, \\(\\mathcal E=-\\cos a\\cos b\\), and \\(|S|=|{-\\cos a}(\\cos b-\\cos b')-\\cos a'(\\cos b+\\cos b')|\\le|\\cos b-\\cos b'|+|\\cos b+\\cos b'|\\le2\\). The whole violation is \\(-\\sin a\\sin b\\) in (5): the term \\(E\\) of VI.1. (The even combination, \\(D+E\\), gives \\(+\\sin a\\sin b\\) and the same \\(|S|\\): the exchange term is the entanglement; its sign is exclusion, not Bell.) (`bell_half_angle_calculation.py`: (4)–(8) symbolically.)

**(B′) The scattered pair in closed form.** (A)’s second preparation leaves at the detectors the pair of columns \\(M=r\\,\\chi_1\\chi_2^{\\mathsf T}-t\\,\\chi_2\\chi_1^{\\mathsf T}\\) — \\(\\chi_1,\\chi_2\\) the launched axes, \\(\\theta\_{12}\\) the angle between them, \\(r\\) the bounce and \\(t=|t|e^{i\\delta}\\) the pass, \\(\\delta\\) the phase the barrier puts between them — and (2) reads it as before, \\(W(+\_a,+\_b)=\\|P_aMP_b^{\\mathsf T}\\|^2\\). Three lines of algebra on \\(M\\), then the correlation.

\\\[M=r\\,\\chi_1\\chi_2^{\\mathsf T}-t\\,\\chi_2\\chi_1^{\\mathsf T}=(r-t)\\,\\tfrac12\\big(\\chi_1\\chi_2^{\\mathsf T}+\\chi_2\\chi_1^{\\mathsf T}\\big)+(r+t)\\,\\tfrac12\\big(\\chi_1^{\\mathsf T}T\\chi_2\\big)\\,T\\tag{9}\\\]

9. 0. on the odd half: the pair is its symmetric part weighted by \\(r-t\\) and the turn \\(T\\) weighted by \\(r+t\\). It is \\(T\\) times a number — the shared-packet pair of (B) — only when \\(r=t\\) exactly; but that is not what the detectors are sensitive to.

\\\[M^{\\mathsf T}TM=(\\det M)\\,T,\\qquad\\det M=r\\,t\\,\\big(\\chi_1^{\\mathsf T}T\\chi_2\\big)^2=r\\,|t|\\,e^{i\\delta}\\sin^2\\tfrac{\\theta\_{12}}2\\tag{10}\\\]

10. the pair with \\(T\\) between its two columns, read against itself, is \\(T\\) again, times \\(\\det M\\); each outer product has no determinant of its own, and for two of them \\(\\det(uv^{\\mathsf T}+u'v'^{\\mathsf T})=(u\\wedge u')(v\\wedge v')\\), which here is \\((r\\chi_1\\wedge(-t)\\chi_2)(\\chi_2\\wedge\\chi_1)=rt\\,(\\chi_1^{\\mathsf T}T\\chi_2)^2\\). A turn of either detector’s axis is a unitary on one side of \\(M\\), which leaves \\(|\\det M|\\) and \\(\\|M\\|^2\\) alone; the number the two detectors can read from the pair, whatever their axes, is therefore the ratio

\\\[\\|M\\|^2=|r|^2+|t|^2-2|rt|\\cos\\delta\\,\\cos^2\\tfrac{\\theta\_{12}}2,\\qquad C\\equiv\\frac{2|\\det M|}{\\|M\\|^2}=\\frac{2|rt|\\sin^2\\tfrac{\\theta\_{12}}2}{|r|^2+|t|^2-2|rt|\\cos\\delta\\cos^2\\tfrac{\\theta\_{12}}2}\\tag{11}\\\]

11. \\(\\|M\\|^2\\) is III.15a’s weight of the pair summed over all four outcomes; the cross term is twice the real part of \\(r\\bar t\\,(\\chi_1^{\\mathsf T}\\chi_2)^2\\). \\(C\\) is \\(0\\) for equal axes (\\(\\theta\_{12}=0\\): \\(M=(r-t)\\chi\\chi^{\\mathsf T}\\), one column twice, nothing to read) and largest for opposite ones; the dependence is \\(\\sin^2\\tfrac{\\theta\_{12}}2\\), the square of (0)’s number, not \\(\\sin\\theta\_{12}\\).

\\\[\\theta\_{12}=\\pi:\\quad C=\\frac{2|rt|}{|r|^2+|t|^2};\\qquad\\delta=\\pm\\tfrac\\pi2:\\quad C=\\frac{2|rt|\\sin^2\\tfrac{\\theta\_{12}}2}{|r|^2+|t|^2}\\tag{12}\\\]

12. opposite launched axes: \\(C=1\\) whenever bounce and pass are equal in size, *for any* \\(\\delta\\); the symmetric barrier’s \\(\\delta=\\pm\\pi/2\\) (the lattice’s \\(\\arg\\sum F(L,R)F(R,L)^\*=-\\pi/2\\) at every \\(V_0\\)) removes the \\(\\cos\\delta\\) term. On the lattice \\(r\\), \\(t\\) are packets, \\(|r|^2\\to w\_{\\rm b}=\\sum|F(L,R)|^2\\) (bounce), \\(|t|^2\\to w\_{\\rm p}\\) (pass), \\(r\\bar t\\to X=\\sum F(L,R)\\overline{F(R,L)}\\), so \\(C=2|X|/(w\_{\\rm b}+w\_{\\rm p})\\) — weights, written \\(w\\) so as not to collide with the turn \\(T\\) and the region \\(R\\).

\\\[\\chi(\\theta,\\phi)=\\big(\\cos\\tfrac\\theta2,\\ e^{i\\phi}\\sin\\tfrac\\theta2\\big)^{\\mathsf T},\\quad\\theta\_{12}=\\pi:\\qquad \\mathcal E(a,b)=-\\cos\\theta_a\\cos\\theta_b-C\\,\\sin\\theta_a\\sin\\theta_b\\cos(\\phi_a-\\phi_b-\\delta)\\tag{13}\\\]

13. the correlation (6) for detector axes anywhere on the sphere — IV.3’s columns, \\(\\theta\\) from the launched axis \\(\\hat{\\mathbf x}\\), \\(\\phi\\) the turn about it, \\(P(\\hat{\\mathbf n})=\\chi\\chi^\\dagger\\) — by (2)–(3) with these columns, the four weights as in (5). With \\(C=1\\), \\(\\delta=0\\) it is \\(-\\hat{\\mathbf n}\_a\\cdot\\hat{\\mathbf n}\_b\\), (6) again. The bounce–pass phase sits in the cross term as \\(\\cos(\\phi_a-\\phi_b-\\delta)\\): for the symmetric barrier the two detectors must sort in planes through \\(\\hat{\\mathbf x}\\) a quarter-turn apart; sorting in one plane (\\(\\phi_a=\\phi_b\\), the real columns of (B)) reads the real part of an imaginary cross term, \\(E=-\\cos\\theta_a\\cos\\theta_b\\), and \\(|S|\\le2\\) by (8) — the lattice’s \\(S=2.000\\) in the plane. The third half-turn of IV.3 is what reads the pair.

\\\[\\phi_a-\\phi_b=\\delta,\\ \\ \\theta_a=0,\\ \\theta_a'=\\tfrac\\pi2,\\ \\tan\\theta_b=C,\\ \\theta_b'=\\pi-\\theta_b:\\qquad|S|=\\frac{2}{\\sqrt{1+C^2}}+\\frac{2C^2}{\\sqrt{1+C^2}}=2\\sqrt{1+C^2}\\tag{14}\\\]

14. 13. at these angles: \\(\\mathcal E(a,b)-\\mathcal E(a,b')=-2\\cos\\theta_b\\), \\(\\mathcal E(a',b)+\\mathcal E(a',b')=-2C\\sin\\theta_b\\), and \\(\\cos\\theta_b=1/\\sqrt{1+C^2}\\), \\(\\sin\\theta_b=C/\\sqrt{1+C^2}\\). \\(|S|>2\\) for every \\(C>0\\): any bounce at all entangles, and \\(2\\sqrt2\\) needs \\(|r|=|t|\\). That (14) is the largest value over all axes is Horodecki’s theorem; the lattice optimizer finds it: \\(S=2.194,\\ 2.489,\\ 2.150\\) at \\(V_0=0.15,\\ 0.2,\\ 0.3\\), equal to \\(2\\sqrt{1+C^2}\\) with \\(C=2|X|/(w\_{\\rm b}+w\_{\\rm p})\\) to the optimizer’s tolerance. What keeps the lattice below \\(2\\sqrt2\\) is \\(|X|/\\sqrt{w\_{\\rm b}w\_{\\rm p}}=0.76\\): \\(r\\) and \\(t\\) carry different phases for the different \\(k\\) across the packet. *This is regular quantum mechanics’ answer*: for two identical fermions scattered by a spin-independent potential, the outgoing axis state is the two-qubit state \\(r\\chi_1\\chi_2-t\\chi_2\\chi_1\\), (11) is its concurrence (Wootters 1998) and (14) its maximal CHSH value (Horodecki, Horodecki, Horodecki 1995). Here the same numbers come from (0), (2) and III.15a’s square on two walkers’ columns, with the pair made by VI.2’s field. (`bell_scattering_closed_form.py` 9–14; `bell_from_scattering.py` 1–3.)

**(C) In \\(1{+}1\\) the sense cannot do it (this setup).** The same construction with the sense as the two-valued thing, the detectors on the two sides of the source and the settings turns of the sense pair (V.1’s first knob) gives \\(S=1.81\\) for the odd state, the even state and the product state alike (\\(1.807\\), \\(1.811\\), \\(1.810\\), the last being the product state’s closed form \\(2|r_L||r_R|\\); a point source, coin angle \\(0.2\\), \\(20\\) ticks, detectors beyond \\(8\\) sites). The reason: in \\(1{+}1\\) the sense *is* the direction of motion, so a walker that reached the left detector is the backward sense with weight \\(0.96\\) — nothing is left to turn, and (4)’s overlaps are fixed by the side. The single-site sense matrix is of rank two, so the exchange term is there, but the degree of freedom it entangles has been read by the geometry. Bell needs a two-valued thing the motion does not fix: the axis, which exists from \\(2{+}1\\) up and is free to turn from \\(3{+}1\\) up (in \\(2{+}1\\) it can only point out of the plane, IV.8 (C)). The statement is about an internal two-valued reading at fixed detectors; an entanglement of *arrival times* along the line, read with path delays and the pair angle of V.2 (Franson’s arrangement), is a different setup, not done here. (`bell_two_walkers_1plus1.py` 1, 2.)

**(D) The run: the walk’s own step, and the relativistic corrections to (6).** (B) assumed the shared packet and the fixed axes, which hold for slow walkers. The lattice checks both and finds the correction: one Gaussian packet \\(e^{-r^2/2w^2}\\) at rest, axes \\(\\pm\\hat{\\mathbf x}\\), odd; the exact Dirac evolution for the time \\(t=1.5\\,mw^2\\) (and, beside it, the three-coin step of IV.10, which converges to it); detectors \\(z\<-w\\) and \\(z>w\\), the edge placed midway between two grid planes so that the sums over a detector are second order in the grid spacing; (2) evaluated with the evolved columns and \\(S\\) maximized over the four angles. The numbers depend on that time and on where the detectors sit — at \\(1/(mw)=0.5\\) with the detectors at \\(d=0\\) the coefficient of \\(\\cos b\\) below is \\(0.56\\), \\(0.69\\), \\(0.76\\), \\(0.69\\), \\(0.69\\), \\(0.70\\) at \\(t=1,\\ 1.5,\\ 2,\\ 3,\\ 4,\\ 8\\ mw^2\\), not monotonic, against the table’s \\(0.61\\) with \\(d=w\\) at \\(1.5\\,mw^2\\) — so the setup is part of the table. (The run is written in the mass-diagonal basis \\(\\beta=\\mathrm{diag}(1,1,-1,-1)\\), \\(\\alpha_k=\\begin{pmatrix}0&\\sigma_k\\\\\\sigma_k&0\\end{pmatrix}\\), reached from IV.4’s chiral basis by \\(W=\\tfrac1{\\sqrt2}\\begin{pmatrix}1&1\\\\1&-1\\end{pmatrix}\\otimes1\\) and from IV.4b’s real counters by that after the chiral map of IV.5; the axis is sorted with \\(\\Sigma_k=\\mathrm{diag}(\\sigma_k,\\sigma_k)\\), the same \\(\\sigma_k\\) on both pairs, not with the direction matrices \\(\\alpha_k\\).)

| \\(1/(mw)\\), the packet’s wavenumber scale over the mass | odd | even | product | \\(\\mathcal E(0,b)\\) |
| --- | --- | --- | --- | --- |
| \\(0.5\\) | \\(1.713\\) | \\(1.713\\) | \\(1.211\\) | \\(-0.606\\cos b\\) |
| \\(0.25\\) | \\(2.486\\) | \\(2.486\\) | \\(1.758\\) | \\(-0.879\\cos b\\) |
| \\(0.08\\) | \\(2.785\\) | \\(2.785\\) | \\(1.970\\) | \\(-0.985\\cos b\\) |
| \\(0.04\\) | \\(2.814\\) | \\(2.814\\) | \\(1.990\\) | \\(-0.995\\cos b\\) |

The table is regular quantum mechanics — two free Dirac particles in the \\(J=0\\) state, evolved by the Dirac equation, the axis read with \\(\\tfrac12(1+\\hat{\\mathbf n}\\cdot\\boldsymbol\\Sigma)\\) on the full column — and the walk’s step is to be compared with it. In every row the odd state exceeds the product state by exactly \\(\\sqrt2\\), which is the exchange term; it exceeds \\(2\\) in the three slower rows and not in the fastest; even and odd agree, as (8) said; \\(\\mathcal E(0,b)\\) is \\(-\\kappa\\cos b\\) with \\(S=2\\sqrt2\\,\\kappa\\), to round-off. The walk’s step converges to this limit at first order in the step: for \\(1/(mw)=0.5\\) at \\(t=3.2/m\\) the two columns differ by \\(0.255\\), \\(0.127\\), \\(0.063\\) at \\(m\\varepsilon=0.4,\\ 0.2,\\ 0.1\\) and the step gives \\(S=1.455\\), \\(1.642\\), \\(1.695\\) against \\(1.718\\), \\(1.714\\), \\(1.714\\) on the same lattices. A step of \\(m\\varepsilon=0.4\\) is too coarse for the table’s own times (\\(60\\) ticks at \\(1/(mw)=0.25\\) give \\(1.18\\) against \\(2.49\\)), which is why the table is the exact evolution. The fall-off with \\(1/(mw)\\) is physics: the two-walker state stays \\(J=0\\), but as the walkers get faster the small components carry axis and motion coupled (IV.8 (C)), the shared-packet form of (B) is lost in part, and \\(E\\) drops below \\(-\\cos(a-b)\\). \\(S\\to2\\sqrt2\\) only as \\(1/(mw)\\to0\\). One caveat on the relativistic rows: they are quantum mechanics’ numbers *for this detector*, which sorts the full column by \\(\\boldsymbol\\Sigma\\) — IV.7a’s half-turn, the document’s reading of “the axis”; a detector reading the mean spin of Foldy and Wouthuysen instead would give a different fall-off, and the two agree only as \\(1/(mw)\\to0\\). The nonrelativistic row is the one every reading shares. No signalling: the marginal at \\(L\\) is \\((0.5,0.5)\\) for every setting at \\(R\\), to round-off. The table is the same on \\(64^3\\) and \\(96^3\\) points to \\(6\\times10^{-4}\\). (`bell_two_walkers_3plus1.py` 1, 1a, 1\\('\\), 2–4.)

**(E) Bell’s three, in the walk’s words.** Bell’s bound follows from writing the joint weight as

\\\[P(A,B\\,|\\,a,b)=\\sum\_\\lambda p(\\lambda)\\,P(A\\,|\\,a,\\lambda)\\,P(B\\,|\\,b,\\lambda),\\tag{15}\\\]

with \\(\\lambda\\) fixed at the source and independent of \\(a,b\\). Three things are assumed. *Locality of the step:* nothing at one detector changes the other walker’s step. The walk has it — every coin is drawn at its own event (VI.1 (B)), and the three-coin step is local. *Independence of the settings:* \\(p(\\lambda)\\) does not depend on \\(a,b\\). The walk has it — no coin and nothing at the source knows \\(\\hat{\\mathbf n}\_a\\) or \\(\\hat{\\mathbf n}\_b\\); the settings sort columns in the detectors only, and the marginals say so. *Outcome independence:* given what left the source, \\(\\lambda\\), the two readings are independent — \\(P(A,B|a,b,\\lambda)=P(A|a,\\lambda)P(B|b,\\lambda)\\). Note what this does not say: it does not say the outcomes are fixed by \\(\\lambda\\); \\(P(A|a,\\lambda)\\) may be any number in \\(\[0,1\]\\), and a local model with random outcomes obeys (15) and the bound just the same, so “the outcome was not decided at the source” is no way around Bell. What the walk lacks is the factorization itself: in III.8c (A)‘s words, the weight belongs to an outward route and a return route that agree at *both* detectors, and it is never source data times two local readings. Take \\(\\lambda\\) to be the pair that left the source, (0): at equal settings (5) gives \\(W(+\_a,+\_a)=0\\), while the product of the two marginals is \\(\\tfrac14\\). The reason is the event: it has two endpoints, and its weight (2) is built from the return legs, which start at the absorbers with the sorting on them. So (15) has nothing to factor: the product \\(P(A|a,\\lambda)P(B|b,\\lambda)\\) is replaced by the square of \\(D-E\\), whose cross term \\(-\\sin a\\sin b\\) pairs both walkers’ out- and return legs. (Which event closes is the only random thing, III.15a (D)3; the step itself is orthogonal.) This is the transactional way around Bell (Cramer), here with everything real and the two legs derived rather than posited; it is the price of inputs 9 and 10 of the one list.

What this section does not claim: that the violation is derived from the coin — it is derived from the exchange term and the shared packet, both already on the page — nor that the walk is nonlocal in any signalling sense: \\(\\sum_BW\\) in (2) sums \\(P(\\pm\\hat{\\mathbf n}\_b)\\) to \\(1\\) and the far setting drops out.

*What went in.* No new input. (A) is VI.1, IV.10, IV.2a and V.5, and (0) is two lines of algebra on the odd combination; (B) is III.15a’s rule on VI.1’s state with III.8a’s half-angle columns, and (B′) the same rule on the scattered pair with IV.3’s columns; (C)–(D) are runs. The reading of (E) — which assumption is absent, outcome independence — is the document’s definition of an event (input 9) read against (15).

(`scripts/bell_half_angle_calculation.py` (B); `scripts/bell_scattering_closed_form.py` (B′); `scripts/bell_from_scattering.py` the state from an interaction; `scripts/bell_two_walkers_1plus1.py` (C); `scripts/bell_two_walkers_3plus1.py` (D), with the automaton and the exact evolution.)

## VI.4 Many walkers on a cut: what creation and annihilation are, read off the walk

*The picture first.* Quantum field theory speaks of a space of states with any number of particles, and of operators that create and annihilate them. None of that is put under the walk here. A state on a cut is a list of which columns hold a walker; its amplitude is what exclusion makes of the sum over zigzags; “create” and “annihilate” are the names of two things one can do to that list. The vacuum is the one configuration nothing happens in, and from it Feynman’s rule, the antiparticle as a hole, pair creation as a lifted walker, and the vacuum’s energy as a sum over closed loops are computed, not stated. A missing walker turns out to be, in every amplitude, one walker of the opposite orientation and charge, so nothing has to be absent. The section ends with what is taken — that the walker is read by all its detectors, and that the vacuum is the configuration the long-clock count selects — and with what the ribbon’s own backward crossings have to do with the loops a vacuum needs.

**(A) A walker is a column; two walkers are a two-index column; the step makes a determinant.** III.1’s walker read in \\(t\\) (III.5a (H)) has at each site \\(x\\) two counters, \\(Q\_+(x)\\) for the forward sense (a right-mover) and \\(Q\_-(x)\\) for the backward sense (a left-mover, III.3’s antiparticle crossing). On three sites the column is

\\\[q=\\begin{pmatrix}Q\_+(1)\\\\Q\_-(1)\\\\Q\_+(2)\\\\Q\_-(2)\\\\Q\_+(3)\\\\Q\_-(3)\\end{pmatrix},\\\]

and one tick (III.15a (A): shift, then the coin, \\(c=\\cos a\\,dt\\), \\(s=\\sin a\\,dt\\)) is \\(Q'\_+(y)=c\\,Q\_+(y-1)+s\\,Q\_-(y+1)\\), \\(Q'\_-(y)=-s\\,Q\_+(y-1)+c\\,Q\_-(y+1)\\):

\\\[U=\\begin{pmatrix}0&0&0&s&c&0\\\\0&0&0&c&-s&0\\\\c&0&0&0&0&s\\\\-s&0&0&0&0&c\\\\0&s&c&0&0&0\\\\0&c&-s&0&0&0\\end{pmatrix},\\qquad q'=Uq,\\qquad U^{\\mathsf T}U=1 .\\tag{1}\\\]

1. row 1: the forward counter at site 1 came from site 3 (one step right, wrapping), the backward one from site 2, mixed by the coin. \\(U^t\_{XA}\\), \\(X=(x,\\varepsilon)\\), is the sum over all zigzags from \\(A\\) to \\(X\\) in \\(t\\) ticks — III.15a’s reading.

*Which column.* From here on a walker’s column has complex entries: it is the walker read by all its detectors, \\(u+iv\\) (III.8c (B), III.9) — the original detectors \\(u\\) and the added ones \\(v\\), with the walker’s own turn \\(J\\) written as \\(i\\) (IV.5a). \\(U\\) is real and acts alike on the original and the added detectors. This fixes the count: the slots of (C)–(D) are complex directions. A walker in the original detectors alone is its own conjugate, has half as many independent slots and half the vacuum energy (at rest \\(-\\tfrac12\\hbar a\\) against \\(-\\hbar a\\); III.15e (B), `audit_2026_11_09.py` D).

Two walkers are a number for every pair of (site, sense), \\(C\_{XY}\\), \\(6\\times6\\) on three sites; \\(C\_{(1,+),(2,-)}\\) is “one walker at site 1 forward, the other at site 2 backward”. Exclusion in every state (VI.1 (D)) makes it antisymmetric, \\(C\_{XY}=-C\_{YX}\\), \\(C\_{XX}=0\\): two walkers at one site with one sense, never; with opposite senses, allowed (VI.1 (E)). The step is one \\(U\\) per walker:

\\\[C'\_{XY}=\\sum\_{K,M}U\_{XK}\\,U\_{YM}\\,C\_{KM}\\qquad(C'=UCU^{\\mathsf T}).\\tag{2}\\\]

Launch walker 1 at \\(A\\) and walker 2 at \\(B\\), \\(C=e_Ae_B^{\\mathsf T}-e_Be_A^{\\mathsf T}\\). After \\(t\\) ticks,

\\\[C\_{XY}=U^t\_{XA}\\,U^t\_{YB}-U^t\_{XB}\\,U^t\_{YA}:\\tag{3}\\\]

3. (the zigzags from \\(A\\) to \\(X\\))(the zigzags from \\(B\\) to \\(Y\\)) minus the same with the destinations swapped, VI.1’s \\(D-E\\). That two-term expression is what the word *determinant* names; with three walkers it is the six pairings with their signs. Nothing was added: \\(U\\) is the walk and the minus is exclusion. (`fock_space_from_the_walk.py` 1: the two- and three-walker columns evolved by (2), against (3) and its \\(3\\times3\\) form, at every configuration, to \\(10^{-15}\\).)

**(B) Adding and removing a walker.** The list of configuration amplitudes for \\(n\\) walkers is an antisymmetric \\(n\\)-index column \\(T\\). “Add a walker in column \\(v\\)” is a definite thing to do to it — put \\(v\\) in a new index and antisymmetrize, which is the Laplace expansion of the determinant along its new row:

\\\[(c^\\dagger_vT)\_{X_0X_1\\ldots X_n}=\\sum\_{j=0}^{n}(-1)^j\\,v\_{X_j}\\,T\_{X_0\\ldots\\widehat{X_j}\\ldots X_n};\\qquad(c_wT)\_{X_1\\ldots X_n}=\\sum\_{X_0}\\overline{w\_{X_0}}\\,T\_{X_0X_1\\ldots X_n}.\\tag{4}\\\]

4. the second is “remove the walker that overlaps \\(w\\)”. From (4) alone,

\\\[c_w\\,c^\\dagger_v+c^\\dagger_v\\,c_w=\\langle w,v\\rangle,\\qquad c^\\dagger_v\\,c^\\dagger_v=0 .\\tag{5}\\\]

5. the second equation is exclusion itself; the first is the polarization identity of VI.1 (D). These are the creation and annihilation operators of the textbook, and their algebra is the determinant’s antisymmetry — not an input (ibid. 2, to \\(10^{-15}\\)). The free step keeps the walker count: (2) has one \\(U\\) per index. Nothing in the walk takes a list of \\(n\\) walkers to a list of \\(n+1\\); that is III.3’s conserved charge.

**(C) The selected configuration, and Feynman’s rule from it.** \\(U\\) is real and orthogonal, so its standing patterns come in pairs: columns the tick turns one way, \\(U\\varphi=e^{-i\\varepsilon}\\varphi\\), and their conjugates turning the other, \\(e^{+i\\varepsilon}\\) — the walker whose clock turns one way, and the walker whose clock turns the other way (III.15e (C)). Half of all columns have \\(\\varepsilon\<0\\). A single walker in any one of them is a solution; nothing in the free walk in real time prefers a configuration (the sum over configurations is \\(\\prod(1+e^{-i\\varepsilon_jt})\\), phases). Two things do, and they agree.

*The count.* In the walk’s own clock the same sum has weights, not phases, \\(\\prod_j(1+e^{-\\varepsilon_j\\tau})\\), and at long clock it is carried by one configuration: every column with \\(\\varepsilon_j\<0\\) occupied, the others empty (VI.5 (C); frame by frame, III.15d). That the vacuum is what this limit selects is taken (VI.5).

*The lamp.* With V.2’s turn and a carrier that takes energy away, a walker in a positive column can step down into an empty negative one. The plate run measured it: a lone walker has \\(93\\%\\) of its line weight downward, at a rate \\(\\approx2mc^2\\), and there is no bottom. The only configuration with no line going down anywhere is the same one, every negative column occupied (`plate_of_walkers_sum_and_difference.py` 3, 6: the downward step of each walker is cancelled line by line by the exchange weights of the walkers below it).

Call it the sea: \\(N\\) walkers, one in each \\(\\varphi_j\\) with \\(\\varepsilon_j\<0\\), the antisymmetric \\(N\\)-index column \\(\\varphi_1\\wedge\\cdots\\wedge\\varphi_N\\).

Add one walker at \\(A\\) and ask for it at \\(X\\) after \\(t\\) ticks. The \\((N+1)\\)-walker determinant is a bordered one, and its Schur complement is one line:

\\\[\\frac{\\langle\\text{sea}+X\\,|\\,U^t\\,|\\,\\text{sea}+A\\rangle}{\\langle\\text{sea}\\,|\\,U^t\\,|\\,\\text{sea}\\rangle}=\\frac{\\det\\begin{pmatrix}\\Phi^\\dagger U^t\\Phi&\\Phi^\\dagger U^te_A\\\\e_X^\\dagger U^t\\Phi&e_X^\\dagger U^te_A\\end{pmatrix}}{\\det\\Phi^\\dagger U^t\\Phi}=(U^t)\_{XA}-\\sum_j e^{-i\\varepsilon_jt}\\varphi_j(X)\\overline{\\varphi_j(A)}=\\big(U^t\\,\\Pi\_+\\big)\_{XA}.\\tag{6}\\\]

6. \\(\\Phi\\) the matrix of sea columns; the subtracted term is the walker’s amplitude to have entered an occupied column, which the determinant removes; what is left is the positive half, forward. Remove a sea walker instead (a hole at \\(A\\), \\(c\_{e_A}\\) on the sea) and the same algebra gives

\\\[\\frac{\\langle\\text{hole at }X\\,|\\,U^t\\,|\\,\\text{hole at }A\\rangle}{\\langle\\text{sea}\\,|\\,U^t\\,|\\,\\text{sea}\\rangle}=\\overline{(U^t\\Pi\_-)\_{XA}}=\\big(\\Pi\_-\\,U^{-t}\\big)\_{AX}:\\tag{7}\\\]

7. the negative column, run from \\(X\\) back to \\(A\\) — and, because \\(U\\) is real, \\(\\overline{\\Pi\_-}=\\Pi\_+\\), the same number is \\((U^t\\Pi\_+)\_{XA}\\): **the amplitude of a hole from \\(A\\) to \\(X\\) is the amplitude of a walker in a positive column from \\(A\\) to \\(X\\).** Conjugation is counting the added detectors with the other sign, \\(v\\to-v\\) of V.2 (H): the clock turned the other way, \\(k\\to-k\\), \\(e\\to-e\\). With a field on, \\(W\_\\theta\\) in place of \\(U^t\\), the statement is exact: the hole’s amplitude in the field \\(\\theta\\) is the amplitude of a walker above the sea in the field \\(-\\theta\\) (VI.6 (D), four lines from \\(U\\) real and orthogonal and the field a turn of the pair; on a ring to \\(10^{-15}\\), against the same-charge walker off by \\(0.6\\); two holes against two walkers, and a pair made from the vacuum against the pair with walker and hole exchanged, the same — `hole_many_body_check.py` 1–4, plain lists of configuration amplitudes). (6)–(7) are IV.10’s bottom row — the positive half forward, the negative half backward — derived from exclusion plus the configuration, not chosen; the amplitude onto a sea column is \\(0\\) (ibid. 3, to \\(10^{-15}\\); and by brute force on 4 sites, the rank-5 antisymmetric column evolved index by index — the plain list of 56 configuration amplitudes, no determinant — to \\(10^{-15}\\)). The antiparticle is a missing walker; III.3’s backward crossing, read on the cut, is a hole; and a missing walker is, in every amplitude run, one walker of the opposite orientation and charge. The sea then enters in two ways only: through the split \\(\\Pi\_\\pm\\) on each line — which is a property of one walker’s column, the sense in which its clock turns (III.8c (15)–(16)) — and through the common factor \\(\\langle\\text{sea}|W|\\text{sea}\\rangle\\). Nothing has to be absent.

**(D) The vacuum’s energy is a sum over closed loops.** \\(\\langle\\text{sea}|U^t|\\text{sea}\\rangle=\\det\\Phi^\\dagger U^t\\Phi=e^{-iE_0t}\\) with \\(E_0=\\sum\_{\\varepsilon_j\<0}\\varepsilon_j=-\\tfrac12\\operatorname{Tr}|H|\\), and with \\(|\\varepsilon|=\\frac\\pi2-\\frac4\\pi\\sum\_{k\\ \\rm odd}\\frac{\\cos k\\varepsilon}{k^2}\\),

\\\[E_0=-\\frac{\\dim}{4}\\pi+\\frac2\\pi\\sum\_{k\\ \\rm odd}\\frac{\\operatorname{Tr}U^k}{k^2},\\tag{8}\\\]

8. \\(\\operatorname{Tr}U^k\\) the signed count of closed \\(k\\)-tick zigzags. The companion document’s “coin’s version”, line 2, is this, with the sea as its reason. On this tick an odd number of ticks moves a walker by an odd number of sites, so the closed zigzags of odd \\(k\\) are the ones that wind round the ring. On a ring with an even number of sites there are none: the sum is zero and \\(E_0=-\\tfrac{\\dim}4\\pi\\) exactly (\\(-62.831853\\) on 40 sites, ibid. 4, which therefore checks the first term only). On 7 sites the first is at \\(k=7\\) and \\(E_0=-10.950479=-10.995574+0.045095\\), the last number the loop sum (`review_v33_loop_sum_on_odd_ring.py`). The first term is \\(-\\tfrac\\pi4\\) per state, whatever the ring; the loop sum is the part that knows the ring’s length. The sign of a closed history is in the same count: around the selected configuration a history that winds \\(n\\) times before closing — of a walker above the sea, or of a hole — carries the sign of the \\(n\\)-cycle it closes, \\((-1)^{n-1}\\) (VI.5 (A)): \\(+1\\) for a single turn, \\(-1\\) for two exchanging places. The smallest loop a field closes is such a 2-cycle, the two-vertex bubble of (F) with its exclusion \\(-1\\), and that case is machine-checked there. This is where the \\(-1\\) per loop comes from; that it is Feynman’s \\(-1\\) for a loop of any number of vertices, in his convention for the backward leg (VI.6 (7)), is not checked here beyond two.

**(E) Pair creation is the turn, between configurations.** V.2’s turn, \\(e^{i\\theta}\\) on one site for one tick, takes the sea to “the sea with column \\(o\\) replaced by an empty column \\(e\\)” with amplitude

\\\[i\\theta\\,N\_{eo}\\,e^{-i(\\varepsilon_e-\\varepsilon_o)},\\qquad N\_{eo}=\\sum_c\\overline{e_c(x_0)}\\,o_c(x_0),\\tag{9}\\\]

9. to first order in \\(\\theta\\) and relative to the vacuum’s own factor (\\(3\\times10^{-5}\\) of the largest amplitude at \\(\\theta=0.02\\)) — these are the lines of the plates runs; and \\(\\sum\_{\\rm pairs}|{\\rm amp}|^2=1-|\\text{vacuum persistence}|^2\\) to relative order \\(\\theta^2\\) (\\(2\\times10^{-4}\\); ibid. 5). In the continuum, with III.15e (C)’s projectors, the element is \\(\\mathbb P\_+(k)\\,\\mathbb J\\,\\mathbb P\_-(k')\\) — lift a lagging loop at \\(k'\\) to a leading loop at \\(k\\); zero for a uniform turn; its transpose, \\(\\mathbb P\_-\\mathbb J\\mathbb P\_+\\), is the reverse process (`audit_2026_11_09.py` C). By (7) the loop that was removed is a leading loop \\(\\lambda\_{-k'}\\) of charge \\(-e\\): one loop re-oriented, counted by the forward-time observer as two. A pair is never a new walker: one walker steps up from the sea; the forward-time observer sees a particle above and a hole below, and the walker count is unchanged — III.3, charge conserved while crossings are not.

**(F) Why a sea at all, which half, and what it is needed for.** Three statements, each with its check.

*Half, not both and not neither.* Neither filled: the lamp. Both filled: nothing can move, and a hole in it shines — it is the empty case with the added detectors counted with the other sign (V.2 (H): \\(J\\to-J\\) exchanges filled with empty and \\(e\\to-e\\)). Half: the one dark configuration, the one that operation takes to itself, and the one the long-clock count selects (C).

*Which half is a name.* The count selects the columns that grow along the walk’s clock. That they are called “negative” is the convention \\(\\psi\\propto e^{-iEt/\\hbar}\\) together with an orientation already in the document: III.15b’s cut has a future side — the source fires before the counter clicks — and the carrier takes energy away in that direction. Fill the positive columns and run the walk with \\(t\\) reversed (\\(U\\to U^{\\mathsf T}\\)): every line of (6)–(9) reads the same with particle and hole exchanged and \\(e\\to-e\\). “Lower” is the direction the record defines; the selection does not depend on the name.

*The sea is needed only for its own weight.* By (7), every statement about a hole is a statement about one walker of the other orientation and charge, and the sea enters the amplitudes of walkers and holes in two ways only: the split \\(\\Pi\_\\pm\\) on each line, and the common factor \\(\\langle\\text{sea}|W|\\text{sea}\\rangle\\). The split is a property of one walker — the sense in which its clock turns (III.8c: \\(u\\pm iJu=2u\_{(\\pm)}\\)) — and with a field on it is applied at every vertex: (11)’s loop term. So for every process of walkers that exist — emitted, scattered, their antiparticles — nothing has to be stored: each ribbon carries the sense of its own clock, and the rule is (6)–(7). The one place the sea itself is measured is where that common factor is: a field pushing on nothing, and walls. A lone charge at the origin, its potential by VI.2’s time-link angle, three vacua, on the document’s continuum lines (`vacuum_polarization_sea_or_nothing.py`): (i) the sea, whose charge response is the sum over its lines; (ii) nothing plus the one-walker propagator with Feynman’s rule closed into a loop through two vertices, with the exclusion \\(-1\\); (iii) nothing. Then

\\\[\\chi\_{\\rm sea}=\\chi\_{\\rm loop}\\ \\ (3\\times10^{-9}),\\qquad\\chi\_{\\rm nothing}=0,\\qquad\\chi(k)=-\\alpha k^2+O(k^4),\\ \\chi(0)=0,\\qquad\\alpha a^2\\to\\frac1{6\\pi}\\ \\ (0.0647,\\,0.0572,\\,0.0547,\\,0.0537),\\tag{10}\\\]

10. the vacuum polarizes and does not screen a uniform charge; around the point charge an opposite-sign cloud, half of it within two Compton lengths, totalling \\(-(1-1/\\epsilon)\\) with \\(\\epsilon=1+q^2\\alpha\\); and \\(\\alpha\\to1/6\\pi m^2\\) is the one-loop vacuum polarization of \\(1{+}1\\) QED, \\(\\Pi(0)=q^2/6\\pi m^2\\). Nothing screens nothing; nothing plus the closed loop is the sea; and the closed loop is what nature has (the Uehling part of the Lamb shift). (The lattice tick’s own filled reference gives zero here; that is the tick’s sublattice — odd powers of \\(U\\) never return to their site — and not a property of the lines. `return_as_pair_response_test.py` 4.)

**(G) The ribbon’s Z and the closed loop are exchange partners.** III.1’s ribbon, cut at constant \\(t\\), crosses several times; a Z crosses three — forward, backward, forward. In the \\(t\\)-reading that is one column with a negative-energy part, and III.8c (E) splits it: the forward branch has \\(v=+Ju\\), the backward branch \\(v=-Ju\\) — a particle and, counting the added detectors with the other sign (V.2 (H): \\(J\\to-J\\), \\(e\\to-e\\)), an antiparticle of positive energy and opposite charge; each a complete walker, read by all its detectors. Three operations are kept apart here: the sign of the added detectors, which fixes charge and energy and touches no route; the measurement return of III.8c (A), \\(F(t_1,t_2)=F(t_2,t_1)^{\\mathsf T}\\), which sums the return routes of an event for either charge and changes no charge; and the exchange below, which acts between two crossings on a cut. Now exclusion: the Z’s forward branch and the original ribbon are two forward crossings on one cut. Antisymmetrize them — let the Z’s branch continue to the detector and the original forward part join the backward branch, which closes — and the exchange partner of *a ribbon with a Z* is *a plain ribbon plus a closed ribbon*. On the lattice this is one identity: one walker in the sea with a static field on,

\\\[\\frac{\\langle\\text{sea}+X\\,|\\,U\_\\Phi^t\\,|\\,\\text{sea}+A\\rangle}{\\langle\\text{sea}\\,|\\,U\_\\Phi^t\\,|\\,\\text{sea}\\rangle}=\\underbrace{(U\_\\Phi^t)\_{XA}}\_{\\text{the ribbon alone, Z's in}}-\\underbrace{e_X^\\dagger U\_\\Phi^t\\Phi\\,(\\Phi^\\dagger U\_\\Phi^t\\Phi)^{-1}\\Phi^\\dagger U\_\\Phi^te_A}\_{\\text{its exchange with the sea: the loop}},\\tag{11}\\\]

11. (`z_graph_and_closed_loop_are_exchange_partners.py` 1–3.) Field off, the loop term is \\((U^t\\Pi\_-)\_{XA}\\), (6)‘s blocked half. Field on, the ribbon alone ends partly in occupied columns (\\(5\\times10^{-5}\\) to \\(2\\times10^{-3}\\) for \\(q=0.1\\) to \\(0.4\\)) — exclusion violated — and the loop term removes that weight exactly (\\(10^{-16}\\)); its field-dependent part has the Z’s size at every \\(q\\) (ratio \\(1.000\\), \\(1.001\\), \\(1.004\\)). Two different ribbons stay two: \\(\\langle\\text{sea}+X+Y|U\_\\Phi^t|\\text{sea}+A+B\\rangle/\\langle\\text{sea}|U\_\\Phi^t|\\text{sea}\\rangle=\\det\\begin{pmatrix}G\_{XA}&G\_{XB}\\\\G\_{YA}&G\_{YB}\\end{pmatrix}\\) of the dressed amplitudes (11), to \\(10^{-17}\\); their mutual exchange is VI.1’s \\(E\\), each one’s loop is inside its own \\(G\\), and neither is the other’s loop (the undressed determinant is off by \\(0.016\\)). So the closed loops are not a second kind of thing: they are the ribbons’ own Z’s, seen across an exchange, with the exchange sign; and the sea is the cut-reading of all Z’s with all their partners. The smallest case, by hand, is VI.1a (C): a Z whose trough creates into the state a second ribbon occupies has \\(D-E=0\\).

**(H) What went in.**

*Taken:* the walker is read by a detector for each time sense and each orientation — eight in three dimensions, four in one (III.8c (B), IV.5a; input 11). The lines do not force the added detectors: the original ones are a closed set (IV.5a (C)). They are what rest and charge need. V.1–V.2 show that a walker a field can *turn* must be read by all of them; one in the original detectors alone can still be coupled through its coin, so this is “charged \\(\\Rightarrow\\) all the detectors”, not “detectable \\(\\Rightarrow\\) all the detectors”.

*Input:* exclusion (input 15), one process per state, for every column.

*Taken:* that the vacuum is the configuration the long-clock count selects ((C); VI.5).

*Derived from these:* the determinant and its two operations (A)–(B), with their algebra; Feynman’s rule (6)–(7); the hole as one walker of the opposite orientation and charge, exactly and with a field on (7), one hole by hand and two holes and a pair by machine; pair creation as one loop re-oriented (E); the vacuum’s energy as a sum over closed loops (8), the pair’s; its polarization (10); the Z and the closed loop as exchange partners (11).

*Not derived:* the name “negative” (an orientation, (F)); the carrier’s own state (VI.2 (F′)), which the near-field plates did without and the far field has not been run for.

This replaces, as a consequence and not as an input, IV.10’s bottom-row rule and the companion document’s vacuum prescription.

(`scripts/fock_space_from_the_walk.py` (A)–(E), 6 checks; `scripts/hole_is_opposite_loop.py` and `scripts/hole_many_body_check.py` for (7); `scripts/vacuum_polarization_sea_or_nothing.py` (F), 4 checks; `scripts/z_graph_and_closed_loop_are_exchange_partners.py` (G), 3 checks; `scripts/plate_of_walkers_sum_and_difference.py` 3, 6 for the lamp; `scripts/return_as_pair_response_test.py` 1–3: the measurement return of III.8c (A) pairs every out with every back whatever the added detectors hold, and makes no pairs; `scripts/audit_2026_11_09.py` C, D.)

## VI.5 The vacuum: what the walk’s own clock selects

*The picture first.* The walk has no sea. One walker has a kernel and a count (III.15d), and neither has a vacuum in it. Many walkers with exclusion have a count of configurations, and that count, taken along the walk’s own clock, does something the one-walker count cannot: it puts all its weight on one configuration. That configuration, described on a time cut, is “every negative-energy column occupied”. So the sea is not put under the walk and not chosen for being lowest: it is what the long-clock count of exclusive configurations is carried by. What is taken is that the vacuum *is* that configuration. The same count gives the sign of a closed history, and shows that without exclusion there would be no reference at all.

*Symbols.* \\(U\\) a one-walker evolution over some duration; \\(\\Lambda^N\\) the antisymmetric \\(N\\)-walker columns of VI.1 (D); \\(\\tau=x/c\\) the walk’s clock; \\(M_k\\), \\(L_k(\\tau)=e^{\\tau M_k}\\) and \\(\\omega_k\\) from III.15d. Nothing else is new.

**(A) A closed history is an exchange.** One walker’s closed histories are \\(\\operatorname{Tr}U=\\sum_sU\_{ss}\\): leave the state \\(s\\), come back to \\(s\\). For \\(N\\) identical walkers the closed histories are the trace over the antisymmetric columns, which by VI.4 (4)’s expansion is the sum of the \\(N\\times N\\) principal minors of \\(U\\),

\\\[\\operatorname{Tr}\_{\\Lambda^N}U=\\sum\_{|S|=N}\\det U\_{SS},\\qquad Z\\equiv\\sum\_{N\\ge0}\\operatorname{Tr}\_{\\Lambda^N}U=\\det(1+U).\\tag{1}\\\]

1. no walker, one, two, …, every slot: a finite sum over configurations, always defined (`gluing_sign_and_the_selected_reference.py` 1, built minor by minor). Its logarithm turns “all configurations” into “single closed histories”:

\\\[\\log\\det(1+U)=\\sum\_{n\\ge1}\\frac{(-1)^{n+1}}{n}\\operatorname{Tr}U^n .\\tag{2}\\\]

2. \\(\\operatorname{Tr}U^n\\) is one possible walker going round \\(n\\) times before closing. The \\(1/n\\) is where on the loop one starts. *The \\((-1)^{n-1}\\) is the sign of an \\(n\\)-cycle*: in the antisymmetric column, one walker going round \\(n\\) times is \\(n\\) walkers each stepping into the next one’s place. Gluing a history to its own start costs the exchange sign of the cycle it closes — VI.1 (D)’s \\(-1\\), on closed histories. For Part II’s undirected walker (the symmetric column, VI.1 (F)) the same sum is \\(1/\\det(1-U)\\) and every coefficient is \\(+1/n\\) (ibid. 2, by permanents).

*Where (2) holds.* It is the logarithm of (1) only while every eigenvalue of \\(U\\) is smaller than \\(1\\) in size. In (C) one slot has an eigenvalue \\(e^{a\\tau}>1\\) and the series as written diverges. The form that converges takes that slot out first:

\\\[\\log\\big(1+e^{a\\tau}\\big)=a\\tau+\\sum\_{n\\ge1}\\frac{(-1)^{n+1}}{n}e^{-na\\tau},\\qquad\\log Z\_{\\rm rest}=a\\tau+2\\sum\_{n\\ge1}\\frac{(-1)^{n+1}}{n}e^{-na\\tau}.\\tag{3}\\\]

3. the slot is occupied, which gives the \\(a\\tau\\); what winds after that is its *absence*, with weight \\(e^{-a\\tau}\\) per turn, exactly as an occupant of the other slot winds. The two series are a hole and a particle, each carrying the cycle sign. So the selection below is a statement about the finite sum (1), and (2) is the expansion around the selected configuration (`audit_gluing_draft_2026_11_09.py` 1).

**(B) In the \\(t\\)-reading nothing is selected.** At rest the walker turns in its own plane as \\(e^{\\mp at\\mathbb J}\\) in the two classes (III.15e (3)); with \\(\\mathbb J\\) written as \\(i\\) the two slots carry \\(e^{\\mp iat}\\) over a time \\(t\\), and

\\\[Z(t)=(1+e^{-iat})(1+e^{+iat})=2+2\\cos at:\\tag{4}\\\]

4. four configurations — none, one slot, the other, both — each a phase, none larger than the others (`gluing_sign_and_the_selected_reference.py` 3a). Every closed history in \\(t\\) is a turn, and a sum of turns prefers nothing.

**(C) Along the walk’s clock one configuration carries the count.** III.15d (3)–(4): along its own clock the walker’s weighted count is \\(L_k(\\tau)=e^{\\tau M_k}\\), with one direction growing as \\(e^{+\\omega_k\\tau}\\) and one shrinking as \\(e^{-\\omega_k\\tau}\\). Count the closed configurations of any number of walkers with it:

\\\[Z_k(\\tau)=\\det\\big(1+L_k(\\tau)\\big)=\\underbrace{1}\_{\\text{none}}+\\underbrace{e^{-\\omega_k\\tau}}\_{\\text{decaying slot}}+\\underbrace{e^{+\\omega_k\\tau}}\_{\\text{growing slot}}+\\underbrace{1}\_{\\text{both}}=2+2\\cosh\\omega_k\\tau .\\tag{5}\\\]

5. now the four configurations are weights, and one grows. Its share of \\(Z\\) at rest is \\(0.64\\) at \\(a\\tau=1.4\\), \\(0.993\\) at \\(5.6\\), \\(0.999998\\) at \\(14\\) (`sea_in_coordinates_check.py` 4), and

\\\[-\\frac1\\tau\\log Z_k(\\tau)\\ \\longrightarrow\\ -\\omega_k,\\qquad -\\frac\\hbar\\tau\\log\\prod_kZ_k(\\tau)\\ \\longrightarrow\\ -\\sum_k\\hbar\\omega_k=-\\tfrac12\\operatorname{Tr}|H| .\\tag{6}\\\]

6. (\\(-\\sum_k\\hbar\\omega_k\\) for five \\(k\\)’s to \\(10^{-6}\\) at \\(a\\tau=20\\), `gluing_sign_and_the_selected_reference.py` 3d.) The configuration that carries the count is “the growing slot occupied, the decaying one empty, for every \\(k\\)”. By III.15e (5) the growing direction of the walk’s count is, read in \\(t\\), the lagging loop — the class of negative energy. So in the words of the time cut the configuration is *every negative-energy column occupied, every positive one empty*: the sea of VI.4 (C), and (6) is its energy, VI.4 (8). The slots are those of a walker read by all four detectors; one in the original two alone has half of them (III.15e (B)). And by III.15d (C) every \\(k\\) is the rest configuration seen from one frame, \\(\\hbar k=mc\\sinh\\eta\\): the vacuum is the rest count’s growing slot, occupied once in every frame, with the box’s own measure over \\(k\\), \\(dk=(a/c)\\cosh\\eta\\,d\\eta\\). (This is not III.16’s sum over frames, which is uniform in \\(\\eta\\): per unit length \\(-\\sum_k\\hbar\\omega_k=-\\frac{\\hbar a^2}{2\\pi c}\\int\\cosh^2\\eta\\,d\\eta\\), an energy in a box of fixed length.)

No “lowest energy” was said. “Lowest” is what the long-clock limit of the count is called on the time cut.

**(D) Without the exchange sign there is no reference at all.** The same count for the undirected walker is \\(Z=1/\\det(1-L_k)\\). The growing slot contributes \\(1/(1-e^{+\\omega_k\\tau})\<0\\) for every \\(\\tau>0\\): the geometric series over how many walkers sit in that slot does not converge (ibid. 3c). A walker without exclusion and with a growing slot has no vacuum; one can always add another. It is one-per-slot that makes the count of closed configurations settle on a reference. The fermion’s negative reference energy and the existence of its vacuum are one fact.

**(E) What the count is a count of, and what the vacuum looks like.**

*The clock as an inverse temperature.* Many walkers with exclusion in contact with anything at temperature \\(T\\) have the weights \\(e^{-E\_{\\rm config}/k_BT}\\) and the partition function \\(\\prod_k(1+e^{-\\beta\\hbar\\omega_k})(1+e^{+\\beta\\hbar\\omega_k})\\). That is \\(\\prod_kZ_k(\\tau)\\) with \\(\\beta\\hbar\\leftrightarrow\\tau\\): an identity between two formulas, the clock in the place of \\(\\beta\\hbar\\). Cooling and running the clock long are the same limit.

*What the product over \\(k\\) is.* On a *time* cut a walker’s states are functions of \\(x\\), the modes labelled by \\(k\\) are a basis there, and VI.1 (D) gives one occupant per mode: that is the determinant. Each mode’s weight is then, by III.15d, a positive count along the walk’s clock of one ribbon in one frame. So the product is one ribbon per frame, each counted along its own clock, at most one occupant per frame and slot; ribbons of different frames do not meet. It is not one ensemble of non-crossing lines in the \\((x,t)\\) plane. A determinant with positive weights counts non-crossing lines when its rows and columns are *positions* on the cut; here they are the label \\(k\\), and the two are different objects:

|  | lines in the real \\((x,t)\\) plane | the product of (6) |
| --- | --- | --- |
| what diagonalises the clock step | \\(e^{-i\\nu t}\\), \\(\\nu\\) real | \\(e^{-ckt}\\) |
| step matrix | \\(\\begin{pmatrix}i\\nu&a\\\\a&-i\\nu\\end{pmatrix}\\) | \\(M_k=\\begin{pmatrix}-ck&a\\\\a&ck\\end{pmatrix}\\) |
| rates along the clock | \\(\\pm\\sqrt{a^2-\\nu^2}\\): growth only for \\(\\lvert\\nu\\rvert\<a\\) | \\(\\pm\\sqrt{a^2+c^2k^2}\\): every \\(k\\) grows |
| growth per unit length | \\(\\int\_{\\lvert\\nu\\rvert\<a}\\sqrt{a^2-\\nu^2}\\,\\frac{d\\nu}{2\\pi}=\\frac{a^2}4\\), bounded | \\(\\int\\sqrt{a^2+c^2k^2}\\,\\frac{dk}{2\\pi}\\), unbounded |

Real exponentials do not span the functions on a clock cut, so the right column is not the left one in another basis (`audit_gluing_draft_2026_11_09.py` 3). Nothing is packed in \\(t\\), and the walls of the companion document, which sit at two values of \\(x\\), change which \\(k\\) occur and confine nothing sideways.

*Two pictures of the one vacuum.* VI.1 (D) is exclusion in every state, and a state is any column; read in two bases it gives two pictures (`vacuum_picture.py`).

- *Frames, in the walk’s clock.* A slot is a sheet: the rest zigzag seen from one frame, spread evenly over \\(t\\). Sheets of different frames pass through one another. At most one unit per sheet; the vacuum is every growing sheet holding its unit. The weights are positive here, and this is where “occupied” is derived.
- *Positions, on a time cut.* Each site has two states, moving right and moving left. Two walkers of opposite sense may share a site (VI.1 (E)); two of the same sense may not. In the vacuum exactly half of all (site, sense) states are occupied in every sample, each with probability exactly \\(\\tfrac12\\) (from \\(\\overline{\\Pi\_-}=\\Pi\_+\\) and \\(\\Pi\_-+\\Pi\_+=1\\)), and walkers of one sense keep apart: on a ring of 12 sites with III.15a’s tick the same sense is present on two neighbouring sites with probability \\(0.197\\) against the \\(0.25\\) of independence. From the continuum lines, \\(h(k)=ck\\sigma_z-a\\sigma_y\\) with exact \\(k\\) on a ring, opposite senses are both present on one site with probability \\(0.199\\) at \\(a=0.6\\) and \\(0.242\\) at \\(a=0.15\\): the mass ties the two senses at a point, and a site holds exactly one walker more often than chance (\\(0.60\\), \\(0.52\\) against \\(0.5\\)) (`vacuum_snapshot_continuum_lines.py`). (On the lattice tick that last probability is exactly \\(0.2500\\); that is the tick’s sublattice — odd powers of \\(U\\) never return to their site — and not a property of the lines.) These are probabilities of the determinant squared; between cuts the walkers evolve by amplitudes that interfere, so this is what is there on a cut, not a classical gas to be followed.

**(F) Where the sea came from.** Step by step:

| step | one walker? | which reading | has a sea? |
| --- | --- | --- | --- |
| the kernel, III.15d (1) | yes | \\(x\\) (outside the cone) | no — the distribution of one zigzag |
| the weighted count, III.15d (4) | yes | \\(x\\) | no — a growth rate, \\(\\omega_k\\) |
| the configuration count (5) | **many**, with exclusion | \\(x\\) | not yet — a determinant |
| the configuration that carries it (6) | many | \\(x\\), then named on the time cut by III.15e (5) | **yes** |

The sea appears at one step: when the many-walker count with exclusion is asked which configuration grows fastest along the walker’s own clock. It is not in the one-walker equations in either reading, and not in the \\(t\\)-reading of many walkers (there the count is (4), four phases). It is in the reading that counts, because the other reading only turns. In coordinates no continuation is made anywhere: III.14’s kernel already lives outside the cone, and the weight is a Laplace transform with bounded support. The one quarter turn used is III.15e’s \\(\\mathcal D\\), between the two backward detectors, and it only renames the result.

*What went in.* Exclusion (input 15) and the completed detectors (III.8c (B), input 11). **Taken:** that the vacuum is the configuration the long-clock count selects, with the box’s measure over \\(k\\). Derived from these: the sign \\((-1)^{n-1}\\) of an \\(n\\)-fold closed history; that one configuration carries the count; that it is the lagging class in every frame; its energy \\(-\\tfrac12\\operatorname{Tr}|H|\\); that without exclusion there is no reference. “Energy” for a growth rate along the clock is III.8c (D)’s definition of energy as \\(\\hbar\\) times a turn rate, the two rates being one number by III.15e (4); the definition, with \\(\\hbar\\) the scale, is not derived. As a kernel on the whole line \\(e^{\\tau M_k}\\) has growing modes and its integral over \\(k\\) does not converge; every count here is on a finite system, as the companion document’s round trips are.

(`scripts/master_v33_new_text_checks.py` 6 for (3) and the shares; `scripts/gluing_sign_and_the_selected_reference.py`, 15 checks: (1)–(2), the symmetric count, (4), the shares and limits, (D), five \\(k\\)’s, \\(M_k^2=\\omega_k^2\\) and (5) at \\(k=0,\\,0.7,\\,2\\); `scripts/sea_in_coordinates_check.py` 4–5; `scripts/audit_gluing_draft_2026_11_09.py` 1 for (3), 3 for the table, 4 for the measure; `scripts/vacuum_picture.py`, `scripts/vacuum_snapshot_continuum_lines.py` for (E).)

## VI.6 The hole is the opposite loop

*The picture first.* A hole needs something to make a hole in, and the walk has only things that zigzag. VI.5 supplies the something: the configuration the long-clock count selects holds, in every frame, the loop of one orientation — the lagging loop — once. Take one of those loops out. What is left behaves, to the last digit and with a field on, as *one loop of the other orientation with the opposite charge*. Nothing is absent; one loop is turned over. This section does the removal as an operation, reads what it costs, and proves the statement free and with a field. It then reads Feynman’s rule at a source in the same words.

*Symbols.* \\(\\lambda\\), \\(\\ell\\) the leading and lagging loops, \\(\\mathbb P\_\\pm\\), \\(\\mathbb G\\), \\(\\mathbb J\\) of III.15e; \\(c^\\dagger_w\\), \\(c_w\\) the adding and removing of VI.4 (4); \\(U\\) the tick of VI.4 (1), real and orthogonal; \\(\\Pi\_\\pm\\) the projectors onto its columns of positive and negative \\(\\varepsilon\\); \\(W\_\\theta\\) the evolution with V.2’s turn \\(e^{i\\theta(x)}\\) of the pair on each site at each tick. Columns are complex: the original and the added detectors together, \\(\\mathbb J\\) written as \\(i\\) (IV.5a).

**(A) The count, and what “occupied” means.** Exclusion (VI.1 (D)): a slot holds \\(0\\) or \\(1\\). Per \\(k\\) there are two slots and four configurations, with weights along the clock \\(1\\), \\(e^{-\\omega_k\\tau}\\), \\(e^{+\\omega_k\\tau}\\), \\(1\\) (VI.5 (5)). At long clock one configuration carries the count: the growing slot occupied, the decaying one empty, for every \\(k\\); and every \\(k\\) is the rest configuration seen from one frame (III.15d (C)). By III.15e (5) the occupant is a lagging loop. So

\\\[|\\text{vac}\\rangle=\\ell\_{k_1}\\wedge\\ell\_{k_2}\\wedge\\cdots,\\qquad E_0=-\\sum_k\\hbar\\omega_k,\\tag{1}\\\]

1. “occupied” meaning *present in the configuration that carries the count*. Two slots per \\(k\\), each a complex direction, is the count of the walker read by all four detectors; a self-conjugate one, in the original two alone, has half as many and half the energy (\\(-0.5\\) against \\(-1.0\\) in units of \\(\\hbar a\\) at rest, with its own lines reproduced as its equations of motion; `audit_2026_11_09.py` D). The companion document’s \\(-\\pi\\hbar c/24L\\) is the pair’s number.

**(B) Removing a loop, and what it costs.** With VI.4 (4): two loops, \\(T\_{XY}=\\ell_Xm_Y-m_X\\ell_Y\\) with \\(\\ell\\perp m\\); removing \\(\\ell\\) gives \\(\\sum_X\\overline{\\ell_X}(\\ell_Xm_Y-m_X\\ell_Y)=m_Y\\). Removing is taking a column out of the determinant. Relative to the vacuum’s own factor \\(\\langle\\text{vac}|U^t|\\text{vac}\\rangle\\), on a ring of four sites with the Fock space built index by index,

\\\[c^\\dagger\_\\lambda|\\text{vac}\\rangle\\ \\to\\ e^{-i\\omega t},\\qquad c\_\\ell|\\text{vac}\\rangle\\ \\to\\ e^{-i\\omega t}\\qquad(\\text{both of modulus }1.000000000000).\\tag{2}\\\]

2. adding a leading loop costs \\(+\\hbar\\omega\\). Removing a lagging loop also costs \\(+\\hbar\\omega\\): the vacuum’s factor is a product that contained that loop’s \\(e^{+i\\omega t}\\), and taking the loop out takes the factor out. On the walk’s side the same two changes each multiply the count by \\(e^{-\\omega\\tau}\\) (`hole_is_opposite_loop.py` 3).

**(C) The hole propagates as the opposite loop: the free walk.** Remove at \\(A\\), run \\(t\\) ticks, ask for the removal at \\(X\\):

\\\[\\frac{\\langle c\_{e_X}\\text{vac}\\,|\\,U^t\\,|\\,c\_{e_A}\\text{vac}\\rangle}{\\langle\\text{vac}|U^t|\\text{vac}\\rangle}=\\overline{(U^t\\Pi\_-)\_{XA}}=(U^t\\Pi\_+)\_{XA}.\\tag{3}\\\]

3. the first equality is VI.4 (7), here from the plain list of configuration amplitudes. The second is one line: \\(U\\) is real, so the conjugate of a column with \\(U\\varphi=e^{-i\\varepsilon}\\varphi\\) has \\(U\\overline\\varphi=e^{+i\\varepsilon}\\overline\\varphi\\), and \\(\\overline{\\Pi\_-}=\\Pi\_+\\). In the walk’s words conjugation is \\(u+iv\\to u-iv\\), the added detectors counted with the other sign, \\(v\\to-v\\) of III.15e (3): it turns a lagging loop into a leading one, and it reverses the winding, \\(k\\to-k\\). *The amplitude for a hole to go from \\(A\\) to \\(X\\) is the amplitude for a leading loop to go from \\(A\\) to \\(X\\)* (\\(-0.056801\\,i\\) three ways; ibid. 2).

**(D) With a field: the opposite charge, exactly.** \\(W\_\\theta\\) is the tick with V.2’s turn, so \\(W\_{-\\theta}=\\overline{W\_\\theta}\\). \\(\\Phi\\) the matrix whose columns are the sea’s, \\(P=\\Phi\\Phi^\\dagger\\) the projector onto them, \\(Q=1-P=\\overline P\\).

*The hole.* The many-walker evolution carries a removal along, \\(\\mathcal W\\,c_w\\,\\mathcal W^{-1}=c\_{W\_\\theta w}\\), and the one-walker transition density between two determinants gives

\\\[H\_{AX}\\equiv\\frac{\\langle c\_{e_X}\\text{vac}\\,|\\,W\_\\theta\\,|\\,c\_{e_A}\\text{vac}\\rangle}{\\langle\\text{vac}|W\_\\theta|\\text{vac}\\rangle}=\\big\[\\Phi\\,(\\Phi^\\dagger W\_\\theta\\Phi)^{-1}\\Phi^\\dagger\\big\]\_{AX}:\\qquad H=(PW\_\\theta P)^{-1}\\ \\text{on the sea's columns}.\\tag{4}\\\]

*A walker of the other charge, above the same sea* (VI.4 (11) with \\(-\\theta\\)):

\\\[G'=W\_{-\\theta}-W\_{-\\theta}\\Phi\\,(\\Phi^\\dagger W\_{-\\theta}\\Phi)^{-1}\\Phi^\\dagger W\_{-\\theta}=\\big\[\\big(W\_{-\\theta}^{-1}\\big)\_{QQ}\\big\]^{-1}\\tag{5}\\\]

5. the Schur complement of the sea block is the inverse of the other block of the inverse. Now \\(W\_{-\\theta}\\) is unitary and equals \\(\\overline{W\_\\theta}\\), so \\(W\_{-\\theta}^{-1}=W\_\\theta^{\\mathsf T}\\); and \\(Q=P^{\\mathsf T}\\), so \\((W\_\\theta^{\\mathsf T})\_{QQ}=(PW\_\\theta P)^{\\mathsf T}\\). Hence

\\\[\\boxed{\\;H\_{AX}(\\theta)=G'\_{XA}(-\\theta)\\;}\\tag{6}\\\]

6. the hole in the field \\(\\theta\\) is the walker in the field \\(-\\theta\\). Three things went in: the free step is real (III.15a), the field is a turn of the pair (V.2), and the step is orthogonal. On the ring, \\(\\theta=(0.3,-0.2,0.5,0.1)\\), \\(t=3\\): the two sides agree to \\(10^{-15}\\) in every entry; against the walker in \\(+\\theta\\) they differ by \\(0.60\\); the plain-list amplitude equals the determinant formula (\\(-0.086441-0.004988\\,i\\)); and the vacuum’s own factor has modulus \\(0.902\\) — the field is making pairs, and the identity is exact all the same (`audit_2026_11_09.py` E). With a field that varies from site to site and from tick to tick the same holds for one hole at all 64 \\((A,X)\\), for two holes against two walkers, and for a pair made from the vacuum against the pair with walker and hole exchanged, \\(10^{-15}\\) each (`hole_many_body_check.py` 1–4). This is V.2 (H)’s \\(J\\to-J\\), \\(e\\to-e\\), now for the removed loop.

**So \\(c\_{\\ell_k}|\\text{vac}\\rangle\\) is, in every amplitude run, one leading loop \\(\\overline{\\ell_k}=\\lambda\_{-k}\\) of charge \\(-q\\).** The sea enters in two ways only: through the split \\(\\Pi\_\\pm\\) on each line, which is one walker’s own property — the sense in which its clock turns — and with a field on is applied at every vertex (VI.4 (11)’s loop term); and through the common factor \\(\\langle\\text{vac}|W|\\text{vac}\\rangle\\).

**(E) What the field does, in these words.** A field with a gradient has \\(\\mathbb P\_+(k)\\,\\mathbb J\\,\\mathbb P\_-(k')\\neq0\\) for \\(k\\neq k'\\) (III.15e (C); `audit_2026_11_09.py` C). On the vacuum that element is \\(c^\\dagger\_{\\lambda_k}c\_{\\ell\_{k'}}\\): remove a lagging loop at \\(k'\\), add a leading loop at \\(k\\). By (D) the removed loop is a leading loop \\(\\lambda\_{-k'}\\) of charge \\(-q\\). One loop is re-oriented; the forward-time observer counts two leading loops, of charges \\(q\\) and \\(-q\\), the field supplying the momentum \\(k-k'\\). The count’s weight changes by \\(e^{-(\\omega_k+\\omega\_{k'})\\tau}\\). A uniform turn does nothing: \\(\\mathbb P\_+\\mathbb J\\mathbb P\_-=0\\) at \\(k=k'\\), and a constant turn commutes with the step and turns every pair by the same angle (VI.2 (A)). The transpose, \\(\\mathbb P\_-(k)\\,\\mathbb J\\,\\mathbb P\_+(k')\\), is the reverse process.

**(F) The same rule read at a source.** The evolution question and the source question have two kernels (III.15, IV.10). With \\(\\Pi\_\\pm(k)=\\tfrac12\\big(1\\pm iG(k)/\\omega_k\\big)\\) the projectors of III.15e (C) in complex form,

\\\[K=\\theta(t)\\big\[e^{-i\\omega t}\\Pi\_++e^{+i\\omega t}\\Pi\_-\\big\],\\qquad K_F=\\theta(t)\\,e^{-i\\omega t}\\Pi\_+-\\theta(-t)\\,e^{+i\\omega t}\\Pi\_- .\\tag{7}\\\]

7. both solve the same equation with a point source — the jump \\(\\Pi\_++\\Pi\_-=1\\) at \\(t=0\\) fixes the relative sign — and they differ by a free solution: where the \\(\\Pi\_-\\) leg lies. In \\(K\\) both legs leave the source upward. In \\(K_F\\) the \\(\\Pi\_-\\) leg points down, and \\(e^{ikx}e^{+i\\omega t}=e^{-i\[(-k)x-\\omega t\]}\\) is a leading loop of momentum \\(-k\\) run the other way in \\(t\\). By (C)–(D) that leg, read forward, is a leading loop of charge \\(-q\\) arriving at the source. III.3’s charge then says which kernel a source of one walker has. One walker injected, charge \\(q\\); weights \\(w\_\\pm=\\|\\Pi\_\\pm\\Psi(0)\\|^2\\), \\(w\_++w\_-=1\\):

- both legs leaving: the charge that leaves is \\(q\\,(w\_+-w\_-)\\neq q\\);
- the leading leg leaving, the other arriving and absorbed: the change is \\(q\\,w\_++q\\,w\_-=q\\).

So with the removed loop read as a walker of charge \\(-q\\), conservation of III.3’s charge at the source puts the \\(\\Pi\_-\\) leg in the past. That is Feynman’s rule at a source, and it agrees with VI.4 (6)–(7), where it came from exclusion and the vacuum. In the forward-time picture: the antiparticle’s leg comes from the past into the source and is absorbed, the particle’s leg leaves, and the charge \\(q\\) is injected once. A Z in a field is the same reading of (E): a pair at an earlier point of the field, an annihilation at a later one, three lines between — III.3’s trough and crest, with the trough before the crest in \\(t\\).

*What went in.*

| item | status |
| --- | --- |
| the walker is read by all its detectors | taken (III.8c (B), IV.5a); needed for a walker a field can turn (V.1–V.2), not for every walker that can be detected |
| the vacuum is the configuration that carries the count | taken (VI.5); the count is the pair’s |
| adding, removing, their algebra | VI.4 (B) |
| the hole is one leading loop of charge \\(-q\\), free and with a field on | **derived**, (C)–(D) |
| pair creation as one loop re-oriented | derived, (E) |
| the \\(\\Pi\_-\\) leg in the past at a source | derived from III.3’s charge, (F); and from exclusion, VI.4 (C) |

Not derived: whether the orientation of the cut, the name “negative”, can be traced to the carrier.

(`scripts/hole_is_opposite_loop.py` 0–4: \\(\\overline{\\Pi\_-}=\\Pi\_+\\), (3), (2) with the moduli, the conjugate of a lagging loop; `scripts/audit_2026_11_09.py` C, D, E; `scripts/hole_many_body_check.py` 1–4; `scripts/master_v33_new_text_checks.py` 7–8 for (F). Finite matrices and exact enumeration only.)

## VI.7 Leading loop, lagging loop, hole: the settled count

*The picture first.* A walker flips between its two senses at the rate \\(a\\). Follow it along its own clock and one thing fades and one thing stays: which sense it is in is forgotten, at the rate \\(2a\\), and what stays is the walker with equal weight in both senses. The two are the two combinations of III.1’s two counters, \\(Q\_+\\) and \\(Q\_-\\), in the walk’s clock — the difference and the sum. (These are not the sums and differences of IV.5b, which are of the forward and the backward detector of one orientation.) Read in \\(t\\) the same two combinations do not fade or stay; they turn, in opposite senses, and are the leading loop and the lagging loop of III.15e. So the lagging loop is the flip process at equilibrium, and the leading loop is the memory of the sense. With many walkers and exclusion every slot settles the same way, once: that configuration is VI.5’s vacuum, and a gap in it is VI.6’s hole, the antiparticle. The names are the document’s: *leading loop* (the walker whose clock turns one way: at rest, the same amount in the forward and the backward detector of each orientation; positive energy), *lagging loop* (the walker whose clock turns the other way: at rest, opposite amounts; negative energy — the leading loop reversed), *hole* (a lagging loop missing from the settled configuration; it moves as the lagging loop traversed in reverse, VI.6 (C)).

*Symbols.* \\(P\_\\pm\\) III.2’s two probabilities, \\(Q\_\\pm\\) the same with the loss out (III.5a (F)); \\(\\tau=x/c\\) the walk’s clock; \\(\\mathbb A,\\mathbb S,\\mathbb T,\\mathbb J\\), \\(\\mathbb G(k)\\), \\(\\mathbb P\_\\pm(k)\\), \\(\\lambda\\), \\(\\ell\\) of III.15e, order \\((u\_+,u\_-,v\_+,v\_-)\\); \\(M_k\\), \\(\\eta\\) of III.15d; \\(\\varepsilon\\) a one-walker energy on a ring (VI.4). \\(\\mathbb K\\) counts the added detectors with the other sign, \\(v\\to-v\\) of III.15e (3) and V.2 (H), written as a matrix.

**(A) Two reversals of a loop.** A leading loop can be reversed in two ways: count the added detectors with the other sign, or run the process backward.

\\\[\\mathbb K=\\begin{pmatrix}1&0&0&0\\\\0&1&0&0\\\\0&0&-1&0\\\\0&0&0&-1\\end{pmatrix}\\quad(v\\to-v:\\ \\ v=Tu\\ \\text{becomes}\\ v=T^{-1}u),\\qquad\\mathbb S=\\begin{pmatrix}0&1&0&0\\\\1&0&0&0\\\\0&0&0&1\\\\0&0&1&0\\end{pmatrix}\\ \\text{together with}\\ t\\to-t .\\tag{1}\\\]

1. \\(\\mathbb K\\) commutes with \\(\\mathbb A\\), \\(\\mathbb S\\), \\(\\mathbb T\\) and anticommutes with \\(\\mathbb J\\); \\(\\mathbb S\\) exchanges the forward and the backward counter and commutes with \\(\\mathbb J\\). On a mode,

\\\[\\mathbb K\\,\\mathbb G(k)\\,\\mathbb K=\\mathbb G(-k),\\qquad\\mathbb S\\,\\mathbb G(k)\\,\\mathbb S=-\\mathbb G(k);\\qquad\\mathbb K\\,\\mathbb P\_+(k)\\,\\mathbb K=\\mathbb P\_-(-k),\\qquad\\mathbb S\\,\\mathbb P\_+(k)\\,\\mathbb S=\\mathbb P\_-(k).\\tag{2}\\\]

2. each carries the leading class onto the lagging class: \\(\\mathbb K\\lambda=\\ell\\). \\(\\mathbb K\\) works because the step in \\(x\\), \\(-c\\,\\mathbb A\\,\\partial_x+a\\,\\mathbb T\\), contains no \\(\\mathbb J\\). The product \\(\\mathbb K\\mathbb S\\) keeps the class and reverses the motion, \\(\\mathbb P\_\\pm(k)\\to\\mathbb P\_\\pm(-k)\\): the two reversals are one, up to turning the motion round. With V.2’s field, a turn \\(\\varphi\\,\\mathbb J\\) of the pair, both reverse the field: \\(\\mathbb K\\,\\mathbb G\_\\varphi(k)\\,\\mathbb K=\\mathbb G\_{-\\varphi}(-k)\\). Two consequences. *The lagging loop is the leading loop of the opposite charge.* And *the energies pair off*: the energy matrix \\(\\hbar\\,\\mathbb J\\,\\mathbb G(k)\\) has trace \\(0\\), and on a ring the tick’s energies come in pairs \\(\\pm\\varepsilon\\) and sum to \\(0\\) (4, 7, 12 sites, to round-off; with a field on, the energies in \\(\\theta\\) are minus those in \\(-\\theta\\)).

**(B) The flip process: what stays and what fades.** III.2’s lines at rest, as probabilities, loss in:

\\\[c\\,\\partial_x\\begin{pmatrix}P\_+\\\\P\_-\\end{pmatrix}=a\\begin{pmatrix}-1&1\\\\1&-1\\end{pmatrix}\\begin{pmatrix}P\_+\\\\P\_-\\end{pmatrix},\\qquad\\begin{pmatrix}P\_+\\\\P\_-\\end{pmatrix}(\\tau)=\\underbrace{\\frac12\\begin{pmatrix}1\\\\1\\end{pmatrix}}\_{\\text{the result}}+\\underbrace{\\frac12\\,e^{-2a\\tau}\\begin{pmatrix}1\\\\-1\\end{pmatrix}}\_{\\text{what is still left to go}}\\quad\\text{from }\\begin{pmatrix}1\\\\0\\end{pmatrix}.\\tag{3}\\\]

3. the sum has rate \\(0\\): it is the equilibrium, the walker whose sense is no longer known. The difference has rate \\(-2a\\): it is the memory of the sense. “The flip process at equilibrium” and “the memory of the sense fades” are one statement, the result and the process, as \\(1\\) and \\(0.9+0.09+0.009+\\dots\\) are: the first term of (3) is the sum, and the second is what remains to be added, \\(0.1,\\ 0.01,\\ 0.001,\\dots\\) Taking the loss out with the one factor \\(e^{a\\tau}\\) (III.5a (F)) moves both rates by \\(a\\), to \\(+a\\) and \\(-a\\), and changes no difference: the gap is \\(2a=2mc^2/\\hbar\\) either way. Reading in \\(t\\) changes the sign of one turn (III.15a (A)), which is one multiplication:

\\\[\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix}\\begin{pmatrix}0&1\\\\1&0\\end{pmatrix}=\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}:\\qquad\\cosh a\\tau,\\ \\sinh a\\tau\\ \\ \\longrightarrow\\ \\ \\cos at,\\ \\sin at .\\tag{4}\\\]

4. the same flips at the same rate; along the clock they make a count grow or shrink, in \\(t\\) they turn the pair. By III.15e (5) the two combinations are the two loops:

|  | along the walk’s clock | in \\(t\\) |
| --- | --- | --- |
| sum \\((1,1)\\) | stays (rate \\(0\\); \\(+a\\) with the loss out) | lagging loop, \\(E=-\\hbar a\\) |
| difference \\((1,-1)\\) | fades (rate \\(-2a\\); \\(-a\\) with the loss out) | leading loop, \\(E=+\\hbar a\\) |
| long run | settles on the sum | never settles: keeps turning |

Every row but the last is one fact in two readings. “Negative energy” is the plain count of all paths, which has no cancellation in it; “positive energy” is the count taken with a sign on every flip, even minus odd. The weight along the clock is \\(e^{-E\\tau/\\hbar}\\) for both.

**(C) Blocks: the visible flip rate.** Cut the clock into blocks of length \\(\\tau\\) and record only whether the walker leaves a block in the sense opposite to the one it entered with. Two flips inside a block cancel. From (3),

\\\[P\_{\\rm flip}(\\tau)=\\frac{1-e^{-2a\\tau}}2,\\qquad a\_{\\rm block}(\\tau)=\\frac{1-e^{-2a\\tau}}{2\\tau}\\ \\longrightarrow\\ a\\ \\ (\\tau\\to0),\\qquad\\longrightarrow\\ \\frac1{2\\tau}\\ \\ (\\tau\\to\\infty).\\tag{5}\\\]

5. an interval holding \\(aT=20\\) flips on the average, cut into \\(N\\) blocks, \\(6000\\) walks:

| blocks \\(N\\) | 2000 | 100 | 20 | 10 | 4 | 1 |
| --- | --- | --- | --- | --- | --- | --- |
| blocks that come out reversed, simulated | 19.79 | 16.49 | 8.61 | 4.91 | 2.01 | 0.51 |
| \\(N\\,(1-e^{-2aT/N})/2\\) | 19.80 | 16.48 | 8.65 | 4.91 | 2.00 | 0.50 |

Zoomed in, all twenty are seen; zoomed out, half a flip — the sense is forgotten. The fall of the visible flip rate is the second term of (3) dying; what it leaves is the first term. So the falling flip rate is not the lagging loop itself: it is the leading loop’s part fading, and the lagging loop is what is left. The flips have not gone: the count of paths in the block is still \\(e^{a\\tau}\\), every flip in it. This is a statement in the clock, where the weights are probabilities. In \\(t\\) the same block has the flip weight \\(\\sin^2a\\tau\\), which does not settle.

**(D) In motion: the same equilibrium, and the speed.** A moving observer sees the flip rate \\(a\\,e^{\\eta}\\) while the walker runs forward and \\(a\\,e^{-\\eta}\\) while it runs backward (III.15d (D)), a process again:

\\\[c\\,\\partial\_{x'}\\begin{pmatrix}P\_+\\\\P\_-\\end{pmatrix}=\\begin{pmatrix}-a\\,e^{\\eta}&a\\,e^{-\\eta}\\\\a\\,e^{\\eta}&-a\\,e^{-\\eta}\\end{pmatrix}\\begin{pmatrix}P\_+\\\\P\_-\\end{pmatrix}:\\qquad\\text{rates }0\\text{ and }-2\\omega_k;\\qquad\\text{equilibrium }\\frac1{2\\cosh\\eta}\\begin{pmatrix}e^{-\\eta}\\\\e^{\\eta}\\end{pmatrix}.\\tag{6}\\\]

6. III.15d (12). The gap is \\(2\\omega_k\\). At equilibrium the two senses are not equally weighted, and their imbalance is \\(\\tanh\\eta=v/c\\): *the mean speed of the walker is the imbalance of senses of the settled state*. Every segment is at \\(c\\); a block long against \\(1/\\omega_k\\) moves at \\(v\\). The flips per unit clock at equilibrium are \\(a/\\cosh\\eta=a/\\gamma\\), while the energy is \\(\\hbar a\\gamma\\): the energy is the flip rate only at rest; in general it is the rate of the count. The equilibrium column is \\(M_k\\)’s growing column, both entries of one sign at every \\(k\\); the fading one has a change of sign. Blocks: reversed with probability \\(2\\pi\_+\\pi\_-(1-e^{-2\\omega_k\\tau})\\) (\\(\\eta=0.8\\), \\(400\\) to \\(2\\) blocks: \\(13.98\\), \\(8.25\\), \\(2.77\\), \\(0.56\\) simulated against \\(14.00\\), \\(8.25\\), \\(2.78\\), \\(0.56\\)).

**(E) Many walkers: the settled configuration, and the hole.** A slot holds \\(0\\) or \\(1\\) (VI.1 (D)). The many-walker step is built from the minors of the tick (VI.5 (1)), so the factor of a configuration is the product of the factors of the loops in it, and its energy is the sum of theirs. With no walker the product is empty: factor \\(1\\), energy \\(0\\). On four sites, eight slots, the \\(256\\) configurations built minor by minor:

\\\[E(\\text{configuration})=\\sum\_{\\text{slots occupied}}\\varepsilon,\\qquad E(\\text{none})=0,\\qquad E(\\text{all})=0\\ \\ (\\det U=1),\\qquad E_0=\\sum\_{\\varepsilon\<0}\\varepsilon=-\\tfrac12\\sum\_{\\text{all slots}}|\\varepsilon| .\\tag{7}\\\]

7. the lowest of the \\(256\\) is unique: every lagging loop present, no leading loop. Its energy is minus half the sum of the sizes of *all* the slots’ energies, present or not, because by (A) they pair off. The \\(256\\) energies are symmetric about \\(0\\): a configuration and its complement have opposite energies. For one pair of slots \\(\\pm\\varepsilon\\), along the clock (VI.5 (5)) and in energy:

| filling | weight along the clock | energy from the empty cut | from the settled one |
| --- | --- | --- | --- |
| neither | \\(1\\) | \\(0\\) | \\(+\|\\varepsilon\|\\) — a hole |
| leading loop only | \\(e^{-\\omega\\tau}\\) | \\(+\|\\varepsilon\|\\) | \\(+2\|\\varepsilon\|\\) — a leading loop and a hole |
| lagging loop only | \\(e^{+\\omega\\tau}\\) | \\(-\|\\varepsilon\|\\) | \\(0\\) |
| both | \\(1\\) (\\(\\det L_k=1\\)) | \\(0\\) | \\(+\|\\varepsilon\|\\) — a leading loop |

This is (B) once per slot: the long clock leaves the filling in which each slot holds its settled combination. Counted from it every other configuration costs \\(\\sum|\\varepsilon|\\) over the leading loops present and the lagging loops absent, never less than \\(0\\); and a lagging loop missing costs what a leading loop added costs. The missing lagging loop is the hole, the antiparticle: by (A) the lagging loop has energy \\(-\\omega-\\varphi\\) in a field, so its absence costs \\(+\\omega+\\varphi\\), the energy of a leading loop in the field \\(-\\varphi\\) (VI.6 has the propagation). A field is not needed to have the two; it is needed to tell them apart.

*What is settled and what is taken.* A constant added to all \\(256\\) energies cannot be seen on one ring. Two things single out the zero of (7): it is the one for which complements have opposite energies, and it is the empty cut, the one configuration with nothing on it to change when the ring’s length does — which is what the companion document’s force needs. That the vacuum *is* the settled configuration is VI.5’s step, taken there.

*The same number without a snapshot.* In a snapshot — one time, the list of slots, which of them hold a walker — there is one walker for every pair of slots. They are not put in by a separate rule: they are how closed histories look in a snapshot, a history that goes round \\(n\\) times showing as \\(n\\) walkers (VI.5 (A)):

\\\[\\sum\_{\\text{all fillings}}(\\text{weight})=\\det(1+U)=\\exp\\Big\[\\sum\_{n\\ge1}\\frac{(-1)^{n+1}}n\\operatorname{Tr}U^n\\Big\].\\tag{7a}\\\]

7a. on the left, every filling of the list and every number of walkers; on the right, one possible walker followed round a closed history \\(n\\) times, with the exchange sign of the cycle. For one pair of slots, \\((1+e^{a\\tau})(1+e^{-a\\tau})=1+e^{-a\\tau}+e^{+a\\tau}+1\\): the four fillings of the table. The rule sits on each path, and the filled list follows from it.

**(F) Three dimensions: the reversal reverses the axis.** IV.5a’s eight counters, the same \\(\\mathbb K=\\operatorname{diag}(1_4,-1_4)\\). The axis, like the energy, is a sense of turning read against \\(\\mathbb J\\):

\\\[\\Sigma_z=-\\mathbb J\\begin{pmatrix}A_xA_y&0\\\\0&A_xA_y\\end{pmatrix}\\ \\text{and cyclic},\\qquad\\mathbb K\\,\\Sigma_j\\,\\mathbb K=-\\Sigma_j,\\qquad(E,\\ \\mathbf k,\\ \\hat{\\mathbf s})\\ \\xrightarrow{\\ \\mathbb K\\ }\\ (-E,\\ -\\mathbf k,\\ -\\hat{\\mathbf s}).\\tag{8}\\\]

8. so a missing lagging loop with \\((-E,-\\mathbf k,-\\hat{\\mathbf s})\\) is a leading loop with \\((+E,+\\mathbf k,+\\hat{\\mathbf s})\\) and the opposite charge. In IV.4’s chiral notation \\(\\mathbb K\\) reads \\(\\psi\\to C\\,\\overline{\\psi}\\) with

\\\[C=(\\text{a phase})\\times\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}\_{\\text{the two senses}}\\times\\begin{pmatrix}0&1\\\\-1&0\\end{pmatrix}\_{\\text{the two axis components}}:\\tag{9}\\\]

9. conjugate, and apply the turn matrix twice — once between the blocks, once inside each. On the axis pair alone, \\(T\\overline{\\sigma_k}T^{-1}=-\\sigma_k\\): the turn matrix of the walker at rest, \\(v=Tu\\), is the one that reverses the axis.

*What went in.* No new input. (A) is products of III.15e’s matrices, \\(\\mathbb K\\) among them. (B)–(D) are III.2’s lines and III.15d (12) read as the flip processes they are. (E) is VI.1 (D)’s exclusion and VI.5 (1)’s minors; the empty cut as the reference when the ring changes is taken, and VI.5’s identification of the vacuum is used. (F) is IV.4b’s matrices. What the pair made at a vertex is, leading loop against hole, is left for the section on interactions.

(`scripts/loop_reversal_and_the_sea.py`: 1.0–1i for (A); 4c for (3), (5) and the table; 1f, 2d for (4) and the clock weights; 4a, 4b for (D); 2a–2c, 2f for (E); 3a–3e for (F). Exact matrices, minors written out, and two simulations of the flip process.)

## VI.8 The count as a partition function: free energy, and three entropies

*The picture first.* VI.5 and VI.7 counted the fillings of the slots along the walk’s own clock, each with the weight \\(e^{-E\\tau/\\hbar}\\), and let the clock run long. Stop the clock at a finite length instead. Then several fillings share the count, each slot holds a walker with a probability, and that probability is the Fermi function with the clock length standing where \\(\\hbar/k_BT\\) would. This is not a resemblance: the Boltzmann factor of the time reading *is* the walk’s own count over a clock of that length, by III.15e (4). So a temperature is a clock closed on itself after the length \\(\\hbar c/k_BT\\); the settled energy is the free energy of a long clock; and “settling” is an entropy going to zero. Two other entropies sit close by and are different things: the entropy of one walker’s sense, which rises as the sense is forgotten, and the entropy of the settled configuration seen from inside a region, which counts what the region shares with the rest and stops growing at the Compton length.

*Symbols.* The slots are those of a walker read by all four detectors (VI.5 (C)). Slots, fillings, \\(L_k(\\tau)=e^{\\tau M_k}\\), \\(\\omega_k\\) as in VI.5 and VI.7; \\(\\mathbb J\\), \\(\\mathbb G(k)\\), \\(\\mathbb M_k\\), \\(\\mathcal D\\) of III.15e, order \\((u\_+,u\_-,v\_+,v\_-)\\); \\(\\Pi\_-\\) the settled configuration on a ring (VI.4); \\(\\beta=1/k_BT\\); \\(h(p)=-p\\ln p-(1-p)\\ln(1-p)\\), the entropy of a yes/no with probability \\(p\\). In the formulas \\(\\hbar=c=k_B=1\\).

**(A) The count of fillings is a partition function.** One pair of slots at \\(k\\), clock length \\(\\tau\\). VI.7’s four fillings, with their weights, added:

\\\[Z_k(\\tau)=\\underbrace{1}\_{\\text{neither}}+\\underbrace{e^{-\\omega_k\\tau}}\_{\\text{leading only}}+\\underbrace{e^{+\\omega_k\\tau}}\_{\\text{lagging only}}+\\underbrace{1}\_{\\text{both}}=\\det\\big(1+L_k(\\tau)\\big)=\\big(1+e^{\\omega_k\\tau}\\big)\\big(1+e^{-\\omega_k\\tau}\\big).\\tag{1}\\\]

1. VI.5 (1) for one pair; for many pairs \\(Z\\) is the product. The share of the count in which a slot of energy \\(E\\) holds a walker is

\\\[f(E)=\\frac1{e^{E\\tau}+1},\\qquad E=-\\omega_k\\ \\text{(lagging slot)},\\quad E=+\\omega_k\\ \\text{(leading slot)},\\tag{2}\\\]

2. and the slots are filled independently. It is the Fermi function, with \\(\\tau\\) where \\(1/T\\) stands. The \\(+1\\) is exclusion: a slot holds \\(0\\) or \\(1\\). The long clock, \\(f\\to1\\) and \\(f\\to0\\), is VI.5’s settled configuration.

\\\[F(\\tau)=-\\frac1\\tau\\ln Z=\\sum_k\\Big\[-\\omega_k-\\frac2\\tau\\ln\\big(1+e^{-\\omega_k\\tau}\\big)\\Big\]\\ \\xrightarrow{\\ \\tau\\to\\infty\\ }\\ E_0,\\qquad\\langle E\\rangle=-\\frac{\\partial\\ln Z}{\\partial\\tau}=-\\sum_k\\omega_k\\tanh\\frac{\\omega_k\\tau}2 .\\tag{3}\\\]

3. *the settled energy is the free energy of the long clock.* At a finite clock the count has not finished settling, and how far it has is an entropy:

\\\[S\_{\\rm count}=-\\sum\_{\\text{fillings}}p\\ln p=\\ln Z+\\tau\\,\\langle E\\rangle=\\sum\_{\\text{slots}}h\\big(f(E)\\big).\\tag{4}\\\]

4. three expressions, equal (64 fillings built minor by minor). Per pair of slots it is \\(2\\ln2=1.386\\) when the clock is short — the four fillings share the count — and \\(1.364,\\ 1.164,\\ 0.382,\\ 0.001\\) at \\(\\omega\\tau=0.3,\\ 1,\\ 3,\\ 10\\): it falls to \\(0\\) as one filling takes the whole count. *Settling is this entropy going to zero.*

**(B) The clock length is the inverse temperature.** The energy matrix of the time reading on a mode, and the clock’s step matrix, written out:

\\\[\\mathbb J\\,\\mathbb G(k)=\\begin{pmatrix}ck&0&0&-a\\\\0&-ck&a&0\\\\0&a&ck&0\\\\-a&0&0&-ck\\end{pmatrix},\\qquad\\mathbb M_k=\\begin{pmatrix}-ck&a&0&0\\\\a&ck&0&0\\\\0&0&-ck&a\\\\0&0&a&ck\\end{pmatrix},\\qquad\\mathcal D=\\begin{pmatrix}1&0&0&0\\\\0&0&0&-1\\\\0&0&1&0\\\\0&1&0&0\\end{pmatrix}.\\\]

Both are real. \\(\\mathcal D\\) is III.15e (4)’s relabelling: it turns the two backward detectors into each other, \\((u\_-,v\_-)\\to(-v\_-,u\_-)\\), and leaves the two forward ones alone; it commutes with \\(\\mathbb J\\), and (4) there, multiplied by \\(\\mathbb J\\), reads

\\\[\\mathcal D^{-1}\\,\\big(\\mathbb J\\,\\mathbb G(k)\\big)\\,\\mathcal D=-\\,\\mathbb M_k,\\qquad\\text{so}\\qquad e^{-\\beta\\,\\mathbb J\\,\\mathbb G(k)}=\\mathcal D\\;e^{\\beta\\,\\mathbb M_k}\\;\\mathcal D^{-1}.\\tag{5}\\\]

5. on the left, the Boltzmann factor of the mode \\(k\\) at inverse temperature \\(\\beta\\); on the right, the walk’s own count over a clock of length \\(\\tau=\\beta\\), in relabelled counters. The same \\(k\\) is the wavenumber in \\(t\\) and the tilt in the clock (III.15e (D)). So (A) is not an analogy: *a thermal state of the time reading is the count of fillings of a clock glued to its own start after the length \\(\\hbar c/k_BT\\)* — VI.5 (A)’s gluing, with its exchange sign. Nothing is continued: where one usually replaces \\(t\\) by \\(-i\\beta\\), the walk already has the reading in which the same generator grows and shrinks. What (5) does not say is that a walker which has advanced some distance is hot. The temperature belongs to the closed clock.

**(C) Size and temperature exchange.** Put the walk on a ring of length \\(L\\) in the time reading, closed with the exchange sign, at inverse temperature \\(\\beta\\): wavenumbers \\(p_m=(2m+1)\\pi/L\\), energies \\(\\pm\\varepsilon_m\\), \\(\\varepsilon_m^2=a^2+p_m^2\\). Or run the clock over the length \\(L\\), with the tilts \\(k_n=(2n+1)\\pi/\\beta\\), rates \\(\\pm\\omega_n\\), \\(\\omega_n^2=a^2+k_n^2\\). The two counts are one double sum, done in either order:

\\\[\\sum_m\\ln\\Big\[\\big(1+e^{\\beta\\varepsilon_m}\\big)\\big(1+e^{-\\beta\\varepsilon_m}\\big)\\Big\]\\ \\ \\text{and}\\ \\ \\sum_n\\ln\\Big\[\\big(1+e^{L\\omega_n}\\big)\\big(1+e^{-L\\omega_n}\\big)\\Big\]\\ \\ \\text{are both}\\ \\ \\sum\_{m,n}\\ln\\big(a^2+p_m^2+k_n^2\\big)+\\text{const},\\tag{6}\\\]

6. by \\((1+e^x)(1+e^{-x})=4\\prod_n\\big\[1+x^2/((2n+1)\\pi)^2\\big\]\\). The \\(+\\) of exclusion puts the zeros at odd multiples of \\(i\\pi\\): that is why the partner’s numbers are half-integers, the companion document’s \\((-1)^m\\) on the \\(m\\)-th round trip. Each sum in (6) grows without bound, by a term proportional to \\(L\\beta\\) that is the same on both sides; what is equal, and checked, is everything else. With mass: the change of \\(\\ln Z\\) between two masses, from either side, to \\(10^{-15}\\). Without mass:

\\\[\\frac{\\beta\\,\\pi}{6L}+4\\sum\_{m\\ge0}\\ln\\big(1+e^{-\\beta(2m+1)\\pi/L}\\big)\\qquad\\text{is unchanged when }\\beta\\text{ and }L\\text{ are exchanged.}\\tag{7}\\\]

7. the first term is \\(-\\beta E_0(L)\\) with \\(E_0=-\\pi\\hbar c/6L\\), the settled energy of this ring; the sum is what a finite temperature adds. Read from the other side the same number is \\(L\\pi/6\\beta=-\\beta\\,L\\,\\mathsf f\\) with \\(\\mathsf f=-\\pi/6\\beta^2\\): *the settled energy of a short ring is the thermal free energy of a long hot one.* Per unit length the entropy of the count is \\(\\pi/3\\beta\\). (This ring closes with the exchange sign, hence \\(-\\pi/6L\\); the companion document’s \\(-\\pi/24L\\) is between walls. One number, \\(c=1\\), in both.)

**(D) Three entropies, kept apart.**

| entropy of | what it measures | what it does |
| --- | --- | --- |
| the count of fillings, (4) | how many fillings still share the count at clock length \\(\\tau\\); by (5), the thermal entropy at \\(T=1/\\tau\\) | falls from \\(2\\ln2\\) per pair to \\(0\\) |
| one walker’s sense, VI.7 (3) | how far the sense is forgotten | rises from \\(0\\) to \\(\\ln2\\) at rest; its distance from equilibrium only falls |
| the settled configuration inside a region, (E) | what the region shares with the rest | grows with the logarithm of the region, and stops at the Compton length |

The second, in motion: the equilibrium of VI.7 (6) is lopsided, so the sense is not forgotten entirely. Its entropy at equilibrium is \\(h\\big(\\tfrac12(1-v/c)\\big)\\) — \\(0.582,\\ 0.251,\\ 0.090\\) at rapidity \\(0.5,\\ 1.3,\\ 2\\) — equal to \\(\\ln2\\) at rest and \\(0\\) at the speed of light, where the sense is known.

**(E) The settled configuration seen from inside a region.** A configuration of VI.4’s kind is fixed by \\(\\Pi\_-\\). Seen from a region of \\(\\ell\\) sites only, what remains is \\(\\Pi\_-\\) with its rows and columns kept to that region; if its eigenvalues are \\(\\nu_i\\), the region’s entropy is \\(\\sum_ih(\\nu_i)\\) (the standard result for such configurations: I. Peschel, “Calculation of reduced density matrices from correlation functions”, arXiv:cond-mat/0212631, doi:10.1088/0305-4470/36/14/101; cited, not derived here). It is \\(0\\) if the region shares nothing with the rest. Computed:

| the settled configuration of | region | entropy | reading |
| --- | --- | --- | --- |
| the continuum lines, 601 sites, almost no mass | \\(\\ell=4\\) to \\(200\\) | a straight line in \\(\\ln\\ell\\), slope \\(0.666\\) | two edges in \\(k\\), one third each |
| the same with mass, Compton length 20 sites | well inside / well outside it | slope \\(0.645\\) / \\(0.306\\) | the mass closes the edge at \\(k=0\\) |
| the same, \\(\\ell=200\\), mass varied |  | slope \\(0.335\\) in \\(\\ln(1/a)\\) | one third from that edge |
| the lattice tick, 1200 sites, coin angle \\(\\theta\\) | \\(\\ell\\gg1/\\theta\\) | level; slope \\(0.666,\\ 0.663\\) in \\(\\ln(1/\\theta)\\) | two gaps, one third each |

One third of a logarithm for every place where the filling changes abruptly from one \\(k\\) to the next. The physical one is at \\(k=0\\), where the leading and lagging slots meet; the mass rounds it off, and then the region shares with the rest \\(\\tfrac13\\ln(\\text{Compton length})\\) and no more, however large the region is. The other edge in each row belongs to the lattice, not to the walk: for the continuum lines the sharp stop at the largest \\(k\\); for the tick a second gap at quasi-energy \\(\\pi\\), the sublattice of VI.4 (F). So what the region shares is set by its two ends and not by its size, and its coefficient is the same \\(c=1\\) as in (7), in \\(\\pi/3\\beta\\) and in the companion document’s \\(-\\pi/24L\\).

*What went in.* No new input in (A)–(D): (A) is VI.5 (1) at a finite clock; (B) is III.15e (4); (C) is a product formula and the companion document’s exchange of modes and round trips; (D) is VI.7. In (E) the formula for the region’s entropy is the textbook one, cited. Not done: (E) in three dimensions, where the two ends become a surface; and why the one number \\(c\\) appears in all four places, which is used here as an observation.

(`scripts/count_entropy_free_energy.py`, 15 checks: 1a–1d and 2 for (A); 6a for (B); 6b–6d and 3 for (C); 4a, 4b for (D); 5a–5c for (E).)

### What went in (Part VI, so far)

| the coin is the event’s, not the walker’s | **the input of VI.1**: a reading of II.1’s noise, of III.15b’s kind |
| --- | --- |
| one process per point and direction (exclusion) | from that reading (VI.1 (B)) |
| the odd combination; \\(-1\\) under exchange; conservation | derived from exclusion in every state, by polarization, in \\(1{+}1\\) and \\(3{+}1\\) (VI.1 (D), (E)); not from the free step or conservation, which keep both sectors — exclusion in the site–sense basis alone admits a symmetric column (VI.1 (D)) |
| a constraint on paths leaks; an elastic collision is a third statistics | derived, with the numbers (VI.1 (C)) |
| statistics follows the presence of a direction at a point: directed walker fermion, Part II’s walker boson; the link angle is a relation between points and its carrier is Part II’s walker at \\(m=0\\) | a reading (VI.1 (F), VI.2 (E)) |
| the coin does not read the clock hand: a constant turn of the pair passes through the step; a local turn is a potential, \\(\\tilde\\Phi=\\Phi-\\frac\\hbar e\\partial_t\\theta\\), \\(\\tilde A_x=A_x+\\frac\\hbar e\\partial_x\\theta\\), with zero field | derived from III.1’s rule on V.2’s step (VI.2 (A), (B)) |
| an event reads the pair only through the loop angle \\(\\theta_i-\\theta_j\\); flux periodic in \\(2\\pi\\) | derived (VI.2 (C)) |
| the force is the gradient of the turn angle per tick; \\(d\\boldsymbol\\pi/dt=e(\\mathbf E+c\\boldsymbol\\alpha\\times\\mathbf B)\\) exactly | derived (VI.2 (D)) |
| the carrier is Part II’s walker at \\(m=0\\) with a one-form index, \\(\\eta\_{\\mu\\nu}/(4\\pi^2s^2)\\); Coulomb’s law; like walkers repel by \\(\\eta\_{00}=-1\\) | the Casimir document’s sentence read with VI.1 (F); the index from what is carried; then III.16, derived (VI.2 (E)) |
| the link angle is a noise of the event with a local quadratic weight | **the one input of VI.2**, of II.1’s kind; given it, \\(F^2\\) (textbook, cited) and \\(\\partial\_\\mu F^{\\nu\\mu}=j^\\nu\\) follow, with the source conserved (VI.2 (F)) |
| two identical walkers in one packet have the turn \\(T\\) for their axis part, whatever their axes (identity (0)); the entanglement is VI.1’s exchange term \\(-\\sin a\\sin b\\); \\(\\mathcal E=-\\cos(a-b)\\), \\(\|S\|=2\\sqrt2\\) | derived, no new input (VI.3 (A), (B)) |
| two strangers scattered by VI.2’s field leave the pair \\(r\\chi_1\\chi_2^{\\mathsf T}-t\\chi_2\\chi_1^{\\mathsf T}\\); \\(C=2\|\\det M\|/\\\|M\\\|^2\\), \\(\|S\|\_{\\max}=2\\sqrt{1+C^2}\\); the bounce–pass phase picks the sorting planes | derived; regular QM’s numbers (VI.3 (B′)) |
| no violation in \\(1{+}1\\); in \\(3{+}1\\) the relativistic fall-off of \\(S\\) is QM’s for the \\(\\Sigma\\)-reading detector; no signalling | runs (VI.3 (C), (D)) |
| of Bell’s three, the walk lacks outcome-independence: an event has two endpoints (input 9) | a reading (VI.3 (E)) |
| the rule between two walkers is Maxwell solved for one walker’s current and read by the other’s segments; which inverse of \\(\\Box\\) is the carrier’s state | derived from V.2 and (F); the inverse is a choice, named (VI.2 (F′)) |
| Coulomb’s law, the push \\(-e\\langle\\partial_x\\Phi\\rangle\\), and the magnetic sign from \\(\\eta\_{00}\\), run as one process on the lattice | runs (VI.2 (G)) |
| a Z created into an occupied state has weight zero, \\(D-E=0\\) | derived by hand (VI.1a) |
| a configuration of walkers is a determinant; adding and removing a walker and their algebra | derived from exclusion (VI.4 (A)–(B)) |
| the vacuum is the configuration the long-clock count of exclusive configurations selects | **taken** (VI.4 (C), VI.5); that one configuration carries the count, that it is the lagging class in every frame, and its energy \\(-\\tfrac12\\operatorname{Tr}\\lvert H\\rvert\\) are derived; without exclusion there is no reference |
| the sign of a closed history is the exchange sign of the cycle it closes, \\((-1)^{n-1}\\) | derived (VI.5 (A)) |
| Feynman’s rule; pair creation as the turn between configurations; the vacuum’s energy as a sum over closed loops; its polarization, \\(\\alpha a^2\\to1/6\\pi\\); the Z and the closed loop as exchange partners | derived, given the vacuum (VI.4 (C)–(G)) |
| a missing walker is one walker of the opposite orientation and charge, free and with a field on; the \\(\\Pi\_-\\) leg lies in the past at a source | derived (VI.6) |
| a leading loop has two reversals, counting the added detectors with the other sign, and the process run backward, and they pair every energy with its opposite; the lagging loop is the flip process at equilibrium and the leading loop the memory of the sense, in every frame; the mean speed is the imbalance of senses of that equilibrium; \\(E_0=-\\tfrac12\\sum\\lvert\\varepsilon\\rvert\\) is the lowest of the configurations counted from the empty cut, and a hole costs what a leading loop costs | derived (VI.7 (A)–(E)); the empty cut as the reference when the ring’s length changes is **taken** |
| the reversal that reverses the energy reverses the axis; in chiral notation it is conjugation and the turn matrix twice | derived (VI.7 (F)) |
| the count of fillings stopped at a finite clock is a partition function: the Fermi function with the clock’s length where \\(\\hbar/k_BT\\) stands; the settled energy is the free energy of the long clock; settling is an entropy going to zero | derived (VI.8 (A)), from VI.5 (1) |
| a thermal state of the time reading is the count of a clock closed on itself after the length \\(\\hbar c/k_BT\\); size and temperature exchange, \\(-\\pi/6L\\leftrightarrow\\pi/3\\beta\\) | derived (VI.8 (B)–(C)), from III.15e (4) and a product formula |
| three entropies kept apart: of the count, of one walker’s sense, of the settled configuration inside a region; the last is a third of a logarithm per edge and stops at the Compton length | computed (VI.8 (D)–(E)); the formula for a region’s entropy is cited; three dimensions not done; why one number \\(c=1\\) in four places is an observation |
| spin–statistics as a theorem; bosons beyond Part II | **not yet** |

---

# The three routes side by side

|  | Part I: diffusing particle | Part II: walk in spacetime | Part III: time-sense particle (telegraph in remarks) |
| --- | --- | --- | --- |
| clock | \\(t\\) | step label \\(\\tau\\); \\(t\\) is a coordinate | \\(x\\) (the telegraph: \\(t\\)) |
| random element | displacement in \\(x\\), size \\(\\sqrt{\\Delta t}\\) | displacement in \\((t,x)\\), variance \\(\\propto\\eta^{\\mu\\nu}\\): the time kick is imaginary | the time sense \\(\\varepsilon=dt/d\\tau\\) flips at rate \\(a\\); \\(dx/d\\tau=c\\) exactly (telegraph: the space direction flips, \\(t\\) exact) |
| Taylor step | second order survives → \\(D\\partial_x^2\\) | second order survives → \\(D\\Box\\) | first order only → \\(\\pm\\partial_t\\) (telegraph: \\(\\pm c\\partial_x\\)) |
| clock reversal | backward drift; \\(u=D\\partial_x\\ln\\rho\\) | backward drift; \\(u^\\mu=D\\partial^\\mu\\ln\\rho\\) | the two time senses swap; rates \\(aP\_\\mp/P\_\\pm\\); no \\(u\\) |
| kinetic term from the rule | \\(\\frac m2(v^2+u^2)\\) | \\(\\frac m2(v\\cdot v+u\\cdot u)\\) | constant: no \\(\\rho\\)-only term |
| the Hamiltonian | built from the kinetic energy | built from the kinetic term | read off from the equations |
| the \\(i\\) | heat → Schrödinger, \\(D\\to iD\\) | same | already present: the time senses mix by a quarter-turn (telegraph: a half-turn, needs \\(a\\to-ia\\)) |
| mass | — | \\(\\tau\\)-independence, \\(K=-mc^2/2\\) | reversal rate of the time sense, \\(mc^2/\\hbar\\), real |
| output | Schrödinger | Klein–Gordon; Schrödinger as slow limit | Klein–Gordon with the right sign from a real rate (III.6); Dirac; each component a Klein–Gordon fluid of Part II (III.13); Schrödinger, \\(u\\), \\(D=\\hbar/2m\\) and the osmotic term as its slow limit |
| three dimensions | \\(\\partial_x^2\\to\\nabla^2\\) | \\(\\Box\\) with \\(\\nabla^2\\) | motion at \\(c\\) along an axis; the axis must be held as a column; the two time senses are the two handednesses, turned into each other by the mass |
| what a point does | the Gaussian of I.1a | the kernel summed over \\(\\tau\\): the source answer, \\(K_0\\) in \\(1{+}1\\), \\(\\frac{m}{4\\pi^2 s}K_1\\) in \\(3{+}1\\) (III.16, IV.10) | the evolution answer: \\(\\delta\\) on the light line plus \\(I_0,I_1\\) outside the cone (III.14); continued across the cone it is the Dirac amplitude (III.15); in three dimensions no one-clock kernel exists (IV.2a); the process is III.1 axis by axis, in turn, whose limit is Dirac (IV.10) |

---

## What is put in, and where (the one list)

Every Part ends with a ledger; this is the one list. An *input* is a sentence the document asserts and does not derive; a *reading* is a number already on the page called by a physical name; a *choice* selects one of several processes the same inputs allow. Two rows carry the word *taken*: 11 says which detectors read the walker, which the lines do not force; 18 is an identification — the configuration the count singles out is called the vacuum. Everything else in the document is algebra on these, each line checked by a script.

|  | what | kind | where it enters | what it buys |
| --- | --- | --- | --- | --- |
| 1 | the step is a coin: a drift plus an independent kick, scale set by \\(D\\) | input | I.1 | Part I’s kinematics (\\(\\mathbf b\_\\pm\\), \\(\\mathbf u\\), continuity) |
| 2 | the energy is the average of the forward and backward kinetic energies | input (a rule) | I.5 | \\(T=\\frac m2(v^2+u^2)\\), hence the quantum potential; I.5 says what the other two choices give |
| 3 | \\(\\mathbf v=\\nabla S/m\\) (irrotational flow) and \\(\\dot S=-\\delta H/\\delta\\rho\\) | input (the one dynamical law) | I.6–I.7 | Hamilton–Jacobi with the quantum potential; Schrödinger with \\(D=\\hbar/2m\\) |
| 4 | \\(\\psi\\) single-valued | input (Wallstrom) in Part I; from III.8c on a consequence of 11 | I.9 ledger, IV.5a (H) | circulation quantization, stated once: each component pair \\((u_j,v_j)\\) is two single-valued real fields, so its angle winds by \\(2\\pi n\\) around any loop avoiding its zeros; for a uniform axis the four windings are equal and that is Part I’s single \\(S\\); for a varying axis the current’s circulation adds IV.8 (B)’s \\((h/m)k\\sin^2\\frac{\\theta_0}2\\). In Part I, with one density and a real \\(S\\), it is an input; once the walker is read by all its detectors it is derived (IV.5a (H)) |
| 5 | no preferred frame: the kick’s covariance is \\(\\eta^{\\mu\\nu}\\) | input | II.1 | the Lorentzian sign; the time kick imaginary, so Part II is a formal walk, real in Euclidean signature (III.16) |
| 6 | nothing observable depends on the step label, with \\(K=-mc^2/2\\) | input | II.7 | the mass shell; Klein–Gordon |
| 7 | the clock is \\(x\\) and the coin is in \\(t\\) (the time-sense walker S, not the telegraph T) | choice | III.0–III.1 | Klein–Gordon with the right sign from a real rate; the sign the continuation supplies elsewhere is the choice of clock (III.5a (J)). Thesis: McKeon–Ord 1992 |
| 8 | a record is data on a spacelike cut | input, of 5’s kind | III.15b | the \\(t\\)-reading is the physical one; the return leg exists on that cut (III.15a (B), III.15c) |
| 9 | an event is an outward route and a return route that agree on the endpoint; every such pair counts | definition | III.15a (C) | the weight of an event is the product of the two readings, \\(Q^2\\); the degree is then forced by 5 and conservation (III.15c); the assumption of Bell’s the walk does not have (VI.3 (E)) |
| 10 | the weight of an event is its frequency | reading, at the place III.1 read step factors as probabilities | III.15a, V.5 | the Born rule as used; nothing in the document derives it; on two walkers it is the weight (2) of VI.3 |
| 11 | the walker is read by a detector for each time sense and each orientation: eight in three dimensions, four in one (time sense and position of the clock hand) | taken: the lines do not force the added detectors (IV.5a (C)); needed for a walker a field can turn (V.1–V.2), not for every walker that can be detected; the Hamiltonian form of III.9–III.10 is a reading of equations III.8c (C) derives | III.8c (B), IV.5a | the walker at rest (\\(v=Tu\\): the same amount forward and backward); Dirac on complex columns, the complex column being the original and the added detectors with the walker’s own turn \\(J\\) written as \\(i\\) (IV.5a); III.10’s fluid; the sign of energy as the sense in which the walker’s clock turns (III.15e); sums, differences and the frame (IV.5b) |
| 12 | one speed \\(c\\) in every direction, for every component | input | IV.1, IV.4a | the four anticommuting matrices; four components; no positive one-clock process in \\(3{+}1\\) (IV.2a), the automaton instead (IV.10) |
| 13 | the field enters as a phase per step, \\(e^{\\frac{ie}\\hbar(-\\Phi dt+A_xdx)}\\), i.e. a turn of the walker in its own plane, its time sense untouched | input | V.2 | minimal coupling; gauge freedom as where the clock hand is set; \\(g=2\\) |
| 14 | a detector is a walker coupled through the walker’s own turn, turned by \\(g\\) while the particle is on its site | model | V.5 | the detector registers a return, \\(Q^2+Q'^2\\), never anything linear |
| 15 | the coin is drawn at the event, once — two walkers in one state, any state (any column, not only a site and a sense), draw it together | reading of 5, stated for every state | VI.1 (B), (D) | one process per state (exclusion); by polarization, the odd combination and the \\(-1\\) under exchange; the free step keeps the sign. Exclusion in every state is the input; antisymmetry is what follows |
| 16 | the link angle is a noise of the event with a local quadratic weight in two derivatives | input, of 5’s kind | VI.2 (F) | \\(F^2\\) (textbook from here), Maxwell with the walker’s current; the carrier is Part II’s walker at \\(m=0\\) with a one-form index. Rows 15 and 16 are one idea read twice — everything random belongs to the event, not to the walker — applied to the coin (15) and to the link (16) |
| 17 | \\(c\\), \\(\\hbar\\), \\(m\\) (\\(=\\hbar a/c^2\\)), \\(e\\) | constants | — | scales |
| 18 | the vacuum is the configuration that the count of exclusive configurations along the walk’s own clock selects | taken | VI.4 (C), VI.5 | Feynman’s rule; the hole as one walker of the opposite orientation and charge (VI.6); \\(E\_{\\rm ref}=-\\tfrac12\\operatorname{Tr}\\lvert H\\rvert\\) and the companion document’s vacuum; the name “lowest” is the orientation of the record’s cut (VI.4 (F)) |

What is claimed as derived, then, is derived *from these*; where a section says “no new input” it means no line beyond this table. What is not derived anywhere: 10 (why a weight is a frequency), 4 in Part I, the scale of 16, 18 (that the vacuum is that selection), 11 (that the walker is read by all the detectors: the lines do not force it), the carrier’s own state (VI.2 (F′)), spin–statistics as a theorem, bosons beyond Part II’s walker.

---

## Next questions

*What this section is.* The one list says what is put in and what is not derived. This section lists questions the document raises and has not worked out. Each says where it attaches. For none of them is it known yet whether the answer follows from what is already here, or needs one more thing put in.

1. **A flip rate without a coin.** Every stretch of the walk is light-like, and the mass is the rate at which the time sense flips, \\(a=c/\\ell\\) with \\(\\ell\\) the mean length of a stretch. So mass is already a confinement of light-like motion: the time sense cannot run on for longer than \\(1/a\\) on average. For the walker at rest with the same amount forward and backward, a potential and the flip rate are one turn and their rates add (V.4 (2)). But a field’s turn never flips the time sense (V.4 (C)): with \\(a=0\\) the forward and the backward detectors exchange nothing, and a potential cannot make them. The question: with the coin at \\(a=0\\), what supplies the exchange of forward and backward at a nonzero rate? A mass gap would be such a rate. Where it attaches: several sets of detectors give more turns than one (V.1, the door to non-abelian fields); the carrier is Part II’s walker at \\(m=0\\) (VI.2 (E)); a stretch ended by an exchange with the vacuum’s loops and not by the walker’s own coin (VI.4 (G)).
2. **The walker that is not one clock.** IV.5b’s \\(\\gamma\\neq0\\): a walker that has neither sums only nor differences only in any running detectors. Not treated.
3. **The neutral walker.** A walker in the original detectors alone cannot be turned by a field (III.15e (B)). What it is as a particle is not treated.
4. **The particle at a vertex.** What the walker of VI.7 is at a vertex. Not treated.
5. **Gravity.** Two places in the document touch it: a change of speed by sense is a change of the light cone (V.1, third row), and a tilt of the step matrices is a local frame (V.3). Neither is followed. A third route is thermodynamic (T. Jacobson, *Phys. Rev. Lett.* 75, 1260 (1995), arXiv:gr-qc/9504004). If the settled configuration has an entropy across any edge that goes as the edge’s area, and an accelerated detector reads a temperature, then “heat = temperature × change of entropy”, asked at every edge, gives Einstein’s equation. What that needs from the walk: the entropy of a region of the settled configuration, which VI.8 (E) has in one dimension (a third of a logarithm per edge, stopping at the Compton length); the same in three dimensions, where the two ends become a surface, which is not done; and the temperature an accelerated detector reads, for which VI.8 (B) has only the closed clock.
6. **How many slots are filled.** The vacuum’s energy is \\(-\\tfrac12\\operatorname{Tr}|H|\\), a sum over every slot (VI.4 (8)), and that the count selects this configuration is taken (input 18). Summed down to a shortest length, such a sum exceeds the energy density observed in empty space by a factor of about \\(10^{55}\\) (cut at \\(100\\) GeV) to about \\(10^{123}\\) (cut at the Planck length). What fixes the number of settled walkers is not treated. One back-of-envelope estimate, recorded for how it fails: one settled walker per particle of matter, \\(|E|=mc^2\\) each, gives the rest-energy density of matter, \\(2.4\\times10^{-10}\\ \\mathrm{J/m^3}\\) (\\(3.8\\times10^{-11}\\) for baryons alone, \\(0.25\\) per cubic metre), against \\(5.3\\times10^{-10}\\ \\mathrm{J/m^3}\\) observed: within a factor of \\(2\\) to \\(14\\). It is not an answer, for three reasons. It restates that the two densities are comparable today, and they were not earlier: matter was \\(3.7\\) times larger at redshift \\(1\\) and \\(29\\) times at redshift \\(3\\). The settled energy is negative. And matter is the leading loops, not the settled ones (VI.7). (Densities from the parameters of Planck Collaboration, “Planck 2018 results. VI. Cosmological parameters”, arXiv:1807.06209: \\(H_0=67.4\\) km/s/Mpc, \\(\\Omega_m=0.315\\), \\(\\Omega_bh^2=0.0224\\). `scripts/napkin_vacuum_estimate.py`: arithmetic only.)
7. ** Non abelian include reflection not just rotation types of preservation, refleciton -> mass making
