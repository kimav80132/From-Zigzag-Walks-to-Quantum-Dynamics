# From Zigzag Walks to Quantum Dynamics

A derivation of quantum mechanics, and of the first tools of quantum field theory, from random walks
that are local in spacetime and run on real numbers. Work in progress.

## The idea in one paragraph

A walker moves at the speed *c* along a line that serves as its clock. What flips at random is whether
its time runs forward or backward. Its two counters are ordinary probabilities. Read with time as the
clock, the flip becomes a quarter-turn of the pair of counters, and that quarter-turn is the *i* of the
Dirac equation: nothing imaginary is put in by hand. Mass is the rate of the flips. The same construction
is then carried to three space dimensions on eight real counters, to the electromagnetic field, to two
walkers and Bell's inequality, and to the vacuum.

## Where to start

| If you have | Read |
|---|---|
| ten minutes | `notes/walk_to_dirac_propagator_1plus1.pdf`: seven pages, one space dimension, from the walk to the Dirac propagator. Every step is on real counters, with the usual formula quoted beside it |
| ten minutes, in the usual notation | `notes/walk_to_dirac_propagator_1plus1_usual_symbols.pdf`: six pages, the same note with its last part written with *i*, as physicists usually read it |
| an afternoon | `Diffusion and Telegraph(13).md`: the main document, in six parts |
| an interest in sums | `Go-and-come-back sums Euler–Maclaurin, zeta and the Casimir energy(13).md`: the companion document |

## What the main document covers

- **Part I.** The Schrödinger equation from a diffusing particle (Nelson's route), every step written out.
- **Part II.** The Klein–Gordon equation from the same walk taken in spacetime.
- **Part III.** The Dirac equation in one space dimension from a walk whose time sense flips at a real rate;
  its kernels and propagators as sums over zigzag paths; what a measurement is in the walk.
- **Part IV.** Three space dimensions: the Dirac equation on eight real counters, spin, and the propagators,
  including Feynman's rule.
- **Part V.** The electromagnetic field as a turn of the counters.
- **Part VI.** Two and many walkers: exclusion, Bell's inequality (the value 2√2, and which of Bell's
  assumptions the walk does not satisfy), the vacuum and its energy, vacuum polarization in one space
  dimension.

The companion document gives the Casimir energy as a sum over round trips.

## Where it stands

- It is a reformulation. It rederives known equations and does not yet predict a new number.
- It has not been peer reviewed.
- It was written with publicly available AI models, in their consumer versions, with no special access.
  Results that cite a script are checked by that script.
- The inputs it assumes are listed in one table near the end of the main document
  ("What is put in, and where"), followed by the open questions.
- Not done: the interacting theory in three space dimensions. Vacuum polarization is worked out in one
  space dimension only. References are still being completed.

## The checks

`scripts/` holds the 121 scripts that the two documents cite by name. Each one checks the statements
made at the place where it is cited.

```
python3 run_all_checks.py          # one at a time
python3 run_all_checks.py -j 4     # four at a time
```

It needs Python 3 with numpy, scipy, sympy, mpmath and matplotlib (tested with Python 3.13, numpy 2.5,
scipy 1.18, sympy 1.14, mpmath 1.3, matplotlib 3.11). Each script runs in its own process. A script is clean when it exits with code 0 and prints no line starting with `FAIL`. A
report is written to `CHECK_REPORT.md`.

The last full run, on 25 November 2026, from a fresh copy of this repository: **121 scripts, 121 clean,
805 PASS lines, 0 FAIL lines.** Eight of the clean scripts print numbers or figures without PASS lines.
The whole run takes about ten minutes. Its report is in `CHECK_REPORT.md`.

What a check is, and is not:

- Algebraic identities are checked symbolically or to rounding error.
- Comparisons of a lattice walk with a continuum formula are checked to a tolerance stated in the script,
  usually between one part in a hundred and one part in ten thousand. They are not exact.
- A few scripts print numbers or draw a figure without a pass line. `run_all_checks.py` counts them
  separately.
- A passing script shows that the computation it contains is right. Whether that computation is the one
  the text claims is a separate question, and it is the one an audit should ask.

Four scripts in `notes/` check the note: `walk_to_propagator_page_checks.py` (its formulas, in both versions: 25 checks,
the steps on four real counters each against the usual complex form), `note_figures_three_walks_and_lattice_sum.py`
(what its two figures show; with `--write` it redraws them), `walk_lines_against_ord_mann_eq6.py` (its equations against
those of Ord and Mann, 2003) and `real_corners_against_feynmans_i.py` (its real corner weights against Feynman's imaginary one).

## For an automated audit

This repository is meant to be gone through by a model that can run code. A suggested instruction:

> Read the main document. For each statement that cites a script, run the script and judge whether what
> the script checks is what the text claims. Report: (1) scripts that fail; (2) scripts that pass but do
> not test the claim made in the text; (3) steps in the text that have no check behind them; (4) results
> that are already known under another name, with the reference.

Errors found are welcome. Please open an issue.

## License

See `LICENSE`: the text is under CC BY 4.0 and the code under the MIT terms given there.

© 2026 kimav80132
