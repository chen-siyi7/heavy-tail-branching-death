# heavy-tail-branching-death

Reproducible R code for the numerical experiments accompanying

> *Fixed Cell Death and Annealed Moment Explosion in Heavy-Tailed
> Branching Models of Secondary Tumors.*

The script reproduces six figures and four LaTeX tables that appear in the
Numerical Verification section of the manuscript and its Appendix C.

## What this code does

The paper studies a continuous-time, infinite-type branching process in
which clone-specific birth rates `a` are drawn from `Exp(λ)` and a fixed
per-cell death rate `d` is added under the *independent-death convention*.
The script verifies the paper's main analytical claims by direct
simulation:

| Experiment | What it shows | Output |
|---|---|---|
| 1 | Annealed mean `E[X_1(t)]` matches `e^{-dt} λ/(λ - (1-μ)t)` and diverges at `t* = λ/(1-μ)` for every fixed `d` | `fig1_annealed_mean.pdf`, `table2_annealed.tex` |
| 2 | Eventual survival probability matches `exp(-λc) - λc E_1(λc)` with `c = d/(1-μ)`, strictly below the naive bound `exp(-λc)` | `fig2_survival.pdf`, `table1_survival.tex` |
| 3 | Bounded birth-rate support `a ∈ [0, a_max]` eliminates the finite-time moment explosion; `d* = (1-μ) a_max` separates long-time growth from decay | `fig3_bounded.pdf`, `table3_bounded.tex` |
| 4 | Random independent death rate damps the annealed mean by the Laplace transform `E[exp(-dt)]` but leaves the explosion time unchanged | `fig4_random_d.pdf` |
| 5 | Empirical means systematically underestimate the analytic mean near `t*`; the bias scales as `n^{-(ν-1)/ν}` with `ν → 1` | `fig5_bias_scaling.pdf` |
| 6 | Leading small-`d` expansion `1 - λc(1 - γ - log(λc))` matches the exact survival probability up to `O((λc)^2)` | `fig6_smalld.pdf`, `table4_smalld.tex` |

## Requirements

R ≥ 4.0 (tested on 4.4.1) with the following packages (auto-installed
on first run if missing):

```
ggplot2, scales, expint, dplyr, tidyr, xtable, gridExtra
```

## Running

From the command line:

```sh
Rscript simulations.R
```

This runs all six experiments and writes all outputs to the working
directory. Wall time is roughly 5–15 minutes on a recent laptop; the
dominant costs are Experiments 2 and 5 (20,000 and up to 250,000 total
replicates, respectively).

To run interactively, e.g. to inspect intermediate results:

```r
source("simulations.R")
res <- run_all("results")     # any output directory
res$e2$data                   # data frame for Experiment 2
```

## Reproducibility notes

- The seed `set.seed(20260509)` at the top of the script reproduces
  the empirical values reported in the manuscript tables (within R 4.4
  and the listed package versions).
- Figures use ASCII-safe labels for portability across PDF viewers and
  embed TrueType fonts via the `pdf()` device defaults.
- `Experiments 1, 3, 4, 5` cap each trajectory at `max_pop = 2e5`
  (`1e6` for Experiment 5) to bound memory; the script reports the
  fraction of capped trajectories per parameter setting so the
  user can see whether the heavy right tail is being truncated.
- `Experiment 2` uses a Rao-Blackwellized estimator for eventual
  survival: it walks the embedded discrete-time chain until extinction
  or until the population reaches `N_escape = 200`, then returns the
  exact conditional survival probability `1 - q(a)^{N_escape}` rather
  than treating every escape as certain survival. This removes the
  small upward bias that a binary estimator would introduce near
  criticality.

## Files

```
simulations.R    # the script (single file, all experiments)
README.md
LICENSE
.gitignore
```

## License

MIT. See `LICENSE`.

## Citation

If you use this code, please cite the manuscript. A `CITATION.cff` will
be added once the paper is published.
