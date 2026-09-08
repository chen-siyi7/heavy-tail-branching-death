# Cell death in a heavy-tailed branching process

Reproducibility code accompanying Siyi Chen's manuscript, *Fixed Cell Death and Annealed Moment Explosion in a Heavy-Tailed Branching Model of Secondary Tumors*.

The code reproduces the numerical experiments for the ancestral type-1 lineage under independent cell death. It generates all six figures, the data underlying four tables, LaTeX table fragments, and numerical checks. It does not simulate the full infinite-type tumor population or analyze empirical patient data.

## Run

From the repository directory, run:

```sh
Rscript --vanilla run.R
```

The output directory is `results/simulation_audit_output/`. The runner records the R session, random-number generator settings, elapsed time, and file checksums. It refuses to overwrite an existing audit directory. For another run, choose a fresh output parent:

```sh
Rscript --vanilla run.R --output-dir results/run-02
```

An explicit output path is relative to the directory from which the command is invoked. The default `results/` directory is relative to `run.R`. Paths containing spaces should be enclosed in quotes. Help is available with `Rscript --vanilla run.R --help`.

The complete original workflow can also be run directly:

```sh
Rscript --vanilla simulation_audit.R
```

Direct execution writes to `simulation_audit_output/` in the current working directory and can overwrite earlier outputs there. The original script retains an optional `--publish` mode that copies generated figures and table fragments into the current directory with replacement. The recommended `run.R` entry point does not accept that option and does not copy files into a manuscript folder.

## Interpretation and limitations

Near the moment boundary, clone sizes are heavy-tailed. A typical sample mean can substantially underestimate the analytical mean even though the estimator is unbiased. Large empirical deviations near that boundary should not be interpreted as a change in the theoretical explosion time.

The survival experiment samples Bernoulli indicators using the exact conditional survival formula. It checks averaging and implementation, not an independent proof of that formula. A separate embedded birth-death walk checks the probability of reaching population size 200 before extinction. That finite threshold is an approximation to eventual survival.

The tail check compares a finite count with an asymptotic limit. It is not expected to agree to machine precision. The leading small-death approximation is used only near zero and becomes inaccurate at larger death rates. The transition sampler contains a floating-point safeguard for very small geometric success probabilities; it is not an arbitrary-precision sampler for extreme parameter regimes.
