# Cell death in a heavy-tailed branching process

Reproducibility code accompanying Siyi Chen's manuscript, *Fixed Cell Death and Annealed Moment Explosion in a Heavy-Tailed Branching Model of Secondary Tumors*.

The code reproduces the numerical experiments for the ancestral type-1 lineage under independent cell death. It generates all six figures, the data underlying four tables, LaTeX table fragments, and numerical checks. It does not simulate the full infinite-type tumor population or analyze empirical patient data.

## Requirements

The reference computation used R 4.4.1. Only packages supplied with the standard R distribution are used. No CRAN packages, input datasets, network connection, LaTeX installation, or GitHub token are needed to run the code.

Use R 4.4.1 when comparing closely with the reference results. Other R versions and operating systems have not been validated for this release. Floating-point calculations and PDF metadata can differ across platforms.

See [VALIDATION.md](VALIDATION.md) for the completed local reproduction check and source checksum.

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

## Code and reproducibility

`simulation_audit.R` is the unchanged manuscript audit. `run.R` provides path handling, overwrite protection, output checks, and session recording. The master seed is 20260902. Separate deterministic seeds are assigned to each experimental setting and repeated trial. The runner explicitly uses Mersenne-Twister, Inversion, and Rejection, matching the reference R session.

The main parameters are `lambda = 1` and `mu = 0.3`. The moment boundary is `t_star = lambda / (1 - mu)`. The script is a reproduction workflow for these settings, not a general-purpose R package. Changing parameters requires reviewing the time grids, reference constants, captions, and numerical assertions together.

The script first checks transition-law moments and the mean, survival, bounded-support, random-death, tail, and small-death calculations. It then runs the simulations. Figures 1 and 2 use 20,000 replicates per setting. The bounded-support experiment uses 3,000, and the random-death experiment uses 4,000. The finite-sample experiment uses 200 independent trials at each of eight sample sizes.

## Generated figures

The figure filenames retain the six-figure workflow identifiers. The final manuscript places the last three figures in the supplement.

| File | Content | Placement |
| --- | --- | --- |
| `fig1_annealed_mean.pdf` | Annealed mean under fixed death | Main Figure 1 |
| `fig2_survival.pdf` | Eventual type-1 survival | Main Figure 2 |
| `fig3_bounded.pdf` | Bounded birth-rate support | Main Figure 3 |
| `fig4_finite_sample_instability.pdf` | Distribution of sample means near the boundary | Supplement |
| `fig5_random_death.pdf` | Independent random death rates | Supplement |
| `fig6_smalld.pdf` | Leading small-death approximation | Supplement |

The output folder also contains 15 CSV files, eight LaTeX table fragments, `audit_summary.txt`, `session_info.txt`, `run_info.txt`, and `checksums.md5`. The LaTeX fragments are not standalone documents. Their colored variants require the manuscript's revision macros; they are not full manuscript change tracking.

## Interpretation and limitations

Near the moment boundary, clone sizes are heavy-tailed. A typical sample mean can substantially underestimate the analytical mean even though the estimator is unbiased. Large empirical deviations near that boundary should not be interpreted as a change in the theoretical explosion time.

The survival experiment samples Bernoulli indicators using the exact conditional survival formula. It checks averaging and implementation, not an independent proof of that formula. A separate embedded birth-death walk checks the probability of reaching population size 200 before extinction. That finite threshold is an approximation to eventual survival.

The tail check compares a finite count with an asymptotic limit. It is not expected to agree to machine precision. The leading small-death approximation is used only near zero and becomes inaccurate at larger death rates. The transition sampler contains a floating-point safeguard for very small geometric success probabilities; it is not an arbitrary-precision sampler for extreme parameter regimes.

The original PDF export uses Helvetica and preserves the manuscript-era graphical settings. It does not embed all figure fonts. Embed and check the fonts before using these PDFs for journal production. The code package does not change the existing manuscript figures or resolve that production-format issue.

## Citation

`CITATION.cff` identifies this software and its manuscript author. No repository URL, publication DOI, or release DOI has been invented. Add those identifiers when available. GitHub can display citation information from a root-level citation file. [GitHub citation documentation](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files).

## License

No software license has been selected for this package. The author should choose the intended reuse terms before presenting the repository as open source. A public repository is not by itself an open-source license. [GitHub licensing guidance](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/licensing-a-repository).

## Upload to GitHub

Upload the contents of this folder as the repository root, including `.gitignore`, `.gitattributes`, and `CITATION.cff`. Keep generated outputs separate unless you intentionally decide to archive them. The supplied package contains no reviewer reports, response letters, submitted manuscript PDFs, or legacy simulation files. Nothing has been uploaded automatically.
