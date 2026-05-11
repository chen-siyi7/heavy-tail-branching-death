## =========================================================================
##  simulations.R
##
##  Reproducible simulations for the paper
##    "Fixed Cell Death and Annealed Moment Explosion in Heavy-Tailed
##     Branching Models of Secondary Tumors."
## -------------------------------------------------------------------------
## Dependencies. Auto-installs anything missing on first run.
## -------------------------------------------------------------------------
.required <- c("ggplot2", "scales", "expint", "dplyr",
               "tidyr", "xtable", "gridExtra")
.missing <- .required[!vapply(.required, requireNamespace,
                              logical(1), quietly = TRUE)]
if (length(.missing) > 0L) {
  message("Installing missing packages: ",
          paste(.missing, collapse = ", "))
  install.packages(.missing, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(expint)     # provides expint::expint_E1(x) = E_1(x)
  library(dplyr)
  library(tidyr)
  library(xtable)
})

## Reproducibility: this seed produces the empirical values reported in
## the manuscript tables (within R 4.4 and the listed package versions).
set.seed(20260509)

## -------------------------------------------------------------------------
## Plot theme and palette: matches the polished manuscript figures.
## -------------------------------------------------------------------------
PAL <- list(
  blue   = "#2b6cb0",
  orange = "#dd6b20",
  green  = "#2f855a",
  red    = "#c53030",
  gray   = "#718096",
  ink    = "#1a202c"
)
SERIES <- c(PAL$blue, PAL$orange, PAL$green, PAL$red)

theme_paper <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "serif") +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(colour = "#cbd5e0",
                                        linewidth = 0.3),
      axis.line          = element_line(colour = "#4a5568",
                                        linewidth = 0.4),
      axis.ticks         = element_line(colour = "#4a5568",
                                        linewidth = 0.3),
      axis.text          = element_text(colour = "#4a5568"),
      legend.key         = element_blank(),
      legend.background  = element_blank(),
      legend.position    = "right",
      legend.title       = element_text(size = base_size - 1),
      legend.text        = element_text(size = base_size - 1),
      plot.title         = element_text(face = "plain",
                                        size = base_size,
                                        hjust = 0)
    )
}

## =========================================================================
## CORE SIMULATOR
##
## The type-1 clone under the independent-death convention is a linear
## birth-death process with effective birth rate b = (1 - mu) * a and
## death rate d, started from a single cell. The next event time is
## Exp((b + d) X), and the event is a birth with probability b/(b + d),
## a death with probability d/(b + d).
## =========================================================================

## -------------------------------------------------------------------------
## simulate_clone_grid()
##
## Exact Gillespie simulation of a single linear birth-death trajectory.
## Records the population size at each time in t_grid.
##
## Args:
##   a, mu, d  : birth rate, mutation probability, death rate
##   t_grid    : sorted vector of times in (0, t_max] at which to record X
##   t_max     : truncation time
##   max_pop   : memory cap; trajectory truncated when X >= max_pop
##   X0        : initial population (default 1)
##
## Returns a list:
##   x       : numeric vector of length(t_grid), the population at each grid time
##   capped  : logical, TRUE if the cap was hit (right tail truncated)
## -------------------------------------------------------------------------
simulate_clone_grid <- function(a, mu, d, t_grid,
                                t_max, max_pop = 1e6, X0 = 1L) {
  b      <- (1 - mu) * a
  X      <- X0
  t      <- 0
  n_grid <- length(t_grid)
  out    <- numeric(n_grid)
  idx    <- 1L
  capped <- FALSE

  while (idx <= n_grid) {
    if (X == 0) {
      out[idx:n_grid] <- 0
      return(list(x = out, capped = capped))
    }
    if (X >= max_pop) {
      out[idx:n_grid] <- X
      capped <- TRUE
      return(list(x = out, capped = capped))
    }
    rate_total <- (b + d) * X
    if (rate_total <= 0) {
      out[idx:n_grid] <- X
      return(list(x = out, capped = capped))
    }
    dt     <- rexp(1, rate = rate_total)
    t_next <- t + dt
    while (idx <= n_grid && t_grid[idx] <= t_next) {
      out[idx] <- X
      idx <- idx + 1L
    }
    if (idx > n_grid) return(list(x = out, capped = capped))
    if (t_next > t_max) {
      out[idx:n_grid] <- X
      return(list(x = out, capped = capped))
    }
    if (runif(1) < b / (b + d)) X <- X + 1L else X <- X - 1L
    t <- t_next
  }
  list(x = out, capped = capped)
}

## -------------------------------------------------------------------------
## annealed_mean(), annealed_mean_truncated()
##
## Sample a ~ Exp(lambda) (or truncated exponential) once per replicate,
## simulate the corresponding clone, and average X_1 at each grid time.
##
## d may be a scalar, or `d_sampler` may supply a per-replicate sample
## (used in Experiment 4 for random independent death).
##
## Returns a list:
##   mean             : numeric vector, empirical annealed mean at each grid time
##   capped_fraction  : fraction of replicates whose trajectory hit max_pop;
##                      a non-negligible value near t_* indicates that the
##                      empirical mean is biased downward by truncation of
##                      the heavy right tail.
## -------------------------------------------------------------------------
annealed_mean <- function(lam, mu, d, t_grid, t_max,
                          n_reps, max_pop = 1e6, d_sampler = NULL) {
  sums     <- numeric(length(t_grid))
  n_capped <- 0L
  for (k in seq_len(n_reps)) {
    a     <- rexp(1, rate = lam)
    d_use <- if (is.null(d_sampler)) d else d_sampler()
    sim   <- simulate_clone_grid(a, mu, d_use, t_grid, t_max, max_pop)
    sums  <- sums + sim$x
    if (sim$capped) n_capped <- n_capped + 1L
  }
  list(mean = sums / n_reps,
       capped_fraction = n_capped / n_reps)
}

## Inverse-CDF sampler for a ~ Exp(lambda) truncated to [0, a_max].
sample_truncated_exp <- function(lam, a_max, n = 1L) {
  u <- runif(n)
  -log(1 - u * (1 - exp(-lam * a_max))) / lam
}

annealed_mean_truncated <- function(lam, mu, d, a_max, t_grid, t_max,
                                    n_reps, max_pop = 1e6) {
  sums     <- numeric(length(t_grid))
  n_capped <- 0L
  for (k in seq_len(n_reps)) {
    a    <- sample_truncated_exp(lam, a_max, 1L)
    sim  <- simulate_clone_grid(a, mu, d, t_grid, t_max, max_pop)
    sums <- sums + sim$x
    if (sim$capped) n_capped <- n_capped + 1L
  }
  list(mean = sums / n_reps,
       capped_fraction = n_capped / n_reps)
}

## -------------------------------------------------------------------------
## eventual_survival_weight()
##
## Rao-Blackwellized estimator of eventual survival for a single clone
## with birth rate a. Walks the embedded discrete-time chain (Bernoulli
## random walk on Z_{>=0} with p_birth = b/(b+d)) until extinction or
## until X >= N_escape. On hitting the escape threshold, returns the
## *exact* conditional eventual-survival probability of a linear
## birth-death process started from X cells:
##
##   1 - q(a)^X,    q(a) = d/((1-mu) a)   if (1-mu) a > d
##   0                                    if (1-mu) a <= d
##
## By the tower property this is unbiased for Pr(eventual survival).
## A naive estimator that simply returns 1 on hitting N_escape is biased
## upward by E[1{escape} q(a)^N_escape], which is small for N_escape = 200
## but nonzero near criticality.
## -------------------------------------------------------------------------
eventual_survival_weight <- function(a, mu, d, N_escape = 200,
                                     max_steps = 1e6) {
  b          <- (1 - mu) * a
  rate_total <- b + d
  if (rate_total <= 0) return(0)
  p_birth <- b / rate_total
  X       <- 1L

  for (i in seq_len(max_steps)) {
    if (X == 0) return(0)
    if (X >= N_escape) {
      if (b <= d) return(0)               # subcritical: extinction certain
      return(1 - (d / b) ^ X)
    }
    if (runif(1) < p_birth) X <- X + 1L else X <- X - 1L
  }
  ## Fallthrough: max_steps reached without extinction or escape.
  if (X == 0 || b <= d) return(0)
  1 - (d / b) ^ X
}

empirical_eventual_survival <- function(lam, mu, d, n_reps,
                                        N_escape = 200) {
  vals <- numeric(n_reps)
  for (k in seq_len(n_reps)) {
    a <- rexp(1, rate = lam)
    vals[k] <- eventual_survival_weight(a, mu, d, N_escape = N_escape)
  }
  list(emp = mean(vals),
       se  = sd(vals) / sqrt(n_reps))
}

## =========================================================================
## ANALYTIC FORMULAS
## =========================================================================

## E[X_1(t)] = e^{-d t} * lam / (lam - (1-mu) t),  for t < lam/(1-mu).
## Manuscript Proposition 3.1.
analytic_mean_fixed_d <- function(lam, mu, d, t) {
  denom <- lam - (1 - mu) * t
  ifelse(denom > 0, exp(-d * t) * lam / pmax(denom, 1e-300), Inf)
}

## E[X_1(t)] under random independent d (Laplace-transform damping):
##   E[X_1(t)] = lam / (lam - (1-mu) t) * L_G(t),
## where L_G(t) = E[exp(-d t)] is the Laplace transform of the
## death-rate distribution G.
analytic_mean_random_d <- function(lam, mu, t, L_G) {
  denom <- lam - (1 - mu) * t
  ifelse(denom > 0,
         (lam / pmax(denom, 1e-300)) * L_G(t),
         Inf)
}

## Eventual survival probability of a clone:
##   Pr(survival) = exp(-lam c) - lam c E_1(lam c),  c = d/(1-mu).
## Manuscript Theorem 4.3.
analytic_survival <- function(lam, mu, d) {
  if (d <= 0) return(1)
  cc <- d / (1 - mu)
  exp(-lam * cc) - lam * cc * expint::expint_E1(lam * cc)
}

## Naive bound: Pr(supercritical) = Pr((1-mu) a > d) = exp(-lam d/(1-mu)).
## Discussed in Section 4.3 as the looser upper bound.
analytic_naive_bound <- function(lam, mu, d) {
  cc <- d / (1 - mu)
  exp(-lam * cc)
}

## Truncated-exponential annealed mean (closed form):
##   E[X_1(t)] = e^{-d t} * lam / (1 - e^{-lam a_max})
##              * integral_0^{a_max} e^{((1-mu) t - lam) a} da.
## Logical-indexed assignment is used instead of ifelse() to avoid
## evaluating (e^{s a_max} - 1) / s at s = 0 (which gives NaN even
## though the s = 0 branch is selected).
truncated_exp_mean <- function(lam, mu, d, a_max, t) {
  Z <- 1 - exp(-lam * a_max)
  s <- (1 - mu) * t - lam

  integ             <- numeric(length(s))
  near_zero         <- abs(s) < 1e-12
  integ[near_zero]  <- a_max
  integ[!near_zero] <- (exp(s[!near_zero] * a_max) - 1) / s[!near_zero]

  exp(-d * t) * lam / Z * integ
}

## =========================================================================
## EXPERIMENTS
##
## Each experiment_k() returns a named list with simulation parameters,
## empirical/analytic data, and (where applicable) the capped fraction.
## The corresponding plot_exp_k() and table_exp_k() consume that list
## to produce a single PDF or .tex file.
## =========================================================================

## -------------------------------------------------------------------------
## Experiment 1: annealed mean at fixed death rate (Figure 1, Table 2).
## -------------------------------------------------------------------------
experiment_1 <- function(lam = 1, mu = 0.3, n_reps = 4000,
                         max_pop = 2e5) {
  t_star <- lam / (1 - mu)
  t_grid <- seq(0.05, 0.95 * t_star, length.out = 20)
  d_vals <- c(0, 0.5, 1, 2)

  results    <- list()
  cap_report <- setNames(numeric(length(d_vals)), as.character(d_vals))

  for (d in d_vals) {
    cat(sprintf("  exp1: d = %.2f\n", d))
    sim <- annealed_mean(lam, mu, d, t_grid,
                         t_max = 0.95 * t_star,
                         n_reps = n_reps, max_pop = max_pop)
    cap_report[as.character(d)] <- sim$capped_fraction
    ana <- analytic_mean_fixed_d(lam, mu, d, t_grid)
    results[[as.character(d)]] <- data.frame(
      t = t_grid, empirical = sim$mean, analytic = ana, d = d
    )
  }
  if (any(cap_report > 0)) {
    cat("  exp1 capped fraction (max_pop reached):\n")
    print(round(cap_report, 4))
  }
  list(t_star = t_star, n_reps = n_reps, d_vals = d_vals,
       lam = lam, mu = mu,
       capped_fraction = cap_report,
       data = do.call(rbind, results))
}

plot_exp1 <- function(res) {
  df <- res$data
  base <- ggplot() +
    geom_line(data = df,
              aes(t, analytic, colour = factor(d), group = d),
              linewidth = 0.7) +
    geom_point(data = df,
               aes(t, empirical, colour = factor(d)),
               shape = 21, fill = "white", size = 2) +
    geom_vline(xintercept = res$t_star, linetype = "dashed",
               colour = PAL$gray) +
    scale_colour_manual(values = SERIES,
                        name   = "death rate",
                        labels = paste0("d = ", res$d_vals)) +
    labs(x = "t", y = expression(E*"["*X[1](t)*"]")) +
    theme_paper()

  list(linear = base + ggtitle("(a) Linear scale"),
       log    = base + scale_y_log10() + ggtitle("(b) Semilog scale"))
}

table_exp1 <- function(res, file_out = "table2_annealed.tex") {
  df       <- res$data
  df$ratio <- df$empirical / df$analytic
  ## Pick four time-indices per d (early, mid, late, near-explosion).
  picks <- df %>%
    group_by(d) %>%
    arrange(t) %>%
    slice(c(1, 8, 14, 19)) %>%
    ungroup() %>%
    select(t, d, empirical, analytic, ratio)

  cap <- sprintf(
    paste0(
      "Annealed mean $\\mathbb{E}[X_1(t)]$ at selected times for ",
      "four fixed death rates $d$ ($\\lambda = 1$, $\\mu = 0.3$, ",
      "$n_{\\mathrm{rep}} = %d$). Empirical averages match analytic ",
      "values closely well below $t_*$ but systematically underestimate ",
      "near $t_*$, where the distribution becomes heavy-tailed."),
    res$n_reps
  )
  xt <- xtable(picks,
               digits  = c(0, 3, 2, 3, 3, 3),
               caption = cap, label = "tab:annealed",
               align   = c("l", "r", "r", "r", "r", "r"))
  print(xt, file = file_out, include.rownames = FALSE,
        booktabs = TRUE, sanitize.text.function = identity,
        caption.placement = "top")
}

## -------------------------------------------------------------------------
## Experiment 2: eventual survival vs. naive bound (Figure 2, Table 1).
## -------------------------------------------------------------------------
experiment_2 <- function(lam = 1, mu = 0.3, n_reps = 20000,
                         N_escape = 200) {
  d_vals <- seq(0, 2, length.out = 9)
  rows   <- vector("list", length(d_vals))
  for (i in seq_along(d_vals)) {
    d   <- d_vals[i]
    cat(sprintf("  exp2: d = %.2f\n", d))
    est <- empirical_eventual_survival(lam, mu, d, n_reps,
                                       N_escape = N_escape)
    ana <- analytic_survival(lam, mu, d)
    nv  <- analytic_naive_bound(lam, mu, d)
    z   <- (est$emp - ana) / max(est$se, 1e-12)
    rows[[i]] <- data.frame(d = d, empirical = est$emp, se = est$se,
                            analytic = ana, naive = nv, z = z)
  }
  list(lam = lam, mu = mu, n_reps = n_reps, N_escape = N_escape,
       data = do.call(rbind, rows))
}

plot_exp2 <- function(res) {
  df          <- res$data
  d_dense     <- seq(0.001, 2, length.out = 200)
  c_dense     <- d_dense / (1 - res$mu)
  ana_dense   <- exp(-res$lam * c_dense) -
                 res$lam * c_dense * expint::expint_E1(res$lam * c_dense)
  naive_dense <- exp(-res$lam * c_dense)

  curves <- data.frame(
    d = c(d_dense, d_dense),
    survival = c(ana_dense, naive_dense),
    type = factor(rep(c("analytic", "naive"), each = length(d_dense)),
                  levels = c("analytic", "naive"))
  )

  ggplot() +
    geom_line(data = curves,
              aes(d, survival, colour = type, linetype = type),
              linewidth = 0.7) +
    geom_errorbar(data = df,
                  aes(d, ymin = empirical - 1.96 * se,
                          ymax = empirical + 1.96 * se),
                  width = 0.04, colour = PAL$blue) +
    geom_point(data = df, aes(d, empirical),
               shape = 21, fill = "white", size = 2.2,
               colour = PAL$blue) +
    scale_colour_manual(
      values = c(analytic = PAL$blue, naive = PAL$gray),
      labels = c(expression(e^{-lambda*c} - lambda*c*E[1](lambda*c)),
                 expression("naive bound " * e^{-lambda*c}))
    ) +
    scale_linetype_manual(
      values = c(analytic = "solid", naive = "solid"),
      guide  = "none"
    ) +
    scale_y_log10() +
    labs(x = "d", y = "Pr(eventual survival)", colour = NULL) +
    theme_paper()
}

table_exp2 <- function(res, file_out = "table1_survival.tex") {
  df     <- res$data
  df_out <- data.frame(
    d         = sprintf("%.2f", df$d),
    Empirical = sprintf("$%.4f\\;(\\pm %.4f)$",
                        df$empirical, 1.96 * df$se),
    Analytic  = sprintf("$%.4f$", df$analytic),
    Naive     = sprintf("$%.4f$", df$naive),
    z         = sprintf("$%+.2f$", df$z),
    stringsAsFactors = FALSE
  )
  cap <- sprintf(
    paste0(
      "Eventual survival probability of the type-1 clone lineage as a ",
      "function of $d$ ($\\lambda = 1$, $\\mu = 0.3$, ",
      "$n_{\\mathrm{rep}} = %d$, $N_{\\mathrm{escape}} = %d$). ",
      "Each replicate uses the Rao-Blackwellized estimator that ",
      "returns $1 - q(a)^{N_{\\mathrm{escape}}}$ on hitting the ",
      "escape threshold, where $q(a) = d/((1-\\mu)a)$. ",
      "$z = (\\text{empirical} - \\text{analytic})/\\mathrm{s.e.}$; ",
      "the observed $z$ values are consistent with sampling ",
      "fluctuation around the analytic curve."),
    res$n_reps, res$N_escape
  )
  xt <- xtable(df_out, caption = cap, label = "tab:survival",
               align = c("l", "r", "c", "c", "c", "c"))
  print(xt, file = file_out, include.rownames = FALSE,
        booktabs = TRUE, sanitize.text.function = identity,
        caption.placement = "top")
}

## -------------------------------------------------------------------------
## Experiment 3: bounded birth-rate support (Figure 3, Table 3).
## -------------------------------------------------------------------------
experiment_3 <- function(lam = 1, mu = 0.3, a_max = 2,
                         n_reps = 3000, max_pop = 2e5) {
  d_star     <- (1 - mu) * a_max
  d_vals     <- c(0.5, d_star, 2.0)             # below, at, above d*
  t_star_unb <- lam / (1 - mu)
  t_grid     <- seq(0.1, 2.5 * t_star_unb, length.out = 20)

  rows       <- list()
  cap_report <- setNames(numeric(length(d_vals)), as.character(d_vals))

  for (d in d_vals) {
    cat(sprintf("  exp3: d = %.2f\n", d))
    sim <- annealed_mean_truncated(lam, mu, d, a_max, t_grid,
                                   t_max  = max(t_grid),
                                   n_reps = n_reps,
                                   max_pop = max_pop)
    cap_report[as.character(d)] <- sim$capped_fraction
    ana <- truncated_exp_mean(lam, mu, d, a_max, t_grid)
    rows[[as.character(d)]] <- data.frame(t = t_grid,
                                          empirical = sim$mean,
                                          analytic  = ana, d = d)
  }
  if (any(cap_report > 0)) {
    cat("  exp3 capped fraction (max_pop reached):\n")
    print(round(cap_report, 4))
  }
  list(lam = lam, mu = mu, a_max = a_max, d_star = d_star,
       d_vals = d_vals, n_reps = n_reps,
       t_star_unbounded = t_star_unb,
       capped_fraction  = cap_report,
       data             = do.call(rbind, rows))
}

plot_exp3 <- function(res) {
  df       <- res$data
  d_colors <- setNames(c(PAL$green, PAL$orange, PAL$red),
                       as.character(res$d_vals))

  ggplot() +
    geom_line(data = df, aes(t, analytic, colour = factor(d)),
              linewidth = 0.7) +
    geom_point(data = df, aes(t, empirical, colour = factor(d)),
               shape = 21, fill = "white", size = 2) +
    geom_vline(xintercept = res$t_star_unbounded, linetype = "dashed",
               colour = PAL$gray) +
    scale_colour_manual(values = d_colors,
                        labels = paste0("d = ", res$d_vals),
                        name   = "death rate") +
    scale_y_log10() +
    labs(x = "t", y = expression(E*"["*X[1](t)*"]"),
         title = sprintf("Bounded support: a_* = %g, d* = %.2f",
                         res$a_max, res$d_star)) +
    theme_paper()
}

table_exp3 <- function(res, file_out = "table3_bounded.tex") {
  df    <- res$data
  picks <- df %>%
    group_by(d) %>%
    arrange(t) %>%
    slice(c(1, 5, 10, 14, 19)) %>%
    ungroup() %>%
    mutate(regime = case_when(
      d < res$d_star               ~ "$d < d^*$",
      abs(d - res$d_star) < 1e-6   ~ "$d = d^*$",
      TRUE                         ~ "$d > d^*$"
    )) %>%
    select(t, d, regime, empirical, analytic)

  cap <- sprintf(
    paste0(
      "Annealed mean $\\mathbb{E}[X_1(t)]$ under truncated-exponential ",
      "birth rate, $a \\in [0, a_{\\max}]$ with $a_{\\max} = %g$ ",
      "($\\lambda = 1$, $\\mu = 0.3$, $n_{\\mathrm{rep}} = %d$). ",
      "Critical death rate is $d^* = (1-\\mu)a_* = %.2f$."),
    res$a_max, res$n_reps, res$d_star
  )
  xt <- xtable(picks,
               digits  = c(0, 3, 2, 0, 3, 3),
               caption = cap, label = "tab:bounded",
               align   = c("l", "r", "r", "l", "r", "r"))
  print(xt, file = file_out, include.rownames = FALSE,
        booktabs = TRUE, sanitize.text.function = identity,
        caption.placement = "top")
}

## -------------------------------------------------------------------------
## Experiment 4: random independent death rate (Figure 4).
## -------------------------------------------------------------------------
experiment_4 <- function(lam = 1, mu = 0.3, n_reps = 4000,
                         max_pop = 2e5) {
  delta_vals <- c(0.5, 1, 2)
  t_star     <- lam / (1 - mu)
  t_grid     <- seq(0.05, 0.95 * t_star, length.out = 20)

  rows       <- list()
  cap_report <- setNames(numeric(length(delta_vals)),
                         as.character(delta_vals))

  for (delta in delta_vals) {
    cat(sprintf("  exp4: delta = %.2f\n", delta))
    ## local() captures `delta` by value (avoids R's lazy-evaluation
    ## binding to the loop variable's final value).
    sampler <- local({ del <- delta; function() rexp(1, rate = del) })
    L_G     <- local({ del <- delta; function(t) del / (del + t) })

    sim <- annealed_mean(lam, mu, d = NA, t_grid,
                         t_max = 0.95 * t_star,
                         n_reps = n_reps, max_pop = max_pop,
                         d_sampler = sampler)
    cap_report[as.character(delta)] <- sim$capped_fraction

    ana <- analytic_mean_random_d(lam, mu, t_grid, L_G)
    rows[[as.character(delta)]] <- data.frame(t = t_grid,
                                              empirical = sim$mean,
                                              analytic  = ana,
                                              delta     = delta)
  }
  if (any(cap_report > 0)) {
    cat("  exp4 capped fraction (max_pop reached):\n")
    print(round(cap_report, 4))
  }
  list(lam = lam, mu = mu, t_star = t_star, n_reps = n_reps,
       delta_vals       = delta_vals,
       capped_fraction  = cap_report,
       data             = do.call(rbind, rows))
}

plot_exp4 <- function(res) {
  df       <- res$data
  d_colors <- setNames(c(PAL$blue, PAL$orange, PAL$green),
                       as.character(res$delta_vals))
  ggplot() +
    geom_line(data = df, aes(t, analytic, colour = factor(delta)),
              linewidth = 0.7) +
    geom_point(data = df, aes(t, empirical, colour = factor(delta)),
               shape = 21, fill = "white", size = 2) +
    geom_vline(xintercept = res$t_star, linetype = "dashed",
               colour = PAL$gray) +
    scale_colour_manual(values = d_colors,
                        labels = paste0("delta = ", res$delta_vals),
                        name   = NULL) +
    scale_y_log10() +
    labs(x = "t", y = expression(E*"["*X[1](t)*"]"),
         title = "Random independent death") +
    theme_paper()
}

## -------------------------------------------------------------------------
## Experiment 5: heavy-tail bias scaling (Figure 5).
##   At t = 0.9 t_*, sweep n_rep over a logarithmic range and record
##   the empirical/analytic ratio across multiple trials per n.
## -------------------------------------------------------------------------
experiment_5 <- function(lam = 1, mu = 0.3, d = 0,
                         n_vals = c(200, 500, 1000, 2000, 5000,
                                    10000, 20000, 50000),
                         n_trials = 5, max_pop = 1e6) {
  t_star <- lam / (1 - mu)
  t_eval <- 0.9 * t_star
  nu     <- lam / ((1 - mu) * t_eval)
  ana    <- analytic_mean_fixed_d(lam, mu, d, t_eval)

  rows         <- list()
  k            <- 0L
  total_capped <- 0
  total_reps   <- 0L

  for (n in n_vals) {
    cat(sprintf("  exp5: n = %d\n", n))
    for (trial in seq_len(n_trials)) {
      sim <- annealed_mean(lam, mu, d, t_eval, t_max = t_eval,
                           n_reps = n, max_pop = max_pop)
      total_capped <- total_capped + sim$capped_fraction * n
      total_reps   <- total_reps   + n
      k            <- k + 1L
      rows[[k]]    <- data.frame(n = n, trial = trial,
                                 ratio = sim$mean / ana)
    }
  }
  cap_frac <- if (total_reps > 0) total_capped / total_reps else 0
  if (cap_frac > 0) {
    cat(sprintf("  exp5 overall capped fraction: %.4f\n", cap_frac))
  }
  list(lam = lam, mu = mu, d = d, t_eval = t_eval, nu = nu,
       analytic = ana, n_trials = n_trials,
       capped_fraction = cap_frac,
       data = do.call(rbind, rows))
}

plot_exp5 <- function(res) {
  df   <- res$data
  summ <- df %>%
    group_by(n) %>%
    summarise(mean_ratio = mean(ratio),
              se = sd(ratio) / sqrt(n_distinct(trial)),
              .groups = "drop")

  ggplot() +
    geom_point(data = df, aes(n, ratio),
               colour = PAL$gray, size = 1.6, alpha = 0.6) +
    geom_errorbar(data = summ,
                  aes(n,
                      ymin = mean_ratio - 1.96 * se,
                      ymax = mean_ratio + 1.96 * se),
                  width = 0.06, colour = PAL$blue) +
    geom_line(data = summ, aes(n, mean_ratio),
              colour = PAL$blue, linewidth = 0.7) +
    geom_point(data = summ, aes(n, mean_ratio),
               colour = PAL$blue, size = 2.2) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    scale_x_log10(breaks = unique(df$n),
                  labels = trans_format("log10",
                                        math_format(10^.x))) +
    labs(x = expression(n[rep]),
         y = expression("empirical / analytic at " ~ t == 0.9 * t["*"]),
         title = sprintf(
           "Heavy-tail bias scaling at t = 0.9 t_*  (nu = %.2f)",
           res$nu)) +
    theme_paper()
}

## -------------------------------------------------------------------------
## Experiment 6: small-d asymptotic verification (Figure 6, Table 4).
##   Compare exact survival to leading expansion
##   1 - lam c (1 - gamma - log(lam c)),  c = d/(1-mu).
## -------------------------------------------------------------------------
EULER_GAMMA <- 0.5772156649015329

survival_small_d <- function(lam, mu, d) {
  if (d <= 0) return(1)
  cc <- d / (1 - mu)
  lc <- lam * cc
  1 - lc * (1 - EULER_GAMMA - log(lc))
}

experiment_6 <- function(lam = 1, mu = 0.3) {
  d_vals <- 10 ^ seq(-4, 0, length.out = 25)
  ex     <- vapply(d_vals, function(d) analytic_survival(lam, mu, d),
                   numeric(1))
  sd_    <- vapply(d_vals, function(d) survival_small_d(lam, mu, d),
                   numeric(1))
  data.frame(d = d_vals, exact = ex, small_d = sd_,
             abs_err = abs(ex - sd_),
             rel_err = abs(ex - sd_) / pmax(ex, 1e-15))
}

plot_exp6 <- function(df) {
  long <- pivot_longer(df, c(exact, small_d),
                       names_to  = "curve",
                       values_to = "p")
  long$curve <- factor(long$curve,
                       levels = c("exact", "small_d"),
                       labels = c("exact", "small-d expansion"))

  ggplot(long, aes(d, p, colour = curve, linetype = curve)) +
    geom_line(linewidth = 0.7) +
    scale_x_log10(labels = trans_format("log10",
                                        math_format(10^.x))) +
    scale_colour_manual(values = c("exact"             = PAL$blue,
                                   "small-d expansion" = PAL$orange),
                        name = NULL) +
    scale_linetype_manual(values = c("exact"             = "solid",
                                     "small-d expansion" = "dashed"),
                          name = NULL) +
    coord_cartesian(ylim = c(0, 1.05)) +
    labs(x = "d", y = "Pr(eventual survival)",
         title = "Small-d expansion vs exact") +
    theme_paper()
}

table_exp6 <- function(df, file_out = "table4_smalld.tex") {
  picks <- df[seq(1, nrow(df), by = 3), ]
  out <- data.frame(
    d      = sprintf("%.2e",  picks$d),
    Exact  = sprintf("$%.6f$", picks$exact),
    SmallD = sprintf("$%.6f$", picks$small_d),
    RelErr = sprintf("$%.2e$", picks$rel_err),
    stringsAsFactors = FALSE
  )
  cap <- paste0(
    "Validation of the small-$d$ asymptotic expansion in ",
    "Remark~\\ref{rem:small-d}: the leading-order approximation ",
    "$1 - \\lambda c\\,(1 - \\gamma - \\log(\\lambda c))$ ",
    "with $c = d/(1-\\mu)$ approaches the exact survival ",
    "probability $e^{-\\lambda c} - \\lambda c\\,E_1(\\lambda c)$ ",
    "as $d \\to 0^+$. $\\lambda = 1$, $\\mu = 0.3$."
  )
  xt <- xtable(out, caption = cap, label = "tab:smalld",
               align = c("l", "r", "c", "c", "c"))
  print(xt, file = file_out, include.rownames = FALSE,
        booktabs = TRUE, sanitize.text.function = identity,
        caption.placement = "top")
}

## =========================================================================
## DRIVER
## =========================================================================
run_all <- function(out_dir = ".") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  old_wd <- setwd(out_dir); on.exit(setwd(old_wd))

  cat("Running Experiment 1 (annealed mean, fixed d)...\n")
  e1 <- experiment_1()
  p1 <- plot_exp1(e1)
  ## Stack the two panels vertically (gridExtra is a hard dependency
  ## auto-installed at the top of the file).
  pdf("fig1_annealed_mean.pdf", width = 6.5, height = 7)
  gridExtra::grid.arrange(p1$linear, p1$log, ncol = 1)
  dev.off()
  table_exp1(e1)

  cat("Running Experiment 2 (eventual survival)...\n")
  e2 <- experiment_2()
  ggsave("fig2_survival.pdf", plot_exp2(e2), width = 5.5, height = 4)
  table_exp2(e2)

  cat("Running Experiment 3 (bounded support)...\n")
  e3 <- experiment_3()
  ggsave("fig3_bounded.pdf", plot_exp3(e3), width = 6, height = 4)
  table_exp3(e3)

  cat("Running Experiment 4 (random independent death)...\n")
  e4 <- experiment_4()
  ggsave("fig4_random_d.pdf", plot_exp4(e4), width = 6, height = 4)

  cat("Running Experiment 5 (heavy-tail bias scaling)...\n")
  e5 <- experiment_5()
  ggsave("fig5_bias_scaling.pdf", plot_exp5(e5),
         width = 5.5, height = 4)

  cat("Running Experiment 6 (small-d expansion)...\n")
  e6 <- experiment_6()
  ggsave("fig6_smalld.pdf", plot_exp6(e6), width = 5.5, height = 4)
  table_exp6(e6)

  cat("\nAll figures and tables written to ",
      normalizePath(out_dir), "\n", sep = "")
  invisible(list(e1 = e1, e2 = e2, e3 = e3, e4 = e4, e5 = e5, e6 = e6))
}

## -------------------------------------------------------------------------
## Run all experiments only when invoked as a script (Rscript simulations.R),
## not when sourced interactively.
## -------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  run_all()
}
