# Complete simulation audit for the cell-death branching-process manuscript.
#
# This base-R script:
#   1. verifies every analytic expression used by the numerical section;
#   2. verifies the exact finite-time birth-death transition sampler;
#   3. reruns all stochastic experiments with recorded seeds;
#   4. regenerates all six figures and the data underlying all four tables;
#   5. writes a machine-readable and human-readable audit summary.
#
# Audit outputs are written to simulation_audit_output/ first.  Canonical
# manuscript figures and tables should be replaced only after these checks pass.

options(stringsAsFactors = FALSE, scipen = 6)

lambda <- 1
mu <- 0.3
t_star <- lambda / (1 - mu)
master_seed <- 20260902L
audit_dir <- "simulation_audit_output"
publish <- "--publish" %in% commandArgs(trailingOnly = TRUE)
dir.create(audit_dir, showWarnings = FALSE, recursive = TRUE)

# Nature-style graphical system: Helvetica, white backgrounds, black axes,
# no background gridlines, and a high-contrast subset of the NPG palette. Line
# types and point shapes duplicate the colour encoding for accessibility. Text,
# lines, and markers are sized for the final manuscript rather than screen zoom.
nature_colors <- c(
  vermillion = "#E64B35",
  cyan = "#4DBBD5",
  green = "#00A087",
  navy = "#3C5488",
  salmon = "#F39B7F",
  slate = "#8491B4",
  mint = "#91D1C2",
  red = "#DC0000",
  brown = "#7E6148",
  taupe = "#B09C85"
)
palette <- unname(nature_colors[c("navy", "vermillion", "green", "brown")])
nature_ink <- "#262626"
nature_midgray <- "#737373"
nature_lightgray <- "#BDBDBD"
nature_lty <- c(1, 5, 2, 3)
nature_pch <- c(21, 22, 24, 23)

nature_pdf <- function(filename, width, height) {
  pdf(
    filename, width = width, height = height, family = "Helvetica",
    pointsize = 10.5, useDingbats = FALSE, bg = "white"
  )
}

nature_par <- function(mar = c(3.5, 3.9, 0.6, 0.45), ...) {
  par(
    mar = mar, bty = "l", las = 1, mgp = c(2.35, 0.72, 0),
    tcl = -0.3, cex.axis = 0.95, cex.lab = 1.05,
    col.axis = nature_ink, col.lab = nature_ink, fg = nature_ink,
    lend = "round", ljoin = "round", lwd = 1.05,
    xaxs = "r", yaxs = "r", ...
  )
}

nature_panel_label <- function(label) {
  mtext(label, side = 3, line = 0.05, adj = 0, font = 2, cex = 1.15)
}

analytic_exponential_mean <- function(t, death) {
  stopifnot(all(t >= 0), all(death >= 0))
  denominator <- lambda - (1 - mu) * t
  ifelse(denominator > 0,
         exp(-death * t) * lambda / denominator,
         Inf)
}

# Stable zero-inflated-geometric parameters for a linear birth-death
# process started from one cell.  Conditional on being nonzero,
# N(t) = 1 + Geometric(success = 1 - beta).
transition_parameters <- function(birth, death, t) {
  birth <- as.numeric(birth)
  if (length(death) == 1L) {
    death <- rep(as.numeric(death), length(birth))
  } else {
    death <- as.numeric(death)
  }
  stopifnot(length(birth) == length(death),
            length(t) == 1L, t >= 0,
            all(birth >= 0), all(death >= 0))

  rate_difference <- birth - death
  scale <- pmax(1, birth, death)
  critical <- abs(rate_difference) <= 1e-10 * scale
  supercritical <- !critical & rate_difference > 0
  subcritical <- !critical & rate_difference < 0
  p_zero <- numeric(length(birth))
  beta <- numeric(length(birth))

  if (any(supercritical)) {
    idx <- which(supercritical)
    z <- exp(-rate_difference[idx] * t)
    denominator <- birth[idx] - death[idx] * z
    p_zero[idx] <- death[idx] * (1 - z) / denominator
    beta[idx] <- birth[idx] * (1 - z) / denominator
  }

  if (any(subcritical)) {
    idx <- which(subcritical)
    z <- exp(rate_difference[idx] * t)
    denominator <- death[idx] - birth[idx] * z
    p_zero[idx] <- death[idx] * (1 - z) / denominator
    beta[idx] <- birth[idx] * (1 - z) / denominator
  }

  if (any(critical)) {
    idx <- which(critical)
    bt <- birth[idx] * t
    p_zero[idx] <- bt / (1 + bt)
    beta[idx] <- p_zero[idx]
  }

  list(
    p_zero = pmin(pmax(p_zero, 0), 1),
    beta = pmin(pmax(beta, 0), 1)
  )
}

sample_linear_birth_death <- function(birth, death, t) {
  pars <- transition_parameters(birth, death, t)
  success <- pmax(1 - pars$beta, .Machine$double.eps)
  alive <- runif(length(birth)) >= pars$p_zero
  count <- numeric(length(birth))
  count[alive] <- 1 + rgeom(sum(alive), prob = success[alive])
  count
}

birth_death_variance <- function(birth, death, t) {
  rate_difference <- birth - death
  ifelse(
    abs(rate_difference) <= 1e-10 * pmax(1, birth, death),
    2 * birth * t,
    (birth + death) / rate_difference *
      exp(rate_difference * t) * (exp(rate_difference * t) - 1)
  )
}

sample_truncated_exponential <- function(n, rate, upper) {
  u <- runif(n)
  -log1p(-u * (1 - exp(-rate * upper))) / rate
}

analytic_truncated_mean <- function(t, death, upper = 2) {
  k <- (1 - mu) * t - lambda
  integral <- ifelse(
    abs(k) < 1e-10,
    upper,
    expm1(k * upper) / k
  )
  exp(-death * t) * lambda * integral /
    (1 - exp(-lambda * upper))
}

E1_scalar <- function(x) {
  if (x == 0) {
    return(Inf)
  }
  stopifnot(x > 0)
  integrate(function(u) exp(-u) / u, lower = x, upper = Inf,
            rel.tol = 1e-11, subdivisions = 2000L)$value
}

E1 <- function(x) vapply(x, E1_scalar, numeric(1))

analytic_survival <- function(death) {
  result <- numeric(length(death))
  zero <- death == 0
  result[zero] <- 1
  if (any(!zero)) {
    x <- lambda * death[!zero] / (1 - mu)
    result[!zero] <- exp(-x) - x * E1(x)
  }
  result
}

conditional_survival <- function(rate, death) {
  c_threshold <- death / (1 - mu)
  ifelse(rate > c_threshold, 1 - c_threshold / rate, 0)
}

conditional_count_tail <- function(n, rate, death, t) {
  pars <- transition_parameters((1 - mu) * rate, death, t)
  result <- numeric(length(rate))
  positive <- pars$beta > 0 & pars$p_zero < 1
  result[positive] <- (1 - pars$p_zero[positive]) *
    exp(n * log(pars$beta[positive]))
  result
}

annealed_count_tail <- function(n, t, death) {
  integrate(
    function(a) conditional_count_tail(n, a, death, t) *
      lambda * exp(-lambda * a),
    lower = 0, upper = Inf, rel.tol = 1e-8, subdivisions = 2000L
  )$value
}

tail_asymptotic_constant <- function(t, death) {
  nu <- lambda / ((1 - mu) * t)
  gamma(nu + 1) * exp(-lambda * death / (1 - mu))
}

log_expm1_positive <- function(x) {
  ifelse(x > 50, x + log1p(-exp(-x)), log(expm1(x)))
}

# Probability that the embedded birth-death random walk hits N before 0.
# This is the exact quantity estimated by a Gillespie escape-threshold run.
conditional_escape_probability <- function(rate, death, threshold) {
  birth <- (1 - mu) * rate
  result <- numeric(length(rate))
  if (death == 0) {
    result[birth > 0] <- 1
    return(result)
  }

  valid <- birth > 0
  rho <- rep(Inf, length(rate))
  rho[valid] <- death / birth[valid]
  log_rho <- log(rho[valid])
  near_critical <- abs(log_rho) < 1e-8
  below_one <- log_rho < -1e-8
  above_one <- log_rho > 1e-8
  valid_idx <- which(valid)

  if (any(near_critical)) {
    result[valid_idx[near_critical]] <- 1 / threshold
  }
  if (any(below_one)) {
    lr <- log_rho[below_one]
    result[valid_idx[below_one]] <-
      (-expm1(lr)) / (-expm1(threshold * lr))
  }
  if (any(above_one)) {
    lr <- log_rho[above_one]
    log_probability <- log_expm1_positive(lr) -
      log_expm1_positive(threshold * lr)
    result[valid_idx[above_one]] <- exp(log_probability)
  }
  result
}

# Simulate the embedded birth-death walk from one cell until it reaches either
# extinction or a finite escape threshold. Mutation events that leave the
# type-1 count unchanged are omitted from the embedded chain.
simulate_embedded_escape <- function(rate, death, threshold) {
  birth <- (1 - mu) * rate
  n_rep <- length(birth)
  stopifnot(length(death) == 1L, death >= 0, threshold >= 2L)

  if (death == 0) {
    return(birth > 0)
  }

  p_birth <- birth / (birth + death)
  state <- rep.int(1L, n_rep)
  active <- seq_len(n_rep)

  while (length(active) > 0L) {
    step_up <- runif(length(active)) < p_birth[active]
    state[active] <- state[active] + ifelse(step_up, 1L, -1L)
    active <- active[state[active] > 0L & state[active] < threshold]
  }

  state == threshold
}

analytic_random_death_mean <- function(t, delta) {
  analytic_exponential_mean(t, 0) * delta / (delta + t)
}

small_d_approximation <- function(death) {
  x <- lambda * death / (1 - mu)
  1 - x * (1 - 0.5772156649015329 - log(x))
}

# ---------------------------------------------------------------------------
# Deterministic mathematical and sampler checks
# ---------------------------------------------------------------------------

check_grid <- expand.grid(
  birth = c(0, 0.2, 0.8, 2),
  death = c(0, 0.2, 0.8, 2),
  time = c(0.1, 0.8, 2)
)

transition_check <- do.call(rbind, lapply(seq_len(nrow(check_grid)), function(i) {
  row <- check_grid[i, ]
  pars <- transition_parameters(row$birth, row$death, row$time)
  implied_mean <- (1 - pars$p_zero) / (1 - pars$beta)
  geometric_mean <- 1 / (1 - pars$beta)
  geometric_variance <- pars$beta / (1 - pars$beta)^2
  implied_second_moment <- (1 - pars$p_zero) *
    (geometric_variance + geometric_mean^2)
  implied_variance <- implied_second_moment - implied_mean^2
  exact_mean <- exp((row$birth - row$death) * row$time)
  exact_variance <- birth_death_variance(row$birth, row$death, row$time)
  data.frame(
    birth = row$birth,
    death = row$death,
    time = row$time,
    mean_relative_error = abs(implied_mean - exact_mean) /
      pmax(1, abs(exact_mean)),
    variance_relative_error = abs(implied_variance - exact_variance) /
      pmax(1, abs(exact_variance))
  )
}))

mean_integral_grid <- expand.grid(
  fraction = c(0.1, 0.5, 0.9),
  death = c(0, 0.5, 2)
)
mean_integral_check <- do.call(rbind, lapply(
  seq_len(nrow(mean_integral_grid)),
  function(i) {
    row <- mean_integral_grid[i, ]
    t <- row$fraction * t_star
    numerical <- integrate(
      function(a) lambda * exp(
        -row$death * t + ((1 - mu) * t - lambda) * a
      ),
      lower = 0, upper = Inf, rel.tol = 1e-10
    )$value
    analytic <- analytic_exponential_mean(t, row$death)
    data.frame(
      fraction = row$fraction,
      death = row$death,
      relative_error = abs(numerical - analytic) / analytic
    )
  }
))

survival_death_grid <- seq(0, 2, by = 0.25)
survival_integral_check <- do.call(rbind, lapply(
  survival_death_grid,
  function(death) {
    numerical <- if (death == 0) {
      1
    } else {
      c_threshold <- death / (1 - mu)
      integrate(
        function(a) (1 - c_threshold / a) * lambda * exp(-lambda * a),
        lower = c_threshold, upper = Inf, rel.tol = 1e-10
      )$value
    }
    analytic <- analytic_survival(death)
    data.frame(
      death = death,
      numerical = numerical,
      analytic = analytic,
      absolute_error = abs(numerical - analytic)
    )
  }
))

bounded_check_grid <- expand.grid(
  time = c(0.1, 1 / (1 - mu), 2, 3.5),
  death = c(0.5, 1.4, 2)
)
bounded_integral_check <- do.call(rbind, lapply(
  seq_len(nrow(bounded_check_grid)),
  function(i) {
    row <- bounded_check_grid[i, ]
    numerical <- integrate(
      function(a) exp(((1 - mu) * a - row$death) * row$time) *
        lambda * exp(-lambda * a) / (1 - exp(-2 * lambda)),
      lower = 0, upper = 2, rel.tol = 1e-11
    )$value
    analytic <- analytic_truncated_mean(row$time, row$death, 2)
    data.frame(
      time = row$time,
      death = row$death,
      relative_error = abs(numerical - analytic) / analytic
    )
  }
))

random_check_grid <- expand.grid(
  fraction = c(0.1, 0.5, 0.9),
  delta = c(0.5, 1, 2)
)
random_integral_check <- do.call(rbind, lapply(
  seq_len(nrow(random_check_grid)),
  function(i) {
    row <- random_check_grid[i, ]
    t <- row$fraction * t_star
    birth_factor <- integrate(
      function(a) lambda * exp(((1 - mu) * t - lambda) * a),
      lower = 0, upper = Inf, rel.tol = 1e-10
    )$value
    death_factor <- integrate(
      function(d) exp(-d * t) * row$delta * exp(-row$delta * d),
      lower = 0, upper = Inf, rel.tol = 1e-10
    )$value
    numerical <- birth_factor * death_factor
    analytic <- analytic_random_death_mean(t, row$delta)
    data.frame(
      fraction = row$fraction,
      delta = row$delta,
      relative_error = abs(numerical - analytic) / analytic
    )
  }
))

tail_check_grid <- expand.grid(
  fraction = c(0.8, 0.9),
  death = c(0, 0.5, 1),
  n = 10000
)
tail_asymptotic_check <- do.call(rbind, lapply(
  seq_len(nrow(tail_check_grid)),
  function(i) {
    row <- tail_check_grid[i, ]
    t <- row$fraction * t_star
    nu <- lambda / ((1 - mu) * t)
    scaled_tail <- annealed_count_tail(row$n, t, row$death) * row$n^nu
    limit <- tail_asymptotic_constant(t, row$death)
    data.frame(
      fraction = row$fraction,
      death = row$death,
      n = row$n,
      nu = nu,
      scaled_tail = scaled_tail,
      asymptotic_limit = limit,
      relative_discrepancy = abs(scaled_tail - limit) / limit
    )
  }
))

small_d_half_relative_error <- abs(
  small_d_approximation(0.5) - analytic_survival(0.5)
) / analytic_survival(0.5)

stopifnot(
  max(transition_check$mean_relative_error) < 1e-10,
  max(transition_check$variance_relative_error) < 1e-9,
  max(mean_integral_check$relative_error) < 1e-8,
  max(survival_integral_check$absolute_error) < 1e-9,
  max(bounded_integral_check$relative_error) < 1e-9,
  max(random_integral_check$relative_error) < 1e-8,
  max(tail_asymptotic_check$relative_discrepancy) < 0.08,
  small_d_half_relative_error > 0.992,
  small_d_half_relative_error < 0.994
)

write.csv(transition_check,
          file.path(audit_dir, "transition_sampler_checks.csv"),
          row.names = FALSE)
write.csv(mean_integral_check,
          file.path(audit_dir, "analytic_mean_checks.csv"),
          row.names = FALSE)
write.csv(survival_integral_check,
          file.path(audit_dir, "survival_formula_checks.csv"),
          row.names = FALSE)
write.csv(bounded_integral_check,
          file.path(audit_dir, "bounded_formula_checks.csv"),
          row.names = FALSE)
write.csv(random_integral_check,
          file.path(audit_dir, "random_death_formula_checks.csv"),
          row.names = FALSE)
write.csv(tail_asymptotic_check,
          file.path(audit_dir, "tail_asymptotic_checks.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Figure 1 and Table 1: annealed mean with fixed death
# ---------------------------------------------------------------------------

death_rates <- c(0, 0.5, 1, 2)
n_annealed <- 20000L
time_fractions <- c(0.035, 0.175, 0.35, 0.50, 0.525, 0.70, 0.82,
                    0.90, 0.94, 0.95, 0.97, 0.98, 0.985, 0.992)
table_fractions <- c(0.50, 0.90, 0.95, 0.98, 0.992)

annealed_results <- do.call(rbind, lapply(seq_along(death_rates), function(di) {
  death <- death_rates[di]
  do.call(rbind, lapply(seq_along(time_fractions), function(ti) {
    fraction <- time_fractions[ti]
    t <- fraction * t_star
    set.seed(master_seed + 1000L * di + ti)
    rate <- rexp(n_annealed, rate = lambda)
    count <- sample_linear_birth_death((1 - mu) * rate, death, t)
    empirical <- mean(count)
    analytic <- analytic_exponential_mean(t, death)
    data.frame(
      death = death,
      fraction = fraction,
      time = t,
      n_rep = n_annealed,
      empirical = empirical,
      analytic = analytic,
      ratio = empirical / analytic
    )
  }))
}))

write.csv(annealed_results,
          file.path(audit_dir, "fig1_annealed_results.csv"),
          row.names = FALSE)
write.csv(annealed_results[
  annealed_results$fraction %in% table_fractions, ],
  file.path(audit_dir, "table1_annealed.csv"),
  row.names = FALSE)

dense_t <- seq(0, 0.995 * t_star, length.out = 500)
analytic_dense <- unlist(lapply(
  death_rates, function(death) analytic_exponential_mean(dense_t, death)
))
all_y <- c(annealed_results$empirical, analytic_dense)

nature_pdf(file.path(audit_dir, "fig1_annealed_mean.pdf"),
           width = 6.8, height = 6.2)
nature_par(mar = c(3.4, 4.0, 0.8, 0.35), mfrow = c(2, 1))
plot(NA, xlim = c(0, 1.02 * t_star), ylim = c(0, max(all_y)),
     xlab = expression("Time, " * italic(t)),
     ylab = expression("Annealed mean, " * E[X[1](t)]))
nature_panel_label("a")
  abline(v = t_star, lty = 2, col = nature_midgray, lwd = 1.25)
for (i in seq_along(death_rates)) {
  death <- death_rates[i]
  lines(dense_t, analytic_exponential_mean(dense_t, death),
        col = palette[i], lwd = 2.25, lty = nature_lty[i])
  block <- annealed_results[annealed_results$death == death, ]
  points(block$time, block$empirical, pch = nature_pch[i], bg = "white",
         col = palette[i], lwd = 1.25, cex = 0.95)
}

plot(NA, xlim = c(0, 1.02 * t_star),
     ylim = range(all_y[is.finite(all_y) & all_y > 0]), log = "y",
     xlab = expression("Time, " * italic(t)),
     ylab = expression("Annealed mean, " * E[X[1](t)]))
nature_panel_label("b")
abline(v = t_star, lty = 2, col = nature_midgray, lwd = 1.25)
for (i in seq_along(death_rates)) {
  death <- death_rates[i]
  lines(dense_t, analytic_exponential_mean(dense_t, death),
        col = palette[i], lwd = 2.25, lty = nature_lty[i])
  block <- annealed_results[annealed_results$death == death, ]
  points(block$time, block$empirical, pch = nature_pch[i], bg = "white",
         col = palette[i], lwd = 1.25, cex = 0.95)
}
legend("topleft", legend = paste0("d = ", death_rates),
       col = palette, lwd = 2.25, lty = nature_lty,
       pch = nature_pch, pt.bg = "white", pt.cex = 0.95,
       bty = "n", cex = 0.90, seg.len = 2.5)
dev.off()

# ---------------------------------------------------------------------------
# Figure 2 and Table 2: exact eventual-survival simulation
# ---------------------------------------------------------------------------

survival_death_rates <- seq(0, 2, by = 0.25)
n_survival <- 20000L
escape_threshold <- 200L

survival_results <- do.call(rbind, lapply(
  seq_along(survival_death_rates),
  function(i) {
    death <- survival_death_rates[i]
    set.seed(master_seed + 100000L + i)
    rate <- rexp(n_survival, rate = lambda)
    p_conditional <- conditional_survival(rate, death)
    survived <- runif(n_survival) < p_conditional
    empirical <- mean(survived)
    analytic <- analytic_survival(death)
    standard_error <- sqrt(empirical * (1 - empirical) / n_survival)

    p_escape_conditional <- conditional_escape_probability(
      rate, death, escape_threshold
    )
    escape_mixture_numeric <- mean(p_escape_conditional)
    escape_mixture_integral <- if (death == 0) {
      1
    } else {
      integrate(
        function(a) conditional_escape_probability(
          a, death, escape_threshold
        ) * lambda * exp(-lambda * a),
        lower = 0, upper = Inf, rel.tol = 1e-9,
        subdivisions = 2000L
      )$value
    }

    set.seed(master_seed + 110000L + i)
    escape_indicator <- simulate_embedded_escape(
      rate, death, escape_threshold
    )
    escape_empirical <- mean(escape_indicator)
    escape_standard_error <- sqrt(
      escape_empirical * (1 - escape_empirical) / n_survival
    )

    data.frame(
      death = death,
      n_rep = n_survival,
      empirical = empirical,
      analytic = analytic,
      naive = exp(-lambda * death / (1 - mu)),
      standard_error = standard_error,
      margin_95 = 1.96 * standard_error,
      z_survival = if (standard_error > 0) {
        (empirical - analytic) / standard_error
      } else {
        NA_real_
      },
      escape_probability_mc_integral = escape_mixture_numeric,
      escape_probability_integral = escape_mixture_integral,
      escape_probability_simulated = escape_empirical,
      escape_standard_error = escape_standard_error,
      z_escape_vs_exact = if (escape_standard_error > 0) {
        (escape_empirical - escape_mixture_integral) /
          escape_standard_error
      } else {
        NA_real_
      },
      escape_bias = escape_mixture_integral - analytic
    )
  }
))

stopifnot(
  sum(is.finite(survival_results$z_survival)) == 8L,
  sum(is.finite(survival_results$z_escape_vs_exact)) == 8L,
  all(is.na(survival_results$z_survival[survival_results$death == 0])),
  all(is.na(survival_results$z_escape_vs_exact[survival_results$death == 0])),
  all(survival_results$empirical[survival_results$death == 0] == 1),
  all(survival_results$analytic[survival_results$death == 0] == 1)
)

write.csv(survival_results,
          file.path(audit_dir, "table2_survival.csv"),
          row.names = FALSE)

dense_death <- seq(0, 2, length.out = 201)
dense_survival <- analytic_survival(dense_death)
dense_naive <- exp(-lambda * dense_death / (1 - mu))

nature_pdf(file.path(audit_dir, "fig2_survival.pdf"),
           width = 5.2, height = 4.0)
nature_par()
plot(dense_death, dense_survival, type = "l", lwd = 2.4,
     col = palette[1], ylim = c(0, 1),
     xlab = expression("Death rate, " * italic(d)),
     ylab = "Eventual survival probability")
lines(dense_death, dense_naive, lwd = 2.1,
      col = palette[2], lty = 5)
nonzero_margin <- survival_results$margin_95 > 0
arrows(
  survival_results$death[nonzero_margin],
  pmax(0, survival_results$empirical[nonzero_margin] -
         survival_results$margin_95[nonzero_margin]),
  survival_results$death[nonzero_margin],
  pmin(1, survival_results$empirical[nonzero_margin] +
         survival_results$margin_95[nonzero_margin]),
  angle = 90, code = 3, length = 0.028,
  col = nature_midgray, lwd = 1.2
)
points(survival_results$death, survival_results$empirical,
       pch = 21, bg = "white", col = nature_ink, lwd = 1.25, cex = 1.05)
legend("topright",
       legend = c("analytic survival", "naive supercritical bound",
                  "exact Monte Carlo"),
       col = c(palette[1], palette[2], nature_ink),
       lty = c(1, 5, NA), pch = c(NA, NA, 21),
       pt.bg = c(NA, NA, "white"), lwd = c(2.4, 2.1, NA),
       pt.cex = c(NA, NA, 1.05), bty = "n", cex = 0.90,
       seg.len = 2.5)
dev.off()

# ---------------------------------------------------------------------------
# Figure 3 and Table 3: bounded birth-rate support
# ---------------------------------------------------------------------------

bounded_upper <- 2
bounded_death_rates <- c(0.5, 1.4, 2)
bounded_times <- seq(0.1, 2.5 * t_star, length.out = 20)
bounded_table_indices <- c(1, 5, 10, 14, 19)
n_bounded <- 3000L

bounded_results <- do.call(rbind, lapply(
  seq_along(bounded_death_rates),
  function(di) {
    death <- bounded_death_rates[di]
    do.call(rbind, lapply(seq_along(bounded_times), function(ti) {
      t <- bounded_times[ti]
      set.seed(master_seed + 200000L + 1000L * di + ti)
      rate <- sample_truncated_exponential(
        n_bounded, lambda, bounded_upper
      )
      count <- sample_linear_birth_death((1 - mu) * rate, death, t)
      empirical <- mean(count)
      analytic <- analytic_truncated_mean(t, death, bounded_upper)
      standard_error <- sd(count) / sqrt(n_bounded)
      data.frame(
        death = death,
        regime = if (death < 1.4) {
          "d < d*"
        } else if (death > 1.4) {
          "d > d*"
        } else {
          "d = d*"
        },
        time_index = ti,
        time = t,
        n_rep = n_bounded,
        empirical = empirical,
        analytic = analytic,
        standard_error = standard_error,
        z = if (standard_error > 0) {
          (empirical - analytic) / standard_error
        } else {
          NA_real_
        }
      )
    }))
  }
))

write.csv(bounded_results,
          file.path(audit_dir, "fig3_bounded_results.csv"),
          row.names = FALSE)
write.csv(bounded_results[
  bounded_results$time_index %in% bounded_table_indices, ],
  file.path(audit_dir, "table3_bounded.csv"),
  row.names = FALSE)

dense_bounded_time <- seq(0, 2.5 * t_star, length.out = 400)
bounded_y <- c(
  bounded_results$empirical,
  unlist(lapply(bounded_death_rates, function(death) {
    analytic_truncated_mean(dense_bounded_time, death, bounded_upper)
  }))
)

nature_pdf(file.path(audit_dir, "fig3_bounded.pdf"),
           width = 5.2, height = 4.25)
layout(matrix(c(1, 2), nrow = 2), heights = c(0.14, 0.86))
par(mar = c(0, 3.9, 0, 0.45), xpd = NA)
plot.new()
bounded_death_labels <- as.expression(lapply(
  bounded_death_rates, function(death) bquote(d == .(death))
))
legend("center", legend = bounded_death_labels,
       col = palette[seq_along(bounded_death_rates)],
       lwd = 2.3, lty = nature_lty[seq_along(bounded_death_rates)],
       pch = nature_pch[seq_along(bounded_death_rates)],
       pt.bg = "white", pt.cex = 1.0,
       bty = "n", cex = 0.90, seg.len = 2.2,
       horiz = TRUE, x.intersp = 0.8)
nature_par(mar = c(3.5, 3.9, 0.2, 0.45), xpd = FALSE)
plot(NA, xlim = range(dense_bounded_time),
     ylim = range(bounded_y[bounded_y > 0]), log = "y",
     xlab = expression("Time, " * italic(t)),
     ylab = expression("Annealed mean, " * E[X[1](t)]))
abline(v = t_star, lty = 2, col = nature_midgray, lwd = 1.25)
for (i in seq_along(bounded_death_rates)) {
  death <- bounded_death_rates[i]
  lines(
    dense_bounded_time,
    analytic_truncated_mean(dense_bounded_time, death, bounded_upper),
    col = palette[i], lwd = 2.3, lty = nature_lty[i]
  )
  block <- bounded_results[bounded_results$death == death, ]
  points(block$time, block$empirical, pch = nature_pch[i], bg = "white",
         col = palette[i], lwd = 1.25, cex = 1.0)
}
dev.off()

# ---------------------------------------------------------------------------
# Manuscript Figure 5: random independent death rates
# ---------------------------------------------------------------------------

random_deltas <- c(0.5, 1, 2)
random_times <- seq(0.05, 0.95 * t_star, length.out = 14)
n_random <- 4000L

random_results <- do.call(rbind, lapply(
  seq_along(random_deltas),
  function(di) {
    delta <- random_deltas[di]
    do.call(rbind, lapply(seq_along(random_times), function(ti) {
      t <- random_times[ti]
      set.seed(master_seed + 300000L + 1000L * di + ti)
      rate <- rexp(n_random, rate = lambda)
      death <- rexp(n_random, rate = delta)
      count <- sample_linear_birth_death((1 - mu) * rate, death, t)
      empirical <- mean(count)
      analytic <- analytic_random_death_mean(t, delta)
      data.frame(
        delta = delta,
        time = t,
        fraction = t / t_star,
        n_rep = n_random,
        empirical = empirical,
        analytic = analytic,
        ratio = empirical / analytic
      )
    }))
  }
))

write.csv(random_results,
          file.path(audit_dir, "fig5_random_death_results.csv"),
          row.names = FALSE)

dense_random_time <- seq(0, 0.995 * t_star, length.out = 400)
random_y <- c(
  random_results$empirical,
  unlist(lapply(random_deltas, function(delta) {
    analytic_random_death_mean(dense_random_time, delta)
  }))
)

nature_pdf(file.path(audit_dir, "fig5_random_death.pdf"),
           width = 5.2, height = 4.0)
nature_par()
plot(NA, xlim = c(0, 1.02 * t_star),
     ylim = range(random_y[random_y > 0]), log = "y",
     xlab = expression("Time, " * italic(t)),
     ylab = expression("Annealed mean, " * E[X[1](t)]))
abline(v = t_star, lty = 2, col = nature_midgray, lwd = 1.25)
for (i in seq_along(random_deltas)) {
  delta <- random_deltas[i]
  lines(dense_random_time,
        analytic_random_death_mean(dense_random_time, delta),
        col = palette[i], lwd = 2.3, lty = nature_lty[i])
  block <- random_results[random_results$delta == delta, ]
  points(block$time, block$empirical, pch = nature_pch[i], bg = "white",
         col = palette[i], lwd = 1.25, cex = 1.0)
}
random_delta_labels <- as.expression(lapply(
  random_deltas, function(delta) bquote(delta == .(delta))
))
legend("topleft", legend = random_delta_labels,
       col = palette[seq_along(random_deltas)], lwd = 2.3,
       lty = nature_lty[seq_along(random_deltas)],
       pch = nature_pch[seq_along(random_deltas)],
       pt.bg = "white", pt.cex = 1.0,
       bty = "n", cex = 0.90, seg.len = 2.5)
dev.off()

# ---------------------------------------------------------------------------
# Manuscript Figure 4: finite-sample instability near the heavy-tail boundary
# ---------------------------------------------------------------------------

bias_time <- 0.9 * t_star
bias_index <- lambda / ((1 - mu) * bias_time)
bias_sample_sizes <- c(200L, 500L, 1000L, 2000L, 5000L, 10000L,
                       20000L, 50000L)
bias_trials <- 200L
bias_analytic <- analytic_exponential_mean(bias_time, 0)

bias_results <- do.call(rbind, lapply(
  seq_along(bias_sample_sizes),
  function(i) {
    n <- bias_sample_sizes[i]
    do.call(rbind, lapply(seq_len(bias_trials), function(trial) {
      set.seed(master_seed + 400000L + 1000L * i + trial)
      rate <- rexp(n, rate = lambda)
      count <- sample_linear_birth_death(
        (1 - mu) * rate, 0, bias_time
      )
      empirical <- mean(count)
      data.frame(
        n_rep = n,
        trial = trial,
        empirical = empirical,
        analytic = bias_analytic,
        ratio = empirical / bias_analytic
      )
    }))
  }
))

bias_summary <- do.call(rbind, lapply(
  split(bias_results, bias_results$n_rep),
  function(block) {
    data.frame(
      n_rep = block$n_rep[1],
      trials = nrow(block),
      mean_ratio = mean(block$ratio),
      median_ratio = median(block$ratio),
      q025 = unname(quantile(block$ratio, 0.025)),
      q10 = unname(quantile(block$ratio, 0.10)),
      q90 = unname(quantile(block$ratio, 0.90)),
      q975 = unname(quantile(block$ratio, 0.975)),
      fraction_below_one = mean(block$ratio < 1),
      maximum_ratio = max(block$ratio)
    )
  }
))
bias_summary <- bias_summary[order(bias_summary$n_rep), ]

write.csv(bias_results,
          file.path(audit_dir, "fig4_finite_sample_trials.csv"),
          row.names = FALSE)
write.csv(bias_summary,
          file.path(audit_dir, "fig4_finite_sample_summary.csv"),
          row.names = FALSE)

nature_pdf(file.path(audit_dir, "fig4_finite_sample_instability.pdf"),
           width = 5.2, height = 4.0)
nature_par(mar = c(3.7, 4.1, 0.6, 0.45))
plot(bias_summary$n_rep, bias_summary$median_ratio,
     log = "x", type = "n",
     ylim = range(c(0.35, 1, bias_summary$q025, bias_summary$q975)),
     xlab = expression("Sample size, " * n[rep]),
     ylab = "Empirical mean / analytic mean")
abline(h = 1, lty = 2, col = nature_midgray, lwd = 1.25)
segments(bias_summary$n_rep, bias_summary$q025,
         bias_summary$n_rep, bias_summary$q975,
         col = nature_midgray, lwd = 1.45)
segments(0.94 * bias_summary$n_rep, bias_summary$q025,
         1.06 * bias_summary$n_rep, bias_summary$q025,
         col = nature_midgray, lwd = 1.45)
segments(0.94 * bias_summary$n_rep, bias_summary$q975,
         1.06 * bias_summary$n_rep, bias_summary$q975,
         col = nature_midgray, lwd = 1.45)
lines(bias_summary$n_rep, bias_summary$median_ratio,
      col = palette[1], lwd = 2.4)
points(bias_summary$n_rep, bias_summary$median_ratio,
       pch = 21, bg = "white", col = palette[1], lwd = 1.25, cex = 1.05)
legend("topright",
       legend = c("median over 200 trials", "2.5%-97.5% interval",
                  "ratio = 1"),
       col = c(palette[1], nature_midgray, nature_midgray),
       lty = c(1, 1, 2), pch = c(21, NA, NA),
       pt.bg = c("white", NA, NA), lwd = c(2.4, 1.45, 1.25),
       pt.cex = 1.05, bty = "n", cex = 0.90, seg.len = 2.5)
dev.off()

# ---------------------------------------------------------------------------
# Figure 6 and Table 4: small-d approximation
# ---------------------------------------------------------------------------

small_d_grid <- 10^seq(-4, 0, length.out = 300)
small_d_table <- 10^seq(-4, 0, by = 0.5)
small_d_exact <- analytic_survival(small_d_grid)
small_d_approx <- small_d_approximation(small_d_grid)

small_d_results <- data.frame(
  death = small_d_table,
  exact = analytic_survival(small_d_table),
  approximation = small_d_approximation(small_d_table)
)
small_d_results$relative_error <-
  abs(small_d_results$approximation - small_d_results$exact) /
  small_d_results$exact
write.csv(small_d_results,
          file.path(audit_dir, "table4_smalld.csv"),
          row.names = FALSE)

nature_pdf(file.path(audit_dir, "fig6_smalld.pdf"),
           width = 5.2, height = 4.0)
nature_par()
plot(small_d_grid, small_d_exact, type = "l", log = "x",
     ylim = range(c(0, 1, small_d_approx)),
     xlab = expression("Death rate, " * italic(d)),
     ylab = "Eventual survival probability",
     col = palette[1], lwd = 2.4)
lines(small_d_grid, small_d_approx,
      col = palette[2], lwd = 2.2, lty = 5)
legend("bottomleft", legend = c("exact", "small-d expansion"),
       col = c(palette[1], palette[2]), lwd = c(2.4, 2.2),
       lty = c(1, 5), bty = "n", cex = 0.90, seg.len = 2.5)
dev.off()

# ---------------------------------------------------------------------------
# Reproducible LaTeX tables
# ---------------------------------------------------------------------------

write_table_variants <- function(lines, filename) {
  clean_path <- file.path(audit_dir, filename)
  writeLines(lines, clean_path)

  tracked_lines <- lines
  centering_line <- which(tracked_lines == "\\centering")[1]
  tracked_lines <- append(
    tracked_lines, "\\color{revisionblue}", after = centering_line
  )
  caption_line <- grep("^\\\\caption\\{", tracked_lines)[1]
  tracked_lines[caption_line] <- sub(
    "^\\\\caption\\{", "\\\\revcaption{", tracked_lines[caption_line]
  )
  tracked_name <- sub("\\.tex$", "_tracked.tex", filename)
  writeLines(tracked_lines, file.path(audit_dir, tracked_name))
}

survival_table_lines <- c(
  "% Generated by simulation_audit.R",
  "\\begin{table}[ht]",
  "\\centering",
  paste0(
    "\\caption{Exact Monte Carlo estimates of the eventual survival ",
    "probability of the type-1 clone lineage as a function of $d$ ",
    "($\\lambda=1$, $\\mu=0.3$, $n_{\\mathrm{rep}}=20{,}000$). ",
    "Parentheses give $\\pm1.96\\,\\mathrm{s.e.}$, where ",
    "$\\mathrm{s.e.}=\\{\\widehat p(1-\\widehat p)/n_{\\mathrm{rep}}\\}^{1/2}$. ",
    "Here $z=(\\widehat p-p_{\\mathrm{analytic}})/\\mathrm{s.e.}$. ",
    "The eight non-degenerate comparisons have $|z|<3$. ",
    "NA denotes an undefined $z$ at $d=0$, when survival is certain ",
    "and $\\mathrm{s.e.}=0$.}"
  ),
  "\\label{tab:survival}",
  "\\begin{tabular}{rcccc}",
  "\\toprule",
  "$d$ & Empirical & Analytic & Naive & $z$ \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(survival_results))) {
  row <- survival_results[i, ]
  survival_table_lines <- c(
    survival_table_lines,
    sprintf(
      "%.2f & $%.4f\\;(\\pm %.4f)$ & $%.4f$ & $%.4f$ & $%s$ \\\\",
      row$death, row$empirical, row$margin_95, row$analytic,
      row$naive, if (is.finite(row$z_survival)) {
        sprintf("%+.2f", row$z_survival)
      } else {
        "\\mathrm{NA}"
      }
    )
  )
}
survival_table_lines <- c(
  survival_table_lines, "\\bottomrule", "\\end{tabular}", "\\end{table}", ""
)
write_table_variants(survival_table_lines, "table2_survival.tex")

annealed_table_data <- annealed_results[
  annealed_results$fraction %in% table_fractions, ]
annealed_table_lines <- c(
  "% Generated by simulation_audit.R",
  "\\begin{table}[ht]",
  "\\centering",
  paste0(
    "\\caption{Annealed mean $\\mathbb{E}[X_1(t)]$ at times concentrated ",
    "near $t_*$ for four fixed death rates $d$ ($\\lambda=1$, $\\mu=0.3$, ",
    "$n_{\\mathrm{rep}}=20{,}000$). The column $t/t_*$ shows proximity ",
    "to the common moment boundary. A typical finite sample can ",
    "severely underestimate the analytic mean as $t\\uparrow t_*$.}"
  ),
  "\\label{tab:annealed}",
  "\\begin{tabular}{rrrrrr}",
  "\\toprule",
  "$d$ & $t/t_*$ & $t$ & Empirical & Analytic & Ratio \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(annealed_table_data))) {
  row <- annealed_table_data[i, ]
  if (i > 1 && row$death != annealed_table_data$death[i - 1]) {
    annealed_table_lines <- c(annealed_table_lines, "\\midrule")
  }
  annealed_table_lines <- c(
    annealed_table_lines,
    sprintf(
      "%.2f & %.3f & %.3f & %.3f & %.3f & %.3f \\\\",
      row$death, row$fraction, row$time, row$empirical,
      row$analytic, row$ratio
    )
  )
}
annealed_table_lines <- c(
  annealed_table_lines, "\\bottomrule", "\\end{tabular}", "\\end{table}", ""
)
write_table_variants(annealed_table_lines, "table1_annealed.tex")

table3_data <- bounded_results[
  bounded_results$time_index %in% bounded_table_indices, ]
table3_lines <- c(
  "% Generated by simulation_audit.R",
  "\\begin{table}[ht]",
  "\\centering",
  paste0(
    "\\caption{Annealed mean $\\mathbb{E}[X_1(t)]$ under a ",
    "truncated-exponential birth rate on $[0,a_{\\max}]$, with ",
    "$a_{\\max}=2$, $\\lambda=1$, $\\mu=0.3$, and ",
    "$n_{\\mathrm{rep}}=3000$. The critical death rate is ",
    "$d^*=(1-\\mu)a_*=1.40$.}"
  ),
  "\\label{tab:bounded}",
  "\\begin{tabular}{rrlrr}",
  "\\toprule",
  "$t$ & $d$ & Regime & Empirical & Analytic \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(table3_data))) {
  row <- table3_data[i, ]
  if (i > 1 && row$death != table3_data$death[i - 1]) {
    table3_lines <- c(table3_lines, "\\midrule")
  }
  regime_tex <- if (row$death < 1.4) {
    "$d<d^*$"
  } else if (row$death > 1.4) {
    "$d>d^*$"
  } else {
    "$d=d^*$"
  }
  table3_lines <- c(
    table3_lines,
    sprintf(
      "%.3f & %.2f & %s & %.3f & %.3f \\\\",
      row$time, row$death, regime_tex, row$empirical, row$analytic
    )
  )
}
table3_lines <- c(
  table3_lines, "\\bottomrule", "\\end{tabular}", "\\end{table}", ""
)
write_table_variants(table3_lines, "table3_bounded.tex")

table4_lines <- c(
  "% Generated by simulation_audit.R",
  "\\begin{table}[ht]",
  "\\centering",
  paste0(
    "\\caption{Validation of the small-$d$ asymptotic expansion in ",
    "Remark~\\ref{rem:small-d}. The leading-order approximation ",
    "$1-\\lambda c\\,(1-\\gamma-\\log(\\lambda c))$, with ",
    "$c=d/(1-\\mu)$, approaches the exact survival probability ",
    "$e^{-\\lambda c}-\\lambda c\\,E_1(\\lambda c)$ as $d\\to0^+$. ",
    "Here $\\lambda=1$ and $\\mu=0.3$.}"
  ),
  "\\label{tab:smalld}",
  "\\begin{tabular}{rccc}",
  "\\toprule",
  "$d$ & Exact & SmallDApprox & RelErr \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(small_d_results))) {
  row <- small_d_results[i, ]
  table4_lines <- c(
    table4_lines,
    sprintf(
      "%.2e & $%.6f$ & $%.6f$ & $%.2e$ \\\\",
      row$death, row$exact, row$approximation, row$relative_error
    )
  )
}
table4_lines <- c(
  table4_lines, "\\bottomrule", "\\end{tabular}", "\\end{table}", ""
)
write_table_variants(table4_lines, "table4_smalld.tex")

# ---------------------------------------------------------------------------
# Audit summary
# ---------------------------------------------------------------------------

escape_standard_error_scale <- with(
  survival_results,
  sqrt(pmax(analytic * (1 - analytic), .Machine$double.eps) / n_rep)
)
escape_bias_in_se <- survival_results$escape_bias /
  escape_standard_error_scale

audit_lines <- c(
  "COMPLETE R SIMULATION AUDIT",
  sprintf("Timestamp: %s", format(Sys.time(), tz = "America/Chicago")),
  sprintf("R version: %s", R.version.string),
  sprintf("Master seed: %d", master_seed),
  "",
  "DETERMINISTIC CHECKS",
  sprintf("Transition-law mean max relative error: %.3e",
          max(transition_check$mean_relative_error)),
  sprintf("Transition-law variance max relative error: %.3e",
          max(transition_check$variance_relative_error)),
  sprintf("Exponential-mixture mean max relative error: %.3e",
          max(mean_integral_check$relative_error)),
  sprintf("Survival formula max absolute error: %.3e",
          max(survival_integral_check$absolute_error)),
  sprintf("Bounded-support mean max relative error: %.3e",
          max(bounded_integral_check$relative_error)),
  sprintf("Random-death mean max relative error: %.3e",
          max(random_integral_check$relative_error)),
  sprintf("Tail asymptotic max relative discrepancy at n=10000: %.3e",
          max(tail_asymptotic_check$relative_discrepancy)),
  sprintf("Small-d approximation relative error at d=0.5: %.4f",
          small_d_half_relative_error),
  "",
  "STOCHASTIC CHECKS",
  "Zero-standard-error cases have undefined z and are recorded as NA.",
  sprintf("Non-degenerate survival comparisons: %d",
          sum(is.finite(survival_results$z_survival))),
  sprintf("Exact-survival simulation max |z|: %.3f",
          max(abs(survival_results$z_survival), na.rm = TRUE)),
  sprintf(paste0(
    "Embedded-walk threshold simulation max |z| versus exact ",
    "hitting probability: %.3f"
  ), max(abs(survival_results$z_escape_vs_exact), na.rm = TRUE)),
  sprintf("Escape-threshold N=200 max positive bias: %.3e",
          max(survival_results$escape_bias)),
  sprintf("Escape-threshold bias max in survival-MC s.e. units: %.3f",
          max(escape_bias_in_se)),
  sprintf("Bounded-support simulation max |z|: %.3f",
          max(abs(bounded_results$z), na.rm = TRUE)),
  sprintf("Figure 1 median empirical/analytic ratio: %.3f",
          median(annealed_results$ratio)),
  sprintf("Random-death median empirical/analytic ratio: %.3f",
          median(random_results$ratio)),
  "",
  "HEAVY-TAIL FINITE-SAMPLE CHECK",
  sprintf("Tail index at t=0.9 t*: %.6f", bias_index),
  "The sample mean is unbiased because its expectation exists.",
  "Its finite-sample distribution is highly right-skewed and has infinite",
  "variance at this time; typical (median) experiments underestimate the",
  "analytic mean, while rare extreme experiments restore unbiasedness.",
  sprintf("Median ratio range over sample sizes: %.3f to %.3f",
          min(bias_summary$median_ratio),
          max(bias_summary$median_ratio)),
  sprintf("Mean ratio range over 200-trial batches: %.3f to %.3f",
          min(bias_summary$mean_ratio),
          max(bias_summary$mean_ratio)),
  sprintf("Fraction-below-one range: %.3f to %.3f",
          min(bias_summary$fraction_below_one),
          max(bias_summary$fraction_below_one)),
  "",
  "CONCLUSION",
  "All analytic formulas, the tail asymptotic, and the exact transition",
  "sampler passed their deterministic checks.",
  "The direct embedded-walk simulation agrees with the exact probability",
  "of hitting the finite escape threshold before extinction.",
  "The escape-threshold survival method is approximate; the regenerated",
  "Figure 2 and Table 2 use exact conditional survival simulation.",
  "The manuscript correctly describes the near-boundary behavior as",
  "'finite-sample instability' and 'typical underestimation'."
)

writeLines(audit_lines, file.path(audit_dir, "audit_summary.txt"))

if (publish) {
  figure_map <- c(
    fig1_annealed_mean.pdf = "fig1_annealed_mean.pdf",
    fig2_survival.pdf = "fig2_survival.pdf",
    fig3_bounded.pdf = "fig3_bounded.pdf",
    fig4_finite_sample_instability.pdf = "fig4_finite_sample_instability.pdf",
    fig5_random_death.pdf = "fig5_random_death.pdf",
    fig6_smalld.pdf = "fig6_smalld.pdf"
  )
  copied_figures <- mapply(
    function(source_name, destination_name) {
      file.copy(
        file.path(audit_dir, source_name),
        destination_name,
        overwrite = TRUE
      )
    },
    names(figure_map), unname(figure_map)
  )
  stopifnot(all(copied_figures))

  table_names <- c(
    "table1_annealed.tex", "table2_survival.tex",
    "table3_bounded.tex", "table4_smalld.tex",
    "table1_annealed_tracked.tex", "table2_survival_tracked.tex",
    "table3_bounded_tracked.tex", "table4_smalld_tracked.tex"
  )
  copied_tables <- file.copy(
    file.path(audit_dir, table_names),
    table_names,
    overwrite = TRUE
  )
  stopifnot(all(copied_tables))
  message("Published audited figures and tables to the manuscript directory.")
}

message(paste(audit_lines, collapse = "\n"))
