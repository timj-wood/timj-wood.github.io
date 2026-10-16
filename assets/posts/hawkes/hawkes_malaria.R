# hawkes_malaria.R
# Simulate and fit a near-elimination malaria outbreak with the epihawkes
# package from Unwin et al. (2021): https://github.com/mrc-ide/epihawkes
#   remotes::install_github("mrc-ide/epihawkes")
library(epihawkes)
library(ggplot2)
set.seed(42)


# ---- 1. Scenario -----------------------------------------------------------
T_max  <- 5 * 365
delay  <- 15          # days with no onward transmission
m      <- 0.10        # mean importation rate (cases / day)
a      <- 0.80        # seasonal amplitude
t_peak <- 200         # day of year importations peak
eta    <- 0.60        # branching ratio (R_c)
sigma  <- 20          # Rayleigh spread (days)

# epihawkes parameterises the same model differently:
#   kernel      alpha * (u - delay) * exp(-delta (u - delay)^2 / 2)   => eta = alpha/delta, sigma^2 = 1/delta
#   background  A + M cos(2 pi t / 365.25) + N sin(2 pi t / 365.25)   => A = m, (M, N) = m a (cos, sin)(w t_peak)
w <- 2 * pi / 365.25
true <- list(alpha = eta / sigma^2, delta = 1 / sigma^2, delay = delay,
             A = m, M = m * a * cos(w * t_peak), N = m * a * sin(w * t_peak))

# ---- 2a. Simulate with epihawkes (Ogata thinning on the full intensity) -----
events_thin <- hawkes_simulation(events = c(0), T_max = T_max, kernel = ray_kernel,
                                 parameters = true, mu_fn = mu_sinusoidal_constant,
                                 mu_fn_diff = mu_diff_sinusoidal_constant, print_level = 0)
cat(sprintf("epihawkes thinning simulation: %d cases\n", length(events_thin)))

# ---- 2b. Simulate as a family tree (records who infected whom) -------------
# Same model, built from the Hawkes & Oakes cluster representation. Uses the
# epihawkes background function; children arrive delay + Rayleigh(sigma) later.
simulate_tree <- function(T_max, p) {
  mu_max <- p$A + sqrt(p$M^2 + p$N^2)
  cand   <- sort(runif(rpois(1, mu_max * T_max), 0, T_max))
  times  <- cand[runif(length(cand)) < mu_sinusoidal_constant(cand, p) / mu_max]
  parent <- rep(0L, length(times))                 # 0 = imported
  frontier <- seq_along(times)
  while (length(frontier) > 0) {
    nxt <- integer(0)
    for (i in frontier) {
      n_child <- rpois(1, p$alpha / p$delta)
      if (n_child == 0) next
      t_child <- times[i] + p$delay + sqrt(-2 * log(runif(n_child)) / p$delta)  # Rayleigh
      t_child <- t_child[t_child < T_max]
      times  <- c(times, t_child); parent <- c(parent, rep(i, length(t_child)))
      nxt <- c(nxt, length(times) - rev(seq_along(t_child)) + 1)
    }
    frontier <- nxt
  }
  o <- order(times); rank <- match(seq_along(times), o)
  list(times = times[o], parent = ifelse(parent[o] > 0, rank[parent[o]], 0L))
}
tree <- simulate_tree(T_max, true)
events   <- tree$times
imported <- tree$parent == 0
cat(sprintf("family-tree simulation:        %d cases (%d imported, %d local, local fraction %.2f)\n",
            length(events), sum(imported), sum(!imported), mean(!imported)))

# ---- 3. Fit with the epihawkes likelihood and analytic gradient -------------
fit_epihawkes <- function(events, n_starts = 10) {
  best <- NULL
  for (k in seq_len(n_starts)) {
    start <- list(alpha = runif(1, 1e-4, 5e-3), delta = runif(1, 1e-3, 1e-2),
                  A = runif(1, 0.02, 0.3), M = runif(1, -0.1, 0.1), N = runif(1, -0.1, 0.1))
    out <- tryCatch(suppressWarnings(optim(
      par = start, fn = neg_log_likelihood, gr = ray_derivatives, method = "BFGS",
      events = events, delay = delay, kernel = ray_kernel,
      mu_fn = mu_sinusoidal_constant, mu_diff_fn = mu_diff_sinusoidal_constant,
      mu_int_fn = mu_int_sinusoidal_constant, print_level = 0,
      control = list(maxit = 500))), error = function(e) NULL)
    if (!is.null(out) && (is.null(best) || out$value < best$value)) best <- out
  }
  p <- as.list(best$par); p$delay <- delay; p$negLL <- best$value; p
}
report <- function(p) {
  data.frame(parameter = c("m", "a", "t_peak", "eta", "sigma"),
             true = c(m, a, t_peak, eta, sigma),
             estimate = c(p$A, sqrt(p$M^2 + p$N^2) / p$A, (atan2(p$N, p$M) / w) %% 365.25,
                          p$alpha / p$delta, 1 / sqrt(p$delta)))
}
fit <- fit_epihawkes(events)
cat("\nMaximum-likelihood fit (epihawkes::neg_log_likelihood, BFGS, best of 10 starts):\n")
print(report(fit), row.names = FALSE, digits = 3)

# ---- 4. Imported or local? ---------------------------------------------------
# P(case i imported) = mu(t_i) / lambda*(t_i); the rest is shared among earlier cases.
branching_probabilities <- function(events, p) {
  lags <- outer(events, events, "-")                        # lags[i, j] = t_i - t_j
  contrib <- ray_kernel(lags, p) * (lags > p$delay)
  contrib[!is.finite(contrib)] <- 0
  mu_i  <- mu_sinusoidal_constant(events, p)
  total <- mu_i + rowSums(contrib)
  list(p_imported = mu_i / total, P = contrib / total)
}
bp <- branching_probabilities(events, fit)
cat(sprintf("\nClassification accuracy at P > 0.5: %.2f\n", mean((bp$p_imported > 0.5) == imported)))
cat(sprintf("Expected number imported: %.1f  (true %d)\n", sum(bp$p_imported), sum(imported)))

# ---- 5. Figures --------------------------------------------------------------
# Fig 2: epihawkes' own intensity plot, with the true labels
p2 <- plot_intensity_imported(events = events, imported_events = events[imported], T_max = T_max,
                              kernel = ray_kernel, parameters = true, mu_fn = mu_sinusoidal_constant) +
  labs(x = "day", y = "cases / day") + theme_minimal()
ggsave("simulated_outbreak_r.png", p2, width = 10, height = 4, dpi = 150)

# Fig 3: true vs fitted kernel
u <- seq(0, 120, by = 0.5)
df3 <- rbind(data.frame(u, phi = ray_kernel(u, true) * (u > delay), which = "true"),
             data.frame(u, phi = ray_kernel(u, fit)  * (u > delay), which = "fitted"))
p3 <- ggplot(df3, aes(u, phi, linetype = which)) + geom_line() +
  labs(x = "days since a case", y = expression(phi(u)), linetype = NULL) + theme_minimal()
ggsave("kernel_r.png", p3, width = 7, height = 3, dpi = 150)

# Fig 4: estimated P(imported), split by truth
df4 <- data.frame(p = bp$p_imported, truth = ifelse(imported, "truly imported", "truly local"))
p4 <- ggplot(df4, aes(p, fill = truth)) + geom_histogram(bins = 20, alpha = 0.6, position = "identity") +
  labs(x = "estimated P(imported)", y = "number of cases", fill = NULL) + theme_minimal()
ggsave("classification_r.png", p4, width = 8, height = 4, dpi = 150)
cat("\nSaved simulated_outbreak_r.png, kernel_r.png, classification_r.png\n")
