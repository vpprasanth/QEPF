## =====================================================================
##  qepf_functions.R
##  Core functions for the quantile-based effectiveness persistence
##  function (QEPF),  P(u) = V(u) / Q(u),  V(u) = (1-u)^{-1} int_u^1 Q(p) dp.
##
##  Sankaran, P.G., Prasanth, V.P. and Midhu, N.N.
##  "Quantile-Based Effectiveness Persistence Function: A Tail-Focused
##   Metric with Theory, Estimation, and Application to Biosimilar
##   Evaluation."
##
##  Use:  source("R/qepf_functions.R")
##  The analysis scripts in R/ are self-contained and do not need this
##  file; it collects the methods in one place for use on new data.
## =====================================================================

## ---------------------------------------------------------------------
## 1. Estimation
## ---------------------------------------------------------------------

## k = ceiling(n u), kept in 1..n-1. The small offset avoids the
## floating-point trap ceiling(100 * 0.6) = 61.
kidx <- function(n, u) pmin(pmax(as.integer(ceiling(n * u - 1e-9)), 1L), n - 1L)

## P-hat on a vector of indices k, from a SORTED sample xs (fast path).
Pgrid <- function(xs, k) {
  n  <- length(xs)
  cs <- cumsum(xs)
  (cs[n] - cs[k]) / ((n - k) * xs[k])
}

## Nonparametric estimator, eq. (5.2) of the paper (vectorised in u):
##   P_n(u) = [ (n-k)^{-1} sum_{i=k+1}^n X_(i) ] / X_(k),  k = ceiling(n u)
qepf_hat <- function(x, u) {
  stopifnot(all(x > 0), all(u > 0 & u < 1))
  Pgrid(sort(x), kidx(length(x), u))
}

## Empirical upper-tail mean V_n(u) and quantile Q_n(u)
vitality_hat <- function(x, u) {
  xs <- sort(x); n <- length(xs); k <- kidx(n, u); cs <- cumsum(xs)
  (cs[n] - cs[k]) / (n - k)
}
quantile_hat <- function(x, u) sort(x)[kidx(length(x), u)]

## ---------------------------------------------------------------------
## 2. Closed-form P(u) (parameterizations as in Table 1 of the paper)
## ---------------------------------------------------------------------
P_uniform     <- function(u) (1 + u) / (2 * u)
P_exponential <- function(u) 1 - 1 / log1p(-u)
P_pareto      <- function(u, alpha) rep_len(alpha / (alpha - 1), length(u))   # alpha > 1
P_power       <- function(u, beta) beta * (1 - u^(1 + 1 / beta)) / ((beta + 1) * (1 - u) * u^(1 / beta))
P_weibull     <- function(u, k) {                       # shape k; scale-free
  t <- -log1p(-u); a <- 1 + 1 / k
  gamma(a) * pgamma(t, shape = a, lower.tail = FALSE) / ((1 - u) * t^(1 / k))
}
P_gamma       <- function(u, k) {                       # shape k; scale-free
  t <- qgamma(u, shape = k)                             # Gamma(k+1,t)/Gamma(k) = k * Q(k+1, t)
  k * pgamma(t, shape = k + 1, lower.tail = FALSE) / ((1 - u) * t)
}
P_lmrqd       <- function(u, alpha, mu) {               # Q(u) = -(alpha+mu) log(1-u) - 2 alpha u
  1 + (mu + alpha * u) / (-(alpha + mu) * log1p(-u) - 2 * alpha * u)
}

## ---------------------------------------------------------------------
## 3. Test of equal persistence on U = [u1, u2]  (Section 6.1-6.2)
##    T = sqrt(n1 n2 / (n1 + n2)) * sup_U |P1-hat(u) - P2-hat(u)|
##    calibrated by permutation. normalize = TRUE divides each arm by its
##    own Q-hat(u1) before pooling (scale-normalized permutation).
## ---------------------------------------------------------------------
qepf_perm_test <- function(x, y, U = c(0.60, 0.90), L = 101, B = 2000,
                           normalize = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n1 <- length(x); n2 <- length(y); N <- n1 + n2
  ug <- seq(U[1], U[2], length.out = L)
  k1 <- kidx(n1, ug); k2 <- kidx(n2, ug)
  xs <- sort(x); ys <- sort(y)
  sc <- sqrt(n1 * n2 / N)
  D0 <- Pgrid(xs, k1) - Pgrid(ys, k2)
  T0 <- sc * max(abs(D0))
  z  <- if (normalize) sort(c(xs / xs[k1[1]], ys / ys[k2[1]])) else sort(c(xs, ys))
  Tb <- numeric(B)
  for (b in seq_len(B)) {
    p <- sample.int(N)
    Tb[b] <- sc * max(abs(Pgrid(z[sort.int(p[seq_len(n1)], method = "radix")], k1) -
                          Pgrid(z[sort.int(p[(n1 + 1L):N], method = "radix")], k2)))
  }
  list(T = T0, p.value = (1 + sum(Tb >= T0)) / (B + 1),
       crit95 = unname(quantile(Tb, 0.95)), u = ug, D = D0)
}

## ---------------------------------------------------------------------
## 4. Equivalence test for tail shape  (Section 6.3)
##    H0: sup_U |D(u)| >= Delta  vs  H1: sup_U |D(u)| < Delta,
##    D(u) = P_y(u) - P_x(u). Within-arm bootstrap; equivalence is declared
##    if Delta_min = sup_U |D-hat| + c_{1-alpha} < Delta.
## ---------------------------------------------------------------------
qepf_equiv_test <- function(x, y, Delta, U = c(0.60, 0.90), L = 31, B = 2000,
                            alpha = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n1 <- length(x); n2 <- length(y)
  ug <- seq(U[1], U[2], length.out = L)
  k1 <- kidx(n1, ug); k2 <- kidx(n2, ug)
  xs <- sort(x); ys <- sort(y)
  D0 <- Pgrid(ys, k2) - Pgrid(xs, k1)
  supD <- numeric(B)
  for (b in seq_len(B)) {
    bx <- xs[sort.int(sample.int(n1, n1, replace = TRUE), method = "radix")]
    by <- ys[sort.int(sample.int(n2, n2, replace = TRUE), method = "radix")]
    supD[b] <- max(abs((Pgrid(by, k2) - Pgrid(bx, k1)) - D0))
  }
  c_hat <- unname(quantile(supD, 1 - alpha))
  Delta_min <- max(abs(D0)) + c_hat
  list(sup_D = max(abs(D0)), c = c_hat, Delta_min = Delta_min,
       equivalent = Delta_min < Delta, u = ug, D = D0)
}

## ---------------------------------------------------------------------
## 5. Level component: tail-mean ratio r = V_y(u1) / V_x(u1)  (Section 6.4)
##    Percentile bootstrap CI on the log scale. If cluster ids are given,
##    clusters are resampled within each arm.
## ---------------------------------------------------------------------
qepf_tail_mean_ratio <- function(x, y, u1 = 0.60, B = 2000, level = 0.95,
                                 cluster_x = NULL, cluster_y = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  r0 <- vitality_hat(y, u1) / vitality_hat(x, u1)
  resample <- function(v, cl) {
    if (is.null(cl)) return(v[sample.int(length(v), replace = TRUE)])
    groups <- split(v, cl)
    unlist(groups[sample.int(length(groups), replace = TRUE)], use.names = FALSE)
  }
  rb <- replicate(B, vitality_hat(resample(y, cluster_y), u1) /
                     vitality_hat(resample(x, cluster_x), u1))
  a <- (1 - level) / 2
  ci <- exp(quantile(log(rb), c(a, 1 - a), names = FALSE))
  list(ratio = r0, ci = ci)
}
