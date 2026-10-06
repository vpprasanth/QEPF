## =====================================================================
##  08_asymptotic_variance_check.R  --  numerical check of Proposition 5.1
##
##  Compares the asymptotic variance sigma^2_P(u) of Proposition 5.1 with
##  n * Var{P_n(u)} from Monte Carlo samples, for the Weibull distribution
##  with k = 2 (u = 0.9, 0.8) and the exponential distribution (k = 1, u = 0.9).
##
##    sigma^2_P(u) = [ sigma_V^2 - 2 P(u) sigma_VQ + P(u)^2 sigma_Q^2 ] / Q(u)^2
##    sigma_Q^2  = u (1-u) q(u)^2,   sigma_VQ = u q(u) M(u),
##    sigma_V^2  = Var{(X - t)_+} / (1-u)^2,   t = Q(u),  q(u) = Q'(u)
##
##  Run from the repository root:  source("R/08_asymptotic_variance_check.R")
##  Monte Carlo design: R = 50000 samples of size n = 10000.
##  Output: results/asymptotic_variance_check.csv     Run time: about 1.5-2 hours
## =====================================================================

## ---- find the repository root -------------------------------------------
if (!file.exists("R/qepf_functions.R")) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
if (!dir.exists("results")) stop("Please set the working directory to the repository root (the folder containing R/ and results/).")

## ---- Weibull(shape k, scale 1): Q(p) = s^(1/k) with s = -log(1-p) ----------
## Integrals over p in (u, 1) are computed on the s-scale, dp = e^{-s} ds,
## which removes the singularity at p = 1.
sigma2_P_weibull <- function(u, k) {
  s0 <- -log1p(-u)
  t  <- s0^(1 / k)                                        # Q(u)
  q  <- (1 / k) * s0^(1 / k - 1) / (1 - u)                # Q'(u)
  I1 <- integrate(function(s) s^(1 / k) * exp(-s), s0, Inf, rel.tol = 1e-10)$value
  I2 <- integrate(function(s) (s^(1 / k) - t)^2 * exp(-s), s0, Inf, rel.tol = 1e-10)$value
  V  <- I1 / (1 - u); M <- V - t; P <- V / t
  var_plus <- I2 - ((1 - u) * M)^2                        # Var{(X - t)_+}
  sV  <- var_plus / (1 - u)^2
  sVQ <- u * q * M
  sQ  <- u * (1 - u) * q^2
  c(P = P, sigma2 = (sV - 2 * P * sVQ + P^2 * sQ) / t^2)
}

## ---- estimator, eq. (5.2) ----------------------------------------------------
Phat <- function(x, u) {
  xs <- sort(x); n <- length(xs); k <- as.integer(ceiling(n * u - 1e-9))
  mean(xs[(k + 1):n]) / xs[k]
}

## ---- Monte Carlo -------------------------------------------------------------
set.seed(20261006)
n <- 10000; R <- 50000
cases <- data.frame(dist = c("Weibull (k = 2)", "Weibull (k = 2)", "Exponential (k = 1)"),
                    k = c(2, 2, 1), u = c(0.9, 0.8, 0.9))

res <- do.call(rbind, lapply(seq_len(nrow(cases)), function(i) {
  k <- cases$k[i]; u <- cases$u[i]
  th <- sigma2_P_weibull(u, k)
  est <- replicate(R, Phat(rweibull(n, shape = k, scale = 1), u))
  nv  <- n * var(est)
  data.frame(dist = cases$dist[i], u = u, P_true = round(th[["P"]], 4),
             sigma2_formula = round(th[["sigma2"]], 3),
             n_var_MC = round(nv, 3),
             MC_se_of_n_var = round(nv * sqrt(2 / (R - 1)), 3),   # approx. SE of the MC variance
             sqrt_n_bias = round(sqrt(n) * (mean(est) - th[["P"]]), 3))
}))

print(res, row.names = FALSE)
write.csv(res, "results/asymptotic_variance_check.csv", row.names = FALSE)

## The formula and the Monte Carlo values should agree to within about two
## Monte Carlo standard errors (column MC_se_of_n_var).
