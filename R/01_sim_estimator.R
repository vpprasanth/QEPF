## =====================================================================
##  01_sim_estimator.R  --  Table 2: bias and MSE of the estimator P_n(u)
##
##  LMRQD(alpha = 0.5, mu = 5), Weibull(k = 2), Pareto I(alpha = 2.5);
##  u in {0.85, 0.90, 0.95}; n in {25, 50, 100, 1000}; 1000 replications.
##  (The exponential distribution is also computed; it is not in Table 2.)
##
##  Run from the repository root:  source("R/01_sim_estimator.R")
##  Output: results/sim_estimator_bias_mse.csv      Run time: a few minutes
## =====================================================================

## ---- find the repository root (works from QEPF.Rproj, from the root, or
## ---- with source("<full path>/R/<script>.R")) ----------------------------
if (!file.exists("R/qepf_functions.R")) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
if (!dir.exists("results")) stop("Please set the working directory to the repository root (the folder containing R/ and results/).")
set.seed(12345)

## --- True QEPF -----------------------------------------------------------
qepf_true_exp <- function(u) 1 + 1 / (-log(1 - u))
qepf_true_weibull <- function(u, k = 2, scale = 1) {
  a <- -log(1 - u); s <- 1 + 1 / k
  Gamma_upper <- gamma(s) * pgamma(a, shape = s, lower.tail = FALSE)
  (exp(a) * Gamma_upper) / (a^(1 / k))
}
qepf_true_pareto <- function(u, alpha = 2.5, xm = 1) alpha / (alpha - 1)
qepf_true_lmrqd  <- function(u, alpha = 0.5, mu = 5) {
  denom <- -(alpha + mu) * log(1 - u) - 2 * alpha * u
  1 + (mu + alpha * u) / denom
}

## --- Estimator, eq. (5.2) ----------------------------------------------------
qepf_hat <- function(x, u) {
  x <- sort(x); n <- length(x); k <- ceiling(n * u)
  if (k >= n) return(NA_real_)
  mean(x[(k + 1):n]) / x[k]
}

## --- Random number generators (quantile inversion) ---------------------------
rpareto <- function(n, alpha = 2.5, xm = 1) { u <- runif(n); xm * (1 - u)^(-1 / alpha) }
rlmrqd  <- function(n, alpha = 0.5, mu = 5) { u <- runif(n); -(alpha + mu) * log(1 - u) - 2 * alpha * u }

## --- One design cell ---------------------------------------------------------
run_cell <- function(dist, n, u, reps = 1000, pars) {
  hats <- numeric(reps)
  for (r in seq_len(reps)) {
    x <- switch(dist,
                exp     = rexp(n, rate = pars$rate),
                weibull = rweibull(n, shape = pars$shape, scale = 1),
                pareto  = rpareto(n, alpha = pars$alpha, xm = pars$xm),
                lmrqd   = rlmrqd(n, alpha = pars$alpha_lmrqd, mu = pars$mu_lmrqd))
    hats[r] <- qepf_hat(x, u)
  }
  hats <- hats[is.finite(hats)]
  true <- switch(dist,
                 exp     = qepf_true_exp(u),
                 weibull = qepf_true_weibull(u, k = pars$shape),
                 pareto  = qepf_true_pareto(u, alpha = pars$alpha, xm = pars$xm),
                 lmrqd   = qepf_true_lmrqd(u, alpha = pars$alpha_lmrqd, mu = pars$mu_lmrqd))
  data.frame(dist = dist, u = u, n = n, true = true, mean = mean(hats),
             bias = mean(hats) - true, mse = mean((hats - true)^2), reps_eff = length(hats))
}

## --- Grid and run --------------------------------------------------------------
ns    <- c(25, 50, 100, 1000)
us    <- c(0.85, 0.90, 0.95)
dists <- c("exp", "weibull", "pareto", "lmrqd")
pars  <- list(shape = 2, rate = 1, alpha = 2.5, xm = 1, alpha_lmrqd = 0.5, mu_lmrqd = 5)

res <- do.call(rbind, lapply(dists, function(d)
  do.call(rbind, lapply(us, function(u)
    do.call(rbind, lapply(ns, function(n) run_cell(d, n, u, reps = 1000, pars = pars)))))))

res_print <- transform(res, true = round(true, 4), bias = round(bias, 4), mse = round(mse, 4))
print(res_print[, c("dist", "u", "n", "true", "bias", "mse")], row.names = FALSE)
write.csv(res, "results/sim_estimator_bias_mse.csv", row.names = FALSE)
