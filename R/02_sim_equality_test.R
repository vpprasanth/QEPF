## =====================================================================
##  02_sim_equality_test.R  --  Tables 3 and 4: size and power of the
##  test of equal persistence under several calibrations
##
##  Run from the repository root:  source("R/02_sim_equality_test.R")
##  Output: results/sim_equality_size.csv, results/sim_equality_power.csv
##  Run time: several hours for the full run (QUICK <- FALSE) on a laptop;
##  set QUICK <- TRUE for a check that takes about a minute.
##
##  Calibrations, all computed from the same simulated data sets:
##    boot      : pooled bootstrap,  sup |D|  (calibration of the first submission)
##    perm      : permutation,       sup |D|  (default in the revised paper)
##    perm_norm : permutation after dividing each arm by its own Q-hat(u1)
##                (scale-normalized permutation)
##    boot_stud : pooled bootstrap,  sup |D(u)| / sd*(u)  (studentized)
##    perm_stud : permutation,       sup |D(u)| / sd*(u)
##    boot_mn   : m-out-of-n pooled bootstrap (power setting of the first
##                submission; kept for comparison)
##  Scenario gamma vs 0.7*gamma has equal persistence functions but a 30%
##  lower level (Proposition 3.3 / Remark on scale blindness).
## =====================================================================

suppressPackageStartupMessages(library(parallel))

## ---- find the repository root (works from QEPF.Rproj, from the root, or
## ---- with source("<full path>/R/<script>.R")) ----------------------------
if (!file.exists("R/qepf_functions.R")) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
if (!dir.exists("results")) stop("Please set the working directory to the repository root (the folder containing R/ and results/).")
out_dir <- "results"

## ---- 1. Data generators ---------------------------------------------------
rpareto <- function(n, alpha = 2.5, xm = 1) xm * runif(n)^(-1 / alpha)
rLMRQD  <- function(n, alpha = 3.5, mu = 5) {
  U <- runif(n); -(alpha + mu) * log1p(-U) - 2 * alpha * U
}
gen <- list(
  gamma     = function(n) rgamma(n, shape = 4, scale = 1),
  LMRQD     = function(n) rLMRQD(n, alpha = 3.5, mu = 5),
  lognormal = function(n) rlnorm(n, meanlog = 0, sdlog = 2),
  weibull   = function(n) rweibull(n, shape = 5, scale = 20),
  pareto    = function(n) rpareto(n, alpha = 2.5, xm = 1),
  gamma_x07 = function(n) 0.7 * rgamma(n, shape = 4, scale = 1)  # uniform 30% lower level
)

## ---- 2. Grid, k-index and vectorized estimator -----------------------------
U_RANGE <- c(0.60, 0.90)
L_GRID  <- 101

ugrid_for <- function(n1, n2 = n1, U = U_RANGE, L = L_GRID, eps = 1e-8) {
  u_max <- min(U[2], min((n1 - 1) / n1, (n2 - 1) / n2) - eps)
  seq(max(U[1], eps), u_max, length.out = L)
}
kidx <- function(n, u) pmin(pmax(as.integer(ceiling(n * u - 1e-9)), 1L), n - 1L)

## P-hat(u) = [ sum_{i>k} X_(i) / (n-k) ] / X_(k), from a SORTED sample
Pgrid <- function(xs, k) {
  n  <- length(xs)
  cs <- cumsum(xs)
  (cs[n] - cs[k]) / ((n - k) * xs[k])
}

rowMaxAbs <- function(M) {
  r <- abs(M[, 1])
  for (j in 2:ncol(M)) r <- pmax(r, abs(M[, j]))
  r
}
colSD <- function(M) {
  mu <- colMeans(M)
  sqrt(colSums((M - rep(mu, each = nrow(M)))^2) / (nrow(M) - 1))
}

## ---- 3. One data set -> p-values under each calibration --------------------
qepf_test_all <- function(x, y, ug, B = 499, m_frac = 0.85) {
  n1 <- length(x); n2 <- length(y); N <- n1 + n2
  k1 <- kidx(n1, ug); k2 <- kidx(n2, ug)
  sc <- sqrt(n1 * n2 / N)
  D0 <- Pgrid(sort.int(x), k1) - Pgrid(sort.int(y), k2)
  T0 <- sc * max(abs(D0))
  z  <- sort.int(c(x, y))                      # pooled sample, sorted once

  m1 <- floor(m_frac * n1); m2 <- floor(m_frac * n2)
  km1 <- kidx(m1, ug); km2 <- kidx(m2, ug)
  scm <- sqrt(m1 * m2 / (m1 + m2))

  ## scale-normalized pooled sample: each arm divided by its own Q-hat(u1)
  xs <- sort.int(x); ys <- sort.int(y)
  zn <- sort.int(c(xs / xs[k1[1]], ys / ys[k2[1]]))

  G  <- length(ug)
  Db <- matrix(0, B, G); Dp <- matrix(0, B, G); Dn <- matrix(0, B, G); Tmn <- numeric(B)
  for (b in seq_len(B)) {
    ## pooled bootstrap (with replacement from the pooled sample)
    i1 <- sort.int(sample.int(N, n1, replace = TRUE), method = "radix")
    i2 <- sort.int(sample.int(N, n2, replace = TRUE), method = "radix")
    Db[b, ] <- Pgrid(z[i1], k1) - Pgrid(z[i2], k2)
    ## permutation (split the pooled sample without replacement)
    p <- sample.int(N)
    j1p <- sort.int(p[seq_len(n1)], method = "radix")
    j2p <- sort.int(p[(n1 + 1L):N], method = "radix")
    Dp[b, ] <- Pgrid(z[j1p], k1) - Pgrid(z[j2p], k2)
    ## scale-normalized permutation (same split, normalized pooled sample)
    Dn[b, ] <- Pgrid(zn[j1p], k1) - Pgrid(zn[j2p], k2)
    ## m-out-of-n pooled bootstrap
    j1 <- sort.int(sample.int(N, m1, replace = TRUE), method = "radix")
    j2 <- sort.int(sample.int(N, m2, replace = TRUE), method = "radix")
    Tmn[b] <- scm * max(abs(Pgrid(z[j1], km1) - Pgrid(z[j2], km2)))
  }

  p_raw  <- function(Dstar) (1 + sum(sc * rowMaxAbs(Dstar) >= T0)) / (B + 1)
  p_stud <- function(Dstar) {
    s <- colSD(Dstar)
    (1 + sum(rowMaxAbs(sweep(Dstar, 2, s, "/")) >= max(abs(D0) / s))) / (B + 1)
  }
  c(boot = p_raw(Db), perm = p_raw(Dp), perm_norm = p_raw(Dn),
    boot_stud = p_stud(Db), perm_stud = p_stud(Dp),
    boot_mn = (1 + sum(Tmn >= T0)) / (B + 1))
}

## ---- 4. Parallel driver (works on Windows, macOS and Linux) ---------------
sim_one <- function(r, d1, d2, n, B, m_frac) {
  qepf_test_all(gen[[d1]](n), gen[[d2]](n), ug = ugrid_for(n), B = B, m_frac = m_frac)
}

run_grid <- function(scen, ns, R = 2000, B = 499, alpha = 0.05,
                     m_frac = 0.85, cores = max(1, detectCores() - 1), seed = 20251123) {
  cl <- makeCluster(cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c("gen", "rpareto", "rLMRQD", "ugrid_for", "kidx", "Pgrid",
                      "rowMaxAbs", "colSD", "qepf_test_all", "sim_one", "U_RANGE", "L_GRID"),
                envir = globalenv())
  clusterSetRNGStream(cl, seed)

  res <- list()
  for (s in seq_len(nrow(scen))) for (n in ns) {
    t0 <- proc.time()[["elapsed"]]
    pv <- parSapply(cl, seq_len(R), sim_one,
                    d1 = scen$d1[s], d2 = scen$d2[s], n = n, B = B, m_frac = m_frac)
    rate <- rowMeans(pv <= alpha)
    row  <- data.frame(scenario = paste0(scen$d1[s], "_vs_", scen$d2[s]),
                       truth = scen$truth[s], n = n,
                       t(round(100 * rate, 2)),
                       mc_se = round(100 * sqrt(alpha * (1 - alpha) / R), 2),
                       secs = round(proc.time()[["elapsed"]] - t0, 1))
    print(row, row.names = FALSE)
    res[[length(res) + 1]] <- row
  }
  do.call(rbind, res)
}

## =====================================================================
##  5. RUN
## =====================================================================
QUICK <- FALSE      # FALSE: full run used in the paper; TRUE: quick check

if (QUICK) {
  R_rep <- 100; B_rep <- 199; ns_size <- c(100, 1000); ns_pow <- c(100)
} else {
  R_rep <- 2000; B_rep <- 499
  ns_size <- c(50, 100, 300, 500, 1000, 2000)
  ns_pow  <- c(50, 100, 300, 500, 1000)
}

size_scen <- data.frame(
  d1 = c("gamma", "LMRQD", "lognormal", "gamma"),
  d2 = c("gamma", "LMRQD", "lognormal", "gamma_x07"),
  truth = c("H0", "H0", "H0", "H0: equal P, 30% lower level")
)
pow_scen <- data.frame(
  d1 = c("gamma", "gamma", "LMRQD"),
  d2 = c("LMRQD", "lognormal", "lognormal"),
  truth = "H1"
)

size_res <- run_grid(size_scen, ns_size, R = R_rep, B = B_rep)
pow_res  <- run_grid(pow_scen,  ns_pow,  R = R_rep, B = B_rep)

write.csv(size_res, file.path(out_dir, "sim_equality_size.csv"),  row.names = FALSE)
write.csv(pow_res,  file.path(out_dir, "sim_equality_power.csv"), row.names = FALSE)

## Reading the results: with R = 2000 the Monte Carlo SE of a 5% rate is
## about 0.49 percentage points, so rates of about 4.0-6.0% match the
## nominal level.
