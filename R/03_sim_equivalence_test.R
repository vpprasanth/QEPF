## =====================================================================
##  03_sim_equivalence_test.R  --  Table 5: equivalence test for tail shape
##
##  H0: sup_U |P_bio(u) - P_ref(u)| >= Delta  vs  H1: sup_U |D(u)| < Delta
##  Decision: within-arm bootstrap (not pooled). Let c be the (1-alpha)
##  quantile of sup_U |D*(u) - D(u)|. Declare equivalence if
##  sup_U |D(u)| + c < Delta (simultaneous band inside (-Delta, Delta)).
##  Also reports the relative margin on the excess persistence,
##  log{(P_bio-1)/(P_ref-1)} within +/- log(theta).
##
##  Scenarios
##   power    : Gamma(4) vs Gamma(4)
##   boundary : Gamma(4) vs Gamma(2.461), for which sup_U |D(u)| = 0.15
##              (found numerically), so with Delta = 0.15 the rate of
##              declaring equivalence is the type I error at the boundary.
##
##  Run from the repository root:  source("R/03_sim_equivalence_test.R")
##  Output: results/sim_equivalence.csv     Run time: about 15 minutes
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

ugrid <- seq(0.60, 0.90, length.out = 31)
kidx  <- function(n, u = ugrid) pmin(pmax(as.integer(ceiling(n * u - 1e-9)), 1L), n - 1L)
Pgrid <- function(xs, k) { n <- length(xs); cs <- cumsum(xs); (cs[n] - cs[k]) / ((n - k) * xs[k]) }

equiv_one <- function(x, y, B = 499, alpha = 0.05,
                      Deltas = c(0.10, 0.15, 0.20), thetas = c(1.25, 1.5)) {
  n1 <- length(x); n2 <- length(y); k1 <- kidx(n1); k2 <- kidx(n2)
  xs <- sort.int(x); ys <- sort.int(y)
  px <- Pgrid(xs, k1); py <- Pgrid(ys, k2)
  D0 <- py - px; L0 <- log((py - 1) / (px - 1))
  supD <- supL <- numeric(B)
  for (b in seq_len(B)) {
    bx <- xs[sort.int(sample.int(n1, n1, replace = TRUE), method = "radix")]
    by <- ys[sort.int(sample.int(n2, n2, replace = TRUE), method = "radix")]
    qx <- Pgrid(bx, k1); qy <- Pgrid(by, k2)
    supD[b] <- max(abs((qy - qx) - D0))
    supL[b] <- max(abs(log((qy - 1) / (qx - 1)) - L0))
  }
  cD <- quantile(supD, 1 - alpha, names = FALSE)
  cL <- quantile(supL, 1 - alpha, names = FALSE)
  c(setNames(max(abs(D0)) + cD < Deltas, paste0("abs_", Deltas)),
    setNames(max(abs(L0)) + cL < log(thetas), paste0("rel_", thetas)),
    Delta_min = max(abs(D0)) + cD)
}

gen <- list(g4 = function(n) rgamma(n, shape = 4), g2461 = function(n) rgamma(n, shape = 2.461))
sim_one <- function(r, d1, d2, n, B) equiv_one(gen[[d1]](n), gen[[d2]](n), B = B)

run_equiv <- function(scen, ns, R = 2000, B = 499, cores = max(1, detectCores() - 1), seed = 2026) {
  cl <- makeCluster(cores); on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c("ugrid", "kidx", "Pgrid", "equiv_one", "gen", "sim_one"), envir = globalenv())
  clusterSetRNGStream(cl, seed)
  res <- list()
  for (s in seq_len(nrow(scen))) for (n in ns) {
    t0 <- proc.time()[["elapsed"]]
    m  <- parSapply(cl, seq_len(R), sim_one, d1 = scen$d1[s], d2 = scen$d2[s], n = n, B = B)
    dec <- m[rownames(m) != "Delta_min", , drop = FALSE]
    row <- data.frame(scenario = scen$label[s], n = n, t(round(100 * rowMeans(dec), 2)),
                      median_Delta_min = round(median(m["Delta_min", ]), 3),
                      secs = round(proc.time()[["elapsed"]] - t0, 1))
    print(row, row.names = FALSE); res[[length(res) + 1]] <- row
  }
  do.call(rbind, res)
}

QUICK <- FALSE   # FALSE: full run used in the paper; TRUE: quick check
scen <- data.frame(d1 = c("g4", "g4"), d2 = c("g4", "g2461"),
                   label = c("power: Gamma(4) vs Gamma(4)",
                             "boundary: sup|D| = 0.15 (read column abs_0.15)"))
if (QUICK) {
  out <- run_equiv(scen, ns = c(500, 1000), R = 100, B = 199)
} else {
  out <- run_equiv(scen, ns = c(300, 500, 1000, 2000, 5000), R = 2000, B = 499)
}
write.csv(out, file.path(out_dir, "sim_equivalence.csv"), row.names = FALSE)
