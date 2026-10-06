## =====================================================================
##  04_fig3_estimator_curves.R  --  data for Figure 3 (Section 7.1)
##
##  LMRQD(alpha = 0.5, mu = 5), Weibull(k = 2), Pareto I(alpha = 2.5).
##  For each distribution: theoretical P(u), the estimate P_n(u) from one
##  sample of size n = 100, and the Monte Carlo mean with a pointwise 95%
##  band from R = 1000 samples, on u in [0.50, 0.95].
##
##  Output (written to figures/, read by figures/fig3.tex):
##    qepf_<Dist>_Theoretical.csv  (u, P)
##    qepf_<Dist>_Empirical.csv    (u, P)
##    qepf_<Dist>_MCband.csv       (u, mean, lo, hi)
##  plus a quick ggplot preview, figures/fig3_preview.png (if ggplot2 is
##  installed). Then compile figures/fig3.tex from the figures/ folder.
##
##  Run from the repository root:  source("R/04_fig3_estimator_curves.R")
## =====================================================================

## ---- find the repository root (works from QEPF.Rproj, from the root, or
## ---- with source("<full path>/R/<script>.R")) ----------------------------
if (!file.exists("R/qepf_functions.R")) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
if (!dir.exists("figures")) stop("Please set the working directory to the repository root (the folder containing R/ and figures/).")

set.seed(20251015)
n      <- 100                         # sample size of the single displayed sample
R      <- 1000                        # Monte Carlo samples for the band
u_grid <- seq(0.50, 0.95, by = 0.01)  # upper half; k = ceiling(n*u) < n for n >= 21

## ---- estimator, eq. (5.2), vectorised over u ------------------------------
Phat <- function(x, u) {
  xs <- sort.int(x); n <- length(xs)
  k  <- pmin(pmax(as.integer(ceiling(n * u - 1e-9)), 1L), n - 1L)
  cs <- cumsum(xs)
  (cs[n] - cs[k]) / ((n - k) * xs[k])
}

## ---- distributions: random generator + closed-form P(u) --------------------
dists <- list(
  LMRQD = list(
    label = "LMRQD (alpha = 0.5, mu = 5)",
    r = function(n) { U <- runif(n); -(0.5 + 5) * log1p(-U) - 2 * 0.5 * U },
    P = function(u) 1 + (5 + 0.5 * u) / (-(0.5 + 5) * log1p(-u) - 2 * 0.5 * u)
  ),
  Weibull = list(
    label = "Weibull (k = 2)",
    r = function(n) rweibull(n, shape = 2, scale = 1),
    P = function(u) {                      # Gamma(a, x) / [(1-u) x^(1/k)]
      k <- 2; x <- -log1p(-u); a <- 1 + 1 / k
      gamma(a) * pgamma(x, shape = a, lower.tail = FALSE) / ((1 - u) * x^(1 / k))
    }
  ),
  Pareto_I = list(
    label = "Pareto I (alpha = 2.5)",
    r = function(n) runif(n)^(-1 / 2.5),   # sigma = 1
    P = function(u) rep(2.5 / 1.5, length(u))
  )
)

## ---- compute and write CSVs (unquoted headers, as pgfplots requires) --------
summ <- list()
for (nm in names(dists)) {
  d     <- dists[[nm]]
  theo  <- d$P(u_grid)
  emp   <- Phat(d$r(n), u_grid)
  M     <- t(replicate(R, Phat(d$r(n), u_grid)))          # R x length(u_grid)
  band  <- data.frame(u    = u_grid,
                      mean = colMeans(M),
                      lo   = apply(M, 2, quantile, 0.025),
                      hi   = apply(M, 2, quantile, 0.975))
  write.csv(data.frame(u = u_grid, P = theo), paste0("figures/qepf_", nm, "_Theoretical.csv"), row.names = FALSE, quote = FALSE)
  write.csv(data.frame(u = u_grid, P = emp),  paste0("figures/qepf_", nm, "_Empirical.csv"),   row.names = FALSE, quote = FALSE)
  write.csv(band,                             paste0("figures/qepf_", nm, "_MCband.csv"),      row.names = FALSE, quote = FALSE)
  summ[[nm]] <- data.frame(Distribution = d$label,
                           max_abs_bias_MCmean = round(max(abs(band$mean - theo)), 4),
                           band_hi_max = round(max(band$hi), 3))
}
print(do.call(rbind, summ), row.names = FALSE)

## ---- quick preview with ggplot2 (optional) -----------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  pdat <- do.call(rbind, lapply(names(dists), function(nm) {
    th <- read.csv(paste0("figures/qepf_", nm, "_Theoretical.csv"))
    em <- read.csv(paste0("figures/qepf_", nm, "_Empirical.csv"))
    bd <- read.csv(paste0("figures/qepf_", nm, "_MCband.csv"))
    data.frame(Distribution = dists[[nm]]$label, u = th$u,
               Theoretical = th$P, Empirical = em$P, MCmean = bd$mean, lo = bd$lo, hi = bd$hi)
  }))
  pdat$Distribution <- factor(pdat$Distribution, levels = sapply(dists, `[[`, "label"))
  pp <- ggplot(pdat, aes(u)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#56B4E9", alpha = 0.25) +
    geom_line(aes(y = Theoretical), colour = "#009E73", linewidth = 1) +
    geom_line(aes(y = MCmean), colour = "#0072B2", linetype = "dotted", linewidth = 0.9) +
    geom_line(aes(y = Empirical), colour = "#D55E00", linetype = "dashed", linewidth = 0.8) +
    facet_wrap(~ Distribution, scales = "free_y", nrow = 1) +
    labs(x = "u", y = "P(u)") + theme_bw(base_size = 11)
  ggsave(file.path("figures", "fig3_preview.png"), pp, width = 10, height = 3.6, dpi = 200)
  print(pp)
}
