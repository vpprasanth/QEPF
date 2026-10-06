## =====================================================================
##  06_smarter_level_component.R  --  Section 8: level component and
##  Delta_min for the SMARTER trial
##
##  1. Upper-tail mean ratio r = V_int(0.6) / V_ctl(0.6) with a 95%
##     village-level bootstrap CI (villages resampled within each arm).
##  2. Entry thresholds Q(0.6) in each arm.
##  3. Delta_min = sup_U |D-hat(u)| + c_0.95 for the shape comparison on
##     U = [0.60, 0.90], with c from a village-level bootstrap within arms.
##
##  These are the values reported in the manuscript (r = 1.161, 95% CI
##  1.077-1.255; Delta_min = 0.261), computed with seed 2025 and B = 2000.
##
##  Run from the repository root:  source("R/06_smarter_level_component.R")
##  Output: results/smarter_level_component.csv     Run time: a few minutes
## =====================================================================

## ---- find the repository root (works from QEPF.Rproj, from the root, or
## ---- with source("<full path>/R/<script>.R")) ----------------------------
if (!file.exists("R/qepf_functions.R")) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
if (!file.exists("R/qepf_functions.R")) stop("Please set the working directory to the repository root (the folder containing R/, data/ and results/).")
source("R/qepf_functions.R")

dat <- read.csv("data/smarter_anonymised_data.csv")
dat$Y <- dat$SBP0 - dat$SBP4
dat <- dat[is.finite(dat$Y) & dat$Y > 0, ]
int <- dat[dat$group == 1, ]   # intervention
ctl <- dat[dat$group == 2, ]   # control
cat(sprintf("Positive reductions: intervention n = %d (%d villages), control n = %d (%d villages)\n",
            nrow(int), length(unique(int$village_id)), nrow(ctl), length(unique(ctl$village_id))))

u1 <- 0.60
## 1-2. level component (y = intervention, x = control)
lev <- qepf_tail_mean_ratio(x = ctl$Y, y = int$Y, u1 = u1, B = 2000,
                            cluster_x = ctl$village_id, cluster_y = int$village_id, seed = 2025)
cat(sprintf("V(0.6): intervention %.2f, control %.2f\n", vitality_hat(int$Y, u1), vitality_hat(ctl$Y, u1)))
cat(sprintf("Q(0.6): intervention %.1f, control %.1f\n", quantile_hat(int$Y, u1), quantile_hat(ctl$Y, u1)))
cat(sprintf("Tail-mean ratio r = %.3f (95%% CI %.3f-%.3f)\n", lev$ratio, lev$ci[1], lev$ci[2]))

## 3. Delta_min with a village-level bootstrap within arms
set.seed(2025)
ug <- seq(0.60, 0.90, length.out = 30)
P_int <- qepf_hat(int$Y, ug); P_ctl <- qepf_hat(ctl$Y, ug)
D0 <- P_int - P_ctl
res_int <- function() { g <- split(int$Y, int$village_id); unlist(g[sample.int(length(g), replace = TRUE)]) }
res_ctl <- function() { g <- split(ctl$Y, ctl$village_id); unlist(g[sample.int(length(g), replace = TRUE)]) }
supdev <- replicate(2000, max(abs((qepf_hat(res_int(), ug) - qepf_hat(res_ctl(), ug)) - D0)))
c95 <- unname(quantile(supdev, 0.95))
cat(sprintf("sup|D| = %.3f, c_0.95 = %.3f, Delta_min = %.3f\n", max(abs(D0)), c95, max(abs(D0)) + c95))

write.csv(data.frame(
  quantity = c("V_int(0.6)", "V_ctl(0.6)", "Q_int(0.6)", "Q_ctl(0.6)", "r", "r_lo95", "r_hi95",
               "P_int(0.6)", "P_ctl(0.6)", "sup_abs_D", "c95", "Delta_min"),
  value = c(vitality_hat(int$Y, u1), vitality_hat(ctl$Y, u1), quantile_hat(int$Y, u1), quantile_hat(ctl$Y, u1),
            lev$ratio, lev$ci, P_int[1], P_ctl[1], max(abs(D0)), c95, max(abs(D0)) + c95)),
  "results/smarter_level_component.csv", row.names = FALSE)
