## =====================================================================
##  05_smarter_analysis.R  --  Section 8: SMARTER trial
##  Figure 4 (P-hat curves), test of equal persistence with village-level
##  permutation and bootstrap calibration, and Table 6 (sensitivity to the
##  tail interval).
##
##  Data: data/smarter_anonymised_data.csv
##        (Dryad, doi:10.5061/dryad.tmpg4f58w; group 1 = intervention, 2 = control)
##  Outcome: SBP reduction at 12 months, SBP0 - SBP4, restricted to
##           participants with a positive reduction.
##
##  Run from the repository root:  source("R/05_smarter_analysis.R")
##  Output: figures/fig4_smarter_curves.pdf/.png, figures/smarter_*.png,
##          results/smarter_T_perm.csv, results/smarter_T_boot.csv,
##          results/smarter_sensitivity.csv
##  Run time: the baseline test (B = 20000) and the sensitivity analysis
##  (16 intervals, B = 10000 each) take several hours; reduce B_perm,
##  B_boot and B_sens for a quick check.
## =====================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(latex2exp)
})

## ---------------------------
## Paths
## ---------------------------
## find the repository root (works from QEPF.Rproj, from the root, or with
## source("<full path>/R/<script>.R"))
if (!file.exists("R/qepf_functions.R")) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) setwd(dirname(dirname(this_file)))
}
if (!dir.exists("data")) stop("Please set the working directory to the repository root (the folder containing R/ and data/).")
in_csv      <- "data/smarter_anonymised_data.csv"
out_dir_fig <- "figures"
out_dir_tab <- "results"
stopifnot(file.exists(in_csv))

## ---------------------------
## Load & basic prep
## ---------------------------
dat <- readr::read_csv(in_csv, show_col_types = FALSE)
needed <- c("ID","village_id","group","SBP0","SBP4","risk0","risk4")
if (!all(needed %in% names(dat))) {
  stop("Some required columns are missing. Check: ", paste(setdiff(needed, names(dat)), collapse = ", "))
}

df0 <- dat %>%
  mutate(
    arm = factor(group, levels = c(1, 2), labels = c("Intervention","Control")),
    village_id = as.integer(village_id),
    SBP0 = as.numeric(SBP0), SBP4 = as.numeric(SBP4),
    risk0 = as.numeric(risk0), risk4 = as.numeric(risk4)
  ) %>%
  filter(!is.na(arm), !is.na(village_id))

## ---------------------------
## Analysis settings
## ---------------------------
analysis_outcome <- "SBP"      # "SBP" or "ASCVD"
seed    <- 2025
L_grid  <- 30                  # grid points on U
B_perm  <- 20000               # village-level permutations (baseline test)
B_boot  <- 20000               # village-level pooled bootstrap (baseline test)
B_sens  <- 10000               # per interval in the sensitivity analysis (Table 6)
RUN_SENSITIVITY <- TRUE

## ---------------------------
## Outcome: positive reductions only
## ---------------------------
df <- df0 %>%
  mutate(
    base = if (analysis_outcome == "SBP") SBP0 else risk0,
    foll = if (analysis_outcome == "SBP") SBP4 else risk4,
    Y    = pmax(0, base - foll)
  ) %>%
  filter(is.finite(Y), is.finite(base), is.finite(foll)) %>%
  filter(Y > 0)
benefit_var <- "Y"

n_int <- sum(df$arm == "Intervention"); n_ctl <- sum(df$arm == "Control")
cat(sprintf("Participants with a positive reduction: Intervention = %d, Control = %d\n", n_int, n_ctl))
cl_tab <- df %>% distinct(village_id, arm)
cat(sprintf("Villages: Intervention = %d, Control = %d\n",
            sum(cl_tab$arm == "Intervention"), sum(cl_tab$arm == "Control")))

## ============================================================
## QEPF estimator and test statistic
## ============================================================
qepf_hat_vec_fast <- function(x, u_grid) {
  x <- sort(x); n <- length(x)
  k_idx <- pmin(pmax(ceiling(n * u_grid), 1L), n - 1L)
  cs <- c(0, cumsum(x))
  ((cs[n + 1L] - cs[k_idx + 1L]) / (n - k_idx)) / x[k_idx]
}

build_u_grid <- function(U, n_ref, n_bio, L = 201, eps = 1e-8) {
  u_max_allowed <- min((n_ref - 1) / n_ref, (n_bio - 1) / n_bio) - eps
  u_min <- max(U[1], eps); u_max <- min(U[2], u_max_allowed)
  if (u_max <= u_min) stop("U too close to 1 for given n.")
  seq(u_min, u_max, length.out = L)
}

qepf_diff_and_T <- function(x_ref, x_bio, U = c(0.60, 0.90), L = 201) {
  n_r <- length(x_ref); n_b <- length(x_bio)
  u_grid <- build_u_grid(U, n_r, n_b, L)
  p_ref <- qepf_hat_vec_fast(x_ref, u_grid)
  p_bio <- qepf_hat_vec_fast(x_bio, u_grid)
  ok <- is.finite(p_ref) & is.finite(p_bio)
  diff <- p_ref[ok] - p_bio[ok]
  ne <- (n_r * n_b) / (n_r + n_b)          # scaling factor of the statistic
  list(u = u_grid[ok], diff = diff, T_EQ = sqrt(ne) * max(abs(diff)), ne = ne,
       p_ref = p_ref[ok], p_bio = p_bio[ok])
}

## ============================================================
## Village-level null distributions
## ============================================================
## 1) Permutation: shuffle village arm labels, preserving the number of villages per arm
permute_T_EQ_cluster <- function(df, benefit_var, U = c(0.60, 0.90), L = 201, B = 5000, seed = 2025) {
  set.seed(seed)
  cl_tab <- df %>% distinct(village_id, arm)
  K <- nrow(cl_tab)
  n_int_cl <- sum(cl_tab$arm == "Intervention"); n_ctl_cl <- sum(cl_tab$arm == "Control")
  T_perm <- numeric(B)
  for (b in seq_len(B)) {
    perm_arms <- c(rep("Intervention", n_int_cl), rep("Control", n_ctl_cl))[sample(K)]
    cl_perm <- cl_tab %>% mutate(arm_perm = perm_arms)
    df_perm <- df %>% left_join(cl_perm %>% select(village_id, arm_perm), by = "village_id")
    x_ref_b <- df_perm %>% filter(arm_perm == "Intervention") %>% pull(benefit_var)
    x_bio_b <- df_perm %>% filter(arm_perm == "Control")      %>% pull(benefit_var)
    T_perm[b] <- qepf_diff_and_T(x_ref_b, x_bio_b, U = U, L = L)$T_EQ
  }
  T_perm[is.finite(T_perm)]
}

## 2) Pooled village bootstrap: resample villages from the pooled set and reassign arms
bootstrap_null_T_cluster_pooled <- function(df, benefit_var, U = c(0.60, 0.90), L = 201, B = 5000, seed = 2025) {
  set.seed(seed)
  cl_tab <- df %>% distinct(village_id, arm)
  K <- nrow(cl_tab)
  n_int_cl <- sum(cl_tab$arm == "Intervention")
  T_boot <- numeric(B)
  for (b in seq_len(B)) {
    sampled <- sample(cl_tab$village_id, size = K, replace = TRUE)
    ref_cl <- sampled[seq_len(n_int_cl)]; bio_cl <- sampled[(n_int_cl + 1):K]
    df_b <- df %>%
      filter(village_id %in% sampled) %>%
      mutate(arm_boot = case_when(village_id %in% ref_cl ~ "Intervention",
                                  village_id %in% bio_cl ~ "Control",
                                  TRUE ~ NA_character_)) %>%
      filter(!is.na(arm_boot))
    x_ref_b <- df_b %>% filter(arm_boot == "Intervention") %>% pull(benefit_var)
    x_bio_b <- df_b %>% filter(arm_boot == "Control")      %>% pull(benefit_var)
    if (length(x_ref_b) < 2 || length(x_bio_b) < 2) { T_boot[b] <- NA_real_; next }
    T_boot[b] <- qepf_diff_and_T(x_ref_b, x_bio_b, U = U, L = L)$T_EQ
  }
  T_boot[is.finite(T_boot)]
}

## ============================================================
## Baseline analysis on U = [0.60, 0.90] and Figure 4
## ============================================================
x_int <- df %>% filter(arm == "Intervention") %>% pull(benefit_var)
x_ctl <- df %>% filter(arm == "Control")      %>% pull(benefit_var)
base  <- qepf_diff_and_T(x_int, x_ctl, U = c(0.60, 0.90), L = L_grid)
T_obs <- base$T_EQ

df_plot <- data.frame(u = base$u, Intervention = base$p_ref, Control = base$p_bio, diff = base$diff)
p_curves <- ggplot(df_plot, aes(x = u)) +
  geom_point(aes(y = Intervention, color = "Intervention"), size = 2) +
  geom_point(aes(y = Control,      color = "Control"),      size = 2) +
  geom_line(aes(y = Intervention, color = "Intervention"), linewidth = 0.6, linetype = "dashed", alpha = 0.6) +
  geom_line(aes(y = Control,      color = "Control"),      linewidth = 0.6, linetype = "dashed", alpha = 0.6) +
  scale_color_manual(values = c("Intervention" = "#0B5F88", "Control" = "#C0392B"), name = NULL) +
  labs(x = latex2exp::TeX("$u$"), y = latex2exp::TeX("$\\hat{P}(u)$")) +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
print(p_curves)
ggsave(file.path(out_dir_fig, "fig4_smarter_curves.png"), p_curves, width = 7.5, height = 4.5, dpi = 300)
ggsave(file.path(out_dir_fig, "fig4_smarter_curves.pdf"), p_curves, width = 7.5, height = 4.5)

p_diff <- ggplot(df_plot, aes(x = u, y = diff)) +
  geom_line(color = "#0B5F88", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(title = "D(u) = Intervention - Control", x = "u", y = "D(u)") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave(file.path(out_dir_fig, "smarter_qepf_diff.png"), p_diff, width = 7.5, height = 4.5, dpi = 300)

T_perm <- permute_T_EQ_cluster(df, benefit_var, U = c(0.60, 0.90), L = L_grid, B = B_perm, seed = seed)
T_boot <- bootstrap_null_T_cluster_pooled(df, benefit_var, U = c(0.60, 0.90), L = L_grid, B = B_boot, seed = seed)
cat(sprintf("\nT = %.3f (scaling factor %.1f)\n", T_obs, base$ne))
cat(sprintf("Village permutation: p = %.5f, crit95 = %.3f (B = %d)\n",
            mean(T_perm >= T_obs), quantile(T_perm, 0.95), length(T_perm)))
cat(sprintf("Village pooled bootstrap: p = %.5f, crit95 = %.3f (B = %d)\n",
            mean(T_boot >= T_obs), quantile(T_boot, 0.95), length(T_boot)))
readr::write_csv(data.frame(T_perm = T_perm), file.path(out_dir_tab, "smarter_T_perm.csv"))
readr::write_csv(data.frame(T_boot = T_boot), file.path(out_dir_tab, "smarter_T_boot.csv"))

## ============================================================
## Table 6: sensitivity to the tail interval
## ============================================================
if (RUN_SENSITIVITY) {
  intervals <- list(c(0.10,0.90), c(0.15,0.90), c(0.20,0.90), c(0.25,0.90),
                    c(0.30,0.90), c(0.35,0.90), c(0.40,0.90), c(0.45,0.90),
                    c(0.50,0.90), c(0.55,0.90), c(0.60,0.90), c(0.65,0.90),
                    c(0.70,0.90), c(0.75,0.90), c(0.80,0.90), c(0.85,0.975))
  res_sens <- purrr::map_dfr(intervals, function(U) {
    b  <- qepf_diff_and_T(x_int, x_ctl, U = U, L = L_grid)
    Tp <- permute_T_EQ_cluster(df, benefit_var, U = U, L = L_grid, B = B_sens, seed = seed)
    Tb <- bootstrap_null_T_cluster_pooled(df, benefit_var, U = U, L = L_grid, B = B_sens, seed = seed)
    data.frame(U_lower = U[1], U_upper = U[2], L = L_grid, T_EQ = b$T_EQ, n_eff = b$ne,
               p_perm = mean(Tp >= b$T_EQ), p_boot = mean(Tb >= b$T_EQ),
               crit95_perm = unname(quantile(Tp, 0.95)), crit95_boot = unname(quantile(Tb, 0.95)),
               width = U[2] - U[1])
  })
  print(res_sens)
  readr::write_csv(res_sens, file.path(out_dir_tab, "smarter_sensitivity.csv"))
}
