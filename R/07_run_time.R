## =====================================================================
##  07_run_time.R  --  Section 6.5: time one permutation test of equal
##  persistence (n = 500 per arm, B = 2000, 101 grid points on U).
##  Reported in the paper: under one second (elapsed) on a standard laptop
##  (0.49-0.70 s in our runs).
##  Run from anywhere:  source("R/07_run_time.R")
## =====================================================================
set.seed(1)
x  <- rgamma(500, shape = 4)                 # arm 1
y  <- rgamma(500, shape = 4)                 # arm 2
ug <- seq(0.60, 0.90, length.out = 101)      # grid on U = [0.60, 0.90]

kidx  <- function(n, u) pmin(pmax(as.integer(ceiling(n * u - 1e-9)), 1L), n - 1L)
Pgrid <- function(xs, k) { n <- length(xs); cs <- cumsum(xs); (cs[n] - cs[k]) / ((n - k) * xs[k]) }

perm_test <- function(x, y, ug, B = 2000) {
  n1 <- length(x); n2 <- length(y); N <- n1 + n2
  k1 <- kidx(n1, ug); k2 <- kidx(n2, ug)
  T0 <- max(abs(Pgrid(sort(x), k1) - Pgrid(sort(y), k2)))
  z  <- sort(c(x, y)); Tb <- numeric(B)
  for (b in seq_len(B)) {
    p <- sample.int(N)
    Tb[b] <- max(abs(Pgrid(z[sort.int(p[1:n1], method = "radix")], k1) -
                       Pgrid(z[sort.int(p[(n1 + 1):N], method = "radix")], k2)))
  }
  (1 + sum(Tb >= T0)) / (B + 1)              # permutation p-value
}

system.time(perm_test(x, y, ug, B = 2000))