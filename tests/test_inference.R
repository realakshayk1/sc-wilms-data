library(testthat)

# Load the patient-level permutation machinery (09) and composition helpers (12).
# Source from inside the script dir so their internal relative source() calls
# resolve; main() is not run (guarded by sys.nframe()). Functions land in globalenv.
local({
  phase1 <- if (dir.exists("phase1_mechanotypes")) "phase1_mechanotypes" else file.path("..", "phase1_mechanotypes")
  old <- setwd(phase1); on.exit(setwd(old))
  suppressWarnings(suppressMessages({
    source("09_distributional_validation.R")
    source("12_composition_analysis.R")
  }))
})

skip_if_no_transport <- function() {
  if (!requireNamespace("transport", quietly = TRUE)) skip("transport not installed")
}

# Build synthetic (cell x sample) data: n_samp samples per group, n_cell cells each,
# with a controllable true between-group shift at the SAMPLE level.
make_nested <- function(n_samp = 8, n_cell = 200, shift = 0, seed = 1) {
  set.seed(seed)
  rows <- list()
  for (g in c("pos", "neg")) {
    for (i in seq_len(n_samp)) {
      sid <- paste(g, i, sep = "_")
      mu <- if (g == "pos") shift else 0
      # sample-level random intercept so cells within a sample are correlated
      mu <- mu + rnorm(1, 0, 0.2)
      rows[[length(rows) + 1]] <- data.frame(
        score = rnorm(n_cell, mu, 1), sample_id = sid, grp = g,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

test_that("patient-level permutation is NOT anti-conservative under the null", {
  skip_if_no_transport()
  # No true effect -> p should be roughly uniform, definitely not pinned to the floor.
  ps <- vapply(1:20, function(s) {
    d <- make_nested(shift = 0, seed = s)
    sample_perm_p_w1(d$score, d$sample_id, d$grp, n_perm = 199, seed = s)$p_perm
  }, numeric(1))
  # The old (buggy) version returned ~0.001 every time. Guard against that.
  expect_gt(mean(ps), 0.2)
  expect_gt(max(ps), 0.5)
})

test_that("patient-level permutation detects a strong true sample-level effect", {
  skip_if_no_transport()
  d <- make_nested(shift = 2.0, seed = 42)
  res <- sample_perm_p_w1(d$score, d$sample_id, d$grp, n_perm = 199, seed = 42)
  expect_lt(res$p_perm, 0.05)
  expect_true(is.finite(res$w1) && res$w1 > 0)
})

test_that("permutation requires >= 2 samples per group", {
  skip_if_no_transport()
  d <- make_nested(n_samp = 1, shift = 1, seed = 3)
  res <- sample_perm_p_w1(d$score, d$sample_id, d$grp, n_perm = 99, seed = 3)
  expect_true(is.na(res$p_perm))
})

test_that("CLR transform rows are centred (sum to ~0)", {
  m <- matrix(c(0.6, 0.3, 0.1, 0.2, 0.2, 0.6), nrow = 2, byrow = TRUE)
  z <- clr(m)
  expect_true(all(abs(rowSums(z)) < 1e-8))
})
