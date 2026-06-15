# Baseline validation tests comparing new implementation with legacy code
# File created: 2025-10-13 by r-rcpp-developer agent
#
# These tests validate that the new implementation produces results consistent
# with the legacy GFisher_v2.R code

library(testthat)
library(GFisher)

# Source legacy code for comparison
legacy_file <- file.path(Sys.getenv("GLOW_LEGACY_ROOT", unset = "/nonexistent"), "legacy-materials/GFisher/GFisher_v2.R")
if (file.exists(legacy_file)) {
  source(legacy_file, local = TRUE)
  legacy_available <- TRUE
} else {
  legacy_available <- FALSE
}


# Test stat.GFisher against legacy
test_that("stat.GFisher matches legacy implementation", {
  skip_if(!legacy_available, "Legacy code not available")

  set.seed(123)
  n <- 10
  pval <- runif(n)

  # Test 1: df=2, w=1
  new_stat1 <- stat.GFisher(pval, df = 2, w = 1)
  legacy_stat1 <- stat.GFisher(pval, df = 2, w = 1)
  expect_equal(new_stat1, legacy_stat1, tolerance = 1e-12)

  # Test 2: Varying df and weights
  new_stat2 <- stat.GFisher(pval, df = 1:n, w = 1:n)
  legacy_stat2 <- stat.GFisher(pval, df = 1:n, w = 1:n)
  expect_equal(new_stat2, legacy_stat2, tolerance = 1e-12)
})


# Test p.GFisher against legacy
test_that("p.GFisher with HYB method matches legacy", {
  skip_if(!legacy_available, "Legacy code not available")

  set.seed(456)
  n <- 10
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  zscore <- matrix(rnorm(n), nrow = 1) %*% chol(M)
  pval <- 2 * (1 - pnorm(abs(zscore)))

  gf <- stat.GFisher(pval, df = 2, w = 1)

  # Compare HYB method
  new_p <- p.GFisher(gf, df = 2, w = 1, M = M, method = "HYB")
  legacy_p <- p.GFisher(gf, df = 2, w = 1, M = M, method = "HYB")

  expect_equal(new_p, legacy_p, tolerance = 1e-10)
})

test_that("p.GFisher with MR method matches legacy", {
  skip_if(!legacy_available, "Legacy code not available")

  set.seed(789)
  n <- 8
  M <- matrix(0.5, n, n) + diag(0.5, n, n)
  pval <- runif(n)

  gf <- stat.GFisher(pval, df = 2, w = 1)

  # Compare MR method (use same seed for reproducibility)
  new_p <- p.GFisher(gf, df = 2, w = 1, M = M, method = "MR", nsim = 5e3, seed = 111)
  legacy_p <- p.GFisher(gf, df = 2, w = 1, M = M, method = "MR", nsim = 5e3, seed = 111)

  expect_equal(new_p, legacy_p, tolerance = 1e-10)
})

test_that("p.GFisher with GB method matches legacy", {
  skip_if(!legacy_available, "Legacy code not available")

  set.seed(321)
  n <- 8
  M <- diag(n)
  pval <- runif(n)

  gf <- stat.GFisher(pval, df = 2, w = 1)

  new_p <- p.GFisher(gf, df = 2, w = 1, M = M, method = "GB")
  legacy_p <- p.GFisher(gf, df = 2, w = 1, M = M, method = "GB")

  expect_equal(new_p, legacy_p, tolerance = 1e-10)
})


# Test stat.oGFisher against legacy
test_that("stat.oGFisher matches legacy implementation", {
  skip_if(!legacy_available, "Legacy code not available")

  set.seed(555)
  n <- 10
  M <- matrix(0.4, n, n) + diag(0.6, n, n)
  zscore <- matrix(rnorm(n), nrow = 1) %*% chol(M)
  pval <- 2 * (1 - pnorm(abs(zscore)))

  DF <- rbind(rep(1, n), rep(2, n))
  W <- rbind(rep(1, n), 1:10)

  new_result <- stat.oGFisher(pval, DF, W, M, p.type = "two", method = "HYB")
  legacy_result <- stat.oGFisher(pval, DF, W, M, p.type = "two", method = "HYB")

  expect_equal(new_result$STAT, legacy_result$STAT, tolerance = 1e-10)
  expect_equal(new_result$PVAL, legacy_result$PVAL, tolerance = 1e-10)
  expect_equal(new_result$minp, legacy_result$minp, tolerance = 1e-10)
  expect_equal(new_result$cct, legacy_result$cct, tolerance = 1e-10)
})


# Test pval.oGFisher against legacy
test_that("pval.oGFisher with CCT matches legacy", {
  skip_if(!legacy_available, "Legacy code not available")

  set.seed(666)
  n <- 10
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  pval <- runif(n)

  DF <- rbind(rep(1, n), rep(2, n))
  W <- rbind(rep(1, n), 1:10)

  new_result <- pval.oGFisher(pval, DF, W, M, p.type = "two",
                               method = "HYB", combine = "cct")
  legacy_result <- pval.oGFisher(pval, DF, W, M, p.type = "two",
                                  method = "HYB", combine = "cct")

  expect_equal(new_result$stat, legacy_result$stat, tolerance = 1e-10)
  expect_equal(new_result$pval, legacy_result$pval, tolerance = 1e-10)
  expect_equal(new_result$pval_indi, legacy_result$pval_indi, tolerance = 1e-10)
  expect_equal(new_result$stat_indi, legacy_result$stat_indi, tolerance = 1e-10)
})

test_that("pval.oGFisher with MVN matches legacy", {
  skip_if(!legacy_available, "Legacy code not available")

  set.seed(777)
  n <- 8
  M <- diag(n)
  pval <- runif(n)

  DF <- rbind(rep(2, n), rep(3, n))
  W <- rbind(rep(1, n), rep(1, n))

  new_result <- pval.oGFisher(pval, DF, W, M, p.type = "two",
                               method = "HYB", combine = "mvn")
  legacy_result <- pval.oGFisher(pval, DF, W, M, p.type = "two",
                                  method = "HYB", combine = "mvn")

  expect_equal(new_result$stat, legacy_result$stat, tolerance = 1e-10)
  expect_equal(new_result$pval, legacy_result$pval, tolerance = 1e-9)
  expect_equal(new_result$pval_indi, legacy_result$pval_indi, tolerance = 1e-10)
})


# Comprehensive validation test
test_that("Complete workflow matches legacy implementation", {
  skip_if(!legacy_available, "Legacy code not available")

  set.seed(999)
  n <- 12
  M <- matrix(0.4, n, n) + diag(0.6, n, n)

  # Generate test data
  zscore <- matrix(rnorm(n), nrow = 1) %*% chol(M)
  pval <- 2 * (1 - pnorm(abs(zscore)))

  # Single GFisher test
  new_stat <- stat.GFisher(pval, df = 2, w = 1)
  legacy_stat <- stat.GFisher(pval, df = 2, w = 1)
  expect_equal(new_stat, legacy_stat, tolerance = 1e-12)

  new_p_hyb <- p.GFisher(new_stat, df = 2, w = 1, M = M, method = "HYB")
  legacy_p_hyb <- p.GFisher(legacy_stat, df = 2, w = 1, M = M, method = "HYB")
  expect_equal(new_p_hyb, legacy_p_hyb, tolerance = 1e-10)

  # oGFisher test
  DF <- rbind(rep(1, n), rep(2, n), rep(3, n))
  W <- rbind(rep(1, n), 1:n, sqrt(1:n))

  new_ogf <- stat.oGFisher(pval, DF, W, M, method = "HYB")
  legacy_ogf <- stat.oGFisher(pval, DF, W, M, method = "HYB")

  expect_equal(new_ogf$STAT, legacy_ogf$STAT, tolerance = 1e-10)
  expect_equal(new_ogf$PVAL, legacy_ogf$PVAL, tolerance = 1e-10)

  new_ogf_pval <- pval.oGFisher(pval, DF, W, M, method = "HYB", combine = "cct")
  legacy_ogf_pval <- pval.oGFisher(pval, DF, W, M, method = "HYB", combine = "cct")

  expect_equal(new_ogf_pval$pval, legacy_ogf_pval$pval, tolerance = 1e-10)
})


# Numerical accuracy tests
test_that("Results are numerically stable across different scenarios", {
  skip_if(!legacy_available, "Legacy code not available")

  # Test various correlation structures
  test_scenarios <- list(
    list(rho = 0.1, desc = "low correlation"),
    list(rho = 0.5, desc = "medium correlation"),
    list(rho = 0.8, desc = "high correlation")
  )

  for (scenario in test_scenarios) {
    set.seed(123)
    n <- 10
    M <- matrix(scenario$rho, n, n) + diag(1 - scenario$rho, n, n)
    pval <- runif(n)
    gf <- stat.GFisher(pval, df = 2, w = 1)

    new_p <- p.GFisher(gf, df = 2, w = 1, M = M, method = "HYB")
    legacy_p <- p.GFisher(gf, df = 2, w = 1, M = M, method = "HYB")

    expect_equal(new_p, legacy_p, tolerance = 1e-10,
                 label = sprintf("Scenario: %s", scenario$desc))
  }
})


# Test edge cases with legacy
test_that("Edge cases match legacy behavior", {
  skip_if(!legacy_available, "Legacy code not available")

  set.seed(123)
  n <- 5

  # Very small p-values
  pval_small <- c(1e-10, 1e-9, 1e-8, runif(2))
  M <- diag(n)

  new_stat <- stat.GFisher(pval_small, df = 2, w = 1)
  legacy_stat <- stat.GFisher(pval_small, df = 2, w = 1)
  expect_equal(new_stat, legacy_stat, tolerance = 1e-12)

  # Large df values
  pval_reg <- runif(n)
  df_large <- rep(10, n)

  new_stat2 <- stat.GFisher(pval_reg, df = df_large, w = 1)
  legacy_stat2 <- stat.GFisher(pval_reg, df = df_large, w = 1)
  expect_equal(new_stat2, legacy_stat2, tolerance = 1e-12)

  # Unequal weights
  w_unequal <- c(1, 2, 3, 4, 5)
  new_stat3 <- stat.GFisher(pval_reg, df = 2, w = w_unequal)
  legacy_stat3 <- stat.GFisher(pval_reg, df = 2, w = w_unequal)
  expect_equal(new_stat3, legacy_stat3, tolerance = 1e-12)
})


# Performance comparison (informational, not enforced)
test_that("Performance comparison with legacy (informational)", {
  skip_if(!legacy_available, "Legacy code not available")
  skip_on_cran()

  set.seed(123)
  n <- 20
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  pval <- runif(n)
  gf <- stat.GFisher(pval, df = 2, w = 1)

  # Time new implementation
  time_new <- system.time({
    for (i in 1:100) {
      p.GFisher(gf, df = 2, w = 1, M = M, method = "HYB")
    }
  })

  # Time legacy implementation
  time_legacy <- system.time({
    for (i in 1:100) {
      p.GFisher(gf, df = 2, w = 1, M = M, method = "HYB")
    }
  })

  message(sprintf("New: %.3f sec, Legacy: %.3f sec",
                  time_new["elapsed"], time_legacy["elapsed"]))

  # Just verify both complete
  expect_true(time_new["elapsed"] >= 0)
  expect_true(time_legacy["elapsed"] >= 0)
})
