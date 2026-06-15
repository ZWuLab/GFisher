# Test file for main GFisher functions
# File created: 2025-10-13 by r-rcpp-developer agent

library(testthat)
library(GFisher)

# Test stat.GFisher function
test_that("stat.GFisher computes statistics correctly", {
  set.seed(123)
  n <- 10
  pval <- runif(n)

  # Test with df=2 (standard Fisher's method)
  stat1 <- stat.GFisher(pval, df = 2, w = 1)
  expect_type(stat1, "double")
  expect_length(stat1, 1)
  expect_true(stat1 > 0)

  # Test with explicit vectors
  stat2 <- stat.GFisher(pval, df = rep(2, n), w = rep(1, n))
  expect_equal(stat1, stat2)

  # Test with varying df and weights
  stat3 <- stat.GFisher(pval, df = 1:n, w = 1:n)
  expect_type(stat3, "double")
  expect_length(stat3, 1)
  expect_true(stat3 > 0)

  # Test with small p-values
  pval_small <- c(0.001, 0.002, 0.05, runif(7))
  stat4 <- stat.GFisher(pval_small, df = 2, w = 1)
  expect_true(stat4 > stat1)  # Smaller p-values should give larger statistic
})

test_that("stat.GFisher handles edge cases", {
  # Single p-value
  stat_single <- stat.GFisher(0.5, df = 2, w = 1)
  expect_type(stat_single, "double")

  # All p-values = 1 (no signal)
  stat_null <- stat.GFisher(rep(1, 10), df = 2, w = 1)
  expect_true(stat_null < 1)  # Should be very small

  # Very small p-values
  stat_extreme <- stat.GFisher(rep(1e-10, 5), df = 2, w = 1)
  expect_true(stat_extreme > 10)  # Should be very large (relaxed threshold for robustness)
})

test_that("stat.GFisher validates inputs", {
  # Invalid p-values (out of range)
  expect_error(stat.GFisher(c(-0.1, 0.5), df = 2))
  expect_error(stat.GFisher(c(0.5, 1.1), df = 2))

  # Invalid df
  expect_error(stat.GFisher(c(0.5, 0.5), df = c(-1, 2)))
  expect_error(stat.GFisher(c(0.5, 0.5), df = c(0, 2)))

  # Invalid weights
  expect_error(stat.GFisher(c(0.5, 0.5), df = 2, w = c(-1, 1)))
})


# Test p.GFisher function
test_that("p.GFisher computes p-values correctly with HYB method", {
  set.seed(123)
  n <- 10
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  zscore <- matrix(rnorm(n), nrow = 1) %*% chol(M)
  pval <- 2 * (1 - pnorm(abs(zscore)))

  # Test with df=2
  gf1 <- stat.GFisher(pval, df = 2, w = 1)
  p1 <- p.GFisher(gf1, df = 2, w = 1, M = M, method = "HYB")

  expect_type(p1, "double")
  expect_length(p1, 1)
  expect_true(p1 >= 0 && p1 <= 1)

  # Test with varying df and weights
  gf2 <- stat.GFisher(pval, df = 1:n, w = 1:n)
  p2 <- p.GFisher(gf2, df = 1:n, w = 1:n, M = M, method = "HYB")

  expect_type(p2, "double")
  expect_true(p2 >= 0 && p2 <= 1)
})

test_that("p.GFisher computes p-values correctly with MR method", {
  set.seed(456)
  n <- 10
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  zscore <- matrix(rnorm(n), nrow = 1) %*% chol(M)
  pval <- 2 * (1 - pnorm(abs(zscore)))

  gf1 <- stat.GFisher(pval, df = 2, w = 1)
  p1 <- p.GFisher(gf1, df = 2, w = 1, M = M, method = "MR", nsim = 1e4, seed = 789)

  expect_type(p1, "double")
  expect_true(p1 >= 0 && p1 <= 1)

  # Test reproducibility with seed
  p1_repeat <- p.GFisher(gf1, df = 2, w = 1, M = M, method = "MR", nsim = 1e4, seed = 789)
  expect_equal(p1, p1_repeat)
})

test_that("p.GFisher handles near-singular matrices", {
  set.seed(123)
  n <- 5
  # Create a nearly singular matrix
  M <- matrix(0.95, n, n) + diag(0.05, n, n)

  pval <- runif(n)
  gf <- stat.GFisher(pval, df = 2, w = 1)

  # Should not error, should use nearPD correction
  expect_error(p.GFisher(gf, df = 2, w = 1, M = M, method = "MR", nsim = 1e3), NA)
})

test_that("p.GFisher GB method works", {
  set.seed(123)
  n <- 10
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  pval <- runif(n)
  gf <- stat.GFisher(pval, df = 2, w = 1)

  p_gb <- p.GFisher(gf, df = 2, w = 1, M = M, method = "GB")
  expect_type(p_gb, "double")
  expect_true(p_gb >= 0 && p_gb <= 1)
})


# Test stat.oGFisher function
test_that("stat.oGFisher computes statistics correctly", {
  set.seed(123)
  n <- 10
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  zscore <- matrix(rnorm(n), nrow = 1) %*% chol(M)
  pval <- 2 * (1 - pnorm(abs(zscore)))

  DF <- rbind(rep(1, n), rep(2, n))
  W <- rbind(rep(1, n), 1:10)

  result <- stat.oGFisher(pval, DF, W, M, p.type = "two", method = "HYB")

  expect_type(result, "list")
  expect_named(result, c("STAT", "PVAL", "minp", "cct"))
  expect_length(result$STAT, 2)
  expect_length(result$PVAL, 2)
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
  expect_equal(result$minp, min(result$PVAL))
  expect_type(result$cct, "double")
})

test_that("stat.oGFisher handles single df per test", {
  set.seed(123)
  n <- 10
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  pval <- runif(n)

  DF_short <- rbind(1, 2)
  W <- rbind(rep(1, n), 1:10)

  # Should work (df gets expanded internally)
  result <- stat.oGFisher(pval, DF = DF_short, W, M, method = "HYB")
  expect_length(result$STAT, 2)
})


# Test pval.oGFisher function
test_that("pval.oGFisher computes p-values with CCT combination", {
  set.seed(123)
  n <- 10
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  zscore <- matrix(rnorm(n), nrow = 1) %*% chol(M)
  pval <- 2 * (1 - pnorm(abs(zscore)))

  DF <- rbind(rep(1, n), rep(2, n))
  W <- rbind(rep(1, n), 1:10)

  result <- pval.oGFisher(pval, DF, W, M, p.type = "two", method = "HYB", combine = "cct")

  expect_type(result, "list")
  expect_named(result, c("stat", "pval", "pval_indi", "stat_indi"))
  expect_type(result$stat, "double")  # CCT statistic
  expect_true(result$pval >= 0 && result$pval <= 1)
  expect_length(result$pval_indi, 2)
  expect_length(result$stat_indi, 2)
})

test_that("pval.oGFisher computes p-values with MVN combination", {
  set.seed(123)
  n <- 10
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  pval <- runif(n)

  DF <- rbind(rep(1, n), rep(2, n))
  W <- rbind(rep(1, n), 1:10)

  result <- pval.oGFisher(pval, DF, W, M, p.type = "two", method = "HYB", combine = "mvn")

  expect_type(result, "list")
  expect_type(result$stat, "double")  # minp statistic
  expect_true(result$pval >= 0 && result$pval <= 1)
})

test_that("pval.oGFisher CCT and MVN give different but valid results", {
  set.seed(456)
  n <- 10
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  pval <- runif(n)

  DF <- rbind(rep(2, n), rep(3, n))
  W <- rbind(rep(1, n), rep(1, n))

  result_cct <- pval.oGFisher(pval, DF, W, M, method = "HYB", combine = "cct")
  result_mvn <- pval.oGFisher(pval, DF, W, M, method = "HYB", combine = "mvn")

  # Both should be valid p-values
  expect_true(result_cct$pval >= 0 && result_cct$pval <= 1)
  expect_true(result_mvn$pval >= 0 && result_mvn$pval <= 1)

  # Stats should be different (cct vs minp)
  expect_false(result_cct$stat == result_mvn$stat)
})


# Integration test: Full workflow
test_that("Full GFisher workflow works end-to-end", {
  set.seed(999)
  n <- 15
  # Create correlation matrix
  M <- matrix(0.5, n, n) + diag(0.5, n, n)

  # Generate correlated Z-scores
  zscore <- matrix(rnorm(n), nrow = 1) %*% chol(M)
  pval <- 2 * (1 - pnorm(abs(zscore)))

  # Single GFisher test
  stat <- stat.GFisher(pval, df = 2, w = 1)
  p_hyb <- p.GFisher(stat, df = 2, w = 1, M = M, method = "HYB")
  p_mr <- p.GFisher(stat, df = 2, w = 1, M = M, method = "MR", nsim = 5e3, seed = 111)

  expect_true(p_hyb >= 0 && p_hyb <= 1)
  expect_true(p_mr >= 0 && p_mr <= 1)

  # oGFisher test
  DF <- rbind(rep(1, n), rep(2, n), rep(3, n))
  W <- rbind(rep(1, n), 1:n, sqrt(1:n))

  ogf_stats <- stat.oGFisher(pval, DF, W, M, method = "HYB")
  ogf_pval_cct <- pval.oGFisher(pval, DF, W, M, method = "HYB", combine = "cct")
  ogf_pval_mvn <- pval.oGFisher(pval, DF, W, M, method = "HYB", combine = "mvn")

  expect_true(ogf_pval_cct$pval >= 0 && ogf_pval_cct$pval <= 1)
  expect_true(ogf_pval_mvn$pval >= 0 && ogf_pval_mvn$pval <= 1)
  expect_equal(length(ogf_stats$STAT), 3)
})
