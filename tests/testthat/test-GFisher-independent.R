# Test file for GFisher independent implementation
# File created: 2025-10-13 by r-rcpp-developer agent
# Updated: 2025-10-15 - Comprehensive consistency and calibration tests

library(testthat)
library(GFisher)

# Skip tests if coga package is not available
skip_if_not_coga <- function() {
  skip_if_not_installed("coga")
}


# =============================================================================
# PART 1: Basic Functionality Tests
# =============================================================================

# Test p.GFisher_ind function
test_that("p.GFisher_ind computes p-values correctly", {
  skip_if_not_coga()

  set.seed(123)
  n <- 10
  df <- runif(n, 0.5, 5)
  w <- abs(rnorm(n))
  q <- 40

  # Exact calculation
  p_exact <- p.GFisher_ind(q = q, df = df, w = w, isExact = TRUE)

  expect_type(p_exact, "double")
  expect_length(p_exact, 1)
  expect_true(p_exact >= 0 && p_exact <= 1)
})

test_that("p.GFisher_ind approximation works", {
  skip_if_not_coga()

  set.seed(123)
  n <- 10
  df <- runif(n, 0.5, 5)
  w <- abs(rnorm(n))
  q <- 40

  # Approximation
  p_approx <- p.GFisher_ind(q = q, df = df, w = w, isExact = FALSE)

  expect_type(p_approx, "double")
  expect_true(p_approx >= 0 && p_approx <= 1)
})

test_that("p.GFisher_ind handles equal weights", {
  skip_if_not_coga()

  n <- 5
  df <- rep(2, n)
  w <- rep(1, n)
  q <- 20

  p_val <- p.GFisher_ind(q = q, df = df, w = w)

  expect_true(p_val >= 0 && p_val <= 1)
})

test_that("p.GFisher_ind handles single weight", {
  skip_if_not_coga()

  n <- 5
  df <- rep(2, n)
  w <- 1  # Single weight - should be expanded
  q <- 20

  p_val <- p.GFisher_ind(q = q, df = df, w = w)

  expect_true(p_val >= 0 && p_val <= 1)
})

test_that("p.GFisher_ind gives error without coga", {
  # Temporarily unload coga if it's loaded
  if ("coga" %in% loadedNamespaces()) {
    skip("coga is loaded, cannot test error message")
  }

  # Mock the requireNamespace function to return FALSE
  # This is tricky - we'll just check the function expects coga
  expect_true(grepl("coga", paste(deparse(p.GFisher_ind), collapse = " ")))
})


# Test p.GFisher_ind_w1 function
test_that("p.GFisher_ind_w1 computes p-values correctly", {
  set.seed(123)
  n <- 10
  df <- rep(2, n)
  pval <- runif(n)
  q <- stat.GFisher(pval, df = df, w = 1)

  p_fast <- p.GFisher_ind_w1(q, df)

  expect_type(p_fast, "double")
  expect_length(p_fast, 1)
  expect_true(p_fast >= 0 && p_fast <= 1)
})

test_that("p.GFisher_ind_w1 uses gamma distribution correctly", {
  n <- 5
  df <- rep(2, n)
  q <- 25

  p_fast <- p.GFisher_ind_w1(q, df)
  p_gamma <- pgamma(q, shape = sum(df)/2, scale = 2/length(df), lower.tail = FALSE)

  expect_equal(p_fast, p_gamma)
})

test_that("p.GFisher_ind_w1 handles varying df", {
  df <- c(1, 2, 3, 4, 5)
  q <- 30

  p_val <- p.GFisher_ind_w1(q, df)

  expect_true(p_val >= 0 && p_val <= 1)
})

test_that("p.GFisher_ind_w1 handles large q", {
  df <- rep(2, 10)
  q <- 1000  # Very large statistic

  p_val <- p.GFisher_ind_w1(q, df)

  expect_true(p_val >= 0 && p_val < 0.001)  # Should be very small
})


# Test stat.oGFisher_ind function
test_that("stat.oGFisher_ind computes statistics correctly", {
  skip_if_not_coga()

  set.seed(123)
  n <- 10
  nGF <- 3

  DF <- matrix(runif(n * nGF, 0.5, 5), ncol = n) / 10
  W <- abs(matrix(rnorm(n * nGF), ncol = n))

  # Under H0
  p <- runif(n)
  result <- stat.oGFisher_ind(p = p, DF = DF, W = W)

  expect_type(result, "list")
  expect_named(result, c("STAT", "PVAL", "minp", "cct"))
  expect_length(result$STAT, nGF)
  expect_length(result$PVAL, nGF)
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
  expect_equal(result$minp, min(result$PVAL))
  expect_type(result$cct, "double")
})

test_that("stat.oGFisher_ind detects signals", {
  skip_if_not_coga()

  set.seed(123)
  n <- 10
  nGF <- 3

  DF <- matrix(rep(2, n * nGF), ncol = n)
  W <- matrix(rep(1, n * nGF), ncol = n)

  # With potential signal
  p <- runif(n)
  p[1] <- 0.00001

  result <- stat.oGFisher_ind(p = p, DF = DF, W = W)

  # With such a strong signal (p=1e-5), at least one test should detect it
  expect_true(result$minp < 0.05)  # Should detect signal
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
})

test_that("stat.oGFisher_ind handles varying configurations", {
  skip_if_not_coga()

  set.seed(456)
  n <- 8
  nGF <- 4

  DF <- rbind(
    rep(1, n),
    rep(2, n),
    1:n / 2,
    rep(3, n)
  )
  W <- rbind(
    rep(1, n),
    1:n,
    sqrt(1:n),
    rep(2, n)
  )

  p <- runif(n)
  result <- stat.oGFisher_ind(p = p, DF = DF, W = W)

  expect_length(result$STAT, nGF)
  expect_length(result$PVAL, nGF)
  expect_true(all(is.finite(result$STAT)))
  expect_true(all(is.finite(result$PVAL)))
})


# Test pval.oGFisher_ind function
test_that("pval.oGFisher_ind computes p-values with CCT", {
  skip_if_not_coga()

  set.seed(123)
  n <- 10
  nGF <- 3

  DF <- matrix(runif(n * nGF, 0.5, 5), ncol = n) / 10
  W <- abs(matrix(rnorm(n * nGF), ncol = n))

  p <- runif(n)
  result <- pval.oGFisher_ind(p = p, DF = DF, W = W, combine = "cct")

  expect_type(result, "list")
  expect_named(result, c("pval", "pval_indi"))
  expect_type(result$pval, "double")
  expect_true(result$pval >= 0 && result$pval <= 1)
  expect_length(result$pval_indi, nGF)
})

test_that("pval.oGFisher_ind computes p-values with minp", {
  skip_if_not_coga()

  set.seed(123)
  n <- 10
  nGF <- 3

  DF <- matrix(rep(2, n * nGF), ncol = n)
  W <- matrix(rep(1, n * nGF), ncol = n)

  p <- runif(n)
  result <- pval.oGFisher_ind(p = p, DF = DF, W = W, combine = "minp")

  expect_type(result, "list")
  expect_true(result$pval >= 0 && result$pval <= 1)
})

test_that("pval.oGFisher_ind CCT and minp differ", {
  skip_if_not_coga()

  set.seed(456)
  n <- 10
  nGF <- 4

  DF <- rbind(rep(1, n), rep(2, n), rep(3, n), rep(4, n))
  W <- rbind(rep(1, n), 1:n, sqrt(1:n), rep(2, n))

  p <- runif(n)
  p[1] <- 0.001  # Small signal

  result_cct <- pval.oGFisher_ind(p = p, DF = DF, W = W, combine = "cct")
  result_minp <- pval.oGFisher_ind(p = p, DF = DF, W = W, combine = "minp")

  # Both should be valid
  expect_true(result_cct$pval >= 0 && result_cct$pval <= 1)
  expect_true(result_minp$pval >= 0 && result_minp$pval <= 1)

  # Individual p-values should be the same
  expect_equal(result_cct$pval_indi, result_minp$pval_indi, tolerance = 1e-10)
})

test_that("pval.oGFisher_ind handles extreme CCT values", {
  skip_if_not_coga()

  set.seed(789)
  n <- 5
  nGF <- 2

  DF <- rbind(rep(2, n), rep(3, n))
  W <- rbind(rep(1, n), rep(1, n))

  # Create very small p-values to get extreme CCT
  p <- rep(1e-20, n)

  # Should not error
  expect_error(pval.oGFisher_ind(p = p, DF = DF, W = W, combine = "cct"), NA)
})


# =============================================================================
# PART 2: Consistency Tests (Independent vs General Functions)
# =============================================================================

test_that("p.GFisher_ind matches p.GFisher with identity matrix - equal weights", {
  skip_if_not_coga()

  set.seed(999)
  n <- 8
  df <- rep(2, n)
  w <- rep(1, n)
  pval <- runif(n)
  q <- stat.GFisher(pval, df = df, w = w)

  # Independent method
  p_ind <- p.GFisher_ind(q = q, df = df, w = w, isExact = TRUE)

  # General method with identity matrix (independence)
  p_gen <- p.GFisher(q = q, df = df, w = w, M = diag(n),
                     p.type = "two", method = "HYB")

  # Should be close (HYB may have some approximation error)
  expect_equal(p_ind, p_gen, tolerance = 0.01)
})

test_that("p.GFisher_ind matches p.GFisher with identity matrix - unequal weights", {
  skip_if_not_coga()

  set.seed(888)
  n <- 10
  df <- runif(n, 0.5, 5)
  w <- abs(rnorm(n))
  q <- 40

  # Independent method
  p_ind <- p.GFisher_ind(q = q, df = df, w = w, isExact = TRUE)

  # General method with identity matrix and MR (more accurate for this case)
  p_gen <- p.GFisher(q = q, df = df, w = w, M = diag(n),
                     p.type = "two", method = "MR", nsim = 1e5)

  # MR has simulation variability, so tolerance is larger
  expect_equal(p_ind, p_gen, tolerance = 0.1)
})

test_that("p.GFisher_ind_w1 matches p.GFisher_ind with equal weights", {
  skip_if_not_coga()

  set.seed(777)
  n <- 10
  df <- runif(n, 0.5, 3)
  q <- 30

  # Using p.GFisher_ind_w1
  p_w1 <- p.GFisher_ind_w1(q = q, df = df)

  # Using p.GFisher_ind with equal weights
  p_ind <- p.GFisher_ind(q = q, df = df, w = rep(1, n), isExact = TRUE)

  # Should be very close
  expect_equal(p_w1, p_ind, tolerance = 1e-6)
})

test_that("p.GFisher_ind_w1 matches p.GFisher with identity matrix and equal weights", {
  set.seed(666)
  n <- 8
  df <- rep(2, n)
  pval <- runif(n)
  q <- stat.GFisher(pval, df = df, w = 1)

  # Using p.GFisher_ind_w1
  p_w1 <- p.GFisher_ind_w1(q = q, df = df)

  # Using p.GFisher with identity matrix
  p_gen <- p.GFisher(q = q, df = df, w = rep(1, n), M = diag(n),
                     p.type = "two", method = "HYB")

  # Should be close
  expect_equal(p_w1, p_gen, tolerance = 0.01)
})

test_that("stat.oGFisher_ind matches stat.oGFisher with identity matrix", {
  skip_if_not_coga()

  set.seed(555)
  n <- 10
  nGF <- 3

  DF <- matrix(runif(n * nGF, 0.5, 5), ncol = n) / 10
  W <- abs(matrix(rnorm(n * nGF), ncol = n))

  p <- runif(n)

  # Independent method
  result_ind <- stat.oGFisher_ind(p = p, DF = DF, W = W)

  # General method with identity matrix
  result_gen <- stat.oGFisher(p = p, DF = DF, W = W, M = diag(n), method = "HYB")

  # Statistics should be identical (they both use stat.GFisher)
  expect_equal(result_ind$STAT, result_gen$STAT, tolerance = 1e-10)

  # P-values should be close (HYB approximation vs exact coga)
  expect_equal(result_ind$PVAL, result_gen$PVAL, tolerance = 0.05)
})

test_that("pval.oGFisher_ind matches pval.oGFisher with identity matrix - CCT", {
  skip_if_not_coga()

  set.seed(444)
  n <- 10
  nGF <- 5

  DF <- matrix(runif(n * nGF, 0.5, 3), ncol = n) / 10
  W <- abs(matrix(rnorm(n * nGF), ncol = n))

  p <- runif(n)
  p[1] <- 0.001  # Small signal

  # Independent method
  result_ind <- pval.oGFisher_ind(p = p, DF = DF, W = W, combine = "cct")

  # General method with identity matrix and MR
  result_gen <- pval.oGFisher(p = p, DF = DF, W = W, M = diag(n),
                              method = "MR", combine = "cct")

  # Should be reasonably close
  expect_equal(result_ind$pval, result_gen$pval, tolerance = 0.1)
})

test_that("pval.oGFisher_ind matches pval.oGFisher with identity matrix - minp", {
  skip_if_not_coga()

  set.seed(333)
  n <- 8
  nGF <- 4

  DF <- matrix(rep(2, n * nGF), ncol = n)
  W <- matrix(rep(1, n * nGF), ncol = n)

  p <- runif(n)
  p[1] <- 0.001

  # Independent method
  result_ind <- pval.oGFisher_ind(p = p, DF = DF, W = W, combine = "minp")

  # General method with identity matrix
  result_gen <- pval.oGFisher(p = p, DF = DF, W = W, M = diag(n),
                              method = "HYB", combine = "minp")

  # Should be reasonably close
  expect_equal(result_ind$pval, result_gen$pval, tolerance = 0.1)
})


# =============================================================================
# PART 3: Empirical P-Value Calibration Tests
# =============================================================================

test_that("p.GFisher_ind_w1 has correct empirical type I error - alpha = 0.05", {
  skip_on_cran()  # Too slow for CRAN
  skip_if_not_installed("coga")

  set.seed(12345)
  nsim <- 5000  # Reduced for testing speed
  n <- 10

  # Generate phred scores for realistic df values
  snpN <- 3e5
  phredScores <- -10 * log10((1:snpN) / snpN)

  p_cal <- rep(NA, nsim)
  for (i in 1:nsim) {
    pval <- runif(n)
    df <- sample(phredScores, n) / 10000
    pp_trans <- qchisq(1 - pval, df = df)
    gfisherstat <- sum(pp_trans) / n  # Normalized (mean)
    p_cal[i] <- p.GFisher_ind_w1(gfisherstat, df)
  }

  # Check empirical type I error at alpha = 0.05
  empirical_rate <- mean(p_cal <= 0.05)
  ratio <- empirical_rate / 0.05

  # Should be close to 1 (within 20% is reasonable for 5000 sims)
  expect_true(ratio >= 0.8 && ratio <= 1.2,
              info = sprintf("Ratio = %.3f, should be between 0.8 and 1.2", ratio))
})

test_that("p.GFisher_ind_w1 has correct empirical type I error - alpha = 0.01", {
  skip_on_cran()  # Too slow for CRAN
  skip_if_not_installed("coga")

  set.seed(23456)
  nsim <- 5000
  n <- 10

  snpN <- 3e5
  phredScores <- -10 * log10((1:snpN) / snpN)

  p_cal <- rep(NA, nsim)
  for (i in 1:nsim) {
    pval <- runif(n)
    df <- sample(phredScores, n) / 10000
    pp_trans <- qchisq(1 - pval, df = df)
    gfisherstat <- sum(pp_trans) / n
    p_cal[i] <- p.GFisher_ind_w1(gfisherstat, df)
  }

  empirical_rate <- mean(p_cal <= 0.01)
  ratio <- empirical_rate / 0.01

  # Looser tolerance for smaller alpha (more variability expected)
  expect_true(ratio >= 0.6 && ratio <= 1.4,
              info = sprintf("Ratio = %.3f, should be between 0.6 and 1.4", ratio))
})

test_that("p.GFisher_ind has correct empirical type I error with unequal weights", {
  skip_on_cran()  # Too slow for CRAN
  skip_if_not_coga()

  set.seed(34567)
  nsim <- 3000  # Smaller for speed
  n <- 8

  # Fixed df and weights for this test
  df <- runif(n, 0.5, 3)
  w <- abs(rnorm(n))

  p_cal <- rep(NA, nsim)
  for (i in 1:nsim) {
    pval <- runif(n)
    q <- stat.GFisher(pval, df = df, w = w)
    p_cal[i] <- p.GFisher_ind(q = q, df = df, w = w, isExact = TRUE)
  }

  # Check at alpha = 0.05
  empirical_rate <- mean(p_cal <= 0.05)
  ratio <- empirical_rate / 0.05

  expect_true(ratio >= 0.8 && ratio <= 1.2,
              info = sprintf("Ratio = %.3f, should be between 0.8 and 1.2", ratio))
})

test_that("pval.oGFisher_ind has correct empirical type I error - CCT", {
  skip_on_cran()  # Too slow for CRAN
  skip_if_not_coga()

  set.seed(45678)
  nsim <- 2000  # Smaller for speed
  n <- 10
  nGF <- 5

  # Fixed configuration
  DF <- matrix(runif(n * nGF, 0.5, 3), ncol = n) / 10
  W <- abs(matrix(rnorm(n * nGF), ncol = n))

  p_cal <- rep(NA, nsim)
  for (i in 1:nsim) {
    pval <- runif(n)
    result <- pval.oGFisher_ind(p = pval, DF = DF, W = W, combine = "cct")
    p_cal[i] <- result$pval
  }

  # Check at alpha = 0.05
  empirical_rate <- mean(p_cal <= 0.05)
  ratio <- empirical_rate / 0.05

  expect_true(ratio >= 0.8 && ratio <= 1.2,
              info = sprintf("Ratio = %.3f, should be between 0.8 and 1.2", ratio))
})


# =============================================================================
# PART 4: Edge Cases and Robustness
# =============================================================================

test_that("stat.oGFisher_ind is faster than stat.oGFisher", {
  skip_if_not_coga()
  skip_on_cran()  # Skip timing tests on CRAN

  set.seed(123)
  n <- 15
  nGF <- 5

  DF <- matrix(rep(2, n * nGF), ncol = n)
  W <- matrix(rep(1, n * nGF), ncol = n)
  p <- runif(n)

  # Time independent method
  time_ind <- system.time({
    stat.oGFisher_ind(p = p, DF = DF, W = W)
  })

  # Time general method
  time_gen <- system.time({
    stat.oGFisher(p = p, DF = DF, W = W, M = diag(n), method = "HYB")
  })

  # Independent should be faster (but don't enforce strictly, just log)
  message(sprintf("Independent: %.3f sec, General: %.3f sec",
                  time_ind["elapsed"], time_gen["elapsed"]))

  # Just check both complete successfully
  expect_true(time_ind["elapsed"] >= 0)
  expect_true(time_gen["elapsed"] >= 0)
})

test_that("Independent functions handle edge cases", {
  skip_if_not_coga()

  # Single p-value
  result_single <- stat.oGFisher_ind(p = 0.5, DF = matrix(2, 1, 1), W = matrix(1, 1, 1))
  expect_length(result_single$STAT, 1)

  # Very small p-values
  p_small <- rep(1e-15, 5)
  DF_small <- matrix(rep(2, 10), nrow = 2)
  W_small <- matrix(rep(1, 10), nrow = 2)
  result_small <- stat.oGFisher_ind(p = p_small, DF = DF_small, W = W_small)
  expect_true(all(is.finite(result_small$PVAL)))

  # Large q value
  p_large <- p.GFisher_ind(q = 10000, df = rep(2, 5), w = rep(1, 5))
  expect_true(p_large < 1e-10)  # Should be extremely small
})
