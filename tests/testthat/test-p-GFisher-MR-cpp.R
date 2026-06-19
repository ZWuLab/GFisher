# File: test-p-GFisher-MR-cpp.R
# Purpose: Test C++ implementation of p.GFisher method="MR" simulation
# Date: 2025-10-13
# Author: r-rcpp-developer agent

test_that("C++ p_GFisher_MR matches R implementation with same seed", {
  skip_if_not(requireNamespace("GFisher", quietly = TRUE))

  set.seed(123)
  n <- 10
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  df <- rep(2, n)
  w <- rep(1, n)

  # Generate test statistic
  pval <- runif(n)
  q <- stat.GFisher(pval, df, w)

  # R version
  set.seed(456)
  p_r <- p.GFisher(q, df, w, M, method = "MR", nsim = 1e4, seed = 456, use.cpp = FALSE)

  # C++ version (seed must be set in R, pass 0 to C++)
  set.seed(456)
  p_cpp <- p_GFisher_MR_cpp(q, df, w, M, "two", 1e4, 0)

  # Should match exactly (both use same seed)
  expect_equal(p_cpp, p_r, tolerance = 1e-10)
})

test_that("C++ version handles varying df", {
  skip_if_not(requireNamespace("GFisher", quietly = TRUE))

  set.seed(789)
  n <- 8
  M <- matrix(0.5, n, n) + diag(0.5, n, n)
  df <- rep(c(1, 2, 3), length.out = n)
  w <- rep(1, n)

  pval <- runif(n)
  q <- stat.GFisher(pval, df, w)

  set.seed(111)
  p_r <- p.GFisher(q, df, w, M, method = "MR", nsim = 5e3, seed = 111, use.cpp = FALSE)
  set.seed(111)
  p_cpp <- p_GFisher_MR_cpp(q, df, w, M, "two", 5e3, 0)

  expect_equal(p_cpp, p_r, tolerance = 1e-10)
})

test_that("C++ version handles one-sided p-values", {
  skip_if_not(requireNamespace("GFisher", quietly = TRUE))

  set.seed(222)
  n <- 5
  M <- diag(n)
  df <- rep(2, n)
  w <- rep(1, n)

  pval <- runif(n)
  q <- stat.GFisher(pval, df, w)

  set.seed(333)
  p_r <- p.GFisher(q, df, w, M, p.type = "one", method = "MR", nsim = 5e3, seed = 333, use.cpp = FALSE)
  set.seed(333)
  p_cpp <- p_GFisher_MR_cpp(q, df, w, M, "one", 5e3, 0)

  expect_equal(p_cpp, p_r, tolerance = 1e-10)
})

test_that("C++ version handles different weights", {
  skip_if_not(requireNamespace("GFisher", quietly = TRUE))

  set.seed(444)
  n <- 10
  M <- matrix(0.4, n, n) + diag(0.6, n, n)
  df <- rep(2, n)
  w <- 1:n

  pval <- runif(n)
  q <- stat.GFisher(pval, df, w)

  set.seed(555)
  p_r <- p.GFisher(q, df, w, M, method = "MR", nsim = 5e3, seed = 555, use.cpp = FALSE)
  set.seed(555)
  p_cpp <- p_GFisher_MR_cpp(q, df, w, M, "two", 5e3, 0)

  expect_equal(p_cpp, p_r, tolerance = 1e-10)
})

test_that("C++ version produces valid p-values", {
  skip_if_not(requireNamespace("GFisher", quietly = TRUE))

  set.seed(666)
  n <- 10
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  df <- rep(2, n)
  w <- rep(1, n)

  pval <- runif(n)
  q <- stat.GFisher(pval, df, w)

  set.seed(777)
  p_cpp <- p_GFisher_MR_cpp(q, df, w, M, "two", 1e4, 0)

  # P-value should be between 0 and 1
  expect_true(p_cpp >= 0 && p_cpp <= 1)
})

test_that("Wrapper function with use.cpp=TRUE matches direct C++ call", {
  skip_if_not(requireNamespace("GFisher", quietly = TRUE))

  set.seed(888)
  n <- 10
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  df <- rep(2, n)
  w <- rep(1, n)

  pval <- runif(n)
  q <- stat.GFisher(pval, df, w)

  # Wrapper with use.cpp=TRUE
  set.seed(999)
  p_wrapper <- p.GFisher(q, df, w, M, method = "MR", nsim = 1e4, seed = 999, use.cpp = TRUE)

  # Direct C++ call
  set.seed(999)
  p_cpp <- p_GFisher_MR_cpp(q, df, w, M, "two", 1e4, 0)

  expect_equal(p_wrapper, p_cpp, tolerance = 1e-10)
})

test_that("C++ version handles df=2 special case", {
  skip_if_not(requireNamespace("GFisher", quietly = TRUE))

  set.seed(100)
  n <- 15
  M <- matrix(0.2, n, n) + diag(0.8, n, n)
  df <- rep(2, n)  # Special case: df=2 uses -2*log(p) transformation
  w <- rep(1, n)

  pval <- runif(n)
  q <- stat.GFisher(pval, df, w)

  set.seed(200)
  p_r <- p.GFisher(q, df, w, M, method = "MR", nsim = 5e3, seed = 200, use.cpp = FALSE)
  set.seed(200)
  p_cpp <- p_GFisher_MR_cpp(q, df, w, M, "two", 5e3, 0)

  expect_equal(p_cpp, p_r, tolerance = 1e-10)
})

test_that("C++ version handles single df value with replication", {
  skip_if_not(requireNamespace("GFisher", quietly = TRUE))

  set.seed(300)
  n <- 10
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  df <- rep(3, n)  # All same non-2 df
  w <- rep(1, n)

  pval <- runif(n)
  q <- stat.GFisher(pval, df, w)

  set.seed(400)
  p_r <- p.GFisher(q, df, w, M, method = "MR", nsim = 5e3, seed = 400, use.cpp = FALSE)
  set.seed(400)
  p_cpp <- p_GFisher_MR_cpp(q, df, w, M, "two", 5e3, 0)

  expect_equal(p_cpp, p_r, tolerance = 1e-10)
})

test_that("C++ version works with complex correlation structure", {
  skip_if_not(requireNamespace("GFisher", quietly = TRUE))

  set.seed(500)
  n <- 12
  # Create a more complex correlation structure
  M <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      M[i, j] <- 0.7^abs(i - j)  # AR(1) structure
    }
  }

  df <- rep(c(1, 2, 3), each = 4)
  w <- 1:n

  pval <- runif(n)
  q <- stat.GFisher(pval, df, w)

  set.seed(600)
  p_r <- p.GFisher(q, df, w, M, method = "MR", nsim = 5e3, seed = 600, use.cpp = FALSE)
  set.seed(600)
  p_cpp <- p_GFisher_MR_cpp(q, df, w, M, "two", 5e3, 0)

  expect_equal(p_cpp, p_r, tolerance = 1e-10)
})

test_that("C++ version handles small p-values correctly", {
  skip_if_not(requireNamespace("GFisher", quietly = TRUE))

  set.seed(700)
  n <- 10
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  df <- rep(2, n)
  w <- rep(1, n)

  # Generate small p-values to test tail behavior
  pval <- c(0.001, 0.002, 0.005, runif(7))
  q <- stat.GFisher(pval, df, w)

  set.seed(800)
  p_r <- p.GFisher(q, df, w, M, method = "MR", nsim = 5e3, seed = 800, use.cpp = FALSE)
  set.seed(800)
  p_cpp <- p_GFisher_MR_cpp(q, df, w, M, "two", 5e3, 0)

  expect_equal(p_cpp, p_r, tolerance = 1e-10)
})

test_that("C++ version handles near-singular correlation matrix", {
  skip_if_not(requireNamespace("GFisher", quietly = TRUE))

  set.seed(900)
  n <- 8
  # Create a nearly singular correlation matrix
  M <- matrix(0.9, n, n) + diag(0.1, n, n)

  df <- rep(2, n)
  w <- rep(1, n)

  pval <- runif(n)
  q <- stat.GFisher(pval, df, w)

  # Both should handle this with nearPD correction
  set.seed(1000)
  p_r <- p.GFisher(q, df, w, M, method = "MR", nsim = 5e3, seed = 1000, use.cpp = FALSE)
  set.seed(1000)
  p_cpp_wrapper <- p.GFisher(q, df, w, M, method = "MR", nsim = 5e3, seed = 1000, use.cpp = TRUE)

  # Should match exactly
  expect_equal(p_cpp_wrapper, p_r, tolerance = 1e-10)
})
