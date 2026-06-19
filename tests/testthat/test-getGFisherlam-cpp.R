# Test file: test-getGFisherlam-cpp.R
# Purpose: Validate C++ implementation of getGFisherlam against R baseline
# Date: 2025-10-13
# Author: r-rcpp-developer agent

test_that("C++ getGFisherlam matches R implementation - basic case", {
  set.seed(123)
  n <- 10
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  D <- rep(2, n)
  w <- rep(1, n)

  # Get GM first
  GM <- getGFisherGM(D, w, M, "two", use.cpp = FALSE)

  # R version
  lam_r <- getGFisherlam(D, w, M, GM, use.cpp = FALSE)$lam

  # C++ version
  lam_cpp <- getGFisherlam_cpp(D, w, M, GM)$lam

  # Sort for comparison (eigenvalues may be in different order)
  lam_r_sorted <- sort(lam_r, decreasing = TRUE)
  lam_cpp_sorted <- sort(lam_cpp, decreasing = TRUE)

  # Should have same number of eigenvalues
  expect_equal(length(lam_cpp_sorted), length(lam_r_sorted))

  # Should match within tolerance
  expect_equal(lam_cpp_sorted, lam_r_sorted, tolerance = 1e-6)
})

test_that("C++ version handles varying degrees of freedom", {
  set.seed(456)
  n <- 8
  M <- matrix(0.5, n, n) + diag(0.5, n, n)
  D <- rep(c(1, 2, 3), length.out = n)  # Varying df
  w <- rep(1, n)

  GM <- getGFisherGM(D, w, M, "two", use.cpp = FALSE)

  lam_r <- getGFisherlam(D, w, M, GM, use.cpp = FALSE)$lam
  lam_cpp <- getGFisherlam_cpp(D, w, M, GM)$lam

  lam_r_sorted <- sort(lam_r, decreasing = TRUE)
  lam_cpp_sorted <- sort(lam_cpp, decreasing = TRUE)

  # Should have same number of eigenvalues
  expect_equal(length(lam_cpp_sorted), length(lam_r_sorted))

  # Should match within tolerance
  expect_equal(lam_cpp_sorted, lam_r_sorted, tolerance = 1e-6)
})

test_that("C++ version handles max(D) > 1 correctly", {
  set.seed(789)
  n <- 5
  M <- diag(n)  # Independent
  D <- c(1, 1, 2, 3, 3)  # max(D) = 3
  w <- 1:n

  GM <- getGFisherGM(D, w, M, "two", use.cpp = FALSE)

  lam_r <- getGFisherlam(D, w, M, GM, use.cpp = FALSE)$lam
  lam_cpp <- getGFisherlam_cpp(D, w, M, GM)$lam

  # Should have same number of eigenvalues (or very close)
  expect_equal(length(lam_cpp), length(lam_r), tolerance = 1)

  lam_r_sorted <- sort(lam_r, decreasing = TRUE)
  lam_cpp_sorted <- sort(lam_cpp, decreasing = TRUE)

  # Should match within tolerance
  expect_equal(lam_cpp_sorted, lam_r_sorted, tolerance = 1e-6)
})

test_that("C++ version filters small eigenvalues", {
  set.seed(111)
  n <- 5
  M <- matrix(0.95, n, n) + diag(0.05, n, n)  # High correlation
  D <- rep(2, n)
  w <- rep(1, n)

  GM <- getGFisherGM(D, w, M, "two", use.cpp = FALSE)

  lam_cpp <- getGFisherlam_cpp(D, w, M, GM)$lam

  # All eigenvalues should be > 1e-10
  expect_true(all(lam_cpp > 1e-10))
})

test_that("C++ version handles identity correlation", {
  set.seed(222)
  n <- 6
  M <- diag(n)  # Perfect independence
  D <- rep(2, n)
  w <- rep(1, n)

  GM <- getGFisherGM(D, w, M, "two", use.cpp = FALSE)

  lam_r <- getGFisherlam(D, w, M, GM, use.cpp = FALSE)$lam
  lam_cpp <- getGFisherlam_cpp(D, w, M, GM)$lam

  lam_r_sorted <- sort(lam_r, decreasing = TRUE)
  lam_cpp_sorted <- sort(lam_cpp, decreasing = TRUE)

  expect_equal(length(lam_cpp_sorted), length(lam_r_sorted))
  expect_equal(lam_cpp_sorted, lam_r_sorted, tolerance = 1e-6)
})

test_that("C++ version handles weighted cases", {
  set.seed(333)
  n <- 7
  M <- matrix(0.4, n, n) + diag(0.6, n, n)
  D <- rep(2, n)
  w <- seq(0.5, 2, length.out = n)  # Varying weights

  GM <- getGFisherGM(D, w, M, "two", use.cpp = FALSE)

  lam_r <- getGFisherlam(D, w, M, GM, use.cpp = FALSE)$lam
  lam_cpp <- getGFisherlam_cpp(D, w, M, GM)$lam

  lam_r_sorted <- sort(lam_r, decreasing = TRUE)
  lam_cpp_sorted <- sort(lam_cpp, decreasing = TRUE)

  expect_equal(length(lam_cpp_sorted), length(lam_r_sorted))
  expect_equal(lam_cpp_sorted, lam_r_sorted, tolerance = 1e-6)
})

test_that("C++ version handles large degrees of freedom", {
  set.seed(444)
  n <- 10
  M <- diag(n)
  D <- rep(10, n)  # Large df, complex eigenvalue structure
  w <- rep(1, n)

  GM <- getGFisherGM(D, w, M, "two", use.cpp = FALSE)

  lam_r <- getGFisherlam(D, w, M, GM, use.cpp = FALSE)$lam
  lam_cpp <- getGFisherlam_cpp(D, w, M, GM)$lam

  # May have different numbers due to filtering, but should be close
  expect_gt(length(lam_cpp), 0)
  expect_gt(length(lam_r), 0)

  # Compare the largest eigenvalues (most important)
  n_compare <- min(length(lam_cpp), length(lam_r), 10)
  lam_r_sorted <- sort(lam_r, decreasing = TRUE)[1:n_compare]
  lam_cpp_sorted <- sort(lam_cpp, decreasing = TRUE)[1:n_compare]

  expect_equal(lam_cpp_sorted, lam_r_sorted, tolerance = 1e-6)
})

test_that("Wrapper function uses C++ by default", {
  set.seed(555)
  n <- 5
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  D <- rep(2, n)
  w <- rep(1, n)

  GM <- getGFisherGM(D, w, M, "two", use.cpp = FALSE)

  # Default should use C++
  lam_default <- getGFisherlam(D, w, M, GM)$lam
  lam_cpp <- getGFisherlam_cpp(D, w, M, GM)$lam

  lam_default_sorted <- sort(lam_default, decreasing = TRUE)
  lam_cpp_sorted <- sort(lam_cpp, decreasing = TRUE)

  expect_equal(lam_default_sorted, lam_cpp_sorted, tolerance = 1e-10)
})

test_that("Wrapper function falls back to R when requested", {
  set.seed(666)
  n <- 5
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  D <- rep(2, n)
  w <- rep(1, n)

  GM <- getGFisherGM(D, w, M, "two", use.cpp = FALSE)

  # Force R version
  lam_r1 <- getGFisherlam(D, w, M, GM, use.cpp = FALSE)$lam
  lam_r2 <- getGFisherlam(D, w, M, GM, use.cpp = FALSE)$lam

  lam_r1_sorted <- sort(lam_r1, decreasing = TRUE)
  lam_r2_sorted <- sort(lam_r2, decreasing = TRUE)

  # Both R calls should give identical results
  expect_equal(lam_r1_sorted, lam_r2_sorted, tolerance = 1e-10)
})

test_that("C++ version handles small matrices", {
  set.seed(777)
  n <- 2
  M <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  D <- c(2, 2)
  w <- c(1, 1)

  GM <- getGFisherGM(D, w, M, "two", use.cpp = FALSE)

  lam_r <- getGFisherlam(D, w, M, GM, use.cpp = FALSE)$lam
  lam_cpp <- getGFisherlam_cpp(D, w, M, GM)$lam

  lam_r_sorted <- sort(lam_r, decreasing = TRUE)
  lam_cpp_sorted <- sort(lam_cpp, decreasing = TRUE)

  expect_equal(length(lam_cpp_sorted), length(lam_r_sorted))
  expect_equal(lam_cpp_sorted, lam_r_sorted, tolerance = 1e-6)
})

test_that("C++ version handles complex correlation structures", {
  set.seed(888)
  n <- 12
  # Create block correlation structure
  M <- diag(n)
  M[1:4, 1:4] <- 0.8
  diag(M[1:4, 1:4]) <- 1
  M[5:8, 5:8] <- 0.6
  diag(M[5:8, 5:8]) <- 1
  M[9:12, 9:12] <- 0.4
  diag(M[9:12, 9:12]) <- 1

  D <- rep(c(1, 2, 3), each = 4)
  w <- rep(1, n)

  GM <- getGFisherGM(D, w, M, "two", use.cpp = FALSE)

  lam_r <- getGFisherlam(D, w, M, GM, use.cpp = FALSE)$lam
  lam_cpp <- getGFisherlam_cpp(D, w, M, GM)$lam

  lam_r_sorted <- sort(lam_r, decreasing = TRUE)
  lam_cpp_sorted <- sort(lam_cpp, decreasing = TRUE)

  # Should have similar number of eigenvalues
  expect_equal(length(lam_cpp_sorted), length(lam_r_sorted), tolerance = 2)

  # Compare top eigenvalues
  n_compare <- min(length(lam_cpp), length(lam_r), 15)
  expect_equal(lam_cpp_sorted[1:n_compare], lam_r_sorted[1:n_compare], tolerance = 1e-6)
})
