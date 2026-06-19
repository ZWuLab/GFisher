# Test file: test-getGFisherGM-cpp.R
# Purpose: Validate C++ implementation of getGFisherGM against R baseline
# Date: 2025-10-13
# Author: r-rcpp-developer agent

test_that("C++ getGFisherGM matches R implementation for two-sided p-values", {
  set.seed(123)
  n <- 10
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  D <- rep(2, n)
  w <- rep(1, n)

  # R version (force use.cpp = FALSE)
  GM_r <- getGFisherGM(D, w, M, p.type = "two", use.cpp = FALSE)

  # C++ version (direct call)
  GM_cpp <- getGFisherGM_cpp(D, w, M, "two")

  # Should match within numerical tolerance
  expect_equal(GM_cpp, GM_r, tolerance = 1e-4)
})

test_that("C++ getGFisherGM matches R implementation for one-sided p-values", {
  set.seed(456)
  n <- 8
  M <- matrix(0.5, n, n) + diag(0.5, n, n)
  D <- rep(3, n)
  w <- rep(1, n)

  # R version
  GM_r <- getGFisherGM(D, w, M, p.type = "one", use.cpp = FALSE)

  # C++ version
  GM_cpp <- getGFisherGM_cpp(D, w, M, "one")

  # Should match within numerical tolerance
  expect_equal(GM_cpp, GM_r, tolerance = 1e-4)
})

test_that("C++ version handles different degrees of freedom", {
  set.seed(789)
  n <- 5
  M <- matrix(0.4, n, n) + diag(0.6, n, n)
  D <- 1:5  # Different df for each
  w <- rep(1, n)

  # R version
  GM_r <- getGFisherGM(D, w, M, p.type = "two", use.cpp = FALSE)

  # C++ version
  GM_cpp <- getGFisherGM_cpp(D, w, M, "two")

  expect_equal(GM_cpp, GM_r, tolerance = 1e-4)
})

test_that("C++ version handles different weights", {
  set.seed(111)
  n <- 6
  M <- diag(n)  # Independent case
  D <- rep(2, n)
  w <- 1:n  # Different weights

  # R version
  GM_r <- getGFisherGM(D, w, M, p.type = "two", use.cpp = FALSE)

  # C++ version
  GM_cpp <- getGFisherGM_cpp(D, w, M, "two")

  expect_equal(GM_cpp, GM_r, tolerance = 1e-4)
})

test_that("C++ version handles non-identical correlation matrices", {
  set.seed(222)
  n <- 7
  # Generate a valid correlation matrix
  M <- matrix(runif(n * n, -0.5, 0.5), n, n)
  M <- (M + t(M)) / 2  # Make symmetric
  diag(M) <- 1
  # Ensure positive definite
  M <- Matrix::nearPD(M, corr = TRUE)$mat
  M <- as.matrix(M)

  D <- rep(2, n)
  w <- rep(1, n)

  # R version
  GM_r <- getGFisherGM(D, w, M, p.type = "two", use.cpp = FALSE)

  # C++ version
  GM_cpp <- getGFisherGM_cpp(D, w, M, "two")

  expect_equal(GM_cpp, GM_r, tolerance = 1e-5)  # Slightly larger tolerance for complex matrix
})

test_that("Wrapper function uses C++ by default", {
  set.seed(333)
  n <- 5
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  D <- rep(2, n)
  w <- rep(1, n)

  # Wrapper should use C++ by default
  GM_wrapper <- getGFisherGM(D, w, M, p.type = "two")

  # Direct C++ call
  GM_cpp <- getGFisherGM_cpp(D, w, M, "two")

  expect_equal(GM_wrapper, GM_cpp)
})

test_that("Wrapper fallback works when use.cpp = FALSE", {
  set.seed(444)
  n <- 5
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  D <- rep(2, n)
  w <- rep(1, n)

  # Force R implementation
  GM_r <- getGFisherGM(D, w, M, p.type = "two", use.cpp = FALSE)

  # Should still return valid matrix
  expect_true(is.matrix(GM_r))
  expect_equal(dim(GM_r), c(n, n))
  expect_equal(GM_r, t(GM_r))  # Should be symmetric
})

test_that("C++ version handles edge case: single df value", {
  set.seed(555)
  n <- 4
  M <- matrix(0.2, n, n) + diag(0.8, n, n)
  D <- rep(1, n)  # Minimum df
  w <- rep(1, n)

  # R version
  GM_r <- getGFisherGM(D, w, M, p.type = "two", use.cpp = FALSE)

  # C++ version
  GM_cpp <- getGFisherGM_cpp(D, w, M, "two")

  expect_equal(GM_cpp, GM_r, tolerance = 1e-4)
})

test_that("C++ version handles large df values", {
  set.seed(666)
  n <- 5
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  D <- rep(10, n)  # Larger df
  w <- rep(1, n)

  # R version
  GM_r <- getGFisherGM(D, w, M, p.type = "two", use.cpp = FALSE)

  # C++ version
  GM_cpp <- getGFisherGM_cpp(D, w, M, "two")

  # Large df requires slightly more tolerance due to integration complexity
  expect_equal(GM_cpp, GM_r, tolerance = 3e-4)
})

test_that("C++ version produces symmetric matrix", {
  set.seed(777)
  n <- 6
  M <- matrix(0.4, n, n) + diag(0.6, n, n)
  D <- rep(2, n)
  w <- rep(1, n)

  # C++ version
  GM_cpp <- getGFisherGM_cpp(D, w, M, "two")

  # Check symmetry
  expect_equal(GM_cpp, t(GM_cpp), tolerance = 1e-10)
})

test_that("C++ version produces correct diagonal elements", {
  set.seed(888)
  n <- 5
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  D <- rep(2, n)
  w <- rep(1, n)

  # C++ version
  GM_cpp <- getGFisherGM_cpp(D, w, M, "two")

  # Diagonal should be sqrt(2*D)*w squared
  expected_diag <- (sqrt(2 * D) * w)^2
  expect_equal(diag(GM_cpp), expected_diag, tolerance = 1e-6)
})
