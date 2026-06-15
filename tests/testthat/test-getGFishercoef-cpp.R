## File: test-getGFishercoef-cpp.R
## Purpose: Test C++ implementation of getGFishercoef against pure R version
## Date: 2025-10-13
## Author: r-rcpp-developer agent

test_that("getGFishercoef_cpp matches R for two-sided p-values with single df", {
  # Setup: simple correlation matrix
  M <- matrix(c(1, 0.5, 0.3,
                0.5, 1, 0.4,
                0.3, 0.4, 1), nrow = 3)
  D <- 1  # Single df

  # Pure R version
  result_r <- getGFishercoef(D, M, p.type = "two", use.cpp = FALSE)

  # C++ version
  result_cpp <- getGFishercoef_cpp(D, M, p_type = "two")

  # Check that all coefficients match
  expect_equal(result_cpp$coeff2, result_r$coeff2, tolerance = 1e-4)
  expect_equal(result_cpp$coeff4, result_r$coeff4, tolerance = 1e-4)
  expect_equal(result_cpp$coeff6, result_r$coeff6, tolerance = 1e-4)
  expect_equal(result_cpp$coeff8, result_r$coeff8, tolerance = 1e-4)

  # Check that single df is replicated correctly
  expect_equal(length(result_cpp$coeff2), 3)
  expect_equal(length(result_cpp$coeff4), 3)
  expect_equal(length(result_cpp$coeff6), 3)
  expect_equal(length(result_cpp$coeff8), 3)
})


test_that("getGFishercoef_cpp matches R for two-sided p-values with varying df", {
  # Setup: simple correlation matrix
  M <- matrix(c(1, 0.6, 0.2,
                0.6, 1, 0.5,
                0.2, 0.5, 1), nrow = 3)
  D <- c(1, 2, 3)  # Varying df

  # Pure R version
  result_r <- getGFishercoef(D, M, p.type = "two", use.cpp = FALSE)

  # C++ version
  result_cpp <- getGFishercoef_cpp(D, M, p_type = "two")

  # Check that all coefficients match
  expect_equal(result_cpp$coeff2, result_r$coeff2, tolerance = 1e-4)
  expect_equal(result_cpp$coeff4, result_r$coeff4, tolerance = 1e-4)
  expect_equal(result_cpp$coeff6, result_r$coeff6, tolerance = 1e-4)
  expect_equal(result_cpp$coeff8, result_r$coeff8, tolerance = 1e-4)
})


test_that("getGFishercoef_cpp matches R for one-sided p-values with single df", {
  # Setup: simple correlation matrix
  M <- matrix(c(1, 0.4, 0.3, 0.2,
                0.4, 1, 0.5, 0.1,
                0.3, 0.5, 1, 0.6,
                0.2, 0.1, 0.6, 1), nrow = 4)
  D <- 1  # Single df

  # Pure R version
  result_r <- getGFishercoef(D, M, p.type = "one", use.cpp = FALSE)

  # C++ version
  result_cpp <- getGFishercoef_cpp(D, M, p_type = "one")

  # Check that all coefficients match
  expect_equal(result_cpp$coeff1, result_r$coeff1, tolerance = 1e-4)
  expect_equal(result_cpp$coeff2, result_r$coeff2, tolerance = 1e-4)
  expect_equal(result_cpp$coeff3, result_r$coeff3, tolerance = 1e-4)
  expect_equal(result_cpp$coeff4, result_r$coeff4, tolerance = 1e-4)

  # Check that single df is replicated correctly
  expect_equal(length(result_cpp$coeff1), 4)
  expect_equal(length(result_cpp$coeff2), 4)
  expect_equal(length(result_cpp$coeff3), 4)
  expect_equal(length(result_cpp$coeff4), 4)
})


test_that("getGFishercoef_cpp matches R for one-sided p-values with varying df", {
  # Setup: simple correlation matrix
  M <- matrix(c(1, 0.7, 0.3,
                0.7, 1, 0.4,
                0.3, 0.4, 1), nrow = 3)
  D <- c(1, 2, 4)  # Varying df

  # Pure R version
  result_r <- getGFishercoef(D, M, p.type = "one", use.cpp = FALSE)

  # C++ version
  result_cpp <- getGFishercoef_cpp(D, M, p_type = "one")

  # Check that all coefficients match
  expect_equal(result_cpp$coeff1, result_r$coeff1, tolerance = 1e-4)
  expect_equal(result_cpp$coeff2, result_r$coeff2, tolerance = 1e-4)
  expect_equal(result_cpp$coeff3, result_r$coeff3, tolerance = 1e-4)
  expect_equal(result_cpp$coeff4, result_r$coeff4, tolerance = 1e-4)
})


test_that("R wrapper getGFishercoef uses C++ by default", {
  # Setup
  M <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  D <- 1

  # Call wrapper with default use.cpp = TRUE
  result_default <- getGFishercoef(D, M, p.type = "two")

  # Call C++ directly
  result_cpp <- getGFishercoef_cpp(D, M, p_type = "two")

  # Should be identical (no fallback warning)
  expect_equal(result_default$coeff2, result_cpp$coeff2, tolerance = 1e-10)
  expect_equal(result_default$coeff4, result_cpp$coeff4, tolerance = 1e-10)
  expect_equal(result_default$coeff6, result_cpp$coeff6, tolerance = 1e-10)
  expect_equal(result_default$coeff8, result_cpp$coeff8, tolerance = 1e-10)
})


test_that("R wrapper getGFishercoef can force pure R implementation", {
  # Setup
  M <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  D <- c(1, 2)

  # Force pure R
  result_r <- getGFishercoef(D, M, p.type = "two", use.cpp = FALSE)

  # Call C++ directly
  result_cpp <- getGFishercoef_cpp(D, M, p_type = "two")

  # Should match within integration tolerance
  expect_equal(result_r$coeff2, result_cpp$coeff2, tolerance = 1e-4)
  expect_equal(result_r$coeff4, result_cpp$coeff4, tolerance = 1e-4)
  expect_equal(result_r$coeff6, result_cpp$coeff6, tolerance = 1e-4)
  expect_equal(result_r$coeff8, result_cpp$coeff8, tolerance = 1e-4)
})


test_that("getGFishercoef_cpp handles edge cases correctly", {
  # Test with identity matrix (no correlation)
  M_identity <- diag(3)
  D <- c(1, 2, 3)

  result_cpp <- getGFishercoef_cpp(D, M_identity, p_type = "two")
  result_r <- getGFishercoef(D, M_identity, p.type = "two", use.cpp = FALSE)

  expect_equal(result_cpp$coeff2, result_r$coeff2, tolerance = 1e-4)
  expect_equal(result_cpp$coeff4, result_r$coeff4, tolerance = 1e-4)
  expect_equal(result_cpp$coeff6, result_r$coeff6, tolerance = 1e-4)
  expect_equal(result_cpp$coeff8, result_r$coeff8, tolerance = 1e-4)
})


test_that("getGFishercoef_cpp works with larger matrices", {
  # Test with larger matrix (n=10)
  set.seed(123)
  n <- 10
  # Create a valid correlation matrix
  A <- matrix(rnorm(n * n), n, n)
  M <- cov2cor(A %*% t(A))

  D <- rep(1, n)

  result_cpp <- getGFishercoef_cpp(D, M, p_type = "two")
  result_r <- getGFishercoef(D, M, p.type = "two", use.cpp = FALSE)

  # Check that results match
  expect_equal(result_cpp$coeff2, result_r$coeff2, tolerance = 1e-4)
  expect_equal(result_cpp$coeff4, result_r$coeff4, tolerance = 1e-4)
  expect_equal(result_cpp$coeff6, result_r$coeff6, tolerance = 1e-4)
  expect_equal(result_cpp$coeff8, result_r$coeff8, tolerance = 1e-4)

  # Check correct dimensions
  expect_equal(length(result_cpp$coeff2), n)
  expect_equal(length(result_cpp$coeff4), n)
  expect_equal(length(result_cpp$coeff6), n)
  expect_equal(length(result_cpp$coeff8), n)
})


test_that("getGFishercoef_cpp error handling works", {
  # Test with mismatched dimensions
  M <- diag(3)
  D <- c(1, 2)  # Wrong length

  expect_error(getGFishercoef_cpp(D, M, p_type = "two"),
               "Length of D must be 1 or match dimension of M")

  # Test with invalid p_type
  D <- c(1, 1, 1)
  expect_error(getGFishercoef_cpp(D, M, p_type = "invalid"),
               "p_type must be 'two' or 'one'")
})


test_that("getGFishercoef_cpp one-sided with large df", {
  # Test one-sided with larger degrees of freedom
  M <- matrix(c(1, 0.5, 0.3,
                0.5, 1, 0.6,
                0.3, 0.6, 1), nrow = 3)
  D <- c(5, 10, 15)

  result_cpp <- getGFishercoef_cpp(D, M, p_type = "one")
  result_r <- getGFishercoef(D, M, p.type = "one", use.cpp = FALSE)

  expect_equal(result_cpp$coeff1, result_r$coeff1, tolerance = 1e-4)
  expect_equal(result_cpp$coeff2, result_r$coeff2, tolerance = 1e-4)
  expect_equal(result_cpp$coeff3, result_r$coeff3, tolerance = 1e-4)
  expect_equal(result_cpp$coeff4, result_r$coeff4, tolerance = 1e-4)
})
