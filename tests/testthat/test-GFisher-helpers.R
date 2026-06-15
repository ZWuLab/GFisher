# Test file for GFisher helper functions
# File created: 2025-10-13 by r-rcpp-developer agent

library(testthat)
library(GFisher)

# Test getGFishercoef function
test_that("getGFishercoef computes coefficients for two-sided p-values", {
  set.seed(123)
  n <- 5
  M <- diag(n)
  D <- rep(2, n)

  coeff <- GFisher:::getGFishercoef(D, M, p.type = "two")

  expect_type(coeff, "list")
  expect_named(coeff, c("coeff2", "coeff4", "coeff6", "coeff8"))
  expect_length(coeff$coeff2, n)
  expect_length(coeff$coeff4, n)
  expect_length(coeff$coeff6, n)
  expect_length(coeff$coeff8, n)

  # All coefficients should be numeric
  expect_true(all(is.finite(coeff$coeff2)))
  expect_true(all(is.finite(coeff$coeff4)))
  expect_true(all(is.finite(coeff$coeff6)))
  expect_true(all(is.finite(coeff$coeff8)))
})

test_that("getGFishercoef computes coefficients for one-sided p-values", {
  set.seed(123)
  n <- 5
  M <- diag(n)
  D <- rep(2, n)

  coeff <- GFisher:::getGFishercoef(D, M, p.type = "one")

  expect_type(coeff, "list")
  expect_named(coeff, c("coeff1", "coeff2", "coeff3", "coeff4"))
  expect_length(coeff$coeff1, n)
  expect_length(coeff$coeff2, n)
  expect_length(coeff$coeff3, n)
  expect_length(coeff$coeff4, n)

  # All coefficients should be numeric
  expect_true(all(is.finite(coeff$coeff1)))
  expect_true(all(is.finite(coeff$coeff2)))
  expect_true(all(is.finite(coeff$coeff3)))
  expect_true(all(is.finite(coeff$coeff4)))
})

test_that("getGFishercoef handles single df", {
  n <- 10
  M <- diag(n)
  D <- 2  # Single df

  coeff <- GFisher:::getGFishercoef(D, M, p.type = "two")

  # Should replicate to length n
  expect_length(coeff$coeff2, n)
  expect_true(all(coeff$coeff2 == coeff$coeff2[1]))  # All should be equal
})

test_that("getGFishercoef handles varying df", {
  n <- 5
  M <- diag(n)
  D <- 1:5  # Varying df

  coeff <- GFisher:::getGFishercoef(D, M, p.type = "two")

  expect_length(coeff$coeff2, n)
  # Coefficients should be different for different df
  expect_true(length(unique(coeff$coeff2)) > 1)
})


# Test getGFishercov function
test_that("getGFishercov computes covariance correctly", {
  set.seed(123)
  n <- 5
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  D1 <- rep(2, n)
  D2 <- rep(2, n)
  W1 <- rep(1/n, n)
  W2 <- rep(1/n, n)

  cov_val <- GFisher:::getGFishercov(D1, D2, W1, W2, M, p.type = "two", var.correct = TRUE)

  expect_type(cov_val, "double")
  expect_length(cov_val, 1)
  expect_true(is.finite(cov_val))
})

test_that("getGFishercov gives larger covariance for higher correlation", {
  n <- 5
  D <- rep(2, n)
  W <- rep(1/n, n)

  # Low correlation
  M_low <- matrix(0.1, n, n) + diag(0.9, n, n)
  cov_low <- GFisher:::getGFishercov(D, D, W, W, M_low, p.type = "two")

  # High correlation
  M_high <- matrix(0.7, n, n) + diag(0.3, n, n)
  cov_high <- GFisher:::getGFishercov(D, D, W, W, M_high, p.type = "two")

  expect_true(cov_high > cov_low)
})

test_that("getGFishercov works for one-sided p-values", {
  n <- 5
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  D1 <- rep(2, n)
  D2 <- rep(2, n)
  W1 <- rep(1/n, n)
  W2 <- rep(1/n, n)

  cov_val <- GFisher:::getGFishercov(D1, D2, W1, W2, M, p.type = "one", var.correct = TRUE)

  expect_type(cov_val, "double")
  expect_true(is.finite(cov_val))
})


# Test getGFisherCOR function
test_that("getGFisherCOR computes correlation matrix correctly", {
  set.seed(123)
  n <- 5
  m <- 3  # Number of GFisher tests
  M <- matrix(0.3, n, n) + diag(0.7, n, n)

  DD <- matrix(rep(1:m, each = n), nrow = m, byrow = TRUE)
  W <- matrix(rep(1/n, m * n), nrow = m)

  COR <- GFisher:::getGFisherCOR(DD, W, M, var.correct = TRUE, p.type = "two")

  expect_true(is.matrix(COR))
  expect_equal(dim(COR), c(m, m))

  # Should be symmetric
  expect_equal(COR, t(COR))

  # Diagonal should be 1 (correlation matrix)
  expect_equal(diag(COR), rep(1, m), tolerance = 1e-10)

  # Off-diagonal should be between -1 and 1
  expect_true(all(COR >= -1 & COR <= 1))
})

test_that("getGFisherCOR handles different configurations", {
  n <- 10
  m <- 4
  M <- diag(n)

  # Different df and weights for each test
  DD <- rbind(rep(1, n), rep(2, n), 1:n, rep(3, n))
  W <- rbind(rep(1, n), 1:n, sqrt(1:n), rep(2, n))

  COR <- GFisher:::getGFisherCOR(DD, W, M, p.type = "two")

  expect_equal(dim(COR), c(m, m))
  expect_true(all(is.finite(COR)))
})


# Test getGFisherGM function
test_that("getGFisherGM computes covariance matrix correctly", {
  set.seed(123)
  n <- 5
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  D <- rep(2, n)
  w <- rep(1/n, n)

  GM <- GFisher:::getGFisherGM(D, w, M, p.type = "two")

  expect_true(is.matrix(GM))
  expect_equal(dim(GM), c(n, n))

  # Should be symmetric
  expect_equal(GM, t(GM), tolerance = 1e-10)

  # Should be positive semi-definite (all eigenvalues >= 0 or very small negative)
  eig_vals <- eigen(GM, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig_vals > -1e-10))
})

test_that("getGFisherGM handles varying df and weights", {
  n <- 8
  M <- matrix(0.5, n, n) + diag(0.5, n, n)
  D <- 1:n
  w <- sqrt(1:n) / sum(sqrt(1:n))

  GM <- GFisher:::getGFisherGM(D, w, M, p.type = "two")

  expect_equal(dim(GM), c(n, n))
  expect_true(all(is.finite(GM)))
})

test_that("getGFisherGM works for one-sided p-values", {
  n <- 5
  M <- diag(n)
  D <- rep(2, n)
  w <- rep(1/n, n)

  GM <- GFisher:::getGFisherGM(D, w, M, p.type = "one")

  expect_true(is.matrix(GM))
  expect_equal(dim(GM), c(n, n))
})


# Test getGFisherlam function
test_that("getGFisherlam computes eigenvalues correctly", {
  set.seed(123)
  n <- 5
  M <- matrix(0.3, n, n) + diag(0.7, n, n)
  D <- rep(2, n)
  w <- rep(1/n, n)

  # First compute GM
  GM <- GFisher:::getGFisherGM(D, w, M, p.type = "two")

  # Then compute eigenvalues
  result <- GFisher:::getGFisherlam(D, w, M, GM)

  expect_type(result, "list")
  expect_named(result, "lam")
  expect_true(length(result$lam) > 0)

  # All eigenvalues should be positive (> 1e-10)
  expect_true(all(result$lam > 1e-10))
})

test_that("getGFisherlam handles varying df", {
  n <- 6
  M <- diag(n)
  D <- c(1, 1, 2, 2, 3, 3)  # Mix of df values
  w <- rep(1/n, n)

  GM <- GFisher:::getGFisherGM(D, w, M, p.type = "two")
  result <- GFisher:::getGFisherlam(D, w, M, GM)

  # Should have more eigenvalues when df > 1
  expect_true(length(result$lam) >= n)
})

test_that("getGFisherlam handles large df", {
  n <- 4
  M <- diag(n)
  D <- rep(10, n)  # Large df
  w <- rep(1/n, n)

  GM <- GFisher:::getGFisherGM(D, w, M, p.type = "two")
  result <- GFisher:::getGFisherlam(D, w, M, GM)

  # Should have many more eigenvalues
  expect_true(length(result$lam) > n * 5)
})

test_that("getGFisherlam handles near-singular matrices", {
  n <- 4
  # Create a matrix that needs nearPD correction
  M <- matrix(0.99, n, n) + diag(0.01, n, n)
  D <- rep(2, n)
  w <- rep(1/n, n)

  GM <- GFisher:::getGFisherGM(D, w, M, p.type = "two")

  # Should not error even with near-singular matrix
  expect_error(GFisher:::getGFisherlam(D, w, M, GM), NA)
})


# Integration test: Helper functions work together
test_that("Helper functions integrate correctly", {
  set.seed(777)
  n <- 8
  m <- 3
  M <- matrix(0.4, n, n) + diag(0.6, n, n)

  # Test workflow: compute correlation matrix among multiple GFisher tests
  DD <- rbind(rep(1, n), rep(2, n), rep(3, n))
  W <- rbind(rep(1, n), 1:n, sqrt(1:n))
  W <- W / rowSums(W)  # Normalize

  # This uses getGFishercov internally
  COR <- GFisher:::getGFisherCOR(DD, W, M, p.type = "two")

  expect_equal(dim(COR), c(m, m))
  expect_true(all(is.finite(COR)))
  expect_equal(diag(COR), rep(1, m), tolerance = 1e-10)

  # Test workflow: eigenvalue calculation
  D <- rep(2, n)
  w <- rep(1/n, n)

  GM <- GFisher:::getGFisherGM(D, w, M, p.type = "two")
  lam_result <- GFisher:::getGFisherlam(D, w, M, GM)

  expect_true(length(lam_result$lam) >= n)
  expect_true(all(lam_result$lam > 0))
})
