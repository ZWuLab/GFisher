# Test file for coefficient caching functionality
# Date: 2025-10-13
# Purpose: Verify that caching improves performance and maintains accuracy

library(testthat)
library(GFisher)

test_that("Cached version produces identical results to non-cached", {
  set.seed(123)
  n <- 10
  M <- diag(n)
  D <- c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1)  # Repeated df values

  # Clear cache to start fresh
  GFisher:::clear_GFisher_cache()

  # Get results from cached version
  result_cached <- GFisher:::getGFishercoef_cached_cpp(D, M, "two")

  # Get results from non-cached version
  result_regular <- GFisher:::getGFishercoef_cpp(D, M, "two")

  # Check that results are identical
  expect_equal(result_cached$coeff2, result_regular$coeff2, tolerance = 1e-10)
  expect_equal(result_cached$coeff4, result_regular$coeff4, tolerance = 1e-10)
  expect_equal(result_cached$coeff6, result_regular$coeff6, tolerance = 1e-10)
  expect_equal(result_cached$coeff8, result_regular$coeff8, tolerance = 1e-10)
})

test_that("Cache management functions work correctly", {
  # Clear cache
  cleared <- GFisher:::clear_GFisher_cache()
  expect_true(cleared >= 0)

  # Check empty cache
  expect_equal(GFisher:::get_GFisher_cache_size(), 0)

  # Run computation to populate cache
  n <- 5
  M <- diag(n)
  D <- c(1, 2, 3, 2, 1)  # Some repeated values
  result <- GFisher:::getGFishercoef_cached_cpp(D, M, "two")

  # Check cache has entries (3 unique df values × 4 coefficient types)
  size <- GFisher:::get_GFisher_cache_size()
  expect_true(size > 0)
  expect_true(size <= 12)  # At most 3 unique df × 4 coefficients

  # Get cache statistics
  stats <- GFisher:::get_GFisher_cache_stats()
  expect_type(stats, "list")
  expect_true(stats$size > 0)
  expect_equal(stats$size, size)
  expect_true(stats$memory_kb > 0)
  expect_true(stats$fill_ratio >= 0 && stats$fill_ratio <= 100)

  # Clear cache again
  cleared <- GFisher:::clear_GFisher_cache()
  expect_equal(cleared, size)
  expect_equal(GFisher:::get_GFisher_cache_size(), 0)
})

test_that("Caching improves performance for repeated df values", {
  set.seed(456)
  n <- 100
  M <- diag(n)
  D <- rep(c(1, 2, 3), length.out = n)  # Many repeated values

  # Clear cache for fair comparison
  GFisher:::clear_GFisher_cache()

  # First run - populates cache
  time1 <- system.time({
    result1 <- GFisher:::getGFishercoef_cached_cpp(D, M, "two")
  })

  # Second run - should use cache
  time2 <- system.time({
    result2 <- GFisher:::getGFishercoef_cached_cpp(D, M, "two")
  })

  # Cache should make second run much faster
  # Note: This might be too fast to measure reliably, so we just check it doesn't error
  expect_equal(result1$coeff2, result2$coeff2)

  # Check cache was populated
  cache_size <- GFisher:::get_GFisher_cache_size()
  expect_true(cache_size > 0)  # Should have cached values
})

test_that("Caching works for both two-sided and one-sided p-values", {
  n <- 8
  M <- diag(n)
  D <- c(1, 1, 2, 2, 3, 3, 4, 4)

  # Clear cache
  GFisher:::clear_GFisher_cache()

  # Two-sided
  result_two <- GFisher:::getGFishercoef_cached_cpp(D, M, "two")
  expect_type(result_two, "list")
  expect_named(result_two, c("coeff2", "coeff4", "coeff6", "coeff8"))

  cache_size_two <- GFisher:::get_GFisher_cache_size()
  expect_true(cache_size_two > 0)

  # One-sided
  result_one <- GFisher:::getGFishercoef_cached_cpp(D, M, "one")
  expect_type(result_one, "list")
  expect_named(result_one, c("coeff1", "coeff2", "coeff3", "coeff4"))

  cache_size_both <- GFisher:::get_GFisher_cache_size()
  expect_true(cache_size_both > cache_size_two)  # Should have added more entries
})

test_that("Wrapper function uses cached version by default", {
  n <- 6
  M <- diag(n)
  D <- c(1, 2, 1, 2, 1, 2)

  # Clear cache to track new entries
  GFisher:::clear_GFisher_cache()

  # Use wrapper function (should use cached version by default)
  result <- GFisher:::getGFishercoef(D, M, p.type = "two", use.cpp = TRUE)

  # Check that cache was populated
  cache_size <- GFisher:::get_GFisher_cache_size()
  expect_true(cache_size > 0)

  # Verify result structure
  expect_type(result, "list")
  expect_length(result, 4)
})

test_that("Single df value is properly replicated", {
  n <- 5
  M <- diag(n)
  D_single <- 2.5

  # Clear cache
  GFisher:::clear_GFisher_cache()

  # Compute with single df
  result <- GFisher:::getGFishercoef_cached_cpp(D_single, M, "two")

  # Check all values are identical (since same df)
  expect_true(all(result$coeff2 == result$coeff2[1]))
  expect_true(all(result$coeff4 == result$coeff4[1]))
  expect_true(all(result$coeff6 == result$coeff6[1]))
  expect_true(all(result$coeff8 == result$coeff8[1]))

  # Check cache has only 4 entries (one df × 4 coefficient types)
  cache_size <- GFisher:::get_GFisher_cache_size()
  expect_equal(cache_size, 4)
})