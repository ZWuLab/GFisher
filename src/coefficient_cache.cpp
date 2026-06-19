// File: coefficient_cache.cpp
// Purpose: Cached version of getGFishercoef_cpp with memory management
// Date: 2025-10-13
// Author: Claude Code

#include <RcppArmadillo.h>
#include <unordered_map>
#include <string>
#include <sstream>
#include <mutex>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Global cache for coefficient results
// Key format: "df_p_type_coeff_type" (e.g., "2.5_two_2" for df=2.5, two-sided, coeff2)
static std::unordered_map<std::string, double> coefficient_cache;
static std::mutex cache_mutex;  // Thread safety for cache access

// Maximum cache size (prevent unbounded memory growth)
static const size_t MAX_CACHE_SIZE = 10000;  // ~80KB max memory usage

// Forward declarations of integration functions from integration.cpp
arma::vec compute_coefficients_two_cpp(arma::vec D, int coeff_type);
arma::vec compute_coefficients_one_cpp(arma::vec D, int coeff_type);

// Helper function to create cache key
std::string make_cache_key(double df, const std::string& p_type, int coeff_type) {
  std::ostringstream oss;
  oss << df << "_" << p_type << "_" << coeff_type;
  return oss.str();
}

// Helper function to compute single coefficient with caching
double get_cached_coefficient(double df, const std::string& p_type, int coeff_type) {
  // Create cache key
  std::string key = make_cache_key(df, p_type, coeff_type);

  // Thread-safe cache lookup
  {
    std::lock_guard<std::mutex> lock(cache_mutex);

    // Check if value exists in cache
    auto it = coefficient_cache.find(key);
    if (it != coefficient_cache.end()) {
      return it->second;  // Cache hit!
    }
  }

  // Cache miss - compute the value
  arma::vec D_single = {df};
  arma::vec result;

  if (p_type == "two") {
    result = compute_coefficients_two_cpp(D_single, coeff_type);
  } else {
    result = compute_coefficients_one_cpp(D_single, coeff_type);
  }

  double value = result(0);

  // Store in cache (thread-safe)
  {
    std::lock_guard<std::mutex> lock(cache_mutex);

    // Check cache size and clear if necessary
    if (coefficient_cache.size() >= MAX_CACHE_SIZE) {
      // Clear cache when it gets too large (simple strategy)
      // In production, could use LRU eviction instead
      coefficient_cache.clear();
      Rcpp::Rcout << "GFisher coefficient cache cleared (size limit reached)" << std::endl;
    }

    coefficient_cache[key] = value;
  }

  return value;
}

//' Cached Compute Hermite Polynomial Coefficients for GFisher
//'
//' Fast C++ implementation with caching for repeated df values.
//' This can significantly speed up analyses where the same df values are used repeatedly
//' (e.g., df=1,2,3 across thousands of genes).
//'
//' @param D Vector of degrees of freedom for a GFisher statistic.
//' @param M Correlation matrix of the input Z-scores.
//' @param p_type String: "two" for two-sided, "one" for one-sided input p-values.
//'
//' @return A list of 4 vectors of integrals:
//' \itemize{
//'   \item For two-sided: list(coeff2, coeff4, coeff6, coeff8)
//'   \item For one-sided: list(coeff1, coeff2, coeff3, coeff4)
//' }
//'
//' @details
//' This function provides cached computation of coefficients using the same integration
//' routines as getGFishercoef_cpp but with an efficient caching layer at the individual
//' coefficient level. When the same df value is requested multiple times, the cached result is
//' returned immediately without recomputing expensive integrals. The cache has
//' a size limit to prevent unbounded memory growth.
//'
//' Cache efficiency: For genome-wide analyses where df values are reused across
//' many tests, this caching can provide significant speedups by avoiding redundant
//' numerical integration.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List getGFishercoef_cached_cpp(arma::vec D, arma::mat M, std::string p_type) {
  int n = M.n_rows;

  // Handle single df case: replicate to match dimension of M
  arma::vec D_vec = D;
  if (D.n_elem == 1) {
    D_vec = arma::vec(n);
    D_vec.fill(D(0));
  }

  // Validate inputs
  if (D_vec.n_elem != (unsigned)n) {
    Rcpp::stop("Length of D must be 1 or match dimension of M");
  }
  if (M.n_cols != (unsigned)n) {
    Rcpp::stop("M must be a square matrix");
  }

  if (p_type == "two") {
    // Two-sided: compute coefficients for Hermite polynomials H_2, H_4, H_6, H_8
    arma::vec coeff2(n), coeff4(n), coeff6(n), coeff8(n);

    // Process each df value with caching
    for (int i = 0; i < n; i++) {
      coeff2(i) = get_cached_coefficient(D_vec(i), p_type, 2);
      coeff4(i) = get_cached_coefficient(D_vec(i), p_type, 4);
      coeff6(i) = get_cached_coefficient(D_vec(i), p_type, 6);
      coeff8(i) = get_cached_coefficient(D_vec(i), p_type, 8);
    }

    // Convert arma::vec to NumericVector to ensure plain vectors in R
    return Rcpp::List::create(
      Rcpp::Named("coeff2") = Rcpp::NumericVector(coeff2.begin(), coeff2.end()),
      Rcpp::Named("coeff4") = Rcpp::NumericVector(coeff4.begin(), coeff4.end()),
      Rcpp::Named("coeff6") = Rcpp::NumericVector(coeff6.begin(), coeff6.end()),
      Rcpp::Named("coeff8") = Rcpp::NumericVector(coeff8.begin(), coeff8.end())
    );

  } else if (p_type == "one") {
    // One-sided: compute coefficients for Hermite polynomials H_1, H_2, H_3, H_4
    arma::vec coeff1(n), coeff2(n), coeff3(n), coeff4(n);

    // Process each df value with caching
    for (int i = 0; i < n; i++) {
      coeff1(i) = get_cached_coefficient(D_vec(i), p_type, 1);
      coeff2(i) = get_cached_coefficient(D_vec(i), p_type, 2);
      coeff3(i) = get_cached_coefficient(D_vec(i), p_type, 3);
      coeff4(i) = get_cached_coefficient(D_vec(i), p_type, 4);
    }

    // Convert arma::vec to NumericVector to ensure plain vectors in R
    return Rcpp::List::create(
      Rcpp::Named("coeff1") = Rcpp::NumericVector(coeff1.begin(), coeff1.end()),
      Rcpp::Named("coeff2") = Rcpp::NumericVector(coeff2.begin(), coeff2.end()),
      Rcpp::Named("coeff3") = Rcpp::NumericVector(coeff3.begin(), coeff3.end()),
      Rcpp::Named("coeff4") = Rcpp::NumericVector(coeff4.begin(), coeff4.end())
    );

  } else {
    Rcpp::stop("p_type must be 'two' or 'one'");
  }
}

//' Clear GFisher Coefficient Cache
//'
//' Clears the internal cache used by getGFishercoef_cached_cpp.
//' This can be useful to free memory or force recomputation.
//'
//' @return The number of entries that were cleared from the cache.
//'
//' @export
// [[Rcpp::export]]
int clear_GFisher_cache() {
  std::lock_guard<std::mutex> lock(cache_mutex);
  int size = coefficient_cache.size();
  coefficient_cache.clear();
  return size;
}

//' Get GFisher Coefficient Cache Size
//'
//' Returns the current number of entries in the coefficient cache.
//'
//' @return The number of cached coefficient values.
//'
//' @export
// [[Rcpp::export]]
int get_GFisher_cache_size() {
  std::lock_guard<std::mutex> lock(cache_mutex);
  return coefficient_cache.size();
}

//' Get GFisher Cache Statistics
//'
//' Returns detailed statistics about the coefficient cache.
//'
//' @return A list containing cache statistics:
//' \itemize{
//'   \item size: Current number of entries
//'   \item max_size: Maximum cache size before automatic clearing
//'   \item memory_kb: Approximate memory usage in KB
//'   \item fill_ratio: Percentage of max capacity used
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List get_GFisher_cache_stats() {
  std::lock_guard<std::mutex> lock(cache_mutex);

  int size = coefficient_cache.size();
  double memory_kb = size * (sizeof(std::string) * 20 + sizeof(double)) / 1024.0;  // Approximate
  double fill_ratio = (double)size / MAX_CACHE_SIZE * 100.0;

  return Rcpp::List::create(
    Rcpp::Named("size") = size,
    Rcpp::Named("max_size") = (int)MAX_CACHE_SIZE,
    Rcpp::Named("memory_kb") = memory_kb,
    Rcpp::Named("fill_ratio") = fill_ratio
  );
}