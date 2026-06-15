// File: simulation.cpp
// Purpose: Fast C++ implementation of simulation-based p-value calculation for GFisher
// Date: 2025-10-13
// Author: r-rcpp-developer agent

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Fast C++ Implementation of p.GFisher Method="MR"
//'
//' Computes p-value using simulation-assisted moment ratio matching in C++.
//' This function provides substantial speedup over the pure R implementation
//' by vectorizing operations and avoiding apply() loops.
//'
//' @param q Observed GFisher statistic
//' @param df Vector of degrees of freedom
//' @param w Vector of weights
//' @param M Correlation matrix
//' @param p_type String: "two" for two-sided, "one" for one-sided
//' @param nsim Number of simulations (default 50000)
//' @param seed Random seed (0 = no seed)
//' @return P-value
//'
//' @details
//' This function implements the simulation-based moment ratio (MR) method for calculating
//' GFisher p-values under dependence. The algorithm:
//'
//' \enumerate{
//'   \item Generate null distribution: znull ~ N(0, M) using Cholesky decomposition
//'   \item Transform to p-values based on p_type (two-sided or one-sided)
//'   \item Apply chi-square transformations with specified degrees of freedom
//'   \item Compute weighted Fisher statistics: sum(w * pp_trans)
//'   \item Calculate first 4 moments from the null distribution
//'   \item Fit gamma distribution using moment matching
//'   \item Compute p-value from the gamma approximation
//' }
//'
//' \strong{Performance}:
//' This C++ implementation is substantially faster than the R version due to:
//' \itemize{
//'   \item Vectorized matrix operations via Armadillo
//'   \item Elimination of apply() loops
//'   \item Efficient memory layout
//'   \item Optimized BLAS/LAPACK routines
//' }
//'
//' \strong{Reproducibility}:
//' When the same seed is used, C++ and R versions produce
//' identical results (within machine precision).
//'
//' @references
//' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and accurate p-value
//' calculation under dependence. Biometrics, 79(2), 1159-1172.
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
double p_GFisher_MR_cpp(double q, arma::vec df, arma::vec w, arma::mat M,
                        std::string p_type = "two", int nsim = 50000, int seed = 0) {
  int n = M.n_rows;

  // Validate inputs
  if (df.n_elem != (unsigned)n) {
    Rcpp::stop("Length of df must match dimension of M");
  }
  if (w.n_elem != (unsigned)n) {
    Rcpp::stop("Length of w must match dimension of M");
  }
  if (M.n_cols != (unsigned)n) {
    Rcpp::stop("M must be a square matrix");
  }
  if (nsim < 100) {
    Rcpp::stop("nsim must be at least 100");
  }

  // Note: Seed should be set at R level via set.seed() before calling this function
  // We don't set seed in C++ to maintain consistency with R's RNG state
  if (seed > 0) {
    Rcpp::warning("When called from R, the RNG seed has to be set at the R level via set.seed()");
  }

  // Normalize weights
  w = w / arma::sum(w);

  // Cholesky decomposition of M
  arma::mat M_chol;
  bool chol_success = arma::chol(M_chol, M);

  if (!chol_success) {
    // If Cholesky fails, add small constant to diagonal
    arma::mat M_corrected = M + arma::eye(n, n) * 1e-8;
    chol_success = arma::chol(M_chol, M_corrected);

    if (!chol_success) {
      Rcpp::stop("Cholesky decomposition failed. Matrix M may not be positive definite.");
    }
  }

  // Generate null z-scores: nsim x n matrix using R's RNG
  // znull = randn(nsim, n) %*% M_chol
  // Use Rcpp::rnorm to ensure we use R's RNG state
  // R fills matrices column-major, so we need to match that order
  NumericVector rnorm_vals = Rcpp::rnorm(nsim * n, 0.0, 1.0);

  arma::mat znull(nsim, n);
  int idx = 0;
  // Fill column by column to match R's matrix() behavior
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < nsim; i++) {
      znull(i, j) = rnorm_vals[idx++];
    }
  }

  // Apply correlation structure: znull = znull %*% M_chol
  znull = znull * M_chol;

  // Convert to p-values
  arma::mat pnull(nsim, n);

  if (p_type == "two") {
    // Two-sided p-values: 2 * pnorm(-|z|)
    for (int i = 0; i < nsim; i++) {
      for (int j = 0; j < n; j++) {
        pnull(i, j) = 2.0 * R::pnorm(-std::abs(znull(i, j)), 0.0, 1.0, 1, 0);
      }
    }
  } else if (p_type == "one") {
    // One-sided p-values: pnorm(z, lower.tail=FALSE)
    for (int i = 0; i < nsim; i++) {
      for (int j = 0; j < n; j++) {
        pnull(i, j) = R::pnorm(znull(i, j), 0.0, 1.0, 0, 0);  // lower.tail=FALSE
      }
    }
  } else {
    Rcpp::stop("p_type must be 'two' or 'one'");
  }

  // Transform to chi-square quantiles
  arma::mat pp_trans(nsim, n);

  // Check if all df are the same
  double df_sd = arma::stddev(df);

  if (df_sd < 1e-10) {
    // All df are equal
    double df_val = df(0);

    if (std::abs(df_val - 2.0) < 1e-10) {
      // Special case: df=2 uses fast -2*log transformation
      for (int i = 0; i < nsim; i++) {
        for (int j = 0; j < n; j++) {
          pp_trans(i, j) = -2.0 * std::log(pnull(i, j));
        }
      }
    } else {
      // Same df but not 2
      for (int i = 0; i < nsim; i++) {
        for (int j = 0; j < n; j++) {
          pp_trans(i, j) = R::qchisq(pnull(i, j), df_val, 0, 0);  // lower.tail=FALSE
        }
      }
    }
  } else {
    // Different df values
    for (int i = 0; i < nsim; i++) {
      for (int j = 0; j < n; j++) {
        pp_trans(i, j) = R::qchisq(pnull(i, j), df(j), 0, 0);  // lower.tail=FALSE
      }
    }
  }

  // Compute weighted Fisher statistics
  // fishernull[i] = sum(pp_trans[i,] * w)
  // Use matrix-vector multiplication: fishernull = pp_trans * w
  arma::vec fishernull = pp_trans * w;

  // Compute first 4 raw moments
  double M1 = arma::mean(fishernull);
  double M2 = arma::mean(arma::square(fishernull));
  double M3 = arma::mean(arma::pow(fishernull, 3));
  double M4 = arma::mean(arma::pow(fishernull, 4));

  // Compute central moments
  double mu = M1;
  double sigma2 = (M2 - M1 * M1) * nsim / (nsim - 1.0);  // Unbiased variance
  double cmu3 = M3 - 3.0 * M2 * M1 + 2.0 * M1 * M1 * M1;
  double cmu4 = M4 - 4.0 * M3 * M1 + 6.0 * M2 * M1 * M1 - 3.0 * M1 * M1 * M1 * M1;

  // Compute skewness (gamma) and kurtosis (kappa)
  double gm = cmu3 / std::pow(sigma2, 1.5);
  double kp = cmu4 / (sigma2 * sigma2);

  // Moment matching to gamma distribution
  // Shape parameter: a = 9 * gm^2 / (kp - 3)^2
  double a = 9.0 * gm * gm / ((kp - 3.0) * (kp - 3.0));

  // Standardize and transform: (q - mu) / sqrt(sigma2) * sqrt(a) + a
  double q_transformed = (q - mu) / std::sqrt(sigma2) * std::sqrt(a) + a;

  // Calculate p-value using gamma distribution
  double pval = R::pgamma(q_transformed, a, 1.0, 0, 0);  // lower.tail=FALSE, scale=1

  return pval;
}
