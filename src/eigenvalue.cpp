// File: eigenvalue.cpp
// Purpose: Fast C++ implementation of eigenvalue calculations for GFisher
// Date: 2025-10-13
// Author: r-rcpp-developer agent

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Fast C++ Implementation of getGFisherlam
//'
//' Computes eigenvalues needed for quadratic approximation using RcppArmadillo's
//' optimized LAPACK routines for eigenvalue decomposition.
//'
//' @param D Vector of degrees of freedom
//' @param w Vector of weights
//' @param M Correlation matrix of input Z-scores
//' @param GM Covariance matrix from getGFisherGM
//' @return List with eigenvalues (lam)
//'
//' @details
//' This function implements Algorithm from Section 3.4 of the GFisher paper:
//'
//' 1. Construct M_tilde matrix using formula (13):
//'    M_tilde_{ij} = sign(M_{ij}) * sqrt(|GM_{ij}| / (min(d_i, d_j) * 2))
//'
//' 2. Ensure M_tilde is a valid correlation matrix (values <= 1, positive definite)
//'
//' 3. Compute weighted Cholesky: WM_chol = chol(M_tilde) * diag(sqrt(w))
//'
//' 4. Extract eigenvalues from WM_chol^T * WM_chol
//'
//' 5. For max(D) > 1, add additional eigenvalues for higher degrees of freedom
//'
//' 6. Filter out small eigenvalues (< 1e-10)
//'
//' The C++ implementation uses Armadillo's eig_sym() which directly calls
//' optimized LAPACK routines, providing significant speedup over R's eigen().
//'
//' @references
//' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and accurate p-value
//' calculation under dependence. Biometrics, 79(2), 1159-1172. See Section 3.4.
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
List getGFisherlam_cpp(arma::vec D, arma::vec w, arma::mat M, arma::mat GM) {
  int n = M.n_rows;

  // Validate inputs
  if (D.n_elem != (unsigned)n) {
    Rcpp::stop("Length of D must match dimension of M");
  }
  if (w.n_elem != (unsigned)n) {
    Rcpp::stop("Length of w must match dimension of M");
  }
  if (M.n_cols != (unsigned)n || GM.n_rows != (unsigned)n || GM.n_cols != (unsigned)n) {
    Rcpp::stop("M and GM must be n x n matrices");
  }

  // Step 1: Construct DM matrix - min(D_i, D_j) for all pairs
  arma::mat DM(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      DM(i, j) = std::min(D(i), D(j));
    }
  }

  // Step 2: Construct M_tilde using formula (13)
  // M_tilde = sign(M) * sqrt(|GM| / DM / 2)
  arma::mat M_tilde = arma::sign(M) % arma::sqrt(arma::abs(GM) / DM / 2.0);

  // Step 3: Cap values > 1 at 0.999 (ensure valid correlation)
  M_tilde.elem(arma::find(M_tilde > 1.0)).fill(0.999);
  M_tilde.elem(arma::find(M_tilde < -1.0)).fill(-0.999);

  // Step 4: Check if M_tilde is positive definite
  // If not, use iterative correction for numerical stability
  arma::vec eigval_check;
  bool eig_success = arma::eig_sym(eigval_check, M_tilde);

  if (!eig_success || arma::any(eigval_check < 1e-10)) {
    // Use iterative correction: add to diagonal until positive definite
    double min_eigval = arma::min(eigval_check);
    double correction = std::max(1e-8, -min_eigval + 1e-6);

    // Try increasing corrections until it works (up to 5 attempts)
    int attempts = 0;
    while (attempts < 5 && (!eig_success || arma::any(eigval_check < 1e-10))) {
      M_tilde.diag() += correction;
      eig_success = arma::eig_sym(eigval_check, M_tilde);
      correction *= 2.0;  // Double the correction each time
      attempts++;
    }

    if (attempts >= 5 && (!eig_success || arma::any(eigval_check < 1e-10))) {
      Rcpp::warning("M_tilde required significant correction for positive definiteness. Results may be less accurate.");
    }
  }

  // Step 5: Cholesky decomposition with fallback
  arma::mat M_tilde_chol;
  bool chol_success = arma::chol(M_tilde_chol, M_tilde);

  if (!chol_success) {
    // Last resort: add larger diagonal correction and try again
    M_tilde.diag() += 1e-4;
    chol_success = arma::chol(M_tilde_chol, M_tilde);

    if (!chol_success) {
      Rcpp::stop("Cholesky decomposition failed despite corrections. Matrix may be ill-conditioned.");
    } else {
      Rcpp::warning("Required large diagonal correction for Cholesky decomposition.");
    }
  }

  // Step 6: Compute WM_chol = M_tilde_chol * diag(sqrt(w))
  arma::vec sqrt_w = arma::sqrt(arma::abs(w));  // Take abs to handle any negative weights
  arma::mat WM_chol = M_tilde_chol * arma::diagmat(sqrt_w);

  // Step 7: First set of eigenvalues from WM_chol^T * WM_chol
  arma::mat temp = WM_chol.t() * WM_chol;
  arma::vec lam;
  arma::eig_sym(lam, temp);

  // Step 8: Additional eigenvalues for max(D) > 1
  double max_D = arma::max(D);
  if (max_D > 1) {
    for (int i = 2; i <= static_cast<int>(max_D); i++) {
      // Create Ai matrix (indicator for D >= i)
      arma::mat Ai = arma::eye<arma::mat>(n, n);

      // Find indices where D < i and set diagonal to 0
      arma::uvec Did = arma::find(D < i);
      if (Did.n_elem > 0) {
        for (unsigned int j = 0; j < Did.n_elem; j++) {
          Ai(Did(j), Did(j)) = 0.0;
        }
      }

      // Compute additional eigenvalues from WM_chol * Ai * WM_chol^T
      arma::vec lam_i;
      arma::mat temp_i = WM_chol * Ai * WM_chol.t();
      arma::eig_sym(lam_i, temp_i);

      // Append to lam vector
      lam = arma::join_cols(lam, lam_i);
    }
  }

  // Step 9: Filter small eigenvalues (keep only > 1e-10)
  arma::uvec large_idx = arma::find(lam > 1e-10);
  if (large_idx.n_elem > 0) {
    lam = lam.elem(large_idx);
  } else {
    // If all eigenvalues are small, keep at least one
    lam = arma::vec(1);
    lam(0) = 1e-10;
    Rcpp::warning("All eigenvalues were very small. Keeping one small eigenvalue.");
  }

  return List::create(Named("lam") = lam);
}
