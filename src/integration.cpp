// File: integration.cpp
// Purpose: Fast C++ implementation of numerical integration for GFisher
// Date: 2025-10-13
// Author: r-rcpp-developer agent

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Forward declarations
// Cached coefficient function from coefficient_cache.cpp
Rcpp::List getGFishercoef_cached_cpp(arma::vec D, arma::mat M, std::string p_type);
// Non-cached coefficient function (defined later in this file)
Rcpp::List getGFishercoef_cpp(arma::vec D, arma::mat M, std::string p_type);

// Composite Gauss-Legendre quadrature for 1D numerical integration
// Uses multiple subintervals for better accuracy
// Using 15-point Gauss-Legendre on 8 subintervals = 120 evaluation points total
class CompositeGaussLegendre {
private:
  // 15-point Gauss-Legendre quadrature nodes and weights
  static constexpr int N = 15;
  static constexpr int NUM_INTERVALS = 8;  // Divide interval into 8 parts for better accuracy

  // Pre-computed nodes (zeros of Legendre polynomial P_15)
  static constexpr double nodes[15] = {
    -0.9879925180204854, -0.9372733924007060, -0.8482065834104272,
    -0.7244177313601701, -0.5709721726085388, -0.3941513470775634,
    -0.2011940939974345, -0.0000000000000000,  0.2011940939974345,
     0.3941513470775634,  0.5709721726085388,  0.7244177313601701,
     0.8482065834104272,  0.9372733924007060,  0.9879925180204854
  };

  // Pre-computed weights
  static constexpr double weights[15] = {
    0.0307532419961173, 0.0703660474881081, 0.1071592204671720,
    0.1395706779261544, 0.1662692058169939, 0.1861610000155622,
    0.1984314853271116, 0.2025782419255613, 0.1984314853271116,
    0.1861610000155622, 0.1662692058169939, 0.1395706779261544,
    0.1071592204671720, 0.0703660474881081, 0.0307532419961173
  };

public:
  // Integrate function f from a to b using composite Gauss-Legendre quadrature
  double integrate(std::function<double(double)> f, double a, double b) {
    double total_result = 0.0;
    double interval_length = (b - a) / NUM_INTERVALS;

    // Integrate over each subinterval
    for (int interval = 0; interval < NUM_INTERVALS; interval++) {
      double sub_a = a + interval * interval_length;
      double sub_b = sub_a + interval_length;

      // Transform from [-1, 1] to [sub_a, sub_b]
      double mid = (sub_b + sub_a) / 2.0;
      double half_length = (sub_b - sub_a) / 2.0;

      double interval_result = 0.0;
      for (int i = 0; i < N; i++) {
        double x = mid + half_length * nodes[i];
        interval_result += weights[i] * f(x);
      }

      total_result += interval_result * half_length;
    }

    return total_result;
  }
};


// Helper function: Hermite polynomial H_2(x) = x^2 - 1
inline double hermite_2(double x) {
  return x * x - 1.0;
}

// Helper function: Hermite polynomial H_4(x) = x^4 - 6*x^2 + 3
inline double hermite_4(double x) {
  double x2 = x * x;
  return x2 * x2 - 6.0 * x2 + 3.0;
}

// Helper function: Hermite polynomial H_6(x) = x^6 - 15*x^4 + 45*x^2 - 15
inline double hermite_6(double x) {
  double x2 = x * x;
  double x4 = x2 * x2;
  return x2 * x4 - 15.0 * x4 + 45.0 * x2 - 15.0;
}

// Helper function: Hermite polynomial H_8(x) = x^8 - 28*x^6 + 210*x^4 - 420*x^2 + 105
inline double hermite_8(double x) {
  double x2 = x * x;
  double x4 = x2 * x2;
  double x6 = x2 * x4;
  return x4 * x4 - 28.0 * x6 + 210.0 * x4 - 420.0 * x2 + 105.0;
}

// Helper function: Hermite polynomial H_1(x) = x
inline double hermite_1(double x) {
  return x;
}

// Helper function: Hermite polynomial H_3(x) = x^3 - 3*x
inline double hermite_3(double x) {
  return x * x * x - 3.0 * x;
}


//' Compute Coefficient Integrals for Two-Sided P-Values (C++)
//'
//' Fast C++ implementation of coefficient calculations using Gauss-Legendre quadrature.
//' For two-sided p-values, computes integrals of the form:
//' \deqn{I_j(k) = \int_{-8}^{8} F^{-1}_{d_j}(\Phi(x^2)) \phi(x) H_k(x) dx}
//'
//' @param D Vector of degrees of freedom
//' @param coeff_type Integer: 2, 4, 6, or 8 for the Hermite polynomial order
//' @return Vector of coefficient integrals
//' @keywords internal
// [[Rcpp::export]]
arma::vec compute_coefficients_two_cpp(arma::vec D, int coeff_type) {
  int n = D.n_elem;
  arma::vec coeffs(n);

  CompositeGaussLegendre integrator;

  for (int i = 0; i < n; i++) {
    double df = D(i);

    // Define integrand based on coefficient type
    std::function<double(double)> integrand;

    if (coeff_type == 2) {
      integrand = [df](double x) {
        double chi_val = R::pchisq(x * x, 1.0, 1, 0);  // pchisq(x^2, df=1)
        double chi_inv = R::qchisq(chi_val, df, 1, 0);  // qchisq(..., df=df)
        double norm_dens = R::dnorm(x, 0.0, 1.0, 0);
        return chi_inv * norm_dens * hermite_2(x);
      };
    } else if (coeff_type == 4) {
      integrand = [df](double x) {
        double chi_val = R::pchisq(x * x, 1.0, 1, 0);
        double chi_inv = R::qchisq(chi_val, df, 1, 0);
        double norm_dens = R::dnorm(x, 0.0, 1.0, 0);
        return chi_inv * norm_dens * hermite_4(x);
      };
    } else if (coeff_type == 6) {
      integrand = [df](double x) {
        double chi_val = R::pchisq(x * x, 1.0, 1, 0);
        double chi_inv = R::qchisq(chi_val, df, 1, 0);
        double norm_dens = R::dnorm(x, 0.0, 1.0, 0);
        return chi_inv * norm_dens * hermite_6(x);
      };
    } else if (coeff_type == 8) {
      integrand = [df](double x) {
        double chi_val = R::pchisq(x * x, 1.0, 1, 0);
        double chi_inv = R::qchisq(chi_val, df, 1, 0);
        double norm_dens = R::dnorm(x, 0.0, 1.0, 0);
        return chi_inv * norm_dens * hermite_8(x);
      };
    } else {
      Rcpp::stop("Invalid coeff_type. Must be 2, 4, 6, or 8.");
    }

    // Integrate from -8 to 8
    coeffs(i) = integrator.integrate(integrand, -8.0, 8.0);
  }

  return coeffs;
}


//' Compute Coefficient Integrals for One-Sided P-Values (C++)
//'
//' Fast C++ implementation for one-sided p-values.
//' Computes integrals of the form:
//' \deqn{I_j(k) = \int_{-8}^{8} F^{-1}_{d_j}(\Phi(x)) \phi(x) H_k(x) dx}
//'
//' @param D Vector of degrees of freedom
//' @param coeff_type Integer: 1, 2, 3, or 4 for the Hermite polynomial order
//' @return Vector of coefficient integrals
//' @keywords internal
// [[Rcpp::export]]
arma::vec compute_coefficients_one_cpp(arma::vec D, int coeff_type) {
  int n = D.n_elem;
  arma::vec coeffs(n);

  CompositeGaussLegendre integrator;

  for (int i = 0; i < n; i++) {
    double df = D(i);

    // Define integrand based on coefficient type
    std::function<double(double)> integrand;

    if (coeff_type == 1) {
      integrand = [df](double x) {
        double norm_cdf = R::pnorm(x, 0.0, 1.0, 1, 0);  // pnorm(x)
        double chi_inv = R::qchisq(norm_cdf, df, 1, 0);  // qchisq(..., df=df)
        double norm_dens = R::dnorm(x, 0.0, 1.0, 0);
        return chi_inv * norm_dens * hermite_1(x);
      };
    } else if (coeff_type == 2) {
      integrand = [df](double x) {
        double norm_cdf = R::pnorm(x, 0.0, 1.0, 1, 0);
        double chi_inv = R::qchisq(norm_cdf, df, 1, 0);
        double norm_dens = R::dnorm(x, 0.0, 1.0, 0);
        return chi_inv * norm_dens * hermite_2(x);
      };
    } else if (coeff_type == 3) {
      integrand = [df](double x) {
        double norm_cdf = R::pnorm(x, 0.0, 1.0, 1, 0);
        double chi_inv = R::qchisq(norm_cdf, df, 1, 0);
        double norm_dens = R::dnorm(x, 0.0, 1.0, 0);
        return chi_inv * norm_dens * hermite_3(x);
      };
    } else if (coeff_type == 4) {
      integrand = [df](double x) {
        double norm_cdf = R::pnorm(x, 0.0, 1.0, 1, 0);
        double chi_inv = R::qchisq(norm_cdf, df, 1, 0);
        double norm_dens = R::dnorm(x, 0.0, 1.0, 0);
        return chi_inv * norm_dens * hermite_4(x);
      };
    } else {
      Rcpp::stop("Invalid coeff_type. Must be 1, 2, 3, or 4.");
    }

    // Integrate from -8 to 8
    coeffs(i) = integrator.integrate(integrand, -8.0, 8.0);
  }

  return coeffs;
}


//' Cached C++ Implementation of getGFisherGM
//'
//' Computes covariance matrix for weighted GFisher statistics using optimized
//' numerical integration with caching for repeated df values.
//'
//' @param D Vector of degrees of freedom
//' @param w Vector of weights
//' @param M Correlation matrix
//' @param p_type String: "two" for two-sided, "one" for one-sided
//' @return Covariance matrix GM
//'
//' @details
//' This function implements the covariance matrix calculation with automatic
//' caching. It retrieves coefficients using the cached version
//' \code{getGFishercoef_cached_cpp}, which checks the cache first before
//' computing. For analyses with repeated df values (e.g., analyzing thousands
//' of gene sets with df=2), this provides significant speedup through caching.
//'
//' For two-sided p-values, uses cached coefficients coeff2, coeff4, coeff6, coeff8.
//' For one-sided p-values, uses cached coeff1, coeff2, coeff3, coeff4.
//'
//' The covariance matrix is built using Hermite polynomial expansions as described
//' in the GFisher paper.
//'
//' @references
//' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and accurate p-value
//' calculation under dependence. Biometrics, 79(2), 1159-1172.
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
arma::mat getGFisherGM_cached_cpp(arma::vec D, arma::vec w, arma::mat M, std::string p_type) {
  int n = M.n_rows;

  // Validate inputs
  if (D.n_elem != (unsigned)n) {
    Rcpp::stop("Length of D must match dimension of M");
  }
  if (w.n_elem != (unsigned)n) {
    Rcpp::stop("Length of w must match dimension of M");
  }
  if (M.n_cols != (unsigned)n) {
    Rcpp::stop("M must be a square matrix");
  }

  arma::mat GM;

  if (p_type == "two") {
    // Two-sided: use cached coefficients from getGFishercoef_cached_cpp
    // This eliminates redundant integration calculations and leverages caching
    Rcpp::List coeff_list = getGFishercoef_cached_cpp(D, M, p_type);
    arma::vec coeff2 = Rcpp::as<arma::vec>(coeff_list["coeff2"]);
    arma::vec coeff4 = Rcpp::as<arma::vec>(coeff_list["coeff4"]);
    arma::vec coeff6 = Rcpp::as<arma::vec>(coeff_list["coeff6"]);
    arma::vec coeff8 = Rcpp::as<arma::vec>(coeff_list["coeff8"]);

    // Compute covariance matrix using Hermite expansion
    // GM = M^2/2 * coeff2 %*% t(coeff2) + M^4/24 * coeff4 %*% t(coeff4) + ...
    arma::mat M2 = arma::square(M);  // Element-wise square
    arma::mat M4 = arma::square(M2);
    arma::mat M6 = M2 % M4;  // Element-wise multiplication
    arma::mat M8 = arma::square(M4);

    GM = M2 / 2.0 % (coeff2 * coeff2.t()) +
         M4 / 24.0 % (coeff4 * coeff4.t()) +
         M6 / 720.0 % (coeff6 * coeff6.t()) +
         M8 / 40320.0 % (coeff8 * coeff8.t());

    // Apply weights and convert to correlation, then back to covariance
    // GM = diag(sqrt(2*D)*w) %*% cov2cor(GM) %*% diag(sqrt(2*D)*w)

    // First, compute correlation matrix from GM
    arma::vec diag_GM = arma::diagvec(GM);
    arma::vec sqrt_diag_GM = arma::sqrt(diag_GM);

    // Avoid division by zero
    for (int i = 0; i < n; i++) {
      if (sqrt_diag_GM(i) < 1e-10) {
        sqrt_diag_GM(i) = 1.0;
      }
    }

    // Compute correlation matrix
    arma::mat inv_diag = arma::diagmat(1.0 / sqrt_diag_GM);
    arma::mat GM_cor = inv_diag * GM * inv_diag;

    // Ensure diagonal is exactly 1 (numerical stability)
    GM_cor.diag().ones();

    // Apply final scaling with weights
    arma::vec scaling = arma::sqrt(2.0 * D) % w;
    arma::mat diag_scale = arma::diagmat(scaling);
    GM = diag_scale * GM_cor * diag_scale;

  } else if (p_type == "one") {
    // One-sided: use cached coefficients from getGFishercoef_cached_cpp
    // This eliminates redundant integration calculations and leverages caching
    Rcpp::List coeff_list = getGFishercoef_cached_cpp(D, M, p_type);
    arma::vec coeff1 = Rcpp::as<arma::vec>(coeff_list["coeff1"]);
    arma::vec coeff2 = Rcpp::as<arma::vec>(coeff_list["coeff2"]);
    arma::vec coeff3 = Rcpp::as<arma::vec>(coeff_list["coeff3"]);
    arma::vec coeff4 = Rcpp::as<arma::vec>(coeff_list["coeff4"]);

    // Compute covariance matrix
    arma::mat M2 = arma::square(M);
    arma::mat M3 = M % M2;
    arma::mat M4 = arma::square(M2);

    GM = M % (coeff1 * coeff1.t()) +
         M2 / 2.0 % (coeff2 * coeff2.t()) +
         M3 / 6.0 % (coeff3 * coeff3.t()) +
         M4 / 24.0 % (coeff4 * coeff4.t());

    // Apply weights and convert to correlation, then back to covariance
    arma::vec diag_GM = arma::diagvec(GM);
    arma::vec sqrt_diag_GM = arma::sqrt(diag_GM);

    for (int i = 0; i < n; i++) {
      if (sqrt_diag_GM(i) < 1e-10) {
        sqrt_diag_GM(i) = 1.0;
      }
    }

    arma::mat inv_diag = arma::diagmat(1.0 / sqrt_diag_GM);
    arma::mat GM_cor = inv_diag * GM * inv_diag;
    GM_cor.diag().ones();

    arma::vec scaling = arma::sqrt(2.0 * D) % w;
    arma::mat diag_scale = arma::diagmat(scaling);
    GM = diag_scale * GM_cor * diag_scale;

  } else {
    Rcpp::stop("p_type must be 'two' or 'one'");
  }

  return GM;
}


//' Non-Cached C++ Implementation of getGFisherGM
//'
//' Computes covariance matrix for weighted GFisher statistics using optimized
//' numerical integration WITHOUT caching (robust fallback version).
//'
//' @param D Vector of degrees of freedom
//' @param w Vector of weights
//' @param M Correlation matrix
//' @param p_type String: "two" for two-sided, "one" for one-sided
//' @return Covariance matrix GM
//'
//' @details
//' This function implements the covariance matrix calculation without caching.
//' It calls \code{getGFishercoef_cpp()} which always computes coefficients fresh
//' without checking the cache. This version is used as a fallback if the cached
//' version encounters issues (extremely rare).
//'
//' This function shares the same coefficient computation logic with the cached
//' version, ensuring consistency. The only difference is the absence of caching,
//' making it simpler and more robust for fallback scenarios.
//'
//' For two-sided p-values, computes coefficients coeff2, coeff4, coeff6, coeff8.
//' For one-sided p-values, computes coeff1, coeff2, coeff3, coeff4.
//'
//' @references
//' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and accurate p-value
//' calculation under dependence. Biometrics, 79(2), 1159-1172.
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
arma::mat getGFisherGM_cpp(arma::vec D, arma::vec w, arma::mat M, std::string p_type) {
  int n = M.n_rows;

  // Validate inputs
  if (D.n_elem != (unsigned)n) {
    Rcpp::stop("Length of D must match dimension of M");
  }
  if (w.n_elem != (unsigned)n) {
    Rcpp::stop("Length of w must match dimension of M");
  }
  if (M.n_cols != (unsigned)n) {
    Rcpp::stop("M must be a square matrix");
  }

  arma::mat GM;

  if (p_type == "two") {
    // Two-sided: use non-cached coefficients from getGFishercoef_cpp
    // Always computes fresh without caching (robust fallback)
    Rcpp::List coeff_list = getGFishercoef_cpp(D, M, p_type);
    arma::vec coeff2 = Rcpp::as<arma::vec>(coeff_list["coeff2"]);
    arma::vec coeff4 = Rcpp::as<arma::vec>(coeff_list["coeff4"]);
    arma::vec coeff6 = Rcpp::as<arma::vec>(coeff_list["coeff6"]);
    arma::vec coeff8 = Rcpp::as<arma::vec>(coeff_list["coeff8"]);

    // Compute covariance matrix using Hermite expansion
    // GM = M^2/2 * coeff2 %*% t(coeff2) + M^4/24 * coeff4 %*% t(coeff4) + ...
    arma::mat M2 = arma::square(M);  // Element-wise square
    arma::mat M4 = arma::square(M2);
    arma::mat M6 = M2 % M4;  // Element-wise multiplication
    arma::mat M8 = arma::square(M4);

    GM = M2 / 2.0 % (coeff2 * coeff2.t()) +
         M4 / 24.0 % (coeff4 * coeff4.t()) +
         M6 / 720.0 % (coeff6 * coeff6.t()) +
         M8 / 40320.0 % (coeff8 * coeff8.t());

    // Apply weights and convert to correlation, then back to covariance
    // GM = diag(sqrt(2*D)*w) %*% cov2cor(GM) %*% diag(sqrt(2*D)*w)

    // First, compute correlation matrix from GM
    arma::vec diag_GM = arma::diagvec(GM);
    arma::vec sqrt_diag_GM = arma::sqrt(diag_GM);

    // Avoid division by zero
    for (int i = 0; i < n; i++) {
      if (sqrt_diag_GM(i) < 1e-10) {
        sqrt_diag_GM(i) = 1.0;
      }
    }

    // Compute correlation matrix
    arma::mat inv_diag = arma::diagmat(1.0 / sqrt_diag_GM);
    arma::mat GM_cor = inv_diag * GM * inv_diag;

    // Ensure diagonal is exactly 1 (numerical stability)
    GM_cor.diag().ones();

    // Apply final scaling with weights
    arma::vec scaling = arma::sqrt(2.0 * D) % w;
    arma::mat diag_scale = arma::diagmat(scaling);
    GM = diag_scale * GM_cor * diag_scale;

  } else if (p_type == "one") {
    // One-sided: use non-cached coefficients from getGFishercoef_cpp
    // Always computes fresh without caching (robust fallback)
    Rcpp::List coeff_list = getGFishercoef_cpp(D, M, p_type);
    arma::vec coeff1 = Rcpp::as<arma::vec>(coeff_list["coeff1"]);
    arma::vec coeff2 = Rcpp::as<arma::vec>(coeff_list["coeff2"]);
    arma::vec coeff3 = Rcpp::as<arma::vec>(coeff_list["coeff3"]);
    arma::vec coeff4 = Rcpp::as<arma::vec>(coeff_list["coeff4"]);

    // Compute covariance matrix
    arma::mat M2 = arma::square(M);
    arma::mat M3 = M % M2;
    arma::mat M4 = arma::square(M2);

    GM = M % (coeff1 * coeff1.t()) +
         M2 / 2.0 % (coeff2 * coeff2.t()) +
         M3 / 6.0 % (coeff3 * coeff3.t()) +
         M4 / 24.0 % (coeff4 * coeff4.t());

    // Apply weights and convert to correlation, then back to covariance
    arma::vec diag_GM = arma::diagvec(GM);
    arma::vec sqrt_diag_GM = arma::sqrt(diag_GM);

    for (int i = 0; i < n; i++) {
      if (sqrt_diag_GM(i) < 1e-10) {
        sqrt_diag_GM(i) = 1.0;
      }
    }

    arma::mat inv_diag = arma::diagmat(1.0 / sqrt_diag_GM);
    arma::mat GM_cor = inv_diag * GM * inv_diag;
    GM_cor.diag().ones();

    arma::vec scaling = arma::sqrt(2.0 * D) % w;
    arma::mat diag_scale = arma::diagmat(scaling);
    GM = diag_scale * GM_cor * diag_scale;

  } else {
    Rcpp::stop("p_type must be 'two' or 'one'");
  }

  return GM;
}


//' Compute Hermite Polynomial Coefficients for GFisher (C++)
//'
//' Fast C++ implementation of coefficient calculations for GFisher.
//' This function computes the raw coefficient integrals needed for
//' computing covariances between GFisher statistics.
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
//' This function implements the literal calculation of integrals in formula (7)
//' of the GFisher paper. Uses composite Gauss-Legendre quadrature for numerical
//' integration, providing substantial speedup over R's integrate() function.
//'
//' For two-sided p-values, computes:
//' \deqn{I_j(k) = \int_{-8}^{8} F^{-1}_{d_j}(\Phi(x^2)) \phi(x) H_k(x) dx}
//'
//' For one-sided p-values, computes:
//' \deqn{I_j(k) = \int_{-8}^{8} F^{-1}_{d_j}(\Phi(x)) \phi(x) H_k(x) dx}
//'
//' If a single degree of freedom is provided, it is replicated to match
//' the dimension of M.
//'
//' @references
//' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and accurate p-value
//' calculation under dependence. Biometrics, 79(2), 1159-1172.
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List getGFishercoef_cpp(arma::vec D, arma::mat M, std::string p_type) {
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
    arma::vec coeff2 = compute_coefficients_two_cpp(D_vec, 2);
    arma::vec coeff4 = compute_coefficients_two_cpp(D_vec, 4);
    arma::vec coeff6 = compute_coefficients_two_cpp(D_vec, 6);
    arma::vec coeff8 = compute_coefficients_two_cpp(D_vec, 8);

    // Convert arma::vec to NumericVector to ensure plain vectors in R
    return Rcpp::List::create(
      Rcpp::Named("coeff2") = Rcpp::NumericVector(coeff2.begin(), coeff2.end()),
      Rcpp::Named("coeff4") = Rcpp::NumericVector(coeff4.begin(), coeff4.end()),
      Rcpp::Named("coeff6") = Rcpp::NumericVector(coeff6.begin(), coeff6.end()),
      Rcpp::Named("coeff8") = Rcpp::NumericVector(coeff8.begin(), coeff8.end())
    );

  } else if (p_type == "one") {
    // One-sided: compute coefficients for Hermite polynomials H_1, H_2, H_3, H_4
    arma::vec coeff1 = compute_coefficients_one_cpp(D_vec, 1);
    arma::vec coeff2 = compute_coefficients_one_cpp(D_vec, 2);
    arma::vec coeff3 = compute_coefficients_one_cpp(D_vec, 3);
    arma::vec coeff4 = compute_coefficients_one_cpp(D_vec, 4);

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
