# GFisher-helpers.R
# Helper functions for computing covariance, correlation, and eigenvalues
# needed for GFisher p-value calculations under dependence.

# Helper functions for computing algorithms


#' Compute Hermite Polynomial Coefficients for GFisher
#'
#' @title Compute Necessary Integrals for Theorem 1
#' @description Calculate the integrals (up to 4 non-zero terms) used in Theorem 1 (Formula 7)
#'   to approximate the covariance between transformed statistics.
#'
#' @param D A vector of degrees of freedom for a GFisher statistic.
#' @param M Correlation matrix of the input Z-scores from which the input p-values were obtained.
#' @param p.type Character string: \code{"two"} for two-sided (default), \code{"one"} for one-sided input p-values.
#' @param use.cpp Logical, default \code{TRUE}. If \code{TRUE}, uses fast C++ implementation.
#'   Set to \code{FALSE} to force pure R implementation.
#'
#' @return A list of 4 vectors of integrals used in formula (7):
#'   \item{coeff2}{Vector of \eqn{I_1(2), I_2(2), ..., I_n(2)} for two-sided, or \eqn{I_1(1), ..., I_n(1)} for one-sided}
#'   \item{coeff4}{Vector of \eqn{I_1(4), I_2(4), ..., I_n(4)} for two-sided, or \eqn{I_1(2), ..., I_n(2)} for one-sided}
#'   \item{coeff6}{Vector of \eqn{I_1(6), I_2(6), ..., I_n(6)} for two-sided, or \eqn{I_1(3), ..., I_n(3)} for one-sided}
#'   \item{coeff8}{Vector of \eqn{I_1(8), I_2(8), ..., I_n(8)} for two-sided, or \eqn{I_1(4), ..., I_n(4)} for one-sided}
#'
#' @details
#' This function implements the literal calculation of integrals in formula (7) of the GFisher paper.
#'
#' \strong{For two-sided p-values:}
#'
#' The function computes:
#' \deqn{I_j(k) = \int_{-8}^{8} F^{-1}_{d_j}(\Phi(x^2)) \phi(x) H_k(x) dx}
#' where \eqn{F^{-1}_{d_j}} is the inverse chi-square CDF, \eqn{\Phi} is the standard normal CDF,
#' \eqn{\phi} is the standard normal PDF, and \eqn{H_k(x)} are Hermite polynomials:
#' \itemize{
#'   \item \eqn{H_2(x) = x^2 - 1}
#'   \item \eqn{H_4(x) = x^4 - 6x^2 + 3}
#'   \item \eqn{H_6(x) = x^6 - 15x^4 + 45x^2 - 15}
#'   \item \eqn{H_8(x) = x^8 - 28x^6 + 210x^4 - 420x^2 + 105}
#' }
#'
#' \strong{For one-sided p-values:}
#'
#' Similar integrals are computed but with \eqn{\Phi(x)} instead of \eqn{\Phi(x^2)} and
#' different Hermite polynomials (\eqn{H_1, H_2, H_3, H_4}).
#'
#' If a single degree of freedom is provided, it is replicated to match the dimension of \code{M}.
#'
#' \strong{Performance Note:}
#' When \code{use.cpp = TRUE} (default), this function uses a C++ implementation with
#' composite Gauss-Legendre quadrature for numerical integration, providing substantial
#' speedup over the pure R implementation. The C++ version automatically falls back to
#' R if unavailable.
#'
#' @references
#' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and accurate p-value
#' calculation under dependence. \emph{Biometrics}, 79(2), 1159-1172. See Theorem 1 and Formula (7).
#'
#' @keywords internal
#' @author Hong Zhang
getGFishercoef <- function(D, M, p.type = "two", use.cpp = TRUE) {
  # Try C++ version first if requested
  if (use.cpp) {
    tryCatch({
      # Use cached version for better performance with repeated df values
      return(getGFishercoef_cached_cpp(D, M, p.type))
    }, error = function(e) {
      # If cached version fails, try regular C++ version
      tryCatch({
        return(getGFishercoef_cpp(D, M, p.type))
      }, error = function(e2) {
        warning("C++ version failed, falling back to R implementation: ", e2$message)
      })
    })
  }

  # Pure R implementation (fallback)
  # Get dimension from correlation matrix
  n <- dim(M)[1]

  if (p.type == "two") {
    # Two-sided p-values: integrate with Hermite polynomials H_2, H_4, H_6, H_8
    # Integrate from -8 to 8 for numerical stability
    coeff2 <- sapply(D, function(y) {
      stats::integrate(function(x) {
        stats::qchisq(stats::pchisq(x^2, df = 1), df = y) *
          stats::dnorm(x) * (x^2 - 1)
      }, -8, 8)$value
    })

    coeff4 <- sapply(D, function(y) {
      stats::integrate(function(x) {
        stats::qchisq(stats::pchisq(x^2, df = 1), df = y) *
          stats::dnorm(x) * (x^4 - 6 * x^2 + 3)
      }, -8, 8)$value
    })

    coeff6 <- sapply(D, function(y) {
      stats::integrate(function(x) {
        stats::qchisq(stats::pchisq(x^2, df = 1), df = y) *
          stats::dnorm(x) * (x^6 - 15 * x^4 + 45 * x^2 - 15)
      }, -8, 8)$value
    })

    coeff8 <- sapply(D, function(y) {
      stats::integrate(function(x) {
        stats::qchisq(stats::pchisq(x^2, df = 1), df = y) *
          stats::dnorm(x) * (x^8 - 28 * x^6 + 210 * x^4 - 420 * x^2 + 105)
      }, -8, 8)$value
    })

    # Replicate if single df provided
    if (length(D) == 1) {
      coeff2 <- rep(coeff2, n)
      coeff4 <- rep(coeff4, n)
      coeff6 <- rep(coeff6, n)
      coeff8 <- rep(coeff8, n)
    }

    return(list(coeff2 = coeff2, coeff4 = coeff4, coeff6 = coeff6, coeff8 = coeff8))

  } else {
    # One-sided p-values: integrate with Hermite polynomials H_1, H_2, H_3, H_4
    coeff1 <- sapply(D, function(y) {
      stats::integrate(function(x) {
        stats::qchisq(stats::pnorm(x), df = y) * stats::dnorm(x) * (x)
      }, -8, 8)$value
    })

    coeff2 <- sapply(D, function(y) {
      stats::integrate(function(x) {
        stats::qchisq(stats::pnorm(x), df = y) * stats::dnorm(x) * (x^2 - 1)
      }, -8, 8)$value
    })

    coeff3 <- sapply(D, function(y) {
      stats::integrate(function(x) {
        stats::qchisq(stats::pnorm(x), df = y) * stats::dnorm(x) * (x^3 - 3 * x)
      }, -8, 8)$value
    })

    coeff4 <- sapply(D, function(y) {
      stats::integrate(function(x) {
        stats::qchisq(stats::pnorm(x), df = y) * stats::dnorm(x) * (x^4 - 6 * x^2 + 3)
      }, -8, 8)$value
    })

    # Replicate if single df provided
    if (length(D) == 1) {
      coeff1 <- rep(coeff1, n)
      coeff2 <- rep(coeff2, n)
      coeff3 <- rep(coeff3, n)
      coeff4 <- rep(coeff4, n)
    }

    return(list(coeff1 = coeff1, coeff2 = coeff2, coeff3 = coeff3, coeff4 = coeff4))
  }
}


#' Compute Covariance Between Two GFisher Statistics
#'
#' @title Calculate Covariance Between GFisher Statistics (Corollary 2)
#' @description Calculate the covariance between two GFisher statistics with potentially
#'   different degrees of freedom and weights.
#'
#' @param D1 A vector of degrees of freedom for the first GFisher statistic.
#' @param D2 A vector of degrees of freedom for the second GFisher statistic.
#' @param W1 A vector of weights for the first GFisher statistic.
#' @param W2 A vector of weights for the second GFisher statistic.
#' @param M Correlation matrix of the input Z-scores from which the input p-values were obtained.
#' @param p.type Character string: \code{"two"} for two-sided (default), \code{"one"} for one-sided input p-values.
#' @param var.correct Logical, default \code{TRUE}. If \code{TRUE}, ensures exact variance is used
#'   by applying a variance correction factor.
#'
#' @return A numeric value representing the covariance between the two GFisher statistics \eqn{T^{(l)}} and \eqn{T^{(r)}}.
#'
#' @details
#' This function implements Corollary 2 from the GFisher paper, which provides a formula
#' for computing the covariance between two GFisher statistics:
#'
#' \deqn{\text{Cov}(T^{(l)}, T^{(r)}) = \sum_{i,j} w_i^{(l)} w_j^{(r)} \text{Cov}(T_i^{(l)}, T_j^{(r)})}
#'
#' The component covariances are approximated using Hermite polynomial expansions based on
#' the correlation structure \code{M}.
#'
#' \strong{For two-sided p-values:}
#'
#' The cross-covariance matrix is computed as:
#' \deqn{GM_{ij} = \sum_{k=1}^{4} \frac{\rho_{ij}^{2k}}{(2k)!} I_i^{(l)}(2k) I_j^{(r)}(2k)}
#'
#' \strong{For one-sided p-values:}
#'
#' A similar formula is used but with odd-order terms included.
#'
#' If \code{var.correct = TRUE}, the function rescales the covariance to ensure exact marginal variances.
#'
#' This includes Corollary 1 as a special case (when \eqn{l = r}).
#'
#' @references
#' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and accurate p-value
#' calculation under dependence. \emph{Biometrics}, 79(2), 1159-1172. See Corollary 2.
#'
#' @keywords internal
#' @author Hong Zhang
getGFishercov <- function(D1, D2, W1, W2, M, p.type = "two", var.correct = TRUE) {
  # Get dimension
  n <- dim(M)[1]

  # Compute coefficients for both statistics
  res1.coeff <- getGFishercoef(D1, M, p.type = p.type, use.cpp = TRUE)
  res2.coeff <- getGFishercoef(D2, M, p.type = p.type, use.cpp = TRUE)

  if (p.type == "two") {
    # Two-sided: use M^2, M^4, M^6, M^8 terms
    GM_cross <- (M^2 / 2 * res1.coeff$coeff2 %*% t(res2.coeff$coeff2) +
                 M^4 / 24 * res1.coeff$coeff4 %*% t(res2.coeff$coeff4) +
                 M^6 / 720 * res1.coeff$coeff6 %*% t(res2.coeff$coeff6) +
                 M^8 / 40320 * res1.coeff$coeff8 %*% t(res2.coeff$coeff8))

    if (var.correct == TRUE) {
      # Compute exact marginal variances
      v1 <- (res1.coeff$coeff2)^2 / 2 +
            (res1.coeff$coeff4)^2 / 24 +
            (res1.coeff$coeff6)^2 / 720 +
            (res1.coeff$coeff8)^2 / 40320

      v2 <- (res2.coeff$coeff2)^2 / 2 +
            (res2.coeff$coeff4)^2 / 24 +
            (res2.coeff$coeff6)^2 / 720 +
            (res2.coeff$coeff8)^2 / 40320

      # Apply variance correction
      GM_cross <- diag(sqrt(2 * D1 / v1), n, n) %*% GM_cross %*% diag(sqrt(2 * D2 / v2), n, n)
    }

  } else {
    # One-sided: use M, M^2, M^3, M^4 terms
    GM_cross <- (M * res1.coeff$coeff1 %*% t(res2.coeff$coeff1) +
                 M^2 / 2 * res1.coeff$coeff2 %*% t(res2.coeff$coeff2) +
                 M^3 / 6 * res1.coeff$coeff3 %*% t(res2.coeff$coeff3) +
                 M^4 / 24 * res1.coeff$coeff4 %*% t(res2.coeff$coeff4))

    if (var.correct == TRUE) {
      # Compute exact marginal variances
      v1 <- (res1.coeff$coeff1)^2 +
            (res1.coeff$coeff2)^2 / 2 +
            (res1.coeff$coeff3)^2 / 6 +
            (res1.coeff$coeff4)^2 / 24

      v2 <- (res2.coeff$coeff1)^2 +
            (res2.coeff$coeff2)^2 / 2 +
            (res2.coeff$coeff3)^2 / 6 +
            (res2.coeff$coeff4)^2 / 24

      # Apply variance correction
      GM_cross <- diag(sqrt(2 * D1 / v1), n, n) %*% GM_cross %*% diag(sqrt(2 * D2 / v2), n, n)
    }
  }

  # Return weighted sum of covariances
  return(sum(GM_cross * (W1 %*% t(W2))))
}


#' Calculate Correlation Matrix Between Multiple GFisher Statistics
#'
#' @title Correlation Matrix for Multiple GFisher Statistics (Corollary 2)
#' @description Calculate the correlation matrix between multiple GFisher statistics,
#'   each potentially having different degrees of freedom and weights.
#'
#' @param DD An \eqn{m \times n} matrix of degrees of freedom, where \eqn{m} is the number
#'   of GFisher statistics and \eqn{n} is the number of p-values combined by each GFisher test.
#' @param W An \eqn{m \times n} matrix of weights, where \eqn{m} is the number of GFisher
#'   statistics and \eqn{n} is the number of p-values combined by each GFisher test.
#' @param M Correlation matrix of the input Z-scores from which the input p-values were obtained.
#' @param p.type Character string: \code{"two"} for two-sided (default), \code{"one"} for one-sided input p-values.
#' @param var.correct Logical, passed to \code{getGFishercov()}. Default \code{TRUE}.
#'
#' @return An \eqn{m \times m} correlation matrix between the GFisher statistics \eqn{T^{(1)}, T^{(2)}, ..., T^{(m)}}.
#'
#' @details
#' This function computes the correlation matrix among multiple GFisher statistics by:
#' \enumerate{
#'   \item Computing all pairwise covariances using \code{\link{getGFishercov}}
#'   \item Converting the covariance matrix to a correlation matrix
#' }
#'
#' The correlation matrix is used in the minimum p-value approach of the oGFisher test
#' to compute the p-value via multivariate normal distribution.
#'
#' Each row of \code{DD} and \code{W} represents the configuration (degrees of freedom and weights)
#' for one GFisher test.
#'
#' @references
#' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and accurate p-value
#' calculation under dependence. \emph{Biometrics}, 79(2), 1159-1172. See Corollary 2.
#'
#' @keywords internal
#' @author Hong Zhang
getGFisherCOR <- function(DD, W, M, var.correct = TRUE, p.type = "two") {
  # Get number of GFisher tests
  m <- dim(DD)[1]

  # Initialize covariance matrix
  COV <- matrix(NA, ncol = m, nrow = m)

  # Compute upper triangle of covariance matrix
  for (i in 1:m) {
    for (j in i:m) {
      COV[i, j] <- getGFishercov(DD[i, ], DD[j, ], W[i, ], W[j, ],
                                  M, var.correct = var.correct, p.type = p.type)
    }
  }

  # Fill lower triangle (symmetric matrix)
  COV[lower.tri(COV)] <- t(COV)[lower.tri(COV)]

  # Convert covariance to correlation matrix
  return(stats::cov2cor(COV))
}


#' Calculate Covariance Matrix for Weighted Components
#'
#' @title Covariance Matrix for Weighted Components (Theorem 1 / Corollary 1)
#' @description Calculate the covariance matrix for the weighted components \eqn{w_1 T_1, ..., w_n T_n}
#'   of a single GFisher statistic. Automatically uses C++ implementation when available.
#'
#' @param D An n-dimensional vector of degrees of freedom.
#' @param w An n-dimensional vector of weights.
#' @param M Correlation matrix of the input Z-scores from which the input p-values were obtained.
#' @param p.type Character string: \code{"two"} for two-sided (default), \code{"one"} for one-sided input p-values.
#' @param use.cpp Logical, default \code{TRUE}. If \code{TRUE}, uses fast C++ implementation.
#'   Set to \code{FALSE} to force pure R implementation.
#'
#' @return An \eqn{n \times n} covariance matrix for the weighted components \eqn{w_1 T_1, ..., w_n T_n}.
#'
#' @details
#' This function calculates the covariance matrix for the components of a GFisher statistic,
#' which is needed for computing the eigenvalues in the quadratic approximation method.
#'
#' The covariance matrix is computed using Hermite polynomial expansions:
#'
#' \strong{For two-sided p-values:}
#' \deqn{GM_{ij} = \text{Corr}(T_i, T_j) \cdot \sqrt{2d_i} w_i \cdot \sqrt{2d_j} w_j}
#'
#' where the correlation is approximated using the coefficients from \code{\link{getGFishercoef}}.
#'
#' \strong{For one-sided p-values:}
#'
#' Similar formula with additional odd-order terms.
#'
#' The output of this function is the target covariance matrix \eqn{M} mentioned in
#' Section 3.4 (Quadratic approximation) of the GFisher paper.
#'
#' \strong{Performance Note:}
#' When \code{use.cpp = TRUE} (default), this function uses a C++ implementation with
#' composite Gauss-Legendre quadrature for numerical integration, providing substantial
#' speedup over the pure R implementation. The C++ version automatically falls back to
#' R if unavailable.
#'
#' @references
#' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and accurate p-value
#' calculation under dependence. \emph{Biometrics}, 79(2), 1159-1172. See Theorem 1 and Section 3.4.
#'
#' @keywords internal
#' @author Hong Zhang
getGFisherGM <- function(D, w, M, p.type = "two", use.cpp = TRUE) {
  # Try C++ version first if requested
  if (use.cpp) {
    tryCatch({
      # Tier 1: Use cached version for better performance with repeated df values
      return(getGFisherGM_cached_cpp(D, w, M, p.type))
    }, error = function(e) {
      # Tier 2: If cached version fails, try regular C++ version
      tryCatch({
        return(getGFisherGM_cpp(D, w, M, p.type))
      }, error = function(e2) {
        warning("C++ version failed, falling back to R implementation: ", e2$message)
      })
    })
  }

  # Tier 3: Pure R implementation (fallback)
  # Get dimension
  n <- dim(M)[1]

  # Reuse cached coefficients from getGFishercoef
  # This eliminates redundant integration calculations
  coeff <- getGFishercoef(D, M, p.type = p.type, use.cpp = FALSE)

  if (p.type == "two") {
    # Two-sided: use cached coefficients
    # Compute covariance matrix
    GM <- M^2 / 2 * coeff$coeff2 %*% t(coeff$coeff2) +
          M^4 / 24 * coeff$coeff4 %*% t(coeff$coeff4) +
          M^6 / 720 * coeff$coeff6 %*% t(coeff$coeff6) +
          M^8 / 40320 * coeff$coeff8 %*% t(coeff$coeff8)

    # Apply weights and scale to get covariance
    GM <- diag(sqrt(2 * D) * w, n, n) %*% stats::cov2cor(GM) %*% diag(sqrt(2 * D) * w, n, n)

  } else {
    # One-sided: use cached coefficients
    # Compute covariance matrix
    GM <- M * coeff$coeff1 %*% t(coeff$coeff1) +
          M^2 / 2 * coeff$coeff2 %*% t(coeff$coeff2) +
          M^3 / 6 * coeff$coeff3 %*% t(coeff$coeff3) +
          M^4 / 24 * coeff$coeff4 %*% t(coeff$coeff4)

    # Apply weights and scale to get covariance
    GM <- diag(sqrt(2 * D) * w, n, n) %*% stats::cov2cor(GM) %*% diag(sqrt(2 * D) * w, n, n)
  }

  return(GM)
}


#' Calculate Eigenvalues for Quadratic Approximation
#'
#' @title Eigenvalues for Quadratic Approximation (HYB Method)
#' @description Calculate eigenvalues needed for the quadratic approximation method,
#'   which only works for two-sided p-values. Automatically uses C++ implementation
#'   when available for 2-5x speedup.
#'
#' @param D An n-dimensional vector of degrees of freedom.
#' @param w An n-dimensional vector of weights.
#' @param M Correlation matrix of the input Z-scores from which the input p-values were obtained.
#' @param GM Covariance matrix between \eqn{w_1 T_1, ..., w_n T_n}, typically the output
#'   from \code{\link{getGFisherGM}}.
#' @param use.cpp Logical, default \code{TRUE}. If \code{TRUE}, uses fast C++ implementation.
#'   Set to \code{FALSE} to force pure R implementation.
#'
#' @return A list with:
#'   \item{lam}{Vector of eigenvalues (positive values > 1e-10) used in the quadratic approximation}
#'
#' @details
#' This function implements the eigenvalue calculation described in Section 3.4 of the GFisher paper.
#'
#' \strong{Algorithm:}
#' \enumerate{
#'   \item Construct the tilde correlation matrix \eqn{\tilde{M}} from formula (13):
#'     \deqn{\tilde{M}_{ij} = \text{sign}(M_{ij}) \sqrt{|GM_{ij}| / (\min(d_i, d_j) \cdot 2)}}
#'   \item Ensure \eqn{\tilde{M}} is a valid correlation matrix (all elements <= 1, positive definite)
#'   \item Compute weighted Cholesky: \eqn{WM = \text{chol}(\tilde{M}) \cdot \text{diag}(\sqrt{w})}
#'   \item Extract eigenvalues from \eqn{WM^T WM}
#'   \item For \eqn{d > 1}, add additional eigenvalues corresponding to higher degrees of freedom
#' }
#'
#' When weights are negative, eigenvalues can be negative. However, we keep the restriction
#' to positive eigenvalues to be consistent with the GFisher publication.
#'
#' \strong{Note:} This method only works for two-sided p-values. The quadratic approximation
#' leverages the chi-square representation under the null hypothesis.
#'
#' \strong{Performance Note:}
#' When \code{use.cpp = TRUE} (default), this function uses a C++ implementation with
#' Armadillo's optimized LAPACK routines for eigenvalue decomposition, providing substantial
#' speedup over the pure R implementation. The C++ version automatically falls back to
#' R if unavailable.
#'
#' @references
#' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and accurate p-value
#' calculation under dependence. \emph{Biometrics}, 79(2), 1159-1172. See Section 3.4 and Formula (13).
#'
#' @keywords internal
#' @author Hong Zhang
#' @importFrom Matrix nearPD
getGFisherlam <- function(D, w, M, GM, use.cpp = TRUE) {
  # Try C++ version first if requested
  if (use.cpp) {
    tryCatch({
      return(getGFisherlam_cpp(D, w, M, GM))
    }, error = function(e) {
      warning("C++ version failed, falling back to R implementation: ", e$message)
    })
  }

  # Pure R implementation (fallback)
  # Get dimension
  n <- dim(M)[1]

  # Compute min(d_i, d_j) in matrix form
  DM <- pmin(matrix(rep(D, n), nrow = n, byrow = FALSE),
             matrix(rep(D, n), nrow = n, byrow = TRUE))

  # Formula (13): construct tilde M
  M_tilde <- sqrt(abs(GM) / DM / 2) * sign(M)

  # Ensure all correlations are valid (<= 1)
  if (any(M_tilde > 1)) {
    M_tilde[M_tilde > 1] <- 0.999
  }

  # Ensure positive definiteness
  if (any(eigen(M_tilde)$values < 1e-10)) {
    M_tilde <- as.matrix(Matrix::nearPD(M_tilde, corr = TRUE)$mat)
  }

  # Weighted Cholesky decomposition
  M_tilde_chol <- chol(M_tilde)
  WM_chol <- M_tilde_chol %*% diag(sqrt(w))

  # Get eigenvalues from WM^T WM
  lam <- eigen(t(WM_chol) %*% WM_chol, symmetric = TRUE, only.values = TRUE)$values

  # For df > 1, add additional eigenvalues
  if (max(D) > 1) {
    for (i in 2:max(D)) {
      # Create indicator matrix for df >= i
      Ai <- diag(n)
      Did <- which(D < i)
      if (length(Did) > 0) {
        diag(Ai)[Did] <- 0
      }

      # Add eigenvalues from this configuration
      lam <- c(lam, eigen((WM_chol) %*% Ai %*% t(WM_chol),
                          symmetric = TRUE, only.values = TRUE)$values)
    }
  }

  # Keep only positive eigenvalues (> 1e-10)
  # This restriction is consistent with the GFisher publication
  lam <- lam[lam > 1e-10]

  return(list(lam = lam))
}
