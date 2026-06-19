# GFisher-main.R
# Main GFisher test functions for calculating statistics and p-values

# Main functions for users


#' Compute the GFisher Test Statistic
#'
#' @title GFisher Test Statistic
#' @description Compute the GFisher test statistics \eqn{S = \sum_i w_i F^{-1}_{d_i}(1-P_i)},
#'   based on a vector of p-values \eqn{P_i}'s, degrees of freedom \eqn{d_i}'s, and weights \eqn{w_i}'s.
#'   \eqn{F^{-1}_{d_i}} is the inverse CDF of the chi-square distribution with \eqn{d_i} degrees of freedom.
#'
#' @param p A numeric vector of input p-values for the GFisher test. Must be between 0 and 1.
#' @param df Degrees of freedom for inverse chi-square transformation for each p-value.
#'   It can be a vector of the same length as \code{p}, indicating each transformation function
#'   might have a different df, or a single number, indicating the same degrees of freedom for all.
#'   Default is 2 (standard Fisher's method).
#' @param w A numeric vector of non-negative weights for each p-value. Default is 1 (equal weights).
#'   If a single value is provided, equal weights are used for all p-values.
#'   \strong{Note:} Weights are automatically normalized to sum to 1 (i.e., \eqn{w \leftarrow w/\sum(w)}).
#'   This normalization improves numerical stability and is statistically equivalent to other
#'   normalizations such as \eqn{\sum w_i = n}.
#'
#' @return A numeric value representing the GFisher test statistic.
#'
#' @details
#' The statistic is calculated based on normalized weights, i.e., \eqn{w/\sum(w)}, ensuring that
#' the weights always sum to 1.
#'
#' When all degrees of freedom equal 2, the function uses a faster implementation equivalent to
#' Fisher's combination method: \eqn{-2\log(p)}.
#'
#' The GFisher test generalizes Fisher's combination method by allowing:
#' \itemize{
#'   \item Different degrees of freedom for each p-value transformation
#'   \item Flexible weighting schemes
#'   \item Handling of dependent test statistics through correlation structure (see \code{\link{p.GFisher}})
#' }
#'
#' @references
#' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and accurate p-value
#' calculation under dependence. \emph{Biometrics}, 79(2), 1159-1172.
#'
#' @examples
#' # Example 1: Equal weights with df=2 (standard Fisher's method)
#' set.seed(123)
#' n <- 10
#' pval <- runif(n)
#' stat.GFisher(pval, df = 2, w = 1)
#'
#' # Example 2: Equal weights with explicit vectors
#' stat.GFisher(pval, df = rep(2, n), w = rep(1, n))
#'
#' # Example 3: Varying degrees of freedom and weights
#' stat.GFisher(pval, df = 1:n, w = 1:n)
#'
#' # Example 4: Small p-values (signal detection)
#' pval_signal <- c(0.001, 0.002, 0.05, runif(7))
#' stat.GFisher(pval_signal, df = 2, w = 1)
#'
#' @seealso \code{\link{p.GFisher}}, \code{\link{stat.oGFisher}}, \code{\link{pval.oGFisher}}
#'
#' @author Hong Zhang, Zheyang Wu
#' @export
stat.GFisher <- function(p, df = 2, w = 1) {
  # Input validation
  if (any(p < 0 | p > 1)) {
    stop("All p-values must be between 0 and 1")
  }
  if (any(df <= 0)) {
    stop("All degrees of freedom must be positive")
  }
  if (any(w < 0)) {
    stop("All weights must be non-negative")
  }

  # Transform p-values using inverse chi-square CDF
  # Special case: df=2 (standard Fisher's method) has faster implementation
  if (all(df == 2)) {
    # Faster implementation for Fisher's combination: -2*log(p)
    pp_trans <- -2 * log(p)
  } else {
    # General case: use qchisq for arbitrary degrees of freedom
    pp_trans <- stats::qchisq(p, df = df, lower.tail = FALSE)
  }

  # Apply weights
  if (length(w) > 1) {
    # Normalize weights to sum to 1
    w <- w / sum(w)  # Require non-negative weights; sum(w) should be > 0
    fisherstat <- sum(w * pp_trans)
  } else {
    # Equal weights (mean)
    fisherstat <- sum(pp_trans) / length(pp_trans)
  }

  return(fisherstat)
}


#' Compute the P-Value of a GFisher Test
#'
#' @title GFisher Test P-Value Calculation
#' @description Calculate the p-value of a GFisher test statistic under dependence,
#'   using moment matching methods to account for correlation structure among input statistics.
#'
#' @param q Numeric value of the observed GFisher statistic (typically from \code{\link{stat.GFisher}}).
#' @param df Vector of degrees of freedom for inverse chi-square transformation for each
#'   p-value. If all dfs are equal, it can be defined by a single constant.
#' @param w Vector of non-negative weights. If a single value, equal weights are used.
#'   \strong{Note:} Weights are automatically normalized to sum to 1.
#' @param M Correlation matrix of the Z-scores from which the input p-values were obtained.
#'   Must be a square matrix with dimension equal to the length of \code{df}.
#' @param p.type Character string specifying p-value type:
#'   \itemize{
#'     \item \code{"two"}: Two-sided input p-values (default)
#'     \item \code{"one"}: One-sided input p-values
#'   }
#' @param method Character string specifying the calculation method:
#'   \itemize{
#'     \item \code{"MR"}: Simulation-assisted moment ratio matching (recommended for tail probabilities)
#'     \item \code{"HYB"}: Moment ratio matching by quadratic approximation (default, faster but only for two-sided p-values)
#'     \item \code{"GB"}: Brown's method with calculated variance
#'   }
#' @param nsim Number of simulations used in the \code{"MR"} method. Default is \code{5e4}.
#'   Larger values provide more accurate estimates but take longer to compute.
#' @param seed Optional seed for random number generation in simulation-based methods.
#'   Default is \code{NULL} (no seed set).
#' @param use.cpp Logical indicating whether to use C++ implementation for method="MR".
#'   Default is \code{TRUE}. Set to \code{FALSE} to use pure R implementation.
#'
#' @return A numeric value representing the p-value of the GFisher test.
#'
#' @details
#' This function calculates the p-value for a GFisher test statistic under dependence
#' by approximating the null distribution using moment matching.
#'
#' \strong{Method Selection:}
#' \itemize{
#'   \item \strong{HYB (Hybrid)}: Uses quadratic approximation based on the first four moments.
#'     Faster but only works for two-sided p-values. Best for moderate test statistics.
#'   \item \strong{MR (Moment Ratio)}: Uses simulation to estimate moments, then applies
#'     moment matching. More accurate for extreme tail probabilities and works for both
#'     one-sided and two-sided p-values.
#'   \item \strong{GB (Generalized Brown)}: Uses Brown's method with calculated variance.
#'     Simplest approximation but may be less accurate.
#' }
#'
#' The function automatically handles near-singular correlation matrices by finding
#' the nearest positive definite matrix when necessary.
#'
#' Weights are automatically normalized to sum to 1.
#'
#' @references
#' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and accurate p-value
#' calculation under dependence. \emph{Biometrics}, 79(2), 1159-1172.
#'
#' @examples
#' # Example 1: Simple case with equal correlation
#' set.seed(123)
#' n <- 10
#' M <- matrix(0.3, n, n) + diag(0.7, n, n)  # Correlation matrix
#' zscore <- matrix(rnorm(n), nrow = 1) %*% chol(M)
#' pval <- 2 * (1 - pnorm(abs(zscore)))  # Two-sided p-values
#'
#' # Calculate statistic and p-value with df=2
#' gf1 <- stat.GFisher(pval, df = 2, w = 1)
#' p.GFisher(gf1, df = 2, w = 1, M = M, method = "HYB")
#' p.GFisher(gf1, df = 2, w = 1, M = M, method = "MR", nsim = 5e4)
#'
#' # Example 2: Varying df and weights
#' gf2 <- stat.GFisher(pval, df = 1:n, w = 1:n)
#' p.GFisher(gf2, df = 1:n, w = 1:n, M = M, method = "HYB")
#' p.GFisher(gf2, df = 1:n, w = 1:n, M = M, method = "MR", nsim = 5e4)
#'
#' @seealso \code{\link{stat.GFisher}}, \code{\link{stat.oGFisher}}, \code{\link{pval.oGFisher}}
#'
#' @author Hong Zhang, Zheyang Wu
#' @import mvtnorm
#' @importFrom Matrix nearPD
#' @export
p.GFisher <- function(q, df, w, M, p.type = "two", method = "HYB", nsim = NULL, seed = NULL, use.cpp = TRUE) {
  # Input validation
  if (q < 0) {
    stop("q must be non-negative")
  }
  if (any(df <= 0)) {
    stop("All degrees of freedom must be positive")
  }
  if (any(w < 0)) {
    stop("All weights must be non-negative")
  }
  
  # Get dimension
  n <- dim(M)[1]

  # Standardize df and w to vectors
  if (length(df) == 1) {
    df <- rep(df, n)
  }
  if (length(w) == 1) {
    w <- rep(1 / n, n)
  }

  # Normalize weights
  w <- w / sum(w)

  # Method: MR (Moment Ratio with simulation)
  if (method == "MR") {
    # Set number of simulations
    if (!is.numeric(nsim)) {
      nsim <- 5e4
    }

    # Try C++ version first if enabled
    if (use.cpp) {
      tryCatch({
        # Ensure M is positive definite by Cholesky decomposition
        M_chol_test <- try(chol(M), silent = TRUE)

        # If Cholesky fails, find nearest positive definite matrix
        if ("try-error" %in% class(M_chol_test)) {
          M <- as.matrix(Matrix::nearPD(M, corr = TRUE)$mat)
        }

        # Set seed if provided (must be done at R level for C++ to use R's RNG)
        if (is.numeric(seed)) {
          set.seed(seed)
        }

        # Call C++ version (pass 0 for seed since it's already set in R)
        return(p_GFisher_MR_cpp(q, df, w, M, p.type, as.integer(nsim), 0L))

      }, error = function(e) {
        warning("C++ MR method failed, using R implementation: ", e$message)
      })
    }

    # R implementation (fallback or if use.cpp=FALSE)
    # Ensure M is positive definite by Cholesky decomposition
    M_chol <- try(chol(M), silent = TRUE)

    # If Cholesky fails, find nearest positive definite matrix
    if ("try-error" %in% class(M_chol)) {
      M <- as.matrix(Matrix::nearPD(M, corr = TRUE)$mat)
      M_chol <- chol(M)
    }

    # Set seed if provided
    if (is.numeric(seed)) {
      set.seed(seed)
    }

    # Generate null distribution via simulation
    znull <- matrix(stats::rnorm(n * nsim), ncol = n, nrow = nsim) %*% M_chol

    # Convert to p-values based on type
    if (p.type == "two") {
      pnull <- 2 * stats::pnorm(-abs(znull), 0, 1)
    } else {
      pnull <- stats::pnorm(znull, 0, 1, lower.tail = FALSE)
    }

    # Transform p-values
    if (stats::sd(df) == 0) {
      # All df are the same
      if (all(df == 2)) {
        pp_trans <- -2 * log(pnull)
      } else {
        pp_trans <- stats::qchisq(pnull, df = df[1], lower.tail = FALSE)
      }
    } else {
      # Different df values
      pp_trans <- t(apply(pnull, 1, function(x) stats::qchisq(x, df = df, lower.tail = FALSE)))
    }

    # Calculate weighted statistics
    fishernull <- apply(pp_trans, 1, function(x) sum(x * w))

    # Compute moments
    MM <- sapply(1:4, function(x) mean(fishernull^x))
    mu <- MM[1]
    sigma2 <- (MM[2] - MM[1]^2) * nsim / (nsim - 1)
    cmu3 <- MM[3] - 3 * MM[2] * MM[1] + 2 * MM[1]^3
    cmu4 <- MM[4] - 4 * MM[3] * MM[1] + 6 * MM[2] * MM[1]^2 - 3 * MM[1]^4

    # Calculate gamma and kappa
    gm <- cmu3 / sigma2^(3 / 2)
    kp <- cmu4 / sigma2^2

    # Moment matching to gamma distribution
    a <- 9 * gm^2 / (kp - 3)^2
    pval <- stats::pgamma((q - mu) / sqrt(sigma2) * sqrt(a) + a,
                          shape = a, scale = 1, lower.tail = FALSE)

  } else {
    # Methods: HYB or GB (use theoretical covariance)
    GM <- getGFisherGM(df, w, M, p.type)
    mu <- sum(w * df)
    sigma2 <- sum(GM)

    if (method == "HYB") {
      # Hybrid method (quadratic approximation)
      # Only works for two-sided p-values
      if (p.type == "two") {
        out_GM <- getGFisherlam(df, w, M, GM)
        lam <- out_GM$lam

        # Calculate moments from eigenvalues
        c2 <- sum(lam^2)
        c3 <- sum(lam^3)
        c4 <- sum(lam^4)
        gm <- sqrt(8) * c3 / (c2)^(3 / 2)
        kp <- 12 * c4 / c2^2 + 3
        a <- 9 * gm^2 / (kp - 3)^2
      } else {
        stop("Hybrid method (quadratic approximation) only works for two-sided p-values.")
      }
    } else {
      # Brown's method
      a <- mu^2 / sigma2
    }

    # Calculate p-value using gamma approximation
    pval <- stats::pgamma((q - mu) / sqrt(sigma2) * sqrt(a) + a,
                          shape = a, scale = 1, lower.tail = FALSE)
  }

  return(pval)
}


#' Compute the oGFisher Test Statistics
#'
#' @title Omnibus GFisher (oGFisher) Test Statistics
#' @description Calculate statistics for multiple GFisher tests with different degrees of freedom
#'   and weights, along with their combination via Cauchy combination test (CCT) or minimum p-value.
#'
#' @param p A numeric vector of input p-values for the oGFisher test.
#' @param DF A matrix of degrees of freedom for inverse chi-square transformation for each
#'   p-value. Each row represents a GFisher test. It can be a matrix with one column, indicating
#'   the same degrees of freedom for all p-values's chi-square transformations across different tests.
#' @param W A matrix of non-negative weights. Each row represents a GFisher test.
#'   Must have the same dimensions as \code{DF}.
#' @param M Correlation matrix of the input Z-scores from which the input p-values were obtained.
#' @param p.type Character string: \code{"two"} for two-sided (default), \code{"one"} for one-sided input p-values.
#' @param method Character string specifying calculation method:
#'   \itemize{
#'     \item \code{"MR"}: Simulation-assisted moment ratio matching
#'     \item \code{"HYB"}: Moment ratio matching by quadratic approximation (default)
#'     \item \code{"GB"}: Brown's method with calculated variance
#'   }
#' @param nsim Number of simulations used in the \code{"MR"} method. Default is \code{5e4}.
#' @param seed Optional seed for random number generation. Default is \code{NULL}.
#'
#' @return A list with the following components:
#'   \item{STAT}{Vector of GFisher test statistics, one for each row of \code{DF}}
#'   \item{PVAL}{Vector of individual p-values for each GFisher test}
#'   \item{minp}{Minimum p-value among all GFisher tests}
#'   \item{cct}{Cauchy combination test statistic for combining the p-values}
#'
#' @details
#' The oGFisher (omnibus GFisher) test evaluates multiple GFisher statistics simultaneously,
#' each with potentially different degrees of freedom and weighting schemes. This allows
#' for a flexible and powerful approach to combining p-values.
#'
#' The function computes:
#' \enumerate{
#'   \item Individual GFisher statistics for each row of \code{DF} and \code{W}
#'   \item Corresponding p-values for each statistic
#'   \item Minimum p-value across all tests
#'   \item Cauchy combination statistic (CCT) for robust aggregation
#' }
#'
#' \strong{Cauchy Combination:}
#'
#' The CCT statistic is computed as \eqn{\bar{T} = \frac{1}{m}\sum_{i=1}^m \tan[\pi(0.5 - P_i)]},
#' where \eqn{P_i} are the individual GFisher p-values. This approach is robust to dependence
#' and handles very small p-values through a special transformation.
#'
#' P-values larger than 0.9 are capped at 0.9 to improve stability. Very small p-values
#' (< 1e-15) are handled using a special approximation: \eqn{1/(P_i \cdot \pi)}.
#'
#' @references
#' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and accurate p-value
#' calculation under dependence. \emph{Biometrics}, 79(2), 1159-1172.
#'
#' Liu, Y., & Xie, J. (2020). Cauchy combination test: a powerful test with analytic p-value
#' calculation under arbitrary dependency structures. \emph{Journal of the American Statistical
#' Association}, 115(529), 393-402.
#'
#' @examples
#' # Example: Multiple GFisher tests with different df and weights
#' set.seed(123)
#' n <- 10
#' M <- matrix(0.3, n, n) + diag(0.7, n, n)
#' zscore <- matrix(rnorm(n), nrow = 1) %*% chol(M)
#' pval <- 2 * (1 - pnorm(abs(zscore)))
#'
#' # Define two GFisher tests
#' DF <- rbind(rep(1, n), rep(2, n))
#' W <- rbind(rep(1, n), 1:10)
#'
#' # Calculate oGFisher statistics
#' result <- stat.oGFisher(pval, DF, W, M, p.type = "two", method = "HYB")
#' print(result)
#'
#' # Alternative: single df per test (expanded internally)
#' DF_short <- rbind(1, 2)
#' result2 <- stat.oGFisher(pval, DF = DF_short, W, M, p.type = "two", method = "HYB")
#'
#' @seealso \code{\link{stat.GFisher}}, \code{\link{p.GFisher}}, \code{\link{pval.oGFisher}}
#'
#' @author Hong Zhang, Zheyang Wu
#' @export
stat.oGFisher <- function(p, DF, W, M, p.type = "two", method = "HYB", nsim = NULL, seed = NULL) {
  # Calculate individual GFisher statistics for each row of DF and W
  STAT <- sapply(1:dim(DF)[1], function(x) stat.GFisher(p = p, df = DF[x, ], w = W[x, ]))

  # Calculate p-values for each statistic
  PVAL <- sapply(1:length(STAT), function(x) {
    p.GFisher(STAT[x], df = DF[x, ], w = W[x, ], M, p.type, method, nsim, seed)
  })

  # Calculate minimum p-value
  minp <- min(PVAL)

  # Cap large p-values at 0.9 for stability
  PVAL[PVAL > 0.9] <- 0.9

  # Calculate Cauchy combination statistic (CCT)
  # Handle very small p-values specially
  is.small <- (PVAL < 1e-15)
  CCTSTAT <- PVAL

  if (sum(is.small) == 0) {
    # No very small p-values
    cct <- mean(tan(pi * (0.5 - CCTSTAT)))
  } else {
    # Transform regular and very small p-values separately
    CCTSTAT[!is.small] <- tan((0.5 - CCTSTAT[!is.small]) * pi)
    CCTSTAT[is.small] <- 1 / CCTSTAT[is.small] / pi
    cct <- mean(CCTSTAT)
  }

  return(list(STAT = STAT, PVAL = PVAL, minp = minp, cct = cct))
}


#' Compute the oGFisher Test P-Value
#'
#' @title Omnibus GFisher (oGFisher) Test P-Value
#' @description Calculate the p-value for an omnibus GFisher test, which combines multiple
#'   GFisher statistics using either Cauchy combination or multivariate normal distribution.
#'
#' @param p A numeric vector of input p-values for the oGFisher test.
#' @param DF A matrix of degrees of freedom for inverse chi-square transformation for each
#'   p-value. Each row represents a GFisher test.
#' @param W A matrix of non-negative weights. Each row represents a GFisher test.
#' @param M Correlation matrix of the input Z-scores from which the input p-values were obtained.
#' @param p.type Character string: \code{"two"} for two-sided (default), \code{"one"} for one-sided input p-values.
#' @param method Character string specifying calculation method:
#'   \itemize{
#'     \item \code{"MR"}: Simulation-assisted moment ratio matching
#'     \item \code{"HYB"}: Moment ratio matching by quadratic approximation (default)
#'     \item \code{"GB"}: Brown's method with calculated variance
#'   }
#' @param combine Character string specifying combination method:
#'   \itemize{
#'     \item \code{"cct"}: oGFisher using the Cauchy combination test (default)
#'     \item \code{"mvn"}: oGFisher using multivariate normal distribution (minimum p-value approach)
#'   }
#' @param nsim Number of simulations used in the \code{"MR"} method. Default is \code{5e4}.
#' @param seed Optional seed for random number generation. Default is \code{NULL}.
#'
#' @return A list with the following components:
#'   \item{stat}{The test statistic: CCT statistic if \code{combine = "cct"}, or minimum p-value if \code{combine = "mvn"}}
#'   \item{pval}{The p-value of the oGFisher test}
#'   \item{pval_indi}{Vector of individual p-values for each GFisher test}
#'   \item{stat_indi}{Vector of individual GFisher statistics}
#'
#' @details
#' This function performs the omnibus GFisher test by:
#' \enumerate{
#'   \item Computing individual GFisher statistics and p-values for each test configuration
#'   \item Combining these p-values using either:
#'     \itemize{
#'       \item \strong{CCT (Cauchy Combination Test)}: Robust to dependence, analytically tractable
#'       \item \strong{MVN (Multivariate Normal)}: Uses correlation structure to compute minimum p-value distribution
#'     }
#' }
#'
#' \strong{Cauchy Combination Method (\code{combine = "cct"}):}
#'
#' The CCT approach is preferred when the individual tests may be highly dependent.
#' The p-value is calculated as \eqn{P(\text{Cauchy} > \text{cct})} where cct is the
#' Cauchy combination statistic from \code{\link{stat.oGFisher}}.
#'
#' \strong{Minimum P-Value Method (\code{combine = "mvn"}):}
#'
#' This approach accounts for the correlation among GFisher tests using their theoretical
#' correlation matrix (computed via \code{\link{getGFisherCOR}}). The p-value is:
#' \deqn{P = 1 - P(\text{all normalized p-values} > \Phi^{-1}(1-\text{minp}))}
#' where the probability is computed using the multivariate normal distribution.
#'
#' @references
#' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and accurate p-value
#' calculation under dependence. \emph{Biometrics}, 79(2), 1159-1172.
#'
#' Liu, Y., & Xie, J. (2020). Cauchy combination test: a powerful test with analytic p-value
#' calculation under arbitrary dependency structures. \emph{Journal of the American Statistical
#' Association}, 115(529), 393-402.
#'
#' @examples
#' # Example: oGFisher test with Cauchy combination
#' set.seed(123)
#' n <- 10
#' M <- matrix(0.3, n, n) + diag(0.7, n, n)
#' zscore <- matrix(rnorm(n), nrow = 1) %*% chol(M)
#' pval <- 2 * (1 - pnorm(abs(zscore)))
#'
#' # Define multiple GFisher tests
#' DF <- rbind(rep(1, n), rep(2, n))
#' W <- rbind(rep(1, n), 1:10)
#'
#' # CCT combination
#' result_cct <- pval.oGFisher(pval, DF, W, M, p.type = "two",
#'                              method = "HYB", combine = "cct")
#' print(result_cct)
#'
#' # MVN combination (minimum p-value)
#' result_mvn <- pval.oGFisher(pval, DF, W, M, p.type = "two",
#'                              method = "HYB", combine = "mvn")
#' print(result_mvn)
#'
#' # Alternative: single df per test
#' DF_short <- rbind(1, 2)
#' result <- pval.oGFisher(pval, DF = DF_short, W, M, p.type = "two",
#'                         method = "HYB", combine = "cct")
#'
#' @seealso \code{\link{stat.GFisher}}, \code{\link{p.GFisher}}, \code{\link{stat.oGFisher}}, \code{\link{getGFisherCOR}}
#'
#' @author Hong Zhang, Zheyang Wu
#' @import mvtnorm
#' @export
pval.oGFisher <- function(p, DF, W, M, p.type = "two", method = "HYB",
                          combine = "cct", nsim = NULL, seed = NULL) {
  # Calculate oGFisher statistics
  out <- stat.oGFisher(p = p, DF = DF, W = W, M = M, p.type = "two",
                       method = method, nsim = nsim, seed = seed)

  if (combine == "cct") {
    # Cauchy combination test
    stat <- out$cct
    pval <- stats::pcauchy(stat, lower.tail = FALSE)

  } else {
    # Minimum p-value with MVN distribution
    stat <- out$minp
    nd <- dim(DF)[1]

    # Calculate correlation matrix among GFisher tests
    COR_GFisher <- getGFisherCOR(DD = DF, W = W, M = M, p.type = p.type)

    # Calculate p-value using multivariate normal distribution
    pval <- 1 - mvtnorm::pmvnorm(lower = rep(-Inf, nd),
                                  upper = stats::qnorm(1 - stat),
                                  mean = rep(0, nd),
                                  corr = COR_GFisher)[1]
  }

  return(list(stat = stat, pval = pval, pval_indi = out$PVAL, stat_indi = out$STAT))
}
