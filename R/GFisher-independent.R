# GFisher-independent.R
# Optimized implementations for GFisher tests when input p-values are independent.
# These functions are substantially faster than using M as the identity matrix.

# Functions for independent input p-values


#' Calculate P-Value of GFisher Under Independence
#'
#' @title GFisher P-Value for Independent Inputs
#' @description Calculate the p-value of a GFisher test statistic when input p-values are
#'   independent. This is substantially faster than using \code{\link{p.GFisher}} with
#'   \code{M} as the identity matrix.
#'
#' @param q Numeric value of the observed GFisher statistic.
#' @param df Vector of degrees of freedom for each p-value transformation. Cannot be single value
#'   if the number of independent chi-square variables is greater than 1.
#' @param w Vector of non-negative weights. If single value, it is expanded to a vector of equal weights.
#'   \strong{Note:} Weights are automatically normalized to sum to 1, consistent with all functions in this package.
#' @param isExact Logical. If \code{TRUE} (default), uses exact calculation via the \code{coga} package.
#'   If \code{FALSE}, uses an approximation which may be less accurate but faster.
#'
#' @return A numeric p-value.
#'
#' @details
#' Under independence, the GFisher statistic follows a linear combination of independent
#' chi-square distributions:
#'
#' \deqn{S = \sum_{i=1}^n w_i \chi^2_{d_i}}
#'
#' This function computes \eqn{P(S > q)} using:
#' \itemize{
#'   \item \strong{Exact method (\code{isExact = TRUE})}: Uses the \code{coga} package to compute
#'     the distribution of a weighted sum of gamma random variables (chi-square is a special
#'     case of gamma). This method is accurate even for very small p-values.
#'   \item \strong{Approximation (\code{isExact = FALSE})}: Uses a faster approximation from
#'     the \code{coga} package, which may be less accurate for extreme tail probabilities.
#' }
#'
#' Weights are automatically scaled so that \eqn{\sum w_i = 1}.
#'
#' \strong{Note on dependencies:}
#'
#' This function requires the \code{coga} package, which is listed in \code{Suggests}.
#' If \code{coga} is not installed, the function will fail with an error message.
#'
#' @references
#' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and accurate p-value
#' calculation under dependence. \emph{Biometrics}, 79(2), 1159-1172.
#'
#' Witkovsky, V. (2016). Numerical inversion of a characteristic function: An alternative
#' tool to form the probability distribution of output quantity in linear measurement models.
#' \emph{Acta IMEKO}, 5(3), 32-44.
#' (Describes the characteristic function inversion methodology underlying the \code{coga} package)
#'
#' @examples
#' \dontrun{
#' # Requires coga package
#' if (requireNamespace("coga", quietly = TRUE)) {
#'   set.seed(123)
#'   n <- 10
#'   df <- runif(n, 0.5, 5)
#'   w <- abs(rnorm(n))
#'   q <- 40
#'
#'   # Exact calculation
#'   p_exact <- p.GFisher_ind(q = q, df = df, w = w, isExact = TRUE)
#'   print(p_exact)
#'
#'   # Approximation (faster but less accurate)
#'   p_approx <- p.GFisher_ind(q = q, df = df, w = w, isExact = FALSE)
#'   print(p_approx)
#'
#'   # Compare with general method (should be similar but slower)
#'   p_general <- p.GFisher(q = q, df = df, w = w, M = diag(n),
#'                          p.type = "two", method = "HYB")
#'   print(p_general)
#' }
#' }
#'
#' @seealso \code{\link{p.GFisher}}, \code{\link{p.GFisher_ind_w1}}, \code{\link{stat.oGFisher_ind}}
#'
#' @author Hong Zhang, Zheyang Wu
#' @export
p.GFisher_ind <- function(q, df, w, isExact = TRUE) {
  # Check if coga package is available
  if (!requireNamespace("coga", quietly = TRUE)) {
    stop("Package 'coga' is required for p.GFisher_ind but is not installed.\n",
         "Install it with: install.packages('coga')")
  }

  # Standardize and normalize weights to match stat.GFisher behavior
  # If w is a single value, expand to vector of equal weights
  if (length(w) == 1) {
    w <- rep(1, length(df))
  }

  # Normalize weights so sum(w) = 1 (consistent with stat.GFisher)
  w <- w / sum(w)

  # Convert chi-square to gamma parameters
  # Chi-square(df) = Gamma(shape = df/2, rate = 1/2)
  # For weighted sum: w * Chi-square(df) = Gamma(shape = df/2, rate = 1/(2*w))
  shape <- df / 2
  rate <- 1 / (2 * w)

  # Calculate p-value using coga package
  if (isExact) {
    # Exact calculation (recommended)
    pval <- 1 - coga::pcoga(x = q, shape = shape, rate = rate)
  } else {
    # Approximation (faster but may be less accurate)
    pval <- 1 - coga::pcoga_approx(x = q, shape = shape, rate = rate)
  }

  return(pval)
}


#' Calculate P-Value of GFisher Under Independence with Equal Weights
#'
#' @title GFisher P-Value for Independent Inputs with Equal Weights
#' @description Fast calculation of GFisher p-value when input p-values are independent
#'   and have equal weights (or weight = 1). Uses the fact that the sum of independent
#'   chi-square random variables is also chi-square distributed.
#'
#' @param q Numeric value of the observed GFisher statistic.
#' @param df Vector of degrees of freedom for each p-value transformation.
#'
#' @return A numeric p-value.
#'
#' @details
#' When weights are equal (or all equal to \eqn{1/n}), the GFisher statistic under independence is:
#'
#' \deqn{S = \frac{1}{n} \sum_{i=1}^n \chi^2_{d_i} \sim \frac{1}{n} \chi^2_{\sum d_i}}
#'
#' @references
#' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and accurate p-value
#' calculation under dependence. \emph{Biometrics}, 79(2), 1159-1172.
#'
#' @examples
#' # Example: Fast calculation with equal weights
#' set.seed(123)
#' n <- 10
#' df <- rep(2, n)  # Equal df
#' pval <- runif(n)
#' q <- stat.GFisher(pval, df = df, w = 1)
#'
#' # Fast exact p-value
#' p_fast <- p.GFisher_ind_w1(q, df)
#' print(p_fast)
#'
#' # Compare with standard chi-square test (df = 2n)
#' p_chisq <- pchisq(q, df = sum(df), lower.tail = FALSE)
#' print(p_chisq)
#'
#' @seealso \code{\link{p.GFisher_ind}}, \code{\link{p.GFisher}}
#'
#' @author Hong Zhang, Zheyang Wu
#' @export
p.GFisher_ind_w1 <- function(q, df) {
  # Sum of independent chi-square distributions is chi-square with summed df
  stats::pgamma(q, shape=sum(df)/2, scale=2/length(df), lower.tail=FALSE)
}


#' Calculate oGFisher Statistics Under Independence
#'
#' @title oGFisher Statistics for Independent Inputs
#' @description Calculate oGFisher statistics when input p-values are independent.
#'   This is the gold standard for independent inputs and is substantially faster than
#'   using \code{\link{stat.oGFisher}} with \code{M = diag(n)}.
#'
#' @param p A numeric vector of input p-values.
#' @param DF A matrix of degrees of freedom. Each row is the df vector for a GFisher test.
#' @param W A matrix of weights. Each row is the weight vector for a GFisher test.
#'
#' @return A list with the following components:
#'   \item{STAT}{Vector of GFisher test statistics, one for each row of \code{DF}}
#'   \item{PVAL}{Vector of individual p-values for each GFisher test}
#'   \item{minp}{Minimum p-value among all GFisher tests}
#'   \item{cct}{Cauchy combination test statistic for combining the p-values}
#'
#' @details
#' This function is the gold standard for computing oGFisher statistics when inputs are independent.
#'
#' It computes:
#' \enumerate{
#'   \item Individual GFisher statistics using \code{\link{stat.GFisher}}
#'   \item P-values using \code{\link{p.GFisher_ind}} (fast independent method)
#'   \item Minimum p-value and Cauchy combination statistic
#' }
#'
#' The Cauchy combination statistic handles very small p-values (< 1e-15) using the approximation
#' \eqn{1/(P_i \cdot \pi)}.
#'
#' @references
#' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and accurate p-value
#' calculation under dependence. \emph{Biometrics}, 79(2), 1159-1172.
#'
#' @examples
#' \dontrun{
#' # Requires coga package
#' if (requireNamespace("coga", quietly = TRUE)) {
#'   set.seed(123)
#'   n <- 10
#'   nGF <- 3
#'
#'   # Create test configurations
#'   DF <- matrix(runif(n * nGF, 0.5, 5), ncol = n) / 10
#'   W <- abs(matrix(rnorm(n * nGF), ncol = n))
#'
#'   # Under H0
#'   p <- runif(n)
#'   result_null <- stat.oGFisher_ind(p = p, DF = DF, W = W)
#'   print(result_null)
#'
#'   # With potential signal
#'   p[1] <- 0.00001
#'   result_signal <- stat.oGFisher_ind(p = p, DF = DF, W = W)
#'   print(result_signal)
#'
#'   # Compare with general method (should be similar but slower)
#'   result_general <- stat.oGFisher(p = p, DF = DF, W = W, M = diag(n), method = "HYB")
#'   print(result_general)
#' }
#' }
#'
#' @seealso \code{\link{stat.oGFisher}}, \code{\link{p.GFisher_ind}}, \code{\link{pval.oGFisher_ind}}
#'
#' @author Hong Zhang, Zheyang Wu
#' @export
stat.oGFisher_ind <- function(p, DF, W) {
  # Calculate individual GFisher statistics for each row of DF and W
  STAT <- sapply(1:dim(DF)[1], function(x) stat.GFisher(p = p, df = DF[x, ], w = W[x, ]))

  # Calculate p-values using fast independent method
  PVAL <- sapply(1:length(STAT), function(x) p.GFisher_ind(q = STAT[x], df = DF[x, ], w = W[x, ]))

  # Calculate minimum p-value
  minp <- min(PVAL)

  # Calculate Cauchy combination statistic (CCT)
  # Handle very small p-values specially
  is.small <- (PVAL < 1e-15)
  CCTSTAT <- PVAL

  # Transform p-values for CCT
  CCTSTAT[!is.small] <- tan((0.5 - CCTSTAT[!is.small]) * pi)
  CCTSTAT[is.small] <- 1 / CCTSTAT[is.small] / pi

  # Mean of transformed statistics
  cct <- mean(CCTSTAT)

  return(list(STAT = STAT, PVAL = PVAL, minp = minp, cct = cct))
}


#' Calculate oGFisher P-Value Under Independence
#'
#' @title oGFisher P-Value for Independent Inputs
#' @description Calculate the p-value for an omnibus GFisher test when input p-values are independent.
#'   This is the gold standard for independent inputs.
#'
#' @param p A numeric vector of input p-values.
#' @param DF A matrix of degrees of freedom. Each row is the df vector for a GFisher test.
#' @param W A matrix of weights. Each row is the weight vector for a GFisher test.
#' @param combine Character string: \code{"cct"} for Cauchy combination (default), or
#'   \code{"minp"} for minimum p-value with multivariate normal distribution.
#'
#' @return A list with the following components:
#'   \item{pval}{The p-value of the oGFisher test}
#'   \item{pval_indi}{Vector of individual p-values for each GFisher test}
#'
#' @details
#' This function is the gold standard for oGFisher tests under independence.
#'
#' \strong{Cauchy Combination (\code{combine = "cct"}):}
#'
#' Uses the Cauchy combination test statistic. The p-value is:
#' \deqn{P = P(\text{Cauchy} > \text{cct})}
#'
#' For extremely large CCT values (> 1e15), uses the approximation \eqn{P \approx 1/(\text{cct} \cdot \pi)}.
#'
#' \strong{Minimum P-Value (\code{combine = "minp"}):}
#'
#' Uses the correlation structure among GFisher tests (computed via \code{\link{getGFisherCOR}}
#' with \code{M = I}) to calculate the p-value of the minimum p-value statistic via
#' multivariate normal distribution.
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
#' \dontrun{
#' # Requires coga package
#' if (requireNamespace("coga", quietly = TRUE)) {
#'   set.seed(122)
#'   n <- 10
#'   nGF <- 20
#'
#'   # Create test configurations
#'   DF <- matrix(runif(n * nGF, 0.5, 5), ncol = n) / 10
#'   W <- abs(matrix(rnorm(n * nGF), ncol = n))
#'
#'   # Test data
#'   p <- runif(n)
#'   p[1] <- 0.00001
#'
#'   # CCT combination (default)
#'   result_cct <- pval.oGFisher_ind(p = p, DF = DF, W = W, combine = "cct")
#'   print(result_cct)
#'
#'   # Minimum p-value combination
#'   result_minp <- pval.oGFisher_ind(p = p, DF = DF, W = W, combine = "minp")
#'   print(result_minp)
#'
#'   # Compare with general method
#'   result_general_cct <- pval.oGFisher(p = p, DF = DF, W = W, M = diag(n),
#'                                       method = "MR", combine = "cct")
#'   print(result_general_cct)
#' }
#' }
#'
#' @seealso \code{\link{pval.oGFisher}}, \code{\link{stat.oGFisher_ind}}, \code{\link{p.GFisher_ind}}
#'
#' @author Hong Zhang, Zheyang Wu
#' @import mvtnorm
#' @export
pval.oGFisher_ind <- function(p, DF, W, combine = "cct") {
  # Calculate oGFisher statistics under independence
  out <- stat.oGFisher_ind(p, DF, W)

  if (combine == "cct") {
    # Cauchy combination test
    thr <- out$cct

    # Calculate p-value
    # For extremely large values, use approximation to avoid numerical issues
    if (thr > 1e+15) {
      pval <- (1 / thr) / pi
    } else {
      pval <- stats::pcauchy(thr, lower.tail = FALSE)
    }

  } else {
    # Minimum p-value with MVN distribution
    thr <- out$minp
    nd <- dim(DF)[1]

    # Calculate correlation matrix among GFisher tests under independence
    COR_GFisher <- getGFisherCOR(DD = DF, W = W, M = diag(length(p)))

    # Calculate p-value using multivariate normal distribution
    pval <- 1 - mvtnorm::pmvnorm(lower = rep(-Inf, nd),
                                  upper = stats::qnorm(1 - thr),
                                  mean = rep(0, nd),
                                  corr = COR_GFisher)[1]
  }

  return(list(pval = pval, pval_indi = out$PVAL))
}
