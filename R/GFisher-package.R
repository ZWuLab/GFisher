#' GFisher: Generalized Fisher's Combination Tests Under Dependence
#'
#' This package provides accurate and computationally efficient methods for computing
#' p-values for a general family of Fisher-type statistics (GFisher), including
#' Fisher's combination, Good's statistic, and Lancaster's statistic.
#'
#' @section Weight Normalization Convention:
#' All functions in this package normalize weights so that \eqn{\sum w_i = 1}.
#' This normalization:
#' \itemize{
#'   \item Improves numerical stability in p-value calculations
#'   \item Is statistically equivalent to other normalizations (e.g., \eqn{\sum w_i = n}) in terms of its p-value results
#'   \item Focuses on hypothesis testing rather than arbitrary distributions
#' }
#' When you provide weights, they are automatically rescaled internally. For example,
#' \code{w = c(1, 2, 3)} is equivalent to \code{w = c(1/6, 2/6, 3/6)}.
#'
#' @author Zheyang Wu \email{zheyangwu@@wpi.edu}
#' @author Hong Zhang \email{consistencyzhang@@gmail.com}
#'
#' @references
#' Zhang, H. and Wu, Z. (2022). The generalized Fisher's combination and accurate
#' p-value calculation under dependence. \emph{Biometrics}, 79(2), 1159-1172.
#' \doi{10.1111/biom.13634}
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib GFisher, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
