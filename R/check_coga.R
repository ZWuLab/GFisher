# check_coga.R
# Utility function to check and install coga package

#' Check and Install coga Package for Independent Functions
#'
#' @description Helper function to check if the \code{coga} package is installed
#'   and optionally install it. The \code{coga} package is required for functions
#'   with \code{_ind} suffix (e.g., \code{\link{p.GFisher_ind}}, \code{\link{stat.oGFisher_ind}}).
#'
#' @param install Logical. If \code{TRUE} and \code{coga} is not installed,
#'   attempt to install it from CRAN. Default is \code{FALSE}.
#'
#' @return Logical (invisibly). \code{TRUE} if \code{coga} is installed (or successfully installed),
#'   \code{FALSE} otherwise.
#'
#' @details
#' The GFisher package includes optimized functions for calculating p-values when input
#' p-values are independent. These functions (with \code{_ind} suffix) provide faster and
#' more accurate calculations by leveraging the \code{coga} package for weighted sums of
#' chi-square distributions.
#'
#' Functions requiring \code{coga}:
#' \itemize{
#'   \item \code{\link{p.GFisher_ind}} - Calculate p-value under independence
#'   \item \code{\link{p.GFisher_ind_w1}} - Calculate p-value with equal weights (does not require coga)
#'   \item \code{\link{stat.oGFisher_ind}} - Calculate oGFisher statistics under independence
#'   \item \code{\link{pval.oGFisher_ind}} - Calculate oGFisher p-value under independence
#' }
#'
#' If \code{coga} is not available, you can still use the general functions
#' (e.g., \code{\link{p.GFisher}}) with \code{M = diag(n)} to indicate independence.
#'
#' @examples
#' # Check if coga is available
#' check_coga()
#'
#' # Install coga if not available
#' \dontrun{
#' check_coga(install = TRUE)
#' }
#'
#' @seealso \code{\link{p.GFisher_ind}}, \code{\link{stat.oGFisher_ind}}
#'
#' @export
check_coga <- function(install = FALSE) {
  has_coga <- requireNamespace("coga", quietly = TRUE)

  if (!has_coga) {
    message("The 'coga' package is not installed.")
    message("")
    message("The following functions require 'coga':")
    message("  - p.GFisher_ind")
    message("  - stat.oGFisher_ind")
    message("  - pval.oGFisher_ind")
    message("")
    message("Alternative: Use general functions with M = diag(n) for independence.")

    if (install) {
      message("")
      message("Installing 'coga' package from CRAN...")
      utils::install.packages("coga")
      has_coga <- requireNamespace("coga", quietly = TRUE)
      if (has_coga) {
        message("Successfully installed 'coga' package.")
      } else {
        warning("Failed to install 'coga' package. Please try manually: install.packages('coga')")
      }
    } else {
      message("")
      message("To install, run one of:")
      message("  install.packages('coga')")
      message("  GFisher::check_coga(install = TRUE)")
    }
  } else {
    message("'coga' package is installed and available.")
  }

  invisible(has_coga)
}
