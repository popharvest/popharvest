#' Maximum annual growth rate for long-lived species
#'
#' For internal use only. The function estimates maximum annual growth rate for long-lived species using the demographic invariants approach (Niel & Lebreton 2005).
#'
#' @param alpha number ; age at first breeding
#' @param survival number ; annual survival rate of the adult age class
#'
#' @return A number
#' @export
#' @seealso \code{\link{PEG}}, \code{\link{PTL}}
#' @references
#' Niel, C., & Lebreton, J. D. (2005). Using demographic invariants to detect overharvested bird populations from incomplete data. Conservation Biology, 19(3), 826-835.
#'
lambdaMax.L <- function(alpha, survival){
  ((survival*alpha - survival + alpha+1) + sqrt((survival-survival*alpha-alpha-1)^2-4*survival*(alpha)^2)) / (2*alpha)
}
