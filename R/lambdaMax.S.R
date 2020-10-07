#' Maximum annual growth rate for short-lived species
#'
#' For internal use only. The function estimates maximum annual growth rate for short-lived species using the demographic invariants approach (Niel & Lebreton 2005).
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
#' @importFrom stats uniroot
#'
lambdaMax.S<-function(alpha,survival){
  lambdaMax <- function(x) {(exp((alpha+survival/(x-survival))^-1))-x}
  lambdaMax <- uniroot(lambdaMax, c(1,4)) # Looking for the value that comes closest to 0 (absolute value)
  lambdaMax <- lambdaMax[[1]] # Extracting the first element of the list
}
