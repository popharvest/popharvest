#' Maximum annual growth rate for short-lived species
#'
#' \code{lambdaMax.S} function estimates maximum annual growth rate for short-lived species using the Demographic Invariants Method (DIM) (Niel &
#' Lebreton 2005).
#'
#' @details
#' The estimate of \eqn{lambdaMax} is based on the Demographic Invariants Method (DIM) (Niel & Lebreton, 2005) (see \code{rmax} function).
#' This approach calculates \eqn{lambdaMax} (Dillingham & Fletcher, 2008) on the basis of 2 parameters :
#' the maximum value of adult annual annual survival rate (\code{surv} arguments) and age at first breeding (\code{alpha} arguments). When adult survival and age at first
#' breeding are known from the user, they can be provided by either point estimates (\code{surv.fixed}, \code{alpha.fixed} arguments) or
#' distributions  (\code{surv.beta = TRUE} with \code{mean.surv} and \code{sd.surv}, \code{alpha.unif = TRUE} with \code{min.alpha} and \code{max.alpha},
#' \code{alpha.lognorm = TRUE} with \code{mean.alpha} and \code{sd.alpha} arguments). Also, adult annual annual survival rate can be estimated from body mass following the
#' relationship described in Johnson et al. (2012) (see \code{surv} function).
#' @param alpha number ; age at first breeding
#' @param survival number ; adult annual annual survival rate
#'
#' @return A number
#' @export
#' @seealso \code{\link{PEG}}, \code{\link{PTL}}
#' @references
#' Dillingham, P. W., & Fletcher, D. (2008). Estimating the ability of birds to sustain additional human-caused mortalities using a simple decision
#' rule and allometric relationships. Biological Conservation, 141(7), 1783-1792.
#'
#' Niel, C., & Lebreton, J. D. (2005). Using demographic invariants to detect overharvested bird populations from incomplete data. Conservation Biology, 19(3), 826-835.
#'
#' @importFrom stats uniroot
#'
lambdaMax.S<-function(alpha,survival){
  lambdaMax <- function(x) {(exp((alpha+survival/(x-survival))^-1))-x}
  lambdaMax <- uniroot(lambdaMax, c(1,4)) # Looking for the value that comes closest to 0 (absolute value)
  lambdaMax <- lambdaMax[[1]] # Extracting the first element of the list
}
