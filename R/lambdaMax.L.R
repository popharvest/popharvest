#' Maximum annual growth rate for long-lived species
#'
#' \code{lambdaMax.L} function estimates maximum annual growth rate for long-lived species using the Demographic Invariants Method (DIM) (Niel &
#' Lebreton 2005).
#'
#' @details
#' The estimate of \eqn{lambdaMax} is based on the Demographic Invariants Method (DIM) (Niel & Lebreton, 2005) (see \code{rmax} function).
#' This approach estimates \eqn{lambdaMax} (Dillingham & Fletcher, 2008) on the basis of 2 parameters :
#' the maximum annual survival rate of adults (\code{surv} arguments) and age at first breeding (\code{alpha} arguments). When adult survival and age at first
#' breeding are known, users, can provided corresponding values either as point estimates (\code{surv.fixed}, \code{alpha.fixed} arguments) or
#' as distributions (\code{surv.beta = TRUE} with \code{mean.surv} and \code{sd.surv}, \code{alpha.unif = TRUE} with \code{min.alpha} and \code{max.alpha},
#' \code{alpha.lognorm = TRUE} with \code{mean.alpha} and \code{sd.alpha} arguments). Alternatively, annual survival rate can be estimated from body mass following the
#' relationship described in Johnson \emph{et al.} (2012) (see \code{surv} function).
#'
#' @param alpha number ; age at first breeding
#' @param survival number ; annual survival rate of adults
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
lambdaMax.L <- function(alpha, survival){
  ((survival*alpha - survival + alpha+1) + sqrt((survival-survival*alpha-alpha-1)^2-4*survival*(alpha)^2)) / (2*alpha)
}
