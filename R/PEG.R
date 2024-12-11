#' Potential Excess Growth estimation
#'
#' @description
#' \code{PEG} is a function used to assess whether a known level of harvesting exceeds the biological potential of renewable of a population, see 'Details'.
#'
#' @details
#' The function calculates the Potential Excess Growth (PEG) and returns a Sustainable Harvest Index (SHI), in absence of others mortality source.
#'
#'
#' PEG the maximum number of individuals that can be removed annually from a population given its rate (Rmax). PEG
#' is calculated as the product of the safety factor (Fs), Rmax, and the population size (popsize) : \eqn{Fs.Rmax.popsize} (Niel & Lebreton 2005).
#'
#' Fs (\code{Fs} argument) is a positive number less than 1 that takes into account several aspects such as the uncertainties associated
#' with estimates of population size and maximum annual population growth rate as well as the shape and magnitude of the underlying density dependent processes.
#' This is a "caution coefficient" that takes only a fraction of the potential excess growth.
#'
#'
#' Rmax is the maximum annual recruitment rate estimation. This is equivalent to \eqn{1-lambdaMax}, where lambdaMax is the maximum annual population growth rate,
#' under optimal growth conditions. Rmax and lambdaMax can be provided either as point estimates (\code{Rmax.fixed},
#' \code{lambdaMax.fixed}) or drawn from log-normal distributions (\code{Rmax.lognorm = TRUE} with \code{mean.Rmax} and \code{sd.Rmax},
#' \code{lambdaMax = TRUE} with \code{mean.lambdaMax} and \code{sd.lambdaMax}). When Rmax (or lambdaMax) is unknown,
#' the value can be estimated using the Demographic Invariants Method (DIM) (Niel & Lebreton, 2005) (see \code{rmax} function). This
#' approach calculates \eqn{Rmax = lambdaMax – 1} (Dillingham & Fletcher, 2008) on the basis of 2 parameters: annaul survival rate of adults under ideal condition
#' (maximal survival, \code{surv} argument) and age at first breeding (\code{alpha} arguments). When adult survival and age at first reproduction are
#' known, users can provided corresponding values either as point estimates (\code{surv.fixed}, \code{alpha.fixed} arguments) or as distributions
#' (\code{surv.beta = TRUE} with \code{mean.surv} and \code{sd.surv}, \code{alpha.unif = TRUE} with \code{min.alpha} and \code{max.alpha},
#' \code{alpha.lognorm = TRUE} with \code{mean.alpha} and \code{sd.alpha} arguments). Alternatively, adult survival rate can be estimated from body mass
#' following the relationship described in Johnson \emph{et al.} (2012) : \eqn{(p/surv.j^(1/(exp(3.22+0.24*log(body mass)+e)-alpha))} (see \code{surv}
#' function). Body mass can be provided either as a point estimate (\code{mass.fixed}) or drawn from a log-normal distribution
#' (\code{mass.lognorm = TRUE}, \code{mean.mass} and \code{sd.mass} arguments). The parameter p can be either drawn from a beta
#' distribution B(3.34, 101.24) (\code{type.p = "random"}) or provided as a point estimate equal to 0.032 (i.e. the mean of the beta distribution is
#' given by a/(a + b) with a = 3.34 and b = 101.24, Boitani & Fuller 2000) (\code{type.p = "determinist"}). The parameter e can be
#' either null (set to zero) (\code{type.e = "determinist"}) or residuals from a normal distribution N(0, 0.087) (\code{type.e = "random"}). See equations in Johnson \emph{et al.} (2012).
#'
#'
#' By default, the argument \code{surv.j} is silent and masked to users (the default value for \code{surv.j} is 1, as such \eqn{p/surv.j = p/1 = p} ; see original equation in
#' Johnson \emph{et al.}, 2012). This is correct for any alpha if, and only if, juvenile survival < adult survival solely for birds aged <1 year.  This is a reasonable
#' approximation especially in passerines and some other medium-sized bird species. For species with delayed sexual maturity and when juvenile survival < adult
#' survival for birds aged 1 to alpha, the default \code{surv} function returns biased estimates unless the users specify a value for annual survival for birds aged 1
#' to alpha. This is done by using the argument \code{surv.j} which is then used to estimate  \eqn{p = p/surv.j}.
#'
#'
#' In all cases, the type of living rate (\code{living.rate} argument), i.e. long-lived species (\code{living.rate = "long"}) or short-lived species (\code{living.rate = "short"}) must be quoted.
#' This argument is used to select the right formula from Niel & Lebreton (2005) for the calculation of lambdaMax (see \code{lambdaMax.L} for long-lived species and \code{lambdaMax.S} for short-lived species)
#' and, therefater Rmax.
#'
#'
#' The population size, as the previous terms, can also be provided either as point estimates (\code{pop.fixed}) or can be drawn from distributions (\code{pop.unif = TRUE}
#' with \code{min.pop} and \code{max.pop} or \code{pop.lognorm = TRUE} with \code{mean.pop} and \code{sd.pop}).
#'
#'
#' The Sustainable Harvest Index (SHI) is calculated as the ratio harvest level. The harvest level can
#' be provided either as point estimates (\code{harvest.fixed}) or distributions (\code{harvest.unif = TRUE}
#' with \code{min.harvest} and \code{max.harvest} or \code{harvest.lognorm = TRUE} with \code{mean.harvest} and \code{sd.harvest}).
#' When harvest exceeds PTL (SHI>1), this suggests an overharvesting regime with regard to management objectives.
#'
#'
#' @param full.option boolean, default \code{FALSE}. If \code{TRUE}, provides all parameter values necessary for calculations.
#' @param Nsim number, default 1 ;  number of simulations.
#' @param NSp number ; number of species for which calculations are done.
#' @param Rmax.fixed number ; fixed estimate of maximum annual recruitment rate provided by the user. To be used instead of \code{mean.Rmax} when \code{Rmax.lognorm = FALSE}.
#' @param Rmax.lognorm boolean, default \code{FALSE}. If \code{TRUE}, simulates Rmax from a log-normal distribution with mu = \code{mean.Rmax} and standard deviation = \code{sd.Rmax}.
#' @param mean.Rmax number ; arithmetic mean of maximum annual recruitment rate.
#' @param sd.Rmax number ; standard deviation of maximum annual recruitment rate.
#' @param lambdaMax.fixed number ; fixed estimate of maximum annual multiplicative growth rate provided by the user and used to calculate Rmax as \eqn{lamdaMax-1}.
#' @param lambdaMax.lognorm boolean, default \code{FALSE}. If \code{TRUE}, simulates maximum annual multiplicative growth rate from a log-normal distribution with mu = \code{mean.lambdaMax} and standard deviation = \code{sd.lambdaMax}.
#' @param mean.lambdaMax number ; arithmetic mean of maximum annual multiplicative growth rate. To be used only if \code{lamndaMax.lognorm = TRUE}.
#' @param sd.lambdaMax number ; standard deviation of maximum multiplicative growth rate. To be used only if \code{lambdaMax.lognorm = TRUE}.
#' @param living.rate character string ; living rate : either "\code{long}" or "\code{short}". No other character string allowed.
#' @param surv.fixed number ; fixed estimate of (annual) adult survival rate provided by the user. To be used instead of \code{mean.surv} when \code{surv.beta = FALSE}.
#' @param surv.beta boolean, default \code{FALSE}. If \code{TRUE}, simulates (annual) adult survival rate from a beta distribution with mu = \code{mean.surv} and standard deviation = \code{sd.surv} (see below).
#' @param mean.surv number ; arithmetic mean of (annual) adult survival rate. Only used if \code{surv.beta = TRUE}.
#' @param sd.surv number ; standard deviation of adult annual survival rate. Only used if \code{surv.beta = TRUE}.
#' @param alpha.fixed number ; fixed estimate of age at first breeding (in years). To be used when \code{alpha.unif = FALSE} and \code{alpha.lognorm = FALSE}.
#' @param alpha.unif boolean, default \code{FALSE}. Simulates age at first breeding from a uniform distribution bounded by \code{min.alpha} and \code{max.alpha}.
#' @param min.alpha number ; minimum age at first breeding. To be used only if \code{alpha.unif = TRUE}.
#' @param max.alpha number ; maximum age at first breeding. To be used only if \code{alpha.unif = TRUE}.
#' @param alpha.lognorm boolean, default \code{FALSE}. If \code{TRUE}, simulates age at first breeding from a log-normal distribution with mu = \code{mean.alpha} and standard deviation = \code{sd.alpha}.
#' @param mean.alpha number ; arithmetic mean of age at first breeding. To be used only used if \code{alpha.lognorm = TRUE}.
#' @param sd.alpha number ; age at first breeding standard deviation. Only used if \code{alpha.lognorm = TRUE}.
#' @param mass.fixed number ; fixed estimate of body mass (in kilograms) used to estimate adult survival rate under ideal conditions from body mass.
#' @param mass.lognorm boolean, default \code{FALSE}. If \code{TRUE}, simulates body mass from a log-normal distribution with mu = \code{mean.mass} and standard deviation = \code{sd.mass}.
#' @param mean.mass number ; arithmetic mean of body mass (in kilograms). Only used if \code{mass.lognorm = TRUE}.
#' @param sd.mass number ; standard deviation of body mass. Only used if \code{mass.lognorm = TRUE}.
#' @param type.p character string ; the value of the parameter p can be either "\code{determinist}" or "\code{random}". No other character string allowed.
#' @param type.e character string ; the value of the parameter e can be either "\code{determinist}" or "\code{random}". No other character string allowed.
#' @param surv.j number, default 1 ; point estimate of average annual survival for birds aged 1 to alpha. juvenile annual survival rate The value can be provided by the user . Used to estimate adult survival rate from body mass under certain circumstances (see 'Details').
#' @param Fs number ; safety factor, see 'Details'.
#' @param harvest.fixed number ; number of harvested individuals.
#' @param harvest.unif boolean ; default \code{FALSE}. If \code{TRUE}, simulates harvest level from a uniform distribution bounded by \code{min.harvest} and \code{max.harvest}.
#' @param min.harvest number ; minimum number of harvested individuals. To be used only if \code{harvest.unif = TRUE}.
#' @param max.harvest number ; maximum number of harvested individuals. To be used only if \code{harvest.unif = TRUE}.
#' @param harvest.lognorm boolean, default \code{FALSE}. If \code{TRUE}, simulates harvest level from a log-normal distribution with mu =  \code{mean.harvest} and standard deviation = \code{sd.harvest}.
#' @param mean.harvest number ; arithmetic mean of the number of harvested individuals. To be used only if \code{harvest.lognorm = TRUE}.
#' @param sd.harvest number ; standard deviation of the number of harvested individuals. To be used only if \code{harvest.lognorm = TRUE}.
#' @param pop.fixed number ; fixed estimate of population size (in number of individuals). To be used when \code{pop.unif = FALSE} and \code{pop.lognorm = FALSE}.
#' @param pop.unif boolean, default \code{FALSE}. If \code{TRUE}, simulates a population size from a uniform distribution bounded by \code{min.pop} and \code{max.pop}.
#' @param min.pop number ; minimum population size. To be used only if \code{pop.unif = TRUE}.
#' @param max.pop number ; maximum population size. To be used only if \code{pop.unif = TRUE}.
#' @param pop.lognorm boolean, default \code{FALSE}.  If \code{TRUE}, simulates a population size from a log-normal distribution with mu = \code{mean.pop} and standard deviation = \code{sd.pop}.
#' @param mean.pop number ; arithmetic mean of population size. To be used only if \code{pop.lognorm = TRUE}.
#' @param sd.pop number ; standard deviation of population size. To be used only if \code{pop.lognorm = TRUE}.
#'
#' @return A dataframe with the following information for each of the n simulations :
#' If the argument \code{full.option} is \code{FALSE}, it provides a simple dataframe :
#' \itemize{
#'   \item \code{SP} the identification number assigned to the studied species.
#'   \item \code{Sim} the index of the simulation.
#'   \item \code{Fs} the value of the safety factor.
#'   \item \code{Rmax} the maximal annual recruitment rate (lambda max – 1).
#'   \item \code{popsize} the population size (in number of individuals).
#'   \item \code{harvest}  the number of individuals harvested each year.
#'   \item \code{PEG} the calculated value of the Potential excess growth (in number of individuals).
#'   \item \code{SHI} the calculated value of the Sustainability harvest index (from 0 to +∞) (if > 1, harvest exceeds PEG).
#' }
#' If the argument \code{full.option} is \code{TRUE}, the following informations are added to the previous dataframe :
#' \itemize{
#'   \item \code{e} a value drawn from a normal distribution (0, 0.087). This parameter appears when the (annual) adult survival rate is estimated from body mass.
#'   \item \code{p} a value drawn from a beta distribution (3.34,101.24). This parameter appears when the annual survival rate of adults is estimated from body mass (\code{type.p = "random")}.
#' }
#'
#' @export
#' @seealso \code{\link{input.summary}}, \code{\link{output.summary}}.
#' @references
#' Johnson, F. A., Walters, M. A., & Boomer, G. S. (2012). Allowable levels of take for the trade in Nearctic songbirds. Ecological Applications, 22(4), 1114-1130.

#' Runge, M. C., Sauer, J. R., Avery, M. L., Blackwell, B. F., & Koneff, M. D. (2009). Assessing allowable take of migratory birds. The Journal of Wildlife Management, 73(4), 556-565.
#'
#' Boitani, L., & Fuller, T. (Eds.). (2000). Research techniques in animal ecology: controversies and consequences. Columbia university press.
#'
#'
#' @importFrom stats rlnorm rnorm runif
#'
#'
#' @examples
#' PEG(Nsim = 20, NSp = 2, living.rate = c("long", "short"),
#' surv.fixed = c(0.8, 0.65),
#' alpha.unif = TRUE, min.alpha = c(2, 1), max.alpha = c(3, 2),
#' pop.fixed = c(3605244, 55805898),
#' harvest.fixed = c(107802, 8447950),
#' Fs = c(0.1, 0.3, 0.5))
#'
PEG <- function(full.option = FALSE,
                        Nsim = 1,
                        NSp = NULL,
                        Rmax.fixed = NULL,
                        Rmax.lognorm = FALSE,
                        mean.Rmax = NULL,
                        sd.Rmax = NULL,
                        lambdaMax.fixed = NULL,
                        lambdaMax.lognorm = FALSE,
                        mean.lambdaMax = NULL,
                        sd.lambdaMax = NULL,
                        living.rate = NULL,
                        surv.fixed = NULL,
                        surv.beta = FALSE,
                        mean.surv = NULL,
                        sd.surv = NULL,
                        mass.fixed = NULL,
                        mass.lognorm = FALSE,
                        mean.mass = NULL,
                        sd.mass = NULL,
                        type.p = NULL,
                        type.e = NULL,
                        surv.j = NULL,
                        alpha.fixed = NULL,
                        alpha.unif = FALSE,
                        min.alpha = NULL,
                        max.alpha = NULL,
                        alpha.lognorm = FALSE,
                        mean.alpha = NULL,
                        sd.alpha = NULL,
                        Fs = NULL,
                        harvest.fixed = NULL,
                        harvest.unif = FALSE,
                        min.harvest = NULL,
                        max.harvest = NULL,
                        harvest.lognorm = FALSE,
                        mean.harvest = NULL,
                        sd.harvest = NULL,
                        pop.fixed = NULL,
                        pop.unif = FALSE,
                        min.pop = NULL,
                        max.pop = NULL,
                        pop.lognorm = FALSE,
                        mean.pop = NULL,
                        sd.pop = NULL
){
  
  ##PREAMBULE
  #Default value for surv.j to 1 if not provided by the user
 if (is.null(surv.j)){ surv.j<-rep(1, NSp)}
  
  ## I. Error messages about vector type ----
  # Print error messages if the vector types of the datas are not correct
  if(!is.null(Nsim) && !is.numeric(Nsim)){ # If 'Nsim' argument is not null and it is not a numeric vector then print an error message
    stop("'Nsim' must be a numeric vector")
  }
  if(!is.null(NSp) && !is.numeric(NSp)){
    stop("'NSp' must be a numeric vector")
  }
  if(!is.null(Fs) && !is.numeric(Fs)){
    stop("'Fs' must be a numeric vector")
  }

  ## Rmax argments ----
  if(!is.null(Rmax.fixed) & !is.numeric(Rmax.fixed)){
    stop("'Rmax.fixed' must be a numeric vector")
  }
  if(!is.null(mean.Rmax) & !is.numeric(mean.Rmax)){
    stop("'mean.Rmax' must be a numeric vector")
  }
  if(!is.null(sd.Rmax) & !is.numeric(sd.Rmax)){
    stop("'sd.Rmax' must be a numeric vector")
  }

  ## LambdaMax arguments ----
  if(!is.null(lambdaMax.fixed) & !is.numeric(lambdaMax.fixed)){
    stop("'lambdaMax.fixed' must be a numeric vector")
  }
  if(!is.null(mean.lambdaMax) & !is.numeric(mean.lambdaMax)){
    stop("'mean.lambdaMax' must be a numeric vector")
  }
  if(!is.null(sd.lambdaMax) & !is.numeric(sd.lambdaMax)){
    stop("'sd.lambdaMax' must be a numeric vector")
  }

  ## living rate arguments ----
  if(!is.null(living.rate) && !is.character(living.rate)){
    stop("'living.rate' must be a character vector")
  }

  ## survival arguments ----
  if(!is.null(surv.fixed) && !is.numeric(surv.fixed)){
    stop("'surv.fixed' must be a numeric vector")
  }
  if(!is.null(mean.surv) && !is.numeric(mean.surv)){
    stop("'mean.surv' must be a numeric vector")
  }
  if(!is.null(sd.surv) && !is.numeric(sd.surv)){
    stop("'sd.surv' must be a numeric vector")
  }

  ## alpha arguments ----
  if(!is.null(alpha.fixed) && !is.numeric(alpha.fixed)){
    stop("'alpha.fixed' must be a numeric vector")
  }
  if(!is.null(min.alpha) && !is.numeric(min.alpha)){
    stop("'min.alpha' must be a numeric vector")
  }
  if(!is.null(max.alpha) && !is.numeric(max.alpha)){
    stop("'max.alpha' must be a numeric vector")
  }
  if(!is.null(mean.alpha) && !is.numeric(mean.alpha)){
    stop("'mean.alpha' must be a numeric vector")
  }
  if(!is.null(sd.alpha) && !is.numeric(sd.alpha)){
    stop("'sd.alpha' must be a numeric vector")
  }

  ## mass arguments ----
  if(!is.null(mass.fixed) && !is.numeric(mass.fixed)){
    stop("'mass.fixed' must be a numeric vector")
  }
  if(!is.null(mean.mass) && !is.numeric(mean.mass)){
    stop("'mean.mass' must be a numeric vector")
  }
  if(!is.null(sd.mass) && !is.numeric(sd.mass)){
    stop("'sd.mass' must be a numeric vector")
  }
  if(!is.null(type.p) && !is.character(type.p)){
    stop("'type.p' must be a character vector")
  }
  if(!is.null(type.e) && !is.character(type.e)){
    stop("'type.e' must be a character vector")
  }

  ## population size arguments ----
  if(!is.null(pop.fixed) && !is.numeric(pop.fixed)){
    stop("'pop.fixed' must be a numeric vector")
  }
  if(!is.null(min.pop) && !is.numeric(min.pop)){
    stop("'min.pop' must be a numeric vector")
  }
  if(!is.null(max.pop) && !is.numeric(max.pop)){
    stop("'max.pop' must be a numeric vector")
  }
  if(!is.null(mean.pop) && !is.numeric(mean.pop)){
    stop("'mean.pop' must be a numeric vector")
  }
  if(!is.null(sd.pop) && !is.numeric(sd.pop)){
    stop("'sd.pop' must be a numeric vector")
  }

  ## harvest arguments -
  if(!is.null(harvest.fixed) && !is.numeric(harvest.fixed)){
    stop("'harvest.fixed' must be a numeric vector")
  }
  if(!is.null(min.harvest) && !is.numeric(min.harvest)){
    stop("'min.harvest' must be a numeric vector")
  }
  if(!is.null(max.harvest) && !is.numeric(max.harvest)){
    stop("'max.harvest' must be a numeric vector")
  }
  if(!is.null(mean.harvest) && !is.numeric(mean.harvest)){
    stop("'mean.harvest' must be a numeric vector")
  }
  if(!is.null(sd.harvest) && !is.numeric(sd.harvest)){
    stop("'sd.harvest' must be a numeric vector")
  }



  ## II. Error messages about input data ----
  # Print error messages if the input datas are not correct

  ## II.1. Arguments which can not be specified together ----
  ## Rmax arguments ----
  if(!is.null(Rmax.fixed) && (isTRUE(Rmax.lognorm)
                              || !is.null(mean.Rmax)
                              || !is.null(sd.Rmax))){
    stop("'Rmax.fixed' argument and 'Rmax' arguments for normal distribution can not be specified together")
  }
  if((!is.null(Rmax.fixed)
     || isTRUE(Rmax.lognorm)
     || !is.null(mean.Rmax)
     || !is.null(sd.Rmax)) && (!is.null(lambdaMax.fixed)
                              || isTRUE(lambdaMax.lognorm)
                              || !is.null(mean.lambdaMax)
                              || !is.null(sd.lambdaMax))){ #modified (line324-327) on 2022-07-20 : '()' added on the first set conditions
    stop("'Rmax' arguments and 'Lambdamax' arguments can not be specified together")
  }
  if((!is.null(Rmax.fixed)
     || isTRUE(Rmax.lognorm)
     || !is.null(mean.Rmax)
     || !is.null(sd.Rmax)) && !is.null(living.rate)){ #modified (line324-327) on 2022-07-20 : '()' added on the first set conditions
    stop("'Rmax' arguments and 'living.rate' arguments can not be specified together")
  }
  if((!is.null(Rmax.fixed)
      || isTRUE(Rmax.lognorm)
      || !is.null(mean.Rmax)
      || !is.null(sd.Rmax)) && (!is.null(surv.fixed)
                                || isTRUE(surv.beta)
                                || !is.null(mean.surv)
                                || !is.null(sd.surv))){
    stop("'Rmax' arguments and 'survival' arguments can not be specified together")
  }
  if((!is.null(Rmax.fixed)
      || isTRUE(Rmax.lognorm)
      || !is.null(mean.Rmax)
      || !is.null(sd.Rmax)) && (!is.null(alpha.fixed)
                                || isTRUE(alpha.lognorm)
                                || !is.null(mean.alpha)
                                || !is.null(sd.alpha)
                                || isTRUE(alpha.unif)
                                || !is.null(min.alpha)
                                || !is.null(max.alpha))){
    stop("'Rmax' arguments and 'alpha' arguments can not be specified together")
  }
  if((!is.null(Rmax.fixed)
      || isTRUE(Rmax.lognorm)
      || !is.null(mean.Rmax)
      || !is.null(sd.Rmax)) && (!is.null(mass.fixed)
                                || isTRUE(mass.lognorm)
                                || !is.null(mean.mass)
                                || !is.null(sd.mass)
                                || !is.null(type.p)
                                || !is.null(type.e))){
    stop("'Rmax' arguments and arguments to estimate survival based on body mass can not be specified together")
  }

  ## lambdaMax argments ----
  if(!is.null(lambdaMax.fixed) && (isTRUE(lambdaMax.lognorm)
                                   || !is.null(mean.lambdaMax)
                                   || !is.null(sd.lambdaMax))){
    stop("'lambdaMax.fixed' and 'lambdaMax' arguments for normal distribution can not be specified together")
  }
  if((!is.null(lambdaMax.fixed)
     || isTRUE(lambdaMax.lognorm)
     || !is.null(mean.lambdaMax)
     || !is.null(sd.lambdaMax)) && !is.null(living.rate)){
    stop("'lambdaMax' arguments and 'living.rate' arguments can not be specified together")
  }
  if((!is.null(lambdaMax.fixed)
      || isTRUE(lambdaMax.lognorm)
      || !is.null(mean.lambdaMax)
      || !is.null(sd.lambdaMax)) && (!is.null(surv.fixed)
                                     || isTRUE(surv.beta)
                                     || !is.null(mean.surv)
                                     || !is.null(sd.surv))){
    stop("'lambdaMax' arguments and 'survival' arguments can not be specified together")
  }
  if((!is.null(lambdaMax.fixed)
      || isTRUE(lambdaMax.lognorm)
      || !is.null(mean.lambdaMax)
      || !is.null(sd.lambdaMax)) && (!is.null(alpha.fixed)
                                     || isTRUE(alpha.lognorm)
                                     || !is.null(mean.alpha)
                                     || !is.null(sd.alpha)
                                     || isTRUE(alpha.unif)
                                     || !is.null(min.alpha)
                                     || !is.null(max.alpha))){
    stop("'lambdaMax' arguments and 'alpha' arguments can not be specified together")
  }
  if((!is.null(lambdaMax.fixed)
      || isTRUE(lambdaMax.lognorm)
      || !is.null(mean.lambdaMax)
      || !is.null(sd.lambdaMax)) && (!is.null(mass.fixed)
                                     || isTRUE(mass.lognorm)
                                     || !is.null(mean.mass)
                                     || !is.null(sd.mass)
                                     || !is.null(type.p)
                                     || !is.null(type.e))){
    stop("'lambdaMax' arguments and arguments to estimate survival based on body mass can not be specified together")
  }

  ## survival arguments ----
  if(!is.null(surv.fixed) && (isTRUE(surv.beta)
                              || !is.null(mean.surv)
                              || !is.null(sd.surv))){
    stop("'surv.fixed' and 'surv' arguments for normal distribution can not be specified together")
  }
  if((!is.null(surv.fixed)
      || isTRUE(surv.beta)
      || !is.null(mean.surv)
      || !is.null(sd.surv)) && (!is.null(mass.fixed)
                                || isTRUE(mass.lognorm)
                                || !is.null(mean.mass)
                                || !is.null(sd.mass)
                                || !is.null(type.p)
                                || !is.null(type.e))){
    stop("'surv.fixed' and arguments to estimate survival based on body mass can not be specified together")
  }

  if(length(surv.j)!=NSp){
    stop("the length of 'surv.j' argument must equal 'NSp'")
   }
  
  ## mass arguments
  if(!is.null(mass.fixed) && (isTRUE(mass.lognorm)
                              || !is.null(mean.mass)
                              || !is.null(sd.mass))){
    stop("'mass.fixed' and 'mass' arguments for normal distribution can not be specified together")
  }

  ## alpha arguments----
  if(!is.null(alpha.fixed) && (isTRUE(alpha.unif)
                               || !is.null(min.alpha)
                               || !is.null(max.alpha)
                               || isTRUE(alpha.lognorm)
                               || !is.null(mean.alpha)
                               || !is.null(sd.alpha))){
    stop("'alpha.fixed' and 'alpha' arguments for uniform or normal distribution can not be specified together")
  }
  if((isTRUE(alpha.unif)
      || !is.null(min.alpha)
      || !is.null(max.alpha)) && (isTRUE(alpha.lognorm)
                                  || !is.null(mean.alpha)
                                  || !is.null(sd.alpha))){
    stop("'alpha' arguments for uniform distribution and 'alpha' arguments for normal distribution
           can not be specified together")
  }

  ## popsize arguments----
  if(!is.null(pop.fixed) && (isTRUE(pop.unif)
                             || !is.null(min.pop)
                             || !is.null(max.pop)
                             || isTRUE(pop.lognorm)
                             || !is.null(mean.pop)
                             || !is.null(sd.pop))){
    stop("'pop.fixed' and 'pop' arguments for uniform and normal distribution can not be
           specified together")
  }
  if((isTRUE(pop.unif)
      || !is.null(min.pop)
      || !is.null(max.pop)) && (isTRUE(pop.lognorm)
                                || !is.null(mean.pop)
                                || !is.null(sd.pop))){
    stop("'pop' arguments for uniform distribution and 'pop' arguments for normal distribution
           can not be specified together")
  }

  ## harvest arguments----
  if(!is.null(harvest.fixed) && (isTRUE(harvest.lognorm)
                                 || !is.null(min.harvest)
                                 || !is.null(max.harvest)
                                 || isTRUE(harvest.lognorm)
                                 || !is.null(mean.harvest)
                                 || !is.null(sd.harvest))){
    stop("'harvest.fixed' and 'harvest' arguments for uniform and normal distribution can not be
           specified together")
  }
  if((isTRUE(harvest.unif)
      || !is.null(min.harvest)
      || !is.null(max.harvest)) && (isTRUE(harvest.lognorm)
                                    || !is.null(mean.harvest)
                                    || !is.null(sd.harvest))){
    stop("'harvest' arguments for uniform distribution and 'harvest' arguments for normal distribution
           can not be specified together")
  }


  



  ## II.2. Arguments which must be specified (together) ----
  if(is.null(NSp)){
    stop("'NSp' must be specified")
  }
  if(is.null(Fs)){
    stop("'Fs' must be specified")
  }

  ## Length of arguments
  var.storage <- c(length(living.rate), length(surv.fixed), length(mean.surv),
                   length(sd.surv), length(mass.fixed), length(mean.mass), length(sd.mass),
                   length(type.p), length(type.e),
                   length(alpha.fixed),length(min.alpha), length(max.alpha),
                   length(mean.alpha), length(sd.alpha), length(Rmax.fixed),
                   length(mean.Rmax), length(sd.Rmax), length(lambdaMax.fixed),
                   length(mean.lambdaMax), length(sd.lambdaMax), length(harvest.fixed),
                   length(min.harvest), length(max.harvest), length(mean.harvest),
                   length(sd.harvest), length(pop.fixed), length(min.pop),
                   length(max.pop), length(mean.pop), length(sd.pop))
  if(length(unique(var.storage[var.storage!=0])) > 1) {
    stop("arguments do not have the same length")
  }

  ## living.rate arguments----
  if (is.null(living.rate)
      && is.null(Rmax.fixed)
      && is.null(mean.Rmax)
      && is.null(sd.Rmax)
      && is.null(lambdaMax.fixed)
      && is.null(mean.lambdaMax)
      && is.null(sd.lambdaMax)){
    stop("'living.rate' or 'Rmax' arguments or 'lambdamax' arguments must be specified")
  }
  error.living.rate <- living.rate[which(living.rate != "short")]
  error.living.rate <- error.living.rate[which(error.living.rate != "long")]
  if(length(error.living.rate) >= 1 ){stop("'living.rate' must be 'short' or 'long'")}

  ## mass ----
  error.type.p <- type.p[which(type.p != "determinist")]
  error.type.p <- error.type.p[which(error.type.p != "random")]
  if(length(error.type.p) >= 1 ){stop("'type.p' must be 'determinist' or 'random'")}

  error.type.e <- type.e[which(type.e != "determinist")]
  error.type.e <- error.type.e[which(error.type.e != "random")]
  if(length(error.type.e) >= 1 ){stop("'type.e' must be 'determinist' or 'random'")}


  ## II.2.1 If input data are extract from a distribution ----
  ## Rmax arguments ----
  if(isTRUE(Rmax.lognorm) && (is.null(mean.Rmax) || is.null(sd.Rmax))){
    stop("if 'Rmax.lognorm' is 'TRUE', 'mean.Rmax' and 'sd.Rmax' must be specified")
  }

  ## lambdaMax arguments ----
  if(isTRUE(lambdaMax.lognorm) && (is.null(mean.lambdaMax) || is.null(sd.lambdaMax))){
    stop("if 'lambdaMax.lognorm' is 'TRUE', 'mean.lambdaMax' and 'sd.lambdaMax' must be specified")
  }

  ## survival arguments ----
  if(isTRUE(surv.beta) && (is.null(mean.surv) || is.null(sd.surv))){
    stop("if 'surv.beta' is 'TRUE', 'mean.surv' and 'sd.surv' must be specified")
  }

  ## alpha arguments ----
  if(isTRUE(alpha.lognorm) && (is.null(mean.alpha) || is.null(sd.alpha))){
    stop("if 'alpha.lognorm' is 'TRUE', 'mean.alpha' and 'sd.alpha' must be specified")
  }
  if(isTRUE(alpha.unif) && (is.null(min.alpha) || is.null(max.alpha))){
    stop("if 'alpha.unif' is 'TRUE', 'min.alpha' and 'max.alpha' must be specified")
  }

  ## mass arguments ----
  if(isTRUE(mass.lognorm) && (is.null(mean.mass) || is.null(sd.mass))){
    stop("if 'mass.lognorm' is 'TRUE', 'mean.mass' and 'sd.mass' must be specified")
  }
  if((!is.null(mass.fixed) || !is.null(mean.mass)) && is.null(type.p)){
    stop("if 'mass' arguments are specified, 'type.p' must be specified")
  }
  if((!is.null(mass.fixed) || !is.null(mean.mass)) && is.null(type.e)){
    stop("if 'mass' arguments are specified, 'type.e' must be specified")
  }

  ## popsize arguments ----
  if(isTRUE(pop.lognorm) && (is.null(mean.pop) || is.null(sd.pop))){
    stop("if 'pop.lognorm' is 'TRUE', 'mean.pop' and 'sd.pop' must be specified")
  }
  if(isTRUE(pop.unif) && (is.null(min.pop) || is.null(max.pop))){
    stop("if 'pop.unif' is 'TRUE', 'min.pop' and 'max.pop' must be specified")
  }

  ## harvest arguments ----
  if(isTRUE(harvest.lognorm) && (is.null(mean.harvest) || is.null(sd.harvest))){
    stop("if 'harvest.lognorm' is 'TRUE', 'mean.harvest' and 'sd.harvest' must be specified")
  }
  if(isTRUE(harvest.unif) && (is.null(min.harvest) || is.null(max.harvest))){
    stop("if 'harvest.unif' is 'TRUE', 'min.harvest' and 'max.harvest' must be specified")
  }




  ## II.2.2 If some arguments are not NULL ----
  ## Rmax arguments ----
  if((!is.null(mean.Rmax)||!is.null(sd.Rmax)) && !isTRUE(Rmax.lognorm)){
    stop("if 'mean.Rmax' or 'sd.Rmax' specified, 'Rmax.lognorm' must be 'TRUE'")
  }

  ## lambdaMax arguments ----
  if((!is.null(mean.lambdaMax)||!is.null(sd.lambdaMax)) && !isTRUE(lambdaMax.lognorm)){
    stop("if 'mean.lambdaMax' or 'sd.lambdaMax' specified, 'lambdaMax.lognorm' must be 'TRUE'")
  }

  ## survival arguments ----
  if((!is.null(mean.surv)||!is.null(sd.surv)) && !isTRUE(surv.beta)){
    stop("if 'mean.surv' or 'sd.surv' specified, 'surv.beta' must be 'TRUE'")
  }

  ## alpha arguments ----
  if((!is.null(mean.alpha)||!is.null(sd.alpha)) && !isTRUE(alpha.lognorm)){
    stop("if 'mean.alpha' or 'sd.alpha' specified, 'alpha.lognorm' must be 'TRUE'")
  }
  if((!is.null(min.alpha)||!is.null(max.alpha)) && !isTRUE(alpha.unif)){
    stop("if 'min.alpha' or 'max.alpha' specified, 'alpha.unif' must be 'TRUE'")
  }

  ## mass arguments ----
  if((!is.null(mean.mass)||!is.null(sd.mass)) && !isTRUE(mass.lognorm)){
    stop("if 'mean.mass' or 'sd.mass' specified, 'mass.lognorm' must be 'TRUE'")
  }

  ## popsize arguments ----
  if((!is.null(mean.pop)||!is.null(sd.pop)) && !isTRUE(pop.lognorm)){
    stop("if 'mean.pop' or 'sd.pop' specified, 'pop.lognorm' must be 'TRUE'")
  }
  if((!is.null(min.pop)||!is.null(max.pop)) && !isTRUE(pop.unif)){
    stop("if 'min.pop' or 'max.pop' specified, 'pop.unif' must be 'TRUE'")
  }

  ## harvest arguments ----
  if((!is.null(mean.harvest)||!is.null(sd.harvest)) && !isTRUE(harvest.lognorm)){
    stop("if 'mean.harvest' or 'sd.harvest' specified, 'harvest.lognorm' must be 'TRUE'")
  }
  if((!is.null(min.harvest)||!is.null(max.harvest)) && !isTRUE(harvest.unif)){
    stop("if 'min.harvest' or 'max.harvest' specified, 'harvest.unif' must be 'TRUE'")
  }

  ## II.3. Other error messages ----
  ## Rmax arguments ----
  if(isTRUE(mean.Rmax < sd.Rmax)){
    stop("'sd.Rmax' must be lower than 'mean.Rmax'")
  }
  if(isTRUE(Rmax.lognorm) && (isTRUE(lambdaMax.lognorm) || isTRUE(surv.beta) || isTRUE(alpha.lognorm) || isTRUE(alpha.unif) || isTRUE(mass.lognorm)) ){
    stop("'Rmax.lognorm' and 'lambdaMax.lognorm'or 'surv.beta'or 'alpha.unif' or alpha.lognorm' or 'mass.lognorm'  can not be TRUE together")
  }

  if(isTRUE(lambdaMax.lognorm) &&  (isTRUE(surv.beta) || isTRUE(alpha.lognorm) || isTRUE(alpha.unif) || isTRUE(mass.lognorm)) ){
    stop("'lambdaMax.lognorm' and 'surv.beta'or 'alpha.unif' or 'alpha.lognorm' or 'mass.lognorm'  can not be TRUE together")
  }

  if(isTRUE(surv.beta) && isTRUE(mass.lognorm) ){
    stop("'surv.beta' and 'mass.lognorm'  can not be TRUE together")
  }


  ##lambdaMax arguments ----
  if(isTRUE(mean.lambdaMax < sd.lambdaMax)){
    stop("'sd.lambdaMax' must be lower than 'mean.lambdaMax'")
  }

  ## survival arguments ----
  if(isTRUE(mean.surv < sd.surv)){
    stop("'sd.surv' must be lower than 'mean.surv'")
  }

  ## alpha arguments ----
  if(isTRUE(mean.alpha < sd.alpha)){
    stop("'sd.alpha' must be lower than 'mean.alpha'")
  }
  if(isTRUE(max.alpha < min.alpha)){
    stop("'min.alpha' must be lower than 'max.alpha'")
  }

  ## mass arguments
  if(isTRUE(mean.mass < sd.mass)){
    stop("'sd.surv' must be lower than 'mean.surv'")
  }

  ## popsize arguments----
  if(isTRUE(mean.pop < sd.pop)){
    stop("'sd.pop' must be lower than 'mean.pop'")
  }
  if(isTRUE(max.pop < min.pop)){
    stop("'min.pop' must be lower than 'max.pop'")
  }

  ## harvest arguments----
  if(isTRUE(mean.harvest < sd.harvest)){
    stop("'sd.harvest' must be lower than 'mean.harvest'")
  }
  if(isTRUE(max.harvest < min.harvest)){
    stop("'min.harvest' must be lower than 'max.harvest'")
  }


  ## III. Calculations ----
  ## Step 1 : Create the end storage vectors to save results of the calculations and load functions----
  Output1 <- NULL
  Output2 <- NULL

  # lognormal distribution (source : Fred Johnson)
  MOM.lognorm <- function(expect,stdev){
    sig <- log(1+(stdev^2/expect^2))
    se <- sqrt(sig)
    est <- log(expect)-(sig*0.5)
    result <- c(est,se)
    return(result)
  }



  ## Step 2 : Create the first loop calculations ----
  # This step allows to repeat the series of calculations whithin the loop as many times as there are studied species fixing temporary input values
  for (a in 1:NSp){ # Loop for calculations from 1 to the number of species studied

    # Create storage vectors to save results of the calculations for the loop a
    Rmax.calc <- NA

    ## Step 3 : Temporary setting of input values ----
    living.rateT <- living.rate[a] # Temporary setting the value of the object living.rate of the index a
    surv.fixedT <- ifelse(!is.null(surv.fixed), surv.fixed[a], NA) # If the argument 'surv.fixed' is not null, then temporary setting the value of 'surv.fixed' of the index a in the object named 'surv.fixedT', else assign 'NA' to the object 'surv.fixedT'. The same philosophy is applied bellow.
    mean.survT <- ifelse(!is.null(mean.surv), mean.surv[a], NA)
    sd.survT <- ifelse(!is.null(sd.surv), sd.surv[a], NA)
    mass.fixedT <- ifelse(!is.null(mass.fixed), mass.fixed[a], NA)
    mean.massT <- ifelse(!is.null(mean.mass), mean.mass[a], NA)
    sd.massT <- ifelse(!is.null(sd.mass), sd.mass[a], NA)
    type.pT <- ifelse(!is.null(type.p), type.p[a], NA)
    type.eT <- ifelse(!is.null(type.e), type.e[a], NA)
    surv.jT <- surv.j[a]
    alpha.fixedT <- ifelse(!is.null(alpha.fixed), alpha.fixed[a], NA)
    min.alphaT <- ifelse(!is.null(min.alpha), min.alpha[a], NA)
    max.alphaT <- ifelse(!is.null(max.alpha), max.alpha[a], NA)
    mean.alphaT <- ifelse(!is.null(mean.alpha), mean.alpha[a], NA)
    sd.alphaT <- ifelse(!is.null(sd.alpha), sd.alpha[a], NA)
    Rmax.fixedT <- ifelse(!is.null(Rmax.fixed), Rmax.fixed[a], NA)
    mean.RmaxT <- ifelse(!is.null(mean.Rmax), mean.Rmax[a], NA)
    sd.RmaxT <- ifelse(!is.null(sd.Rmax), sd.Rmax[a], NA)
    lambdaMax.fixedT <- ifelse(!is.null(lambdaMax.fixed), lambdaMax.fixed[a], NA)
    mean.lambdaMaxT <- ifelse(!is.null(mean.lambdaMax), mean.lambdaMax[a], NA)
    sd.lambdaMaxT <- ifelse(!is.null(sd.lambdaMax), sd.lambdaMax[a], NA)
    harvest.fixedT <- ifelse(!is.null(harvest.fixed), harvest.fixed[a], NA)
    min.harvestT <- ifelse(!is.null(min.harvest), min.harvest[a], NA)
    max.harvestT <- ifelse(!is.null(max.harvest), max.harvest[a], NA)
    mean.harvestT <- ifelse(!is.null(mean.harvest), mean.harvest[a], NA)
    sd.harvestT <- ifelse(!is.null(sd.harvest), sd.harvest[a], NA)
    pop.fixedT <- ifelse(!is.null(pop.fixed), pop.fixed[a], NA)
    min.popT <- ifelse(!is.null(min.pop), min.pop[a], NA)
    max.popT <- ifelse(!is.null(max.pop), max.pop[a], NA)
    mean.popT <- ifelse(!is.null(mean.pop), mean.pop[a], NA)
    sd.popT <- ifelse(!is.null(sd.pop), sd.pop[a], NA)



    ## Step 4 : Create the second loop calculations in the first one ----
    # This step allows to repeat the series of calculations whithin the loop as many times as there are safety factors
    for (b in 1:length(Fs)){ # Loop for calculations according to the safety factor
      FsT <- Fs[b]



      ## Step 5 : Estimate Rmax ----
      ## Condition "if" : if Rmax and lambdaMax are unknown ----
      # In this case, Rmax must be calculated with the rmax function
      if(is.null(Rmax.fixed) && is.null(mean.Rmax) && is.null(lambdaMax.fixed) && is.null(mean.lambdaMax)){ # If the arguments 'Rmax, 'mean.Rmax, 'lambdaMax', and "mean.lambdaMax" are null then :
        ## Rmax is calculated and datas are save in a temporary data frame named 'OutputT1'
        OutputT1 <- rmax(Nsim = Nsim,
                         living.rate = living.rateT,
                         surv.fixed = surv.fixedT,
                         surv.beta = surv.beta,
                         mean.surv = mean.survT,
                         sd.surv = sd.survT,
                         mass.fixed = mass.fixedT,
                         mass.lognorm = mass.lognorm,
                         mean.mass = mean.massT, #modified the 2022-07-20 : changed mean.massT to sd.massT ## back to mean.massT the 2024-12-11
                         sd.mass = sd.massT,
                         type.p = type.pT,
                         type.e = type.eT,
                         surv.j = surv.jT,
                         alpha.fixed = alpha.fixedT,
                         alpha.unif = alpha.unif,
                         min.alpha = min.alphaT,
                         max.alpha = max.alphaT,
                         alpha.lognorm = alpha.lognorm,
                         mean.alpha = mean.alphaT,
                         sd.alpha = sd.alphaT)



      ## Condition "else" : if Rmax or lambdaMax are known ----
      # In this case, Rmax values are generated without calculations
      }else{
        if(isTRUE(!is.na(Rmax.fixedT))){ # If the temporary value of Rmax.fixedT is different to 'NA' (i.e. the value of Rmax is known) then :
          Rmax.calc <- rep(Rmax.fixedT, Nsim) # The value of Rmax.fixedT is replicated 'Nsim' times and saved in the object 'Rmax'
        }
        if(isTRUE(!is.na(lambdaMax.fixedT))){ # If the temporary value of lambdaMax.fixedT is different to 'NA' (i.e. the value of lambdaMax is known) then :
          Rmax.calc <- lambdaMax.fixedT - 1 # Rmax is the difference between 'lambdaMax.fixed' and 1
        }
        if(isTRUE(Rmax.lognorm)){ # If Rmax is extrated from a normal distribution then :
          Rmax.param <- MOM.lognorm(mean.RmaxT, sd.RmaxT)
          Rmax.calc <- rlnorm(Nsim, Rmax.param[1], Rmax.param[2]) # 'Nsim' values are extracted from a log-normal distribution with with the mean 'mean.Rmax' and the standard deviation 'sd.Rmax'
        }
        if(isTRUE(lambdaMax.lognorm)){ # If lambdaMax is extracted from a normal distribution then :
          lambdaMax.param <- MOM.lognorm(mean.lambdaMaxT, sd.lambdaMaxT)
          lambdaMax <- rlnorm(Nsim, lambdaMax.param[1], lambdaMax.param[2]) # 'Nsim' values are extracted from a log-normal distribution with with the mean 'mean.lambdaMax' and the standard deviation 'sd.lambdaMax'
          Rmax.calc <- lambdaMax - 1 # Rmax is the difference between 'lambdaMax' and 1
        }
      } # End of condition "if else"



      ## Step 6 : Generate estimated variables of 'harvest' and 'popsize'----
      ## Generate the estimated variables of 'harvest' ----
      # This step is necessary to compare the harvest with the PTL
      if(isTRUE(harvest.unif)){ # If harvest is extrated from a uniform distribution then :
        harvest.calc <- runif(Nsim,min.harvestT,max.harvestT) # 'Nsim' values are extracted from a uniform distribution with 'min.harvestT' as lower bound and 'max.harvestT' as upper bound
        harvest.calc <- round(harvest.calc) # Round the values of 'harvest' object to the whole
      }
      if(isTRUE(harvest.lognorm)){ # If harvest is extracted from a normal distribution then :
        harvest.param <-MOM.lognorm(mean.harvestT, sd.harvestT)
        harvest.calc <- rlnorm(Nsim, harvest.param[1], harvest.param[2]) # 'Nsim' values are extracted from a log-normal distribution with with the mean 'mean.harvest' and the standard deviation 'sd.harvest'
        harvest.calc <- round(harvest.calc) # Round the values of 'harvest' object to the whole
      }
      if(isTRUE(!is.na(harvest.fixedT))){ # If the temporary value of harvest.fixedT is different to 'NA' (i.e. the value of harvest is known) then :
        harvest.calc <- rep(harvest.fixedT, Nsim) # The value of harvest.fixedT is replicated 'Nsim' times and saved in the object 'harvest'
      }



      ## Generate the estimated variable of 'popsize' ----
      # This step is necessary to calculate the PTL
      if(isTRUE(pop.unif)){ # If popsize is extrated from a uniform distribution then :
        popsize.calc <- runif(Nsim,min.popT,max.popT) # 'Nsim' values are extracted from a uniform distribution with 'min.popT' as lower bound and 'max.popT' as upper bound
        popsize.calc <- round(popsize.calc) # Round the values of 'popsize' object to the whole
      }
      if(isTRUE(pop.lognorm)){  # If popsize is extracted from a normal distribution then :
        pop.param <-MOM.lognorm(mean.popT, sd.popT)
        popsize.calc <- rlnorm(Nsim, pop.param[1], pop.param[2]) # 'Nsim' values are extracted from a log-normal distribution with with the mean 'mean.pop' and the standard deviation 'sd.pop'
        popsize.calc <- round(popsize.calc) # Round the values of 'popsize' object to the whole
      }
      if(isTRUE(!is.na(pop.fixedT))){ # If the temporary value of pop.fixedT is different to 'NA' (i.e. the value of popsize is known) then :
        popsize.calc <- rep(pop.fixedT, Nsim) # The value of pop.fixedT is replicated 'Nsim' times and saved in the object 'popsize'
      }



      ## Save all datas generated in the step 5 in a data frame Output 1
      if(exists("OutputT1")){ # If the data frame 'OutputT1' exist (i.e. Rmax was calculated from the rmax function) then :
        Output1 <- rbind(Output1, OutputT1) # Combine by rows all temporary data frame 'OutputT1' for each index a of the loop and save in the data frame 'Output1'
      }

      ## Step 7 : Save datas generated in the step 6 in a temporary dataframe named 'OutputT2'----
      NsimT <- 1:Nsim # This object allows to known the value of the index a for each row of the data frame below
      OutputT2 <- data.frame(SP = a, Nsim = NsimT, Fs = FsT, popsize = popsize.calc, harvest = harvest.calc) # Create the OutputT2

      ## Step 8 : Complementary columns for the 'OutputT2' ----
      # This step integrates columns in the final data frame according to the arguments informed
      var.storage <- c("Rmax.calc", "mean.RmaxT", "sd.RmaxT", "lambdaMax.fixedT", "mean.lambdaMaxT", "sd.lambdaMaxT",
                       "min.harvestT", "max.harvestT", "mean.harvestT", "sd.harvestT",
                       "min.popT", "max.popT", "mean.popT", "sd.popT") # Save in an object all sumpplementary possible columns
      not.na <- which(!is.na(mget(var.storage))) # Select columns which values are different from NA. Works only with NA. That is why parameters which are NULL are transformed as NA
      OutputT2[,var.storage[not.na]]<- mget(var.storage)[not.na] # Insert values columns which values are different from NA in the temporary data frame 'OutputT2'

      ## Step 9 : Reorder columns of the 'OutputT2' ----
      var.rmax<- c("Rmax.calc", "mean.RmaxT", "sd.RmaxT", "lambdaMax.fixedT", "mean.lambdaMaxT", "sd.lambdaMaxT") # Order of the columns wanted in the final output
      var.rmax.name <- var.rmax[var.rmax %in% colnames(OutputT2)] # Reorder the columns in the data frame named 'OutputT2'

      var.pop <- c("min.popT", "max.popT", "mean.popT", "sd.popT")
      var.pop.name <- var.pop[var.pop %in% colnames(OutputT2)]

      var.harv <- c("min.harvestT", "max.harvestT", "mean.harvestT", "sd.harvestT")
      var.harv.name <- var.harv[var.harv %in% colnames(OutputT2)]

      reorder <- c("SP", "Nsim", "Fs", var.rmax.name, "popsize", var.pop.name, "harvest", var.harv.name) # Final order of the columns wanted in the final output
      OutputT2 <- data.frame(OutputT2[, reorder]) # Reorder the columns in the data frame named 'OutputT2'

      # Step 10 : Save all the data from the 'OutputT2' in the data frame named 'Output2' ----
      Output2 <- rbind(Output2, OutputT2) # Combin by rows all temporary data frame 'OutputT2' for each index a of the loop and save in the data frame 'Output2'

    } # End of loop b (according to the safety factor)

  } # End of loop a (according to the number of species studied)



  ## Step 11 : Save all the results in a unique data frame named 'full.tab' ----
  if(exists("OutputT1")){ # If the 'OutputT1' exixts then :
    Output1 <- Output1[, !colnames(Output1) == "Nsim"] # Remove 'Nsim' columns in the 'Output1' because it already exists in the Output2
    full.tab <- cbind(Output2[,1:3], Output1, Output2[4:ncol(Output2)] ) # Combine by columns the data frames 'Output1' and 'Output2' and save in the data frame 'full.tab'
  }else{ # If the OutputT1' does not exixts then :
    full.tab <- Output2 # The data frame 'Output2' become the data frame 'full.tab'
  }



  ## Step 12 : Rename colums ----
  # The two lines of code below code removes the 'T' at the end of the columns when it exists
  cond<-which(substr(colnames(full.tab), nchar(colnames(full.tab)),  nchar(colnames(full.tab)))=="T")
  colnames(full.tab)[cond]<-substr(colnames(full.tab)[cond], 1, nchar(colnames(full.tab)[cond])-1)
  colnames(full.tab)[colnames(full.tab)=="Nsim"] <- "Sim" # Rename 'Nsim' as 'Sim'
  colnames(full.tab)[colnames(full.tab)=="Rmax.calc"] <- "Rmax"


  ## Step 13 : Final calculations ----
  # Calculations of PEG ----
  PEG.rate <- full.tab$Rmax * full.tab$Fs
  full.tab$PEG <- PEG.rate * full.tab$popsize

  # Calculate SHI ----
  full.tab$SHI <- full.tab$harvest/full.tab$PEG



  ## Step 14 : Format of the final data frame named 'final.tab'
  if(!isTRUE(full.option)){ # If the argument 'full.option' is FALSE then the final data frame named 'final.tab' is in the form with the minimum of necessary columns
    final.tab <- data.frame(SP = full.tab$SP, Sim = full.tab$Sim, Fs = full.tab$Fs,
                            Rmax = full.tab$Rmax, popsize = full.tab$popsize, harvest = full.tab$harvest,
                            PEG = full.tab$PEG, SHI = full.tab$SHI)
    }else{final.tab <- full.tab} # If the argument 'full.option' is TRUE then the final data frame named 'final.tab' is the form with all values of parameters used for the calculations



  ## Step 15 :  Remove of the temporary vectors
  suppressWarnings(rm(OutputT1, OutputT2, living.rateT, surv.fixedT, mean.survT, sd.survT, mass.fixedT, mean.massT,
     sd.massT, type.pT, type.eT, alpha.fixedT, min.alphaT, max.alphaT, mean.alphaT, sd.alphaT, Rmax.fixedT, mean.RmaxT, sd.RmaxT,
     lambdaMax.fixedT, mean.lambdaMaxT, sd.lambdaMaxT, harvest.fixedT, min.harvestT, max.harvestT, mean.harvestT, sd.harvestT,
     pop.fixedT, min.popT, max.popT, mean.popT, sd.popT, FsT))

  return(final.tab)


  } # End of function
