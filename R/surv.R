#' Survival estimation
#'
#' \code{surv} is a function used to estimate (maximal) annual survival rate of adults under ideal conditions, from body mass (in kilograms), see 'Details'.
#'
#' @details
#' Adult survival rate can be estimated from body mass following the relationship described in Johnson \emph{et al.} (2012) : \eqn{(p/surv.j)^(1/(exp(3.22+0.24*log(body mass)+e)-alpha))}.
#' IBody mass can be provided as either a point estimate (\code{mass.fixed}) or drawn from
#' a log-normal distribution (\code{mass.lognorm = TRUE}, \code{mean.mass} and \code{sd.mass} arguments). The parameter p can be either drawn from a beta distribution
#' B(3.34, 101.24) (\code{type.p = "random"}) or provided as a point estimate equal to 0.032 (i.e. the mean of the beta distribution is given by a/(a + b) with a = 3.34 and b = 101.24, Boitani & Fuller 2000) (\code{type.p = "determinist"})
#' The parameter e (residuals) can be either set to zero (\code{type.e = "determinist"}) or residuals drawn of a normal distribution N(0, 0.087) (\code{type.e = "random"}). See Johnson \emph{et al.} (2012).
#'
#'
#' By default, the argument \code{surv.j} is silent and masked to users (the default value for \code{surv.j} is 1, as such \eqn{p/surv.j = p/1 = p} ; see original equation in
#' Johnson \emph{et al.}, 2012). This is correct for any alpha if, and only if, juvenile survival < adult survival solely for birds aged <1 year.  This is a reasonable
#' approximation especially in passerines and some other medium-sized bird species. For species with delayed sexual maturity and when juvenile survival < adult
#' survival for birds aged 1 to alpha, the default \code{surv} function returns biased estimates unless the users specify a value for annual survival for birds aged 1
#' to alpha. This is done by using the argument \code{surv.j} which is then used to estimate  \eqn{p = p/surv.j}.
#'
#' @param Nsim number ; number of simulations.
#' @param NSp number ; number of species for which calculations are performed.
#' @param mass.fixed number ; point estimate of body mass (in kilograms) used to estimate annual survival rate of adults.
#' @param mass.lognorm boolean, default \code{FALSE}. If \code{TRUE}, simulates body mass from a log-normal distribution with mu = \code{mean.mass} and standard deviation = \code{sd.mass}.
#' @param mean.mass number ; arithmetic mean of body mass (in kilograms). Only used if \code{mass.lognorm = TRUE}.
#' @param sd.mass number ; standard deviation of body mass. Only used if \code{mass.lognorm = TRUE}.
#' @param alpha number ; age at the first reproduction
#' @param type.p character string ; the value of the parameter p can be either "\code{determinist}" or "\code{random}". No other character string allowed.
#' @param type.e character string ; the value of the parameter e can be either "\code{determinist}" or "\code{random}". No other character string allowed.
#' @param surv.j number, default 1 ; point estimate average annual survival for birds aged 1 to alpha. The value can be provided by the user to estimate adult survival rate from body mass under certain circumstances (see 'Details).
#'
#' @return A dataframe with the following information for each of the n simulations :
#' \itemize{
#'   \item \code{Sim} number of the simulation.
#'   \item \code{mass} body mass (in kilograms).
#'   \item \code{e} residuals
#'   \item \code{p}
#'   \item \code{survival} estimate of annual survival rate of adult under ideal conditions.
#' }
#'
#'
#' @export
#' @seealso \code{\link{PEG}}, \code{\link{PTL}}
#' @references
#' Johnson, F. A., Walters, M. A., & Boomer, G. S. (2012). Allowable levels of take for the trade in Nearctic songbirds. Ecological Applications, 22(4), 1114-1130.
#'
#' Boitani, L., & Fuller, T. (Eds.). (2000). Research techniques in animal ecology: controversies and consequences. Columbia university press.
#'
#'
#' @importFrom stats rlnorm rnorm rbeta
#'
surv <- function (Nsim = 1,
                  NSp = 1,
                  mass.fixed = NULL,
                  mass.lognorm = FALSE,
                  mean.mass = NULL, sd.mass = NULL,
                  alpha = NULL,
                  type.p = NULL,
                  type.e = NULL,
                  surv.j = 1) {

  ## Step 1 : Create storage vectors to save the results of the calculations and load functions ----
  Output <- NULL

  # lognormal distribution (source : Fred Johnson)
  MOM.lognorm <- function(expect,stdev){
    sig <- log(1+(stdev^2/expect^2))
    se <- sqrt(sig)
    est <- log(expect)-(sig*0.5)
    result <- c(est,se)
    return(result)
  }



  ## Step 2 : Estimate of adult annual survival rate ----
  for (i in 1:NSp){ # Loop to estimate adult annual survival rate i times

    # Temporary setting of the i-th value of input values
    alphaT <- rep(alpha[i], Nsim)
    type.pT <- type.p[i]
    type.eT <- type.e[i]
    surv.jT <- surv.j[i]

    # Estimate of body mass
    if(isTRUE(mass.lognorm)){ # If mass.lognorm argument is TRUE
      mean.massT <- mean.mass[i] # then temporarly fix the i-th value of the arguments mean.mass in the arguments named mean.massT
      sd.massT <- sd.mass[i] # and temporarly fix the i-th value of the arguments sd.mass in the arguments named sd.massT
      mass.param <- MOM.lognorm(mean.mass[i], sd.mass[i])
      massT <- rlnorm(Nsim, mass.param[1], mass.param[2]) # 'Nsim' values are extracted from a log-normal distribution with with the mean 'mean.mass' and the standard deviation 'sd.mass'
    }else{massT <- rep(mass.fixed[i], Nsim)} # else temporarly fix the i-th value of the argument mass.fixed in the argument named massT


    # Estimate adult annual survival rate and save values of the parameter necessary to this estimate
      if(type.pT == "determinist"){
        p <- rep(3.34/(3.34+101.24),Nsim) # The value p is equal to the formula opposite
      } else {
        p<-rbeta(Nsim,3.34,101.24) # Draw the value of p according to the beta distribution (3.34,101.24)
      }
      if (type.eT == "determinist"){
        e <- rep(0, Nsim)
      } else {
        e<-rnorm(n=Nsim,mean=0, sd=sqrt(0.087)) # Random draw of the residue of the relationship between survival and body mass (see Johnson et al., 2012)
      }

      survival <- (p/surv.jT)^(1/(exp(3.22+0.24*log(massT)+e)-alphaT)) # Survival calculation according body mass


    ## Step 3 : Final output ----
    Sim <- 1:Nsim
    OutputT <- data.frame(Sim = Sim, mass = massT, alpha = alphaT, # Save temporarly values used to estimate adult annual survival rate in a dataframe name OutputT
                          p = p, e = e, surv.j = surv.j, survival = survival)
    if(!is.null(mean.mass)){OutputT$mean.mass <- mean.massT} # If mean.mass argument is not null then save the value from mean.massT argument in the dataframe named OutputT
    if(!is.null(sd.mass)){OutputT$sd.mass <- sd.massT} # If sd.mass argument is not null then save the value from sd.massT argument in the dataframe named OutputT

    Output <- rbind(Output, OutputT) # Combine by rows all temporary data frames

  } # End of loop i

  return(Output) # return the object Output in the console

} # End of function
