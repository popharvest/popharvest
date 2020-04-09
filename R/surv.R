#' Survival estimation
#'
#' For internal use only. Survival estimation based on body mass.
#'
#' @param Nsim number ; number of simulations.
#' @param NSp number ; number of species for which calculations are done.
#' @param mass.fixed number ; fixed estimate of body mass (in kilograms) used to estimate adult survival rate from body mass.
#' @param mass.norm boolean, default \code{FALSE}. If \code{TRUE}, simulates body mass from a normal distribution with mu = \code{mean.mass} and standard deviation = \code{sd.mass}.
#' @param mean.mass number ; arithmetic mean of body mass (in kilograms). Only used if \code{mass.norm = TRUE}.
#' @param sd.mass number ; standard deviation of body mass. Only used if \code{mass.norm = TRUE}.
#' @param type.p character string ; the value of the parameter p can be either "\code{determinist}" or "\code{random}". No other character string allowed.
#' @param type.e character string ; the value of the parameter e can be either "\code{determinist}" or "\code{random}". No other character string allowed.
#'
#' @return A dataframe with the following information for each of the n simulations :
#' \itemize{
#'   \item \code{Sim} number of the simulation.
#'   \item \code{mass} body mass (in kilograms).
#'   \item \code{e}
#'   \item \code{p}
#'   \item \code{survival} estimate of adult annual survival rate.
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
#' @importFrom stats rnorm rbeta
#'
surv <- function (Nsim = 1,
                  NSp = 1,
                  mass.fixed = NULL,
                  mass.norm = FALSE,
                  mean.mass = NULL, sd.mass = NULL,
                  type.p = NULL,
                  type.e = NULL) {

  ## Step 1 : Creating storage vectors to save the results of the calculations ----
  Output <- NULL
  e.storage <- NULL
  p.storage <- NULL
  survival.storage <- NULL
  Nsim.storage <- NULL



  ## Step 2 : Estimate of adult annual survival rate ----
  for (i in 1:NSp){ # Loop to estimate adult annual survival rate i times

    # Temporary setting of the i-th value of input values
    type.pT <- type.p[i]
    type.eT <- type.e[i]

    # Estimate of body mass
    if(isTRUE(mass.norm)){ # If mass.norm argument is TRUE
      mean.massT <- mean.mass[i] # then temporarly fix the i-th value of the arguments mean.mass in the arguments named mean.massT
      sd.massT <- sd.mass[i] # and temporarly fix the i-th value of the arguments sd.mass in the arguments named sd.massT
      massT <- rnorm(1, mean.massT, sd.massT) # 1 value is extracted from a normal distribution with the mean 'mean.massT' and the standard deviation 'sd.massT'
    }else{massT <- mass.fixed[i]} # else temporarly fix the i-th value of the argument mass.fixed in the argument named massT


    for (j in 1:Nsim){ # Loop to estimate adult annual survival rate and save values of the parameter necessary to this estimate

      if(type.pT == "determinist"){
        p <- 3.34/(3.34+101.24) # The value p is equal to the formula opposite
      } else {
        p<-rbeta(1,3.34,101.24) # Draw the value of p according to the beta distribution (3.34,101.24)
      }
      if (type.eT == "determinist"){
        e <- 0
      } else {
        e<-rnorm(n=1,mean=0, sd=sqrt(0.087)) # Random draw of the residue of the relationship between survival and body mass (cf. Johnson et al. 2012)
      }

      survival <- p^(1/(exp(3.22+0.24*log(massT)+e)-1)) # Survival calculation according body mass

      e.storage[j] <- e # saving the j-th value of e
      p.storage[j] <- p # saving the j-th value of p
      survival.storage[j] <-survival # saving the j-th value of survival
      Nsim.storage[j] <- j # saving the value of the object j
    } # End of loop j

    ## Step 3 : Final output ----
    OutputT <- data.frame(Sim = Nsim.storage, mass = massT, # Saving temporarly values used to estimate adult annual survival rate in a dataframe name OutputT
                          p = p.storage, e = e.storage, survival = survival.storage)
    if(!is.null(mean.mass)){OutputT$mean.mass <- mean.massT} # If mean.mass argument is not null then save the value from mean.massT argument in the dataframe named OutputT
    if(!is.null(sd.mass)){OutputT$sd.mass <- sd.massT} # If sd.mass argument is not null then save the value from sd.massT argument in the dataframe named OutputT

    Output <- rbind(Output, OutputT) # Combining by rows all temporary data frames

  } # End of loop i

  return(Output) # return the object Output in the console

} # End of function
