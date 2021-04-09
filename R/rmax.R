#' Maximum annual recruitment rate estimation
#'
#' \code{rmax} is a function used to estimate the maximum annual recruitment rate, see 'Details'.
#'
#' @details
#' When the maximum annual recruitment rate, Rmax, is unknown,
#' the value can be estimated using the Demographic Invariants Method (DIM, Niel & Lebreton, 2005). This approach calculates \eqn{Rmax = lambdaMax – 1} (Dillingham & Fletcher, 2008) on the basis of 2 parameters:
#' annual survival rate of adult under ideal condition (maximal survival, \code{surv} arguments) and age at first breeding \code{alpha} arguments). When adult survival and age at first breeding are known,
#' users can be provide corresponding values either as point estimates (\code{surv.fixed}, \code{alpha.fixed} arguments) or as distributions (\code{surv.beta = TRUE} with \code{mean.surv} and \code{sd.surv},
#' \code{alpha.unif = TRUE} with \code{min.alpha} and \code{max.alpha}, \code{alpha.lognorm = TRUE} with \code{mean.alpha} and \code{sd.alpha} arguments). Alternatively,
#' adult survival rate can be estimated from body mass following the relationship described in Johnson \emph{et al.} (2012) : \eqn{(p/surv.j)^(1/(exp(3.22+0.24*log(body mass)+e)-alpha))} (see \code{surv} function).
#' Body mass can be provide either as a point estimate (\code{mass.fixed})
#' or drawn from a log-normal distribution (\code{mass.lognorm = TRUE}, \code{mean.mass} and \code{sd.mass} arguments). The parameter p can be drawn either from a beta distribution
#' B(3.34, 101.24) (\code{type.p = "random"}) or a point estimate equal to 0.032 (i.e. the mean of the beta distribution is given by a/(a + b) with a = 3.34 and b = 101.24, Boitani & Fuller 2000) (\code{type.p = "determinist"})
#' The parameter e (residuals) can be either set to zero (\code{type.e = "determinist"}) or residuals from a normal distribution N(0, 0.087) (\code{type.e = "random"}). See equations in Johnson \emph{et al.} (2012).
#'
#'
#' By default, the argument \code{surv.j} is silent and masked to users (the default value for \code{surv.j} is 1, as such \eqn{p/surv.j = p/1 = p} ; see original equation in
#' Johnson \emph{et al.}, 2012). This is correct for any alpha if, and only if, juvenile survival < adult survival solely for birds aged <1 year.  This is a
#' reasonable approximation especially in passerines and some other medium-sized bird species. For species with delayed sexual maturity and when juvenile
#' survival < adult survival for birds aged 1 to alpha, the default \code{surv} function returns biased estimates unless the users specify a value for annual survival
#' for birds aged 1 to alpha. This is done by using the argument \code{surv.j} which is then used to estimate \eqn{p = p/surv.j}.
#'
#' @param Nsim number, default 1 ;  number of simulations
#' @param living.rate character string ; living rate : either "\code{long}" or "\code{short}". No other character string allowed.
#' @param surv.fixed number ; fixed estimate of (annual) adult survival rate provided by the user. To be used instead of \code{mean.surv} when \code{surv.beta = FALSE}.
#' @param surv.beta surv.beta boolean, default FALSE. If \code{TRUE}, simulates (annual) adult survival rate from a beta distribution with mu = \code{mean.surv} and standard deviation = \code{sd.surv} (see below).
#' @param mean.surv mean.surv number ; arithmetic mean of (annual) adult survival rate. Only used if \code{surv.beta = TRUE}.
#' @param sd.surv sd.surv number ; standard deviation of adult annual survival rate. Only used if \code{surv.beta = TRUE}.
#' @param alpha.fixed alpha.fixed number ; fixed estimate of age at first breeding (in years). To be used when \code{alpha.unif = FALSE} and \code{alpha.lognorm = FALSE}.
#' @param alpha.unif  alpha.unif boolean, default \code{FALSE}. Simulates age at first breeding from a uniform distribution bounded by \code{min.alpha} and \code{max.alpha}.
#' @param min.alpha  min.alpha number ; minimum age at first breeding. To be used only if \code{alpha.unif = TRUE}.
#' @param max.alpha  max.alpha number ; maximum age at first breeding. To be used only if \code{alpha.unif = TRUE}.
#' @param alpha.lognorm  alpha.lognorm boolean, default \code{FALSE}. If \code{TRUE}, simulates age at first breeding from a normal distribution with mu = \code{mean.alpha} and standard deviation = \code{sd.alpha}.
#' @param mean.alpha mean.alpha number ; arithmetic mean of age at first breeding. To be used only used if \code{alpha.lognorm = TRUE}.
#' @param sd.alpha sd.alpha number ; age at first breeding standard deviation. Only used if \code{alpha.lognorm = TRUE}.
#' @param mass.fixed mass.fixed number ; fixed estimate of body mass (in kilograms) used to estimate adult survival rate from body mass.
#' @param mass.lognorm  mass.lognorm boolean, default \code{FALSE}. If \code{TRUE}, simulates body mass from a normal distribution with mu = \code{mean.mass} and standard deviation = \code{sd.mass}.
#' @param mean.mass  mean.mass number ; arithmetic mean of body mass (in kilograms). Only used if \code{mass.lognorm = TRUE}. #' @param sd.mass number ; standard deviation of body mass. Only used if \code{mass.lognorm = TRUE}.
#' @param sd.mass sd.mass number ; standard deviation of body mass. Only used if \code{mass.lognorm = TRUE}.
#' @param type.p type.p character string ; the value of the parameter p can be either "\code{determinist}" or "\code{random}". No other character string allowed.
#' @param type.e type.e character string ; the value of the parameter e can be either "\code{determinist}" or "\code{random}". No other character string allowed.
#' @param surv.j surv.j number, default 1 ; point estimate of average annual survival for birds aged 1 to alpha. The value can be provided by the user. Used to estimate adult survival rate from body mass under certain circumstances (see Details).

#'
#' @return A dataframe
#' @export
#' @seealso \code{\link{PEG}}, \code{\link{PTL}}.
#' @references
#' Dillingham, P. W., & Fletcher, D. (2008). Estimating the ability of birds to sustain additional human-caused mortalities using a simple decision rule and allometric relationships. Biological Conservation, 141(7), 1783-1792.
#'
#' Johnson, F. A., Walters, M. A., & Boomer, G. S. (2012). Allowable levels of take for the trade in Nearctic songbirds. Ecological Applications, 22(4), 1114-1130.
#'
#' Niel, C., & LEBRETON, J. D. (2005). Using demographic invariants to detect overharvested bird populations from incomplete data. Conservation Biology, 19(3), 826-835.
#'
#' @importFrom stats rlnorm rnorm rbeta
#'
rmax <- function (Nsim = 1,
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
                  surv.j = 1,
                  alpha.fixed = NULL,
                  alpha.unif = FALSE,
                  min.alpha = NULL,
                  max.alpha = NULL,
                  alpha.lognorm = FALSE,
                  mean.alpha = NULL,
                  sd.alpha = NULL) {

  ## Step 1 : Create end storage vector to save the results of the function and load functions ----
  Output <- NULL
  lambdaMax.calc <- NULL

  # lognormal distribution (source : Fred Johnson)
  MOM.lognorm <- function(expect,stdev){
    sig <- log(1+(stdev^2/expect^2))
    se <- sqrt(sig)
    est <- log(expect)-(sig*0.5)
    result <- c(est,se)
    return(result)
  }

  # beta distribution (source : Fred Johnson)
    MOM.beta <- function(mu,var) {
    sum_ab <- (mu*(1-mu)/var)-1
    a <- sum_ab*mu
    b <- sum_ab*(1-mu)
    beta <- c(a, b)
    return(beta)
  }



  # Step 2 : Generate estimated variables ----
  ## alpha (age at the first breeding) estimation ----
  if(isTRUE(alpha.unif)) { # If alpha is extracted from  a uniform distribution then :
    alpha.calc <- runif(Nsim, min.alpha, max.alpha) # 'Nsim' values are extracted from a uniform distribution with 'min.alpha' as lower bound and 'max.alpha' as upper bound
  }
  if(isTRUE(alpha.lognorm)) { # If alpha is extracted from  a normal distribution then :
    alpha.param <- MOM.lognorm(mean.alpha, sd.alpha)
    alpha.calc <- rlnorm(Nsim, alpha.param[1], alpha.param[2]) # 'Nsim' values are extracted from a log-normal distribution with with the mean 'mean.alpha' and the standard deviation 'sd.alpha'
  }
  if(!is.na(alpha.fixed)) {
    alpha.calc <- rep(alpha.fixed, Nsim)
  }
  if(!isTRUE(alpha.unif) && !isTRUE(alpha.lognorm) && is.na(alpha.fixed) && !is.na(max.alpha)){
    alpha.calc <- rep(max.alpha, Nsim)
  }




  ## survival estimation ----
  if(!is.na(surv.fixed)){ # If the value of surv.fixed is different to 'NA' (i.e. the value of surv.fixed is known) then :
    survival.calc <- rep(surv.fixed, Nsim) # The value of surv.fixed is replicated 'Nsim' times and saves in the object 'survival.calc'
  }
  if(!is.na(mass.fixed)){ # If the value of mass.fixed is different to 'NA' (i.e. the value of the body mass is known) then :
    survival.tab <- surv(Nsim, mass.fixed = mass.fixed, alpha = alpha.calc, type.p = type.p, type.e = type.e, surv.j = surv.j) # The values of survival are calculated from the surv function and save in a data frame named 'survival.tab'
    survival.calc <- survival.tab$survival # The values of survival are extract from the data frame survival.tab and save in the object 'survival.calc'
  }
  if(isTRUE(mass.lognorm)){ # If the body mass is extracted from a normal distribution then :
    survival.tab <- surv(Nsim, mass.lognorm = TRUE, mean.mass = mean.mass, sd.mass = sd.mass, alpha = alpha.calc, type.p = type.p, type.e = type.e, surv.j = surv.j) # The values of survival are calculated from the surv function and save in a data frame named 'survival.tab'
    survival.calc <- survival.tab$survival # The values of survival are extract from the data frame survival.tab and save in the object 'survival.calc'
  }
  if(isTRUE(surv.beta)){ # If survival is extracted from a normal distribution then :
    beta.param <- MOM.beta(mean.surv, sd.surv^2)
    survival.calc <- rbeta(Nsim, beta.param[1], beta.param[2]) # 'Nsim' values are extracted from a beta distribution
    survival.calc[survival.calc > 1] <- 0.9999999
    rm(beta.param)
  }




  ## lambdaMax estimation ----
  if(living.rate == "short"){ # If the species studied is a short-lived species then :
    for (b in 1:Nsim){ # Loop for calculations from 1 to the number of iterations wanted
      alpha.calcT <- alpha.calc[b] # Temporary setting the value of the object alpha.calc of the index b
      survival.calcT <- survival.calc[b] # Temporary setting the value of the object survival.calc of the index b
      lambdaMax <- lambdaMax.S(alpha.calcT,survival.calcT) # Use the lambdaMax.S function which makes it possible the lamndaMax estimation for a short-lived species
      lambdaMax.calc <- c(lambdaMax.calc, lambdaMax)
    } # End of loop b
  } else { # If the species is not a short-lived species then it implies it is a long-lived species
    lambdaMax.calc <- lambdaMax.L(alpha.calc, survival.calc) # Use of the lambdaMax.L function which makes it possible the lamndaMax estimation for a long-lived species
  }




  ## Rmax estimation ----
  Rmax <-  lambdaMax.calc - 1



  ## Step 3 : Save datas generated before in a dataframe named 'Output'---- ----
  NsimT = 1:Nsim # This object makes it possible to known the value of the index a for each row of the data frame below
  Output <- data.frame(Nsim = NsimT, Rmax = Rmax, lambdaMax = lambdaMax.calc, living.rate = living.rate,
                        survival = survival.calc, alpha = alpha.calc)  # Create the Output



  ## Step 4 : Complementary columns for the 'Output' ----
  # This step integrates columns in the final data frame according to the arguments informed
  var.storage <- c("min.alpha", "max.alpha","mean.alpha", "sd.alpha", "mass.fixed",
                   "mean.mass","sd.mass", "surv.j") # Save in an object all sumpplementary possible columns
  not.na <- which(!is.na(unlist(mget(var.storage)))) # Select columns whose values are different from NA. Works only with NA. That is why parameters which are NULL are transformed as NA
  Output[,var.storage[not.na]]<- mget(var.storage)[not.na] # Insert values columns whose values are different from NA in the data frame 'Output'
  if(exists("survival.tab")){Output$mass <- survival.tab$mass} # If the data frame 'survival.tab' exists then add the mass values
  if(exists("survival.tab")){Output$e <- survival.tab$e} # If the data frame 'survival.tab' exists then add the e values
  if(exists("survival.tab")){Output$p <- survival.tab$p} # If the data frame 'survival.tab' exists then add the p values



  ## Step 5 : Reorder columns of the 'Output'----
  var.surv <- NULL
  if(exists("survival.tab")){var.surv <- c(var.surv, "mass", "e", "p", "surv.j")} # Order of the columns wanted in the final output
  var.surv.name <- var.surv[var.surv %in% colnames(Output)] # Reorder the columns in the data frame named 'Output'

  var.alpha <- c("min.alpha", "max.alpha","mean.alpha", "sd.alpha") # Order of the columns wanted in the final output
  var.alpha.name <- var.alpha[var.alpha %in% colnames(Output)] # Reorder the columns in the data frame named 'Output'

  reordering <- c("Nsim", "Rmax", "living.rate", "lambdaMax", "survival", var.surv.name, "alpha", var.alpha.name) # Final order of the columns wanted in the final output
  Output <- data.frame(Output[, reordering]) # Reorder the columns in the data frame named 'Output'



} # End of the function
