#' Maximum annual recruitment rate estimation
#'
#' For internal use only.
#'
#' @param Nsim number, default 1 ;  number of simulations
#' @param living.rate character string ; living rate : either "\code{long}" or "\code{short}". No other character string allowed.
#' @param surv.fixed number ; fixed estimate of adult annual survival rate provided by the user. To be used instead of \code{mean.surv} when \code{surv.norm = FALSE}.
#' @param surv.norm surv.norm boolean, default FALSE. If \code{TRUE}, simulates adult annual survival rate from a normal distribution with mu = \code{mean.surv} and standard deviation = \code{sd.surv} (see below).
#' @param mean.surv mean.surv number ; arithmetic mean of adult annual survival rate. Only used if \code{surv.norm = TRUE}.
#' @param sd.surv sd.surv number ; standard deviation of adult annual survival rate. Only used if \code{surv.norm = TRUE}.
#' @param alpha.fixed alpha.fixed number ; fixed estimate of age at first breeding (in years). To be used when \code{alpha.unif = FALSE} and \code{alpha.norm = FALSE}.
#' @param alpha.unif  alpha.unif boolean, default \code{FALSE}. Simulates age at first breeding from a uniform distribution bounded by \code{min.alpha} and \code{max.alpha}.
#' @param min.alpha  min.alpha number ; minimum age at first breeding. To be used only if \code{alpha.unif = TRUE}.
#' @param max.alpha  max.alpha number ; maximum age at first breeding. To be used only if \code{alpha.unif = TRUE}.
#' @param alpha.norm  alpha.norm boolean, default \code{FALSE}. If \code{TRUE}, simulates age at first breeding from a normal distribution with mu = \code{mean.alpha} and standard deviation = \code{sd.alpha}.
#' @param mean.alpha mean.alpha number ; arithmetic mean of age at first breeding. To be used only used if \code{alpha.norm = TRUE}.
#' @param sd.alpha sd.alpha number ; age at first breeding standard deviation. Only used if \code{alpha.norm = TRUE}.
#' @param mass.fixed mass.fixed number ; fixed estimate of body mass (in kilograms) used to estimate adult survival rate from body mass.
#' @param mass.norm  mass.norm boolean, default \code{FALSE}. If \code{TRUE}, simulates body mass from a normal distribution with mu = \code{mean.mass} and standard deviation = \code{sd.mass}.
#' @param mean.mass  mean.mass number ; arithmetic mean of body mass (in kilograms). Only used if \code{mass.norm = TRUE}. #' @param sd.mass number ; standard deviation of body mass. Only used if \code{mass.norm = TRUE}.
#' @param sd.mass sd.mass number ; standard deviation of body mass. Only used if \code{mass.norm = TRUE}.
#' @param type.p type.p character string ; the value of the parameter p can be either "\code{determinist}" or "\code{random}". No other character string allowed.
#' @param type.e type.e character string ; the value of the parameter e can be either "\code{determinist}" or "\code{random}". No other character string allowed.

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
rmax <- function (Nsim = 1,
                  living.rate = NULL,
                  surv.fixed = NULL,
                  surv.norm = FALSE,
                  mean.surv = NULL,
                  sd.surv = NULL,
                  mass.fixed = NULL,
                  mass.norm = FALSE,
                  mean.mass = NULL,
                  sd.mass = NULL,
                  type.p = NULL,
                  type.e = NULL,
                  alpha.fixed = NULL,
                  alpha.unif = FALSE,
                  min.alpha = NULL,
                  max.alpha = NULL,
                  alpha.norm = FALSE,
                  mean.alpha = NULL,
                  sd.alpha = NULL) {

  ## Step 1 : Creating end storage vector to save the results of the function ----
  Output <- NULL
  lambdaMax.calc <- NULL

  # Step 2 : Generating estimated variables ----
  ## survival estimation ----
  if(!is.na(surv.fixed)){ # If the value of surv.fixed is different to 'NA' (means the value of surv.fixed is known) then :
    survival.calc <- rep(surv.fixed, Nsim) # The value of surv.fixed is replicated 'Nsim' times and saving in the object 'survival.calc'
  }
  if(!is.na(mass.fixed)){ # If the value of mass.fixed is different to 'NA' (means the value of the body mass is known) then :
    survival.tab <- surv(Nsim, mass.fixed = mass.fixed, type.p = type.p, type.e = type.e) # The values of survival are calculated from the surv function and save in a data frame named 'survival.tab'
    survival.calc <- survival.tab$survival # The values of survival are extract from the data frame survival.tab and save in the object 'survival.calc'
  }
  if(isTRUE(mass.norm)){ # If the body mass is extracted from a normal distribution then :
    survival.tab <- surv(Nsim, mass.norm = TRUE, mean.mass = mean.mass, sd.mass = sd.mass) # The values of survival are calculated from the surv function and save in a data frame named 'survival.tab'
    survival.calc <- survival.tab$survival # The values of survival are extract from the data frame survival.tab and save in the object 'survival.calc'
  }
  if(isTRUE(surv.norm)){ # If survival is extracted from a normal distribution then :
    survival.calc <- rnorm(Nsim,mean.surv,sd.surv) # 'Nsim' values are extracted from a normal distribution with the mean 'mean.surv' and the standard deviation 'sd.surv'
    survival.calc[survival.calc > 1] <- 0.9999999
  }



  ## alpha (age at the first reproduction) estimation ----
  if(isTRUE(alpha.unif)) { # If alpha is extracted from  a uniform distribution then :
    alpha.calc <- runif(Nsim, min.alpha, max.alpha) # 'Nsim' values are extracted from a uniform distribution with 'min.alpha' as lower bound and 'max.alpha' as upper bound
  }
  if(isTRUE(alpha.norm)) { # If alpha is extracted from  a normal distribution then :
    alpha.calc <- rnorm(Nsim, mean.alpha, sd.alpha) # 'Nsim' values are extracted from a normal distribution with with the mean 'mean.alpha' and the standard deviation 'sd.alpha'
  }
  if(!is.na(alpha.fixed)) {
    alpha.calc <- rep(alpha.fixed, Nsim)
  }
  if(!isTRUE(alpha.unif) && !isTRUE(alpha.norm) && is.na(alpha.fixed) && !is.na(max.alpha)){
    alpha.calc <- rep(max.alpha, Nsim)
  }



  ## lambdaMax estimation ----
  if(living.rate == "short"){ # If the species studied is a short-lived species then :
    for (b in 1:Nsim){ # Loop for calculations from 1 to the number of iterations wanted
      alpha.calcT <- alpha.calc[b] # Temporary setting the value of the object alpha.calc of the index b
      survival.calcT <- survival.calc[b] # Temporary setting the value of the object survival.calc of the index b
      lambdaMax <- lambdaMax.S(alpha.calcT,survival.calcT) # Using of the lambdaMax.S function which makes it possible the lamndaMax estimation for a short-lived species
      lambdaMax.calc <- c(lambdaMax.calc, lambdaMax)
    } # End of loop b
  } else { # If the species is not a short-lived species then it implies it is a long-lived species
    lambdaMax.calc <- lambdaMax.L(alpha.calc, survival.calc) # Using of the lambdaMax.L function which makes it possible the lamndaMax estimation for a long-lived species
  }




  ## Rmax calculation ----
  Rmax <-  lambdaMax.calc - 1



  ## Step 3 : Saving datas generated before in a dataframe named 'Output'---- ----
  NsimT = 1:Nsim # This object makes it possible to known the value of the index a for each row of the data frame below
  Output <- data.frame(Nsim = NsimT, Rmax = Rmax, lambdaMax = lambdaMax.calc, living.rate = living.rate,
                        survival = survival.calc, alpha = alpha.calc)  # Creation of the Output



  ## Step 4 : Complementary columns for the 'Output' ----
  # This step integrates columns in the final data frame according the arguments informed
  var.storage <- c("min.alpha", "max.alpha","mean.alpha", "sd.alpha", "mass.fixed",
                   "mean.mass","sd.mass") # Saving in an object all sumpplementary possible columns
  not.na <- which(!is.na(unlist(mget(var.storage)))) # Select columns whose values are different from NA. Works only with NA. That is why parameters which are NULL are transformed as NA
  Output[,var.storage[not.na]]<- mget(var.storage)[not.na] # Insert values columns whose values are different from NA in the data frame 'Output'
  if(exists("survival.tab")){Output$mass <- survival.tab$mass} # If the data frame 'survival.tab' exists then add the mass values
  if(exists("survival.tab")){Output$e <- survival.tab$e} # If the data frame 'survival.tab' exists then add the e values
  if(exists("survival.tab")){Output$p <- survival.tab$p} # If the data frame 'survival.tab' exists then add the p values



  ## Step 5 : Reordering columns of the 'Output'----
  var.surv <- NULL
  if(exists("survival.tab")){var.surv <- c(var.surv, "mass", "e", "p")} # Order of the columns wanted in the final output
  var.surv.name <- var.surv[var.surv %in% colnames(Output)] # Reordering the columns in the data frame named 'Output'

  var.alpha <- c("min.alpha", "max.alpha","mean.alpha", "sd.alpha") # Order of the columns wanted in the final output
  var.alpha.name <- var.alpha[var.alpha %in% colnames(Output)] # Reordering the columns in the data frame named 'Output'

  reordering <- c("Nsim", "Rmax", "living.rate", "lambdaMax", "survival", var.surv.name, "alpha", var.alpha.name) # Final order of the columns wanted in the final output
  Output <- data.frame(Output[, reordering]) # Reordering the columns in the data frame named 'Output'



} # End of the function
