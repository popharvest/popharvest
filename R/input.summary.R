#' Summary of the input paramters used by \code{PEG} and \code{PTL} functions
#'
#' The function returns descriptive statistics (mean, median, se and 95 \% confidence intervals) for the input parameters used by PEG and PTL functions. These parameters are : annual survival rate of adults (surv), age at first breeding (alpha), maximum multiplicative growth rate (lambdaMax) and maximum annual population growth rate (Rmax). Descriptive statistics are drawn from n simulations called within the \code{PEG} and \code{PTL} functions. When one or several parameters are fixed by the user, the corresponding value is reported for mean and median, without any further measure of dispersion (i.e. se and 95 \% CI are equal to zero). Statistics are reported in columns.
#'
#' @param data dataframe from \code{PEG} and \code{PTL} functions input datas. Other inputs will give an error.
#'
#' @return A dataframe with parameters and dispersion statistics in raw and, species (SP) and safety factor (Fs) or objective factor (Fobj) in column :
#' In columns :
#' \itemize{
#'   \item \code{SP} the species
#' }
#' In rows :
#' \itemize{
#'   \item \code{Fs} safety factor
#'   \item \code{Fobj} objective factor
#'   \item \code{mean.surv} arithmetic mean of adult annual survival rate
#'   \item \code{med.surv} median of adult annual survival rate
#'   \item \code{se.surv} standard error of adult annual survival rate
#'   \item \code{surv.CI95.low } lower bound of the 95 \% confidence interval for annual survival rate
#'   \item \code{surv.CI95.up} upper bound of the 95 \% confidence interval for annual survival rate
#'   \item \code{mean.alpha} arithmetic mean of age at first breeding
#'   \item \code{med.alpha} median of age at first breeding
#'   \item \code{se.alpha} standard error of age at first breeding
#'   \item \code{alpha.CI95.low } lower bound of the 95 \% confidence interval for age at first breeding
#'   \item \code{alpha.CI95.up} upper bound of the 95 \% confidence interval for age at first breeding
#'   \item \code{mean.Rmax} arithmetic mean of maximum annual growth rate
#'   \item \code{med.Rmax} median of the maximum annual growth rate
#'   \item \code{se.Rmax} standard error of the maximum annual growth rate
#'   \item \code{Rmax.CI95.low } lower bound of the 95 \% confidence interval of the maximum annual growth rate
#'   \item \code{Rmax.CI95.up} upper bound of the 95 \% confidence interval of the maximum annual growth rate
#' }
#' @export
#' @seealso \code{\link{PEG}}, \code{\link{PTL}}.
#' @importFrom stats var median t.test
#' @importFrom stringr str_sub
#'
#' @examples
#' df <- PEG(Nsim = 20, NSp = 2, living.rate = c("long", "short"),
#' surv.fixed = c(0.8, 0.65),
#' alpha.unif = TRUE, min.alpha = c(2, 1), max.alpha = c(3, 2),
#' pop.fixed = c(3605244, 55805898),
#' harvest.fixed = c(107802, 8447950),
#' Fs = c(0.1, 0.3, 0.5))
#'
#' input.summary(df)
#'
input.summary <- function(data){

  ## I. Error messages ----
  # Print error messages if the vector types of the datas are not correct

  if(is.null(data)){
    stop("'data' must be specified") # If 'data' argument is null then print error message
  }
  if (!is.data.frame(data)) { # If 'data' argument is not a data frame then print error message
    stop("'data' must be a dataframe")
  }



  ## II. Calculations ----
  ## standard error function ----
  se <- function(x) sqrt(var(x)/length(x))



  ## Step 1 : Create the end storage vector to save results ----
  Output <- NULL



  ## Step 2 : Create the input data filter parameter for the loop of calculations bellow ----
  # This step is necessary to adapt the final data frame according to the PEG or PTL function used previously
  col <- which(colnames(data) == "Fs") # Search for a "Fs" columns
  if(length(col)> 0) {  # If TRUE, the PEG function was used previously
    data$code <- as.factor(paste0("SP", data$SP, " Fs ", data$Fs)) # Create a code for calculations in the step 3
  } else { # If the output comes from PTL function
    data$code <- as.factor(paste0("SP", data$SP, " Fobj ", data$Fobj)) # Create a code for calculations in the step 3
  }



  ## Step 3 : Loop of calulcations  ----
  # Descriptive statistics of the parameters used as input data
  for (j in levels(data$code)){
    newdata <- subset(data, data$code == j) # Extract the data frame for which the value of the 'code' column is equal to the level j



    ## Descriptive statistics of survival ----
    col <- which(colnames(newdata) == "survival") # Search for a column named "survival" in the data frame named 'newdata'
    if(length(col)> 0){ # If TRUE then :
      mean.surv <- round(mean(newdata$survival), 3) # Mean rounded to the thousandth
      med.surv <- round(median(newdata$survival), 3) # Median rounded to the thousandth
      if(mean.surv == round(newdata$survival[1], 3)){ # If the mean is equal to the first row of the column named 'survival' (i.e. all values in the 'survival' column are the same) then :
        se.surv <- 0 # Assign 0 to the standard error object
        surv.CI95.low <- 0 # Assign 0 to the lower bound of the confidence interval at 95 \%
        surv.CI95.up <- 0 # Assign 0 to the upper bound of the confidence interval at 95 \%
        # At this step we must assign 0 manually because otherwise an error message is returned
      }else{ # If the mean is different to the first row of the column named 'survival' then :
        se.surv <- round(se(newdata$survival),3) # Standard error rounded to the thousandth
      CI.list <- t.test(newdata$survival) # Confidence interval at 95 \%
      CI <- CI.list[[4]] # Extract the fourth row from the list which returns the lower and the upper bound of the confidence interval at 95 \%
      surv.CI95.low <- round(CI[1], 3) # Extract the lower bound of the confidence interval at 95 \% rounded to the thousandth
      surv.CI95.up <- round(CI[2], 3) # Exract the upper bound of the confidence interval at 95 \% rounded to the thousandth
      }
    }else{ # If there is no column named "survival" then assign NA to the descriptive statistics
      mean.surv <- NA
      med.surv <- NA
      se.surv <- NA
      surv.CI95.low <- NA
      surv.CI95.up <- NA
    }



    ## Descriptive statistics of alpha ----
    col <- which(colnames(newdata) == "alpha") # Search for a column named "alpha" in the data frame named 'newdata'
    if(length(col)> 0){ # If TRUE, then :
      mean.alpha <- round(mean(newdata$alpha), 3) # Mean rounded to the thousandth
      med.alpha <- round(median(newdata$alpha), 3) # Median rounded to the thousandth
      if(mean.alpha == round(newdata$alpha[1], 3)){ # If the mean is equal to the first row of the column named 'alpha' (i.e. all values in the 'alpha' column are the same) then :
        se.alpha <- 0 # Assign 0 to the standard error
        alpha.CI95.low <- 0 # Assign 0 to the lower bound of confidence interval at 95 \%
        alpha.CI95.up <- 0 # Assign 0 to the upper bound of confidence interval at 95 \%
        # At this step we must assign 0 manually because otherwise an error message is returned
      }else{ # If the mean is different to the first row of the column named 'alpha' then :
        se.alpha <- round(se(newdata$alpha),3) # Standard error rounded to the thousandth
        CI.list <- t.test(newdata$alpha) # Confidence interval at 95 \%
        CI <- CI.list[[4]] # Extract the fourth row from the list which returns the lower and the upper bound of the confidence interval at 95 \%
        alpha.CI95.low <- round(CI[1], 3) # Extract the lower bound of the confidence interval at 95 \% rounded to the thousandth
        alpha.CI95.up <- round(CI[2], 3) # Extract the upper bound of the confidence interval at 95 \% rounded to the thousandth
      }
    }else{ # If there is no column named "alpha" then assign NA to the descriptive statistics
      mean.alpha <- NA
      med.alpha <- NA
      se.alpha <- NA
      alpha.CI95.low <- NA
      alpha.CI95.up <- NA
    }



    ## Descriptive statistics of lambdaMax ----
    col <- which(colnames(newdata) == "lambdaMax") # Search for a column named "lambdaMax" in the data frame named 'newdata'
    if(length(col)> 0){ # If TRUE, then :
      mean.lambdaMax <- round(mean(newdata$lambdaMax), 3) # Mean rounded to the thousandth
      med.lambdaMax <- round(median(newdata$lambdaMax), 3) # Median rounded to the thousandth
      if(mean.lambdaMax == round(newdata$lambdaMax[1], 3)){ # If the mean is equal to the first row of the column named 'lambdaMax' (i.e. all values in the 'lambdaMax' column are the same) then :
        se.lambdaMax <- 0 # Assign 0 to the standard error
        lambdaMax.CI95.low <- 0 # Assign 0 to the to the lower bound of confidence interval at 95 \%
        lambdaMax.CI95.up <- 0 # Assign 0 to the to the lower bound of confidence interval at 95 \%
        # At this step we must assign 0 manually because otherwise an error message is returned
      }else{ # If the mean is different to the first row of the column named 'lambdaMax' then :
        se.lambdaMax <- round(se(newdata$lambdaMax),3) # Standard error rounded to the thousandth
        CI.list <- t.test(newdata$lambdaMax) # Confidence interval at 95 \%
        CI <- CI.list[[4]] # Extract the fourth row from the list which returns the lower and the upper bound of the confidence interval at 95 \%
        lambdaMax.CI95.low <- round(CI[1], 3) # Extract the lower bound of the confidence interval at 95 \% rounded to the thousandth
        lambdaMax.CI95.up <- round(CI[2], 3) # Extract the upper bound of the confidence interval at 95 \% rounded to the thousandth
      }
    }else{ # If there is no column named "alpha" then assign NA to the descriptive statistics
      mean.lambdaMax <- NA
      med.lambdaMax <- NA
      se.lambdaMax <- NA
      lambdaMax.CI95.low <- NA
      lambdaMax.CI95.up <- NA
    }



    ## Descriptive statistics of Rmax ----
    mean.Rmax <- round(mean(newdata$Rmax), 3) # Mean rounded to the thousandth
    med.Rmax <- round(median(newdata$Rmax), 3) # Median rounded to the thousandth
    if(mean.Rmax == round(newdata$Rmax[1], 3)){ # If the mean is equal to the first row of the column named 'Rmax' (i.e. all values in the 'Rmax' column are the same) then :
    se.Rmax <- 0 # Assign 0 to the standard error
      Rmax.CI95.low <- 0 # Assign 0 to the to the lower bound of confidence interval at 95 \%
      Rmax.CI95.up <- 0 # Assign 0 to the to the upper bound of confidence interval at 95 \%
    }else{ # If the mean is different to the first row of the column named 'Rmax' then :
      se.Rmax <- round(se(newdata$Rmax),3) # Standard error rounded to the thousandth
      CI.list <- t.test(newdata$Rmax) # Confidence interval at 95 \%
      CI <- CI.list[[4]] # Extract the fourth row from the list which returns the lower and the upper bound of the confidence interval at 95 \%
      Rmax.CI95.low <- round(CI[1], 3) # Extract the lower bound of the confidence interval at 95 \% rounded to the thousandth
      Rmax.CI95.up <- round(CI[2], 3) # Extract the upper bound of the confidence interval at 95 \% rounded to the thousandth
    }

    ##  Step 4 :  Outputs  ----
    # Save temporary output from the temporary data frame named 'newdata'
    col <- which(colnames(newdata) == "Fs") # Search for a column named "Fs" in the data frame named 'newdata' : if TRUE, it means the function previously used is the PEG function and the vector returns a number

    ## Condition "if" : # If the object 'col' is > 0 ----
    # In this case, there is a column named "Fs" in the data frame 'newdata' and the PEG function was used previously
    if (length(col) > 0){
      Fs <- unique(newdata$Fs) # Delete duplicates in the column 'Fs' in the data frame named 'data'
      OutputT <- data.frame( # Create a data frame named 'OutputT' containing descriptive statistics of Rmax
          code = newdata$code[1],
          Fs = Fs,
          med.Rmax = med.Rmax,
          mean.Rmax = mean.Rmax,
          se.Rmax = se.Rmax,
          Rmax.CI95.low = Rmax.CI95.low,
          Rmax.CI95.up = Rmax.CI95.up)



    ## Condition "else" : If there is no column named 'Fs' in the data frame named 'data' ----
    # In this case, there is a column named 'Fobj' exists
    } else {
      Fobj <- unique(newdata$Fobj) # Delete duplicates in the column 'Fs' in the data frame named 'newdata'
      OutputT <- data.frame( # Create a data frame named 'OutputT' with descriptive statistics of Rmax
          code = newdata$code[1],
          Fobj = Fobj,
          med.Rmax = med.Rmax,
          mean.Rmax = mean.Rmax,
          se.Rmax = se.Rmax,
          Rmax.CI95.low = Rmax.CI95.low,
          Rmax.CI95.up = Rmax.CI95.up
        )
    } # End of condition "if else"



    ## Step 5 : Complementary columns for the data frame named 'OutputT'
    # This step integrates columns in the final data frame according to the arguments informed
    var.storage <- c("med.surv", "mean.surv", "se.surv", "surv.CI95.low", "surv.CI95.up",
                     "med.alpha", "mean.alpha", "se.alpha", "alpha.CI95.low", "alpha.CI95.up",
                     "med.lambdaMax", "mean.lambdaMax", "se.lambdaMax", "lambdaMax.CI95.low", "lambdaMax.CI95.up") # Save in an object all sumpplementary possible columns
    not.na <- which(!is.na(unlist(mget(var.storage)))) # Select columns which do not have NA values. That is why NULL parameters are renamed NA previously.
    OutputT[,var.storage[not.na]]<- mget(var.storage)[not.na] # insert columns which do not have NA values in the temporary data frame named 'OutputT'



    ## Step 6 : Final ouptut ----
    Output <- rbind(Output, OutputT) # Combine by rows all temporary data frame for each level of j

  } # End of loop j



  ## Step 7 : Transpose the dataframe ----
  code <- as.character(Output$code)
  code <- str_sub(code, 1, 3) # Extract a character string which informs the number of the 'SP'
  Output <- t(Output) # Transpose the data frame named 'Output'
  Output <- as.data.frame(Output) # Convert the matrix named 'Output' as a data frame



  ## Step 8 : Rename columns ----
  for(i in 1:ncol(Output)){
    name.col <- code[i]
    colnames(Output)[i] <- name.col
  }



  ## Step 9 : Delete the first row of the table which does not correspond to calculations ----
  Output <- Output[2:nrow(Output),]

  print(Output)



} # End of function
