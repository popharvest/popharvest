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
#'   \item \code{percent2.5.surv} 2.5% pencentile of adult annual survival rate
#'   \item \code{percent97.5.surv} 97.5% percentile of adult annual survival rate
#'   \item \code{mean.alpha} arithmetic mean of age at first breeding
#'   \item \code{med.alpha} median of age at first breeding
#'   \item \code{percent2.5.alpha} 2.5% pencentile of age at first breeding
#'   \item \code{percent97.5.alpha} 97.5% percentile of age at first breeding
#'   \item \code{mean.lambdaMax} arithmetic mean of maximum multiplicative growth rate
#'   \item \code{med.lambdaMax} median of maximum multiplicative growth rate
#'   \item \code{percent2.5.lambdaMax} 2.5% pencentile of maximum multiplicative growth rate
#'   \item \code{percent97.5.lambdaMax} 97.5% percentile of maximum multiplicative growth rate#'
#'   \item \code{mean.Rmax} arithmetic mean of maximum annual growth rate
#'   \item \code{med.Rmax} median of maximum annual growth rate
#'   \item \code{percent2.5.Rmax} 2.5% pencentile of maximum annual growth rate
#'   \item \code{percent97.5.Rmax} 97.5% percentile of maximum annual growth rate
#' }
#' @export
#' @seealso \code{\link{PEG}}, \code{\link{PTL}}.
#' @importFrom stats var median t.test quantile
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
      percent2.5.surv <- round(quantile(newdata$survival, probs=0.025), 3)
      percent97.5.surv <- round(quantile(newdata$survival, probs=0.975),3)
    }else{ # If there is no column named "survival" then assign NA to the descriptive statistics
      mean.surv <- NA
      med.surv <- NA
      percent2.5.surv <- NA
      percent97.5.surv <- NA
    }



    ## Descriptive statistics of alpha ----
    col <- which(colnames(newdata) == "alpha") # Search for a column named "alpha" in the data frame named 'newdata'
    if(length(col)> 0){ # If TRUE, then :
      mean.alpha <- round(mean(newdata$alpha), 3) # Mean rounded to the thousandth
      med.alpha <- round(median(newdata$alpha), 3) # Median rounded to the thousandth
      percent2.5.alpha <- round(quantile(newdata$alpha, probs=0.025),3)
      percent97.5.alpha <- round(quantile(newdata$alpha, probs=0.975),3)
    }else{ # If there is no column named "alpha" then assign NA to the descriptive statistics
      mean.alpha <- NA
      med.alpha <- NA
      percent2.5.alpha <- NA
      percent97.5.alpha <- NA
    }



    ## Descriptive statistics of lambdaMax ----
    col <- which(colnames(newdata) == "lambdaMax") # Search for a column named "lambdaMax" in the data frame named 'newdata'
    if(length(col)> 0){ # If TRUE, then :
      mean.lambdaMax <- round(mean(newdata$lambdaMax), 3) # Mean rounded to the thousandth
      med.lambdaMax <- round(median(newdata$lambdaMax), 3) # Median rounded to the thousandth
      percent2.5.lambdaMax <- round(quantile(newdata$lambdaMax, probs=0.025),3)
      percent97.5.lambdaMax <- round(quantile(newdata$lambdaMax, probs=0.975),3)
    }else{ # If there is no column named "alpha" then assign NA to the descriptive statistics
      mean.lambdaMax <- NA
      med.lambdaMax <- NA
      percent2.5.lambdaMax <- NA
      percent97.5.lambdaMax <- NA
    }



    ## Descriptive statistics of Rmax ----
    mean.Rmax <- round(mean(newdata$Rmax), 3) # Mean rounded to the thousandth
    med.Rmax <- round(median(newdata$Rmax), 3) # Median rounded to the thousandth
    percent2.5.Rmax <- round(quantile(newdata$Rmax, probs=0.025),3)
    percent97.5.Rmax <- round(quantile(newdata$Rmax, probs=0.975),3)

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
          percent2.5.Rmax = percent2.5.Rmax,
          percent97.5.Rmax = percent97.5.Rmax)



    ## Condition "else" : If there is no column named 'Fs' in the data frame named 'data' ----
    # In this case, there is a column named 'Fobj' exists
    } else {
      Fobj <- unique(newdata$Fobj) # Delete duplicates in the column 'Fs' in the data frame named 'newdata'
      OutputT <- data.frame( # Create a data frame named 'OutputT' with descriptive statistics of Rmax
          code = newdata$code[1],
          Fobj = Fobj,
          med.Rmax = med.Rmax,
          mean.Rmax = mean.Rmax,
          percent2.5.Rmax = percent2.5.Rmax,
          percent97.5.Rmax = percent97.5.Rmax)
    } # End of condition "if else"



    ## Step 5 : Complementary columns for the data frame named 'OutputT'
    # This step integrates columns in the final data frame according to the arguments informed
    var.storage <- c("med.surv", "mean.surv", "percent2.5.surv", "percent97.5.surv",
                     "med.alpha", "mean.alpha", "percent2.5.alpha", "percent97.5.alpha",
                     "med.lambdaMax", "mean.lambdaMax", "percent2.5.lambdaMax", "percent97.5.lambdaMax") # Save in an object all sumpplementary possible columns
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
