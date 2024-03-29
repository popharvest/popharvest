% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.summary.R
\name{output.summary}
\alias{output.summary}
\title{Summary of output parameters from the \code{PEG} and the \code{PTL} functions}
\usage{
output.summary(data)
}
\arguments{
\item{data}{dataframe generated by the \code{PEG} function or the the \code{PTL} function detailing calculations for each of the n iterations. Other input will give an error.}
}
\value{
A dataframe containing the following statistics calculated from the n iterations:

In columns :
\itemize{
  \item \code{SP} the species
}
In rows :
\itemize{
  \item \code{Fs} safety factor or \code{Fobj} objective factor
  \item \code{med.Rmax} median of the maximum annual recruitment rate
  \item \code{mean.Rmax} arithmetic mean of the maximum annual recruitment rate
  \item \code{percent2.5.Rmax} 2.5th pencentile of maximum annual recruitment rate
  \item \code{percent97.5.Rmax} 97.5th percentile of maximum annual recruitment rate
  \item \code{med.PEG} median of the potential excess growth
  \item \code{mean.PEG} arithmetic mean of the potential excess growth
  \item \code{percent2.5.PEG} 2.5th pencentile of the potential excess growth
  \item \code{percent97.5.PEG} 97.5th percentile of the potential excess growth
  \item \code{med.PTL} median of the maximum theoretical harvest
  \item \code{mean.PTL} arithmetic mean of the maximum theoretical harvest
  \item \code{percent2.5.PTL} 2.5th pencentile of the maximum theoretical harvest
  \item \code{percent97.5.PTL} 97.5th percentile of the maximum theoretical harvest
  \item \code{med.SHI} median of the sustainability harvest index
  \item \code{mean.SHI} arithmetic mean of the sustainability harvest index
  \item \code{percent2.5.SHI} 2.5th pencentile of the sustainability harvest index
  \item \code{percent97.5.SHI} 97.5th percentile of the sustainability harvest index
  \item \code{unsustain.harvest} frequency (1 = 100\%) of the n iterations for which observed harvest levels exceed the theoretical sustainable harvest levels
}
}
\description{
The function returns descriptive statistics (mean, median, se and 95 \% confidence intervals) for output parameters calculated by PEG or PTL functions. These parameters include : potential excess growth (PEG) or prescribed take level (PTL), sustainable harvest index (SHI), the maximum annual recruitment rate (Rmax) and the number of iterations for which the observed harvest level exceeds the theorical sustainable harvest (unsustain.harvest). Descriptive statistics are calculated from the \emph{n} simulations called within the \code{PEG} and \code{PTL} functions. When one or more parameters are fixed by the user, the corresponding value is reported for both the mean and the median, without any further measure of dispersion (i.e. se and 95 \% CI are set to zero). Statistics are reported in columns.
}
\examples{
df <- PEG(Nsim = 20, NSp = 2, living.rate = c("long", "short"),
surv.fixed = c(0.8, 0.65),
alpha.unif = TRUE, min.alpha = c(2, 1), max.alpha = c(3, 2),
pop.fixed = c(3605244, 55805898),
harvest.fixed = c(107802, 8447950),
Fs = c(0.1, 0.3, 0.5))

output.summary(df)

}
\seealso{
\code{\link{PEG}}, \code{\link{PTL}}
}
