% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stmc.R
\name{stmc}
\alias{stmc}
\title{Single-Test (Single-Reader) Multi-Case ROC Analysis}
\usage{
stmc(response, case, data, cov = jackknife)
}
\arguments{
\item{response}{response metric expressed in terms of a package-supplied
performance \code{\link[=metrics]{metric}}.}

\item{case}{optional variable of case identifiers.}

\item{data}{data frame containing the \code{response}, \code{test},
\code{reader}, and \code{case} variables.}

\item{cov}{function, function call, or character string naming the
\code{\link[=cov_methods]{method}} to use in calculating performance
metric covariances.}
}
\value{
Returns a \code{stmc} class object with the following elements.
\describe{
  \item{\code{metric}}{character name of the performance metric.}
  \item{\code{est}}{estimated mean.}
  \item{\code{se}}{estimated standard error.}
  \item{\code{stmc_data}}{data frame of case-specific ratings.}
}
}
\description{
Estimation of ROC performance metrics for a single test of multiple cases.
}
\examples{
VanDyke11 <- subset(VanDyke, treatment == "1" & reader == "1")
est <- stmc(empirical_auc(truth, rating), data = VanDyke11)
plot(est)
summary(est)

}
\seealso{
\code{\link{metrics}}, \code{\link{cov_methods}},
\code{\link{parameters}}, \code{\link{plot}}, \code{\link{roc_curves}},
\code{\link{summary}}
}
