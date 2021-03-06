% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{deprecated}
\alias{deprecated}
\alias{proproc_auc}
\alias{proproc_sens}
\alias{proproc_spec}
\title{Deprecated Functions}
\usage{
proproc_auc(...)

proproc_sens(...)

proproc_spec(...)
}
\arguments{
\item{...}{arguments passed to non-deprecated equivalent.}
}
\description{
Functions that have been deprecated and will be removed in the future.
}
\details{
Use \code{\link[=binormalLR_auc]{binormalLR_auc()}} instead of
\code{proproc_auc()}.

Use \code{\link[=binormalLR_sens]{binormalLR_sens()}} instead of
\code{proproc_sens()}.

Use \code{\link[=binormalLR_spec]{binormalLR_spec()}} instead of
\code{proproc_spec()}.
}
