% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootPairsCI.R
\name{bootPairsCI}
\alias{bootPairsCI}
\title{Pairs Bootstrap Confidence Intervals}
\usage{
bootPairsCI(
  dependent,
  independents,
  repetitions,
  hypothesisColumn,
  alpha,
  suppressConstant = FALSE
)
}
\arguments{
\item{dependent}{numeric vector of dependent variable}

\item{independents}{matrix of independent variable(s), do not include constant term.}

\item{repetitions}{number of bootstrap repetitions to be performed}

\item{hypothesisColumn}{column index in 'independents' of which coefficient restriction will be tested}

\item{alpha}{Specified confidence level of the interval.}

\item{suppressConstant}{logical. If true, there will be no constant term added to the regression. Defaults to false.}
}
\value{
Numerical pair of the lower and upper bound of the confidence interval.
}
\description{
This function takes a vector or matrix of a dependent variable, a vector or matrix of independent variable(s), a number of bootstrap repetitions, a hypothesis, and a option to rescale the residuals
and creates confidence intervals using the wild bootstrap for the coefficient of a variable (specified by hypothesisColumn) in the context of a ordinary least squares regression.
}
\details{
To learn more about the theory, please read "Econometric Theory and Methods, Oxford University Press, USA (2003).
}
