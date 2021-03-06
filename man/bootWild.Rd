% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootWild.R
\name{bootWild}
\alias{bootWild}
\title{Wild Bootstrap Hypothesis Testing}
\usage{
bootWild(
  dependent,
  independents,
  repetitions,
  hypothesisColumn,
  hypothesisValue,
  rv,
  transform,
  suppressConstant = FALSE
)
}
\arguments{
\item{dependent}{numeric vector of dependent variable}

\item{independents}{matrix of independent variable(s), do not include constant term.}

\item{repetitions}{number of bootstrap repetitions to be performed}

\item{hypothesisColumn}{column index in 'independents' of which coefficient restriction will be tested}

\item{hypothesisValue}{hypothesized value of the coefficient of hypothesisColumn, numeric}

\item{rv}{Distribution of random variable (with mean 0 and variance 1) to multiply the residuals by. Acceptable types are 'mammen', 'rademacher', and 'standard normal'.}

\item{transform}{Type of transform to apply to residuals. Acceptable values are "simple", where transform is f(u) = u(n/(n-k))^.5 (n is the number of observations and
k is the degrees of freedom),  "leverage", where transform is f(u) = u/((1-h)^.5) (h is the diagonal of the hat matrix), or "none".}

\item{suppressConstant}{logical. If true, there will be no constant term added to the regression. Defaults to false.}
}
\value{
A character string detailing the resulting p-value of the test.
}
\description{
This function takes a vector or matrix of a dependent variable, a vector or matrix of independent variable(s), a number of bootstrap repetitions, a hypothesis, and a option to rescale the residuals
and performs the wild bootstrap for the hypothesis test (specified by hypothesisValue) on the coefficient of a variable (specified by hypothesisColumn) against its two tailed alternative
in the context of a ordinary least squares regression.
}
\details{
To learn more about the theory, please read section 3.4 of "Bootstrap Methods
in Econometrics" by James A. MacKinnon, Queen’s Economics Department Working Paper No. 1028, 2006.
}
