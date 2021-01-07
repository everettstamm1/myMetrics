#' Wild Bootstrap Confidence Intervals
#'
#' This function takes a vector or matrix of a dependent variable, a vector or matrix of independent variable(s), a number of bootstrap repetitions, a hypothesis, and a option to rescale the residuals
#' and creates confidence intervals using the wild bootstrap for the coefficient of a variable (specified by hypothesisColumn) in the context of a ordinary least squares regression.
#'
#' To learn more about the theory, please read "Econometric Theory and Methods, Oxford University Press, USA (2003).
#'
#'
#' @param dependent numeric vector of dependent variable
#' @param independents matrix of independent variable(s), do not include constant term.
#' @param repetitions number of bootstrap repetitions to be performed
#' @param hypothesisColumn column index in 'independents' of which coefficient restriction will be tested
#' @param rv Distribution of random variable (with mean 0 and variance 1) to multiply the residuals by. Acceptable types are 'mammen', 'rademacher', and 'standard normal'.
#' @param alpha Specified confidence level of the interval.
#' @param transform Type of transform to apply to residuals. Acceptable values are "simple", where transform is f(u) = u(n/(n-k))^.5 (n is the number of observations and
#' k is the degrees of freedom),  "leverage", where transform is f(u) = u/((1-h)^.5) (h is the diagonal of the hat matrix), or "none".
#' @param suppressConstant logical. If true, there will be no constant term added to the regression. Defaults to false.
#' @return Numerical pair of the lower and upper bound of the confidence interval.
#' @export
bootWildCI <- function(dependent,independents,repetitions,hypothesisColumn,alpha,rv,transform,suppressConstant = FALSE){
  #Making sure all the inputs are cool and if dependent or independent are vectors, making the matrices
  dependent <- checkDependent(dependent)
  hypothesisColumnSkipCons <- hypothesisColumn
  hypothesisColumn <- checkSuppress(suppressConstant,hypothesisColumn)
  independents <- checkIndependent(independents,suppressConstant)
  if(transform != "simple" & transform != "leverage" & transform != "none")
    return(stop("transform must be either 'simple', 'leverage', or 'none'"))
  checkReps(repetitions)
  if(nrow(dependent) != nrow(independents))
    return(stop("dependent and independent matrices have different number of observations."))
  if(ncol(independents) == 1 & suppressConstant){
    return(stop("Cannot have number of independent variables be one and have suppressConstant==TRUE"))
  }

  #Premultiplying (Xt*X)^-1 and ((Xt*X)^-1)*Xt so it doesn't need to be done at every recursion of the helper functions
  invXtX <- solve(t(independents) %*% independents)
  invXtX_Xt <- invXtX %*% t(independents)

  #Creating base sample estimates
  n <- nrow(independents)
  k <- ncol(independents)
  b_ols <- invXtX_Xt %*% dependent
  u_ols <- dependent - independents %*% b_ols
  se_ols <- sqrt(diag(c((t(u_ols) %*% u_ols)/(n-k)) * invXtX))


  #Creating bootstrapped dependent values
  y_u_star <- setupWildCI(dependent,independents,b_ols,u_ols,repetitions,rv,invXtX_Xt, transform)




  #Creating t-values for each bootstrap iteration
  t_vals <- runWildCI(y_u_star[[1]],y_u_star[[2]],independents,hypothesisColumn,b_ols,invXtX,invXtX_Xt,t_vec = c())
  t_vals <- sort(t_vals)

  crit_low <- t_vals[(alpha/2)*(repetitions+1)]
  crit_high <- t_vals[(1-alpha/2)*(repetitions+1)]

  ci <- c(b_ols[hypothesisColumn]-crit_high*se_ols[hypothesisColumn],b_ols[hypothesisColumn]-crit_low*se_ols[hypothesisColumn])


  print(paste0("The alpha = ",alpha," confidence interval is [",ci[1],", ",ci[2],"] for the coefficient of variable ",hypothesisColumnSkipCons," in the matrix of independents using ",repetitions," bootstrap repetitions."))
  return(ci)
}
setupWildCI <- function(dep,indep,b_ols,u_ols,repetitions,rv,invXtX_Xt,transform){
  n <- length(dep)
  k <- ncol(indep)

  v <- randomVariableMatrix(rv, n,repetitions)
  u_star <- v*c(u_ols)
  if(transform == "simple"){
    u_star <- u_star*((n/(n-k))^.5)
  }else if(transform == "leverage"){
    h <- diag(indep %*%invXtX_Xt)
    u_star <- u_star/(1-h)^.5
  }
  Y_star <- c(indep %*% b_ols) + u_star
  return(list(Y_star,u_star))
}
randomVariable <- function(type,n){
  if(type == "mammen" | type == "Mammen"){
    high <- (sqrt(5)+1)/2
    low <- -(sqrt(5)-1)/2
    v <- runif(n)
    v[v > (1-(sqrt(5)-1)/(2*sqrt(5)))] <- high
    v[v <= (1-(sqrt(5)-1)/(2*sqrt(5)))] <- low
    return(v)
  }else if(type == "rademacher"|| type == "Rademacher"){
    high <- 1
    low <- -1
    v <- runif(n)
    v[v > .5] <- high
    v[v <= .5] <- low
    return(v)
  }else if(type =="Standard Normal" | type =="Standard normal" | type =="standard Normal"){
    v <- rnorm(n,0,1)
    return(v)
  }else {
    return(stop("Variable rv ncorrectly specified. Must be either 'mammen', 'rademacher', or 'standard normal"))
  }
}

#Creates an n x m matrix of m repetitions of n drawings from the distribution specified by "type"
randomVariableMatrix <- function(type,n,m){
  dims <- as.matrix(rep(n,m))
  v <- apply(dims,1,function(x){randomVariable(type,x)})
  return(v)
}

runWildCI <- function(Y_stars,u_stars,indep,index,b_ols,invXtX,invXtX_Xt,t_vec){
  if(is.null(ncol(Y_stars))){
    b_star <- invXtX_Xt %*% Y_stars
    u_star <- Y_stars - indep %*% b_star
    se_star <- sqrt(diag(c((t(u_stars) %*% u_stars)/(nrow(indep)-ncol(indep))) * invXtX))
    t_star <- (b_star[index]-b_ols[index])/se_star[index]
    t_vec <- c(t_vec,t_star)
    return(t_vec)
  }else{
    b_star <- invXtX_Xt %*% Y_stars[,1]
    u_star <- Y_stars[,1] - indep %*% b_star
    se_star <- sqrt(diag(c((t(u_stars[,1]) %*% u_stars[,1])/(nrow(indep)-ncol(indep))) * invXtX))
    t_star <- (b_star[index]-b_ols[index])/se_star[index]
    t_vec <- c(t_vec,t_star)
    Y_stars <- Y_stars[,-c(1)]
    u_stars <- u_stars[,-c(1)]
    return(runWildCI(Y_stars,u_stars,indep,index,b_ols,invXtX,invXtX_Xt,t_vec))
  }
}
