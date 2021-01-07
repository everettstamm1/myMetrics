#' Pairs Bootstrap Confidence Intervals
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
#' @param alpha Specified confidence level of the interval.
#' @param suppressConstant logical. If true, there will be no constant term added to the regression. Defaults to false.
#' @return Numerical pair of the lower and upper bound of the confidence interval.
#' @export
bootPairsCI <- function(dependent,independents,repetitions,hypothesisColumn,alpha,suppressConstant = FALSE){
  #Making sure all the inputs are cool and if dependent or independent are vectors, making the matrices
  dependent <- checkDependent(dependent)
  hypothesisColumnSkipCons <- hypothesisColumn
  hypothesisColumn <- checkSuppress(suppressConstant,hypothesisColumn)
  independents <- checkIndependent(independents,suppressConstant)
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
  y_x_star <- setupPairsCI(dependent,independents,repetitions,n)



  #Creating t-values for each bootstrap iteration
  t_vals <- runPairsCI(y_x_star[[1]],y_x_star[[2]],hypothesisColumn,b_ols,c(),n,k)
  t_vals <- sort(t_vals)
  print(t_vals)

  if((alpha/2)*(repetitions+1) <1)
    return(stop("Too large a confidence interval for number of sets. As critical indexes are truncated, the calculated lower critical value has index of zero (index starts at one) in the vector of bootstrap t-values. Increase alpha or number of repetitions."))
  crit_low <- t_vals[(alpha/2)*(repetitions+1)]
  crit_high <- t_vals[(1-alpha/2)*(repetitions+1)]
  ci <- c(b_ols[hypothesisColumn]-crit_high*se_ols[hypothesisColumn],b_ols[hypothesisColumn]-crit_low*se_ols[hypothesisColumn])


  print(paste0("The alpha = ",alpha," confidence interval is [",ci[1],", ",ci[2],"] for the coefficient of variable ",  hypothesisColumnSkipCons," in the matrix of independents using ",repetitions," bootstrap repetitions."))
  return(ci)
}
setupPairsCI <- function(dep,indep,repetitions,n){
  if(is.null(nrow(indep))){
    newIndep <- replicate(repetitions, indep[sample.int(n,replace = TRUE)])
  }else{
    newIndep <- replicate(repetitions, indep[sample.int(n,replace = TRUE),])
  }
  newDep <- replicate(repetitions, dep[sample.int(n,replace = TRUE)])
  return(list(newDep,newIndep))
}

runPairsCI <- function(Y_stars,X_stars,index,b_ols,t_vec,n,k){
  if(is.null(ncol(Y_stars))){
    invXtX <- solve(t(X_stars)%*%X_stars)
    invXtX_Xt <- invXtX %*%t(X_stars)

    b_star <- invXtX_Xt %*% Y_stars
    u_star <- Y_stars - X_stars %*% b_star
    se_star <- sqrt(diag(c((t(u_star) %*% u_star)/(n-k)) * invXtX))
    t_star <- (b_star[index]-b_ols[index])/se_star[index]
    t_vec <- c(t_vec,t_star)
    return(t_vec)
  }else{
    if(length(dim(X_stars)) ==2){
      x_star <- X_stars[,1]
      X_stars <- X_stars[,-c(1)]
    }else{
      x_star <- X_stars[,,1]
      X_stars <- X_stars[,,-c(1)]
    }
    invXtX <- solve(t(x_star)%*%x_star)
    invXtX_Xt <- invXtX %*% t(x_star)

    b_star <- invXtX_Xt %*% Y_stars[,1]
    u_star <- Y_stars[,1] - x_star %*% b_star
    se_star <- sqrt(diag(c((t(u_star) %*% u_star)/(n-k)) * invXtX))

    t_star <- (b_star[index]-b_ols[index])/se_star[index]
    t_vec <- c(t_vec,t_star)
    Y_stars <- Y_stars[,-c(1)]
    return(runPairsCI(Y_stars,X_stars,index,b_ols,t_vec,n,k))
  }
}
