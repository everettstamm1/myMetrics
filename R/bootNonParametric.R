#' Non-Parametric Restricted Residual Resampling Bootstrap
#'
#' This function takes a vector of a dependent variable, a vector of a single independent variable or a list of
#' several independent variables, a number of bootstrap repetitions, a hypothesis, and a option to rescale the residuals
#' and performs restricted residual resampling. To learn more about the theory, please read section 3.1 of "Bootstrap Methods
#' in Econometrics" by James A. MacKinnon, Queen’s Economics Department Working Paper No. 1028, 2006.
#'
#' @importFrom stats lm na.fail rnorm runif
#' @param dependent numeric vector of dependent variable
#' @param independents matrix of independent variable(s), do not include constant term.
#' @param repetitions number of bootstrap repetitions to be performed
#' @param hypothesisColumn column index in 'independents' of which coefficient restriction will be tested
#' @param hypothesisValue hypothesized value of the coefficient of hypothesisColumn, numeric
#' @param transform Type of transform to apply to residuals. Acceptable values are "simple", where transform is f(u) = u(n/(n-k))^.5 (n is the number of observations and
#' k is the degrees of freedom),  "leverage", where transform is f(u) = u/((1-h)^.5) (h is the diagonal of the hat matrix), or "none".
#' @param replace Whether or not to sample residuals with replacement. Logical
#' @param suppressConstant logical. If true, there will be no constant term added to the regression. Defaults to false.
#' @export
bootNonParametric <- function(dependent,independents,repetitions,hypothesisColumn,hypothesisValue,transform,replace, suppressConstant){
  #Making sure all the inputs are cool
  dependent <- checkDependent(dependent)
  hypothesisColumn <- checkSuppress(suppressConstant,hypothesisColumn)
  independents <- checkIndependent(independents,suppressConstant)
  checkHypothesis(independents,hypothesisColumn,hypothesisValue)
  checkReps(repetitions)
  if(transform != "simple" & transform != "leverage" & transform != "none")
    return(stop("transform must be either 'simple', 'leverage', or 'none'"))
  if(nrow(dependent) != nrow(independents))
    return(stop("dependent and independent matrices have different number of observations."))
  if(ncol(independents) == 1 & suppressConstant){
    return(stop("Cannot have number of independent variables be one and have suppressConstant==TRUE"))
  }
  if(replace != TRUE & replace != FALSE){
    return(stop("replace must be of type 'logical'"))
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
  t_ols <- (b_ols[hypothesisColumn]-hypothesisValue)/se_ols[hypothesisColumn]


  #Creating Bootstrap Y and u values
  y_u_star <- bootThis(dependent,independents,hypothesisColumn,hypothesisValue, repetitions,invXtX_Xt,transform,replace)


  #Creating t-values for each bootstrap iteration
  t_vals <- runThis(y_u_star[[1]],y_u_star[[2]],independents,hypothesisColumn,hypothesisValue,invXtX,invXtX_Xt,t_vec = c())

  #Creating Bootstrap P Values
  p <- boot_p(t_vals,t_ols)

  print(paste0("The P-value is ",p," for ",repetitions," bootstrap repetitions of the null hypothesis against the two tailed alternative."))
  return(p)

}
#' @noRd
checkSuppress <- function(s,hypCol){
  if(s == TRUE){
    return(hypCol)
  }else if(s == FALSE){
    hypCol <- hypCol+1
    return(hypCol)
  }else{
    return(stop("suppressConstant must be of type logical"))
  }
}
boot_p <- function(t_boot, t_base){
  B <- length(t_boot)
  F_val <- (1/B)*sum(t_base > t_boot)
  return(2*min(1-F_val,F_val))
}
bootThis <- function(dep,indep,hypInd,hypVal,reps,invXtX_Xt,rescale,replace,tVec){
  n <- length(dep)
  k <- ncol(indep)
  Y_res <- dep - hypVal*indep[,hypInd]
  X_res <- NULL
  b_res <- NULL
  u_res <- NULL
  if(ncol(indep) >1){
    X_res <- indep[,-c(hypInd)]
    b_res <- solve(t(X_res)%*%X_res) %*%t(X_res) %*% Y_res
    u_res <- Y_res - X_res %*% b_res
  }else{
    u_res <- Y_res
  }
  if(rescale == "simple"){
    u_res <- u_res*(((n/(n-k))^.5))
  }else if(rescale == "leverage"){
    h <- diag(indep %*%invXtX_Xt)
    u_res <- u_res/(1-h)^.5
  }
  u_star <- resampleMatrix(u_res,reps,replace)
  if(ncol(indep) >1){
    Y_star <- c(hypVal*indep[,hypInd] + X_res %*% b_res) + u_star
  }else{
    Y_star <- c(hypVal*indep[,hypInd]) + u_star
  }
  return(list(Y_star,u_star))
}
runThis <- function(Y_stars,u_stars,indep,index,h,invXtX,invXtX_Xt,t_vec){
  if(is.null(ncol(Y_stars))){
    b_star <- invXtX_Xt %*% Y_stars
    u_star <- Y_stars - indep %*% b_star
    se_star <- sqrt(diag(c((t(u_stars) %*% u_stars)/(nrow(indep)-ncol(indep))) * invXtX))
    t_star <- (b_star[index]-h)/se_star[index]
    t_vec <- c(t_vec,t_star)
    return(t_vec)
  }else{
    b_star <- invXtX_Xt %*% Y_stars[,1]
    u_star <- Y_stars[,1] - indep %*% b_star
    se_star <- sqrt(diag(c((t(u_stars[,1]) %*% u_stars[,1])/(nrow(indep)-ncol(indep))) * invXtX))
    t_star <- (b_star[index]-h)/se_star[index]
    t_vec <- c(t_vec,t_star)
    Y_stars <- Y_stars[,-c(1)]
    u_stars <- u_stars[,-c(1)]
    return(runWild(Y_stars,u_stars,indep,index,h,invXtX,invXtX_Xt,t_vec))
  }
}
resampleMatrix <- function(u,reps,replacement){
  n <- length(u)
  dims <- as.matrix(rep(n,reps))
  v <- apply(dims,1,function(x){sample(u,x,replacement)})
  return(v)
}





checkReps <- function(a){
  if(class(a) != "numeric" & class(a) != "integer"){
    return(stop("Repetitions must be a number."))
  }else if(a <=0){
    return(stop("Repetitions must be greater than zero"))
  }
}

#Check if dependent is a number and if it has missing values
checkDependent <- function(a){
  if(is.vector(a)){
    a <- as.matrix(a)
  }
  if(typeof(a) != "numeric" & typeof(a) != "integer" &typeof(a) != "double"){
    return(stop("Dependent variable is not a number."))
  }
  na.fail(a)
  return(a)
}
#Check if independent is a number or list of numbers and if it has missing values
checkIndependent <- function(a){
  if(class(a) == "list"){
    for(i in 1:length(a)){
      if(class(a[[i]]) != "numeric" & class(a[[i]]) != "integer"){
        return(stop("One or more independent variables is not a number."))
      }
      na.fail(a[[i]])
    }
  }else{
    if(class(a) != "numeric" & class(a) != "integer"){
      return(stop("One or more independent variables is not a number."))
    }
    na.fail(a)
  }
}
checkHypothesis <- function(ind,col,val){
  k <- ncol(ind)
  if(class(col) != "numeric" & class(col) != "integer"){
    return(stop("hypothesisColumn must be an integer"))
  }else if(class(val) != "numeric" & class(val) != "integer"){
    return(stop("hypothesisValue must be an number"))
  }else if(col>k |col <=0){
    return(stop("hypothesisColumn is not a recognized column of the matrix of independents."))
  }
}

checkRescale <- function(a){
  if(a != TRUE & a != FALSE){
    return(stop("Rescale input must be logical."))
  }
}
