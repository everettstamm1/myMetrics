#' Wild Bootstrap Hypothesis Testing
#'
#' This function takes a vector of a dependent variable, a vector of a single independent variable or a list of
#' several independent variables, a number of bootstrap repetitions, a hypothesis, and a option to rescale the residuals
#' and performs restricted residual resampling. To learn more about the theory, please read section 3.4 of "Bootstrap Methods
#' in Econometrics" by James A. MacKinnon, Queen’s Economics Department Working Paper No. 1028, 2006.
#'
#'
#' @param dependents numeric vector of dependent variable
#' @param independents matrix of independent variable(s), do not include constant term.
#' @param repetitions number of bootstrap repetitions to be performed
#' @param hypothesisColumn column index in 'independents' of which coefficient restriction will be tested
#' @param hypothesisValue hypothesized value of the coefficient of hypothesisColumn, numeric
#' @param rv Distribution of random variable (with mean 0 and variance 1) to multiply the residuals by. Acceptable types are 'mammen', 'rademacher', and 'standard normal'.
#' @param transform Whether or not to transform the t-th residual by 1/(1-ht)^.5, where ht is the t-th diagonal of the hat matrix `X((X^t*X)^-1)*X^t`.
#' This is typically done when some observations have high leverage.
#' @param suppressConstant logical. If true, there will be no constant term added to the regression. Defaults to false.
#' @return A character string detailing the resulting p-value of the test.
#' @export
bootWild <- function(dependent,independents,repetitions,hypothesisColumn,hypothesisValue, rv, transform = TRUE,suppressConstant = FALSE){
  #Making sure all the inputs are cool and if dependent or independent are vectors, making the matrices
  dependent <- checkDependent(dependent)
  hypothesisColumn <- checkSuppress(suppressConstant,hypothesisColumn)
  independents <- checkIndependent(independents,suppressConstant)
  checkHypothesis(independents,hypothesisColumn,hypothesisValue)
  if(transform != TRUE & transform != FALSE)
    return(stop("transform must be of type logical"))
  checkReps(repetitions)
  if(nrow(dependent) != nrow(independents))
    return(stop("dependent and independent matrices have different number of observations."))
  if(ncol(independents) == 1 & suppressConstant){

  }

  #Premultiplying (Xt*X)^-1 and ((Xt*X)^-1)*Xt so it doesn't need to be done at every recursion of the helper functions
  invXtX <- solve(t(independents) %*% independents)
  invXtX_tX <- invXtX %*% t(independents)

  #Creating base sample estimates
  n <- nrow(independents)
  k <- ncol(independents)
  b_ols <- invXtX_tX %*% dependent
  u_ols <- dependent - independents %*% b_ols
  se_ols <- sqrt(diag(c((t(u_ols) %*% u_ols)/(n-k)) * invXtX))
  t_ols <- (b_ols[hypothesisColumn]-hypothesisValue)/se_ols[hypothesisColumn]


  #Creating bootstrapped dependent values
  Y_star <- setupWild(dependent,independents,hypothesisColumn,hypothesisValue, repetitions,rv,invXtX_tX, transform)



  #Creating t-values for each bootstrap iteration
  t_vals <- runWild(Y_star,independents,hypothesisColumn,hypothesisValue,invXtX,invXtX_tX,t_vec = c())


  #Creating Bootstrap P Values
  p <- boot_p(t_vals,t_ols)

  print(paste0("The P-value is ",p," for ",repetitions," bootstrap repetitions of the null hypothesis against the two tailed alternative."))
  return(p)
}

#HELPER FUNCTIONS

#creates n drawings from the distirbution (with mean zero and variance 1) specified by "rv" in the main function.
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

#Creates a matrix of iterations of bootstrapped dependent values. Number of iterations specified by 'repetitions.'
setupWild <- function(dep,indep,hypInd,hypVal,Repetitions,rv,invXtX_tX,transform){
  n <- length(dep)
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

  v <- randomVariableMatrix(rv, n,Repetitions)
  u_star <- v*c(u_res)
  if(transform == TRUE){
    h <- diag(indep %*%invXtX_tX)
    u_star <- u_star/(1-h)^.5
  }
  if(ncol(indep) >1){
    Y_star <- c(hypVal*indep[,hypInd] + X_res %*% b_res) + u_star
    return(Y_star)
  }else{
    Y_star <- c(hypVal*indep[,hypInd]) + u_star
    return(Y_star)
  }

}

#Creates t-values for all iterations of the restriction test.
runWild <- function(Y_stars,indep,index,h,invXtX,invXtX_tX,t_vec){
  if(length(Y_stars) == nrow(indep)){
    b_star <- invXtX_tX %*% Y_stars
    u_star <- Y_stars - indep %*% b_star
    se_star <- sqrt(diag(c((t(u_star) %*% u_star)/(nrow(indep)-ncol(indep))) * invXtX))
    t_star <- (b_star[index]-h)/se_star[index]
    t_vec <- c(t_vec,t_star)
    return(t_vec)
  }else{
    b_star <- invXtX_tX %*% Y_stars[,1]
    u_star <- Y_stars[,1] - indep %*% b_star
    se_star <- sqrt(diag(c((t(u_star) %*% u_star)/(nrow(indep)-ncol(indep))) * invXtX))
    t_star <- (b_star[index]-h)/se_star[index]
    t_vec <- c(t_vec,t_star)
    Y_stars <- Y_stars[,-c(1)]
    return(runWild(Y_stars,indep,index,h,invXtX,invXtX_tX,t_vec))
  }
}

#Calculates p-value of bootstrapped hypothesis test.
boot_p <- function(t_boot, t_base){
  B <- length(t_boot)
  F_val <- (1/B)*sum(t_base > t_boot)
  return(2*min(1-F_val,F_val))
}

checkIndependent <- function(a,s){
  if(is.vector(a)){
    a <- as.matrix(a)
  }
  if(!is.matrix(a)){
    return(stop("Independent variable must be a matrix!"))
  }else if(typeof(a)!= "numeric" & typeof(a) != "integer"&typeof(a) != "double"){
    return(stop("One or more independent variables is not a number."))
  }else if(s == FALSE){
    n <- nrow(a)
    k <- ncol(a)
    a <- matrix(c(rep(1,n),a),ncol = k+1)
  }
  na.fail(a)
  return(a)
}

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

checkReps <- function(a){
  if(class(a) != "numeric" & class(a) != "integer"){
    return(stop("Repetitions must be a number."))
  }else if(a <=0){
    return(stop("Repetitions must be greater than zero"))
  }
}

checkLengths <- function(dep,ind){
  vec <- c()
  vec <- c(vec,length(dep))
  if(typeof(ind) == "list"){
    for(i in 1:length(ind)){
      vec <- c(vec,length(ind[[i]]))
    }
  }else{
    vec <- c(vec,length(ind))
  }
  if(length(unique(vec)) >1){
    return(stop("Differing number of observations in variables"))
  }
}
