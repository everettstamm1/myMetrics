#' Non-Parametric Restricted Residual Resampling Bootstrap
#'
#' This function takes a vector of a dependent variable, a vector of a single independent variable or a list of
#' several independent variables, a number of bootstrap repetitions, a hypothesis, and a option to rescale the residuals
#' and performs restricted residual resampling. To learn more about the theory, please read section 3.1 of "Bootstrap Methods
#' in Econometrics" by James A. MacKinnon, Queen’s Economics Department Working Paper No. 1028, 2006.
#'
#'
#' @param Dependents numeric vector of dependent variable
#' @param Independents numeric vector or list of numeric vectors of independent variable(s)
#' @param Repetitions number of bootstrap repetitions to be performed
#' @param Hypothesis Linear restriction to be tested, must be of the format "Xi = h" where the i-th vector in the Independent
#' list of vectors is the variable you wish to test a restriction on and h is the restriction. Note: if you have multiple digits of independent variables,
#' you must enter the hypothesis with all digits, i.e. if you have 100 independent variables and you wish to test a restriction on
#' the 3-rd vector in the Independent list of vectors, your hypothesis should be written "X003 = h", rather than "X3 = h".
#' @param rescale Whether or not to rescale the residuals by (n/(n-k))^.5. This is typically done to ensure residuals have correct variance.

#' @return A character string detailing the resulting p-value of the test.
#' @export
bootNonParametric <- function(Dependent,Independents,Repetitions,Hypothesis, rescale = FALSE){
  #Making sure all the inputs are cool
  checkDependent(Dependent)
  checkIndependent(Independents)
  deps <- length(Independents)
  dig <-nchar(as.character(deps))
  checkHypothesis(Hypothesis,deps,dig)
  checkRescale(rescale)
  checkReps(Repetitions)
  checkLengths(Dependent,Independents)


  #Creating data frame
  df <- data.frame(matrix(unlist(Independents), ncol = deps))
  colnames(df) <- str_c("X",1:deps)
  df$Y <- Dependent

  #Estimating the base restriction
  mod <- reg(df)

  #Creating the original T Statistic
  xr <- as.numeric(str_sub(Hypothesis,2,1+dig))
  est_b <- mod$coefficients[(1+xr)]
  sd_b <- summary(mod)$coefficients[(1+xr),2]
  h <- as.numeric(str_sub(Hypothesis,5+dig))
  t_base<- ((est_b-h)/sd_b)


  #Creating Bootstrap T Statistics
  booty <- boot_rest(df,Hypothesis, Repetitions,dig,rescale, c())
  #Get it? Bootstrap T Statistics -> Bootstrap T -> Boot T -> Booty

  #Creating Bootstrap P Values
  peas <- boot_p(booty,t_base)

  output_message <- stringr::str_c("The hypothesis ",Hypothesis," has a p-value of ",peas," after ", Repetitions," repetitions of resampling of restricted residuals ")
  if(rescale == TRUE){
    output_message <- stringr::str_c(output_message, "with rescaling factor (n/(n-k))^.5.")
  }else{
    output_message <- stringr::str_c(output_message, "without rescaling factor (n/(n-k))^.5.")
  }
  return(output_message)

}
boot_p <- function(t_boot, t_base){
  B <- length(t_boot)
  F_val <- (1/B)*sum(t_base > t_boot)
  return(2*min(1-F_val,F_val))
}
boot_rest <- function(d_f_1,hyp, R,digits,rescale, t_vec){
  if(R == 1){
      dat <- boot_restrict_rescale(d_f_1,digits,hyp,rescale)
      mod <- reg(dat)
      xr <- as.numeric(stringr::str_sub(hyp,2,1+digits))
      est_b <- mod$coefficients[(1+xr)]
      sd_b <- summary(mod)$coefficients[(1+xr),2]
      h <- as.numeric(stringr::str_sub(hyp,5+digits))
      t_star <- ((est_b-h)/sd_b)
      t_vec <- c(t_vec,t_star)
      return(t_vec)
  }else{
    dat <- boot_restrict_rescale(d_f_1,digits,hyp,rescale)
    mod <- reg(dat)
    xr <- as.numeric(stringr::str_sub(hyp,2,1+digits))
    est_b <- mod$coefficients[(1+xr)]
    sd_b <- summary(mod)$coefficients[(1+xr),2]
    h <- as.numeric(stringr::str_sub(hyp,5+digits))
    t_star <- ((est_b-h)/sd_b)
    t_vec <- c(t_vec,t_star)
    return(boot_rest(dat,hyp,R-1,digits,rescale,t_vec)) #That's correct, I know how to write recursive functions
  }
}


boot_restrict_rescale <- function(data_f,digits, hyp, rescale){
  mod <- reg(data_f)
  n <- length(mod$residuals)
  k <- ncol(data_f)-1
  cons <- mod$coefficients[1]
  B <- mod$coefficients[2:(k+1)]
  h <- as.numeric(stringr::str_sub(hyp,5+digits))
  #creating restricted residual
  u <- data_f$Y - cons
  for(i in 1:k){
    if(i == as.numeric(stringr::str_sub(hyp,2,1+digits))){
      u <- u - h*data_f[[i]]
    }else{
      u <- u - B[i]*data_f[[i]]
    }
  }
  #rescaling the residual
  if(rescale == TRUE){
    u_scale <- ((n/(n-k))^.5)*u
  }
  #resampling the residual
  u_scale_sample <- sample(u_scale, size = n, replace = TRUE)
  #Creating synthetic y
  synthetic_y <- cons + u_scale_sample
  for(i in 1:k){
    if(i == as.numeric(stringr::str_sub(hyp,2,1+digits))){
      synthetic_y<- synthetic_y + h*data_f[[i]]
    }else{
      synthetic_y<- synthetic_y + B[i]*data_f[[i]]
    }
  }
  synth_df <- data_f
  synth_df$Y <- synthetic_y
  return(synth_df)
}

#OLS regression from data frame of form [Y, X1, X2, ..., Xn]
reg <- function(df){
  n <- ncol(df)-1
  formula <- str_c("Y ~ ",str_c("X",1:n,collapse = "+"))
  return(lm(data = df,formula))
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

checkReps <- function(a){
  if(class(a) != "numeric" & class(a) != "integer"){
    return(stop("Repetitions must be a number."))
  }else if(a <=0){
    return(stop("Repetitions must be greater than zero"))
  }
}

#Check if dependent is a number and if it has missing values
checkDependent <- function(a){
  if(class(a) != "numeric" & class(a) != "integer"){
    return(stop("Dependent variable is not a number."))
  }
  na.fail(a)

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

#Check if hypothesis input is correct, n is number of independent variables
checkHypothesis <- function(hyp, n, digits){
  if(class(hyp) != "character"){
    return(stop("Hypothesis is not a character."))
  }else if((stringr::str_sub(hyp,0,1) !="x" & stringr::str_sub(hyp,0,1) !="X") | as.numeric(stringr::str_sub(hyp,2,1+digits)) > n){
      return(stop("Hypothesis is formatted incorrectly. Must be in format \"Xi = h\", where i is the location of the dependent variable you wish to test and h is the hypothesised value. Must include spaces."))
  }else if(is.na(as.numeric(stringr::str_sub(hyp,5+digits)))){
      return(stop("Hypothesis is formatted incorrectly. Must be in format \"Xi = h\", where i is the location of the dependent variable you wish to test and h is the hypothesised value. Must include spaces."))
  }
}

checkRescale <- function(a){
  if(typeof(a) != "logical"){
    return(stop("Rescale input must be logical."))
  }
}
