adjust_estimator <- function(y, a, w, delta = NULL, method = "unadjust"){
  # y is the outcome, and must be a numerical vector. Here we don't deal with
  # categorial datatype;
  # a is the treatment arm, and must be a binary vector;
  # w is the baseline varible matrix or dataframe.
  # method == DR-WLS-U means we deal with missing data while using DR-WLS.
  # delta is the indicator vector of missing outcomes.
  
  if(!is.data.frame(w)){
    w <- data.frame(w)
  }  
  if(method == "unadjust"){
    est <- mean(y[a == 1]) - mean(y[a == 0])
  }else if(method == "IPW"){
    propensity <- glm(a~., family = binomial(), data = cbind(a,w))
    g <- predict(propensity, type = "response")
    est <- sum(y * a / g)/sum(a / g) - 
      sum(y * (1-a) / (1-g) )/sum((1-a) / (1-g))
  }else if(method == "DR-WLS"){
    propensity <- glm(arm~., family = binomial(), data = cbind(arm = a,w))
    g <- predict(propensity, type = "response")
    arm1 <- lm(y~., weights = 1/g[a == 1], 
                data = cbind(y= y[a == 1], w[a == 1, ]))
    arm0 <- lm(y~., weights = 1/(1 - g[a == 0]), 
                data = cbind(y = y[a == 0], w[a == 0, ]))
    est <- mean(predict(arm1, w) - predict(arm0, w))
  }else if(method == "DR-WLS-U"){
    propensity <- glm(arm~., family = binomial(), data = cbind(arm = a,w))
    g <- predict(propensity, type = "response")
    prob_missing <- glm(delta~., family = binomial(), data = cbind(delta, a, w))
    pm <- predict(prob_missing, type = "response")
    indi1 <- (a == 1) & (!is.na(y))
    indi2 <- (a == 0) & (!is.na(y))
    arm1 <- lm(y~., weights = 1/g[indi1]/pm[indi1], 
               data = cbind(y= y[indi1], w[indi1,]))
    arm0 <- lm(y~., weights = 1/(1 - g[indi2])/pm[indi2], 
               data = cbind(y = y[indi2], w[indi2,]))    
    est <- mean(predict(arm1, w) - predict(arm0, w))
  }else{
    stop("Undefined method")
  }
  return(est)
}
