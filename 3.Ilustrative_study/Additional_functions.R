
# 1. Functions to pool estimates first graph ----

logit <- function(x) {
  log(x/(1-x))
}
inv.logit <- function(x) {
  1/(1+exp(-x))
}
inv.logit.SE<- function(logit.se, c) {
  logit.se * (c*(1-c))
}

#Function for calculating absolute risk
f.prevalence <- function(data,outcome_name) {
  a<-summary(data[,get(outcome_name)])
  b<-prop.test(x = a[[2]], n = sum(a), correct = TRUE)
  prop<-b$estimate[[1]]
  ci<-b$conf
  return(list(prevalence= prop,ci.l=ci[1], ci.u=ci[2]))
}


#Function to calculate absolute prevalence of an outcome per study and per imputed dataset
f.abs.perstudy<-function(data, outcome_name,study_name) {
  
  #Calculate absolute prevalence of outcome per study
  inc.outcome<-setDT(data)[, f.prevalence(data=.SD, outcome_name=outcome_name), by = list(get(study_name),.imp)]
  
  #Replace 0's by 0.000001
  inc.outcome$prevalence[inc.outcome$prevalence==0]<-0.000001
  inc.outcome$ci.l[inc.outcome$ci.l==0]<-0.000001
  
  #Logit transformation
  inc.outcome[,logit.prevalence:=logit(prevalence)]
  inc.outcome[,logit.ci.l:=logit(ci.l)]
  inc.outcome[,logit.ci.u:=logit(ci.u)]
  inc.outcome[,logit.se:=(logit.ci.u-logit.ci.l)/(2*1.96)]
  inc.outcome
}

#Function to pool the absolute prevalence per study per imputed dataset using rubins rules
#Results in absolute risk per study

f.pool<-function(data,m){

  alpha <- 0.05
  tr_est <- data$logit.prevalence
  tr_se  <- data$logit.se
  n <- length(tr_est)  #number observations
  mu <- mean(tr_est, na.rm = T) # pool estimate
  w_var <- mean(tr_se^2, na.rm = T) # within variance
  b_var <- var(tr_est, na.rm = T) # between variance
  t_var <- sum(w_var,b_var,b_var/n,na.rm=T) # total variance
  t_se <- sqrt(t_var) # total standard error
  r <-  sum(b_var,(b_var / n),na.rm=T)/ w_var # relative increase variance due to missing values
  v <- (n - 1) * (1 + r^-1)^2 # degrees of freedom
  t <- qt(1-alpha/2, v) #t criteria
  if (is.infinite(v)|is.na(v)){ # t can not be calculated e.g 1 only observation it is aprox to normal 
    t <- qnorm(1-alpha/2)
  }
  
  prevalence <- inv.logit(mu)*100 # mean(est, na.rm = T)
  ci.lb <- inv.logit(mu - t_se*t)*100
  ci.ub <- inv.logit(mu + t_se*t)*100
  return(list(prevalence=prevalence,ci.lb=ci.lb,ci.ub=ci.ub))
}


f.abs.poolrubin <-function(data,study_name) {
  abs.outcome<-setDT(data)[, f.pool(data=.SD, m<-max(data$.imp)), by = list(get)]
  setnames(abs.outcome, "get", study_name)
  return(abs.outcome)
}


# 2. Functions to pool estimates second graph ----
pool_coef <- function(x){
  colMeans(t(sapply(x$analyses, coef)))
}

pool_variance <- function(x){
  m <- length(x$analyses)
  within <- Reduce(`+`, lapply(x$analyses, vcov)) / m
  between <- var(t(sapply(x$analyses, coef)))
  total_var <- within + (1 + 1 / m) * between
}

#Prediction imputed data 
pred<-function(model,source){
  mm   <- as.matrix(model.matrix(model$analyses[[1]]))
  coef <- as.matrix(pool_coef(model),nrow=1)
  mod  <- data.table(cbind(subreg_id=as.character(hdataF0f$subreg_id),age=as.numeric(hdataF0f$age)))
  mod[,subreg_id := as.factor(subreg_id)]
  mod[,age := as.numeric(age)]
  mod[,fit := mm%*%coef]
  mod[,se.fit := sqrt(diag(mm %*% pool_variance(model) %*% t(mm)))]
  mod[,pred := fn(fit)*100]
  mod[,lb:=fn(fit - 1.96 * se.fit)*100]
  mod[,ub:=fn(fit + 1.96 * se.fit)*100]
  mod[,source:=source]
}


