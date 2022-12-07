#rm(list = ls())
library(lme4)
library(data.table)
library(mvtnorm)
library(fMultivar)
library(data.table)
library(mgcv)
library(GJRM)
library(Matrix)
library(mixmeta)
library(mvtnorm)
library(pbivnorm)
library(mice)
library(micemd)
library(miceMNAR)
library(glmmTMB)


#1. Data generation process ----

# 1.1. Define list of parameters ----
par_data <- list ( N = 10, # Number of studies
                   ni = 1000, # Number of patients per study
                   sigmav = 1, # Mean of sigma_e parameter (error terms)
                   sdsigma = 0.05, # Variance of the sigma_e parameter (error terms)
                   het = 1, # Homogeneous(=0) or heterogeneous (=1) variation on error terms
                   rho = 0.6, # Correlation of errors of the outcome and selection equation
                   vrho = 0.0, # Variance of rho parameter
                   vint =  0.4, # Variance of random intercept ( correlation o-s 0.3*rho)
                   vslp1 = 0.4, # Variance of random slope X1 ( correlation o-s 0.5*rho)
                   vslp2 = 0.4, # Variance of random slope X2 ( correlation o-s 0.5*rho)
                   vslp3 = 0.2, # Variance of random slope X3 ( correlation o-s 0.5*rho)
                   skew = 0, # 1= if error terms are skewed t distributed, 0= Normal distributed
                   tau = 0.4, # correlation of random effects of outcome and selection (same for random intercept and slope)
                   syst = 1, # 1= 20% of data are systematically missing otherwise all studies are sporadically missing
                   typey = "normal", # Type continuous variable"probit"
                   betao = c(0.3,1,1), # Beta output model(continuous)
                   betas = c(-0.8,1.3,-0.7,1.2),
                   typem = 0, # Type missing generation 0 = Heckman type 1= truncated
                   valadd=0
)


# Parameters for continuous missing variable (systematic)
par_data_cont1 <- par_data
par_data_cont1$rho= 0.9

par_data_cont2 <- par_data
par_data_cont2$rho= 0.6

par_data_cont3 <- par_data
par_data_cont3$rho= 0.3

par_data_cont4 <- par_data
par_data_cont4$rho= 0

# Parameters size variation (systematic)
par_data_conts1 <- par_data
par_data_conts1$N=100
par_data_conts1$ni=1000
par_data_conts1$rho = 0.6

par_data_conts2 <- par_data
par_data_conts2$N=50
par_data_conts2$ni=1000
par_data_conts2$rho = 0.6

par_data_conts3 <- par_data
par_data_conts3$N=10
par_data_conts3$ni=100
par_data_conts3$rho = 0.6

par_data_conts4 <- par_data
par_data_conts4$N=10
par_data_conts4$ni=50
par_data_conts4$rho = 0.6


# Parameters binomial (systematic)
par_data_bin <- par_data
par_data_bin$typey<-"binomial"

par_data_bin1 <- par_data_bin
par_data_bin1$rho= 0.9

par_data_bin2 <- par_data_bin
par_data_bin2$rho= 0.6

par_data_bin3 <- par_data_bin
par_data_bin3$rho= 0.3

par_data_bin4 <- par_data_bin
par_data_bin4$rho= 0.0

# Parameters of deviations in distribution assumptions
# Parameters for continuous missing variable MNAR truncated (systematic)
par_data_contrun <- par_data
par_data_contrun$typey <- "normal"
par_data_contrun$rho <- 0
par_data_contrun$het <- 1
par_data_contrun$typem<- 1
par_data_contrun$valadd<- -0.96

# Parameters for skewed distribution (systematic)
par_data_skew <- par_data
par_data_skew$skew <- 1
par_data_skew$valadd <- 0.29

###

#par_datal=list(par_data_cont1,par_data_cont4,par_data_bin1,par_data_bin4)
#par_datal=list(par_data_conts1,par_data_conts2,par_data_conts3,par_data_conts4)
#par_datal=list(par_data_bin1,par_data_bin2,par_data_bin3,par_data_bin4)
#par_datal=list(par_data_contrun,par_data_skew)
#par_datal=list(par_data_cont1,par_data_cont2,par_data_cont3,par_data_cont4)


# 1.2. Generate data with missing outcome ---- 
gen_complete_data <- function( DataID, #ID data
                               par_data){ #list of parameters to generate data
  
  set.seed(DataID*3+56)
  N <- par_data$N # Number of studies
  ni <- par_data$ni # Number of patients per study
  sigmav <- par_data$sigmav # Mean of sigma_e parameter (error terms)
  sdsigma <- par_data$sdsigma # Variance of the sigma_e parameter (error terms)
  het <- par_data$het # Homogeneous(=0) or heterogeneous (=1) variation on error terms
  rho <- par_data$rho # Correlation of errors of the outcome and selection equation
  vrho <- par_data$vrho # Variance of rho parameter
  vint <- par_data$vint # Variance of random intercept ( correlation o-s 0.3*rho)
  vslp1 <- par_data$vslp1 # Variance of random slope X1 ( correlation o-s 0.5*rho)
  vslp2 <- par_data$vslp2 # Variance of random slope X1 ( correlation o-s 0.5*rho)
  vslp3 <- par_data$vslp3 # Variance of random slope X1 ( correlation o-s 0.5*rho)
  skew <- par_data$skew # 1= if error terms are skewed t distributed, 0= Normal distributed
  tau <- par_data$tau # correlation of random effects of outcome and selection (same for random intercept and slope)
  syst <- par_data$syst # 1= 20% of data are systematically missing otherwise all studies are sporadically missing
  typey <- par_data$typey # "normal", "probit"
  typem <- par_data$typem # 0= Heckman, 1= Truncated
  betao0 <- par_data$betao[[1]] # Coefficient beta0 output model
  betao1 <- par_data$betao[[2]] # Coefficient beta1 output model
  betao2 <- par_data$betao[[3]] # Coefficient beta2 output model
  betas0 <- par_data$betas[[1]] # Coefficient beta0 selection model
  betas1 <- par_data$betas[[2]] # Coefficient beta1 selection model
  betas2 <- par_data$betas[[3]] # Coefficient beta2 selection model
  betas3 <- par_data$betas[[4]] # Coefficient beta3 selection model
  valadd <- par_data$valadd # Coefficient beta3 selection model
  
  # No correlation between random intercept and random slope at outcome and selection equations
  
  inv.logit <- function(p){
    return(exp(p)/(1 + exp(p)))
  }
  
  # Generate the covariates
  
  d <- matrix(0.0015, 2, 2)
  d[1,1] <- d[2,2] <-0.2
  Mu <- as.vector(mvtnorm::rmvnorm(n = 1, mean = c(0,0), sigma = d))
  
  x1 <- NULL
  x2 <- NULL
  x3 <- NULL
  
  for (i in 1:N){
    x1i <- rbinom(ni,size=1,prob=0.6)
    x2i <- rnorm(n = ni, mean=Mu[[1]], sd = 1)
    x3i <- rnorm(n = ni, mean=Mu[[2]], sd = 0.7)
    x1 <- c(x1,x1i)
    x2 <- c(x2,x2i)
    x3 <- c(x3,x3i)
  }
  

  
  # Generate the error term
  if ( het == 0){ #Sigma_e homogeneous
    d <- diag(2)
    d[2,2] <- sigmav^2
    d[2,1] <- d[1,2] <- rho*(sigmav)
    eps <- mvtnorm::rmvnorm(n = ni* N, mean = rep(0,2), sigma = d)
    
    if(skew == 1){
      #rmvsnorm
      eps <- fMultivar::rmvst(n = ni*N, dim = 2, mu = c(0,0), Omega = d, alpha = c(-5,-5), df = 8)
    }
    
  }else{ #Sigma_e heterogeneous
    
    sigmae <- exp( rnorm( n= N, mean=log(sigmav), sd = sdsigma))
    eps <- NULL
    
    for( i in 1: N){
      d <- diag(2)
      d[2,2] <- sigmae[i]^2
      d[2,1] <- d[1,2] <- sqrt(d[1,1])*sqrt(d[2,2])*rho
      eps_i <- mvtnorm::rmvnorm(n = ni, mean = rep(0,2), sigma = d)
      
      if( skew == 1){
        eps_i <- fMultivar::rmvst(n = ni, dim = 2, mu = c(0,0), Omega = d, alpha = c(2,-6), df = 4)
      }
      eps <- rbind(eps,eps_i)}
  }
  
  eps_s <- eps[,1] #error selection
  eps_o <- eps[,2] #error outcome
  
  # Generate random intercept
  d_re <- matrix(1,2,2)
  d_re[2,1] <- d_re[1,2] <- tau*rho
  d_rei <- vint*d_re
  alphai <- rmvnorm(n = N, mean = rep(0,2), sigma = d_rei)
  alphai_s <- rep(alphai[,1], each=ni) #random intercept selection
  alphai_o <- rep(alphai[,2], each=ni) #random intercept outcome
  
  # Generate random slope x1
  d_res1 <- vslp1*d_re
  alphas1 <- rmvnorm(n = N, mean = rep(0,2), sigma = d_res1)
  alphas1_s <- rep(alphas1[,1], each=ni) #random intercept selection
  alphas1_o <- rep(alphas1[,2], each=ni) #random intercept outcome
  
  # Generate random slope x1
  d_res2 <- vslp2*d_re
  alphas2 <- rmvnorm(n = N, mean = rep(0,2), sigma = d_res2)
  alphas2_s <- rep(alphas2[,1], each=ni) #random intercept selection
  alphas2_o <- rep(alphas2[,2], each=ni) #random intercept outcome
  
  # Generate random slope x3
  alpha3 <- rnorm(n = N, mean=0, sd = vslp3^0.5)
  alphas3_s <- rep(alpha3, each=ni) #random intercept selection
  
  
  # Generate latent outcome and selection variables
  
  Y <- (betao0+alphai_o) + (betao1+alphas1_o)*x1 + (betao2+alphas2_o)*x2 + eps_o
  R <-  valadd+(betas0+alphai_s) + (betas1+alphas1_s)*x1 + (betas2+alphas2_s)*x2 + (betas3+alphas3_s)*x3 +eps_s
  
  # Generate observed outcome and missing indicator variables
  ry <- ifelse(R > 0, 1, 0)
  
  if( par_data$typey == "binomial"){ # If Y is binary
    Y <- ifelse( Y > 0,1,0)
    Y <- factor(Y)}
  
  if ( typem == 1){ # MNAR f(Y and X1)
    R <-  valadd + Y
    ry <- rbinom(n=ni* N,size=1, prob= inv.logit(R))}
  
  
  # Generate missing in Y
  y <- Y
  y[ ry==0 ] <- NA
  
  
  group <- rep(1:N, each=ni)
  data <- data.table(Y,R,y,ry,x1,x2,x3,group)
  
  # Apply systematically missingness to the 20% of studies
  
  if( par_data$syst == 1){
    
    msys <- ceiling(N*0.2)
    data[group%in%c(1:msys), y := NA]
    data[group%in%c(1:msys), ry:= 0]
  }
  
  
  return(as.data.frame(data))
}


#2. Imputation process ----

# 2.1 Define settings for imputation methods-----

Sigma.co <- matrix(0, 5, 5)
diag(Sigma.co) <- exp(rnorm(5))
data0 <- as.data.frame(GJRM::rMVN(100, rep(0,5), Sigma.co))
colnames(data0) <- c('y','x1','x2','x3','group')
ini <- mice(data0, maxit = 0)

#only Y missing
pred_H <- ini$pred
pred_H[,"group"] <- 0
pred_H["group",] <- 0
pred_H["y","group"]<--2
pred_H["y","x3"] <- -3

pred_R<-pred_H
pred_R["y","x3"] <- 1 #-1 at the beginning

pred_list <- list('H'=pred_H,
                  'R'=pred_R,
                  'RB'=pred_R)

meth_list <- list('H'=c("2l.heckman","","","",""),
                  'R'=c("2l.2stage.norm","","","",""),
                  'RB'=c("2l.2stage.bin","","","",""))
dm<-cbind(name=c('HF','R'),
          meth=c('H','R'),
          full=c(TRUE,FALSE))

dm<-data.table(dm)
par_imp1<-list(dm=dm,pred=pred_list,meth=meth_list)


# 2.2. 2l.Heckman imputation method ----
source("/Users/jmunozav/Desktop/Statsmed_Heckman/4.Codes/mice.impute.2l.heckman.R")

#2.3. Get estimates from imputed datasets ----
get_estimates <- function( data = data, # data with missing values
                           DataID = DataID, #data ID
                           ImpID = ImpID, # imputation ID
                           par_imp = par_imp, #parameters imputation
                           M = M){ # number of imputations
  
  datam <- as.data.frame(data[,c("y","x1","x2","x3","group")])
  if (class(data$y) == "numeric"){
    family<-"normal"}else{
      family<-"probit"
    }
  
  
  if (ImpID == 0){ #Complete Case analysis
    
    # Evaluate the model with imputed data
    if (class(data$y) == "numeric"){
      p_model_imp <- with(data, lme4::lmer(y~x1+x2+(1|group)+(0+x1|group)+(0+x2|group)))
    } else if (class(data$y) == "factor" & nlevels(data$y) == 2){
      p_model_imp <- with(data, lme4::glmer(y~x1+x2+(1|group)+(0+x1|group)+(0+x2|group),family=binomial("probit")))
    }
    
    start_time <- Sys.time()
    end_time <- Sys.time()
    Qbar <- p_model_imp@beta # Parameter estimate
    Ubar <- diag(vcov(p_model_imp)) # Variance
    B <- rep(NA,3) #Between imputation variance
    T <- Ubar #Total variance
    v <- VarCorr(p_model_imp) #Random effect Variance
    Re <- c(rint = attr(v$group,"stddev"), rslp1=attr(v$group.1, "stddev"), rslp2=attr(v$group.2, "stddev"), resi = attr(v, "sc"))
    CI <- data.table(confint(p_model_imp,method="Wald"))
    colnames(CI) <- c("LCI","UCI")
    CI <- CI[!is.na(LCI),]
    LCI <- CI[,1]
    UCI <- CI[,2]
    m <- 1
    
  }else {
    if(ImpID != 3){ #HeckmanIPD or 2l.2stage.norm
      
      name <- par_imp$dm[[ImpID,"name"]]
      methv <- ifelse(par_imp$dm[[ImpID,"meth"]] == "R"&family == "probit",
                      'RB',par_imp$dm[[ImpID,"meth"]])
      fullv <- par_imp$dm[[ImpID,"full"]]
      
      
      # Calculates criteria on each imputation set
      start_time <- Sys.time()
      data_imp <- mice( datam, # dataset with missing values
                        m = M, # number of imputations
                        seed = ImpID*34+24*DataID, #seed attached to the dataID
                        meth = par_imp$meth[[methv]], #imputation method vector
                        pred = par_imp$pred[[methv]], #imputation predictors matrix
                        full = as.logical(fullv), #full or separated by parameter set
                        family= family,
                        meta_method = "reml",
                        print = FALSE,
                        maxit = 1)
      
      end_time <- Sys.time()
      
    }else if(ImpID == 3){ # Heckman method & Galimard
      
      datag <- data[,c("y","x1","x2","x3","group")]
      datag$group <- as.factor(datag$group)
      if(family=="probit"){
        datag$y <- as.factor(datag$y)
      }
      JMEy <- generate_JointModelEq(data=datag,varMNAR = "y")
      # This applies for systematically missing pattern as the function does not allow to include the group term
      JMEy[,"y_var_sel"] <- c(0,1,1,1,0)
      JMEy[,"y_var_out"] <- c(0,1,1,0,0)
      arg <- MNARargument(data=datag,varMNAR="y",JointModelEq=JMEy)
      
      start_time <- Sys.time()
      data_imp <- mice(data = arg$data_mod,
                       seed = ImpID*34+24*DataID, #seed attached to the dataID
                       method = arg$method,
                       predictorMatrix = arg$predictorMatrix,
                       JointModelEq = arg$JointModelEq,
                       control = arg$control,
                       m = M,
                       maxit = 1)
      end_time <- Sys.time()
    }
    
    
    # Evaluate the model with imputed data
    if (class(data$y) == "numeric"){
      p_model_imp <- try(with(data_imp,lme4::lmer(y~x1+x2+(1|group)+(0+x1|group)+(0+x2|group))))

    }else if (class(data$y) == "factor" & nlevels(data$y) == 2){
      p_model_imp <- try(with(data_imp,lme4::glmer(y~x1+x2+(1|group)+(0+x1|group)+(0+x2|group),
                                                   family=binomial("probit"))), silent=TRUE)
      if(inherits( p_model_imp, "try-error")){
        p_model_imp<-NULL
        p_model_imp$glmmBB=1
        for (i in 1:10){
          p_model_imp[[i]] <- try(glmmTMB:: glmmTMB(y~x1+x2+(1|group)+(0+x1|group)+(0+x2|group),
                                                    family=binomial("probit"),data=complete(data_imp, i)), silent=TRUE)
        }
        
      }
    }
    
    # Get estimates for each imputation
    get_imputation_estimates <- function(impN, p_model_imp){
      if (is.null(p_model_imp$glmmBB)){
        betas <- p_model_imp$analyses[[impN]]@beta
        var <- diag(vcov(p_model_imp$analyses[[impN]]))
        v <- VarCorr(p_model_imp$analyses[[impN]])
        res <- c(rint=attr(v$group,"stddev"),rslp1=attr(v$group.1, "stddev"),rslp2=attr(v$group.2, "stddev"),resi=attr(v, "sc"))
      }else{
        betas <- summary(p_model_imp[[impN]])$coefficients$cond[,"Estimate"]
        var <- diag(vcov(p_model_imp[[impN]])$cond)
        v <- VarCorr(p_model_imp[[impN]])$cond
        res <- c(rint=attr(v$group,"stddev"),rslp1=attr(v$group.1, "stddev"),rslp2=attr(v$group.2, "stddev"),resi=attr(v, "sc"))
      }
      return(list(betas,var,res))
    }
    
    imp_est <- lapply(1:M,FUN = get_imputation_estimates, p_model_imp = p_model_imp)
    imp_est <- do.call(Map, c(f = rbind, imp_est))
    
    
    #Pool confidence interval according Rubin's rule for betas taken from mice package
    
    Beta <- imp_est[[1]]# Parameter estimate
    Var <- imp_est[[2]] # Variance
    Re <- colMeans(imp_est[[3]]) #Random effects parameters (int,slp,res) in std.
    
    m <- nrow(Beta) # Number of imputed datasets
    Qbar <- colMeans(Beta) # Pooled parameter estimate
    Ubar <- colMeans(Var) # Average within imputation variance
    B <- colSums((Beta - matrix(rep(Qbar,each=m),nrow=m))^2)/(m-1) #Between imputation variance
    T <- Ubar+B+B/m # Total variance
    rm <- (B+B/m)/Ubar # Relative increase in variance that is due to the missing values.
    v <- (m-1)*(1+1/rm)^2 # Reference t distribution degrees of freedom
    #95% Confidence interval
    LCI <- Qbar-qt(0.975,v)*sqrt(T)
    UCI <- Qbar+qt(0.975,v)*sqrt(T)
    
  }
  
  p <- length(Qbar)-1
  outnam <- c("DataID","ImpID", paste0("Qbar_",0:p), paste0("Ubar_",0:p), paste0("Re_",0:(length(Re)-1)), paste0("B_",0:p), paste0("T_",0:p), paste0("LCI_",0:p), paste0("UCI_",0:p), "m", "time", "pmis")
  output <- cbind(DataID = DataID, ImpID = ImpID, t(Qbar), t(Ubar), t(Re), t(B), t(T), t(LCI), t(UCI), m, time = as.numeric(difftime(end_time, start_time, units = "secs")), pmis = sum(is.na(data$y))/nrow(data)*100)
  colnames(output) <- outnam
  return(output)
}

#3. HPC Function ----
fun_HPC <- function(ID,par_data,par_imp,M){
  DataID <- ceiling(ID/(nrow(par_imp$dm)+2))
  ImpID <- ID%%(nrow(par_imp$dm)+2)
  data <- gen_complete_data( DataID, #ID data
                             par_data) #list of parameters
  final_est <- try(get_estimates(data = data, # data with missing values
                                 DataID = DataID, #data ID
                                 ImpID = ImpID, # imputation ID
                                 par_imp = par_imp, #parameters imputation
                                 M = M), silent=TRUE)
  if(inherits( final_est, "try-error")){
    final_est <- cbind(DataID = DataID,ImpID=ImpID, run = 0) }
  else{
    final_est<-cbind(final_est,run = 1)}
  
  return(final_est)
}



