# Synthetic data based on a obesity survey retrieved from
# https://www.kaggle.com/datasets/fabinmndez/obesitydata/

rm(list=ls())

library(data.table)
library(truncnorm)
library(nlme)
library(devtools)
library(mice)
library(mvtnorm)


generate_data <- function(betao,betas){
  
  set.seed(1234)
  N <- 5
  Ns <- c(500,400,450,400,361)
  Gender <- Age <- Height <- BMI <- FamOb <- Time <- Cluster <- bo1 <-bo2<- bs1 <-bs2 <- NULL
  sigmae <- exp( rnorm( n= N, mean=log(1.8), sd = 0.3))
  es <- eo <- NULL
  # Effects
  
  # Random effects
  # Generate random intercept and slope outcome
  d <-diag(2)
  rho= -0.4
  d[1,1] <- 15^2
  d[2,2] <- 6^2
  d[2,1] <- d[1,2] <- rho*sqrt(d[1,1]*d[2,2])
  bo <- rmvnorm(n = N, mean = rep(0,2), sigma = d)
  # Generate random intercept and slope selection
  d <-diag(2)
  rho= 0.4
  d[1,1] <- 0.8^2
  d[2,2] <- 0.05^2
  d[2,1] <- d[1,2] <- rho*sqrt(d[1,1]*d[2,2])
  bs <- rmvnorm(n = N, mean = rep(0,2), sigma = d)
  
  for (i in 1:N){
    
    n <- Ns[i]
    Gender <- c( Gender, rbinom(n = n, size = 1, prob = 0.49))
    Age <- c( Age, round(rtruncnorm(n, a = 14, b = 61, mean = 24, sd = 6.34), 0)) # rounded to integers
    Height <- c( Height, rtruncnorm(n, a = 1.45, b = 1.98, mean = 1.7, sd = 0.09))
    BMI <- c(BMI,rtruncnorm(n, a = 13, b = 50, mean = 29, sd = 8))
    FamOb <- c(FamOb,rbinom(n = n, size = 1, prob = 0.85))
    Time<- c(Time,rtruncnorm(n=n, a = 1, b = 10, mean = 5, sd = 2.5))
    Cluster <- c(Cluster,rep(i,n))
    bo1 <- c(bo1,rep(bo[i,1],n))
    bo2 <- c(bo2,rep(bo[i,2],n))
    bs1 <- c(bs1,rep(bs[i,1],n))
    bs2 <- c(bs2,rep(bs[i,2],n))
    # Error terms eps1<- selection eps2<- outcome
    d <- diag(2)
    rho=  -0.9
    d[2,2] <- sigmae[i]^2
    d[2,1] <- d[1,2] <- sqrt(d[1,1])*sqrt(d[2,2])*rho
    ei <- mvtnorm::rmvnorm(n = n, mean = rep(0,2), sigma = d)
    es <- c(es,ei[,1])
    eo <- c(eo,ei[,2])
  }  
  
  ds <- data.table(Cluster, Gender, Age, Height, BMI, FamOb,Time)
  ds2 <- apply(ds, MARGIN = 2, FUN = function(X) (if(diff(range(X))==0){1} else{(X - min(X))/diff(range(X))})) # scaled
  ds2 <- as.data.table(ds2)
  
  ds[, y.star := (betao[1]+bo1)+(betao[2]+bo2)*ds2$Age+(betao[3])*as.numeric(ds2$Gender)+(betao[4])*ds2$Height+(betao[5])*as.numeric(ds2$FamOb)+(betao[6])*ds2$BMI + eo ]
  ds[, ry.star :=  (betas[1]+bs1)+(betas[2]+bs2)*ds2$Age+(betas[3])*as.numeric(ds2$Gender)+(betas[4])*ds2$Height+(betas[5])*as.numeric(ds2$FamOb)+(betas[6])*ds2$BMI+(betas[7])*ds2$Time + es]
  ds[, ry :=  ifelse( ry.star>1&Cluster!=3,1,0)]
  ds[, Weight :=  ifelse( ry==1, y.star,NA)]

  
  #Ampute the covariates
  for (i in 1:N){
    data_miss <- ds[is.na(Weight)&Cluster==i,c("Gender","Height","Age","FamOb")]
    # Ampute development dataset
    patterns = rbind(c(1,1,0,1),c(1,1,1,0),c(1,1,0,0),c(1,0,0,0)) # (0=missing, 1=observed)
    freq = c(0.3,0.3,0.25,0.15) # ocurrence of patterns
    mech = c("MAR") 
    weights = ampute.default.weights(patterns, mech) # Weights of weighted sum scores
    
    #md.pattern(data_inc)
    data_inc <-   ampute(data = data_miss,
                         patterns = patterns,
                         prop = 0.4, #Percentage of missing cells
                         freq = freq, 
                         mech = mech,
                         weights=weights)$amp
    ds[is.na(Weight)&Cluster==i, Height := data_inc$Height]
    ds[is.na(Weight)&Cluster==i, Age := data_inc$Age]
    ds[is.na(Weight)&Cluster==i, FamOb := data_inc$FamOb]
  }
  
  ds[, Gender:= factor(Gender,labels=c("Female","Male"))]
  ds[, FamOb:= factor(FamOb,labels=c("no","yes"))]
  ds[, BMI:=Weight/(Height)^2]

  
  return(ds)}

ds <- generate_data(betao = c(40.62, 10.5, -0.90, 20.87, 2.41, 80.30),  
                    betas =  c(3.5, 0.8, -0.8,   3.0, -0.7,-3.5,-2.1))
summary(ds)
Obesity <- ds[,c("Cluster","Time","Gender","Height","FamOb","Age","Weight","BMI")]
save(Obesity,file="Obesity.rda")


library(ggmice)
library(data.table)
data(Obesity)
summary(Obesity)
md.pattern(Obesity)
# Count missingness per group
dataNA<-setDT(Obesity)[, .(nNA = sum(is.na(Weight)),n=.N), by = Cluster]
dataNA[, propNA:=nNA/n]
dataNA
#Plot weight
Obesity$Cluster<-as.factor(Obesity$Cluster)
ggplot(Obesity, aes(x = Weight, group=Cluster)) +
  geom_histogram(aes(color = Cluster,fill= Cluster),
                 position = "identity", bins = 30) + facet_grid(Cluster~.)


ggplot(ds,aes(x=Age,y= y.star ,group=Cluster))+
  geom_point(color="black",size=0.1)+
  geom_smooth(aes(color=as.factor(Cluster)))+
  facet_grid(.~Cluster)+
  theme_classic()+
  scale_color_viridis_d()+
  theme(legend.position = "none")


