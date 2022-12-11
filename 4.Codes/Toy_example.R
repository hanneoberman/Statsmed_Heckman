rm(list=ls())
library(data.table)
library(truncnorm)
library(nlme)
library(devtools)

#0. Load original dataset ----

data <- as.data.table(read.csv('/Users/jmunozav/Desktop/Statsmed_Heckman/4.Codes/ObesityDataSet_raw_and_data_sinthetic.csv', stringsAsFactors=FALSE, fileEncoding="latin1"))

#1. Modify dataset in order to get a variable weight with MNAR missing mechanism ----
# Generate countries where survey was taken
data <- data[,c("Weight","Gender","Age","Height","FAVC")] # Select some of the variables of the original dataset
data[,Cluster:=c(rep(1,500),rep(2,400),rep(3,450),rep(4,400),rep(5,361))] # in the imputation method is necessary cluster is included as numeric
N<-length(unique(data$Cluster))
#  Generation of the exclusion restriction (ERV): The time variable describes the time it takes for a person to answer the weight question.
set.seed(12345)
data[,Time:=rtruncnorm(n=nrow(data), a = 1, b = 10, mean = 5, sd = 3)] 
data[,Gender:=(ifelse(Gender=="Female",1,0))]
data[,FAVC:=(ifelse(FAVC=="yes",1,0))]
data[,Heights:=scale(Height)]
data[,Ages:=scale(Age)]
data[,Times:=scale(Time)]

MAR_model <- with(data, lme4::lmer(Weight~Gender+Ages+Heights+FAVC+(1+Heights+Gender+Ages|Cluster)))
summary(MAR_model)
# Generate outcome variable
data[,XOBO:=predict(MAR_model)+15]

# Generate random slope for Gender Age on selection equation
alpha_gender <- rnorm(n = N, mean=0, sd = 0.01)
alpha_age <- rnorm(n = N, mean=0, sd = 0.01)
nobs<-c(500,400,450,400,361)
data[,Genderi:=unlist(mapply(rep,alpha_gender,each=nobs))]
data[,Agei:=unlist(mapply(rep,alpha_age,each=nobs))]
data[,XSBS:=1.5-(0.7+Genderi)*Gender-(0.5+Agei)*Ages-1.2*Times] # Simulated selection model

# Simulate bivariate normal correlated errors
v <- VarCorr(MAR_model) #Get varcov matrix
eps <- NULL
rho=-0.8 # We assume that non-responders are more likely to have a high weight.
set.seed(12345)
sigmae <- exp( rnorm( n=N, mean=log(attr(v, "sc")), sd = 0.02))
for( i in 1: N){
  d <- diag(2)
  d[2,2] <- sigmae[i]^2
  d[2,1] <- d[1,2] <- sqrt(d[1,1])*sqrt(d[2,2])*rho
  eps_i <- mvtnorm::rmvnorm(n = nrow(data[Cluster==i]), mean = rep(0,2), sigma = d)
  eps <- rbind(eps,eps_i)
}

# Generate latent outcome and selection variables
data[,Weight.star:=XOBO+eps[,2]]
data[,ry.star:=XSBS+eps[,1]]
data[,ry:=ifelse(ry.star>0,1,0)]

# Generate observed variables
data[,ry:=ifelse(ry.star>0&Cluster!=5,1,0)] #Include systematic missingness on cluster 4
data[,Weight :=ifelse(ry==1,Weight.star,NA)] 

#Generate missing in other variables
data_other<-data[,c("Gender","Age","Height","FAVC","Weight.star")]
# ampute the complete data once for every mechanism
ampdata0 <- ampute(data_other, patterns = c(1,1,0,0,1), prop = 0.2, mech = "MAR")$amp
ampdata1 <- ampute(data_other, patterns = c(1,0,0,0,1), prop = 0.1, mech = "MAR")$amp
# create a random allocation vector
# use the prob argument to specify how much of each mechanism should be created
# here, 0.5 of the missingness should be MAR and 0.5 should be MCAR
indices <- sample(x = c(1, 2), size = nrow(ampdata0), replace = TRUE, prob = c(0.7, 0.3))
ampdata <- matrix(NA, nrow = nrow(data_other), ncol = ncol(data_other))
ampdata[indices == 1, ] <- as.matrix(ampdata0[indices == 1, ])
ampdata[indices == 2, ] <- as.matrix(ampdata1[indices == 2, ])
colnames(ampdata)<-colnames(data_other)
ampdata<-data.table(ampdata)
ampdata[,Weight:=data$Weight]
ampdata[,Time:=data$Time]
ampdata[,Cluster:=data$Cluster]
ampdata[,Weight.star:=NULL]
ampdata[,FAVC :=as.factor(FAVC)]
ampdata[,Gender :=as.factor(Gender)]
dataObs<-as.data.frame(ampdata)
save(dataObs,file="dataObs.Rdata")


#2. Descriptive analysis simulated dataset ----
# Description: The dataset used here was based on data collected from 2111 individuals
#in an online obesity survey in different locations. 
# The data was simplified and grouped into five clusters. The values and observability 
#of the weight variable were defined according to the Heckman model in a herarchical model.
# In addition, a systematic loss of this variable at one location was assumed. 
# For the other variables included in the data set, missing values were simulated according to a MAR mechanism. 

# Usage
#data(dataObs)

# Format
# A dataframe with 2111 observations with the following variables:
# Gender a factor value with two levels: 1 ("Female"), 0 ("Male").
# Age a numeric value indicating age of subject in years.
# Height a numeric value with Height in meters.
# FAVC a factor value describing the frequent consumption of high caloric food (FAVC) with two levels:
# 1("Yes"), 0("Male").
# Weight a numeric value with Weight in Kilograms.
# Time a numeric value indicating time in responding the questions in minutes.
# Cluster a numeric indexing the cluster. 

# Details
# For more details about the original dataset. 
# Palechor, F. M., & de la Hoz Manotas, A. (2019). Dataset for estimation of obesity levels based on eating habits and physical condition in individuals from Colombia, Peru and Mexico. Data in brief, 25, 104344.
# Simulation data code on Toyexample github repository.

# Source 
# Dataset obtained from "https://www.kaggle.com/datasets/fabinmndez/obesitydata?select=ObesityDataSet_raw_and_data_sinthetic.csv"
# and it is described here https://www.sciencedirect.com/science/article/pii/S2352340919306985?via%3Dihub

# Example
rm(list=ls())
library(ggmice)
library(data.table)
load("dataObs.Rdata")
#data(dataObs)
summary(dataObs)
md.pattern(dataObs)

# Count missingness per group
dataNA<-setDT(dataObs)[, .(nNA = sum(is.na(Weight)),n=.N), by = Cluster]
dataNA[, propNA:=nNA/n]
dataNA

#Plot weight
ggplot(dataObs, aes(x = Weight, group=as.factor(Cluster))) +
geom_histogram(aes(color = as.factor(Cluster),fill= as.factor(Cluster)),
               position = "identity", bins = 30)+facet_grid(Cluster~.)
  

#3. Imputation model
# Imputation considering weight is MAR

seed=123

ini <- mice(dataObs, maxit = 0)
meth <- ini$method
pred <-ini$predictorMatrix
meth[c("Age","Height","Weight")] <- "2l.2stage.norm" 
meth[c("FAVC")] <- "2l.2stage.bin" 
pred[,"Cluster"] <- -2 # Cluster variable
pred[,"Time"] <- 0 # Cluster variable

imp0 <- mice(data = dataObs, meth = meth, pred = pred, seed = seed)
ggmice(imp0, ggplot2::aes(x = Weight, group = .imp)) +
  ggplot2::geom_density() 
ggmice(imp0, ggplot2::aes(x = Weight, group = .imp)) + 
  ggplot2::geom_density() + 
  ggplot2::scale_size_manual(values = c("observed" = 1, "imputed" = 0.5),guide = "none") +
  ggplot2::theme(legend.position = "none")+
  facet_wrap(~Cluster)

(summary(complete(imp0,"long")$Weight))


source("/Users/jmunozav/Desktop/Statsmed_Heckman/4.Codes/mice.impute.2l.2stage.heckman.R")

# Data passed into the mice function has to specify the cluster variable as a numeric variable, and in case that incomplete
# variable is a binary variable, it should specified as factor variable with 2 levels. 

#It is necessary to specify for the MNAR missing variable the method "2l.heckman" in the mice methods vector.     
meth["Weight"] <- "2l.2stage.heckman" 
#Furthermore, in the prediction matrix, the group or cluster variable should be specified as "-2", all predictor variables
#belonging to the selection and outcome as "1", the exclusion restrictions or predictor variables that are only included in the selection equation 
#as "-3" and those that are only included in the outcome equation as "-4".

pred["Weight","Time"]  <- -3 # Variable only affects the selection model (Exclusion restriction)
pred["Weight",c("Height","FAVC")]  <- -4 # Variables only affect the outcome model
# As Age and Gender are predictors in both selection and outcome model are set as 1
head(dataObs)


# Imputation of weight variable without pmm 
imp <- mice(data = dataObs, meth = meth, pred = pred, seed = seed)
plot(imp)
ggmice(imp, ggplot2::aes(x = Weight, group = .imp)) +
  ggplot2::geom_density() 
ggmice(imp, ggplot2::aes(x = Weight, group = .imp)) + 
  ggplot2::geom_density() + 
  ggplot2::scale_size_manual(values = c("observed" = 1, "imputed" = 0.5),guide = "none") +
  ggplot2::theme(legend.position = "none")+
  facet_wrap(~Cluster)


# Imputation of weight variable with (predictive mean matching) pmm: we can also use the ppm approach, by setting pmm= TRUE and providing a vector of donnors ypmm,
# for instance in this case the vector or donnors are weight from 35 to 180 kilograms and the imputed values are given in this range of values.

imp_pmm <- mice(data = dataObs, meth = meth, pred = pred, seed = seed, pmm=TRUE,ypmm=seq(35,180,0.1))
plot(imp_pmm) #Check convergence
ggmice(imp_pmm, ggplot2::aes(x = Weight, group = .imp)) +
  ggplot2::geom_density() 
ggmice(imp_pmm, ggplot2::aes(x = Weight, group = .imp)) + 
  ggplot2::geom_density() + 
  ggplot2::scale_size_manual(values = c("observed" = 1, "imputed" = 0.5),guide = "none") +
  ggplot2::theme(legend.position = "none")+
  facet_wrap(~Cluster)


# Summary three imputation methods
(summary(complete(imp0,"long")$Weight))
(summary(complete(imp,"long")$Weight))
(summary(complete(imp_pmm,"long")$Weight))

# Model
library(broom.mixed)
model_MAR <- with(imp0, lmer(Weight~Gender+Age+Height+FAVC+(1|Cluster)))
model_MNAR <- with(imp, lmer(Weight~Gender+Age+Height+FAVC+(1|Cluster)))
model_MNAR_pmm <- with(imp_pmm, lmer(Weight~Gender+Age+Height+FAVC+(1|Cluster)))
summary(pool(model_MAR))
summary(pool(model_MNAR))
summary(pool(model_MNAR_pmm))
