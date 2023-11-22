rm(list=ls())
require(mice)
require(micemd)
library(tidyverse)
library(broom.mixed)
library(nlme)
library(ggplot2)
library(ggmice)
library(patchwork)

## Load the dataset ----
# data("Obesity",pkg="micemd")

## Explore the dataset ----

### Missing patterns ----
md.pattern(Obesity)
by(Obesity,INDICES = Obesity$Cluster,summary)

myplots <- lapply(1:5, function(i) {
  ggmice::plot_pattern(setDT(Obesity)[Cluster==i])+ 
    ggtitle(paste0("Cluster", i))+ theme(legend.position = ifelse(i==1,"bottom","none"))
})

myplots[[1]]+ myplots[[2]] /myplots[[3]]+ myplots[[4]] /myplots[[5]]


### BMI vs covariates ----
ggplot(Obesity[complete.cases(Obesity),],aes(x=Age,y=BMI,group=Cluster))+
  geom_point(color="black",size=0.1)+
  geom_smooth(aes(color=as.factor(Cluster)),method = lm)+
  facet_grid(.~Cluster)+
  theme_classic()+
  scale_color_viridis_d()+
  theme(legend.position = "none")

ggplot(Obesity[complete.cases(Obesity),],aes(x=Age,y=Weight,group=Cluster))+
  geom_point(color="black",size=0.1)+
  geom_smooth(aes(color=as.factor(Cluster)))+
  facet_grid(.~Cluster)+
  theme_classic()+
  scale_color_viridis_d()+
  theme(legend.position = "none")

ggplot(Obesity[complete.cases(Obesity),], aes(x = FamOb, y = BMI))+ 
  geom_boxplot(aes(fill = FamOb), show.legend = FALSE) +
  facet_grid(.~Cluster)+
  theme_classic()+
  scale_fill_viridis_d()

ggplot(Obesity[complete.cases(Obesity),], aes(x = Gender, y = BMI))+ 
  geom_boxplot(aes(fill = Gender), show.legend = FALSE) +
  facet_grid(.~Cluster)+
  theme_classic()+
  scale_fill_viridis_d()


### Check ICC ----
Nulmodel=lme4::lmer(BMI ~ 1 + (1|Cluster), data = Obesity)
performance::icc(Nulmodel)


## Density plots of covariates ----

ggplot(Obesity, aes(x=Weight, color=as.factor(Cluster), fill=as.factor(Cluster))) +
  geom_density(alpha=0.3)+
  facet_grid(~Cluster)+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_classic()+theme(legend.position = "none")

ggplot(Obesity, aes(x=Age, color=as.factor(Cluster), fill=as.factor(Cluster))) +
  geom_density(alpha=0.3)+
  facet_grid(~Cluster)+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_classic()+theme(legend.position = "none")

ggplot(Obesity, aes(x=Height, color=as.factor(Cluster), fill=as.factor(Cluster))) +
  geom_density(alpha=0.3)+
  facet_grid(~Cluster)+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_classic()+theme(legend.position = "none")



## MAR ----

library(micemd)
meth_mar <- find.defaultMethod(Obesity, ind.clust=1, I.small = 7, ni.small = 100, prop.small = 0.4)
meth_mar["BMI"]<- "~ I(Weight / (Height)^2)"
meth_mar["Age"]<-"2l.2stage.pmm" 

pred_mar <- mice(Obesity, maxit = 0)$pred 
pred_mar[,"Time"]<-0
pred_mar[,"Cluster"]<--2
pred_mar[pred_mar==1]<-2
pred_mar[c("Height", "Weight"), "BMI"] <- 0

### MAR normal ----
imp_mar <- mice(data = Obesity, meth = meth_mar, pred = pred_mar, m=10, seed = 123, pred_std=FALSE)
summary(complete(imp_mar,"long"))
plot(imp_mar)

### MAR pmm ----
meth_mar["Weight"]<-"2l.2stage.pmm" 
imp_mar_pmm <- mice(data = Obesity, meth = meth_mar, pred = pred_mar, m=10, seed = 123, pred_std=FALSE)
summary(complete(imp_mar_pmm,"long"))
plot(imp_mar_pmm)


ggmice(imp_mar, ggplot2::aes(x = Weight, group = .imp)) +
  ggplot2::geom_density() 
ggmice(imp_mar, ggplot2::aes(x = Weight, group = .imp)) + 
  ggplot2::geom_density() + 
  ggplot2::scale_size_manual(values = c("observed" = 1, "imputed" = 0.5),guide = "none") +
  ggplot2::theme(legend.position = "none")+
  ggplot2::facet_wrap(~Cluster)


## MNAR ----

### Check exclusion restriction variable ----
dObesity <- Obesity %>% mutate_if(is.factor, as.numeric)
dObesity$ry<-ifelse(is.na(Obesity$Weight),1,0)
outm <- lm(Weight ~ Gender + Height + FamOb + Age + Time,data=dObesity)
summary(outm) # Time is not asociated with weight
selm <- glm(ry ~ Gender + Height + FamOb + Age+Time,family=binomial(link = "logit"),data=dObesity)
summary(selm) # Time is not asociated with R


### MNAR normal ----
meth_mnar <- meth_mar
meth_mnar["Weight"]<- "2l.2stage.heckman"
pred_mnar<-pred_mar
pred_mnar["Weight","Time"]  <- -3

imp_mnar<- mice(data = Obesity, meth = meth_mnar, pred = pred_mnar, m=10, seed = 123, pred_std=F)
summary(complete(imp_mnar,"long"))
plot(imp_mnar)

ggmice(imp_mnar, ggplot2::aes(x = Weight, group = .imp)) +
  ggplot2::geom_density() 
ggmice(imp_mnar, ggplot2::aes(x = Weight, group = .imp)) + 
  ggplot2::geom_density() + 
  ggplot2::scale_size_manual(values = c("observed" = 1, "imputed" = 0.5),guide = "none") +
  ggplot2::theme(legend.position = "none")+
  ggplot2::facet_wrap(~Cluster)


### MNAR pmm ----
imp_mnarpmm <- mice(data = Obesity, meth = meth_mnar, pred = pred_mnar, m=10, seed = 123, pmm=T, pred_std=F)
summary(complete(imp_mnarpmm,"long"))
plot(imp_mnarpmm)


ggmice(imp_mnarpmm, ggplot2::aes(x = Weight, group = .imp)) +
  ggplot2::geom_density() 
ggmice(imp_mnarpmm, ggplot2::aes(x = Weight, group = .imp)) + 
  ggplot2::geom_density() + 
  ggplot2::scale_size_manual(values = c("observed" = 1, "imputed" = 0.5),guide = "none") +
  ggplot2::theme(legend.position = "none")+
  ggplot2::facet_wrap(~Cluster)



# Fit models ----
###  Random intercept model ----
cc_ri <- with(Obesity[complete.cases(Obesity),],lme( BMI ~ Age + FamOb + Gender,random=~1|Cluster))
mar_ri <- with(imp_mar,lme( BMI ~ Age + FamOb + Gender,random=~1|Cluster))
mar_pmm_ri <- with(imp_marpmm,lme( BMI ~ Age + FamOb + Gender,random=~1|Cluster))
mnar_ri<- with(imp_mnar,lme(BMI ~ Age + FamOb + Gender,random=~1|Cluster))
mnar_pmm_ri<- with(imp_mnarpmm,lme(BMI ~ Age + FamOb + Gender,random=~1|Cluster))
list_models<-list(cc_ri,mar_ri,mnar_ri,mnar_pmm_ri)
plot_models(list_models,mod_name=c("Complete case", "MAR", "MNAR","MNAR_pmm"))

### Random slope model ----

cc_rs<- with(Obesity[complete.cases(Obesity),],lme( BMI ~ Age + FamOb + Gender,random=~1+Age|Cluster))
mar_rs <- with(imp_mar,lme( BMI ~ Age + FamOb + Gender,random=~1+Age|Cluster))
mar_pmm_rs <- with(imp_mar_pmm,lme( BMI ~ Age + FamOb + Gender,random=~1+Age|Cluster))
mnar_rs<- with(imp_mnar,lme(BMI ~ Age + FamOb + Gender,random=~1+Age|Cluster))
mnar_pmm_rs<- with(imp_mnarpmm,lme(BMI ~ Age + FamOb + Gender,random=~1+Age|Cluster))
list_models<-list(cc_rs,mar_rs,mar_pmm_rs,mnar_rs,mnar_pmm_rs)

### Plot coefficients ----

plot_models <- function(list_models, mod_name){
  mod_table<-NULL
  for(i in 1:length(list_models)){
    mod <- list_models[[i]]  
    if(class(mod)[[1]]=="mira"){
      mod_tab <-summary(pool(mod), conf.int = TRUE)[,c("term","estimate","2.5 %","97.5 %","p.value")]
    }else{
      sum_cc <- (summary(mod,conf.int = TRUE))
      mod_tab<-data.frame(intervals(mod, which = "fixed")$fixed)
      mod_tab$term<-rownames(mod_tab)
      mod_tab$pvalue<-sum_cc$tTable[,5]
      mod_tab <-mod_tab[,c("term","est.","lower","upper","pvalue")]
      rownames(mod_tab)<-NULL
      
    }
    colnames(mod_tab)<-c("Coefficient","Estimate","lower","upper","p.value")
    mod_tab$model<-mod_name[i]
    mod_table<-rbind(mod_table,mod_tab)
  }
  
  plot <-ggplot(mod_table, aes(x = Coefficient, y = Estimate, ymin = lower, ymax = upper, group=model, color=model)) + 
    geom_hline(yintercept=0,linetype="dashed",color="black",size=0.1)+
    geom_linerange( position = position_dodge(width = 0.8)) + 
    geom_pointrange(size=0.1, position = position_dodge(width = 0.8))+coord_flip()+
    theme_classic()+theme(legend.position = "bottom")+scale_color_viridis_d()
  return(plot)
}

plot_models(list_models,mod_name=c("Complete case", "MAR", "MAR_pmm", "MNAR","MNAR_pmm"))

<<<<<<< Updated upstream
||||||| Stash base
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
=======
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
library(ggplot2)
library(data.table)
library(lme4)
library(mice)
library(micemd)
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


source("./4.Codes/mice.impute.2l.2stage.heckman.R")

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
>>>>>>> Stashed changes
