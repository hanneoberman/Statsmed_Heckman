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

