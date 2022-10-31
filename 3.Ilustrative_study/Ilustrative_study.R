library(here) 
library(data.table)
library(micemd)
library(mice)
library(gamm4)
library(tidymv)
library(ggplot2)



#0. Additional required function ----
# 0.1 Functions to pool estimates first graph ----
logit <- function(x) {log(x/(1-x))}
inv.logit <- function(x) {1/(1+exp(-x))}
inv.logit.SE<- function(logit.se, c) {logit.se * (c*(1-c))}

#Function for calculating absolute risk
f.prevalence <- function(data,outcome_name) {
  a<-summary(data[,get(outcome_name)])
  b<-prop.test(x = a[[2]], n = sum(a), correct = TRUE)
  prop<-b$estimate[[1]]
  ci<-b$conf
  return(list(prevalence= prop,ci.l=ci[1], ci.u=ci[2]))
}


#Function to calculate absolute risk of an outcome per study and per imputed dataset
f.abs.perstudy<-function(data, outcome_name,study_name) {
  
  #Calculate absolute risk of outcome per study
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

#Function to pool the absolute risks per study per imputed dataset using rubins rules
#Results in absolute risk per study

f.pool<-function(data,m){
  logit.abs <- mean(data$logit.prevalence)
  within <- mean(data$logit.se^2)
  between <- (1 + (1/m)) * var(data$logit.prevalence)
  logit.var <- within + between
  logit.se <- sqrt(logit.var)
  logit.ci.lb <- logit.abs + qnorm(0.05/2)     * logit.se
  logit.ci.ub <- logit.abs + qnorm(1 - 0.05/2) * logit.se
  prevalence<-inv.logit(logit.abs)*100
  ci.lb<-inv.logit(logit.ci.lb)*100
  ci.ub<-inv.logit(logit.ci.ub)*100
  return(list(prevalence=prevalence,ci.lb=ci.lb,ci.ub=ci.ub))
}


f.abs.poolrubin <-function(data,study_name) {
  abs.outcome<-setDT(data)[, f.pool(data=.SD, m<-max(data$.imp)), by = list(get)]
  setnames(abs.outcome, "get", study_name)
  return(abs.outcome)
}


# 0.2. Functions to pool estimates second graph ----
pool_coef <- function(x){
  colMeans(t(sapply(x$analyses, coef)))
}

pool_variance <- function(x){
  m <- length(x$analyses)
  within <- Reduce(`+`, lapply(x$analyses, vcov)) / m
  between <- var(t(sapply(x$analyses, coef)))
  total_var <- within + (1 + 1 / m) * between
}

#Prediction imputated data 
pred_model<-function(model,source){
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




# 1. Subset original dataset ----

data_original <- read.delim(here('3.Ilustrative_study','ISASimple_Gates_LLINE-UP_rct_RSRC.txt')) # Dataset retrieved from https://clinepidb.org/ce/app/workspace/analyses/DS_7c4cd6bba9/new/details#AccessRequest
dist <- readxl::read_xlsx(here('3.Ilustrative_study','Districts.xlsx')) #district specification
source(here('3.Ilustrative_study','Additional_functions.R'))
source(here('4.Codes','mice.impute.2l.heckman.R'))


data <- data.table()
data[, par_id := data_original$Participant_Id ]
data[, dis_id := data_original$District..EUPATH_0000407. ]
data[, clu_id :=data_original$Cluster.ID..EUPATH_0035024.]
data[, p_anopheles := data_original$Log10.cluster.level.mean.female.Anopheles..EUPATH_0044173. ] #prevalence anopheles
data[, test := data_original$Plasmodium..by.thick.smear.microscopy..EUPATH_0024314. ] # result plasmodium test 
data[, age := data_original$Age..OBI_0001169.]
data[, sex := data_original$Sex..PATO_0000047. ]
data[, bednets_pp := data_original$X1.bednet.per.2.people..EUPATH_0020219. ]
data[, wealthin_hou := data_original$Household.wealth.index..numerical..EUPATH_0000014. ]
data[, date_obs := data_original$Observation.date..EUPATH_0004991. ]
data[, study_time := data_original$Study.timepoint..OBI_0001508. ]
data[, consent := data_original$Consent.for.lab.testing..EUPATH_0044111.]

hdata <- merge(data,dist,by="dis_id") #merge name of subregions
hdata <- setDT(hdata)[consent!="No"&age>=2&age<=10,] #restrict to children aged 2-10 years old, with consent for lab testing
hdata[,test:=ifelse(test=="Yes",1,ifelse(test=="No",0,NA))] # Test value
hdata[,ry:=ifelse(is.na(test),0,1)]

#2. Exclusion restriction variable----
hdata[,date_obs:=as.Date(date_obs,"%d-%m-%Y")]
hdata[,day_obs:= weekdays(date_obs)]
hdata[,weekend:=ifelse(day_obs%in%c("Saturday","Sunday"),"yes","no")]
#Information taken from Uganda school schedules for 2017
hdata[,schoolterm:=ifelse(date_obs>="2017-02-06"&date_obs<="2017-05-05","T117",
                   ifelse(date_obs>="2017-05-29"&date_obs<="2017-08-25","T217",
                   ifelse(date_obs>="2017-09-18"&date_obs<="2017-12-08","T317",
                   ifelse(date_obs>="2018-02-05"&date_obs<="2018-05-04","T118",
                   ifelse(date_obs>="2018-05-28"&date_obs<="2018-08-24","T218",
                   ifelse(date_obs>="2018-09-17"&date_obs<="2018-12-07","T318",
                   ifelse(date_obs>="2019-02-04"&date_obs<="2019-05-03","T119",
                   ifelse(date_obs>="2019-05-27"&date_obs<="2019-08-23","T219",
                   ifelse(date_obs>="2019-09-16"&date_obs<="2019-12-29","T319",
                         "Relax")))))))))]  
pubholidays<-as.Date(c("2017-01-01","2017-01-26","2017-02-16","2017-03-08","2017-04-14","2017-04-17","2017-05-01","2017-06-03","2017-06-09","2017-06-26","2017-09-01","2017-10-09","2017-12-25","2017-12-26",
                     "2018-01-01","2018-01-26","2018-02-16","2018-03-08","2018-03-30","2018-04-02","2018-05-01","2018-06-03","2018-06-09","2018-06-15","2018-08-21","2018-10-09","2018-12-25","2018-12-26",
                     "2019-01-01","2019-01-26","2019-02-16","2019-03-08","2019-04-19","2019-04-22","2019-05-01","2019-06-03","2019-06-09","2019-06-05","2019-08-12","2019-10-09","2019-12-25","2019-12-26"))                                   
hdata[,pholiday:=ifelse(date_obs%in%pubholidays,"holiday","workday")] 
hdata[,fundays:=ifelse((schoolterm=="Relax"|pholiday=="holiday"|weekend=="yes"),"yes","no")]

hdata[,wealthin_hou:=scale(wealthin_hou)]
hdata[,sex:=as.factor(sex)]
hdata[,bednets_pp:=as.factor(bednets_pp)]
hdata[,fundays:=as.factor(fundays)]

hdataF0 <- hdata[study_time==0,] #Only data at baseline
hdataF0$dis_name<-hdataF0$dis_id
hdataF0$dis_id<-as.integer(as.factor(hdataF0$dis_id)) #convert as numeric for imputation

sum.table<-setDT(hdataF0)[,.(did=(length(unique(dis_name))),
                             N=.N,
                             age=round(mean(age),2),
                             p_anopheles= round(mean(exp(p_anopheles)),2),
                             LIp_an=round(exp(quantile(p_anopheles,probs=0.025,type=1)),1),
                             UIp_an=round(exp(quantile(p_anopheles,probs=0.975,type=1)),1),
                             wealthin_hou=round(mean(wealthin_hou),2),
                             LIw=round(quantile(wealthin_hou,probs=0.025,type=1),1),
                             UIw=round(quantile(wealthin_hou,probs=0.975,type=1),1),
                             bednets_pp=round(sum(ifelse(bednets_pp=="Yes",1,0))/.N*100,1),
                             sex=round(sum(ifelse(sex=="Male",1,0))/.N*100,1),
                             fundays=round(sum(ifelse(fundays=="yes",1,0))/.N*100,1),
                             ry=round((.N-sum(ry))/.N*100,1)),by=list(subreg_id)]


sum.table[,p_anopheles:=paste0(p_anopheles,"[",LIp_an,",",UIp_an,"]")]
sum.table[,wealthin_hou:=paste0(wealthin_hou,"[",LIw,",",UIw,"]")]

sum.table[,c("LIp_an","UIp_an","LIw","UIw"):=NULL]

# 2.1. Exclusion restriction assesment----
fits <- gam(ry ~ s(age,k=3)+p_anopheles+wealthin_hou+bednets_pp+fundays+sex, family=binomial,data = hdataF0)
summary(fits)
fito <- gam(test ~s(age,k=3)+p_anopheles+wealthin_hou+bednets_pp+fundays+sex, family=binomial,data = hdataF0)
summary(fito)

tfito<-as.data.table(rbind(summary(fito)$p.table,summary(fito)$s.table))
tfits<-as.data.table(rbind(summary(fits)$p.table,summary(fits)$s.table))
ps<-c("***","***","","***","***","","***")
po <-c("***","***","***","***","","","***")
names<-c("(Intercept)","Female Anopheline/house (Log10 mean)","Household wealth index","Bednet for 2 person-Yes","Sample taken in holidays-Yes","Sex-Male", "s(Age)")
esto<-paste0(round(tfito$Estimate,3),"(",round(tfito$`Std. Error`,3),")",po)
ests<-paste0(round(tfits$Estimate,3),"(",round(tfits$`Std. Error`,3),")",ps)
tablexc<-as.data.table(cbind(names,esto,ests))


#3. Multiple Imputation ---- 

#3.1. Transform dataset for imputation----
#Get covariates of age splines
modg<-gam(wealthin_hou~s(age,k=3),data=hdataF0)
model.matrix(modg)
sage<-data.table(model.matrix(modg)[,2:3])
colnames(sage)<-c("sage1","sage2")

# Subset data for imputation
hdataF0f<- setDT(hdataF0[,c("dis_id","dis_name","subreg_id","p_anopheles","test","age","wealthin_hou","bednets_pp","sex","fundays")])
hdataF0f<-cbind(hdataF0f,sage)
hdataF0f[,subreg_id:=as.factor(subreg_id)]
hdataF0f[,sex:=as.factor(sex)]
hdataF0f[,bednets_pp:=as.factor(bednets_pp)]
hdataF0f[,test:=as.factor(test)]
hdataF0f[,fundays:=as.factor(fundays)]


#3.2. Set prediction matrix and methods -----
ini <- mice(hdataF0f, maxit = 0)
pred <- ini$pred
pred[,"dis_id"] <- -2
pred[,"dis_name"] <- 0
pred[,"subreg_id"] <- 0
pred["age",] <- 0
pred[,"age"] <- 0
pred["test","fundays"] <- -3
pred["test","sage2"]<- 1


#3.3. Multiple imputation with heckman model----
meth<-ini$method
meth[c("test")]<-"4l.heckman"

# Full correlation of parameters
# data_heck <- mice(hdataF0f, # dataset with missing values
#                   m = 20,   # number of imputations
#                   seed = 1234, #seed attached to the dataID
#                   meth = meth, #imputation method vector
#                   pred = pred, #imputation predictors matrix
#                   maxit=1,
#                   meta_method="reml",
#                   pmm=FALSE)
# save(data_heck,file=(here('3.Ilustrative_study','MalariaMNAR.RData')))
load(file=here('3.Ilustrative_study','MalariaMNAR.RData'))

#3.4. Multiple imputation with  2l.2stage.bin method----

# meth[c("test")]<-"2l.2stage.bin" 
# pred["test","fundays"] <- 0
# data_mar <- mice( hdataF0f, # dataset with missing values
#                        m = 20,   # number of imputations
#                        seed = 1234, #seed attached to the dataID
#                        meth = meth, #imputation method vector
#                        pred = pred, #imputation predictors matrix
#                        maxit=1)
# save(data_mar,file=(here('3.Ilustrative_study','MalariaMAR.RData')))
load(file=here('3.Ilustrative_study','MalariaMAR.RData'))

# 4.1 Prevalence plot ----

# Prevalence per subregion----

hdataF0f[,pts_num:=ifelse(is.na(test),NA,ifelse(test=="1",1,0))]
data_MCAR<-hdataF0f[, .(count = .N, positives = sum(pts_num,na.rm=TRUE), NAs=sum(is.na(pts_num))), by = subreg_id]
data_MCAR[,Miss:=(NAs)/(count)*100]
data_MCAR[,p_MCAR:=positives/(count-NAs)*100]
data_MCAR[,LC_MCAR:=(positives)/(count)*100]
data_MCAR[,UC_MCAR:=(positives+NAs)/(count)*100]
data_MCAR[,positives:=NULL]
data_MCAR[,NAs:=NULL]

#Under MAR assumption, calculate the prevalence of the outcome in every study separate and in every imputed dataset
prev.outcome_mar<-f.abs.perstudy(data=complete(data_mar, "long"),outcome_name = "test",study_name="subreg_id")
#under MNAR assumption, pool the prevalence over the imputed datasets, resulting in absolute risks per study
abs.outcome_mar<-f.abs.poolrubin(data=prev.outcome_mar,study_name="subreg_id")
colnames(abs.outcome_mar)<-c("subreg_id","p_MAR","LC_MAR","UC_MAR")

prev.outcome_mnar<-f.abs.perstudy(data=complete(data_heck, "long"),outcome_name = "test",study_name="subreg_id")
#under MNAR assumption, pool the prevalence over the imputed datasets, resulting in absolute risks per study
abs.outcome_mnar<-f.abs.poolrubin(data=prev.outcome_mnar,study_name="subreg_id")
colnames(abs.outcome_mnar)<-c("subreg_id","p_MNAR","LC_MNAR","UC_MNAR")

data_prev<-merge(data_MCAR, abs.outcome_mar,by="subreg_id")
data_prev<-merge(data_prev,abs.outcome_mnar,by="subreg_id")

data_long<- setDT(melt(data_prev, 
                  measure.vars = c( "p_MCAR","p_MAR","p_MNAR",
                                    "LC_MCAR","LC_MAR","LC_MNAR",
                                    "UC_MCAR","UC_MAR","UC_MNAR"),
                  variable.name = "Statistic", value.name = "Value"))

data_long[, c("Statistic", "Parameter") := tstrsplit(Statistic, "_", fixed=TRUE)]
data_wide<- setDT(dcast(data_long, subreg_id + count+ Miss+Parameter ~ Statistic, value.var = "Value"))
data_wide[,Parameter:=factor(Parameter,levels=c("MCAR","MAR","MNAR"))]
levels(data_wide$Parameter)<-c("CC","2l.MAR","2l.Heckman")

pd <- position_dodge(width = 0.4)
plot_dist<-ggplot(data_wide, aes(x = subreg_id, y = p, colour = Parameter)) +
                  geom_pointrange(aes(ymax = UC, ymin = LC),position = pd) + 
                  geom_point(position = pd)+ theme_light()+
                  ylab("Parasite prevalence (%)")+xlab("Sub-region")+
                  scale_color_brewer(palette="Dark2")
                  


# 5.Prevalence per subregion and age ----

# Complete case dataset
hdataF0f<-as.data.table(hdataF0f)
m0 <- gam(test ~ subreg_id + s(age,k=3,by = subreg_id), family=binomial,data = hdataF0f)
mod0 <- tidymv::predict_gam(m0)

fn<- function(x){exp(x)/(1+exp(x))}
mod0$pred <- with(mod0, fn(fit)*100)
mod0$lb <- with(mod0, fn(fit - 1.96 * se.fit)*100)
mod0$ub <- with(mod0, fn(fit + 1.96 * se.fit)*100)
mod0$source<-"CC"

# MNAR imputed dataset no full
mod.mnarf <- with(data_heck, gam(test ~ subreg_id + s(age,k=3,by = subreg_id), family=binomial))
mod1<-pred_model(model=mod.mnarf,source="2l.Heckman")  


# MAR imputed dataset
mod.mar <- with(data_mar, gam(test ~ subreg_id + s(age,k=3,by = subreg_id), family=binomial))
mod2<-pred_model(model=mod.mar,source="2l.MAR")  

#Combine datasets   
mod<-rbind(mod0,mod1,mod2)
setDT(mod)[,source:=as.factor(source)]
mod[,source:=factor(source,levels=c("CC","2l.MAR","2l.Heckman"))]

#Prevalence plot
plot_sdist_age <- ggplot(mod, aes(x = age, y = pred, group = source)) + 
  geom_line(aes(colour=source))+ 
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  geom_ribbon(aes(ymin=lb, ymax=ub,fill =source), alpha = 0.3)+
  facet_wrap(~subreg_id)+
  scale_y_continuous(limits=c(0,100)) +
  scale_x_continuous(limits=c(2,10)) +
  xlab("Age (years)")+ylab("Parasite prevalence (%)")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom",
        legend.title = element_blank())


##### END CODE########
# Prevalence per subregion----
# 
# hdataF0f[,pts_num:=ifelse(is.na(test),NA,ifelse(test=="1",1,0))]
# data_MCAR<-hdataF0f[, .(subreg_id=subreg_id[[1]],count = .N, positives = sum(pts_num,na.rm=TRUE), NAs=sum(is.na(pts_num))), by = dis_id]
# data_MCAR[,Miss:=(NAs)/(count)*100]
# data_MCAR[,p_MCAR:=positives/(count-NAs)*100]
# data_MCAR[,LC_MCAR:=(positives)/(count)*100]
# data_MCAR[,UC_MCAR:=(positives+NAs)/(count)*100]
# data_MCAR[,positives:=NULL]
# data_MCAR[,NAs:=NULL]
# 
# #Under MAR assumption, calculate the prevalence of the outcome in every study separate and in every imputed dataset
# prev.outcome_mar<-f.abs.perstudy(data=complete(data_mar, "long"),outcome_name = "test",study_name="dis_id")
# #under MNAR assumption, pool the prevalence over the imputed datasets, resulting in absolute risks per study
# abs.outcome_mar<-f.abs.poolrubin(data=prev.outcome_mar,study_name="dis_id")
# colnames(abs.outcome_mar)<-c("dis_id","p_MAR","LC_MAR","UC_MAR")
# 
# prev.outcome_mnar<-f.abs.perstudy(data=complete(data_heck, "long"),outcome_name = "test",study_name="dis_id")
# #under MNAR assumption, pool the prevalence over the imputed datasets, resulting in absolute risks per study
# abs.outcome_mnar<-f.abs.poolrubin(data=prev.outcome_mnar,study_name="dis_id")
# colnames(abs.outcome_mnar)<-c("dis_id","p_MNAR","LC_MNAR","UC_MNAR")
# 
# data_prev<-merge(data_MCAR, abs.outcome_mar,by="dis_id")
# data_prev<-merge(data_prev,abs.outcome_mnar,by="dis_id")
# 
# data_long<- setDT(melt(data_prev, 
#                        measure.vars = c( "p_MCAR","p_MAR","p_MNAR",
#                                          "LC_MCAR","LC_MAR","LC_MNAR",
#                                          "UC_MCAR","UC_MAR","UC_MNAR"),
#                        variable.name = "Statistic", value.name = "Value"))
# 
# data_long[, c("Statistic", "Parameter") := tstrsplit(Statistic, "_", fixed=TRUE)]
# data_wide<- setDT(dcast(data_long,subreg_id+ dis_id + count+ Miss+Parameter ~ Statistic, value.var = "Value"))
# data_wide[,Parameter:=factor(Parameter,levels=c("MCAR","MAR","MNAR"))]
# levels(data_wide$Parameter)<-c("CC","2l.MAR","2l.Heckman")
# 
# pd <- position_dodge(width = 0.4)
# plot_dist2<-ggplot(data_wide, aes(x = dis_id, y = p, colour = Parameter)) +
#   geom_errorbar(aes(ymax = UC, ymin = LC),position = pd) + 
#   geom_point(position = pd,size=0.5)+ theme_light()+
#   ylab("Parasite prevalence (%)")+xlab("Sub region")+
#   scale_color_brewer(palette="Dark2")+facet_grid(subreg_id~.)

