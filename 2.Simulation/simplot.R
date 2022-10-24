library(here)
library(data.table)
library(ggplot2)

# Arrange dataset
data<-as.data.table(readxl::read_xlsx(here('2.Simulation','Sim_results.xlsx'),sheet="Sheet1"))
head(data)
data<-data[,grep("Scenario|^Ns|^ni|Method|Variation|time|run|^Cov|^Qbar|^Re|^Emp|^Bias|^MSE|^MCMSE|^Mod", names(data)), with = FALSE]
colnames(data)<-c("Scenario","Ns","ni","Variation","Method",paste0("Mean_",0:6),"time", paste0("Cov_",0:2), paste0("Bias_",0:6),
                  paste0("MSE_",0:6),paste0("ModSE_",0:2),"run",paste0("EmpSE_",0:6),paste0("MCMSE_",0:6))
data<- setDT(melt(data, measure.vars = c( paste0("Mean_",0:6), paste0("Cov_",0:2),paste0("Bias_",0:6),paste0("MSE_",0:6),paste0("ModSE_",0:2),paste0("EmpSE_",0:6)
                                          ,paste0("MCMSE_",0:6)),
                  variable.name = "Statistic", value.name = "Value"))
data[, c("Statistic", "Parameter") := tstrsplit(Statistic, "_", fixed=TRUE)]
data<- setDT(dcast(data, Scenario + Ns + ni+ Variation+ Method + time + run + Parameter ~ Statistic, value.var = "Value"))
data[,SeCov:=sqrt(Cov*(1-Cov)/run)]
data[,Coverage:=Cov]
data[,SEBias:=EmpSE/run^0.5]
data[,Per_run:=run/5]
param_val <- c("0"="beta_0", "1"="beta_1", "2"="beta_2", 
               "3"="sigma_b0", "4"="sigma_b1", "5"="sigma_b2","6"= "sigma_e")

data[,Parameter := as.factor(param_val[data$Parameter])]

unique(data$Scenario)
datagraph<-data[Scenario=="Rho",c("Scenario","Variation","Method","Parameter","Bias","Cov","MSE","SEBias","SeCov","MCMSE")]
colnames(datagraph)<-c("Scenario","Variation","Method","Parameter","p_Bias","p_Cov","p_MSE","SE_Bias","SE_Cov","SE_MSE")

datagraph<-data.table::melt(setDT(datagraph), 
     measure = patterns("p_", "SE_"),
     variable.name = 'var', value.name = c('p', 'se'))
setDT(datagraph)[,criteria:=ifelse(var==1,"Bias",ifelse(var==2,"Coverage","MSE"))]

head(datagraph)
pd = position_dodge(.1)
ggplot(datagraph[Parameter%in%c("beta_0","beta_1","beta_2","sigma_b0","sigma_b1","sigma_b2")],aes(x= Variation,
                     y     = p,
                     color = Method)) +
  geom_point(shape = 1,
             size  = 1,
             position = pd) +
  geom_errorbar(aes(ymin  = p - 1.96*se,
                    ymax  = p + 1.96*se),
                width = 0.2,
                size  = 0.5,
                position = pd) +
  facet_grid(criteria~Parameter, scales="free",labeller = label_parsed)

+
  guides(color=guide_legend(title="Imputation method"))+
  theme(axis.text.x = element_text(angle =0,vjust=0.5,size=8),legend.position = "none")+
  ylab("Estimate")+xlab(varlab)+
  scale_color_brewer(palette="Dark2")
#scale_color_viridis(discrete=TRUE)



plot_cont_rho<-ggplot(long,aes(x=Variation,y=value,colour= Method,group=Method))+
  geom_point(aes(shape = Method),size  = 1) +
  facet_grid(criteria~Parameter, scales="free",labeller = label_parsed)+
  theme_light()+
  geom_hline(data = datalines1, aes(yintercept = Z),linetype='dotted')+
  theme(axis.text.x = element_text(angle =0,vjust=0.5,size=8),legend.position ="bottom")+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black'))+
  ylab("Value")+xlab(expression(rho))+
  scale_color_brewer(palette="Dark2")






# rho variation
data1<-data[Scenario=="Rho"&Method!="Heck_IPD_F"&Parameter%in%c("beta_1","sigma_b1","beta_2","sigma_b2"),]
data1[,Parameter := factor(Parameter, levels = c("beta_0","beta_1","beta_2","sigma_b0","sigma_b1","sigma_b2","sigma_e"),
                           ordered = TRUE, labels=c(expression(beta[0]),expression(beta[1]),expression(beta[2]),
                                                    expression(sigma[b*0]),expression(sigma[b*1]),expression(sigma[b*2]),
                                                    expression(sigma[e])))]
datalines1 <- data.table(Parameter = rep(c("beta_1", "beta_2","sigma_b1","sigma_b2"),4),
                         criteria = rep(c("Bias","Coverage","RMSE"),each=4),
                         Z = c(0,0,0,0,0.95,0.95,NA,NA,0,0,0,0))
datalines1[,Parameter := factor(Parameter, levels = c("beta_0","beta_1","beta_2","sigma_b0","sigma_b1","sigma_b2","sigma_e"),
                           ordered = TRUE, labels=c(expression(beta[0]),expression(beta[1]),expression(beta[2]),
                                                    expression(sigma[b*0]),expression(sigma[b*1]),expression(sigma[b*2]),
                                                    expression(sigma[e])))]

long <- melt(setDT(data1[,c("Variation","Method","Parameter","Bias","Coverage","RMSE")]),
             id.vars = c("Variation", "Method","Parameter"), variable.name = "criteria")
long$Variation<-as.factor(long$Variation)
long$Method<-as.factor(long$Method)
levels(long$Method)
long$Method <- factor(long$Method, levels = c("CC","1l.Heckman","2l.MAR","2l.Heckman"))


plot_cont_rho<-ggplot(long,aes(x=Variation,y=value,colour= Method,group=Method))+
  geom_point(aes(shape = Method),size  = 1) +
  facet_grid(criteria~Parameter, scales="free",labeller = label_parsed)+
  theme_light()+
  geom_hline(data = datalines1, aes(yintercept = Z),linetype='dotted')+
  theme(axis.text.x = element_text(angle =0,vjust=0.5,size=8),legend.position ="bottom")+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black'))+
  ylab("Value")+xlab(expression(rho))+
  scale_color_brewer(palette="Dark2")

# size variation
datas<-data[Scenario=="N"&Method!="Heck_IPD_F"&Parameter%in%c("beta_1","sigma_b1","beta_2","sigma_b2"),]
datas[,Parameter := factor(Parameter, levels = c("beta_0","beta_1","beta_2","sigma_b0","sigma_b1","sigma_b2","sigma_e"),
                           ordered = TRUE, labels=c(expression(beta[0]),expression(beta[1]),expression(beta[2]),
                                                    expression(sigma[b*0]),expression(sigma[b*1]),expression(sigma[b*2]),
                                                    expression(sigma[e])))]
datalines1 <- data.table(Parameter = rep(c("beta_1", "beta_2","sigma_b1","sigma_b2"),4),
                         criteria = rep(c("Bias","Coverage","EmpSE","RMSE"),each=4),
                         Z = c(0,0,0,0,0.95,0.95,NA,NA,NA,NA,NA,NA,0,0,0,0))
datalines1[,Parameter := factor(Parameter, levels = c("beta_0","beta_1","beta_2","sigma_b0","sigma_b1","sigma_b2","sigma_e"),
                                ordered = TRUE, labels=c(expression(beta[0]),expression(beta[1]),expression(beta[2]),
                                                         expression(sigma[b*0]),expression(sigma[b*1]),expression(sigma[b*2]),
                                                         expression(sigma[e])))]

longs <- melt(setDT(datas[,c("Variation","Method","Parameter","Bias","EmpSE","Coverage","RMSE")]),
             id.vars = c("Variation", "Method","Parameter"), variable.name = "criteria")
longs$Variation<-as.factor(longs$Variation)
levels(longs$Variation)<-c("N=10 \n ni=1000","N=10 \n ni=50","N=100 \n ni=1000")
longs$Variation <- factor(longs$Variation, levels = c("N=10 \n ni=50","N=10 \n ni=1000","N=100 \n ni=1000"))
longs$Method<-as.factor(longs$Method)
longs$Method <- factor(longs$Method, levels = c("CC","1l.Heckman","2l.MAR","2l.Heckman"))


figure_final_s<-ggplot(longs,aes(x=Variation,y=value,colour= Method,group=Method))+
  geom_point(aes(shape = Method),size  = 1) +
  facet_grid(criteria~Parameter, scales="free",labeller = label_parsed)+
  theme_light()+
  geom_hline(data = datalines1, aes(yintercept = Z),linetype='dotted')+
  theme(axis.text.x = element_text(angle =0,vjust=0.5,size=8),legend.position ="bottom")+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black'))+
  ylab("Value")+xlab("Number of clusters(N) and sample size(ni)")+
  scale_color_brewer(palette="Dark2")
 
  
