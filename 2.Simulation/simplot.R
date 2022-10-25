library(here)
library(data.table)
library(ggplot2)

# Arrange dataset
data<-as.data.table(readxl::read_xlsx(here('2.Simulation','Sim_results.xlsx'),sheet="Sheet3"))
data<-data[,grep("Scenario|^Ns|^ni|Method|Variation|time|run|^Cov|^Qbar|^Re|^Emp|^Bias|^MSE|^MCMSE|^Mod|^SE", names(data)), with = FALSE]
colnames(data)<-c("Scenario","Ns","ni","Variation","Method",paste0("Mean_",0:6),"time", paste0("Cov_",0:2), paste0("Bias_",0:6),
                  paste0("MSE_",0:6),paste0("ModSE_",0:2),"run",paste0("EmpSE_",0:6),paste0("MCMSE_",0:6),paste0("SEEmp_",0:6))
data<- setDT(melt(data, measure.vars = c( paste0("Mean_",0:6), paste0("Cov_",0:2),paste0("Bias_",0:6),paste0("MSE_",0:6),paste0("ModSE_",0:2),paste0("EmpSE_",0:6)
                                          ,paste0("MCMSE_",0:6),paste0("SEEmp_",0:6)),
                  variable.name = "Statistic", value.name = "Value"))
data[, c("Statistic", "Parameter") := tstrsplit(Statistic, "_", fixed=TRUE)]
data<- setDT(dcast(data, Scenario + Ns + ni+ Variation+ Method + time + run + Parameter ~ Statistic, value.var = "Value"))
data[,SeCov:=sqrt(Cov*(1-Cov)/run)]
data[,Coverage:=Cov]
data[,SEBias:=EmpSE/run^0.5]
data[,Per_run:=run/5]
param_val <- c("0"="beta_0", "1"="beta_1", "2"="beta_2", 
               "3"="sigma_b0", "4"="sigma_b1", "5"="sigma_b2","6"= "sigma_e")
data[,Method:=factor(Method, levels = c("CC","1l.Heckman","2l.MAR","2l.Heckman"))]

data[,Parameter := as.factor(param_val[data$Parameter])]
data[,Parameter := factor(Parameter, levels = c("beta_0","beta_1","beta_2","sigma_b0","sigma_b1","sigma_b2","sigma_e"),
                          ordered = TRUE, labels=c(expression(beta[0]),expression(beta[1]),expression(beta[2]),
                                                   expression(sigma[b*0]),expression(sigma[b*1]),expression(sigma[b*2]),
                                                   expression(sigma[e])))]


plotsim<-function(Scenariov,levelsvar){
datagraph<-data[Scenario==Scenariov,c("Scenario","Variation","Method","Parameter","Bias","SEBias","Cov","SeCov","EmpSE","SEEmp","MSE","MCMSE")]
colnames(datagraph)<-c("Scenario","Variation","Method","Parameter","p_Bias","SE_Bias","p_Cov","SE_Cov","p_EmpSE","SE_EmpSE","p_MSE","SE_MSE")

datagraph<-data.table::melt(setDT(datagraph), 
                            measure = patterns("p_", "SE_"),
                            variable.name = 'var', value.name = c('p', 'se'))
setDT(datagraph)[,criteria:=fcase(var==1,"Bias",var==2,"Coverage",var==3,"EmpSE",var==4,"MSE")]

datalines <- data.table(Parameter = c("beta_0","beta_1", "beta_2","sigma_b0", "sigma_b1","sigma_b2","sigma_e"),
                        Z = c(rep(0,7),rep(0.95,7),rep(0,7),rep(0,7)),
                        criteria=c(rep("Bias",7),rep("Coverage",7),rep("EmpSE",7),rep("MSE",7)))
datalines[,Parameter := factor(Parameter, levels = c("beta_0","beta_1","beta_2","sigma_b0","sigma_b1","sigma_b2","sigma_e"),
                               ordered = TRUE, labels=c(expression(beta[0]),expression(beta[1]),expression(beta[2]),
                                                        expression(sigma[b*0]),expression(sigma[b*1]),expression(sigma[b*2]),
                                                        expression(sigma[e])))]
datagraph[, Variation:=as.factor(Variation)]
datagraph[, Variation:=factor(Variation,levels=levelsvar)]
levels(datagraph$Variation)
pd = position_dodge(.1)
datalines1<-datalines[Parameter%in%c("beta[0]","beta[1]","beta[2]","sigma[b * 0]","sigma[b * 1]","sigma[b * 2]")]

plot<-ggplot(datagraph[Parameter%in%c("beta[0]","beta[1]","beta[2]","sigma[b * 0]","sigma[b * 1]","sigma[b * 2]")],
       aes(y= Variation, x     = p, color = Method)) +
  geom_point(shape = 1,size  = 1,position = pd) +
  geom_vline(data = datalines1, aes(xintercept = Z),linetype='dotted')+
  geom_errorbar(aes(xmin  = p - 1.96*se, xmax  = p + 1.96*se),
                width = 0.2,size  = 0.5,position = pd) +
  facet_grid(Parameter~criteria, scales="free",labeller = label_parsed)+
  scale_color_brewer(palette="Dark2")+
  guides(color=guide_legend(title="Imputation method"))+
  theme_light()+  
  theme(strip.background =element_rect(fill="white"),legend.position="bottom")+
  theme(strip.text = element_text(colour = 'black'))+
  xlab("Value")+ylab(expression(rho))
return(plot)
}


plot_rho<-plotsim(Scenariov="Rho3",levelsvar=c("0","0.3","0.6","0.9"))
plot_bin<-plotsim(Scenariov="Bin3",levelsvar=c("0","0.3","0.6","0.9"))
plot_N<-plotsim(Scenariov="N3",levelsvar=c("10-50","10-100","10-1000","50-1000","100-1000"))
plot_S<-plotsim(Scenariov="S3",levelsvar=c("e~BVGamma","e~BVN","R=f(Y)" ))
