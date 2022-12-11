\name{dataObs}
\docType{data}
\alias{dataObs}
\alias{2l.2stage.heckman}

\title{ A two-level incomplete dataset based on an online obesity survey}
\description{
The dataset used here was based on data collected from 2111 individuals in an online obesity survey in different locations. The data was simplified and grouped into five clusters. The values and observability of the weight variable were defined according to the Heckman model in a herarchical model.
In addition, a systematic loss of this variable at one location was assumed. For the other variables included in the data set, missing values were simulated according to a MAR mechanism. 
}
\usage{
data(dataObs)
}

\format{
A dataframe with 2111 observations with the following variables:
 \tabular{rll}{
     \tab Gender \tab a factor value with two levels: 1 ("Female"), 0 ("Male").\cr
     \tab Age \tab a numeric value indicating age of subject in years.\cr
     \tab Height\tab a numeric value with Height in meters.\cr
     \tab FAVC\tab a factor value describing the frequent consumption of high caloric food (FAVC) with two levels:1("Yes"), 0("Male").\cr
     \tab Weight\tab a numeric value with Weight in Kilograms.\cr
     \tab Time\tab a numeric value indicating time in responding the questions in minutes.\cr
     \tab Cluster\tab a numeric indexing the cluster.\cr
  }
}

\source{Dataset obtained from "https://www.kaggle.com/datasets/fabinmndez/obesitydata?select=ObesityDataSet_raw_and_data_sinthetic.csv"}
\references{
Palechor, F. M., & de la Hoz Manotas, A. (2019). Dataset for estimation of obesity levels based on eating habits and physical condition in individuals from Colombia, Peru and Mexico. Data in brief, 25, 104344.
}

\details{
Simulation data code on gen_dataObs.R github repository
}

\examples{
rm(list=ls())
library(ggmice)
library(data.table)
data(dataObs)
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
}

\keyword{datasets}