#### Program information ######################################################
# Source file name: MenA_VaccSims.R                                           #
#_____________________________________________________________________________#
# PUrpose: run MenA simulations under different vaccination scenarios         #
#   specify scenario, country, region, and number of simulations to run       #
#_____________________________________________________________________________#
# Output:                                                                     #
#_____________________________________________________________________________#
# Author: Chris Stewart chris.c.stewartkp.org 3/2018                          #
#    adapted from Mike Jackson's SAS version                                  #
#_____________________________________________________________________________#
library(lubridate)
library(dplyr)
library(data.table)

##parameters to set:
begin<-Sys.time()
mycountry <- "ETH"
start <- as.Date("2001-01-01")
end <- as.Date("2100-12-31")
myregion <- "not_hyper"
PSA <- FALSE
Vaccination<-TRUE
program <- "campaign" ## "campaign" or "routine" or "both" or "none"
phi<-0.2
sd<-456 #seed for random sto
nSims<-10
###end parameters to set 

#script directory contains functions
setwd("\\\\HOME/stewcc1/MenAModel/R_programming")
source("ModelInputUtilities.R")
source("InitializePopulation.R")
source("MenA_OneSim.R")
source("MenA_helper_functions.R")
source("MenA_summarization_functions.R")

#directory containing inputs
inputdir<-"\\\\HOME/stewcc1/MenAModel/download"
#country-specific parameters
myparams<-GetDemographicParameters(path=inputdir,  mycountry=mycountry, start=start, end=end)
#setwd("\\\\HOME/stewcc1/MenAModel/Rdata")
#params<-read.csv("country_params.csv")
#myparams1 <-params[params$country==mycountry & params$year>=year(start)-1 & params$year<=year(end) ,]

#uSE SAME SEED FOR ALL SCENARIOS (for setting up stochastic parameter)
set.seed(sd, kind = NULL, normal.kind = NULL)
if (Vaccination!=FALSE) {
  myvacc<-GetVaccScenario(mycountry=mycountry, scenario=program, directory=inputdir)
  #vaccdf <- read.csv("VaxCover_ETH.csv")
  #myvaccold <- vaccdf[vaccdf$country=='ETH',]
  #make as vector of years where nothing happens (empty except for campaign only) for efficiency
  if (program=="campaign") {
    nodoses<-as.vector(myvacc[is.na(myvacc$DosesCampaign) | myvacc$DosesCampaign==0,"year"])
  }
}

#fixed parameters
paramfixed<-c(0.23334,0.70002,0.25,0.90,0.75,1.00,0.00000096)
names(paramfixed)<-c("rc","rd","lc","ld","hc","hd","foii")
#disease model for rainy and dry
dxrisk<-rbind(c(0.0018, -0.00000021),
              c(0.0019, -0.0000002))
dimnames(dxrisk)[[1]]<-c("rainy", "dry")
#WAIFW matrix setup
wboth<-GetWAIFWmatrix(path=inputdir, region=myregion)
#waifwin<-read.csv("WAIFW_both.csv", stringsAsFactors = FALSE)  #vector
#Rwaifw<-waifwin[waifwin$region==myregion & waifwin$season=='rainy', 4]
#Dwaifw<-waifwin[waifwin$region==myregion & waifwin$season=='dry', 4]
#wboth<-array(c(expandWaifw(waifw=Rwaifw), expandWaifw(waifw=Dwaifw)), dim=c(361,4,2))
#dimnames(wboth)[[3]]<-c("rainy", "dry")

#initialize population
startSize <- myparams[myparams$year==year(start)-1, 5]
initpop<-InitializePopulation(path=inputdir, start=start, end=end, popsize=startSize, country=mycountry, region=myregion)
#begin simulations
my_data <- list()
for (n in 1:nSims) {
  finalpop<-MenASimulation(startdt=start, enddt=end, pop=initpop,
                           fixedparams=paramfixed, countryparams=myparams,
                           WAIFWmx=wboth, dxr=dxrisk)
  #head(finalpop[,"Inc", 200:225])
  summarizeme<-1
  if (summarizeme > 0) {
    #age-specific death rates (for PSA=no, other option not implemented yet)
    cfr <- c(0.106, 0.096, 0.089, 0.086, 0.079, 0.122) #used AFTER simulation macro in SAS
    my_data[[n]] <-summarizeOneSim(finalpop, n, cfr)
    if (n==1) {
      cohortSize<-getCohortSize(finalpop)
      totalPop<-cohortSize%>%group_by(IterYear)%>%summarize(tot=sum(cohortsize))
      #for checking by plotting
    }
  } #end of conditional summarization
} #end of looping nSims simulations
#final summarization
filename = paste0(mycountry, "_", program, "_", Sys.Date(), ".csv")
#detail output - for testing
detail<-1
if (detail > 0) {
  tensims<-rbindlist(my_data[1:10])
  tensimsum<-tensims%>%group_by(simulation, IterYear)%>%summarize(sumCases=sum(Cases))
  sfile = paste0(mycountry, "_tensims_", program, "_", Sys.Date(), ".csv")
  detfile<-paste0("\\\\HOME/stewcc1/MenAModel/data/", sfile)
  write.csv(tensimsum, detfile)
  print(paste("Simulation detail written to", detfile))
  
}
filename = paste0(mycountry, "_", program, "_", Sys.Date(), ".csv")
finalsummary<-summarizeForOutput(my_data, cohort=cohortSize, write=TRUE, filename=filename)
print(begin)
print(Sys.time())