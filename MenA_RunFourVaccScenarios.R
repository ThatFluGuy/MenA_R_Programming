
#### Program information ######################################################
# Source file name: MenA_RunFourVaccScenarios.R                               #
# Package: MenA_VaccSims                                                      #
# Contact: chris.c.stewart@kp.org, michael.l.jackson@kp.org                   #
# Version Date 12/31/2018                                                     #
#_____________________________________________________________________________#
# PUrpose: Adapted from MenA_VaccSims to run 4 vaccination scenarios at once  #
#   specify country, start, end, region, and number of simulations to run     #
#   Please start with MenA_VaccSims.R and a small number of simulations       #
#   when first using this package                                             #  
#_____________________________________________________________________________#
# Output: 2 .csv files of subsetted detail(optional, use detail parameter)    #
#   and summary output, plus two charts (optional, use charts parameter)      #
#_____________________________________________________________________________#
# Author: Chris Stewart chris.c.stewart@kp.org 2018                           #
#    adapted from Mike Jackson's SAS version                                  #
#_____________________________________________________________________________#
library(lubridate)
library(dplyr)
library(data.table)

##parameters to set:
begin<-Sys.time()
mycountry <- "COD"
start <- as.Date("2001-01-01")
end <- as.Date("2100-12-31")
myregion <- "hyper"
PSA <- FALSE

phi<-0.2
sd<-46209 #seed for random sto, use same for all scenarios
nSims<-10
#directory containing inputs from https://montagu.vaccineimpact.org/
inputdir<-"\\\\HOME/stewcc1/MenAModel/download"
outputdir<-"\\\\HOME/stewcc1/MenAModel/data"
#directory containing R scripts
script.dir <- "\\\\home/stewcc1/MenAModel/R_programming"
detail<-1
charts<-TRUE
###end parameters to set 

#script directory contains functions
setwd(script.dir)
source("ModelInputUtilities.R")
source("MenA_OneSim.R")
source("MenA_helper_functions.R")
source("MenA_summarization_functions.R")
source("MenA_plotting_functions.R")

#country-specific parameters
myparams<-GetDemographicParameters(path=inputdir,  mycountry=mycountry, start=start, end=end)
#vaccination scenario loop
vacc_vector <- c("none", "campaign", "routine", "both")
for (x in 1:length(vacc_vector)) {
  vacc_program<-vacc_vector[x]
  if (vacc_program!="none") {
    myvacc<-GetVaccScenario(mycountry=mycountry, scenario=vacc_program, directory=inputdir)
    #vaccdf <- read.csv("VaxCover_ETH.csv")
    #myvaccold <- vaccdf[vaccdf$country=='ETH',]
    #make as vector of years where nothing happens (empty except for campaign only) for efficiency
    if (vacc_program=="campaign") {
      nodoses<-as.vector(myvacc[is.na(myvacc$DosesCampaign) | myvacc$DosesCampaign==0,"year"])
    }
  }
  
  #fixed parameters
  #uSE SAME SEED FOR ALL SCENARIOS (for setting up stochastic parameter)
  set.seed(sd, kind = NULL, normal.kind = NULL)
  paramfixed<-c(0.23334,0.70002,0.25,0.90,0.75,1.00,0.00000096)
  names(paramfixed)<-c("rc","rd","lc","ld","hc","hd","foii")
  #disease model for rainy and dry
  dxrisk<-rbind(c(0.0018, -0.00000021),
                c(0.0019, -0.0000002))
  dimnames(dxrisk)[[1]]<-c("rainy", "dry")
  #WAIFW matrix setup
  wboth<-GetWAIFWmatrix(path=script.dir, region=myregion)
  
  #initialize population
  startSize <- myparams[myparams$year==year(start)-1, "totalpop"]
  initpop<-InitializePopulation(scriptdir=script.dir, inputdir=inputdir, start=start, end=end, popsize=startSize, country=mycountry, region=myregion)
  #begin simulations
  my_data <- list()
  for (n in 1:nSims) {
    finalpop<-MenASimulation(startdt=start, enddt=end, pop=initpop,
                             fixedparams=paramfixed, countryparams=myparams,
                             WAIFWmx=wboth, dxr=dxrisk, vacc_program=vacc_program)
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
  filename = paste0(mycountry, "_", vacc_program, "_", Sys.Date(), ".csv")
  #detail output - for testing
  if (detail > 0) {
    tensims<-rbindlist(my_data[1:10])
    tensimsum<-tensims%>%group_by(simulation, IterYear)%>%summarize(sumCases=sum(Cases))
    sfile = paste0(mycountry, "_tensims_", vacc_program, "_", Sys.Date(), ".csv")
    detfile<-paste0(outputdir, "/", sfile)
    write.csv(tensimsum, detfile)
    print(paste("Simulation detail written to", detfile))
  }
  filename <- paste0(mycountry, "_", vacc_program, "_", Sys.Date(), ".csv")
  filename1<- paste0(outputdir, "/", filename)
  finalsummary<-summarizeForOutput(my_data, cohort=cohortSize, write=TRUE, filename=filename1)
}
if (charts==TRUE) {
  mygrid<-GridofTenSims(datadir=outputdir, mycountry=mycountry, saveit=TRUE)
  myplot<-PlotAllScenarios(datadir=outputdir, mycountry=mycountry, saveit=TRUE) 
  print(paste("Charts written to", outputdir))
}
print(begin)
print(Sys.time())