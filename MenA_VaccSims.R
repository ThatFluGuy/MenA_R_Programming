

#### Program information ######################################################
# Source file name: MenA_VaccSims.R                                           #
#_____________________________________________________________________________#
# PUrpose: run MenA simulations under different vaccination scenarios         #
#   specify scenario, country, region, and number of simulations to run       #
#_____________________________________________________________________________#
# Input: Script builds some inputs from files downloaded                      #
# from vaccineimpact.org/montagu.  Location of the files is an input parameter#
#_____________________________________________________________________________#
# Output: two csv files: detailed output of first ten simulations, and        #
# summarized output of the requested number of simulations                    #
#_____________________________________________________________________________#
# Author: Chris Stewart chris.c.stewart@kp.org 2018                           #
#    adapted from Mike Jackson's SAS version                                  #
#_____________________________________________________________________________#
library(lubridate)
library(dplyr)
library(data.table) #melt

##parameters to set:
begin<-Sys.time()
mycountry <- "ERI"
start <- as.Date("2001-01-01")
end <- as.Date("2100-12-31")
myregion <- "hyper"  #"hyper" or "not_hyper"
PSA <- FALSE
vacc_program <- "routine" ## "campaign" or "routine" or "both" or "none"
phi<-0.2
sd<-4567 #seed for random sto, use same for all scenarios
nSims<-10  #100 takes ~ 3 mon, 1000 takes 45
#directory containing inputs from https://montagu.vaccineimpact.org/
inputdir<-"\\\\HOME/stewcc1/MenAModel/download"
outputdir<-"\\\\HOME/stewcc1/MenAModel/data"
#directory containing R scripts
script.dir <- "\\\\home/stewcc1/MenAModel/R_programming"
###end parameters to set 

#script directory contains functions
#ha, need to check script.dir so we can source parameter-checking script
if (dir.exists(script.dir)) {
  if (file.exists(paste0(script.dir, "/","MenA_paramCheck.R"))==FALSE) {
    msg<(paste0("MenA_paramcheck.R not found in ", script.dir))
    stop(msg)
  }
} else {
  script.dir<- getSrcDirectory(function(dummy) {dummy})
  if (file.exists(paste0(script.dir, "/","MenA_paramCheck.R"))==FALSE) {
    print(paste0("MenA_paramCheck.R not found in ", script.dir))
    stop("This script requires 6 other scripts; please put in same directory as this one or specify script directory.")
  }
}
setwd(script.dir)
source("MenA_paramCheck.R")
source("ModelInputUtilities.R")
source("InitializePopulation.R")
source("MenA_OneSim.R")
source("MenA_helper_functions.R")
source("MenA_summarization_functions.R")
#check parameters set above
setparams<-as.list(c(mycountry, as.character(start), as.character(end), myregion, PSA, vacc_program, phi, sd, nSims, inputdir, outputdir))
names(setparams)<-c("mycountry", "start", "end", "myregion", "PSA", "vacc_program", "phi", "sd", "nSims", "inputdir", "outputdir")
if (CheckSetParameters(setparams)==FALSE) {
    stop(spmessage)
} else {
  if (length(spmessage)>1) {
    print(spmessage)
  }
}

#country-specific parameters
myparams<-GetDemographicParameters(path=inputdir,  mycountry=mycountry, start=start, end=end)

if (vacc_program!="none") {
  myvacc<-GetVaccScenario(mycountry=mycountry, scenario=vacc_program, directory=inputdir)
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
if (!(is.numeric(wboth))) {
  stop(waifwerr)
}

#initialize population
startSize <- myparams[myparams$year==year(start)-1, "totalpop"]
initpop<-InitializePopulation(scriptdir=script.dir, inputdir=inputdir, start=start, end=end, popsize=startSize, country=mycountry, region=myregion)
#check for errors
if (!(is.numeric(initpop))) {
  if (disterr!="") { print(disterr) } 
  if (dxerr!="") { print(dxerr) } 
  stop(initmsg)
}
#begin simulations
my_data <- list()
for (n in 1:nSims) {
  finalpop<-MenASimulation(startdt=start, enddt=end, pop=initpop,
                           fixedparams=paramfixed, countryparams=myparams,
                           WAIFWmx=wboth, dxr=dxrisk, vacc_program=vacc_program)
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
filename = paste0(mycountry, "_", vacc_program, "_", Sys.Date(), ".csv")
#detail output - for testing
detail<-1
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
print(begin)
print(Sys.time())
