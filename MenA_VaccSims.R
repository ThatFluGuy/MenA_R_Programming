

#### Program information ######################################################
# Source file name: MenA_VaccSims.R                                           #
#### Program information ######################################################
# Package: MenA_VaccSims                                                      #
# Source file name: MenA_VaccSims.R                                           #
# Version Date 12/31/2018                                                     #
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
library(doParallel)
library(dplyr)
library(data.table) #melt
# Chloe 2/5/19: Added package to use data.table
library(reshape2)
library(tidyr)

##parameters to set:
begin<-Sys.time()
mycountry <- "NGA"
start <- as.Date("1951-01-01")
end <- as.Date("2100-12-31")
## NOTE MANY FUNCTIONS APPEAR TO DEPEND ON THESE BEING FIRST AND LAST DATES: EDIT LATER.
myregion <- "hyper"  #"hyper" or "not_hyper"
PSA <- FALSE
vacc_program <- "campaign" ## "campaign" or "routine" or "both" or "none"
vacc_subprogram <- "default"  ## "default" or "bestcase" are allowable options in 2019v3
phi<-0.2
sd<-4567 #seed for random sto, use same for all scenarios
nSims<-117  #Update: 100 takes around 12 minutes if using 4 cores.
#directory containing inputs from https://montagu.vaccineimpact.org/
inputdir<-"G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Data/GAVI inputs/2019_12_gavi_v3"
#outputdir<-"C:/Users/krakcx1/Desktop/Sim_output"
outputdir <- "G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Scratch/Temporary output directory"

#directory containing R scripts
#script.dir <- "C:/Users/krakcx1/Desktop/Cloned_meningitis_code"
#script.dir <- "H:/Git/MenA R"
script.dir <- "C:/Users/jackml4/Documents/Link_To_H_Drive/GAVI MenA predictions/R_programming"
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
source("MenA_OneSim.R")
source("MenA_helper_functions.R")
source("MenA_summarization_functions.R")
source("MenA_calibration_plots.R")
#check parameters set above
setparams<-as.list(c(mycountry, as.character(start), as.character(end), myregion, PSA, vacc_program, phi, sd, nSims, inputdir, outputdir))
names(setparams)<-c("mycountry", "start", "end", "myregion", "PSA", "vacc_program", "phi", "sd", "nSims", "inputdir", "outputdir")
if (CheckSetParameters(setparams)==FALSE) {
    stop(spmessage)
} else {
  if (length(spmessage)>1) { print(spmessage) }
}

# Import country-specific parameters
myparams<-GetDemographicParameters(path=inputdir,  mycountry=mycountry, start=start, end=end)
if (CheckDemogParameters(myparams)==FALSE) {
  stop(dpmessage)
} else {
  if (length(dpmessage)>1) { print(dpmessage) }
}

# Import vaccination program details
if (vacc_program!="none") {
  myvacc<-GetVaccScenario(mycountry=mycountry, scenario=vacc_program, sub.scenario=vacc_subprogram, directory=inputdir)
  if (is.data.frame(myvacc)==FALSE) { stop(vaccmsg)}  #check for output
  #make as vector of years where nothing happens (empty except for campaign only) for efficiency
  if (vacc_program=="campaign") {
    nodoses<-as.vector(myvacc[is.na(myvacc$DosesCampaign) | myvacc$DosesCampaign==0,"year"])
  }
}

# Country-specific life expectancy
my.lifex <- GetLifeExp(path=inputdir, mycountry.s=mycountry)

#Read in parameters calculated in ABC, or a row of parameters to be used by ABC.  Points to model_params.csv
#This file should supply all calibrated parameters, removing all hard-coding throughout the model.
paramfixed <- GetModelParams(path=script.dir, region.val=myregion)

#initialize population
startSize <- myparams[myparams$year==year(start)-1, "totalpop"]
initpop<-InitializePopulation(scriptdir=script.dir, inputdir=inputdir, start=start, end=end, country=mycountry, region=myregion, startSize=startSize)
#check for errors
if (!(is.numeric(initpop))) {
  if (disterr!="") { print(disterr) } 
  if (dxerr!="") { print(dxerr) } 
  stop(initmsg)
}

summarizeme<-1
#set the random number seed based on sd, then create a vector of random number seeds (consistent within sd values)
set.seed(sd, kind = NULL, normal.kind = NULL)
seed.vec <- unique(floor(runif(nSims*2, 0, 1000000)))[1:nSims]


#begin simulations
cl <- makeCluster(4)  #scale this upwards if you're on a workstation with >16gb memory
registerDoParallel(cl)
my_data <- foreach(n=1:nSims, .verbose=TRUE, .packages = c("lubridate", "dplyr", "data.table", "reshape2")) %dopar% {
  set.seed(seed.vec[n])
  paramfixed.row <- paramfixed[n,]  
  finalpop<-MenASimulation(startdt=start, enddt=end, fp=paramfixed.row, initpop=initpop, vacc_program=vacc_program,
                           countryparams=myparams, region=myregion, country=mycountry, inputdir=inputdir)
  if (summarizeme > 0) {
    #age-specific death rates (for PSA=no, other option not implemented yet)
    cfr <- c(0.106, 0.096, 0.089, 0.086, 0.079, 0.122) #used AFTER simulation macro in SAS
    summarizeOneSim(finalpop, n, cfr)
    } #end of conditional summarization
} #end of foreach loop
stopCluster(cl)



#for checking by plotting
onerun <- MenASimulation(startdt=start, enddt=end, fp=paramfixed[4,], initpop=initpop, vacc_program=vacc_program,
                         countryparams=myparams, region=myregion, country=mycountry, inputdir=inputdir)
cohortSize <- getCohortSize(onerun)
totalPop <- cohortSize %>% 
  group_by(year) %>% summarize(tot=sum(cohortsize))



mena.counts(onerun, 1, 60)
mena.counts(onerun, 61, 120)
mena.counts(onerun, 121, 180)
mena.counts(onerun, 181, 240)
mena.counts(onerun, 241, 360)
mena.counts(onerun, 361, 1441)
mena.counts(onerun, 961, 1441)
mena.counts(onerun, 1, 1441)



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
write.csv(cohortSize, file=paste(outputdir, "/cohortsize.csv", sep=""))
finalsummary<-summarizeForOutput(my_data, cohort=cohortSize, write=TRUE, filename=filename1)
print(begin)
print(Sys.time())
