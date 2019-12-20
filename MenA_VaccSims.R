#### Program information ######################################################
# Package: MenA_VaccSims                                                      #
# Source file name: MenA_VaccSims.R                                           #
# Version Date 12/17/2019                                                     #
#_____________________________________________________________________________#
# PUrpose: Run multiple iterations of MenA simulations in a single country    #
# under a specific vaccination scenario                                       #
#   specify scenario, country, region, and number of simulations to run       #
#_____________________________________________________________________________#
# Input: Script builds some inputs from files downloaded                      #
# from vaccineimpact.org/montagu.  Location of the files is an input parameter#
#_____________________________________________________________________________#
# Output: two csv files: detailed output of first ten simulations, and        #
# summarized output of the requested number of simulations                    #
#_____________________________________________________________________________#
# Steps in this program:                                                      #
# (1) Set up libraries                                                        #
# (2) User-set parameters                                                     #
# (3) Import and format data/functions                                        #
# (4) Run simulations                                                         #
#_____________________________________________________________________________#
# Author: Chris Stewart chris.c.stewart@kp.org 2018                           #
#    adapted from Mike Jackson's SAS version                                  #
#_____________________________________________________________________________#


### (1) Set up libraries ######################################################

library(lubridate)
library(doParallel)
library(dplyr)
library(data.table)
library(reshape2)
library(tidyr)

### (2) User-set parameters ###################################################
# Edit the code in this section to specify country and scenario.              #
# Edit the code in this section to specify country and scenario.              #
# Update 2019.12.20 - include option to pseudo-automate. If set to true, the  #
# program will find the first country/scenario combination that has not yet   #
# been completed and will run it, from the scenario_tracker.csv file. If set  #
# to false, manually enter the country, program option, and region type.      #

automate <- TRUE

begin <- Sys.time()
start <- as.Date("1951-01-01") # Use 1/1/1951 to give 50 years burn-in
end <- as.Date("2100-12-31")
PSA <- FALSE
sd <- 4567 # Seed for random sto, use same for all scenarios
nSims <- 117  # Update: 100 takes around 12 minutes if using 4 cores.
use.tensims <- FALSE # If desired, output results from 10 sims for debugging

# Directory containing inputs from https://montagu.vaccineimpact.org/
input.dir<-"G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Data/GAVI inputs/2019_12_gavi_v3"
# Directory for simulation outputs
output.dir <- "G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Analysis/Simulation results"
# Directory containing R scripts
script.dir <- "C:/Users/jackml4/documents/Link_to_H_drive/GAVI MenA predictions/R_programming"

if (automate==TRUE){
  scenario_tracker <- read.csv(paste(input.dir, "scenario_tracker.csv", sep="/"),
                               stringsAsFactors = FALSE)
  not_done <- which(scenario_tracker$completed=="FALSE")
  if(length(not_done) >= 1){
    scen_num <- not_done[1]
    mycountry <- scenario_tracker$country_code[scen_num]
    myregion <- scenario_tracker$region[scen_num]
    vacc_program <- scenario_tracker$vacc_program[scen_num]
    vacc_subprogram <- scenario_tracker$vacc_subprogram[scen_num]
  } else {
    print("All scenarios have been completed")
  }
} else {
  mycountry <- "NGA"
  myregion <- "hyper"  #"hyper" or "not_hyper"
  vacc_program <- "campaign" ## "campaign" or "routine" or "both" or "none"
  vacc_subprogram <- "default"  ## "default" or "bestcase" are allowable options in 2019v3
}

### (3) Import and format data/functions ######################################
# Import scripts used in the simulations, country-specific parameters, and    #
# vaccination program details.                                                #

# Set stochastic parameter
phi<-0.2

# (A) Import scripts
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

# (B) Check parameters set above
setparams<-as.list(c(mycountry, as.character(start), as.character(end), myregion, PSA, vacc_program, phi, sd, nSims, input.dir, output.dir))
names(setparams)<-c("mycountry", "start", "end", "myregion", "PSA", "vacc_program", "phi", "sd", "nSims", "input.dir", "output.dir")
if (CheckSetParameters(setparams)==FALSE) {
    stop(spmessage)
} else {
  if (length(spmessage)>1) { print(spmessage) }
}

# (C) Import country-specific parameters
myparams.full <- GetDemographicParameters(path=input.dir,  mycountry=mycountry, start=start, end=end)
if (CheckDemogParameters(myparams.full)==FALSE) {
  stop(dpmessage)
} else {
  if (length(dpmessage)>1) { print(dpmessage) }
}

# (D) Import vaccination program details
if (vacc_program!="none") {
  myvacc<-GetVaccScenario(mycountry=mycountry, scenario=vacc_program, sub.scenario=vacc_subprogram, directory=input.dir)
  if (is.data.frame(myvacc)==FALSE) { stop(vaccmsg)}  #check for output
  #make as vector of years where nothing happens (empty except for campaign only) for efficiency
  if (vacc_program=="campaign") {
    nodoses<-as.vector(myvacc[is.na(myvacc$DosesCampaign) | myvacc$DosesCampaign==0,"year"])
  }
}

# (E) Country-specific life expectancy
my.lifex <- GetLifeExp(path=input.dir, mycountry.s=mycountry)

# (F) Read in parameters calculated in ABC, or a row of parameters to be used by ABC.  
paramfixed <- GetModelParams(path=script.dir, region.val=myregion)

# (G) Initialize full population.
startSize <- myparams.full[myparams.full$year==year(start)-1, "totalpop"]
initpop.full <- InitializePopulation(scriptdir=script.dir, inputdir=input.dir, start=start, end=end, country=mycountry, region=myregion, startSize=startSize)
#check for errors
if (!(is.numeric(initpop.full))) {
  if (disterr!="") { print(disterr) } 
  if (dxerr!="") { print(dxerr) } 
  stop(initmsg)
}

# (H) Scale down the full population to the modeled population
# This means changing both the starting population size and
# the annual number of births.
pct.modeled <- GetModelPct(path=input.dir, mycountry.s=mycountry)
initpop <- initpop.full
initpop[,,1] <- initpop[,,1] * pct.modeled

myparams <- myparams.full
myparams$births <- myparams$births * pct.modeled

### (4) Run simulations #######################################################
# Set up random seed vector and use parallel processing to run simulations.   #

summarizeme <- 1
# Set the random number seed based on sd, then create a vector of random number seeds (consistent within sd values)
set.seed(sd, kind = NULL, normal.kind = NULL)
seed.vec <- unique(floor(runif(nSims*2, 0, 1000000)))[1:nSims]

# Begin simulations
cl <- makeCluster(4)  #scale this upwards if you're on a workstation with >16gb memory
registerDoParallel(cl)
my_data <- foreach(n=1:nSims, .verbose=TRUE, .packages = c("lubridate", "dplyr", "data.table", "reshape2")) %dopar% {
  set.seed(seed.vec[n])
  paramfixed.row <- paramfixed[n,]  
  finalpop<-MenASimulation(startdt=start, enddt=end, fp=paramfixed.row, initpop=initpop, vacc_program=vacc_program,
                           countryparams=myparams, region=myregion, country=mycountry, inputdir=input.dir)
  if (summarizeme > 0) {
    #age-specific death rates (for PSA=no, other option not implemented yet)
    cfr <- c(0.106, 0.096, 0.089, 0.086, 0.079, 0.122) #used AFTER simulation macro in SAS
    summarizeOneSim(finalpop, n, cfr)
    } #end of conditional summarization
} #end of foreach loop
stopCluster(cl)


# Also run one time to get the cohort size
# VIMC wants the size of the full population, not just the vaccine-
# targetted population, so run this without scaling down
onerun <- MenASimulation(startdt=start, enddt=end, fp=paramfixed[4,], initpop=initpop.full, vacc_program=vacc_program,
                         countryparams=myparams.full, region=myregion, country=mycountry, inputdir=input.dir)
cohortSize <- getCohortSize(onerun)
totalPop <- cohortSize %>% 
  group_by(year) %>% summarize(tot=sum(cohortsize))


# Final summarization
filename = paste0(mycountry, "_", vacc_program, "_", Sys.Date(), ".csv")
# Detail output - for testing
detail<-1
if (use.tensims==TRUE) {
  tensims<-rbindlist(my_data[1:10])
  tensimsum<-tensims%>%group_by(simulation, IterYear)%>%summarize(sumCases=sum(Cases))
  sfile = paste0(mycountry, "_tensims_", vacc_program, "_", Sys.Date(), ".csv")
  detfile<-paste0(output.dir, "/", sfile)
  write.csv(tensimsum, detfile)
  print(paste("Simulation detail written to", detfile))
}

# Output summary results for the country/scenario set
filename <- paste0(mycountry, "_", vacc_program, "_", vacc_subprogram, "_", Sys.Date(), ".csv")
filename1<- paste0(output.dir, "/", filename)
finalsummary<-summarizeForOutput(my_data, cohort=cohortSize, write=TRUE, filename=filename1)


# Update scenario tracker
if (automate==TRUE){
  scenario_tracker$completed[scen_num] <- "TRUE"
  write.csv(scenario_tracker, paste(input.dir, "scenario_tracker.csv", sep="/"))
}

print(begin)
print(Sys.time())
