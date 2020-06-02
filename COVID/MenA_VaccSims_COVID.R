#### Program information ######################################################
# Package: MenA_VaccSims                                                      #
# Source file name: MenA_VaccSims_COVID.R                                     #
# Version Date 12/17/2019                                                     #
#_____________________________________________________________________________#
# PUrpose: Run multiple iterations of MenA simulations in a single country    #
# under a specific vaccination scenario                                       #
#   specify scenario, country, region, and number of simulations to run       #
# 2020.05.20: This version is for the COVID modeling, to explore the possible #
# effects of COVID-related interruptions in vaccination programs on NmA.      #
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
# Update 2020.05.20 - for COVID, remove automation                            #

begin <- Sys.time()
start <- as.Date("1951-01-01") # Use 1/1/1951 to give 50 years burn-in
end <- as.Date("2100-12-31")
PSA <- TRUE
seed <- 4567 # Seed for random sto, use same for all scenarios
nSims <- 200  # Update: 100 takes around 12 minutes if using 4 cores.
use.tensims <- FALSE # If desired, output results from 10 sims for debugging

# Directory containing inputs from https://montagu.vaccineimpact.org/
input.dir<-"G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Data/GAVI inputs/202005gavi_v4"
# Directory for simulation outputs
output.dir <- "G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Analysis/Simulation results"
# Directory containing R scripts
script.dir <- "C:/Users/O992928/Documents/GAVI MenA predictions/R_programming/COVID"

automate <- TRUE

if (automate==TRUE){
  scenario_tracker <- read.csv(paste(input.dir, "scenario_tracker.csv", sep="/"),
                               stringsAsFactors = FALSE)
  not_done <- which(scenario_tracker$completed=="FALSE")
  if(length(not_done) >= 1){
    scen_num <- not_done[1]
    mycountry <- scenario_tracker$country_code[scen_num]
    myregion <- scenario_tracker$region[scen_num]
    vacc_program <- scenario_tracker$scenario[scen_num]
  } else {
    print("All scenarios have been completed")
  }
  
} else {
  mycountry <- "NGA"
  myregion <- "hyper"  #"hyper" or "not_hyper"
  vacc_program <- "scenario1" # scenario number (1-9)
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
source("MenA_paramCheck_COVID.R")
source("ModelInputUtilities_COVID.R")
source("MenA_OneSim_COVID.R")
source("MenA_helper_functions_COVID.R")
source("MenA_summarization_functions.R")
source("MenA_calibration_plots.R")

# (B) Check parameters set above
setparams<-as.list(c(mycountry, as.character(start), as.character(end), myregion, PSA, vacc_program, phi, seed, nSims, input.dir, output.dir))
names(setparams)<-c("mycountry", "start", "end", "myregion", "PSA", "vacc_program", "phi", "seed", "nSims", "input.dir", "output.dir")
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
  nodoses<-as.vector(myvacc[is.na(myvacc$DosesCampaign) | myvacc$DosesCampaign==0,"year"])
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
# Set the random number seed based on seed, then create a vector of random number seeds (consistent within seed values)
set.seed(seed, kind = NULL, normal.kind = NULL)
seed.vec <- unique(floor(runif(nSims*2, 0, 1000000)))[1:nSims]

# Begin simulations
cl <- makeCluster(3)  #scale this upwards if you're on a workstation with >16gb memory
registerDoParallel(cl)
my_data <- foreach(n=1:nSims, .verbose=TRUE, .packages = c("lubridate", "dplyr", "data.table", "reshape2")) %dopar% {
  set.seed(seed.vec[n])
  paramfixed.row <- paramfixed[n,]  
  finalpop<-MenASimulation(startdt=start, enddt=end, fp=paramfixed.row, initpop=initpop, vacc_program=vacc_program,
                           countryparams=myparams, region=myregion, country=mycountry, inputdir=input.dir)
  if (summarizeme > 0) {
    # Using PSA option for CFR
    cfr <- as.numeric(paramfixed.row[, c("cfr1", "cfr2", "cfr3", "cfr4", "cfr5", "cfr6")])
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


### (5) Create outputs ########################################################

# (A) Detail output for testing if needed
detail<-1
if (use.tensims==TRUE) {
  tensims<-rbindlist(my_data[1:10])
  tensimsum<-tensims%>%group_by(simulation, IterYear)%>%summarize(sumCases=sum(Cases))
  sfile = paste0(mycountry, "_tensims_", vacc_program, "_", Sys.Date(), ".csv")
  detfile<-paste0(output.dir, "/", sfile)
  write.csv(tensimsum, detfile)
  print(paste("Simulation detail written to", detfile))
}

# (B) Output summary results for the country/scenario set
filename <- paste0(mycountry, "_", vacc_program, "_", Sys.Date(), ".csv")
filename1<- paste0(output.dir, "/", filename)
finalsummary<-summarizeForOutput(my_data, cohort=cohortSize, write=TRUE, filename=filename1)


# (C) Output PSA results if needed
if (PSA==TRUE){
  filename.psa <- paste0("PSA", "_",  mycountry, "_", vacc_program, "_", 
                    Sys.Date(), ".csv")
  filename.psa1 <- paste0(output.dir, "/", filename.psa)
  
  psa.output <- data.frame(disease=character(0), run_id=numeric(0), year=numeric(0),
                           age=numeric(0), country=character(0), 
                           country_name=character(0), cohort_size=numeric(0),
                           cases=numeric(0), dalys=numeric(0), deaths=numeric(0),
                           stringsAsFactors = FALSE)
  
  names.df <- data.frame(
    country_code=c("BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "ERI", "ETH", "GHA", "GIN", "GMB",
                   "GNB", "KEN", "MLI", "MRT", "NER", "NGA", "RWA", "SDN", "SEN", "SSD", "TCD", "TGO",
                   "TZA", "UGA"),
    country=c("Burundi", "Benin", "Burkina Faso", "Central African Republic",
              "Cote d'Ivoire", "Cameroon", "Congo, the Democratic Republic of the",
              "Eritrea", "Ethiopia", "Ghana", "Guinea", "Gambia", "Guinea-Bissau", "Kenya",                                
              "Mali", "Mauritania", "Niger", "Nigeria", "Rwanda", "Sudan", "Senegal",                              
              "South Sudan", "Chad", "Togo", "Tanzania, United Republic of", "Uganda" ))
  
  oldnames=c("AgeInYears", "Cases", "Deaths", "DALYs", "cohortsize", "simulation")
  newnames=c("age", "cases", "deaths", "dalys", "cohort_size", "run_id")
  
  for (s in 1:nSims){
    result.df <- left_join(x=my_data[[s]], y=cohortSize, by=c("year", "AgeInYears"))
    result.df <- result.df %>% rename_at(vars(oldnames), ~newnames)
    psa.output <- bind_rows(psa.output, result.df)
  }
  records <- length(psa.output$disease)
  psa.output$disease <- rep("MenA", times=records)
  psa.output$country <- rep(mycountry, times=records)
  psa.output$country_name <- rep(names.df$country[names.df$country_code==mycountry], times=records)
  
  write.csv(psa.output, filename.psa1, row.names=FALSE)
  print(paste("Simulation detail written to", filename.psa1))
}

# Update scenario tracker
if (automate==TRUE){
  scenario_tracker$completed[scen_num] <- "TRUE"
  write.csv(scenario_tracker, paste(input.dir, "scenario_tracker.csv", sep="/"),
            row.names=FALSE)
}


print(begin)
print(Sys.time())
