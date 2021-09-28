#### Program information ######################################################
# Package: MenA_VaccSims                                                      #
# Source file name: ConstantFOI_VaccSims.R                                    #
# Version Date 09/24/2021                                                     #
#_____________________________________________________________________________#
# Purpose: Run a single interation of the simulation, using modified          #
# paramater sets so that FOI is constant.                                     #
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


### (1) Set up libraries and directories ######################################

library(lubridate)
library(doParallel)
library(dplyr)
library(data.table)
library(reshape2)
library(tidyr)

# Directory containing inputs from https://montagu.vaccineimpact.org/
input.dir<-"G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Data/GAVI inputs/201910gavi_v4"
# Directory for central simulation outputs
output.dir <- "G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Analysis/ConstantFOI"
# Directory for final PSA outputs
deliv.dir <- "G:/CTRHS/Modeling_Infections/GAVI MenA predictions/Deliverables/Deliverables 2021"
# Directory containing R scripts
script.dir <- "C:/Users/O992928/documents/GAVI MenA predictions/R_programming"


### (2) Set up high-level options #############################################

mycountry <- "TCD"
myregion <- "hyper"
vacc_program <- "routine"
vacc_subprogram <- "default"

### (3) Set up cross-scenario options #########################################

# (A) Universal variables regardless of scenario
start <- as.Date("1951-01-01") # Use 1/1/1951 to give 50 years burn-in
end <- as.Date("2100-12-31")
PSA <- FALSE
seed <- 4567 # Seed for random sto, use same for all scenarios
use.tensims <- FALSE # If desired, output results from 10 sims for debugging
nSims <- 200
phi<-0  # Set stochastic parameter


# (B) Import scripts
if (dir.exists(script.dir)) {
  if (file.exists(paste0(script.dir, "/","fxParamChecks.R"))==FALSE) {
    msg<-(paste0("MenA_paramcheck.R not found in ", script.dir))
    stop(msg)
  }
} else {
  script.dir<- getSrcDirectory(function(dummy) {dummy})
  if (file.exists(paste0(script.dir, "/","fxParamChecks.R"))==FALSE) {
    print(paste0("MenA_paramCheck.R not found in ", script.dir))
    stop("This script requires 6 other scripts; please put in same directory as this one or specify script directory.")
  }
}

setwd(script.dir)
source("fxParamChecks.R")
source("fxModelInputs.R")
source("fxMenASimulation.R")
source("fxSimulationSubFunctions.R")
source("fxSummarization.R")

# Variables and functions to keep between iterations when automate==TRUE
keep.v <- c("automate", "input.dir", "output.dir", "deliv.dir", "script.dir", "scenario_tracker",
                      "not_done", "scenario.loops", "start", "end", "PSA", "seed", "nSims", "use.tensims",
                      "phi", "keep.v", lsf.str())

### (4) Set-up and run simulation #############################################
# One iteration of this loop runs all nSims number of simulations for a       #
# specific country/scenario combination.                                      #

## (I) Set-up ##___________________________________________________________

print(paste("Country=", mycountry, "Program=", vacc_program, "Subprogram=", vacc_subprogram))

# (A) Import country-specific parameters vaccination program details.      #
setparams<-as.list(c(mycountry, as.character(start), as.character(end), myregion, PSA, vacc_program, phi, seed, nSims, input.dir, output.dir))
names(setparams)<-c("mycountry", "start", "end", "myregion", "PSA", "vacc_program", "phi", "seed", "nSims", "input.dir", "output.dir")
if (CheckSetParameters(setparams)==FALSE) {
    stop(spmessage)
} else {
  if (length(spmessage)>1) { print(spmessage) }
}

# (B) Import country-specific parameters
myparams.full <- GetDemographicParameters(path=input.dir,  mycountry=mycountry, start=start, end=end)
if (CheckDemogParameters(myparams.full)==FALSE) {
  stop(dpmessage)
} else {
  if (length(dpmessage)>1) { print(dpmessage) }
}

# (C) Import vaccination program details
if (vacc_program!="none") {
  myvacc<-GetVaccScenario(mycountry=mycountry, scenario=vacc_program, sub.scenario=vacc_subprogram, directory=input.dir)
  if (is.data.frame(myvacc)==FALSE) { stop(vaccmsg)}  #check for output
  #make as vector of years where nothing happens (empty except for campaign only) for efficiency
  if (vacc_program=="campaign") {
    nodoses<-as.vector(myvacc[is.na(myvacc$DosesCampaign) | myvacc$DosesCampaign==0,"year"])
  }
}

# (D) Country-specific life expectancy
my.lifex <- GetLifeExp(path=input.dir, mycountry.s=mycountry)

# (E) Read in parameters calculated in ABC, or a row of parameters to be used by ABC.  
paramfixed <- GetModelParams(path=script.dir, region.val=myregion)

# (F) Initialize full population.
startSize <- myparams.full[myparams.full$year==year(start)-1, "totalpop"]
initpop.full <- InitializePopulation(scriptdir=script.dir, inputdir=input.dir, start=start, end=end, country=mycountry, region=myregion, startSize=startSize)
#check for errors
if (!(is.numeric(initpop.full))) {
  if (disterr!="") { print(disterr) } 
  if (dxerr!="") { print(dxerr) } 
  stop(initmsg)
}

# (G) Scale down the full population to the modeled population
# This means changing both the starting population size and
# the annual number of births.
pct.modeled <- GetModelPct(path=input.dir, mycountry.s=mycountry)
initpop <- initpop.full
initpop[,,1] <- initpop[,,1] * pct.modeled

myparams <- myparams.full
myparams$births <- myparams$births * pct.modeled

## (II) Run simulations ##_________________________________________________
set.seed(seed)

onerun <- MenASimulation(startdt=start, enddt=end, fp=paramfixed[4,], initpop=initpop.full, vacc_program=vacc_program,
                         countryparams=myparams.full, region=myregion, country=mycountry, inputdir=input.dir)


## (III) Create outputs ##_________________________________________________

# Get incident cases by age group and time
inc.wide <- as.data.frame(onerun[, 10, ])
inc.wide$AgeMths <- 0:1440
names(inc.wide)[1:7828] <- paste0("Week", 1:7828)

# Convert to tall data frame
inc.df <- pivot_longer(data=inc.wide, cols=starts_with("Week"), names_to="Week",
                       values_to="Cases")

# Collapse to years, drop years before 2000
inc.df$Week.num <- as.numeric(gsub("Week", "", inc.df$Week))
inc.df$Date <- as.Date("1949-12-25") + (inc.df$Week.num * 7)
inc.df <- inc.df %>% filter(Date >= as.Date("2000-01-01"))
inc.df$year <- year(inc.df$Date)
inc.df2 <- aggregate(Cases~AgeMths+year, data=inc.df, FUN=sum)


# Collapse to age in yeasr
inc.df2$AgeInYears <- inc.df2$AgeMths %/% 12
inc.df3 <- aggregate(Cases~AgeInYears+year, data=inc.df2, FUN=sum)

# Get population size
pop.a <- apply(X=onerun[,1:9,], MARGIN=c(1,3), FUN=sum)
pop.wide <- as.data.frame(pop.a)
pop.wide$AgeMths <- 0:1440
names(pop.wide)[1:7828] <- paste0("Week", 1:7828)

# Convert to tall
pop.df <- pivot_longer(data=pop.wide, cols=starts_with("Week"), names_to="Week",
                       values_to="Pop")

# Collapse to years. Don't sum, but take the midpoint. Also drop if earlier than 2000
pop.df$Week.num <- as.numeric(gsub("Week", "", pop.df$Week))
pop.df$Date <- as.Date("1949-12-25") + (pop.df$Week.num * 7)
pop.df <- pop.df %>% 
  filter(month(pop.df$Date)==7 & day(pop.df$Date)<=7 & pop.df$Date >= as.Date("2000-01-01"))
pop.df$year <- year(pop.df$Date)

# Collapse to age in years
pop.df$AgeInYears <- pop.df$AgeMths %/% 12

pop.df2 <- aggregate(Pop~AgeInYears+year, data=pop.df, FUN=sum)

# Merge
output.df <- left_join(pop.df2, inc.df3, by=c("year", "AgeInYears"))

# Check overall incidence by years
by.year <- aggregate(cbind(Pop, Cases)~year, data=output.df, FUN=sum)
by.year$incidence <- 100000 * by.year$Cases / by.year$Pop

plot(by.year$year, by.year$incidence)

filename <- paste0("ConstantFOI-", vacc_program, ".csv")

write.csv(output.df, paste0(output.dir, "/", filename), row.names=FALSE)
