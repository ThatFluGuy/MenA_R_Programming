#### Program information ######################################################
# Source file name: MenA_helper_functions.R                                   #
# Package: MenA_VaccSims                                                      #
# Contact: chris.c.stewart@kp.org, michael.l.jackson@kp.org                   #
# Version Date 12/13/18                                                       #
#_____________________________________________________________________________#
#_____________________________________________________________________________#
# Functions called by MenA_OneSim and calling scripts                         #
# Contents:                                                                   #
# -InitializePopulation: Used by calling scripts before running simulations  #
# -GetInfectiveRatio: Called by MenA_OneSim in calculating force ofinfection #
#`-Vaccinate: Called by MenA_OneSim according to which program is specified  #
#_____________________________________________________________________________#


#_____________________________________________________________________________#
# INITIALIZE POPULATION Input datasets: PopInputAgeDist.csv, dist_both.csv    #
#_____________________________________________________________________________#
# Parameters:                                                                 #
# Start and end dates of simulation - to be specified in calling program      #
# Startpop from country params                                                #
#	country name to get age distributions from popinputagedist.csv              #
#	region = hyper/not hyper to get age-specific disease state distribution.    #
#_____________________________________________________________________________#
# Purpose: fill first "slice" of the simulation matrix with the population at #   
# the starting point, divided by age and disease state. Matrix is 360 monthly #
# ages by 8 disease states.                                                   #
#_____________________________________________________________________________#
# Created as function 3/5/18, by Chris Stewart chris.c.stewart@kp.org         #
# Changes: remove loop, use vector multiplication instead,                    #
# 30 age group only needs one value in input vectors                          #           
#_____________________________________________________________________________#
InitializePopulation<-function(scriptdir, inputdir, start, end, popsize, country="ETH", region="not_hyper") {
  initmsg<<-""
  #create the matrix
  dur<- ceiling(difftime(end, start, units = "weeks"))
  # poparray <- array(data=0, dim=c(361, 9, dur+1)) # Dimensions are [age groups, states, time]
  # Chloe edit 3/22/19: allow for monthly ages through age 120, including age=0 months.
  # Later on, for sake of comparison with Montagu estimates, will only look through age 70.
  # Note this could make code extremely (and unecessarily) slow...maybe just through age 70 if that's the case?
  poparray <- array(data=0, dim=c(1441, 9, dur+1)) # Dimensions are [age groups, states, time]
  dimnames(poparray)[[2]] <- c("Ns","Nc","Ls","Lc","Hs","Hc","Di","Va", "Inc")
  #parameters for initializing
  mypop<-GetPopAgeDist(path=inputdir, mycountry=country, start=start) 
  # Chloe edit 3/22/19: Now, this should start with pop size of age bands between 0-4 and 100-120.
  # Chloe edit 3/29/10: The error message below was originally intended for 7 age bands; 
  # I'm leaving this unedited for now, since I'm not sure what final groups of age bands we will use. 
  if (length(mypop) < 7) {
    initmsg<<-"Incomplete age distribution information.  Please check the downloaded file [touchstone]_qq_pop_both.csv."
    return(FALSE)
  }
  ##get age-specific proportions of each disease state into vectors: 7 ages X 7 disease states
  ## Chloe 3/22/19: now, we have 21 ages and 7 disease states.
  ## The line below gets fraction of each disease state in each of 7 population groups (5-year age bands up to age 30).
  ## To avoid having too many edits in too many functions, I decided to leave the "GetDiseaseStateDist" function alone, and as before,
  # assume fractions in each disease state in the 30-year age band match those in the older age bands, i.e. can just repeat the
  # final 7 values as needed here to provide proportion estimates up to age 120...so, if limiting to age 70 above, will need to edit this.
  statefract<-GetDiseaseStateDist(directory=scriptdir, region=myregion)
  if (length(statefract) < 49) {
    initmsg<<-"Incomplete disease state distribution information.  Please check the file dist_both.csv provided with the scripts."
    return(FALSE)
  }
  #expand age group fraction as vector to match pop matrix dimension 1
  # Chloe edit 3/22/19: changing this to match new dims; see note above.
  # agefract <- c(rep(as.numeric(mypop[1:6]), each=60), as.numeric(mypop[7]))
  # Chloe 3/29: In new age bands, all are 5-year except for last, which is 120 years.
  agefract <- c(rep(as.numeric(mypop[1:20]), each=60), rep(as.numeric(mypop[21]), each=(12*20)+1))
  # chunks <- c(rep(60, each=360), 1)
  # Chloe: Looks as though the starting assumption is that population age bands are evenly divided 
  # among each included age.
  chunks <- c(rep(60, each=1200), rep((12*20)+1,each=(12*20)+1))
  
  # Exanding this to more age brackets.
  # statemx<-rbind(
  #  matrix(rep(statefract[1:7], each=60), nrow=60),
  #  matrix(rep(statefract[8:14], each=60), nrow=60),
  #  matrix(rep(statefract[15:21], each=60), nrow=60),
  #  matrix(rep(statefract[22:28], each=60), nrow=60),
  #  matrix(rep(statefract[29:35], each=60), nrow=60),
  #  matrix(rep(statefract[36:42], each=60), nrow=60),
  #  matrix(rep(statefract[43:49], each=1), nrow=1)
  # )
  
  # Assuming last 7 values of statefract (originally just for 30-34 age band) apply to all older ages (which was effectively 
  # implied before).
  # Separate matrix for each age chunk,
  # separate column for disease state, separate row for each month of age.
   statemx<-rbind(
    matrix(rep(statefract[1:7], each=60), nrow=60),
    matrix(rep(statefract[8:14], each=60), nrow=60),
    matrix(rep(statefract[15:21], each=60), nrow=60),
    matrix(rep(statefract[22:28], each=60), nrow=60),
    matrix(rep(statefract[29:35], each=60), nrow=60),
    matrix(rep(statefract[36:42], each=60), nrow=60),
    matrix(rep(statefract[43:49], each=60), nrow=60),
    matrix(rep(statefract[43:49], each=60), nrow=60),
    matrix(rep(statefract[43:49], each=60), nrow=60),
    matrix(rep(statefract[43:49], each=60), nrow=60),
    matrix(rep(statefract[43:49], each=60), nrow=60),
    matrix(rep(statefract[43:49], each=60), nrow=60),
    matrix(rep(statefract[43:49], each=60), nrow=60),
    matrix(rep(statefract[43:49], each=60), nrow=60),
    matrix(rep(statefract[43:49], each=60), nrow=60),
    matrix(rep(statefract[43:49], each=60), nrow=60),
    matrix(rep(statefract[43:49], each=60), nrow=60),
    matrix(rep(statefract[43:49], each=60), nrow=60),
    matrix(rep(statefract[43:49], each=60), nrow=60),
    matrix(rep(statefract[43:49], each=60), nrow=60),
    matrix(rep(statefract[43:49], each=(12*20)+1), nrow=(12*20)+1)
   )
  #initialize - 3rd dimension stays at 1
  poparray[, "Ns", 1]<-(startSize*agefract/chunks)*statemx[,1]
  poparray[, "Nc", 1]<-(startSize*agefract/chunks)*statemx[,2]
  poparray[, "Ls", 1]<-(startSize*agefract/chunks)*statemx[,3]
  poparray[, "Lc", 1]<-(startSize*agefract/chunks)*statemx[,4]
  poparray[, "Hs", 1]<-(startSize*agefract/chunks)*statemx[,5]
  poparray[, "Hc", 1]<-(startSize*agefract/chunks)*statemx[,6]
  poparray[, "Di", 1]<-(startSize*agefract/chunks)*statemx[,7]
  #leave Vacc, Inc 0
  return(poparray)
}

#_____________________________________________________________________________#
# GET INFECTIVE RATIO                                                         #
#_____________________________________________________________________________#
# Parameters:  inpop = current slice of population array                      # 
#_____________________________________________________________________________#
# Purpose: four-element vector containing the ratio of infectives             #
# in 4 age groups(0-<5, 5-<13, 14-<20, 20+)                                   #
#_____________________________________________________________________________#
GetInfectiveRatio<-function(inpop){
  #sum total pop 
  mosums<-colSums(inpop)
  N1<-(sum(mosums[1:60]))
  N2<-(sum(mosums[61:156]))
  N3<-(sum(mosums[157:240]))
  N4<-(sum(mosums[241:361]))
  #sum infectious pops
  infpop<-inpop[c("Nc","Lc","Hc","Di"),]
  infmosums<-colSums(infpop)
  I1<-(sum(infmosums[1:60]))
  I2<-(sum(infmosums[61:156]))
  I3<-(sum(infmosums[157:240]))
  I4<-(sum(infmosums[241:361]))
  #calc ratios
  infRatio<-c(I1/N1, I2/N2, I3/N3, I4/N4)
  return(infRatio)
}

#_____________________________________________________________________________#
# VACCINATE                                                                   #
#_____________________________________________________________________________#
# Parameters: popslice = current slice of population array                    # 
#             vlookup = vaccination scenario from VIMC input files            #
#             type = vaccination scenario                                     #
#             mydate=current sim date, to supply year of vaccination scenario #
#_____________________________________________________________________________#
# Purpose: moves people into vaccinated state per vacc. input file info       #
#_____________________________________________________________________________#

vaccinate<-function(popslice, vlookup, type, mydate) {
  #for type="both", both ifs should execute
  eligibles<-c("Ns","Nc","Ls","Lc","Hs","Hc")
  if ((type=="campaign" | type=="both") & (month(mydate)==10)) {
    #get parameters
    cDoses <- vlookup[vlookup$year==year(mydate) & vlookup$activity_type=="campaign","DosesCampaign"]
    #zero-length cDoses (not NA, apparently) is blowing things up
    if (length(cDoses)> 0) {
      if (!is.na(cDoses)) {
       # print(cDoses)
        #change ages i=12 to 359, here 13 to 360
        eligN <- sum(popslice[13:360, eligibles]) 
        pcNLH <- ifelse(cDoses<=eligN, cDoses/eligN, 1 )
        #why do the above if were not going to move anybody?
        popslice[13:360,"Va"] <- popslice[13:360,"Va"] + pcNLH * rowSums(popslice[13:360,eligibles])
        popslice[13:360,eligibles] <- (1-pcNLH) * popslice[13:360, eligibles]
      }
    }
  }
  if (type=="routine" | type=="both") {
    pr <- vlookup[vlookup$year==year(mydate) & vlookup$activity_type=="routine","CoverRoutine"]
   # print(pr)
    if (length(pr)> 0 & pr>0) {
      if (!is.na(pr)) {
        #routine vaccinations, when subjects turn 9 months- thats 10 here
        popslice[10,"Va"]<- popslice[10,"Va"] + (pr*sum(popslice[10,eligibles]))
        popslice[10,eligibles]<- (1-pr)*popslice[10,eligibles]
      }
    }
  }
  return(popslice)
}