#### Program information ######################################################
# Source file name: MenA_helper_functions.R                                   #
#_____________________________________________________________________________#
# Input datasets: none                                                        #
#_____________________________________________________________________________#
# Functions called by MenA_OneSim                                             #
# Contents:                                                                   #
# Get Infective Ratio: Called by MenA_OneSim in calculating force of          #
#	    infection                                                               #
#`Vaccinate: Called by MenA_OneSim according to which program is specified    #
#_____________________________________________________________________________#
#_____________________________________________________________________________#
# Created as function 3/7/18, by Chris Stewart chris.c.stewart@KP.org         #
# Changes:                                                                    #
#_____________________________________________________________________________#


GetInfectiveRatio<-function(inpop){
  #returns a four-element vector containing the ratio of infectives in 4 age groups (0-<5, 5-<13, 14-<20, 20+)
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
  infRatio<-c(I1/N1, I2/N2, I3/N3, I4/N4)
  return(infRatio)
}

vaccinate<-function(popslice, vlookup, type, mydate) {
  #for type="both", both ifs should execute
  eligibles<-c("Ns","Nc","Ls","Lc","Hs","Hc")
  if ((type=="campaign" | type=="both") & (month(mydate)==10)) {
    #get parameters
    cDoses <- vlookup[vlookup$year==year(mydate) & vlookup$activity_type=="campaign","DosesCampaign"]
    #zero-length cDoses (not NA, apparently) is blowing things up
    if (length(cDoses)> 0) {
      if (!is.na(cDoses)) {
      #  print(cDoses)
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
    #print(pr)
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