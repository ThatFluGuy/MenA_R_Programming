#### Program information ######################################################
# Source file name: fxSimulationSubFunctions.R                                #
# Package: MenA_VaccSims                                                      #
# Contact: chris.c.stewart@kp.org, michael.l.jackson@kp.org                   #
# Version Date 12/13/18                                                       #
#_____________________________________________________________________________#
#_____________________________________________________________________________#
# Functions called by MenA_OneSim and calling scripts                         #
# Contents:                                                                   #
# -GetInfectiveRatio: Called by MenA_OneSim in calculating force ofinfection #
#`-Vaccinate: Called by MenA_OneSim according to which program is specified  #
#_____________________________________________________________________________#




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
# EJ 10/23: split out into Vs and Vc.  Added params argument to adjust for    #
#           vaccine effectiveness                                             #
#_____________________________________________________________________________#

vaccinate<-function(popslice, vlookup, type, mydate, params) { 

  # For type="both", both ifs should execute
  eligibles.s <- c("Ns", "Ls", "Hs")
  eligibles.c <- c("Nc", "Lc", "Hc")
  if ((month(mydate)==10)) { # Run campaigns (if any) during October only
    # Get parameters
    cDoses <- vlookup[vlookup$year==year(mydate) & vlookup$activity_type=="campaign","DosesCampaign"]
    # Zero-length cDoses (not NA, apparently) is blowing things up
    if (length(cDoses)> 0) {
      if (!is.na(cDoses)) {
        # print(cDoses)
        #change ages i=12 to 359, here 13 to 360
        #EJ 10/23: separating out susceptibles and colonized, adding vaccine effectiveness 
        min.age <- vlookup$age_first[vlookup$year==year(mydate) & vlookup$activity_type=="campaign"]
        max.age <- vlookup$age_last[vlookup$year==year(mydate) & vlookup$activity_type=="campaign"]
        eligN.s <- sum(popslice[min.age:max.age, eligibles.s])
        eligN.c <- sum(popslice[min.age:max.age, eligibles.c]) 
        pcNLH <- ifelse(cDoses<=(eligN.s + eligN.c), cDoses/(eligN.s + eligN.c), 1 )
        #why do the above if were not going to move anybody?  EJ: it allows for fractional vaccination of the eligible population, if there aren't enough doses available.
        #formula change, example:
        #100 people in eligible, 50% pcNLH, 90% ve.  100*.5*.9 = 45 vaccinated
        #100*(1 - .9*.5) = 100*(1-.45)=55 people remain at risk.  
        popslice[min.age:max.age,"Vs"] <- popslice[min.age:max.age,"Vs"] + params$ve * pcNLH * rowSums(popslice[min.age:max.age,eligibles.s])
        popslice[min.age:max.age,eligibles.s] <- (1 - params$ve * pcNLH) * popslice[min.age:max.age, eligibles.s]
        popslice[min.age:max.age,"Vc"] <- popslice[min.age:max.age,"Vc"] + params$ve * pcNLH * rowSums(popslice[min.age:max.age,eligibles.c])
        popslice[min.age:max.age,eligibles.c] <- (1 - params$ve * pcNLH) * popslice[min.age:max.age, eligibles.c]
              }
    }
  }
  
  if (type %in% c("routine", "both", "booster")) { # Don't run for campaign-only scenarios
    pr <- vlookup[vlookup$year==year(mydate) & vlookup$activity_type=="routine","CoverRoutine"]
    # print(pr)
    if (length(pr)> 0) if(pr>0) {
      if (!is.na(pr)) {
        #routine vaccinations, when subjects turn 9 months- thats 10 here
        popslice[10,"Vs"]<- popslice[10,"Vs"] + (pr*sum(popslice[10,eligibles.s]))
        popslice[10,eligibles.s]<- (1-pr)*popslice[10,eligibles.s]
        popslice[10,"Vc"]<- popslice[10,"Vc"] + (pr*sum(popslice[10,eligibles.c]))
        popslice[10,eligibles.c]<- (1-pr)*popslice[10,eligibles.c]
        
      }
    }
  }
  return(popslice)
}