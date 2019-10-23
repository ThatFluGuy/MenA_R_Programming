#### Program information ######################################################
# Package: MenA_VaccSims                                                      #
# Source file name: MenA_OneSim.R                                             #
# Contact: chris.c.stewart@kp.org, michael.l.jackson@kp.org                   #
# Version Date 12/31/2018                                                     #
#_____________________________________________________________________________#
# Input datasets: none.                                                       #
#_____________________________________________________________________________#
# Parameters:                                                                 #
# Start and end dates of simulation - to be specified in calling program      #
# population array to fill                                                    #
#	data frame containing country parameters                                    #
#	WAIFW matrix (3d with dry vs rainy?)                                        #
# dd1&2, dr1&2 per-carrier risk of invasive disease (age-specific)            #
#_____________________________________________________________________________#
# Purpose: One simulation - iterate weekly from start date to end date        #   
# Move people per WAIFW weekly, age monthly, vaccinate monthly as speficiednnn#
# fill pop matrix slices 2 thru n                                             #
#_____________________________________________________________________________#
# Created as function 3/6/18, by Chris Stewart chric.c.stewart@kp.org         #
# Changes:                                                                    #
#   WAIFW and dxrisk as small 3D matrices                                     #           
#_____________________________________________________________________________#

# Chloe note 3/29: input to function from sims code (just so I can name variables while working with function below)
# startdt=start, enddt=end, pop=initpop,
# fixedparams=paramfixed, countryparams=myparams,
# WAIFWmx=wboth, dxr=dxrisk, vacc_program=vacc_program

# Chloe, 3/29: the WAIFW matrix does appear to be age specific, so I'll need to expand it to be
# appropriate for the older age groups.

# Eric, 10/23: Lots of changes.  Moved dxrisk, waifw, population setup to inside function.  Split vaccinated status into Vs, Vc.  
#              Many calls to fp, the fixed parameters generated via model calibration.  Changed syntax from fp["rc"] to fp$rc to
#              ensure that the result is a vector.  I do not think any hard-coded parameters remain in this code or helper functions.
MenASimulation<-function(startdt, enddt, fp, vacc_program, countryparams, region, country, inputdir) { 
  #setup before loop
  #disease model for rainy and dry
  dxrisk<-rbind(c(fp$dr1, fp$dr2),
                c(fp$dd1, fp$dd2))
  dimnames(dxrisk)[[1]]<-c("rainy", "dry")
  
  #WAIFW matrix setup
  wboth<-GetWAIFWmatrix(params=fp)
  if (!(is.numeric(wboth))) {
    stop("WAIFW matrix is not numeric")
  }
  
  #initialize population
  startSize <- countryparams[countryparams$year==year(start)-1, "totalpop"]
  pop<-InitializePopulation(scriptdir=script.dir, inputdir=inputdir, start=startdt, end=enddt, country=country, region=region, startSize=startSize)
  #check for errors
  if (!(is.numeric(pop))) {
    if (disterr!="") { print(disterr) } 
    if (dxerr!="") { print(dxerr) } 
    stop(initmsg)
  }
  
  #setup before loop
  theDate <- start
  births <- countryparams[countryparams$year==year(theDate), "births"]/52.1775
  imr <- countryparams[countryparams$year==year(theDate), "imr"]/(1000*52.1775)
  ind1 <- which(colnames(countryparams)=="dr59")
  ind2 <- which(colnames(countryparams)=="dr7579")
  ages5through79 <- countryparams[countryparams$year==year(theDate), ind1:ind2]/(1000*52.1775)
  ages1through4 <- countryparams[countryparams$year==year(theDate), "dr14"]/(1000*52.1775)
  over80 <- countryparams[countryparams$year==year(theDate), "dr8084"]/(1000*52.1775)
  deathvec <- c(rep(imr,12),rep(ages1through4,(12*4)),unlist(rep(ages5through79,each=(12*5))),rep(over80,(40*12)+1))
  # v <- countryparams[countryparams$year==year(theDate), "v"] / (1000*52.1775)
  # deathvec<-c(rep(imr,12), rep(v, 349))
  # deathvec<-c(rep(imr,12), rep(v, 1429))
  # Chloe: this is where new death rates could be included; at the moment,
  # has a single death rate for infants and a single death rate for everyone through age 30 (and beyond)
  # Could start by expanding here.
  # For now, I'll need to change how often the final value is repeated.
  # waning age-group vectors 0-5mo, 6mo-2y, 3-10y, 11+y
  # Chloe 5/15: age-dependent rates of waning between stages of immunity;
  # expanding the final group to all individuals up to age 120.
  
  # wanev <- c(rep(1-imr, 7), rep(0.000172, 17), rep(0.000096, 107), rep(0.000307, 230))  #waning from vacc to hi ab scaled to wks-confirm wv(1) = NA
  # waneh <- c(rep(0.01092, 6), rep(0.00654, 18), rep(0.00527, 107), rep(0.00096, 230)) #waning from high to low ab scaled to weeks
  # wanel <- c(rep(0.00970,6), rep(0.00487, 18), rep(0.00364, 107), rep(0.00057, 230)) #waning from low 
  

  #wanev <- c(rep(1-imr, 7), rep(0.000172, 17), rep(0.000096, 107), rep(0.000307, 1310)) #Chloe's edits
  #waneh <- c(rep(0.01092, 6), rep(0.00654, 18), rep(0.00527, 107), rep(0.00096, 1310))
  #wanel <- c(rep(0.00970,6), rep(0.00487, 18), rep(0.00364, 107), rep(0.00057, 1310))

  #EJ: replace hard-coded parameters with references to fp, the fixed parameter input to OneSim.
  wanev <- c(rep((1-imr)*52.1775, 7), rep(fp$wvh2, 17), rep(fp$wvh3, 107), rep(fp$whv4, 1310))/52.1775  #waning from vacc to hi ab scaled to wks-confirm wv(1) = NA     This section is weird in the earliest age group.
  waneh <- c(rep(fp$whl1, 6), rep(fp$whl2, 18), rep(fp$whl3, 107), rep(fp$whl4, 1310))/52.1775 #waning from high to low ab scaled to weeks.  Also, check if dividing by 52.1775 is valid.
  wanel <- c(rep(fp$wln1,6), rep(fp$wln1, 18), rep(fp$wln1, 107), rep(fp$wln1, 1310))/52.1775 #waning from low 
  
  j <- 2 #initialize time interval counter - the first position is filled with starting pop
  iterDate<-start-7 #this will make aging happen the first iteration, like SAS - also supply date for initial pop
  LastMonth <- month(iterDate)
  iterYear<-year(start-7)
  random<-runif(1)
  sto <- 1 + (phi * cos(random*pi))
  #date loop
  while (theDate <= enddt) {
    #store the date in a vector to label iterations later
    #try naming them now, then maybe the following won't be necessary / doesn't seem to work
    #dimnames(pop)[[3]][[j]] <- as.character(theDate)
    iterDate<-c(iterDate, theDate)
    iterYear<-c(iterYear, year(theDate))
    if (day(theDate) <= 7 & month(theDate) == 9)  #yearly update to stochastic param
    {
      random<-runif(1)
      sto <- 1 + (phi * cos(random*pi))
      #sto<-1 #print(paste0("new stochastic parameter ", sto))
    }
    #yearly birth, death and imr rates, only needs to when year changes
    if (month(theDate)== 1 & LastMonth==12)
    {
      births <-countryparams[countryparams$year==year(theDate), "births"]/52.1775
      # imr <- countryparams[countryparams$year==year(theDate), "imr"]/(1000*52.1775)
      # v <- countryparams[countryparams$year==year(theDate), "v"] / (1000*52.1775)
      #new imr, need to update wanev and death
      # wanev <- c(rep(1-imr, 7), wanev[8:361]) 
      # Chloe 5/15: Update to imr changes wanev first few items, rest stays the same.
      imr <- countryparams[countryparams$year==year(theDate), "imr"]/(1000*52.1775)
       wanev <- c(rep(1-imr, 7), wanev[8:1441]) 
      # deathvec<-c(rep(imr,12), rep(v, 349))
      # Chloe 5/22: same as above for each year of sim.
      ind1 <- which(colnames(countryparams)=="dr59")
      ind2 <- which(colnames(countryparams)=="dr7579")
      ages5through79 <- countryparams[countryparams$year==year(theDate), ind1:ind2]/(1000*52.1775)
      ages1through4 <- countryparams[countryparams$year==year(theDate), "dr14"]/(1000*52.1775)
      over80 <- countryparams[countryparams$year==year(theDate), "dr8084"]/(1000*52.1775)
      deathvec <- c(rep(imr,12),rep(ages1through4,(12*4)),unlist(rep(ages5through79,each=(12*5))),rep(over80,(40*12)+1))
    }
    
    #  waifw matrix depends on rainy (Mar-Aug) or dry (Sep-Feb) season
    if (month(theDate)!=LastMonth) {
      if (month(theDate) %in% c(3,4,5,6,7,8)) { 
        wmx <- wboth[,,"rainy"]
      } 
      else { wmx <- wboth[,,"dry"]}
      # Determine age-specific per-carrier risk of invasive disease;
      #like infection, this is also seasonal but peaks later, hence different months (rainy here = 6-10 June-Oct)
      # Chloe 5/30: don't follow how this works for now, but replacing with values similar to original for now;
      # 6/10: confirmed with Mike this should hold.
      # assuming the age-specific per-carrier risk of invasive disease is the same for anyone over 30.
      # iagevec<-seq(0, 360)
      iagevec<-c(seq(0, 360),rep(360,1441-361))
      if (month(theDate) %in% c(1,2,3,4,5,11,12)) {
        sigma = sigmavec<-dxrisk["dry",1] + dxrisk["dry",2] *iagevec
      } else {sigmavec<-dxrisk["rainy",1] + dxrisk["rainy",2] *iagevec }
    } #end of updates conditional on month
    
    #infectious ratios by pop age group (for calculating force) this happens every time point
    IR<-GetInfectiveRatio(t(pop[,,j-1]))
    # Includes force of infection from outside the population (foii, was for in SAS code)
    forcevec <- (sto*wmx[,1]*IR[1]) + (sto*wmx[,2]*IR[2]) + (sto*wmx[,3]*IR[3])+ (sto*wmx[,4]*IR[4]) + fp$foii
    
    
    
    #Transitions from previous time period to current time period
    pop[,"Ns",j] <- wanel*pop[,"Ls",j-1] - (deathvec+forcevec)*pop[,"Ns",j-1] + pop[,"Ns",j-1]
    #births - this order is the same as sas (could just add births to line above)
    pop[1,"Ns",j] <- pop[1,"Ns",j] + births
    pop[,"Nc",j] <- forcevec*pop[,"Ns",j-1] - (deathvec+fp$rc+sigmavec) * pop[,"Nc",j-1] + 
      pop[,"Nc",j-1]
    pop[,"Ls",j] <- waneh*pop[,"Hs",j-1] + fp$rc*pop[,"Nc",j-1] - 
      (deathvec+wanel+(1-fp$lc)*forcevec) * pop[,"Ls",j-1] + pop[,"Ls",j-1]
    pop[,"Lc",j] <- (1-fp$lc)*forcevec*pop[,"Ls",j-1] - 
      (deathvec+fp$rc+(1-fp$ld)*sigmavec) * pop[,"Lc",j-1]  + pop[,"Lc",j-1]
    pop[,"Hs",j] <- fp$rc*(pop[,"Lc",j-1] + pop[,"Hc",j-1]) + 
      fp$rd*pop[,"Di",j-1]  + wanev*pop[,"Vs",j-1] - (deathvec+waneh+(1-fp$hc)*forcevec) * pop[,"Hs",j-1] + pop[,"Hs",j-1] #waning from Vs instead of Va
    pop[,"Hc",j] <- (1-fp$hc)*forcevec*pop[,"Hs",j-1] - 
      (deathvec+fp$rc-(1-fp$hd)*sigmavec)*pop[,"Hc",j-1] + wanev*pop[,"Vc",j-1] + pop[,"Hc",j-1] #added waning from Vc
    pop[,"Di",j] <- sigmavec*(pop[,"Nc",j-1] + (1-fp$ld)*pop[,"Lc",j-1] + 
      (1-fp$hd)*pop[,"Hc",j-1]) - (deathvec+fp$rd)*pop[,"Di",j-1] + pop[,"Di",j-1]
    pop[,"Vs",j] <- pop[,"Vs",j-1]  - wanev*pop[,"Vs",j-1] - deathvec*pop[,"Vs",j-1] +  fp$rc*pop[,"Vc",j-1]  #ways of moving out of Vs: wane to Hs, death.  Moving in requires the vaccinate function, or recovering from Vc
    pop[,"Vc",j] <- pop[,"Vc",j-1] -  fp$rc*pop[,"Vc",j-1] - wanev*pop[,"Vc",j-1] - deathvec*pop[,"Vc",j-1] #ways of moving out of Vc: wane to Hc, recover from Vc to Vs, death.  Moving in requires the vaccinate function
    # Count incident cases of invasive disease NOTE THESE ARE before incrementing = use old pop
    pop[,"Inc",j] <- sigmavec*(pop[,"Nc",j-1] + (1-fp$ld)*pop[,"Lc",j-1] + (1-fp$hd)*pop[,"Hc",j-1])
    #Inc[,j] <- sigmavec*(pop[,"Nc",j-1] + (1-ld)*pop[,"Lc",j-1] + (1-hd)*pop[,"Hc",j-1])
    
    #aging handled monthly ** THIS WORKS **aging incident cases too
    if (month(theDate) != LastMonth){
      for (x in 1:9) {
        #if (x==1) { beforeaging<-pop[,x,j]}
        #save old 361 spot to add to last pot later
        # Chloe 5/30: now, in over 120 pot.
        # save361<-pop[[361,x,j]]
        save1441<-pop[[1441,x,j]]
        #vector to shift: pop[m,x,j]
        pop[,x,j]<-(c(0,  pop[,x,j])[1 : length(pop[,x,j])])
        #move last monthly group into 30+ pot
        # pop[[361,x,j]]<-pop[[361,x,j]] + save361
        pop[[1441,x,j]]<-pop[[1441,x,j]] + save1441
      } #end of x loop
      #Vaccinate here, monthly
      if (vacc_program!="none") {
        #there are some cases where theres nothing to do in many years
        if  (vacc_program=="campaign") { 
          if (year(theDate) %in% nodoses==FALSE) {
            pop[,,j] <- vaccinate(popslice=pop[,,j], vlookup=myvacc, type=vacc_program, mydate=theDate, params=fp)
            #print("vaccinating")
          }
        }
        else {
          pop[,,j] <- vaccinate(popslice=pop[,,j], vlookup=myvacc, type=vacc_program, mydate=theDate, params=fp)
          #print("vaccinating")
        }
      }
    } #end of aging IF
    # if (year(theDate) > 2020) {break}
    LastMonth <- month(theDate)
    theDate<-theDate+7
    j <- j + 1
  } #end of date loop
  #label matrix iteration dimension with dates
  dimnames(pop)[[3]] <- as.Date(iterDate)
  return(pop)
}