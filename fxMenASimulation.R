#### Program information ######################################################
# Package: MenA_VaccSims                                                      #
# Source file name: MenA_OneSim.R                                             #
# Contact: chris.c.stewart@kp.org, michael.l.jackson@kp.org                   #
# Version Date 12/17/2019                                                     #
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

MenASimulation<-function(startdt, enddt, fp, initpop, vacc_program, countryparams, region, country, inputdir) { 
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
  
  pop <- initpop
 
  #setup before loop
  theDate <- start
  births <- countryparams[countryparams$year==year(theDate), "births"]/52.1775
  imr <- countryparams[countryparams$year==year(theDate), "imr"]/(52.1775)
  ind1 <- which(colnames(countryparams)=="dr59")
  ind2 <- which(colnames(countryparams)=="dr7579")
  ages5through79 <- countryparams[countryparams$year==year(theDate), ind1:ind2]/(52.1775)
  ages1through4 <- countryparams[countryparams$year==year(theDate), "dr14"]/(52.1775)
  over80 <- countryparams[countryparams$year==year(theDate), "dr8084"]/(52.1775)
  deathvec <- c(rep(imr,12),rep(ages1through4,(12*4)),unlist(rep(ages5through79,each=(12*5))),rep(over80,(40*12)+1))

  wanev <- c(rep((1-imr), 7), rep(fp$wv2, 17), rep(fp$wv3, 107), rep(fp$wv4, 1310))  #waning from vacc to hi     This section is weird in the earliest age group.
  waneh <- c(rep(fp$wh1, 6),  rep(fp$wh2, 18), rep(fp$wh3, 107), rep(fp$wh4, 1310)) #waning from high to low
  wanel <- c(rep(fp$wl1,6),   rep(fp$wl2, 18), rep(fp$wl3, 107), rep(fp$wl4, 1310)) #waning from low 
  
  j <- 2 #initialize time interval counter - the first position is filled with starting pop
  iterDate<-start-7 #this will make aging happen the first iteration, like SAS - also supply date for initial pop
  LastMonth <- month(iterDate)
  iterYear<-year(start-7)
  random<-runif(1)
  sto <- 1 + (phi * cos(random*pi))
  #date loop
  while (theDate <= enddt) {
    #store the date in a vector to label iterations later
    iterDate<-c(iterDate, theDate)
    iterYear<-c(iterYear, year(theDate))
    if (day(theDate) <= 7 & month(theDate) == 9)  #yearly update to stochastic param
    {
      random<-runif(1)
      sto <- 1 + (phi * cos(random*pi))
    }
    #yearly birth, death and imr rates, only needs to when year changes
    if (month(theDate)== 1 & LastMonth==12)
    {
      births <-countryparams[countryparams$year==year(theDate), "births"]/52.1775
      imr <- countryparams[countryparams$year==year(theDate), "imr"]/(1000*52.1775)
      wanev <- c(rep(1-imr, 7), wanev[8:1441]) 
      ind1 <- which(colnames(countryparams)=="dr59")
      ind2 <- which(colnames(countryparams)=="dr7579")
      ages5through79 <- countryparams[countryparams$year==year(theDate), ind1:ind2]/(52.1775)
      ages1through4 <- countryparams[countryparams$year==year(theDate), "dr14"]/(52.1775)
      over80 <- countryparams[countryparams$year==year(theDate), "dr8084"]/(52.1775)
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
      # assuming the age-specific per-carrier risk of invasive disease is the same for anyone over 30.
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

    #aging handled monthly ** THIS WORKS **aging incident cases too
    if (month(theDate) != LastMonth){
      for (x in 1:9) {
        #save 1441 spot to add to last pot later
        save1441<-pop[[1441,x,j]]
        #vector to shift: pop[m,x,j]
        pop[,x,j]<-(c(0,  pop[,x,j])[1 : length(pop[,x,j])])
        #move last monthly group into 30+ pot
        pop[[1441,x,j]]<-pop[[1441,x,j]] + save1441
      } #end of x loop
      
      #Vaccinate here, monthly
      if (vacc_program!="none") {
        #there are some cases where theres nothing to do in many years
        if  (vacc_program=="campaign") { 
          if (year(theDate) %in% nodoses==FALSE) {
            pop[,,j] <- vaccinate(popslice=pop[,,j], vlookup=myvacc, type=vacc_program, mydate=theDate, params=fp)
          }
        }
        else {
          pop[,,j] <- vaccinate(popslice=pop[,,j], vlookup=myvacc, type=vacc_program, mydate=theDate, params=fp)
        }
      }
    } #end of aging IF
    LastMonth <- month(theDate)
    theDate<-theDate+7
    j <- j + 1
  } #end of date loop
  
  #label matrix iteration dimension with dates
  dimnames(pop)[[3]] <- as.Date(iterDate)
  return(pop)
}